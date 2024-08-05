import os
from datetime import date
from typing import Dict, List

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from managing_files.models import ProcessControler
from pathogen_identification.constants_settings import ConstantsSettings
from pathogen_identification.deployment_main import Run_Main_from_Leaf
from pathogen_identification.models import (
    ParameterSet,
    PIProject_Sample,
    RunMain,
    SoftwareTree,
    SoftwareTreeNode,
)
from pathogen_identification.utilities.tree_deployment import TreeProgressGraph
from pathogen_identification.utilities.utilities_pipeline import (
    SoftwareTreeUtils,
    Utils_Manager,
)
from pathogen_identification.utilities.utilities_views import (
    RawReferenceUtils,
    set_control_reports,
)
from utils.process_SGE import ProcessSGE


class Sample_Staging:
    """
    Class to stage samples for a project
    """

    def __init__(self, sample: PIProject_Sample):
        self.sample = sample
        self.is_deleted = self.sample.is_deleted


class Command(BaseCommand):
    help = "deploy run"

    def add_arguments(self, parser):
        parser.add_argument(
            "--user_id",
            type=int,
            help="user deploying the run (pk)",
        )

        parser.add_argument(
            "--sample_id",
            type=int,
            help="televir sample to be run (pk)",
        )

        parser.add_argument(
            "--leaf_id",
            type=int,
            help="televir leaf to be run on sample (pk)",
        )

        parser.add_argument(
            "--combined_analysis",
            action="store_true",
            help="run combined analysis",
            default=False,
        )

        parser.add_argument(
            "--mapping_request",
            action="store_true",
            help="run mapping only",
            default=False,
        )

        parser.add_argument(
            "--mapping_run_id",
            type=int,
            help="mapping run to be used (pk)",
            required=False,
        )

        parser.add_argument(
            "-o",
            "--outdir",
            type=str,
            help="output directory",
        )

    def handle(self, *args, **options):
        ###
        #### SETUP

        user = User.objects.get(pk=options["user_id"])
        target_sample = PIProject_Sample.objects.get(pk=options["sample_id"])
        project = target_sample.project

        metagenomics = False
        mapping_only = False
        screening = False

        leaf_index = options["leaf_id"]
        combined_analysis = options["combined_analysis"]
        mapping_request = options["mapping_request"]
        mapping_run_pk = options["mapping_run_id"]

        if mapping_request:
            mapping_only = True
            if mapping_run_pk is None:
                raise Exception("mapping_run_id is required for mapping request")
        elif combined_analysis:
            metagenomics = True
        else:
            screening = True

        matched_path_node = SoftwareTreeNode.objects.get(pk=leaf_index)

        ### PROCESS CONTROLER
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_project_sample_metagenomics_run(
                sample_pk=target_sample.pk,
                leaf_pk=matched_path_node.pk,
            ),
            ProcessControler.FLAG_RUNNING,
        )

        ### UTILITIES
        utils = Utils_Manager()
        software_utils = SoftwareTreeUtils(user, project, sample=target_sample)

        local_tree = software_utils.generate_software_tree_safe(
            project,
            sample=target_sample,
            metagenomics=metagenomics,
            screening=screening,
            mapping_only=mapping_only,
        )

        pipeline_tree_index = local_tree.software_tree_pk
        pipeline_tree_query = SoftwareTree.objects.get(pk=pipeline_tree_index)

        ### MANAGEMENT
        submission_dict: Dict[PIProject_Sample, List[Run_Main_from_Leaf]] = {
            target_sample: []
        }

        was_run_killed = utils.parameter_util.check_ParameterSet_killed(
            sample=target_sample, leaf=matched_path_node, project=project
        )

        #### Deployment RUn

        mapping_run = RunMain.objects.get(pk=mapping_run_pk)

        ### SUBMISSION
        try:
            if was_run_killed:
                utils.parameter_util.parameterset_update_status(
                    sample=target_sample,
                    leaf=matched_path_node,
                    project=project,
                    status=ParameterSet.STATUS_NOT_STARTED,
                )

            ### to remove condition for production
            elif mapping_run.status not in [
                RunMain.STATUS_FINISHED,
                # RunMain.STATUS_ERROR,
            ]:
                run = Run_Main_from_Leaf(
                    user=user,
                    input_data=target_sample,
                    project=project,
                    pipeline_leaf=matched_path_node,
                    pipeline_tree=pipeline_tree_query,
                    odir=options["outdir"],
                    threads=ConstantsSettings.DEPLOYMENT_THREADS,
                    combined_analysis=combined_analysis,
                    mapping_request=mapping_request,
                    run_pk=mapping_run_pk,
                )

                run.is_available = True
                run.set_to_queued()

                submission_dict[target_sample].append(run)

                for sample, runs in submission_dict.items():
                    for run in runs:
                        run.Submit()

                # graph_progress.generate_graph()
                set_control_reports(project.pk)

            reference_utils = RawReferenceUtils(target_sample)
            _ = reference_utils.create_compound_references()

            _ = process_SGE.set_submit_televir_sort_pisample_reports(
                user=user,
                pisample_pk=target_sample.pk,
            )
            # calculate_reports_overlaps(target_sample, force=True)

            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_project_sample_metagenomics_run(
                    sample_pk=target_sample.pk,
                    leaf_pk=matched_path_node.pk,
                ),
                ProcessControler.FLAG_FINISHED,
            )

        except Exception as e:
            print(e)
            # graph_progress.generate_graph()

            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_project_sample_metagenomics_run(
                    sample_pk=target_sample.pk,
                    leaf_pk=matched_path_node.pk,
                ),
                ProcessControler.FLAG_ERROR,
            )

            reference_utils = RawReferenceUtils(target_sample)
            _ = reference_utils.create_compound_references()

            raise e
