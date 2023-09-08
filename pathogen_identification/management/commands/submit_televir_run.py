import os
from datetime import date
from typing import List

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from managing_files.models import ProcessControler
from pathogen_identification.constants_settings import ConstantsSettings
from pathogen_identification.deployment_main import Run_Main_from_Leaf
from pathogen_identification.models import (ParameterSet, PIProject_Sample,
                                            Projects, SoftwareTree,
                                            SoftwareTreeNode)
from pathogen_identification.utilities.tree_deployment import TreeProgressGraph
from pathogen_identification.utilities.utilities_pipeline import (
    SoftwareTreeUtils, Utils_Manager)
from pathogen_identification.utilities.utilities_views import \
    set_control_reports
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
            "--project_id",
            type=int,
            help="project to be run (pk)",
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
            "-o",
            "--outdir",
            type=str,
            help="output directory",
        )

    def handle(self, *args, **options):
        ###
        #### SETUP

        user = User.objects.get(pk=options["user_id"])
        project = Projects.objects.get(pk=options["project_id"])
        target_sample = PIProject_Sample.objects.get(pk=options["sample_id"])
        leaf_index = options["leaf_id"]
        matched_path_node = SoftwareTreeNode.objects.get(pk=leaf_index)

        technology = project.technology

        ### PROCESS CONTROLER
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_run(
                project_pk=project.pk,
                sample_pk=target_sample.pk,
                leaf_pk=matched_path_node.pk,
            ),
            ProcessControler.FLAG_RUNNING,
        )

        ### UTILITIES
        utils = Utils_Manager()
        software_utils= SoftwareTreeUtils(user, project)

        local_tree = software_utils.generate_project_tree()

        #tree_makeup = local_tree.makeup
        #pipeline_tree= utils.generate_software_tree_extend(local_tree, user)
        pipeline_tree_index = local_tree.software_tree_pk
        pipeline_tree_query = SoftwareTree.objects.get(pk=pipeline_tree_index)

        ### MANAGEMENT
        submission_dict = {target_sample: []}

        was_run_killed = utils.parameter_util.check_ParameterSet_killed(
            sample=target_sample, leaf=matched_path_node, project=project
        )

        print("was_run_killed", was_run_killed)


        ### draw graph
        graph_progress = TreeProgressGraph(target_sample)


        ### SUBMISSION
        try:

            if was_run_killed:

                utils.parameter_util.parameterset_update_status(
                    sample=target_sample,
                    leaf=matched_path_node,
                    project=project,
                    status=ParameterSet.STATUS_NOT_STARTED,
                )

            else:

                if (
                    utils.parameter_util.check_ParameterSet_available_to_run(
                        sample=target_sample, leaf=matched_path_node, project=project
                    )
                    is False
                ):
                    raise Exception("ParameterSet not available")

                run = Run_Main_from_Leaf(
                    user=user,
                    input_data=target_sample,
                    project=project,
                    pipeline_leaf=matched_path_node,
                    pipeline_tree=pipeline_tree_query,
                    odir=options["outdir"],
                    threads=ConstantsSettings.DEPLOYMENT_THREADS,
                )

                if run.is_available:
                    run.get_in_line()
                    submission_dict[target_sample].append(run)

                for sample, runs in submission_dict.items():
                    for run in runs:

                        run.Submit()
                
                graph_progress.generate_graph()
                set_control_reports(project.pk)

            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_run(
                    project_pk=project.pk,
                    sample_pk=target_sample.pk,
                    leaf_pk=matched_path_node.pk,
                ),
                ProcessControler.FLAG_FINISHED,
            )

        except Exception as e:
            print(e)
            graph_progress.generate_graph()

            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_run(
                    project_pk=project.pk,
                    sample_pk=target_sample.pk,
                    leaf_pk=matched_path_node.pk,
                ),
                ProcessControler.FLAG_ERROR,
            )
            raise e
