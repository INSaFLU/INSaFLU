import os
from datetime import date
from typing import List

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from managing_files.models import ProcessControler
from pathogen_identification.models import PIProject_Sample, Projects, SoftwareTreeNode
from pathogen_identification.utilities.tree_deployment import (
    Tree_Progress,
    TreeProgressGraph,
)
from pathogen_identification.utilities.utilities_pipeline import (
    SoftwareTreeUtils,
    Utils_Manager,
)
from pathogen_identification.utilities.utilities_views import (
    RawReferenceUtils,
    set_control_reports,
)
from utils.process_SGE import ProcessSGE


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
            "-o",
            "--outdir",
            type=str,
            help="output directory",
        )

    def handle(self, *args, **options):
        ###
        # SETUP
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        user = User.objects.get(pk=options["user_id"])
        project = Projects.objects.get(pk=options["project_id"])
        technology = project.technology

        samples = PIProject_Sample.objects.filter(
            project=project, is_deleted=False, pk=options["sample_id"]
        )

        sample = samples.first()

        # PROCESS CONTROLER

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_project_sample(
                project_pk=project.pk, sample_pk=sample.pk
            ),
            ProcessControler.FLAG_RUNNING,
        )

        process = ProcessControler.objects.filter(
            owner__id=user.pk,
            name=process_controler.get_name_televir_project_sample(
                project_pk=project.pk, sample_pk=sample.pk
            ),
        )

        if process.exists():
            process = process.first()

            print(f"Process {process.pk} has been submitted")

        else:
            print("Process does not exist")

        # UTILITIES
        utils = Utils_Manager()
        software_utils = SoftwareTreeUtils(user, project)
        ####
        local_tree = software_utils.generate_project_tree()
        local_paths = local_tree.get_all_graph_paths_explicit()

        # pipeline_tree = utils.generate_software_tree(technology, tree_makeup)
        pipeline_tree = software_utils.generate_software_tree_extend(local_tree)
        pipeline_tree_index = local_tree.software_tree_pk

        # MANAGEMENT
        matched_paths = {
            leaf: utils.utility_manager.match_path_to_tree_safe(path, pipeline_tree)
            for leaf, path in local_paths.items()
        }

        matched_paths = {k: v for k, v in matched_paths.items() if v is not None}

        available_path_nodes = {
            leaf: SoftwareTreeNode.objects.get(
                software_tree__pk=pipeline_tree_index, index=path
            )
            for leaf, path in matched_paths.items()
        }

        available_path_nodes = {
            leaf: utils.parameter_util.check_ParameterSet_available_to_run(
                sample=sample, leaf=matched_path_node, project=project
            )
            for leaf, matched_path_node in available_path_nodes.items()
        }

        matched_paths = {
            k: v for k, v in matched_paths.items() if available_path_nodes[k] == True
        }

        # SUBMISSION
        print("MATCHED PATHS")
        print(matched_paths)

        module_tree = utils.module_tree(pipeline_tree, list(matched_paths.values()))

        try:
            for project_sample in samples:
                if project_sample.is_deleted:
                    continue
                if len(matched_paths) > 0:
                    graph_progress = TreeProgressGraph(project_sample)

                    deployment_tree = Tree_Progress(
                        module_tree, project_sample, project
                    )

                    graph_progress.generate_graph()

                    deployment_tree.cycle_process()

                    graph_progress.generate_graph()
                    set_control_reports(project.pk)

                    reference_utils = RawReferenceUtils(project_sample)
                    _ = reference_utils.create_compound_references()

                    _ = process_SGE.set_submit_televir_sort_pisample_reports(
                        user=user,
                        pisample_pk=project_sample.pk,
                    )

                    break

            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_project_sample(
                    project_pk=project.pk, sample_pk=sample.pk
                ),
                ProcessControler.FLAG_FINISHED,
            )

        except Exception as e:
            print(e)
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_project_sample(
                    project_pk=project.pk, sample_pk=sample.pk
                ),
                ProcessControler.FLAG_ERROR,
            )
            reference_utils = RawReferenceUtils(project_sample)
            _ = reference_utils.create_compound_references()

            raise e
