import os
from datetime import date
from typing import List

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from managing_files.models import ProcessControler
from pathogen_identification.models import (
    PIProject_Sample,
    Projects,
    SoftwareTree,
    SoftwareTreeNode,
)
from pathogen_identification.utilities.tree_deployment import (
    Tree_Progress,
    TreeProgressGraph,
)
from pathogen_identification.utilities.utilities_pipeline import (
    Utility_Pipeline_Manager,
    Utils_Manager,
)
from utils.process_SGE import ProcessSGE
from pathogen_identification.utilities.utilities_views import set_control_reports

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
            "-o",
            "--outdir",
            type=str,
            help="output directory",
        )

    def handle(self, *args, **options):
        ###
        # SETUP

        user = User.objects.get(pk=options["user_id"])
        project = Projects.objects.get(pk=options["project_id"])
        technology = project.technology

        # PROCESS CONTROLER
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_project(project_pk=project.pk),
            ProcessControler.FLAG_RUNNING,
        )

        process = ProcessControler.objects.filter(
            owner__id=user.pk,
            name=process_controler.get_name_televir_project(project_pk=project.pk),
        )

        if process.exists():
            process = process.first()

            print(f"Process {process.pk} has been submitted")

        else:
            print("Process does not exist")

        # UTILITIES
        utils = Utils_Manager()
        samples = PIProject_Sample.objects.filter(project=project)
        local_tree = utils.generate_project_tree(technology, project, user)
        local_paths = local_tree.get_all_graph_paths_explicit()

        tree_makeup = local_tree.makeup

        #pipeline_tree = utils.generate_software_tree(technology, tree_makeup)
        pipeline_tree= utils.generate_software_tree_extend(local_tree)
        pipeline_tree_index = utils.get_software_tree_index(technology, tree_makeup)
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

        try:
            for project_sample in samples:
                if not project_sample.is_deleted:
                    available_path_nodes_sample = {
                        leaf: utils.parameter_util.check_ParameterSet_available_to_run(
                            sample=project_sample,
                            leaf=matched_path_node,
                            project=project,
                        )
                        for leaf, matched_path_node in available_path_nodes.items()
                    }

                    matched_paths_sample = {
                        k: v
                        for k, v in matched_paths.items()
                        if available_path_nodes_sample[k]
                    }
                    ###

                    # SUBMISSION

                    pipeline_utils = Utility_Pipeline_Manager()

                    reduced_tree = utils.tree_subset(
                        pipeline_tree, list(matched_paths_sample.values())
                    )

                    module_tree = pipeline_utils.compress_software_tree(reduced_tree)

                    graph_progress = TreeProgressGraph(project_sample)

                    deployment_tree = Tree_Progress(
                        module_tree, project_sample, project
                    )

                    graph_progress.generate_graph()

                    deployment_tree.cycle_process()

                    graph_progress.generate_graph()

                    ### set control reports
                    set_control_reports(project.pk)

            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_project(project_pk=project.pk),
                ProcessControler.FLAG_FINISHED,
            )

        except Exception as e:
            print(e)
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_project(project_pk=project.pk),
                ProcessControler.FLAG_ERROR,
            )
            raise e
