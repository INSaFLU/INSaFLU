import os
from datetime import date

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from pathogen_identification.constants_settings import ConstantsSettings
from pathogen_identification.deployment_main import Run_Main_from_Leaf
from pathogen_identification.models import (
    PIProject_Sample,
    Projects,
    SoftwareTree,
    SoftwareTreeNode,
)
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager


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
        #
        user = User.objects.get(pk=options["user_id"])
        project = Projects.objects.get(pk=options["project_id"])

        utils = Utils_Manager(owner=user)

        technology = project.technology
        samples = PIProject_Sample.objects.filter(project=project)
        print(samples)

        pipeline_tree = utils.generate_software_tree(technology)
        pipeline_tree_index = utils.get_software_tree_index(technology)

        print("user pk: ", user.pk)
        print("project pk: ", project.pk)
        print("sample pk: ", samples[0].pk)
        print("pipeline_tree_index: ", pipeline_tree_index)

        local_tree = utils.generate_project_tree(technology, project)

        local_paths = local_tree.get_all_graph_paths_explicit()
        print("local paths: ", local_paths.keys())
        sample = samples[0]

        for sample in samples:
            print("sample: ", sample)

            for leaf, path in local_paths.items():
                print("leaf: ", leaf)

                matched_path = utils.utility_manager.match_path_to_tree(
                    path, pipeline_tree
                )

                matched_path_node = SoftwareTreeNode.objects.get(
                    software_tree__pk=pipeline_tree_index, index=matched_path
                )

                matched_path_pk = matched_path_node.pk

                print("matched path: ", leaf, matched_path)
                exists = utils.parameter_util.check_ParameterSet_exists(
                    sample=sample, leaf=matched_path_node, project=project
                )
                print(exists)
                if utils.parameter_util.check_ParameterSet_exists(
                    sample=sample, leaf=matched_path_node, project=project
                ):
                    print("parameter set exists")

                    if utils.parameter_util.check_ParameterSet_processed(
                        sample=sample, leaf=leaf, project=project
                    ):

                        continue

                # path_db = pipeline_tree.df_from_path(path)
                print("path_db: ", path)

                pipeline_tree_query = SoftwareTree.objects.get(pk=pipeline_tree_index)

                run = Run_Main_from_Leaf(
                    user=user,
                    input_data=sample,
                    project=project,
                    pipeline_leaf=matched_path_node,
                    pipeline_tree=pipeline_tree_query,
                    odir=options["outdir"],
                )
                run.Submit()
