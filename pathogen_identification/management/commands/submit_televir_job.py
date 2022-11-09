import os
from datetime import date

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from pathogen_identification.constants_settings import ConstantsSettings
from pathogen_identification.deployment_main import Run_Main_from_Leaf
from pathogen_identification.models import (
    ParameterSet,
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

        utils = Utils_Manager()

        technology = project.technology
        samples = PIProject_Sample.objects.filter(project=project)
        local_tree = utils.generate_project_tree(technology, project, user)
        tree_makeup = local_tree.makeup

        print(samples)

        pipeline_tree = utils.generate_software_tree(technology, tree_makeup)
        pipeline_tree_index = utils.get_software_tree_index(technology, tree_makeup)

        print("user pk: ", user.pk)
        print("project pk: ", project.pk)
        print("sample pk: ", samples[0].pk)
        print("pipeline_tree_index: ", pipeline_tree_index)

        local_paths = local_tree.get_all_graph_paths_explicit()
        print("local paths: ", local_paths.keys())
        sample = samples[0]

        submission_dict = {sample: [] for sample in samples}

        for sample in samples:
            print("sample: ", sample)

            if sample.is_deleted:
                continue

            for leaf, path in local_paths.items():
                print("leaf: ", leaf)
                print("path: ", path)

                matched_path = utils.utility_manager.match_path_to_tree(
                    path, pipeline_tree
                )

                matched_path_node = SoftwareTreeNode.objects.get(
                    software_tree__pk=pipeline_tree_index, index=matched_path
                )

                print("matched path: ", leaf, matched_path)
                exists = utils.parameter_util.check_ParameterSet_exists(
                    sample=sample, leaf=matched_path_node, project=project
                )
                print(exists)
                if exists:
                    print("parameter set exists")

                    if utils.parameter_util.check_ParameterSet_processed(
                        sample=sample, leaf=leaf, project=project
                    ):
                        print("parameter set processed")

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

                if run.is_available:
                    run.get_in_line()
                    submission_dict[sample].append(run)

        for sample, runs in submission_dict.items():
            for run in runs:
                # if run.get_status() == ParameterSet.STATUS_ERROR:

                run.Submit()
