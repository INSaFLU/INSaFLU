import os
from datetime import date
from typing import List

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from managing_files.models import ProcessControler
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
        technology = project.technology

        ### PROCESS CONTROLER
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_project(project_pk=project.pk),
            ProcessControler.FLAG_RUNNING,
        )

        ### UTILITIES
        utils = Utils_Manager()
        samples = PIProject_Sample.objects.filter(project=project)
        local_tree = utils.generate_project_tree(technology, project, user)
        local_paths = local_tree.get_all_graph_paths_explicit()

        tree_makeup = local_tree.makeup

        pipeline_tree = utils.generate_software_tree(technology, tree_makeup)
        pipeline_tree_index = utils.get_software_tree_index(technology, tree_makeup)
        pipeline_tree_query = SoftwareTree.objects.get(pk=pipeline_tree_index)

        ### MANAGEMENT
        submission_dict = {sample: [] for sample in samples if not sample.is_deleted}
        matched_paths = {
            leaf: utils.utility_manager.match_path_to_tree_safe(path, pipeline_tree)
            for leaf, path in local_paths.items()
        }
        available_paths = {
            leaf: path for leaf, path in matched_paths.items() if path is not None
        }

        available_path_nodes = {
            leaf: SoftwareTreeNode.objects.get(
                software_tree__pk=pipeline_tree_index, index=path
            )
            for leaf, path in available_paths.items()
        }

        ### SUBMISSION

        try:

            for sample in samples:

                for leaf, matched_path_node in available_path_nodes.items():

                    if not utils.parameter_util.check_ParameterSet_available(
                        sample=sample, leaf=matched_path_node, project=project
                    ):
                        continue

                    run = Run_Main_from_Leaf(
                        user=user,
                        input_data=sample,
                        project=project,
                        pipeline_leaf=matched_path_node,
                        pipeline_tree=pipeline_tree_query,
                        odir=options["outdir"],
                        threads=ConstantsSettings.DEPLOYMENT_THREADS,
                    )

                    if run.is_available:
                        run.get_in_line()
                        submission_dict[sample].append(run)

            for sample, runs in submission_dict.items():
                for run in runs:

                    run.Submit()

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
