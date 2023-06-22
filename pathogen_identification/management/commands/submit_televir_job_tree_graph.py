import os
from datetime import date
from typing import List

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from managing_files.models import ProcessControler
from pathogen_identification.models import (
    PIProject_Sample,
    Projects,
    ParameterSet,
    SoftwareTreeNode,
)
from pathogen_identification.utilities.tree_deployment import Tree_Progress
from pathogen_identification.utilities.utilities_pipeline import (
    Utility_Pipeline_Manager,
    Utils_Manager,
)
from pathogen_identification.utilities.utilities_pipeline import Parameter_DB_Utility
from utils.process_SGE import ProcessSGE
from pathogen_identification.constants_settings import ConstantsSettings as PICS
from constants.constants import Televir_Metadata_Constants


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
        pipeline_utils = Utility_Pipeline_Manager()
        parameter_db_util = Parameter_DB_Utility()
        samples = PIProject_Sample.objects.filter(project=project)

        try:
            for project_sample in samples:
                if not project_sample.is_deleted:
                    # existing parameter sets
                    existing_parameter_sets = ParameterSet.objects.filter(
                        project=project,
                        status=ParameterSet.STATUS_RUNNING,
                        sample__in=samples,
                    )
                    ## create a tree that contains all the parameter sets
                    (
                        combined_tree,
                        additional_leaves,
                    ) = utils.get_parameterset_leaves_list(existing_parameter_sets)
                    print(additional_leaves)

                    # CRUNCH TREE
                    module_tree = utils.module_tree(combined_tree, additional_leaves)

                    print("####################")
                    print(module_tree.node_index)

                    ## setup a deployment and record the progress
                    deployment_tree = Tree_Progress(
                        module_tree, project_sample, project
                    )

                    media_dir = os.path.join(
                        PICS.media_directory,
                        PICS.televir_subdirectory,
                        str(project_sample.project.owner.pk),
                        str(project_sample.project.pk),
                        str(project_sample.sample.pk),
                    )

                    stacked_df_path = media_dir + f"/{project_sample.pk}_stacked_df.tsv"
                    print(stacked_df_path)
                    stacked_df = deployment_tree.stacked_changes_log()
                    stacked_df.to_csv(
                        stacked_df_path,
                        sep="\t",
                        index=False,
                    )

                    ### create a graph
                    output_html_path = (
                        media_dir + f"/{project_sample.pk}_graph_output.html"
                    )

                    Rgraph_cmd = [
                        Televir_Metadata_Constants.BINARIES["software"]["R"]
                        + "/bin/"
                        + "Rscript",
                        "--vanilla",
                        "pathogen_identification/utilities/pipeline_dendrograph.R",
                        stacked_df_path,
                        output_html_path,
                    ]

                    print(" ".join(Rgraph_cmd))

                    result = os.system(" ".join(Rgraph_cmd))
                    print(result)

                    if os.path.exists(output_html_path):
                        print("graph created")
                    else:
                        print("graph not created")

                    ## read html file
                    def extract_graph_data(html_filepath) -> str:
                        """
                        graph data is stored in line cotaining: data-for="""

                        with open(html_filepath, "r") as f:
                            lines = f.readlines()

                        for line in lines:
                            if "data-for=" in line:
                                return line

                        return None

                    print(stacked_df)
                    graph_data = extract_graph_data(output_html_path)
                    print(graph_data)

        except Exception as e:
            print(e)
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_project(project_pk=project.pk),
                ProcessControler.FLAG_ERROR,
            )
            raise e
