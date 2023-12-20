import os
from datetime import date
from typing import List

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from constants.constants import Televir_Metadata_Constants
from managing_files.models import ProcessControler
from pathogen_identification.constants_settings import ConstantsSettings as PICS
from pathogen_identification.models import PIProject_Sample, Projects
from pathogen_identification.utilities.tree_deployment import TreeProgressGraph
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

        samples = PIProject_Sample.objects.filter(project=project)

        try:
            for project_sample in samples:
                print(f"########### {project_sample.sample.name} ###########")
                if not project_sample.is_deleted:
                    graph_progress = TreeProgressGraph(project_sample)
                    graph_progress.generate_graph()

                    graph_data, _ = graph_progress.get_graph_data()
                    import re

                    # print(graph_data)
                    # print(graph_data)

        except Exception as e:
            print(e)
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_project(project_pk=project.pk),
                ProcessControler.FLAG_ERROR,
            )
            raise e
