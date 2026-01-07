import os
from datetime import date
from typing import List

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from pathogen_identification.models import PIProject_Sample
from managing_files.models import ProcessControler

from pathogen_identification.utilities.utilities_views import calculate_reports_overlaps
from utils.process_SGE import ProcessSGE


class Command(BaseCommand):
    help = "deploy televir sample reports sort"

    def add_arguments(self, parser):
        parser.add_argument(
            "--user_id",
            type=int,
            help="user deploying the run (pk)",
        )

        parser.add_argument(
            "--pisample_id",
            type=int,
            help="project to be run (pk)",
        )

        parser.add_argument(
            "-o",
            "--outdir",
            type=str,
            help="output directory",
        )
        parser.add_argument(
            "-f",
            "--force",
            action="store_true",
            help="force the run",
        )

    def handle(self, *args, **options):
        ###
        # SETUP

        user = User.objects.get(pk=options["user_id"])
        project_sample = PIProject_Sample.objects.get(pk=options["pisample_id"])
        force = options["force"]
        # PROCESS CONTROLER
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_project_sample_sort(
                sample_pk=project_sample.pk
            ),
            ProcessControler.FLAG_RUNNING,
        )

        process = ProcessControler.objects.filter(
            owner__id=user.pk,
            name=process_controler.get_name_televir_project_sample_sort(
                sample_pk=project_sample.pk
            ),
        )

        if process.exists():
            process = process.first()

            print(f"Process {process.pk} has been submitted")

        else:
            print("Process does not exist")

        # UTILITIES
        try:
            if not project_sample.is_deleted:
                calculate_reports_overlaps(project_sample, force=True)

        except Exception as e:
            print(e)
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_project_sample_sort(
                    sample_pk=project_sample.pk
                ),
                ProcessControler.FLAG_ERROR,
            )
            raise e

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_project_sample_sort(
                sample_pk=project_sample.pk
            ),
            ProcessControler.FLAG_FINISHED,
        )
