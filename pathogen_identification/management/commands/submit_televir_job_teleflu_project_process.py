import os
from datetime import date
from typing import List

from django.contrib import messages
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from constants.meta_key_and_values import MetaKeyAndValue
from extend_user.models import Profile
from managing_files.manage_database import ManageDatabase
from managing_files.models import ProcessControler
from managing_files.models import Project as InsafluProject
from managing_files.models import ProjectSample
from pathogen_identification.models import TeleFluProject, TeleFluSample
from pathogen_identification.utilities.reference_utils import (
    create_teleflu_igv_report,
    teleflu_to_insaflu_reference,
)
from settings.default_software_project_sample import DefaultProjectSoftware
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
        user_id = int(options["user_id"])
        user = User.objects.get(pk=user_id)
        project_id = int(options["project_id"])

        # PROCESS CONTROLER
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_teleflu_project_process(
                project_id=project_id,
            ),
            ProcessControler.FLAG_RUNNING,
        )

        process = ProcessControler.objects.filter(
            owner__id=user.pk,
            name=process_controler.get_name_televir_teleflu_project_process(
                project_id=project_id,
            ),
            is_finished=False,
        )

        if process.exists():
            process = process.first()

            print(f"Process {process.pk} has been submitted")

        else:
            print("Process does not exist")

        # UTILITIES
        ### test anonymous account
        try:
            profile = Profile.objects.get(user=user)

        except Profile.DoesNotExist:
            raise Exception("User without profile")

        try:
            success = create_teleflu_igv_report(project_id)

            if success is False:
                process_SGE.set_process_controler(
                    user,
                    process_controler.get_name_televir_teleflu_project_process(
                        project_id=project_id,
                    ),
                    ProcessControler.FLAG_ERROR,
                )
                return

        except Exception as e:
            print(e)
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_teleflu_project_process(
                    project_id=project_id,
                ),
                ProcessControler.FLAG_ERROR,
            )
            raise e

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_teleflu_project_process(
                project_id=project_id,
            ),
            ProcessControler.FLAG_FINISHED,
        )
