import os
from datetime import date
from typing import List

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from managing_files.models import ProcessControler
from pathogen_identification.models import TelefluMapping, TeleFluProject
from pathogen_identification.utilities.reference_utils import create_televir_igv_report
from utils.process_SGE import ProcessSGE


class Command(BaseCommand):
    help = "deploy televir sample reports sort"

    def add_arguments(self, parser):

        parser.add_argument(
            "--leaf_id",
            type=int,
            help="mappings (pk)",
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
        project_id = int(options["project_id"])
        project = TeleFluProject.objects.get(pk=project_id)
        user = User.objects.get(pk=project.owner.pk)
        ref_id = int(options["leaf_id"])

        # PROCESS CONTROLER
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_teleflu_igv_stack(
                teleflu_mapping_id=ref_id,
            ),
            ProcessControler.FLAG_RUNNING,
        )

        teleflu_project = TeleFluProject.objects.get(pk=project_id)
        teleflu_mapping = TelefluMapping.objects.get(
            leaf__pk=ref_id, teleflu_project__pk=project_id
        )

        process = ProcessControler.objects.filter(
            owner__id=user.pk,
            name=process_controler.get_name_televir_teleflu_igv_stack(
                teleflu_mapping_id=ref_id,
            ),
            is_finished=False,
        )

        if process.exists():
            process = process.first()

            print(f"Process {process.pk} has been submitted")

        else:
            print("Process does not exist")

        # UTILITIES
        try:

            create_televir_igv_report(
                teleflu_project_pk=teleflu_project.pk,
                leaf_index=teleflu_mapping.leaf.pk,
            )

        except Exception as e:
            print(e)
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_teleflu_igv_stack(
                    teleflu_mapping_id=ref_id,
                ),
                ProcessControler.FLAG_ERROR,
            )
            raise e

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_teleflu_igv_stack(
                teleflu_mapping_id=ref_id,
            ),
            ProcessControler.FLAG_FINISHED,
        )
