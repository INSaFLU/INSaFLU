from typing import List

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from managing_files.models import ProcessControler
from pathogen_identification.utilities.reference_utils import file_reference_to_insaflu
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
            "--ref_id",
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
        ref_id = int(options["ref_id"])

        # PROCESS CONTROLER
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_file_televir_teleflu_ref_create(
                ref_id=ref_id,
            ),
            ProcessControler.FLAG_RUNNING,
        )

        process = ProcessControler.objects.filter(
            owner__id=user.pk,
            name=process_controler.get_name_file_televir_teleflu_ref_create(
                ref_id=ref_id,
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
            success, ref_id = file_reference_to_insaflu(ref_id, user_id)
            print(success, ref_id)
            if success is False:
                process_SGE.set_process_controler(
                    user,
                    process_controler.get_name_file_televir_teleflu_ref_create(
                        ref_id=ref_id,
                    ),
                    ProcessControler.FLAG_ERROR,
                )
                return

        except Exception as e:
            print(e)
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_file_televir_teleflu_ref_create(
                    ref_id=ref_id,
                ),
                ProcessControler.FLAG_ERROR,
            )
            raise e

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_file_televir_teleflu_ref_create(
                ref_id=ref_id,
            ),
            ProcessControler.FLAG_FINISHED,
        )
