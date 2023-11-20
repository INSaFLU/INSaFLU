import os
from typing import List

import pandas as pd
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from managing_files.models import ProcessControler
from pathogen_identification.constants_settings import ConstantsSettings as CS
from pathogen_identification.models import FinalReport, PIProject_Sample, Projects
from pathogen_identification.templatetags.report_colors import flag_false_positive
from pathogen_identification.utilities.explify_merge import (
    get_illumina_found,
    merge_panels,
    process_televir,
    read_panel,
)
from pathogen_identification.utilities.utilities_general import get_services_dir
from utils.process_SGE import ProcessSGE


class Command(BaseCommand):
    help = "deploy run"

    def add_arguments(self, parser):
        parser.add_argument(
            "--user_id",
            type=int,
            help="user id",
        )

        parser.add_argument(
            "--televir",
            type=str,
            help="televir panel",
        )

        parser.add_argument(
            "--rpip",
            type=str,
            help="RPIP panel",
        )

        parser.add_argument(
            "--upip",
            type=str,
            help="UPIP panel",
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
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        user_pk = options["user_id"]
        user = User.objects.get(pk=user_pk)

        output_file_merged = os.path.join(
            get_services_dir(user), "merged_televir_explify.tsv"
        )

        # PROCESS CONTROLER

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_project_merge_explify_external(
                user_pk == user.pk,
            ),
            ProcessControler.FLAG_RUNNING,
        )

        output_dir = options["outdir"]

        try:
            rpip_panel = read_panel(options["rpip"], panel="Microorganisms (RPIP)")
            upip_panel = read_panel(options["upip"], panel="Microorganisms (UPIP)")
            televir_reports = pd.read_csv(options["televir"], sep="\t")

            illumina_found = get_illumina_found(
                [rpip_panel, upip_panel], tmp_dir=output_dir
            )
            telebac_found = process_televir(televir_reports)

            merged_panel = merge_panels(illumina_found, telebac_found)
            merged_panel.to_csv(output_file_merged, sep="\t", index=False)

        except Exception as e:
            print(e)
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_project_merge_explify_external(
                    user_pk=user.pk,
                ),
                ProcessControler.FLAG_ERROR,
            )
            return

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_project_merge_explify_external(
                user_pk=user.pk,
            ),
            ProcessControler.FLAG_FINISHED,
        )
