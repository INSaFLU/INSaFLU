import os
from typing import List

import pandas as pd
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from managing_files.models import ProcessControler
from pathogen_identification.constants_settings import ConstantsSettings as CS
from pathogen_identification.models import FinalReport, Projects
from pathogen_identification.utilities.explify_merge import (
    get_illumina_found,
    merge_panels,
    process_televir,
    read_panel,
)
from pathogen_identification.utilities.utilities_general import get_project_dir
from utils.process_SGE import ProcessSGE


class Command(BaseCommand):
    help = "deploy run"

    def add_arguments(self, parser):
        parser.add_argument(
            "--project_id",
            type=int,
            help="televir project to be merged (pk)",
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

        project = Projects.objects.filter(
            is_deleted=False, pk=options["project_id"]
        ).first()

        user = project.owner

        output_file_merged = os.path.join(
            get_project_dir(project), CS.EXPLIFY_MERGE_SUFFIX + f".{project.pk}.tsv"
        )

        # PROCESS CONTROLER

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_project_merge_explify(
                project_pk=project.pk,
            ),
            ProcessControler.FLAG_RUNNING,
        )

        all_reports = FinalReport.objects.filter(
            run__project__pk=int(project.pk)
        ).order_by("-coverage")
        output_dir = options["outdir"]

        televir_reports = pd.DataFrame.from_records(all_reports.values())

        try:
            print(output_file_merged)
            rpip_panel = read_panel(options["rpip"], panel="Microorganisms (RPIP)")
            upip_panel = read_panel(options["upip"], panel="Microorganisms (UPIP)")
            print("PANELS READ")

            illumina_found = get_illumina_found(
                [rpip_panel, upip_panel], tmp_dir=output_dir
            )
            print("ILLUMINA FOUND")
            telebac_found = process_televir(televir_reports)
            print("TELEBAC FOUND")

            merged_panel = merge_panels(illumina_found, telebac_found)
            merged_panel.to_csv(output_file_merged, sep="\t", index=False)

        except Exception as e:
            print(e)
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_project_merge_explify(
                    project_pk=project.pk,
                ),
                ProcessControler.FLAG_ERROR,
            )
            return

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_project_merge_explify(
                project_pk=project.pk,
            ),
            ProcessControler.FLAG_FINISHED,
        )
