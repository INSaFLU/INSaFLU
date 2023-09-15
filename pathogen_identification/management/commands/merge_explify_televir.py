from typing import List

import pandas as pd
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from managing_files.models import ProcessControler
from pathogen_identification.models import (
    FinalReport,
    Projects,
    SoftwareTree,
    SoftwareTreeNode,
)
from pathogen_identification.utilities.explify_merge import (
    get_illumina_found,
    merge_panels,
    process_televir,
    read_panel,
)
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

        outdir = options["outdir"]
        output_file_merged = outdir + f"merged_explify_project.{project.pk}.tsv"

        # PROCESS CONTROLER

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_project_merge_explify(
                project_pk=project.pk,
            ),
            ProcessControler.FLAG_RUNNING,
        )

        process = ProcessControler.objects.filter(
            owner__id=user.pk,
            name=process_controler.get_name_televir_project_merge_explify(
                project_pk=project.pk,
            ),
        )

        all_reports = FinalReport.objects.filter(
            run__project__pk=int(project.pk)
        ).order_by("-coverage")

        televir_reports = pd.DataFrame.from_records(all_reports.values())

        try:
            rpip_panel = read_panel(
                options["rpip_panel"], panel="Microorganisms (RPIP)"
            )
            upip_panel = read_panel(
                options["upip_panel"], panel="Microorganisms (UPIP)"
            )
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

        illumina_found = get_illumina_found([rpip_panel, upip_panel])
        telebac_found = process_televir(televir_reports)

        merged_panel = merge_panels(illumina_found, telebac_found)
        merged_panel.to_csv(output_file_merged, sep="\t", index=False)

        if process.exists():
            process = process.first()

            print(f"Process {process.pk} has been submitted")

        else:
            print("Process does not exist")
