import os
from datetime import date
from typing import List

import pandas as pd
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from constants.constants import Televir_Metadata_Constants
from pathogen_identification.models import FinalReport, PIProject_Sample
from pathogen_identification.modules.object_classes import RunCMD
from settings.constants_settings import ConstantsSettings


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
            "--sample_id",
            type=int,
            help="televir sample to be run (pk)",
            default=None,
        )

        parser.add_argument(
            "-p",
            "--project_id",
            type=int,
            help="project id",
            default=None,
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

        sample_pk = options["sample_id"]
        project_pk = options["project_id"]
        outdir = options["outdir"]
        os.makedirs(outdir, exist_ok=True)

        if sample_pk is None and project_pk is None:
            raise ValueError("Must provide either sample or project")

        if project_pk is not None:
            reports = FinalReport.objects.filter(run__project_pk=int(project_pk))
        else:
            reports = FinalReport.objects.filter(sample__pk=int(sample_pk))

        env_bin = os.path.join(
            Televir_Metadata_Constants.BINARIES["ROOT"],
            Televir_Metadata_Constants.BINARIES[
                ConstantsSettings.PIPELINE_NAME_remapping
            ]["default"],
            "bin",
        )

        print("env_bin", env_bin)

        cmd_runner = RunCMD(
            bin=env_bin,
            logdir=outdir,
            prefix=f"sample_{sample_pk}",
            task="update_stats",
        )

        for report in reports:
            bam_path = report.bam_path
            stats_report = os.path.join(
                outdir,
                f"sample_{sample_pk}_report_{report.pk}.tsv",
            )

            if not os.path.exists(bam_path):
                print(f"Skipping {bam_path}")
                continue

            cmd = [
                "samtools",
                "stats",
                bam_path,
                "|",
                "grep ^SN",
                "|",
                "cut -f 2-",
                ">",
                stats_report,
            ]

            cmd_runner.run_script_software(cmd)
            stats_df = pd.read_csv(
                stats_report, sep="\t", header=None, index_col=0
            ).rename(columns={0: "stat", 1: "value", 2: "comment"})

            error_rate = stats_df.loc["error rate:", "value"]
            quality_avg = stats_df.loc["average quality:", "value"]

            report.error_rate = int(error_rate)
            report.quality_avg = int(quality_avg)
            report.save()
