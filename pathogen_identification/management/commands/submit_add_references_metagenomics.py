import os
from datetime import date
from typing import List

import pandas as pd
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from constants.constants import Televir_Metadata_Constants
from pathogen_identification.models import PIProject_Sample, RunMain
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
            "-r",
            "--reference",
            type=int,
            help="sample id",
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
        reference_file_path = options["reference"]
        outdir = options["outdir"]
        os.makedirs(outdir, exist_ok=True)

        try:
            reference_df = pd.read_csv(reference_file_path, sep="\t")
        except FileNotFoundError:
            print("Reference file not found")
            return

        assert any(
            [
                "Description" in reference_df.columns,
                "accession" in reference_df.columns,
                "taxid" in reference_df.columns,
            ]
        ), "Reference file must contain at least one of the following columns: Description, accession, taxid"

        if sample_pk is None and project_pk is None:
            raise ValueError("Must provide either sample or project")

        proxy_run = RunMain()

        if project_pk is not None:
            proxy_run.project = PIProject_Sample.objects.get(pk=int(project_pk)).project

        if sample_pk is not None:
            proxy_run.sample = PIProject_Sample.objects.get(pk=int(sample_pk))
