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
        reference_file_path = options["reference"]
        outdir = options["outdir"]
        os.makedirs(outdir, exist_ok=True)
