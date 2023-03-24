"""
Created on January 30, 2023

@author: daniel.sobral
"""
import os
import logging
import pandas as pd
from django.core.management import BaseCommand
from django.contrib.auth.models import User
from django.db import transaction
from managing_files.models import Sample
from pathogen_identification.models import (
    Projects,
    RunMain,
    ParameterSet,
    FinalReport,
    PIProject_Sample,
)
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager
from settings.constants_settings import ConstantsSettings
from pathogen_identification.constants_settings import (
    ConstantsSettings as ConstantsSettingsPI,
)
from fluwebvirus.settings import MEDIA_ROOT
import pandas as pd


class Command(BaseCommand):
    """
    classdocs
    """

    help = "Checks if a given televir project finished."

    # logging
    logger_debug = logging.getLogger("fluWebVirus.debug")
    logger_production = logging.getLogger("fluWebVirus.production")

    def __init__(self, *args, **kwargs):
        super(Command, self).__init__(*args, **kwargs)

    def add_arguments(self, parser):
        parser.add_argument(
            "--project_name", nargs="?", type=str, required=True, help="Project Name"
        )
        parser.add_argument(
            "--user_login",
            nargs="?",
            type=str,
            required=True,
            help="User login of the project and sample owner",
        )
        # parser.add_argument('--report_file', nargs='?', type=str, required=True,
        # 					help='Path for the report file to be saved (if project finished)')

    # A command must define handle()
    def handle(self, *args, **options):

        project_name = options["project_name"]
        account = options["user_login"]
        # report_file = options['report_file']

        try:

            user = User.objects.get(username=account)

            project = Projects.objects.get(
                name__iexact=project_name,
                is_deleted=False,
                owner__username=user.username,
            )

            project_media_dir = os.path.join(
                MEDIA_ROOT,
                ConstantsSettingsPI.televir_subdirectory,
                str(user.pk),
                str(project.pk),
            )

            project_results_path = os.path.join(project_media_dir, "results.tsv")

            runs = RunMain.objects.filter(
                project=project, parameter_set__status=ParameterSet.STATUS_FINISHED
            )

            final_report = FinalReport.objects.filter(run__in=runs)

            reports_to_pandas = pd.DataFrame(final_report.values())
            reports_to_pandas["project_name"] = project.name

            def get_sample_name(sample_id):
                sample = PIProject_Sample.objects.get(id=sample_id)
                return sample.name

            def get_leaf_id(run_id):
                run = RunMain.objects.get(id=run_id)
                return run.parameter_set.leaf.pk

            reports_to_pandas["sample_name"] = reports_to_pandas.apply(
                lambda row: get_sample_name(row["sample_id"]), axis=1
            )

            reports_to_pandas["leaf_id"] = reports_to_pandas.apply(
                lambda row: get_leaf_id(row["run_id"]), axis=1
            )

            # if not reports_to_pandas.empty:
            reports_to_pandas.to_csv(
                project_results_path, sep="\t", index=False, header=True
            )

            self.stdout.write(project_results_path)

        except User.DoesNotExist as e:
            self.stdout.write("Error: User '{}' does not exist.".format(account))
        except Exception as e:
            self.stdout.write("Error: {}.".format(e))
