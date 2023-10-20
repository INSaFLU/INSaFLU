"""
Created on January 30, 2023

@author: daniel.sobral
"""
import logging
import os

from django.contrib.auth.models import User
from django.core.management import BaseCommand
from django.db import transaction

from constants.constants import Constants
from managing_files.models import Sample
from pathogen_identification.models import PIProject_Sample, Projects
from pathogen_identification.utilities.utilities_pipeline import (
    SoftwareTreeUtils,
    Utils_Manager,
)
from settings.constants_settings import ConstantsSettings
from utils.process_SGE import ProcessSGE


class Command(BaseCommand):
    """
    classdocs
    """

    help = "Create a TELEVIR project and add a sample to it. Returns project id."

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
            "--sample_name", nargs="?", type=str, required=True, help="Sample Name"
        )
        parser.add_argument(
            "--user_login",
            nargs="?",
            type=str,
            required=True,
            help="User login of the project and sample owner",
        )
        parser.add_argument(
            "--test",
            action="store_true",
            required=False,
            default=False,
            help="Test if sample is ready for projects message",
        )

    # A command must define handle()

    def success_message(self, project_id: int):
        self.stdout.write(f"Project {project_id} submitted.")

    def handle(self, *args, **options):
        project_name = options["project_name"]
        sample_name = options["sample_name"]
        account = options["user_login"]

        user = User.objects.get(username=account)

        try:
            sample = Sample.objects.get(name=sample_name, owner=user, is_deleted=False)
        except Sample.DoesNotExist:
            self.stdout.write(f"Sample {sample_name} does not exist.")
            return False

        process_SGE = ProcessSGE()

        if options["test"]:
            self.success_message(1)

        try:
            try:
                project = Projects.objects.get(
                    name__iexact=project_name,
                    is_deleted=False,
                    owner__username=user.username,
                )
            except Projects.DoesNotExist:
                project = Projects()
                project.name = project_name
                project.owner = user

                # TODO Check where these constants are, or define them somewhere...
                technology = ConstantsSettings.TECHNOLOGY_minion
                if sample.type_of_fastq == Sample.TYPE_OF_FASTQ_illumina:
                    technology = ConstantsSettings.TECHNOLOGY_illumina
                project.technology = technology
                project.save()

            try:
                project_sample = PIProject_Sample.objects.get(
                    project=project,
                    sample=sample,
                )

            except PIProject_Sample.DoesNotExist:
                with transaction.atomic():
                    project_sample = PIProject_Sample()
                    project_sample.project = project
                    project_sample.sample = sample
                    project_sample.name = sample.name
                    project_sample_input = sample.file_name_1
                    if sample.is_valid_2:
                        project_sample_input += ";" + sample.file_name_2
                    project_sample.input = project_sample_input
                    sample_technology = ConstantsSettings.TECHNOLOGY_minion
                    if sample.type_of_fastq == Sample.TYPE_OF_FASTQ_illumina:
                        sample_technology = ConstantsSettings.TECHNOLOGY_illumina
                    if project.technology != sample_technology:
                        self.stdout.write(
                            "Project has different technology {} from sample technology {}...".format(
                                project.technology, sample_technology
                            )
                        )
                    project_sample.technology = sample.type_of_fastq
                    project_sample.report = "report"
                    project_sample.save()

            utils = Utils_Manager()
            utils = SoftwareTreeUtils(
                user=user,
                project=project,
            )
            runs_to_deploy = utils.check_runs_to_deploy_sample(project_sample)

            if len(runs_to_deploy) > 0:
                for proj_sample, leafs_to_deploy in runs_to_deploy.items():
                    for leaf in leafs_to_deploy:
                        taskID = process_SGE.set_submit_televir_run(
                            user=user,
                            project_pk=project.pk,
                            sample_pk=proj_sample.pk,
                            leaf_pk=leaf.pk,
                        )

            self.success_message(project.pk)

        except User.DoesNotExist as e:
            self.stdout.write("Error: User '{}' does not exist.".format(account))
        except Exception as e:
            self.stdout.write("Error: {}.".format(e))
