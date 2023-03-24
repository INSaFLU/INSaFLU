"""
Created on January 30, 2023

@author: daniel.sobral
"""
import os
import logging
from django.core.management import BaseCommand
from django.contrib.auth.models import User
from django.db import transaction
from managing_files.models import Sample
from pathogen_identification.models import Projects, PIProject_Sample
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager
from settings.constants_settings import ConstantsSettings
from utils.process_SGE import ProcessSGE
import os

from django.conf import settings
from extend_user.models import Profile
from managing_files.models import ProcessControler

from utils.utils import Utils
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

    # A command must define handle()

    def set_submit_televir_sample(self, user: User, project_pk: int):
        """
        submit the job to televir
        """
        utils = Utils()
        process_SGE = ProcessSGE()
        user_pk = user.pk
        process_controler = ProcessControler()
        out_dir = utils.get_temp_dir()

        vect_command = [
            "python3 {} submit_televir_job --user_id {} --project_id {} -o {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                user_pk,
                project_pk,
                out_dir,
            )
        ]

        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        queue_name = user.profile.queue_name_sge
        (job_name_wait, job_name) = user.profile.get_name_sge_seq(
            Profile.SGE_PROCESS_dont_care, Profile.SGE_LINK
        )
        outdir_sge = utils.get_temp_dir()
        path_file = process_SGE.set_script_run_sge(
            outdir_sge,
            queue_name,
            vect_command,
            job_name,
            True,
            [job_name_wait],
            alternative_temp_dir=out_dir,
        )
        try:
            sge_id = process_SGE.submitte_job(path_file)
            print("project submitted, sge_id: " + str(sge_id))
            if sge_id != None:
                process_SGE.set_process_controlers(
                    user, process_controler.get_name_televir_project(project_pk), sge_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return sge_id

    def handle(self, *args, **options):

        project_name = options["project_name"]
        sample_name = options["sample_name"]
        account = options["user_login"]

        try:

            process_SGE = ProcessSGE()
            user = User.objects.get(username=account)
            sample = Sample.objects.get(name=sample_name, owner=user)

            project_count = Projects.objects.filter(
                name__iexact=project_name,
                is_deleted=False,
                owner__username=user.username,
            ).count()

            if project_count == 0:
                with transaction.atomic():
                    project = Projects()
                    project.name = project_name
                    project.owner = user
                    project.owner_id = user.id
                    # TODO Check where these constants are, or define them somewhere...
                    technology = ConstantsSettings.TECHNOLOGY_minion
                    if sample.type_of_fastq == Sample.TYPE_OF_FASTQ_illumina:
                        technology = ConstantsSettings.TECHNOLOGY_illumina
                    project.technology = technology
                    project.save()

            if project_count > 0:
                project = Projects.objects.filter(
                    name__iexact=project_name,
                    is_deleted=False,
                    owner__username=user.username,
                )[0]

                with transaction.atomic():
                    project_sample = PIProject_Sample()
                    project_sample.project = project
                    project_sample.sample = sample
                    project_sample.name = sample.name
                    project_sample_input = sample.file_name_1
                    if sample.is_valid_2:
                        project_sample_input += ";" + sample.file_name_2
                    project_sample.input = project_sample_input
                    sample_technology = "ONT"
                    if sample.type_of_fastq == Sample.TYPE_OF_FASTQ_illumina:
                        sample_technology = "Illumina/IonTorrent"
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
                runs_to_deploy = utils.check_runs_to_deploy(user, project)

                if runs_to_deploy:
                    taskID = process_SGE.set_submit_televir_job(
                        user=user,
                        project_pk=project.pk,
                    )
                    self.stdout.write(
                        "Project submitted as task {}.".format(project.id, taskID)
                    )

            else:

                sample = Sample.objects.get(name=sample_name, owner=user)

                with transaction.atomic():
                    project = Projects()
                    project.name = project_name
                    project.owner = user
                    project.owner_id = user.id
                    # TODO Check where these constants are, or define them somewhere...
                    technology = ConstantsSettings.TECHNOLOGY_minion
                    if sample.type_of_fastq == Sample.TYPE_OF_FASTQ_illumina:
                        technology = ConstantsSettings.TECHNOLOGY_illumina
                    project.technology = technology
                    project.save()
                    project_sample_input = sample.file_name_1
                    if sample.is_valid_2:
                        project_sample_input += ";" + sample.file_name_2
                    project_sample = PIProject_Sample()
                    project_sample.project = project
                    project_sample.sample = sample
                    project_sample.name = sample.name
                    project_sample.input = project_sample_input
                    project_sample.technology = sample.type_of_fastq
                    project_sample.report = "report"
                    project_sample.save()
                    self.stdout.write("Project created with id {}.".format(project.id))

                utils = Utils_Manager()
                runs_to_deploy = utils.check_runs_to_deploy(user, project)

                if runs_to_deploy:
                    taskID = process_SGE.set_submit_televir_job(
                        user=user,
                        project_pk=project.pk,
                    )
                    self.stdout.write(
                        "Project submitted as task {}.".format(project.id, taskID)
                    )

        except User.DoesNotExist as e:
            self.stdout.write("Error: User '{}' does not exist.".format(account))
        except Exception as e:
            self.stdout.write("Error: {}.".format(e))
