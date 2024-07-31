import os
from datetime import date
from typing import List

from django.contrib import messages
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from constants.meta_key_and_values import MetaKeyAndValue
from extend_user.models import Profile
from managing_files.manage_database import ManageDatabase
from managing_files.models import ProcessControler
from managing_files.models import Project as InsafluProject
from managing_files.models import ProjectSample
from pathogen_identification.models import TeleFluProject, TeleFluSample
from pathogen_identification.utilities.reference_utils import (
    teleflu_to_insaflu_reference,
)
from settings.default_software_project_sample import DefaultProjectSoftware
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
            "--project_id",
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
        project_id = int(options["project_id"])

        # PROCESS CONTROLER
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_teleflu_project_create(
                project_id=project_id,
            ),
            ProcessControler.FLAG_RUNNING,
        )

        process = ProcessControler.objects.filter(
            owner__id=user.pk,
            name=process_controler.get_name_televir_teleflu_project_create(
                project_id=project_id,
            ),
            is_finished=False,
        )

        if process.exists():
            process = process.first()

            print(f"Process {process.pk} has been submitted")

        else:
            print("Process does not exist")

        # UTILITIES
        ### test anonymous account
        try:
            profile = Profile.objects.get(user=user)

        except Profile.DoesNotExist:
            raise Exception("User without profile")

        try:
            success, insaflu_reference = teleflu_to_insaflu_reference(
                project_id, user_id
            )

            print("Success", success, insaflu_reference)

            if success is False:
                process_SGE.set_process_controler(
                    user,
                    process_controler.get_name_televir_teleflu_project_create(
                        project_id=project_id,
                    ),
                    ProcessControler.FLAG_ERROR,
                )
                return

            teleflu_project = TeleFluProject.objects.get(pk=project_id)

            insaflu_project = InsafluProject.objects.create(
                owner=user,
                name=teleflu_project.name,
                reference=insaflu_reference,
            )

            teleflu_project.insaflu_project = insaflu_project
            teleflu_project.save()

            default_software = DefaultProjectSoftware()
            default_software.test_all_defaults(
                user, insaflu_project, None, None
            )  ## the user can have defaults yet

            process_SGE = ProcessSGE()
            process_SGE.set_create_project_list_by_user(user)

            ###### SAMPLES
            samples = TeleFluSample.objects.filter(teleflu_project=teleflu_project)
            for sample in samples:
                original_sample = sample.televir_sample.sample

                try:
                    project_sample = ProjectSample.objects.get(
                        project=insaflu_project,
                        sample=original_sample,
                    )
                except ProjectSample.DoesNotExist:
                    project_sample = ProjectSample.objects.create(
                        project=insaflu_project,
                        sample=original_sample,
                    )

                (job_name_wait, job_name) = ("", "")
                ### create a task to perform the analysis of snippy and freebayes
                ### Important, it is necessary to run again because can have some changes in the parameters.
                manageDatabase = ManageDatabase()
                metaKeyAndValue = MetaKeyAndValue()
                try:
                    if len(job_name_wait) == 0:
                        (
                            job_name_wait,
                            job_name,
                        ) = profile.get_name_sge_seq(
                            Profile.SGE_PROCESS_projects, Profile.SGE_GLOBAL
                        )
                    if original_sample.is_type_fastq_gz_sequencing():
                        taskID = process_SGE.set_second_stage_snippy(
                            project_sample, user, job_name, [job_name_wait]
                        )
                    else:
                        taskID = process_SGE.set_second_stage_medaka(
                            project_sample, user, job_name, [job_name_wait]
                        )

                    ### set project sample queue ID
                    manageDatabase.set_project_sample_metakey(
                        project_sample,
                        user,
                        metaKeyAndValue.get_meta_key_queue_by_project_sample_id(
                            project_sample.id
                        ),
                        MetaKeyAndValue.META_VALUE_Queue,
                        taskID,
                    )

                    ### need to collect global files again
                    taskID = process_SGE.set_collect_global_files(insaflu_project, user)
                    manageDatabase.set_project_metakey(
                        insaflu_project,
                        user,
                        metaKeyAndValue.get_meta_key(
                            MetaKeyAndValue.META_KEY_Queue_TaskID_Project,
                            insaflu_project.id,
                        ),
                        MetaKeyAndValue.META_VALUE_Queue,
                        taskID,
                    )

                except Exception as e:
                    print(e)
                    pass

        except Exception as e:
            print(e)
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_teleflu_project_create(
                    project_id=project_id,
                ),
                ProcessControler.FLAG_ERROR,
            )
            raise e

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_teleflu_project_create(
                project_id=project_id,
            ),
            ProcessControler.FLAG_FINISHED,
        )
