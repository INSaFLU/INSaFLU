"""
Created on Jan 5, 2018

@author: mmp
"""

import logging

from django.contrib.auth.models import User
from django.core.management import BaseCommand

from constants.meta_key_and_values import MetaKeyAndValue
from managing_files.manage_database import ManageDatabase
from managing_files.models import ProjectSample
from utils.process_SGE import ProcessSched
from utils.software_minion import SoftwareMinion


class Command(BaseCommand):
    """
    classdocs
    Ex: python3 manage.py second_stage_medaka --project_sample_id 19 --user_id 1

    """

    help = "Run second stage medaka."

    ## logging
    logger_debug = logging.getLogger("fluWebVirus.debug")
    logger_production = logging.getLogger("fluWebVirus.production")

    def __init__(self, *args, **kwargs):
        super(Command, self).__init__(*args, **kwargs)

    ## https://docs.djangoproject.com/en/dev/howto/custom-management-commands/
    def add_arguments(self, parser):
        parser.add_argument(
            "--project_sample_id", type=int, help="Project sample id to process"
        )
        parser.add_argument(
            "--user_id", nargs="?", type=int, help="User id to join to this process"
        )

    # A command must define handle()
    def handle(self, *args, **options):

        software_minion = SoftwareMinion()
        process_SGE = ProcessSched()
        manageDatabase = ManageDatabase()
        metaKeyAndValue = MetaKeyAndValue()

        project_sample_id = options["project_sample_id"]
        user_id = options["user_id"]
        self.stdout.write("Starting for project_sample_id: " + str(project_sample_id))
        self.logger_production.info(
            "Starting for project_sample_id: " + str(project_sample_id)
        )
        self.logger_debug.info(
            "Starting for project_sample_id: " + str(project_sample_id)
        )
        try:
            project_sample = ProjectSample.objects.get(pk=project_sample_id)
            if user_id is None:
                user = project_sample.project.owner
            else:
                user = User.objects.get(pk=user_id)
            software_minion.process_second_stage_medaka(project_sample, user)

            ### WAIT FOR TESTS
            # is_teleflu_project = TeleFluProject.objects.filter(
            #    insaflu_project=project_sample.project,
            # ).exists()
            # if is_teleflu_project:
            #    teleflu_project = TeleFluProject.objects.get(
            #        insaflu_project=project_sample.project
            #    )
            #    create_teleflu_igv_report(teleflu_project.pk)
                                ### need to collect global files again
            

            taskID = process_SGE.set_collect_global_files(
                project_sample.project, project_sample.project.owner
            )
            manageDatabase.set_project_metakey(
                project_sample.project,
                project_sample.project.owner,
                metaKeyAndValue.get_meta_key(
                    MetaKeyAndValue.META_KEY_Queue_TaskID_Project,
                    project_sample.project.id,
                ),
                MetaKeyAndValue.META_VALUE_Queue,
                taskID,
            )
            self.stdout.write("End")
        except ProjectSample.DoesNotExist as e:
            self.stdout.write(
                "Error: ProjectSample id '{}' does not exist.".format(project_sample_id)
            )
        except User.DoesNotExist as e:
            self.stdout.write("Error: User id '{}' does not exist.".format(user_id))

        except Exception as e:
            self.stdout.write("Error: {}".format(str(e)))
