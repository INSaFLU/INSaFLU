"""
Created on Jan 5, 2018

@author: mmp
"""

import logging

from django.contrib.auth.models import User
from django.core.management import BaseCommand

from managing_files.models import Project
from utils.collect_extra_data import CollectExtraData


class Command(BaseCommand):
    """
    classdocs
    Ex: python3 manage.py collect_update_mutation_report --project_id 11 --user_id 1

    """

    help = "Create global files by project sample and update pangolin."

    ## logging
    logger_debug = logging.getLogger("fluWebVirus.debug")
    logger_production = logging.getLogger("fluWebVirus.production")

    def __init__(self, *args, **kwargs):
        super(Command, self).__init__(*args, **kwargs)

    ## https://docs.djangoproject.com/en/dev/howto/custom-management-commands/
    def add_arguments(self, parser):
        parser.add_argument("--project_id", type=int, help="Project id to process")
        parser.add_argument(
            "--user_id", nargs="?", type=int, help="User id to join to this process"
        )

    # A command must define handle()
    def handle(self, *args, **options):

        collect_extra_data = CollectExtraData()
        project_id = options["project_id"]
        user_id = options["user_id"]
        self.stdout.write("Starting for project_id: " + str(project_id))
        self.logger_production.info(
            "Starting update mutation report for project_id: " + str(project_id)
        )
        self.logger_debug.info(
            "Starting update mutation report for project_id: " + str(project_id)
        )

        try:
            project = Project.objects.get(pk=project_id)
            if user_id == None:
                user = project.owner
            else:
                user = User.objects.get(pk=user_id)
            collect_extra_data.__collect_mutation_report(project, user)
            self.stdout.write("End")
        except Project.DoesNotExist as e:
            self.stdout.write("Project id '{}' does not exist.".format(project_id))
        except User.DoesNotExist as e:
            self.stdout.write("Error: User id '{}' does not exist.".format(user_id))
