import os
from datetime import date

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager
from settings.default_software import DefaultSoftware


class Command(BaseCommand):
    help = "deploy run"

    def handle(self, *args, **options):
        ###
        user_system = User.objects.get(username="system")
        default_software = DefaultSoftware()
        default_software.test_all_defaults_pathogen_identification(user_system)

        utils = Utils_Manager()
        utils.generate_default_trees(user_system)

        default_software.remove_all_parameters(user_system)
