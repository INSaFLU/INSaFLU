import os
from datetime import date

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from pathogen_identification.utilities.utilities_pipeline import Utils_Manager
from settings.default_software import DefaultSoftware
from settings.models import Software


class Command(BaseCommand):
    help = "deploy run"

    def handle(self, *args, **options):
        ###
        user_system = User.objects.get(username="system")
        default_software = DefaultSoftware()

        for user in User.objects.all():
            default_software.test_all_defaults(user)

        default_software.remove_all_parameters(user_system)
