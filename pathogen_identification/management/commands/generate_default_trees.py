import os
from datetime import date

from django.core.management.base import BaseCommand
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager


class Command(BaseCommand):
    help = "deploy run"

    def handle(self, *args, **options):
        ###
        utils = Utils_Manager()
        utils.generate_default_trees()
