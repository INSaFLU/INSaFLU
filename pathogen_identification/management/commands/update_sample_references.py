from datetime import date

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from pathogen_identification.models import PIProject_Sample
from pathogen_identification.utilities.utilities_views import RawReferenceUtils


class Command(BaseCommand):
    help = "deploy run"

    def handle(self, *args, **options):

        all_project_samples = PIProject_Sample.objects.filter(
            is_deleted=False, project__is_deleted=False
        )

        print(f"Updating references for {all_project_samples.count()} samples")

        for sample in all_project_samples:
            reference_utils = RawReferenceUtils(sample)

            _ = reference_utils.retrieve_compound_references()
