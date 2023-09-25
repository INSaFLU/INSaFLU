"""
Created on January 30, 2023

@author: daniel.sobral
"""
import logging
import os

from django.contrib.auth.models import User
from django.core.management import BaseCommand

from managing_files.models import Sample


def get_sample_name(user: User, sample_name: str) -> Sample:
    sample = (
        Sample.objects.filter(name=sample_name, owner=user, is_deleted=False)
        .order_by("creation_date")
        .last()
    )

    if sample is None:
        return None

    return sample


class Command(BaseCommand):
    """
    classdocs
    """

    help = "Check status of sample (whether it is ready to be used)."

    # logging
    logger_debug = logging.getLogger("fluWebVirus.debug")
    logger_production = logging.getLogger("fluWebVirus.production")

    def __init__(self, *args, **kwargs):
        super(Command, self).__init__(*args, **kwargs)

    def add_arguments(self, parser):
        parser.add_argument(
            "--name", nargs="?", type=str, required=True, help="Sample Name"
        )
        parser.add_argument(
            "--user_login",
            nargs="?",
            type=str,
            required=True,
            help="User login of the sample owner",
        )
        parser.add_argument(
            "--test_ready",
            action="store_true",
            required=False,
            default=False,
            help="Test if sample is ready for projects message",
        )
        parser.add_argument(
            "--test_missing",
            action="store_true",
            required=False,
            default=False,
            help="Test if sample is ready for projects message",
        )

    # A command must define handle()
    def success_message(self, sample_name: str, is_ready: bool):
        self.stdout.write(f"Sample {sample_name}. Is Ready: {is_ready}")

    def sample_does_not_exist_message(self, sample_name: str):
        self.stdout.write(f"Sample {sample_name} does not exist.")

    def handle(self, *args, **options):
        sample_name = options["name"]
        account = options["user_login"]

        if options["test_ready"]:
            self.success_message(sample_name, True)
            return False

        if options["test_missing"]:
            self.sample_does_not_exist_message(sample_name)
            return False

        try:
            user = User.objects.get(username=account)
            sample = get_sample_name(user, sample_name)

            if sample:
                self.success_message(sample_name, sample.is_ready_for_projects)

            else:
                self.sample_does_not_exist_message(sample_name)

        except User.DoesNotExist as e:
            self.stdout.write("Error: User '{}' does not exist.".format(account))
        except Exception as e:
            self.stdout.write("Error: {}.".format(e))
