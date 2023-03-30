"""
Created on January 30, 2023

@author: daniel.sobral
"""
from django.core.management import BaseCommand
from django.contrib.auth.models import User
from managing_files.models import Sample
import os
import logging


def get_sample_name(user: User, sample_name: str) -> Sample:
    sample = (
        Sample.objects.filter(name=sample_name, owner=user)
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

    # A command must define handle()
    def handle(self, *args, **options):

        sample_name = options["name"]
        account = options["user_login"]

        try:

            user = User.objects.get(username=account)
            sample = get_sample_name(user, sample_name)

            self.stdout.write("Is Ready: {}".format(sample.is_ready_for_projects))

        except User.DoesNotExist as e:
            self.stdout.write("Error: User '{}' does not exist.".format(account))
        except Exception as e:
            self.stdout.write("Error: {}.".format(e))
