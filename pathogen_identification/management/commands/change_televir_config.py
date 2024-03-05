import logging
import os

import pandas as pd
from django.contrib.auth.models import User
from django.core.management import BaseCommand

from constants.software_names import SoftwareNames
from fluwebvirus.settings import MEDIA_ROOT
from settings.constants_settings import ConstantsSettings
from settings.models import Parameter, Software


def remapping_can_be_off():
    remap_software = Software.objects.filter(
        pipeline_step__name=ConstantsSettings.PIPELINE_NAME_remapping,
    )
    for software in remap_software:
        software.can_be_on_off_in_pipeline = True
        software.save()


def update_remap_params():
    from settings.default_software import DefaultSoftware

    default_software = DefaultSoftware()

    software_remap = Software.objects.filter(
        name=SoftwareNames.SOFTWARE_REMAP_PARAMS_name
    )

    for software in software_remap:
        user = software.owner
        if user is None:
            continue

        remap_params = default_software.default_parameters.get_remap_defaults(
            user,
            Software.TYPE_OF_USE_televir_settings,
            ConstantsSettings.TECHNOLOGY_illumina,
        )

        default_software.default_parameters.persist_parameters_update(
            vect_parameters=remap_params,
            software=software,
        )


class Command(BaseCommand):
    """
    classdocs
    """

    help = "Checks if a given televir project finished."

    # logging
    logger_debug = logging.getLogger("fluWebVirus.debug")
    logger_production = logging.getLogger("fluWebVirus.production")

    def __init__(self, *args, **kwargs):
        super(Command, self).__init__(*args, **kwargs)

    def add_arguments(self, parser):
        parser.add_argument(
            "--metagenomics",
            action="store_true",
            required=False,
            default=False,
            help="Metagenomics",
        )

    # A command must define handle()
    def handle(self, *args, **options):
        if options["metagenomics"]:
            remapping_can_be_off()
            return
