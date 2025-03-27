"""
Created on Jan 5, 2018

@author: mmp
"""

import logging

from django.core.management import BaseCommand

from utils.software import SoftwareFlumut


class Command(BaseCommand):
    """
    classdocs
    Ex: python3 manage.py update_flumut

    """

    help = "Update flumut learn."

    ## logging
    logger = logging.getLogger("fluWebVirus.update_flumut")

    def __init__(self, *args, **kwargs):
        super(Command, self).__init__(*args, **kwargs)

    # A command must define handle()
    def handle(self, *args, **options):

        software_flumut = SoftwareFlumut()
        self.stdout.write("Starting update flumut")
        self.logger.info("Starting update flumut")

        try:
            software_flumut.run_flumut_update()
            self.logger.info("End update flumut")
            self.stdout.write("Success update flumut")
        except Exception as e:
            import traceback

            traceback.print_exc()
            message = "Fail to update flumut: " + str(e)
            self.logger.info(message)
            self.stdout.write(message)
