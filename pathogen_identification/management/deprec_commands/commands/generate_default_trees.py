from django.conf import settings
from django.contrib.auth.models import User
from django.core.management import BaseCommand
from django.core.management.base import BaseCommand

from constants.constants import Constants
from extend_user.models import Profile
from managing_files.models import DataSet
from settings.default_software import DefaultSoftware


class Command(BaseCommand):
    help = "deploy run"

    def create_default_user(self):
        """
        create a default user to link the default references...
        """
        self.__create_account(
            Constants.DEFAULT_USER,
            Constants.DEFAULT_USER_PASS,
            settings.DEFAULT_USER_EMAIL,
            False,
        )
        self.__create_account(
            Constants.USER_ANONYMOUS,
            Constants.USER_ANONYMOUS_PASS,
            settings.USER_ANONYMOUS_EMAIL,
            True,
        )
    
    def __create_account(self, user_name, password, email, b_active):
        """
        try to create a default accounts
        """
        try:
            user = User.objects.get(username=user_name)
            ### great, the default user exist
        except User.DoesNotExist:

            ### need to create it
            self.stdout.write("Add user: {}".format(user_name))
            user = User()
            user.username = user_name
            user.set_password(password)
            user.first_name = user_name
            user.email = email
            user.is_active = b_active
            user.is_staff = False
            user.is_superuser = False
            user.save()

        ### create generic dataset
        for user in User.objects.all():
            result = DataSet.objects.filter(owner__id=user.id)
            if len(result) == 0:
                ### need to create it
                dataSet = DataSet()
                dataSet.name = Constants.DATA_SET_GENERIC
                dataSet.owner = user
                dataSet.save()

        ### for security reasons
        ### set true for anonymous user always
        if user_name == Constants.USER_ANONYMOUS:
            try:
                profile = Profile.objects.get(user__username=user_name)
                profile.only_view_project = True
                profile.save()
            except Profile.DoesNotExist:
                pass


    def handle(self, *args, **options):
        ###
        from settings.constants_settings import ConstantsSettings
        default_software = DefaultSoftware()
        self.stdout.write("Set default users...")
        self.create_default_user()
        from pathogen_identification.utilities.utilities_pipeline import \
            SoftwareTreeUtils
        for user in User.objects.all():
            default_software.test_all_defaults(user)
            for technology in [ConstantsSettings.TECHNOLOGY_illumina, ConstantsSettings.TECHNOLOGY_minion]:
                software_utils = SoftwareTreeUtils(user, None, None)
                software_utils.set_technology(technology)
                software_utils.deactivate_all_nodes()
                _ = software_utils.get_sample_pathnodes(
                    screening=False,
                    mapping_only=True,
                )
                _ = software_utils.get_sample_pathnodes(
                    screening=True,
                    mapping_only=False,
                )
                _ = software_utils.get_sample_pathnodes(
                    screening=False,
                    mapping_only=False,
                )