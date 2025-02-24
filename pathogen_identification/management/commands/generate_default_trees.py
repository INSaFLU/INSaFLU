from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from settings.default_software import DefaultSoftware


class Command(BaseCommand):
    help = "deploy run"

    def handle(self, *args, **options):
        ###
        user_system = User.objects.get(username="system")
        default_software = DefaultSoftware()

        for user in User.objects.all():
            default_software.test_all_defaults(user)

        # default_software.remove_all_parameters(user_system)
