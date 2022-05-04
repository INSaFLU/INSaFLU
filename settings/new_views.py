from braces.views import LoginRequiredMixin
from constants.meta_key_and_values import MetaKeyAndValue
from django.contrib import messages
from django.db import transaction
from django.urls import reverse_lazy
from django.views.generic import ListView, TemplateView, UpdateView
from extend_user.models import Profile
from managing_files.manage_database import ManageDatabase
from managing_files.models import Project, ProjectSample, Sample
from utils.process_SGE import ProcessSGE
from utils.utils import ShowInfoMainPage

from settings.constants_settings import ConstantsSettings
from settings.default_software import DefaultSoftware
from settings.forms import SoftwareForm
from settings.models import Parameter, Software
from settings.tables import SoftwaresTable
