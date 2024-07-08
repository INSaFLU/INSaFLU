import logging
import mimetypes
import ntpath
import os
from typing import Any, Dict, Optional

import pandas as pd
from Bio import SeqIO
from braces.views import FormValidMessageMixin, LoginRequiredMixin
from django import forms
from django.conf import settings
from django.contrib import messages
from django.core.files.temp import NamedTemporaryFile
from django.db import transaction
from django.db.models import Q
from django.http import (Http404, HttpResponse, HttpResponseNotFound,
                         HttpResponseRedirect)
from django.http.response import HttpResponse
from django.shortcuts import render
from django.template.defaultfilters import pluralize
from django.urls import reverse, reverse_lazy
from django.utils.safestring import mark_safe
from django.views import generic
from django.views.generic import ListView
from django_tables2 import RequestConfig

from constants.constants import Constants
from extend_user.models import Profile
from fluwebvirus.settings import MEDIA_ROOT, STATICFILES_DIRS
from managing_files.forms import AddSampleProjectForm
from managing_files.models import ProcessControler
from managing_files.models import ProjectSample as InsafluProjectSample
from managing_files.models import Reference
from managing_files.tables import SampleToProjectsTable
from pathogen_identification.constants_settings import ConstantsSettings
from pathogen_identification.constants_settings import \
    ConstantsSettings as PICS
from pathogen_identification.forms import (PanelReferencesUploadForm,
                                           ReferenceForm, UploadFileForm)
from pathogen_identification.models import (ContigClassification, FinalReport,
                                            ParameterSet, PIProject_Sample,
                                            Projects, RawReference,
                                            ReadClassification,
                                            ReferenceContigs,
                                            ReferenceMap_Main, ReferencePanel,
                                            ReferenceSourceFile,
                                            ReferenceSourceFileMap,
                                            RunAssembly, RunDetail, RunMain,
                                            RunRemapMain, Sample,
                                            TelefluMapping, TeleFluProject,
                                            TeleFluSample, TelevirRunQC)
from pathogen_identification.modules.object_classes import RunQC_report
from pathogen_identification.tables import (AddedReferenceTable,
                                            CompoundRefereceScoreWithScreening,
                                            CompoundReferenceScore,
                                            ContigTable, ProjectTable,
                                            RawReferenceTable,
                                            RawReferenceTable_Basic,
                                            ReferenceSourceTable, RunMainTable,
                                            RunMappingTable, SampleTableOne,
                                            TeleFluInsaFLuProjectTable,
                                            TeleFluReferenceTable)
from pathogen_identification.utilities.reference_utils import \
    generate_insaflu_reference
from pathogen_identification.utilities.televir_bioinf import TelevirBioinf
from pathogen_identification.utilities.televir_parameters import \
    TelevirParameters
from pathogen_identification.utilities.tree_deployment import TreeProgressGraph
from pathogen_identification.utilities.utilities_general import (
    get_services_dir, infer_run_media_dir)
from pathogen_identification.utilities.utilities_pipeline import (
    Parameter_DB_Utility, RawReferenceUtils, SoftwareTreeUtils)
from pathogen_identification.utilities.utilities_views import (
    EmptyRemapMain, ReportSorter, final_report_best_cov_by_accid,
    recover_assembly_contigs)
from settings.constants_settings import ConstantsSettings as CS
from utils.process_SGE import ProcessSGE
from utils.software import Software
from utils.support_django_template import get_link_for_dropdown_item
from utils.utils import ShowInfoMainPage, Utils


class UploadPanelReferencesView:
    pass


def process_query_string(query_string: Optional[str]):
    if query_string is None:
        return ""

    query_string = query_string.strip()

    return query_string


class UploadNewReferencesView(
    LoginRequiredMixin, FormValidMessageMixin, generic.FormView
):

    form_class = PanelReferencesUploadForm
    template_name = "datasets/datasets_new_consensus.html"

    ## Other solution to get the Consensus
    ## https://pypi.python.org/pypi?%3aaction=display&name=django-contrib-requestprovider&version=1.0.1
    def get_form_kwargs(self):
        """
        Set the request to pass in the form
        """
        kw = super(UploadNewReferencesView, self).get_form_kwargs()
        kw["request"] = self.request  # the trick!
        kw["pk"] = self.kwargs["pk"]
        return kw

    def get_success_url(self):
        """
        get source_pk from consensus dataset, need to pass it in context
        """
        return reverse_lazy(
            "add-consensus-dataset", kwargs={"pk": self.kwargs.get("pk")}
        )

    def get_context_data(self, **kwargs):
        context = super(UploadNewReferencesView, self).get_context_data(**kwargs)
        context["nav_dataset"] = True
        context["nav_modal"] = True  ## short the size of modal window
        context["pk"] = self.kwargs["pk"]  ## pk of dataset, need to return
        context["show_info_main_page"] = (
            ShowInfoMainPage()
        )  ## show main information about the institute
        return context

    @transaction.atomic
    def form_valid(self, form):
        software = Software()
        utils = Utils()
        panel_pk = int(self.kwargs["pk"])

        ### test anonymous account
        try:
            profile = Profile.objects.get(user=self.request.user)
            if profile.only_view_project:
                messages.warning(
                    self.request,
                    "'{}' account can not add references.".format(
                        self.request.user.username
                    ),
                    fail_silently=True,
                )
                return super(UploadNewReferencesView, self).form_invalid(form)
        except Profile.DoesNotExist:
            pass

        if form.is_valid():
            name = form.cleaned_data["name"]
            reference_fasta = form.cleaned_data["reference_fasta"]

            reference = form.save(commit=False)
            ## set other data
            reference.display_name = reference.name
            reference.owner = self.request.user
            reference.is_obsolete = False
            reference.number_of_locus = self.request.session[
                Constants.NUMBER_LOCUS_FASTA_FILE
            ]
            reference.reference_fasta_name = utils.clean_name(
                ntpath.basename(reference_fasta.name)
            )
            reference.save()

            ## has the names thar need to pass, if empty pass all
            vect_names_to_upload = (
                self.request.session[Constants.SEQUENCES_TO_PASS].split(",")
                if len(self.request.session[Constants.SEQUENCES_TO_PASS]) > 0
                else []
            )

            ## get panel
            panel = ReferencePanel.objects.get(pk=panel_pk)

            ## clean file
            original_file_name = os.path.join(
                settings.MEDIA_ROOT, reference.reference_fasta.name
            )
            software.dos_2_unix(original_file_name)
            ## test if bases all upper
            software.fasta_2_upper(original_file_name)

            ### check if there all seq names are present in the database yet, and create new files
            ### At least has one sequence
            vect_fail, vect_pass = [], []
            with open(
                os.path.join(settings.MEDIA_ROOT, reference.reference_fasta.name)
            ) as handle_in:
                for record in SeqIO.parse(handle_in, "fasta"):

                    ## name
                    seq_name = (
                        "{}_{}".format(name, record.id) if len(name) > 0 else record.id
                    )

                    ## check the select ones
                    if (
                        len(vect_names_to_upload)
                    ) > 0 and not record.id in vect_names_to_upload:
                        vect_fail.append(seq_name)
                        continue

                    ## try to upload
                    reference_fasta_temp_file_name = NamedTemporaryFile(
                        prefix="flu_fa_", suffix=".fasta", delete=False
                    )

                    with open(reference_fasta_temp_file_name.name, "w") as handle_out:
                        SeqIO.write(record, handle_out, "fasta")

                    ## create the reference
                    final_fasta_filename = f"{seq_name}.fasta"
                    success_create, ref_pk = generate_insaflu_reference(
                        reference_fasta_temp_file_name.name,
                        seq_name,
                        final_fasta_filename,
                        self.request.user,
                    )

                    if success_create == False or ref_pk == None:
                        vect_fail.append(seq_name)
                        continue

                    vect_pass.append(seq_name)
                    insaflu_reference = Reference.objects.get(pk=ref_pk)

            utils.remove_file(original_file_name)
            message = (
                "Consensus '" + "', '".join(vect_pass) + "' were created successfully."
            )
            if len(vect_fail) > 0:
                message += " Not uploaded consensus '" + "', '".join(vect_fail) + "'."
            messages.success(self.request, message, fail_silently=True)
            return super(UploadNewReferencesView, self).form_valid(form)
        else:
            return super(UploadNewReferencesView, self).form_invalid(form)

    form_valid_message = ""  ## need to have this


def clean_check_box_in_session(request):
    """
    check all check boxes on samples/references to add samples
    """
    utils = Utils()
    ## clean all check unique
    if Constants.CHECK_BOX_ALL in request.session:
        del request.session[Constants.CHECK_BOX_ALL]
    vect_keys_to_remove = []
    for key in request.session.keys():
        if (
            key.startswith(Constants.CHECK_BOX)
            and len(key.split("_")) == 3
            and utils.is_integer(key.split("_")[2])
        ):
            vect_keys_to_remove.append(key)
    for key in vect_keys_to_remove:
        del request.session[key]


def is_all_check_box_in_session(
    vect_check_to_test, request, prefix: str = Constants.CHECK_BOX
):
    """
    test if all check boxes are in session
    If not remove the ones that are in and create the new ones all False
    """
    utils = Utils()
    dt_data = {}

    ## get the dictonary
    for key in request.session.keys():
        if (
            key.startswith(prefix)
            and len(key.split("_")) == 3
            and utils.is_integer(key.split("_")[2])
        ):
            dt_data[key] = True

    b_different = False
    if len(vect_check_to_test) != len(dt_data):
        b_different = True

    ## test the vector
    if not b_different:
        for key in vect_check_to_test:
            if key not in dt_data:
                b_different = True
                break

    if b_different:
        ## remove all
        for key in dt_data:
            del request.session[key]

        ## create new
        for key in vect_check_to_test:
            request.session[key] = False
        return False
    return True


# Create your views here.


class IGVform(forms.Form):
    sample_pk = forms.CharField(max_length=100)
    run_pk = forms.CharField(max_length=100)
    reference = forms.CharField(max_length=100)
    unique_id = forms.CharField(max_length=100)
    project_pk = forms.CharField(max_length=100)


class download_form(forms.Form):
    file_path = forms.CharField(max_length=300)

    class Meta:
        widgets = {
            "myfield": forms.TextInput(
                attrs={"style": "border-color:darkgoldenrod; border-radius: 10px;"}
            ),
        }


class download_ref_form(forms.Form):
    file = forms.CharField(max_length=300)
    taxid = forms.CharField(max_length=50)
    run = forms.IntegerField()
    accid = forms.CharField(max_length=50)

    class Meta:
        widgets = {
            "myfield": forms.TextInput(
                attrs={"style": "border-color:darkgoldenrod; border-radius: 10px;"}
            ),
        }


################################################################


class Services(LoginRequiredMixin, generic.CreateView):
    """
    Display a series of applications to run, independent of projects.
    """

    template_name = "pathogen_identification/services.html"

    def get_context_data(self, **kwargs):
        context = {}

        context["nav_services"] = True

        user = self.request.user
        services_dir = get_services_dir(user)
        context["services_dir"] = services_dir
        explify_output_file = os.path.join(services_dir, "merged_televir_explify.tsv")
        explify_file_exists = os.path.exists(explify_output_file)
        merge_explify_file_provide = os.path.join(
            "/media/",
            "televir_projects",
            str(user.pk),
            "services",
            "merged_televir_explify" + ".tsv",
        )

        ## check if merging is runnning
        process_controler = ProcessControler()
        merger_running = ProcessControler.objects.filter(
            name=process_controler.get_name_televir_project_merge_explify_external(
                user_pk=user.pk,
            ),
            is_running=True,
        ).exists()

        context["explify_file_exists"] = explify_file_exists
        context["explify_output_file"] = merge_explify_file_provide
        context["merger_running"] = merger_running
        context["tools"] = []

        return context


class PathId_ProjectsView(LoginRequiredMixin, ListView):
    model = Projects
    template_name = "pathogen_identification/projects.html"
    context_object_name = "projects"
    ##	group_required = u'company-user' security related with GroupRequiredMixin

    def get_context_data(self, **kwargs):
        context = super(PathId_ProjectsView, self).get_context_data(**kwargs)
        tag_search = "search_projects"
        query_set = Projects.objects.filter(
            owner__id=self.request.user.id, is_deleted=False
        ).order_by("-creation_date")

        if self.request.GET.get(tag_search) != None and self.request.GET.get(
            tag_search
        ):
            query_string = process_query_string(self.request.GET.get(tag_search))
            query_set = query_set.filter(
                Q(name__icontains=self.request.GET.get(tag_search))
                | Q(project_samples__name__icontains=query_string)
            ).distinct()

        table = ProjectTable(query_set)

        RequestConfig(
            self.request, paginate={"per_page": Constants.PAGINATE_NUMBER}
        ).configure(table)
        if self.request.GET.get(tag_search) != None:
            context[tag_search] = self.request.GET.get(tag_search)

        ### clean check box in the session
        clean_check_box_in_session(self.request)  ## for both samples and references

        ### clean project name session
        if Constants.PROJECT_NAME_SESSION in self.request.session:
            del self.request.session[Constants.PROJECT_NAME_SESSION]

        context["table"] = table

        context["nav_project"] = True
        context["show_paginatior"] = query_set.count() > Constants.PAGINATE_NUMBER
        context["query_set_count"] = query_set.count()
        context["show_info_main_page"] = ShowInfoMainPage()
        context["query_set_count"] = query_set.count()
        return context


class PathID_ProjectCreateView(LoginRequiredMixin, generic.CreateView):
    """
    Create a new Project
    """

    # utils = Utils()
    model = Projects
    fields = ["name", "description", "technology"]
    success_url = reverse_lazy("PIprojects_main")
    template_name = "pathogen_identification/project_add.html"

    def get_form_kwargs(self):
        """
        Set the request to pass in the form
        """
        kw = super(PathID_ProjectCreateView, self).get_form_kwargs()
        return kw

    def get_context_data(self, **kwargs):
        context = super(PathID_ProjectCreateView, self).get_context_data(**kwargs)

        ### get project name
        project_name = (
            self.request.session[Constants.PROJECT_NAME_SESSION]
            if Constants.PROJECT_NAME_SESSION in self.request.session
            else ""
        )

        ## clean the check box in the session
        if Constants.PROJECT_NAME in self.request.session:
            context[Constants.PROJECT_NAME] = self.request.session[
                Constants.PROJECT_NAME
            ]
            del self.request.session[Constants.PROJECT_NAME]
        if Constants.ERROR_PROJECT_NAME in self.request.session:
            context[Constants.ERROR_PROJECT_NAME] = self.request.session[
                Constants.ERROR_PROJECT_NAME
            ]
            del self.request.session[Constants.ERROR_PROJECT_NAME]

        else:
            context[Constants.ERROR_PROJECT_NAME] = ""

        context["project_name"] = project_name
        context["show_paginatior"] = False
        context["nav_project"] = True
        context["show_info_main_page"] = (
            ShowInfoMainPage()
        )  ## show main information about the institute
        return context

    def form_valid(self, form):
        """
        Validate the form
        """
        ### test anonymous account
        try:
            profile = Profile.objects.get(user=self.request.user)
            if profile.only_view_project:
                messages.warning(
                    self.request,
                    "'{}' account can not create projects.".format(
                        self.request.user.username
                    ),
                    fail_silently=True,
                )
                return super(PathID_ProjectCreateView, self).form_invalid(form)
        except Profile.DoesNotExist:
            pass

        context = self.get_context_data()
        name = form.cleaned_data["name"]

        b_error = False
        try:
            Projects.objects.get(
                name__iexact=name,
                is_deleted=False,
                owner__username=self.request.user.username,
            )
            self.request.session[Constants.ERROR_PROJECT_NAME] = (
                "Exists a project with this name."
            )
            self.request.session[Constants.PROJECT_NAME] = name
            b_error = True
        except Projects.DoesNotExist:
            pass
        ###
        if context[Constants.ERROR_PROJECT_NAME] != "":
            b_error = True

        if not form.cleaned_data["name"]:
            self.request.session[Constants.ERROR_PROJECT_NAME] = (
                "The project name can not be empty."
            )
            self.request.session[Constants.PROJECT_NAME] = name
            b_error = True

        if not form.cleaned_data["name"].replace("_", "").isalnum():
            self.request.session[Constants.ERROR_PROJECT_NAME] = (
                "The project name can only contain letters and numbers."
            )
            self.request.session[Constants.PROJECT_NAME] = name
            b_error = True

        ### exists an error
        if b_error:
            return super(PathID_ProjectCreateView, self).form_invalid(form)

        with transaction.atomic():
            form.instance.owner = self.request.user

            project = form.save()
            project.owner = self.request.user
            project.owner_id = self.request.user.id
            project.technology = form.cleaned_data["technology"]
            project.save()

        return super(PathID_ProjectCreateView, self).form_valid(form)

    form_valid_message = ""  ## need to have this, even empty


class AddSamples_PIProjectsView(
    LoginRequiredMixin, FormValidMessageMixin, generic.CreateView
):
    """
    Create a new reference
    """

    utils = Utils()
    model = Sample
    fields = ["name"]
    success_url = reverse_lazy("PIprojects_main")
    template_name = "project_sample/project_sample_add.html"

    logger_debug = logging.getLogger("fluWebVirus.debug")
    logger_production = logging.getLogger("fluWebVirus.production")

    def get_context_data(self, **kwargs):
        context = super(AddSamples_PIProjectsView, self).get_context_data(**kwargs)

        ### test if the user is the same of the page
        project = Projects.objects.get(pk=self.kwargs["pk"])
        technology = None

        if project.technology == CS.TECHNOLOGY_minion:
            technology = Sample.TYPE_OF_FASTQ_minion
        elif project.technology == CS.TECHNOLOGY_illumina:
            technology = Sample.TYPE_OF_FASTQ_illumina

        context["nav_project"] = True
        if project.owner.id != self.request.user.id:
            context["error_cant_see"] = "1"
            return context

        ##
        count_active_projects = PIProject_Sample.objects.filter(
            project=project, is_deleted=False, is_error=False
        ).count()
        samples_out = PIProject_Sample.objects.filter(
            Q(project=project) & ~Q(is_deleted=True) & ~Q(is_error=True)
        ).values("sample__pk")
        query_set = Sample.objects.filter(
            owner__id=self.request.user.id,
            is_obsolete=False,
            is_deleted=False,
            is_deleted_processed_fastq=False,
            is_ready_for_projects=True,
            type_of_fastq=technology,
        ).exclude(pk__in=samples_out)

        tag_search = "search_add_project_sample"
        if self.request.GET.get(tag_search) != None and self.request.GET.get(
            tag_search
        ):
            query_string = process_query_string(self.request.GET.get(tag_search))
            query_set = query_set.filter(
                Q(name__icontains=query_string)
                | Q(type_subtype__icontains=query_string)
                | Q(data_set__name__icontains=query_string)
                | Q(week__icontains=query_string)
            )
            tag_search
        classified_refs_table = SampleToProjectsTable(query_set)

        ### set the check_box
        if Constants.CHECK_BOX_ALL not in self.request.session:
            self.request.session[Constants.CHECK_BOX_ALL] = False
            is_all_check_box_in_session(
                ["{}_{}".format(Constants.CHECK_BOX, key.id) for key in query_set],
                self.request,
            )

        context[Constants.CHECK_BOX_ALL] = self.request.session[Constants.CHECK_BOX_ALL]
        ## need to clean all the others if are reject in filter
        dt_sample_id_add_temp = {}
        if context[Constants.CHECK_BOX_ALL]:
            for sample in query_set:
                dt_sample_id_add_temp[sample.id] = (
                    1  ## add the ids that are in the tables
                )
            for key in self.request.session.keys():
                if (
                    key.startswith(Constants.CHECK_BOX)
                    and len(key.split("_")) == 3
                    and self.utils.is_integer(key.split("_")[2])
                ):
                    ### this is necessary because of the search. Can occur some checked box that are out of filter.
                    if int(key.split("_")[2]) not in dt_sample_id_add_temp:
                        self.request.session[key] = False
                    else:
                        self.request.session[key] = True
        ## END need to clean all the others if are reject in filter

        key_session_name_project_settings = "project_settings_{}".format(project.name)
        if not key_session_name_project_settings in self.request.session:
            self.request.session[key_session_name_project_settings] = True
        else:
            self.request.session[key_session_name_project_settings] = False

        RequestConfig(
            self.request, paginate={"per_page": Constants.PAGINATE_NUMBER}
        ).configure(classified_refs_table)
        if self.request.GET.get(tag_search) != None:
            context[tag_search] = self.request.GET.get(tag_search)

        context["televir_sample"] = True
        context["table"] = classified_refs_table
        context["show_paginatior"] = query_set.count() > Constants.PAGINATE_NUMBER
        context["query_set_count"] = query_set.count()
        context["project_name"] = project.name
        context["nav_modal"] = True  ## short the size of modal window
        context["show_info_main_page"] = (
            ShowInfoMainPage()
        )  ## show main information about the institute
        context["show_message_change_settings"] = (
            count_active_projects == 0
            and self.request.session[key_session_name_project_settings]
        )  ## Show message to change settings to the project
        context["add_all_samples_message"] = "Add {} sample{}".format(
            query_set.count(), pluralize(query_set.count(), ",s")
        )
        if self.request.POST:
            context["project_sample"] = AddSampleProjectForm(self.request.POST)
        else:
            context["project_sample"] = AddSampleProjectForm()
        return context

    def form_valid(self, form):
        """
        Validate the form
        """

        ### test anonymous account
        try:
            profile = Profile.objects.get(user=self.request.user)
            if profile.only_view_project:
                messages.warning(
                    self.request,
                    "'{}' account can not add samples to a project.".format(
                        self.request.user.username
                    ),
                    fail_silently=True,
                )
                return super(AddSamples_PIProjectsView, self).form_invalid(form)
        except Profile.DoesNotExist:
            pass

        if form.is_valid():

            ### get project sample..
            context = self.get_context_data()

            ### get project
            project = Projects.objects.get(pk=self.kwargs["pk"])

            vect_sample_id_add_temp = []
            for sample in context["table"].data:
                vect_sample_id_add_temp.append(sample.id)

            vect_sample_id_add = []
            if "submit_checked" in self.request.POST:
                for key in self.request.session.keys():
                    if (
                        self.request.session[key]
                        and key.startswith(Constants.CHECK_BOX)
                        and len(key.split("_")) == 3
                        and self.utils.is_integer(key.split("_")[2])
                    ):
                        ### this is necessary because of the search. Can occur some checked box that are out of filter.
                        if int(key.split("_")[2]) in vect_sample_id_add_temp:
                            vect_sample_id_add.append(int(key.split("_")[2]))
            elif "submit_all" in self.request.POST:
                vect_sample_id_add = vect_sample_id_add_temp

            ### get samples already out
            dt_sample_out = {}
            for project_sample in PIProject_Sample.objects.filter(
                project__id=project.id,
                is_deleted=False,
                is_error=False,
                is_deleted_in_file_system=False,
            ):
                dt_sample_out[project_sample.sample.id] = 1

            ### start adding...
            (job_name_wait, job_name) = ("", "")
            project_sample_add = 0
            for id_sample in vect_sample_id_add:
                ## keep track of samples out
                if id_sample in dt_sample_out:
                    continue
                dt_sample_out[id_sample] = 1

                try:
                    sample = Sample.objects.get(pk=id_sample)
                except Sample.DoesNotExist:
                    ## log
                    self.logger_production.error(
                        "Fail to get sample_id {} in ProjectSample".format(
                            key.split("_")[2]
                        )
                    )
                    self.logger_debug.error(
                        "Fail to get sample_id {} in ProjectSample".format(
                            key.split("_")[2]
                        )
                    )
                    continue

                ## the sample can be deleted by other session
                if sample.is_deleted:
                    continue

                ## get project sample
                try:
                    project_sample = PIProject_Sample.objects.get(
                        project__id=project.id, sample__id=sample.id
                    )

                    ### if exist can be deleted, pass to active
                    if (
                        project_sample.is_deleted
                        and not project_sample.is_error
                        and not project_sample.is_deleted_in_file_system
                    ):
                        project_sample.is_deleted = False
                        project_sample.is_finished = False
                        project_sample.save()
                        project_sample_add += 1

                except PIProject_Sample.DoesNotExist:
                    project_sample_input = sample.file_name_1
                    if sample.is_valid_2:
                        project_sample_input += ";" + sample.file_name_2
                    project_sample = PIProject_Sample()
                    project_sample.project = project
                    project_sample.sample = sample
                    project_sample.name = sample.name
                    project_sample.input = project_sample_input
                    project_sample.technology = sample.type_of_fastq
                    project_sample.report = "report"
                    project_sample.save()
                    project_sample_add += 1

            if project_sample_add == 0:
                messages.warning(
                    self.request,
                    "No sample was added to the project '{}'".format(project.name),
                )
            else:
                if project_sample_add > 1:
                    messages.success(
                        self.request,
                        "'{}' samples were added to your project.".format(
                            project_sample_add
                        ),
                        fail_silently=True,
                    )
                else:
                    messages.success(
                        self.request,
                        "One sample was added to your project.",
                        fail_silently=True,
                    )
            return HttpResponseRedirect(reverse_lazy("PIprojects_main"))
        else:
            return super(AddSamples_PIProjectsView, self).form_invalid(form)

    form_valid_message = ""  ## need to have this, even empty


class MainPage(LoginRequiredMixin, generic.CreateView):
    """
    Page with samples of a project
    """

    utils = Utils()
    template_name = "pathogen_identification/main_page.html"
    model = PIProject_Sample
    fields = ["name"]

    def get_context_data(self, **kwargs):
        context = super(MainPage, self).get_context_data(**kwargs)

        try:
            project = Projects.objects.get(pk=self.kwargs["pk"])
        except Projects.DoesNotExist:
            project = None
            messages.error(
                self.request,
                "Project with ID {} does not exist".format(self.kwargs["pk"]),
                fail_silently=True,
            )
            raise Http404

        if project.owner == self.request.user:
            query_set = PIProject_Sample.objects.filter(
                project=project, is_deleted=False
            )
            project_name = project.name
            context["project_owner"] = True

        else:
            messages.error(
                self.request,
                "You do not have permission to access this project.",
            )

            query_set = PIProject_Sample.objects.none()
            project_name = "project"
            context["project_owner"] = False

        tag_search = "search_projects"

        if self.request.GET.get(tag_search) != None and self.request.GET.get(
            tag_search
        ):
            query_string = process_query_string(self.request.GET.get(tag_search))
            query_set = query_set.filter(Q(name__icontains=query_string)).distinct()

        samples = SampleTableOne(query_set)

        ### set the check_box
        if Constants.CHECK_BOX_ALL not in self.request.session:
            self.request.session[Constants.CHECK_BOX_ALL] = False
            is_all_check_box_in_session(
                ["{}_{}".format(Constants.CHECK_BOX, key.id) for key in query_set],
                self.request,
            )

        context[Constants.CHECK_BOX_ALL] = self.request.session[Constants.CHECK_BOX_ALL]
        ## need to clean all the others if are reject in filter
        dt_sample_id_add_temp = {}
        if context[Constants.CHECK_BOX_ALL]:
            for sample in query_set:
                dt_sample_id_add_temp[sample.id] = (
                    1  ## add the ids that are in the tables
                )
            for key in self.request.session.keys():
                if (
                    key.startswith(Constants.CHECK_BOX)
                    and len(key.split("_")) == 3
                    and self.utils.is_integer(key.split("_")[2])
                ):
                    ### this is necessary because of the search. Can occur some checked box that are out of filter.
                    if int(key.split("_")[2]) not in dt_sample_id_add_temp:
                        self.request.session[key] = False
                    else:
                        self.request.session[key] = True
        ## END need to clean all the others if are reject in filter
        ### set the check_box
        RequestConfig(
            self.request, paginate={"per_page": Constants.PAGINATE_NUMBER}
        ).configure(samples)

        ### type of deployment
        DEPLOY_TYPE = PICS.DEPLOYMENT_DEFAULT
        DEPLOY_URL = "deploy_ProjectPI"

        if DEPLOY_TYPE == PICS.DEPLOYMENT_TYPE_PIPELINE:
            DEPLOY_URL = "deploy_runs_ProjectPI"

        context["rows_color"] = [
            "combinations",
            "mapping_runs",
            "running_processes",
            "queued_processes",
        ]

        context["metagenomics"] = ConstantsSettings.METAGENOMICS
        context["table"] = samples
        context["deploy_url"] = DEPLOY_URL
        context["user_id"] = project.owner.pk
        context["project_index"] = project.pk
        context["project_name"] = project_name
        context["nav_sample"] = True
        context["total_items"] = query_set.count()
        context["show_paginatior"] = query_set.count() > Constants.PAGINATE_NUMBER
        context["show_info_main_page"] = ShowInfoMainPage()
        context["query_set_count"] = query_set.count()
        context["demo"] = True if self.request.user.username == "demo" else False

        return context


def excise_paths_leaf_last(string_with_paths):
    """
    if the string has paths, find / and return the last
    """
    split_space = string_with_paths.split(" ")
    new_string = ""
    if len(split_space) > 1:
        for word in split_space:
            new_word = word
            if "/" in word:
                new_word = word.split("/")[-1]
            new_string += new_word + " "
        return new_string
    else:
        return string_with_paths


def teleflu_node_info(node, params_df, node_pk):
    node_info = {
        "pk": node_pk,
        "node": node,
        "modules": [],
    }

    for pipeline_step in CS.vect_pipeline_televir_workflows_display:
        acronym = [x[0] for x in pipeline_step.split(" ")]
        acronym = "".join(acronym).upper()
        params = params_df[params_df.module == pipeline_step].to_dict("records")
        if params:  # if there are parameters for this module
            software = params[0].get("software")
            software = software.split("_")[0]

            params = params[0].get("value")
            params = excise_paths_leaf_last(params)
            params = f"{software} {params}"
        else:
            params = ""

        node_info["modules"].append(
            {
                "module": pipeline_step,
                "parameters": params,
                "short_name": acronym,
                "available": (
                    "software_on"
                    if pipeline_step in params_df.module.values
                    else "software_off"
                ),
            }
        )

    return node_info


class TelefluProjectView(LoginRequiredMixin, generic.CreateView):
    """
    Teleflu Project
    """

    template_name = "pathogen_identification/teleflu.html"
    model = TeleFluProject
    fields = ["name"]

    def get_context_data(self, **kwargs):
        context = super(TelefluProjectView, self).get_context_data(**kwargs)

        teleflu_project_pk = int(self.kwargs["pk"])
        ## TeleFlu Projects
        this_project = TeleFluProject.objects.get(pk=teleflu_project_pk)
        teleflu_samples = TeleFluSample.objects.filter(teleflu_project=this_project)

        teleflu_projects = TeleFluProject.objects.filter(
            pk=teleflu_project_pk, is_deleted=False
        ).order_by("-last_change_date")
        televir_project = teleflu_projects[0].televir_project
        user = televir_project.owner

        context["insaflu_table"] = None
        context["insaflu_projects_exist"] = teleflu_projects.exists()

        if teleflu_projects.exists():

            context["insaflu_table"] = TeleFluInsaFLuProjectTable(teleflu_projects)
            RequestConfig(
                self.request, paginate={"per_page": Constants.PAGINATE_NUMBER}
            ).configure(context["insaflu_table"])

        ##################################### Existring mappings
        software_utils = SoftwareTreeUtils(user, televir_project)

        mappings = TelefluMapping.objects.filter(teleflu_project__pk=teleflu_project_pk)
        mapping_workflows = []
        existing_mapping_pks = []
        for mapping in mappings:
            if mapping.leaf is None:
                continue
            existing_mapping_pks.append(mapping.leaf.pk)

        context["mapping_workflows"] = mapping_workflows
        ####################################### get combinations to deploy
        local_tree = software_utils.generate_software_tree_safe(
            software_utils.project,
            software_utils.sample,
            metagenomics=False,
            mapping_only=True,
            screening=False,
        )

        if local_tree.makeup == -1:
            all_paths = {}
            available_path_nodes = {}
        else:
            all_paths = local_tree.get_all_graph_paths()
            available_path_nodes = software_utils.get_available_pathnodes(local_tree)
        ##########################################
        workflows = []
        for node, params_df in all_paths.items():
            if node not in available_path_nodes:
                continue
            if available_path_nodes[node].pk in existing_mapping_pks:
                continue

            node_pk = available_path_nodes[node].pk
            node_info = teleflu_node_info(node, params_df, node_pk)
            workflows.append(node_info)

        context["insaflu_connection_exists"] = (
            False if this_project.insaflu_project is None else True
        )

        context["workflows"] = workflows
        context["project_nsamples"] = teleflu_samples.count()
        context["project_index"] = televir_project.pk
        context["teleflu_project_pk"] = teleflu_project_pk
        context["project_name"] = televir_project.name
        context["focus_teleflu"] = (
            f"Focus: {this_project.raw_reference.description_first}"
        )
        ###

        return context


from fluwebvirus.settings import (BASE_DIR, MEDIA_ROOT, MEDIA_URL, STATIC_ROOT,
                                  STATIC_URL)
from pathogen_identification.utilities.reference_utils import \
    filter_reference_maps_select
from pathogen_identification.utilities.utilities_general import simplify_name


def remove_pre_static(path: str) -> str:
    cwd = os.getcwd()
    if path.startswith(cwd):
        path = path[len(cwd) :]

    path = path.replace(STATIC_ROOT, STATIC_URL)
    path = path.replace(MEDIA_ROOT, MEDIA_URL)
    # path = mark_safe(path)
    return path


class TelefluMappingIGV(LoginRequiredMixin, generic.TemplateView):
    """
    Teleflu Mapping IGV
    """

    template_name = "pathogen_identification/teleflu_mapping_igv.html"

    def get_context_data(self, **kwargs):
        context = super(TelefluMappingIGV, self).get_context_data(**kwargs)
        televir_bioinf = TelevirBioinf()

        teleflu_mapping_pk = int(self.kwargs["pk"])
        teleflu_mapping = TelefluMapping.objects.get(pk=teleflu_mapping_pk)
        leaf_index = teleflu_mapping.leaf.index
        teleflu_project = teleflu_mapping.teleflu_project
        televir_project_index = teleflu_project.televir_project.pk

        insaflu_project = teleflu_project.insaflu_project

        ### get reference
        teleflu_reference = teleflu_project.raw_reference
        if teleflu_reference is None:
            return False

        reference_file = teleflu_reference.file_path
        reference_index = reference_file + ".fai"
        if os.path.exists(reference_index) is False:
            televir_bioinf.index_fasta(reference_file)
        reference_file = remove_pre_static(reference_file)
        reference_index = remove_pre_static(reference_index)
        # televir_reference
        teleflu_refs = teleflu_project.televir_references

        igv_genome_options = {
            "reference": reference_file,
            "reference_index": reference_index,
            "reference_name": teleflu_mapping.teleflu_project.raw_reference.description,
        }

        # samples
        televir_project_samples = teleflu_mapping.mapped_samples
        sample_dict = {}

        ### get sample files
        accid_list = [ref.accid for ref in teleflu_refs if ref.accid]
        accid_list_simple = [simplify_name(accid) for accid in accid_list] + accid_list

        for sample in televir_project_samples:

            ref_select = filter_reference_maps_select(
                sample, leaf_index, accid_list_simple
            )

            if ref_select is None:
                continue

            sample_dict[sample.pk] = {
                "name": sample.name,
                "bam_file": remove_pre_static(ref_select.bam_file_path),
                "bam_file_index": remove_pre_static(ref_select.bai_file_path),
                "vcf_file": remove_pre_static(ref_select.vcf),
                "sample": sample,
            }

        context["igv_genome"] = igv_genome_options
        context["samples"] = sample_dict
        context["project_index"] = televir_project_index

        return context


class Sample_main(LoginRequiredMixin, generic.CreateView):
    """
    sample main page with list runs per sample
    """

    template_name = "pathogen_identification/sample_main.html"
    model = RunMain
    fields = ["name"]

    def get_context_data(self, **kwargs):
        context = super(Sample_main, self).get_context_data(**kwargs)

        project_pk = int(self.kwargs["pk1"])
        sample_pk = int(self.kwargs["pk2"])
        user = self.request.user

        try:
            project = Projects.objects.get(pk=project_pk)
            sample = PIProject_Sample.objects.get(pk=sample_pk)
        except Exception:
            messages.error(
                self.request,
                "Project with ID '{}' does not exist".format(project_pk),
                fail_silently=True,
            )
            raise Http404

        if project.owner == self.request.user:
            runs = RunMain.objects.filter(
                sample__pk=sample_pk,
                project__pk=project_pk,
                project__owner=user,
                parameter_set__status__in=[
                    ParameterSet.STATUS_FINISHED,
                    ParameterSet.STATUS_RUNNING,
                ],
                run_type__in=[
                    RunMain.RUN_TYPE_PIPELINE,
                ],
            ).order_by("-parameter_set__leaf__index")

            run_mapping = RunMain.objects.filter(
                sample__pk=sample_pk,
                project__pk=project_pk,
                project__owner=user,
                parameter_set__status__in=[
                    ParameterSet.STATUS_KILLED,
                    ParameterSet.STATUS_FINISHED,
                    ParameterSet.STATUS_RUNNING,
                    ParameterSet.STATUS_ERROR,
                ],
                run_type__in=[
                    RunMain.RUN_TYPE_MAP_REQUEST,
                    RunMain.RUN_TYPE_COMBINED_MAPPING,
                    RunMain.RUN_TYPE_PANEL_MAPPING,
                ],
                status__in=[
                    RunMain.STATUS_DEFAULT,
                    RunMain.STATUS_RUNNING,
                    RunMain.STATUS_FINISHED,
                ],
            ).order_by("-created_in")
            sample_name = sample.sample.name
            project_name = project.name

        else:
            messages.error(
                self.request,
                "You do not have permission to access this project.",
                fail_silently=True,
            )
            runs = RunMain.objects.none()
            run_mapping = RunMain.objects.none()
            sample_name = "sample"
            project_name = "project"

        runs_table = RunMainTable(runs, exclude=("created", "nmapped", "mapping"))
        rendered_table = ""

        if run_mapping.exists():
            run_mappings_table = RunMappingTable(run_mapping, order_by=("created",))
            small_context = {
                "nav_sample": True,
                "total_items": run_mapping.count(),
                "show_paginatior": run_mapping.count()
                > ConstantsSettings.PAGINATE_NUMBER,
                "show_info_main_page": ShowInfoMainPage(),
                "table": run_mappings_table,
                "query_set_count": run_mapping.count(),
            }

            template_table_html = "pathogen_identification/mapping_runs.html"
            rendered_table = render_to_string(
                template_table_html, small_context, request=self.request
            )

        RequestConfig(
            self.request, paginate={"per_page": ConstantsSettings.PAGINATE_NUMBER}
        ).configure(runs_table)

        context = {
            "nav_sample": True,
            "total_items": runs.count(),
            "show_paginatior": runs.count() > ConstantsSettings.PAGINATE_NUMBER,
            "show_info_main_page": ShowInfoMainPage(),
            "table": runs_table,
            "table_mapping": rendered_table,
            "mappings_exist": run_mapping.exists(),
            "sample_name": sample_name,
            "project_main": True,
            "project_name": project_name,
            "project_index": project_pk,
            "sample_index": sample_pk,
            "query_set_count": runs.count(),
        }

        return context


from django.http import JsonResponse
from django.template.loader import render_to_string

from fluwebvirus.settings import BASE_DIR


def inject_references_filter(request, max_references: int = 30):
    ###
    tag_add_reference = "search_add_project_reference"
    tag_teleflu = "teleflu_reference"
    table_type = "add_reference"
    tag_panel_id = "panel_id"
    panel = None
    panel_id = None
    project_id = None
    project = None
    user = request.user

    if request.GET.get("max_references") and request.GET.get("max_references") != "":
        max_references = int(request.GET.get("max_references"))

    if request.GET.get(tag_teleflu) and request.GET.get(tag_teleflu) != "":
        table_type = "teleflu_reference"
        max_references = 10
        project_id = int(request.GET.get("project_id"))
        project = Projects.objects.get(pk=project_id)

    if request.GET.get(tag_panel_id) and request.GET.get(tag_panel_id) != "":
        panel_id = int(request.GET.get(tag_panel_id))
        panel = ReferencePanel.objects.get(pk=panel_id)

    references = []

    if request.GET.get(tag_add_reference) is not None:
        if table_type == "teleflu_reference":
            request_tag = request.GET.get(tag_add_reference)
            possible_accid = request_tag.split(".")[0]
            try:
                references = RawReference.objects.filter(
                    Q(accid__icontains=request.GET.get(tag_add_reference))
                    | Q(accid__icontains=possible_accid)
                    | Q(description__icontains=request.GET.get(tag_add_reference))
                    & ~Q(description__in=["root", "NA"])
                    & Q(run__project__pk=project_id)
                ).distinct("accid")

            except Exception as e:
                print(e)

            if references.count() == 0:
                references = []
            else:
                reference_utils = RawReferenceUtils(project=project)
                query_string = request.GET.get(tag_add_reference)

                references = reference_utils.retrieve_compound_references(query_string)

        elif request.GET.get(tag_add_reference) != "":

            existing_reference_taxids = []

            if panel is not None:
                existing_reference_taxids = RawReference.objects.filter(
                    panel=panel
                ).values_list("taxid", flat=True)

            try:
                references = (
                    ReferenceSourceFileMap.objects.filter(
                        Q(
                            reference_source__description__icontains=request.GET.get(
                                tag_add_reference
                            )
                        )
                        | Q(
                            reference_source__accid__icontains=request.GET.get(
                                tag_add_reference
                            )
                        )
                        | Q(
                            reference_source__taxid__taxid__icontains=request.GET.get(
                                tag_add_reference
                            )
                        )
                        | Q(
                            reference_source_file__file__icontains=request.GET.get(
                                tag_add_reference
                            ),
                            reference_source_file__owner__in=[
                                None,
                                user,
                            ],  # allow public references only
                        )
                    )
                    .exclude(
                        reference_source__taxid__taxid__in=existing_reference_taxids
                    )
                    .distinct("reference_source__accid")
                )
            except Exception as e:
                print(e)

        # show max 10 references
        references = references[:max_references]

        data = inject_references(references, request, table_type=table_type)

        return JsonResponse(data)

    else:
        data = inject_references([], request)
        return JsonResponse(data)


def inject_references(references: list, request, table_type: str = "add_reference"):
    context = {}
    data = {}

    if table_type == "add_reference":
        table = ReferenceSourceTable(references)
    else:
        table = TeleFluReferenceTable(references, order_by=("global_score",))

    context["references_table"] = table
    context["references_count"] = len(references)

    if table_type == "add_reference":
        template_table_html = "pathogen_identification/references_table_table_only.html"
    else:
        template_table_html = (
            "pathogen_identification/teleflu_references_table_only.html"
        )
    # render tamplate using context
    try:
        rendered_table = render_to_string(template_table_html, context, request=request)
    except Exception as e:
        print(e)
        rendered_table = ""
    data["my_content"] = rendered_table
    data["references_count"] = len(references)

    return data


def inject_references_added_html(request):
    sample_pk = int(request.GET.get("sample_id"))
    query_set_added_manual = RawReference.objects.filter(
        run__sample__pk=sample_pk, run__run_type=RunMain.RUN_TYPE_STORAGE
    )

    added_references_context = inject__added_references(query_set_added_manual, request)
    empty_data_context = inject_references([], request)
    added_references_context["empty_content"] = empty_data_context["my_content"]

    return JsonResponse(added_references_context)


def inject__added_references(references: list, request):
    context = {}
    data = {}

    context["references_table"] = AddedReferenceTable(references)
    context["references_count"] = len(references)

    template_table_html = (
        "pathogen_identification/references_table_table_only_padding.html"
    )

    # render tamplate using context
    rendered_table = render_to_string(template_table_html, context, request=request)
    data["my_content"] = rendered_table
    data["references_count"] = len(references)

    return data


# class TeleFluProjectCreate(LoginRequiredMixin, generic.CreateView):
class ReferencePanelManagement(LoginRequiredMixin, generic.CreateView):
    """
    page to manage and create insaflu references.

    """

    template_name = "pathogen_identification/reference_panel_management.html"
    fields = "__all__"
    utils = Utils()

    def get_queryset(self, **kwargs):
        user_pk = self.request.user.pk

        return (
            ReferencePanel.objects.filter(owner__id=user_pk)
            .exclude(is_deleted=True)
            .order_by("-creation_date")
        )

    def get_context_data(self, **kwargs) -> Dict[str, Any]:
        context = super().get_context_data(**kwargs)
        user = self.request.user

        panels = (
            ReferencePanel.objects.filter(
                project_sample=None,
                owner__id=user.pk,
                is_deleted=False,
                panel_type=ReferencePanel.PANEL_TYPE_MAIN,
            )
            .exclude(is_deleted=True)
            .order_by("-creation_date")
        )
        context["panels"] = panels
        context["user_id"] = user.pk
        context["nav_reference"] = True

        return context


from django.views.generic import ListView, TemplateView

from pathogen_identification.tables import (ReferenceSourceFileTable,
                                            TelevirReferencesTable)


class ReferenceManagementBase(TemplateView):
    """
    page to manage and create insaflu references files, generate panels.
    """

    template_name = "pathogen_identification/televir_references_base.html"

    def get_context_data(self, **kwargs) -> Dict[str, Any]:

        context = super().get_context_data(**kwargs)
        context["nav_reference"] = True
        return context


class ReferenceFileManagement(LoginRequiredMixin, generic.CreateView):
    """
    page to manage and create insaflu references files.

    """

    template_name = "pathogen_identification/televir_files_management.html"
    fields = "__all__"
    utils = Utils()

    def get_queryset(self, **kwargs):
        user_pk = self.request.user.pk

        return (
            ReferenceSourceFile.objects.filter(Q(owner=None) | Q(owner__id=user_pk))
            .exclude(is_deleted=True)
            .order_by("-creation_date")
        )

    def get_context_data(self, **kwargs) -> Dict[str, Any]:
        context = super().get_context_data(**kwargs)
        user = self.request.user

        files = (
            ReferenceSourceFile.objects.filter(Q(owner=None) | Q(owner__id=user.pk))
            .exclude(is_deleted=True)
            .order_by("-owner", "-creation_date")
        )

        files_table = ReferenceSourceFileTable(files)
        RequestConfig(self.request, paginate={"per_page": 15}).configure(files_table)

        context["files_table"] = files_table
        context["nav_reference"] = True
        context["show_paginatior"] = files.count() > 15
        context["query_set_count"] = files.count()
        context["user_id"] = user.pk

        return context


class ReferenceManagement(LoginRequiredMixin, generic.CreateView):
    """
    page to manage and create insaflu references files.

    """

    template_name = "pathogen_identification/televir_references_view.html"
    fields = "__all__"
    utils = Utils()

    def get_queryset(self, **kwargs):
        user_pk = self.request.user.pk

        return ReferenceSourceFileMap.objects.filter(
            Q(reference_source_file__owner=None)
            | Q(reference_source_file__owner__id=user_pk)
        )

    def get_context_data(self, **kwargs) -> Dict[str, Any]:
        context = super().get_context_data(**kwargs)
        user = self.request.user

        tag_search = "search_references"

        references = (
            ReferenceSourceFileMap.objects.filter(
                Q(reference_source_file__owner__id=user.pk)
                # | Q(reference_source_file__owner=None)
            )
            .order_by("reference_source__description", "reference_source__accid")
            .distinct("reference_source__description", "reference_source__accid")
        )

        if self.request.GET.get(tag_search) != None and self.request.GET.get(
            tag_search
        ):
            query_string = self.request.GET.get(tag_search)
            query_string = process_query_string(query_string)

            references = (
                ReferenceSourceFileMap.objects.filter(
                    Q(reference_source_file__owner__id=user.pk)
                    | Q(reference_source_file__owner=None)
                )
                .order_by("reference_source__description", "reference_source__accid")
                .distinct("reference_source__description", "reference_source__accid")
            )

            references = references.filter(
                Q(
                    reference_source__description__icontains=self.request.GET.get(
                        tag_search
                    )
                )
                | Q(reference_source__accid__icontains=query_string)
                | Q(
                    reference_source__taxid__taxid__icontains=self.request.GET.get(
                        tag_search
                    )
                )
                | Q(
                    reference_source_file__file__icontains=self.request.GET.get(
                        tag_search
                    )
                )
            )

        summary = {
            "total": references.count(),
            "accession_id": references.values_list("reference_source__accid", flat=True)
            .distinct()
            .count(),
            "taxID": references.values_list("reference_source__taxid__taxid", flat=True)
            .distinct()
            .count(),
            "description": references.values_list(
                "reference_source__description", flat=True
            )
            .distinct()
            .count(),
            "Source": references.values_list("reference_source_file__file", flat=True)
            .distinct()
            .count(),
        }

        files_table = TelevirReferencesTable(references, user_id=user.id)
        RequestConfig(
            self.request,
            paginate={"per_page": ConstantsSettings.TELEVIR_REFERENCE_PAGINATE_NUMBER},
        ).configure(files_table)

        context["summary"] = summary
        context["files_table"] = files_table
        context["nav_reference"] = True
        context["show_paginatior"] = references.count() > Constants.PAGINATE_NUMBER
        context["query_set_count"] = references.count()
        context["user_id"] = user.pk

        return context


def check_metadata_table_clean(metadata_table_file) -> Optional[pd.DataFrame]:
    """
    check metadata table
    """

    metadata_table = metadata_table_file.read().decode("utf-8")
    metadata_table = metadata_table.split("\n")
    sep = os.path.splitext(metadata_table_file.name)[1]
    sep = "\t" if sep == ".tsv" else ","
    metadata_table = [x.split(sep) for x in metadata_table]
    metadata_table = [x for x in metadata_table if len(x) > 1]

    try:
        metadata_table = pd.DataFrame(metadata_table[1:], columns=metadata_table[0])

    except Exception as e:
        print(e)
        return None

    if len(metadata_table) == 0:
        return None
    return metadata_table


def download_template_view(request):
    file_path = os.path.join(
        BASE_DIR, "templates", "televir_reference_metadata_template.tsv"
    )
    with open(file_path, "rb") as fh:
        response = HttpResponse(fh.read(), content_type="application/vnd.ms-excel")
        response["Content-Disposition"] = "inline; filename=" + os.path.basename(
            file_path
        )
        return response


class UploadReferencePanel(LoginRequiredMixin, FormValidMessageMixin, generic.FormView):
    """
    page to manage and create insaflu references files.

    """

    template_name = "pathogen_identification/televir_upload_panels.html"
    success_url = reverse_lazy("televir_reference_files")
    form_class = UploadFileForm

    def get_form_kwargs(self):
        """ """

        kw = super(UploadReferencePanel, self).get_form_kwargs()
        kw["request"] = self.request
        return kw

    def get_context_data(self, **kwargs):
        context = super(UploadReferencePanel, self).get_context_data(
            **kwargs
        )  # Call the base implementation first to get a context
        if "form" in kwargs and hasattr(kwargs["form"], "error_in_file"):
            context["error_in_file"] = mark_safe(
                kwargs["form"].error_in_file.replace("\n", "<br>")
            )  ## pass a list
        context["nav_reference"] = True
        context["nav_modal"] = True  ## short the size of modal window
        context["show_info_main_page"] = (
            ShowInfoMainPage()
        )  ## show main information about the institute
        return context

    def form_valid(self, form):
        try:
            profile = Profile.objects.get(user=self.request.user)
            if profile.only_view_project:
                messages.warning(
                    self.request,
                    "'{}' account can not add file with samples.".format(
                        self.request.user.username
                    ),
                    fail_silently=True,
                )
                return super(UploadReferencePanel, self).form_invalid(form)
        except Profile.DoesNotExist:
            pass

        data = {"is_error": False, "is_ok": False, "error_message": ""}
        # Handle the uploaded files here

        # name = form.cleaned_data["name"]
        description = form.cleaned_data["description"]
        reference_fasta_file = form.cleaned_data["fasta_file"]
        metadata_file = form.cleaned_data["metadata"]
        ###
        software = Software()

        reference_metadata_table = check_metadata_table_clean(metadata_file)
        reference_fasta_temp_file_name = NamedTemporaryFile(
            prefix="flu_fa_", delete=False
        )
        reference_metadata_temp_file_name = NamedTemporaryFile(
            prefix="flu_fa_", delete=False, suffix=".tsv"
        )

        try:
            file_data = reference_fasta_file.read()
            reference_fasta_temp_file_name.write(file_data)
            reference_fasta_temp_file_name.flush()
            reference_fasta_temp_file_name.close()
            software.dos_2_unix(reference_fasta_temp_file_name.name)
        except Exception as e:

            data["is_error"] = True
            messages.error(
                self.request,
                "Error in the fasta file",
            )
            return super(UploadReferencePanel, self).form_invalid(form)

        try:
            reference_metadata_table.to_csv(
                reference_metadata_temp_file_name.name, sep="\t", index=False
            )
        except Exception as e:
            print(e)

            messages.error(
                self.request,
                "Error in the metadata file",
            )
            return super(UploadReferencePanel, self).form_invalid(form)

        process_SGE = ProcessSGE()

        try:
            # create reference source file
            reference_source_file = ReferenceSourceFile()
            reference_source_file.owner = self.request.user
            reference_source_file.file = reference_fasta_file
            reference_source_file.description = description
            reference_source_file.save()

            taskID = process_SGE.set_submit_upload_reference_televir(
                user=self.request.user,
                file_id=reference_source_file.pk,
                fasta=reference_fasta_temp_file_name.name,
                metadata=reference_metadata_temp_file_name.name,
            )

        except Exception as e:
            data["is_error"] = True
            messages.error(
                self.request,
                "Error saving file",
            )
            return super(UploadReferencePanel, self).form_invalid(form)

        messages.success(
            self.request,
            "Reference file '{}' uploaded successfully".format(description),
        )
        return super(UploadReferencePanel, self).form_valid(form)

    form_valid_message = ""  ## need to have this


class ReferencesManagementSample(LoginRequiredMixin, generic.CreateView):
    """
    page with raw references table for single sample. used to add references and select for remap
    """

    template_name = "pathogen_identification/references_table.html"
    fields = "__all__"
    utils = Utils()

    def get_queryset(self, **kwargs):
        sample_pk = int(self.kwargs["pk1"])

        return (
            RawReference.objects.filter(run__sample__pk=sample_pk)
            .exclude(accid="-")
            .distinct("accid")
        )

    def get_context_data(self, **kwargs):
        context = super(ReferencesManagementSample, self).get_context_data(**kwargs)
        context["nav_project"] = True

        sample_pk = int(self.kwargs["pk1"])

        try:
            sample_main = PIProject_Sample.objects.get(pk=sample_pk)
        except Exception:
            messages.error(
                self.request,
                "Sample with ID '{}' does not exist".format(sample_pk),
                fail_silently=True,
            )
            raise Http404

        project_main = sample_main.project
        project_name = project_main.name

        if project_main.owner == self.request.user:
            query_set = (
                RawReference.objects.filter(run__sample__pk=sample_pk)
                .exclude(accid="-")
                .distinct("accid")
            )

            query_set_added_manual = query_set.filter(
                run__run_type=RunMain.RUN_TYPE_STORAGE
            )
            query_set = query_set.exclude(run__run_type=RunMain.RUN_TYPE_STORAGE)

            sample_name = sample_main.sample.name
        else:
            messages.error(
                self.request,
                "You do not have permission to access this project.",
                fail_silently=True,
            )
            query_set = RawReference.objects.none()
            query_set_added_manual = RawReference.objects.none()
            project_name = "project"
            sample_name = "sample"

        #### get panels
        panels = (
            ReferencePanel.objects.filter(
                owner__id=self.request.user.pk,
                is_deleted=False,
                project_sample=sample_main,
                panel_type=ReferencePanel.PANEL_TYPE_MAIN,
            )
            .exclude(is_deleted=True)
            .order_by("-creation_date")
        )
        context["panels"] = panels
        #### search add reference bar.
        context["references_table"] = None
        context["references_count"] = 0

        if self.request.method == "POST":
            references_form = ReferenceForm(self.request.POST)
        else:
            references_form = ReferenceForm()

        context["references_form"] = references_form

        ##### check Screening performed
        screening_performed = RunMain.objects.filter(
            sample=sample_main, run_type=RunMain.RUN_TYPE_SCREENING
        ).exists()

        reference_table_class = CompoundReferenceScore
        ordered_by = ("-standard_score",)

        if screening_performed:
            reference_table_class = CompoundRefereceScoreWithScreening
            ordered_by = ("-standard_score", "screening_score")

        #### added references table
        added_references_context = inject__added_references(
            query_set_added_manual, self.request
        )
        context["added_references_table"] = added_references_context["my_content"]
        context["added_references_count"] = added_references_context["references_count"]

        ##### search add reference bar

        # raw_references = raw_references
        tag_search = "search_add_project_sample"
        query_string = None

        if self.request.GET.get(tag_search) != None and self.request.GET.get(
            tag_search
        ):
            query_string = self.request.GET.get(tag_search)
            query_string = process_query_string(query_string)
            query_set = query_set.filter(
                Q(description__icontains=query_string)
                | Q(accid__icontains=query_string)
                | Q(taxid__icontains=query_string)
            )
            # tag_search

        #####
        reference_utils = RawReferenceUtils(sample_main)

        raw_reference_compound = reference_utils.retrieve_compound_references(
            query_string
        )

        # pks of classification runs as integer list
        # runs_pks = [run.pk for run in classification_runs]
        classification_runs = RunMain.objects.filter(
            sample=sample_main, run_type=RunMain.RUN_TYPE_PIPELINE
        )
        if classification_runs.exists():
            compound_reference_table = reference_table_class(
                raw_reference_compound, order_by=ordered_by
            )
        else:
            compound_reference_table = reference_table_class(
                raw_reference_compound, order_by=ordered_by
            )

        if Constants.CHECK_BOX_ALL not in self.request.session:
            self.request.session[Constants.CHECK_BOX_ALL] = False
            is_all_check_box_in_session(
                ["{}_{}".format(Constants.CHECK_BOX, key.id) for key in query_set],
                self.request,
            )

        context[Constants.CHECK_BOX_ALL] = self.request.session[Constants.CHECK_BOX_ALL]
        ## need to clean all the others if are reject in filter
        dt_sample_id_add_temp = {}
        if context[Constants.CHECK_BOX_ALL]:
            for sample in query_set:
                dt_sample_id_add_temp[sample.id] = (
                    1  ## add the ids that are in the tables
                )
            for key in self.request.session.keys():
                if (
                    key.startswith(Constants.CHECK_BOX)
                    and len(key.split("_")) == 3
                    and self.utils.is_integer(key.split("_")[2])
                ):
                    ### this is necessary because of the search. Can occur some checked box that are out of filter.
                    if int(key.split("_")[2]) not in dt_sample_id_add_temp:
                        self.request.session[key] = False
                    else:
                        self.request.session[key] = True

        RequestConfig(self.request, paginate={"per_page": 20}).configure(
            compound_reference_table
        )
        if self.request.GET.get(tag_search) != None:
            context[tag_search] = self.request.GET.get(tag_search)

        ### graph
        #### graph
        graph_progress = TreeProgressGraph(sample_main)
        # graph_progress.generate_graph()
        graph_json, graph_id = graph_progress.get_graph_data()
        ####
        # runs = set([run.pk for run in classification_runs])
        # runs = classification_runs
        runs_number = len(classification_runs) > 0
        ### modal buttons
        # white color icons for the buttons
        color = 'style="color:white;"'
        deploy_metagenomics = (
            "<a "
            + 'href="#id_deploy_metagenomics_modal" data-toggle="modal" data-toggle="tooltip" '
            + 'id="deploy_metagenomics_modal" '
            + 'title="Run combined metagenomics"'
            + f"pk={sample_pk} "
            + f"ref_name={sample_name} style='color: #fff;'"
            + f">Map Combined </a>"
        )
        metagenomics_parameters = (
            "<a href="
            + reverse("pathogenID_sample_settings", kwargs={"sample": sample_pk})
            + ' data-toggle="tooltip" data-toggle="modal" title="Manage combine settings" style="color: #fff;">'
            + f'<i class="padding-button-table fa fa-pencil-square padding-button-table" {color}></i> Parameters </a>'
        )

        context["runs"] = classification_runs
        context["runs_number"] = runs_number
        context["graph_json"] = graph_json
        context["graph_id"] = graph_id
        context["meta_parameters"] = mark_safe(metagenomics_parameters)
        context["deploy_metagenomics"] = mark_safe(deploy_metagenomics)
        context["nav_sample"] = True
        context["show_paginatior"] = query_set.count() > Constants.PAGINATE_NUMBER
        context["query_set_count"] = query_set.count()
        context["show_info_main_page"] = ShowInfoMainPage()
        context["nav_modal"] = True  ## short the size of modal window

        context["owner"] = True
        context["references"] = raw_reference_compound
        context["table"] = compound_reference_table
        context["sample_name"] = sample_name
        context["project_name"] = project_name
        context["project_index"] = project_main.pk
        context["sample_index"] = sample_pk

        return context


def Project_reports(requesdst, pk1):
    """
    Projects reports
    """
    template_name = "pathogen_identification/allreports_table.html"

    try:
        project = Projects.objects.get(pk=int(pk1))
    except Projects.DoesNotExist:
        messages.error(
            requesdst,
            "Project does not exist",
            fail_silently=True,
        )
        raise Http404

    if project.owner == requesdst.user:
        all_reports = FinalReport.objects.filter(run__project__pk=int(pk1)).order_by(
            "-coverage"
        )
        project_name = project.name

    else:
        messages.error(
            requesdst,
            "You do not have permission to access this project.",
            fail_silently=True,
        )
        all_reports = FinalReport.objects.none()
        project_name = "project"

    return render(
        requesdst,
        template_name,
        {
            "final_report": all_reports,
            "project": project_name,
            "project_index": project.pk,
        },
    )


def Sample_reports(requesdst, pk1, pk2):
    """
    sample main page
    """
    template_name = "pathogen_identification/allreports_table.html"

    try:
        project = Projects.objects.get(pk=int(pk1))
    except Projects.DoesNotExist:
        messages.error(
            requesdst,
            "Project does not exist.",
            fail_silently=True,
        )
        raise Http404

    if project.owner == requesdst.user:
        all_reports = FinalReport.objects.filter(
            run__project__pk=int(pk1), sample__pk=int(pk2)
        ).order_by("-coverage")
        project_name = project.name
        sample_name = PIProject_Sample.objects.get(pk=int(pk2)).sample.name

    else:
        messages.error(
            requesdst,
            "You do not have permission to access this project.",
            fail_silently=True,
        )
        all_reports = FinalReport.objects.none()
        project_name = "project"
        sample_name = "sample"

    return render(
        requesdst,
        template_name,
        {
            "final_report": all_reports,
            "project": project_name,
            "sample": sample_name,
            "sample_index": pk2,
            "project_index": project.pk,
        },
    )


class Sample_detail(LoginRequiredMixin, generic.CreateView):
    """
    home page
    """

    template_name = "pathogen_identification/sample_detail.html"
    utils = Utils()

    def get_context_data(self, **kwargs):
        project_pk = int(self.kwargs["pk1"])
        sample_pk = int(self.kwargs["pk2"])
        run_pk = int(self.kwargs["pk3"])

        try:
            project_main = Projects.objects.get(pk=project_pk)

        except Exception as e:
            messages.error(self.request, "Project does not exist")
            raise Http404

        try:
            sample = PIProject_Sample.objects.get(pk=sample_pk)
        except Exception as e:
            messages.error(self.request, "Sample does not exist")
            raise Http404

        try:
            run_main_pipeline = RunMain.objects.get(pk=run_pk)

        except Exception as e:
            messages.error(self.request, "Run does not exist")
            raise Http404

        if self.request.user != project_main.owner:
            messages.error(
                self.request,
                "You do not have permission to access this project.",
                fail_silently=True,
            )

            return {"owner": False}

        project_name = project_main.name
        sample_name = sample.name
        run_name = run_main_pipeline.parameter_set.leaf.index
        sample_main = run_main_pipeline.sample
        #
        is_classification = run_main_pipeline.run_type == RunMain.RUN_TYPE_PIPELINE
        #

        raw_references = (
            RawReference.objects.filter(run=run_main_pipeline)
            .exclude(accid="-")
            .exclude(counts=None)
            .order_by("taxid", "status")
            .distinct("taxid")
        )
        raw_references = sorted(
            raw_references,
            key=lambda x: float(x.read_counts if x.read_counts else 0),
            reverse=True,
        )

        ########
        remapping_performed = True
        if run_main_pipeline.run_type == RunMain.RUN_TYPE_PIPELINE:
            parameter_utils = Parameter_DB_Utility()
            remapping_performed = parameter_utils.check_parameter_set_contains_module(
                run_main_pipeline.parameter_set.leaf, CS.PIPELINE_NAME_remapping
            )

        if is_classification is True:
            raw_reference_table = RawReferenceTable(raw_references)

        else:
            raw_reference_table = RawReferenceTable_Basic(raw_references)

        #
        run_detail = RunDetail.objects.get(sample=sample_main, run=run_main_pipeline)

        #
        try:
            run_qc = TelevirRunQC.objects.get(run=run_main_pipeline)
            qc_report = RunQC_report(
                performed=run_qc.performed,
                method=run_qc.method,
                args=run_qc.args,
                input_reads=run_qc.input_reads,
                output_reads=run_qc.output_reads,
                output_reads_percent=run_qc.output_reads_percent,
            )

        except TelevirRunQC.DoesNotExist:
            qc_report = RunQC_report(
                performed=False,
                method="None",
                args="None",
                input_reads=run_detail.input,
                output_reads=run_detail.input,
                output_reads_percent="1",
            )
        #
        try:
            run_assembly = RunAssembly.objects.get(
                sample=sample_main, run=run_main_pipeline
            )
            recover_assembly_contigs(run_main_pipeline, run_assembly)
            assembly_available = run_assembly.performed
        except RunAssembly.DoesNotExist:
            run_assembly = None
            assembly_available = False

        #
        try:
            run_remap = RunRemapMain.objects.get(
                sample=sample_main, run=run_main_pipeline
            )
            remap_available = run_remap.method != "None"

        except RunRemapMain.DoesNotExist:
            run_remap = EmptyRemapMain
            remap_available = False

        #
        read_classification = ReadClassification.objects.get(
            sample=sample_main, run=run_main_pipeline
        )

        ########
        final_report = FinalReport.objects.filter(
            sample=sample_main, run=run_main_pipeline
        ).order_by("-coverage")
        #
        report_layout_params = TelevirParameters.get_report_layout_params(run_pk=run_pk)
        report_sorter = ReportSorter(final_report, report_layout_params)

        sorted_reports = report_sorter.get_reports()
        excluded_reports_exist = report_sorter.check_excluded_exist()
        empty_reports = report_sorter.get_reports_empty()

        sort_performed = True if report_sorter.analysis_empty is False else False

        if excluded_reports_exist and report_sorter.analysis_empty is False:

            if len(empty_reports.group_list) > 0:
                sorted_reports.append(empty_reports)

        # check has control_flag present
        # has_controlled_flag = False if sample_main.is_control else True
        #########
        clade_heatmap_json = report_sorter.clade_heatmap_json(
            to_keep=[report_group.name for report_group in sorted_reports]
        )

        #########
        private_reads_available = False
        for report_group in sorted_reports:
            if report_group.reports_have_private_reads():
                private_reads_available = True
                break

        contig_classification = ContigClassification.objects.get(
            sample=sample_main, run=run_main_pipeline
        )
        #

        reference_remap_main = ReferenceMap_Main.objects.filter(
            sample=sample_main, run=run_main_pipeline
        )

        context = {
            "project": project_name,
            "run_name": run_name,
            "sort_performed": sort_performed,
            "groups_count": len(sorted_reports),
            "min_shared_reads": round(
                report_layout_params.shared_proportion_threshold * 100, 2
            ),
            "clade_heatmap_json_exists": False if clade_heatmap_json is None else True,
            "clade_heatmap_json": clade_heatmap_json,
            "is_classification": is_classification,
            "remapping_performed": remapping_performed,
            "sample": sample_name,
            "run_main": run_main_pipeline,
            "run_detail": run_detail,
            "qc_report": qc_report,
            "assembly": run_assembly,
            "contig_classification": contig_classification,
            "read_classification": read_classification,
            "run_remap": run_remap,
            "remap_available": remap_available,
            "reference_remap_main": reference_remap_main,
            "number_validated": len(final_report),
            "project_index": project_pk,
            "sample_index": sample_pk,
            "run_index": run_pk,
            "reference_table": raw_reference_table,
            "owner": True,
            "in_control": True,  # has_controlled_flag,
            "report_list": sorted_reports,
            "data_exists": True if not run_main_pipeline.data_deleted else False,
            "excluded_exist": excluded_reports_exist,
            "empty_reports": empty_reports,
            "error_rate_available": report_sorter.error_rate_available,
            "max_error_rate": report_sorter.max_error_rate,
            "quality_avg_available": report_sorter.quality_avg_available,
            "max_quality_avg": report_sorter.max_quality_avg,
            "max_mapped_prop": report_sorter.max_mapped_prop,
            "max_coverage": report_sorter.max_coverage,
            "max_windows_covered": report_sorter.max_windows_covered,
            "overlap_heatmap_available": False,
            "overlap_heatmap_path": report_sorter.overlap_heatmap_path,
            "overlap_pca_exists": report_sorter.overlap_pca_exists,
            "overlap_pca_path": report_sorter.overlap_pca_path,
            "private_reads_available": private_reads_available,
        }

        ### downloadable files
        context["files"] = {}
        # 1. parameters
        params_file_path = run_main_pipeline.params_file_path
        if os.path.exists(params_file_path):
            context["files"]["parameters"] = params_file_path
        # intermediate files zip
        intermediate_reports = run_main_pipeline.intermediate_reports_get()
        run_main_dir = infer_run_media_dir(run_main_pipeline)
        if run_main_dir is not None:
            zip_file_name = "{}_intermediate_reports.zip".format(run_main_pipeline.name)
            file_path = get_create_zip(
                intermediate_reports.files, run_main_dir, zip_file_name
            )
            context["files"]["intermediate_reports_zip"] = file_path

            # final report
            reports_df = run_main_pipeline.get_final_reports_df()
            run_main_dir = infer_run_media_dir(run_main_pipeline)
            reports_df.to_csv(
                os.path.join(run_main_dir, "final_reports.csv"), index=False
            )
            file_path = os.path.join(run_main_dir, "final_reports.csv")
            context["files"]["final_reports_csv"] = file_path

            def eliminate_path_before_media(path: str):
                if PICS.televir_subdirectory in path:
                    path = path.split(PICS.televir_subdirectory)[1]
                    televir_sbdir = os.path.join("/media", PICS.televir_subdirectory)
                    return televir_sbdir + path

                return path.replace(MEDIA_ROOT, "/media")

            for fpath in context["files"]:
                context["files"][fpath] = eliminate_path_before_media(
                    context["files"][fpath]
                )
                context["files"][fpath] = get_link_for_dropdown_item(
                    context["files"][fpath]
                )

        return context


class Sample_ReportCombined(LoginRequiredMixin, generic.CreateView):
    """
    home page
    """

    template_name = "pathogen_identification/sample_detail_compound.html"
    utils = Utils()

    def get_context_data(self, **kwargs):
        project_pk = int(self.kwargs["pk1"])
        sample_pk = int(self.kwargs["pk2"])

        try:
            project_main = Projects.objects.get(pk=project_pk)

        except Exception as e:
            messages.error(self.request, "Project does not exist")
            raise Http404

        try:
            sample = PIProject_Sample.objects.get(pk=sample_pk)
        except Exception as e:
            messages.error(self.request, "Sample does not exist")
            raise Http404

        ####
        project_name = project_main.name
        sample_name = sample.name
        has_controlled_flag = False if sample.is_control else True

        #

        final_report = FinalReport.objects.filter(
            sample=sample, run__project=project_main
        ).order_by("-coverage")

        unique_reports = final_report_best_cov_by_accid(final_report)

        #
        report_layout_params = TelevirParameters.get_report_layout_params(
            project_pk=project_main.pk
        )

        report_sorter = ReportSorter(unique_reports, report_layout_params)
        sort_tree_exists = False
        sort_tree_plot_path = None
        if report_sorter.overlap_manager is not None:
            sort_tree_exists = report_sorter.overlap_manager.tree_plot_exists
            sort_tree_plot_path = report_sorter.overlap_manager.tree_plot_path_render

        sorted_reports = report_sorter.get_reports_compound()
        sort_performed = True if report_sorter.analysis_empty is False else False
        private_reads_available = False
        for report_group in sorted_reports:
            if report_group.reports_have_private_reads():
                private_reads_available = True
                break

        #########
        clade_heatmap_json = report_sorter.clade_heatmap_json(
            to_keep=[report_group.name for report_group in sorted_reports]
        )

        #### graph
        graph_progress = TreeProgressGraph(sample)
        # graph_progress.generate_graph()
        graph_json, graph_id = graph_progress.get_graph_data()
        ####
        runs = set([fr.run.pk for fr in final_report])
        runs_pipeline = RunMain.objects.filter(
            pk__in=runs, run_type=RunMain.RUN_TYPE_PIPELINE
        )
        runs_mapping = RunMain.objects.filter(pk__in=runs).exclude(
            run_type=RunMain.RUN_TYPE_PIPELINE
        )
        runs_number = len(runs)
        runs_exist = runs_number > 0

        context = {
            "project": project_name,
            "graph_json": graph_json,
            "sort_performed": sort_performed,
            "groups_count": len(sorted_reports),
            "min_shared_reads": round(
                report_layout_params.shared_proportion_threshold * 100, 2
            ),
            "clade_heatmap_json_exists": False if clade_heatmap_json is None else True,
            "clade_heatmap_json": clade_heatmap_json,
            "graph_id": graph_id,
            "sample": sample_name,
            "tree_plot_exists": False,
            "tree_plot_path": sort_tree_plot_path,
            "project_index": project_pk,
            "sample_index": sample_pk,
            "report_list": sorted_reports,
            "runs_pipeline": runs_pipeline,
            "runs_mapping": runs_mapping,
            "runs_number": runs_exist,
            "graph_height": runs_number * 22 + 100,
            "owner": True,
            "in_control": has_controlled_flag,
            "error_rate_available": report_sorter.error_rate_available,
            "max_error_rate": report_sorter.max_error_rate,
            "quality_avg_available": report_sorter.quality_avg_available,
            "max_quality_avg": report_sorter.max_quality_avg,
            "max_mapped_prop": report_sorter.max_mapped_prop,
            "max_coverage": report_sorter.max_coverage,
            "max_windows_covered": report_sorter.max_windows_covered,
            "overlap_heatmap_available": False,  # report_sorter.overlap_heatmap_exists,
            "overlap_heatmap_path": report_sorter.overlap_heatmap_path,
            "private_reads_available": private_reads_available,
        }

        return context


class Scaffold_Remap(LoginRequiredMixin, generic.CreateView):
    """
    scaffold remap
    """

    template_name = "pathogen_identification/scaffold_remap.html"

    def get_context_data(self, **kwargs):
        """"""
        # context = super().get_context_data(**kwargs)

        project_pk = int(self.kwargs["pk1"])
        sample_pk = int(self.kwargs["pk2"])
        run_pk = int(self.kwargs["pk3"])
        reference = self.kwargs["reference"]
        user = self.request.user

        try:
            project_main = Projects.objects.get(pk=project_pk)
            sample = PIProject_Sample.objects.get(pk=sample_pk)
            run_main = RunMain.objects.get(pk=run_pk)
        except Exception as e:
            messages.error(self.request, "Project does not exist")
            raise Http404

        if project_main.owner == user:
            run_name = run_main.name
            sample_name = sample.name
            project_name = project_main.name
            reference = (
                reference.replace(".", "_")
                .replace(";", "_")
                .replace(":", "_")
                .replace("|", "_")
            )

            ref_main = ReferenceMap_Main.objects.get(
                reference=reference, sample=sample, run=run_main
            )
            map_db = ReferenceContigs.objects.filter(
                reference=ref_main,
                run=run_main,
            )

        else:
            map_db = ReferenceContigs.objects.none()
            run_name = "run"
            sample_name = "sample"
            project_name = "project"
            reference = "reference"

            messages.error(
                self.request,
                "You do not have permission to access this project.",
                fail_silently=True,
            )

        context = {
            "nav_sample": True,
            "total_items": map_db.count(),
            "show_paginatior": map_db.count() > ConstantsSettings.PAGINATE_NUMBER,
            "show_info_main_page": ShowInfoMainPage(),
            "table": ContigTable(map_db),
            "project": project_name,
            "project_index": project_pk,
            "sample": sample_name,
            "run_name": run_name,
            "reference": reference,
            "sample_index": sample.pk,
            "run_index": run_main.pk,
        }

        return context


def download_file_igv(requestdst):
    """download fasta file"""
    if requestdst.method == "POST":
        form = download_form(requestdst.POST)

        if form.is_valid():
            filepath = form.cleaned_data.get("file_path")

            BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            DLDIR = os.path.join(BASE_DIR, STATICFILES_DIRS[0], "igv_files")
            filepath = os.path.join(DLDIR, filepath)

            if not os.path.exists(filepath):
                return HttpResponseNotFound(f"file {filepath} not found")

            path = open(filepath, "rb")
            # Set the mime type
            mime_type, _ = mimetypes.guess_type(filepath)
            # Set the return value of the HttpResponse
            response = HttpResponse(path, content_type=mime_type)
            # Set the HTTP header for sending to browser
            response["Content-Disposition"] = (
                "attachment; filename=%s" % os.path.basename(filepath)
            )
            # Return the response value
            return response
    else:
        return HttpResponseNotFound("file not found")


def download_file(requestdst):
    """download fasta file"""

    if requestdst.method == "POST":
        form = download_form(requestdst.POST)

        if form.is_valid():
            filepath = form.cleaned_data.get("file_path")

            BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

            if "//" in filepath:
                filepath = "/" + filepath.split("//")[1]

            if not os.path.isfile(filepath):
                return HttpResponseNotFound(f"file {filepath} not found")

            path = open(filepath, "rb")
            # Set the mime type
            mime_type, _ = mimetypes.guess_type(filepath)
            # Set the return value of the HttpResponse
            response = HttpResponse(path, content_type=mime_type)
            # Set the HTTP header for sending to browser
            response["Content-Disposition"] = (
                "attachment; filename=%s" % os.path.basename(filepath)
            )
            # Return the response value
            return response


def download_file_ref(requestdst):
    """download fasta file"""

    if requestdst.method == "POST":
        form = download_ref_form(requestdst.POST)

        if form.is_valid():
            file_request = form.cleaned_data.get("file")
            run_index = int(form.cleaned_data.get("run"))
            taxid = form.cleaned_data.get("taxid")
            accid = form.cleaned_data.get("accid")

            accid_simple = (
                accid.replace(".", "_")
                .replace(";", "_")
                .replace(":", "_")
                .replace("|", "_")
            )

            reference = ReferenceMap_Main.objects.get(
                taxid=taxid,
                run__pk=run_index,
                reference=accid_simple,
            )

            if file_request == "mapped_subset_r1":
                filepath = reference.mapped_subset_r1_fasta
            elif file_request == "mapped_subset_r2":
                filepath = reference.mapped_subset_r2_fasta
            elif file_request == "fasta_file_path":
                filepath = reference.fasta_file_path
            elif file_request == "fasta_index_file_path":
                filepath = reference.fai_file_path
            elif file_request == "bam_file_path":
                filepath = reference.bam_file_path
            elif file_request == "bam_index_file_path":
                filepath = reference.bai_file_path

            else:
                return HttpResponseNotFound(f"file {file_request} not found")

            BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

            # filepath = BASE_DIR + filepath

            if "//" in filepath:
                filepath = "/" + filepath.split("//")[1]

            if not os.path.isfile(filepath):
                return HttpResponseNotFound(f"file {filepath} not found")

            path = open(filepath, "rb")
            # Set the mime type
            mime_type, _ = mimetypes.guess_type(filepath)
            # Set the return value of the HttpResponse
            response = HttpResponse(path, content_type=mime_type)
            # Set the HTTP header for sending to browser
            response["Content-Disposition"] = (
                "attachment; filename=%s" % os.path.basename(filepath)
            )
            # Return the response value
            return response


import zipfile


def generate_zip_file(file_list: list, zip_file_path: str) -> str:
    try:
        with zipfile.ZipFile(zip_file_path, "w") as zip_file:
            for file_path in file_list:
                if os.path.isfile(file_path):
                    zip_file.write(file_path, os.path.basename(file_path))

    except Exception as e:
        print(e)
        return None
    return zip_file_path


def get_create_zip(file_list: list, outdir: str, zip_file_name: str) -> str:
    zip_file_path = os.path.join(outdir, zip_file_name)

    if os.path.exists(zip_file_path):
        return zip_file_path

    zip_file_path = generate_zip_file(file_list, zip_file_path)

    return zip_file_path


def download_intermediate_reports_zipfile(request):
    """
    download intermediate report files in zip"""

    if request.method == "POST":
        run_pk = request.POST.get("run_pk")
        run_main = RunMain.objects.get(pk=int(run_pk))

        intermediate_reports = run_main.intermediate_reports_get()
        run_main_dir = infer_run_media_dir(run_main)
        zip_file_name = "{}_intermediate_reports.zip".format(run_main.name)

        zip_file_path = get_create_zip(
            intermediate_reports.files, run_main_dir, zip_file_name
        )

        path = open(zip_file_path, "rb")
        mime_type, _ = mimetypes.guess_type(zip_file_path)
        response = HttpResponse(path, content_type=mime_type)
        response["Content-Disposition"] = "attachment; filename={}".format(
            zip_file_name
        )

        return response
