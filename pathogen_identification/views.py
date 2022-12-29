import logging
import mimetypes
import os

import pandas as pd
from braces.views import FormValidMessageMixin, LoginRequiredMixin
from constants.constants import Constants, FileExtensions, FileType, TypeFile, TypePath
from constants.meta_key_and_values import MetaKeyAndValue
from django import forms
from django.contrib import messages
from django.db import transaction
from django.db.models import Q
from django.forms.models import model_to_dict
from django.http import (
    Http404,
    HttpResponseNotFound,
    HttpResponseRedirect,
    JsonResponse,
)
from django.http.response import HttpResponse
from django.shortcuts import render
from django.template.defaultfilters import filesizeformat, pluralize
from django.urls import reverse_lazy
from django.utils.safestring import mark_safe
from django.views import generic
from django.views.generic import ListView
from django_tables2 import RequestConfig
from extend_user.models import Profile
from fluwebvirus.settings import STATICFILES_DIRS
from managing_files.forms import AddSampleProjectForm
from managing_files.manage_database import ManageDatabase
from managing_files.models import Sample
from managing_files.tables import SampleToProjectsTable
from settings.constants_settings import ConstantsSettings as CS
from settings.default_software_project_sample import DefaultProjectSoftware
from settings.models import Technology
from utils.process_SGE import ProcessSGE
from utils.utils import ShowInfoMainPage, Utils

from pathogen_identification.constants_settings import ConstantsSettings
from pathogen_identification.models import (
    ContigClassification,
    FinalReport,
    PIProject_Sample,
    Projects,
    RawReference,
    ReadClassification,
    ReferenceContigs,
    ReferenceMap_Main,
    RunAssembly,
    RunDetail,
    RunMain,
    RunRemapMain,
    Sample,
)
from pathogen_identification.tables import (
    ContigTable,
    ProjectTable,
    RawReferenceTable,
    RunMainTable,
    SampleTable,
)


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


def is_all_check_box_in_session(vect_check_to_test, request):
    """
    test if all check boxes are in session
    If not remove the ones that are in and create the new ones all False
    """
    utils = Utils()
    dt_data = {}

    ## get the dictonary
    for key in request.session.keys():
        if (
            key.startswith(Constants.CHECK_BOX)
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
            query_set = query_set.filter(
                Q(name__icontains=self.request.GET.get(tag_search))
                | Q(reference__name__icontains=self.request.GET.get(tag_search))
                | Q(
                    project_samples__sample__name__icontains=self.request.GET.get(
                        tag_search
                    )
                )
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

        ###

        ###

        context["project_name"] = project_name
        context["show_paginatior"] = False
        context["nav_project"] = True
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
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
            self.request.session[
                Constants.ERROR_PROJECT_NAME
            ] = "Exists a project with this name."
            self.request.session[Constants.PROJECT_NAME] = name
            b_error = True
        except Projects.DoesNotExist:
            pass
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

        ## catch everything that is not in connection with project
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
            query_set = query_set.filter(
                Q(name__icontains=self.request.GET.get(tag_search))
                | Q(type_subtype__icontains=self.request.GET.get(tag_search))
                | Q(data_set__name__icontains=self.request.GET.get(tag_search))
                | Q(week__icontains=self.request.GET.get(tag_search))
            )
            tag_search
        table = SampleToProjectsTable(query_set)

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
                dt_sample_id_add_temp[
                    sample.id
                ] = 1  ## add the ids that are in the tables
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

        ### check if it show already the settings message
        key_session_name_project_settings = "project_settings_{}".format(project.name)
        if not key_session_name_project_settings in self.request.session:
            self.request.session[key_session_name_project_settings] = True
        else:
            self.request.session[key_session_name_project_settings] = False

        RequestConfig(
            self.request, paginate={"per_page": Constants.PAGINATE_NUMBER}
        ).configure(table)
        if self.request.GET.get(tag_search) != None:
            context[tag_search] = self.request.GET.get(tag_search)
        context["televir_sample"] = True
        context["table"] = table
        context["show_paginatior"] = query_set.count() > Constants.PAGINATE_NUMBER
        context["query_set_count"] = query_set.count()
        context["project_name"] = project.name
        context["nav_modal"] = True  ## short the size of modal window
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
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
            metaKeyAndValue = MetaKeyAndValue()
            manageDatabase = ManageDatabase()
            process_SGE = ProcessSGE()

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

                ### create a task to perform the analysis of snippy and freebayes
                ### Important, it is necessary to run again because can have some changes in the parameters.
                try:
                    if len(job_name_wait) == 0:
                        (
                            job_name_wait,
                            job_name,
                        ) = self.request.user.profile.get_name_sge_seq(
                            Profile.SGE_PROCESS_projects, Profile.SGE_GLOBAL
                        )
                    if sample.is_type_fastq_gz_sequencing():
                        taskID = process_SGE.set_second_stage_snippy(
                            project_sample, self.request.user, job_name, [job_name_wait]
                        )
                    else:
                        taskID = process_SGE.set_second_stage_medaka(
                            project_sample, self.request.user, job_name, [job_name_wait]
                        )

                    ### set project sample queue ID
                    manageDatabase.set_project_sample_metakey(
                        project_sample,
                        self.request.user,
                        metaKeyAndValue.get_meta_key_queue_by_project_sample_id(
                            project_sample.id
                        ),
                        MetaKeyAndValue.META_VALUE_Queue,
                        taskID,
                    )
                except:
                    pass

            ### necessary to calculate the global results again
            if project_sample_add > 0:
                try:
                    taskID = process_SGE.set_collect_global_files(
                        project, self.request.user
                    )
                    manageDatabase.set_project_metakey(
                        project,
                        self.request.user,
                        metaKeyAndValue.get_meta_key(
                            MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id
                        ),
                        MetaKeyAndValue.META_VALUE_Queue,
                        taskID,
                    )
                except:
                    pass

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
    home page
    """

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

        samples = SampleTable(query_set)
        RequestConfig(
            self.request, paginate={"per_page": Constants.PAGINATE_NUMBER}
        ).configure(samples)

        context["table"] = samples
        context["project_index"] = project.pk
        context["project_name"] = project_name
        context["nav_sample"] = True
        context["total_items"] = query_set.count()
        context["show_paginatior"] = query_set.count() > Constants.PAGINATE_NUMBER
        context["show_info_main_page"] = ShowInfoMainPage()
        context["query_set_count"] = query_set.count()
        print(self.request.user.username)
        context["demo"] = True if self.request.user.username == "demo" else False

        return context


class Sample_main(LoginRequiredMixin, generic.CreateView):
    """
    sample main page
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
            )
            sample_name = sample.sample.name
            project_name = project.name

        else:
            messages.error(
                self.request,
                "You do not have permission to access this project.",
                fail_silently=True,
            )
            runs = RunMain.objects.none()
            sample_name = "sample"
            project_name = "project"

        runs_table = RunMainTable(runs)

        RequestConfig(
            self.request, paginate={"per_page": ConstantsSettings.PAGINATE_NUMBER}
        ).configure(runs_table)

        context = {
            "nav_sample": True,
            "total_items": runs.count(),
            "show_paginatior": runs.count() > ConstantsSettings.PAGINATE_NUMBER,
            "show_info_main_page": ShowInfoMainPage(),
            "table": runs_table,
            "sample_name": sample_name,
            "project_main": True,
            "project_name": project_name,
            "project_index": project_pk,
            "sample_index": sample_pk,
            "query_set_count": runs.count(),
        }

        return context


def Project_reports(requesdst, pk1):
    """
    sample main page
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
            run_main = RunMain.objects.get(pk=run_pk)
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
        run_name = run_main.name

        sample_main = run_main.sample
        #

        raw_references = RawReference.objects.filter(run=run_main).order_by("status")

        raw_reference_table = RawReferenceTable(raw_references)

        #
        run_detail = RunDetail.objects.get(sample=sample_main, run=run_main)
        #
        run_assembly = RunAssembly.objects.get(sample=sample_main, run=run_main)
        #
        run_remap = RunRemapMain.objects.get(sample=sample_main, run=run_main)
        #
        read_classification = ReadClassification.objects.get(
            sample=sample_main, run=run_main
        )
        #
        final_report = FinalReport.objects.filter(
            sample=sample_main, run=run_main
        ).order_by("-coverage")
        #
        contig_classification = ContigClassification.objects.get(
            sample=sample_main, run=run_main
        )
        #

        reference_remap_main = ReferenceMap_Main.objects.filter(
            sample=sample_main, run=run_main
        )
        #

        context = {
            "project": project_name,
            "run_name": run_name,
            "sample": sample_name,
            "run_main": run_main,
            "run_detail": run_detail,
            "assembly": run_assembly,
            "contig_classification": contig_classification,
            "read_classification": read_classification,
            "run_remap": run_remap,
            "reference_remap_main": reference_remap_main,
            "final_report": final_report,
            "number_validated": len(final_report),
            "project_index": project_pk,
            "sample_index": sample_pk,
            "run_index": run_pk,
            "reference_table": raw_reference_table,
            "owner": True,
        }

        return context


# def Scaffold_Remap(requesdst, project="", sample="", run="", reference=""):
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
            response[
                "Content-Disposition"
            ] = "attachment; filename=%s" % os.path.basename(filepath)
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
            response[
                "Content-Disposition"
            ] = "attachment; filename=%s" % os.path.basename(filepath)
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
            response[
                "Content-Disposition"
            ] = "attachment; filename=%s" % os.path.basename(filepath)
            # Return the response value
            return response
