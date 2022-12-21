from braces.views import LoginRequiredMixin
from constants.meta_key_and_values import MetaKeyAndValue
from datasets.models import Dataset
from django.contrib import messages
from django.db import transaction
from django.urls import reverse_lazy
from django.views.generic import ListView, TemplateView, UpdateView
from extend_user.models import Profile
from managing_files.manage_database import ManageDatabase
from managing_files.models import Project, ProjectSample, Sample
from pathogen_identification.models import Projects as Televir_Project
from utils.process_SGE import ProcessSGE
from utils.utils import ShowInfoMainPage

from settings.constants_settings import ConstantsSettings
from settings.default_software import DefaultSoftware
from settings.forms import SoftwareForm
from settings.models import Parameter, Software
from settings.tables import SoftwaresTable

# Create your views here.


class index(TemplateView):
    template_name = "settings/index.html"


class Maintenance(TemplateView):
    template_name = "settings/maintenance.html"


class PISettingsView(LoginRequiredMixin, ListView):
    """
    Home page
    """

    #     model = Software
    #     context_object_name = 'software'
    template_name = "settings/settings.html"

    def get_queryset(self):
        """overwrite queryset to not get all software itens available in Software table"""
        return []

    def update_software_params_global_project(self, project):
        """
        update software global to project
        """
        ### get all global software
        query_set = Software.objects.filter(
            owner=self.request.user,
            type_of_use=Software.TYPE_OF_USE_televir_global,
            type_of_software=Software.TYPE_SOFTWARE,
            is_obsolete=False,
        )
        project = Televir_Project.objects.get(pk=project.pk)
        for software in query_set:

            software_parameters = Parameter.objects.filter(
                software=software,
            )

            software.pk = None
            software.type_of_use = Software.TYPE_OF_USE_televir_project

            try:
                Software.objects.get(
                    name=software.name,
                    type_of_use=Software.TYPE_OF_USE_televir_project,
                    parameter__televir_project=project,
                    pipeline_step=software.pipeline_step,
                )

            except Software.MultipleObjectsReturned:
                pass

            except Software.DoesNotExist:
                software.save()
                for parameter in software_parameters:
                    parameter.pk = None
                    parameter.software = software
                    parameter.televir_project = project
                    parameter.save()

    def duplicate_software_params_global_project(self, project):
        """
        duplicate software global to project
        """
        ### get all global software
        query_set = Software.objects.filter(
            owner=self.request.user,
            type_of_use=Software.TYPE_OF_USE_televir_global,
            type_of_software=Software.TYPE_SOFTWARE,
            is_obsolete=False,
        )
        project = Televir_Project.objects.get(pk=project.pk)
        for software in query_set:

            software_parameters = Parameter.objects.filter(
                software=software,
            )

            software.pk = None
            software.type_of_use = Software.TYPE_OF_USE_televir_project

            software.save()

            for parameter in software_parameters:
                parameter.pk = None
                parameter.software = software
                parameter.televir_project = project
                parameter.save()

    def check_project_params_exist(self, project):
        """
        check if project parameters exist
        """

        query_set = Parameter.objects.filter(televir_project=project.pk)
        if query_set.count() == 0:
            return False
        return True

    def get_context_data(self, **kwargs):
        context = super(PISettingsView, self).get_context_data(**kwargs)
        televir_project = None

        if int(self.kwargs["level"]) > 0:
            televir_project = Televir_Project.objects.get(pk=int(self.kwargs["level"]))
        ### test all defaults first, if exist in database
        default_software = DefaultSoftware()
        default_software.test_all_defaults(
            self.request.user
        )  ## the user can have defaults yet

        ### project parameters
        if televir_project:
            if not self.check_project_params_exist(televir_project):
                self.duplicate_software_params_global_project(televir_project)
            else:
                self.update_software_params_global_project(televir_project)

            technologies = [televir_project.technology]

        else:
            technologies = ConstantsSettings.vect_technology

        all_tables = []  ## order by Technology, PipelineStep, table
        ## [ [unique_id, Technology, [ [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], ...],
        ##    [unique_id, Technology, [ [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], ...], etc
        ## Technology goes to NAV-container, PipelineStep goes to NAV-container, then table
        ## Mix parameters with software
        ### IMPORTANT, must have technology__name, because old versions don't
        for technology in technologies:  ## run over all technology
            vect_pipeline_step = []
            for pipeline_step in ConstantsSettings.vect_pipeline_names:
                # print(f"type of use {Software.TYPE_OF_USE_pident}")
                if televir_project is None:
                    query_set = Software.objects.filter(
                        owner=self.request.user,
                        type_of_use=Software.TYPE_OF_USE_televir_global,
                        type_of_software__in=[
                            Software.TYPE_SOFTWARE,
                            Software.TYPE_INSAFLU_PARAMETER,
                        ],
                        technology__name=technology,
                        pipeline_step__name=pipeline_step,
                        is_obsolete=False,
                    )

                else:
                    query_set = Software.objects.filter(
                        owner=self.request.user,
                        type_of_use=Software.TYPE_OF_USE_televir_project,
                        type_of_software__in=[
                            Software.TYPE_SOFTWARE,
                            Software.TYPE_INSAFLU_PARAMETER,
                        ],
                        technology__name=technology,
                        pipeline_step__name=pipeline_step,
                        parameter__televir_project=televir_project,
                        is_obsolete=False,
                    ).distinct()

                ### if there are software
                if query_set.count() > 0:
                    vect_pipeline_step.append(
                        [
                            "{}_{}".format(
                                pipeline_step.replace(" ", "").replace("/", ""),
                                technology.replace(" ", "").replace("/", ""),
                            ),
                            pipeline_step,
                            SoftwaresTable(query_set, televir_project=televir_project),
                        ]
                    )
            ## if there is software for the pipeline step
            if len(vect_pipeline_step) > 0:
                all_tables.append(
                    [
                        technology.replace(" ", "").replace("/", ""),
                        technology,
                        vect_pipeline_step,
                    ]
                )

        context["all_softwares"] = all_tables
        context["nav_settings"] = True
        ## True for global softwares
        if televir_project:
            context["settings_pathid_project"] = True  ## True for project softwares
            context["project_name"] = televir_project.name
            context["project_id"] = televir_project.pk
        else:
            context["settings_pathogenid"] = True
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        return context


class QCSettingsView(LoginRequiredMixin, ListView):
    """
    Home page
    """

    #     model = Software
    #     context_object_name = 'software'
    template_name = "settings/settings.html"

    def get_queryset(self):
        """overwrite queryset to not get all software itens available in Software table"""
        return []

    def get_context_data(self, **kwargs):
        context = super(QCSettingsView, self).get_context_data(**kwargs)

        ### test all defaults first, if exist in database
        default_software = DefaultSoftware()
        default_software.test_all_defaults(
            self.request.user
        )  ## the user can have defaults yet

        all_tables = []  ## order by Technology, PipelineStep, table
        ## [ [unique_id, Technology, [ [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], ...],
        ##    [unique_id, Technology, [ [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], ...], etc
        ## Technology goes to NAV-container, PipelineStep goes to NAV-container, then table
        ## Mix parameters with software
        ### IMPORTANT, must have technology__name, because old versions don't
        for technology in ConstantsSettings.vect_technology:  ## run over all technology
            vect_pipeline_step = []
            for pipeline_step in ConstantsSettings.vect_pipeline_names:
                query_set = Software.objects.filter(
                    owner=self.request.user,
                    type_of_use=Software.TYPE_OF_USE_qc,
                    type_of_software__in=[
                        Software.TYPE_SOFTWARE,
                        Software.TYPE_INSAFLU_PARAMETER,
                    ],
                    technology__name=technology,
                    pipeline_step__name=pipeline_step,
                    is_obsolete=False,
                )

                ### if there are software
                if query_set.count() > 0:
                    vect_pipeline_step.append(
                        [
                            "{}_{}".format(
                                pipeline_step.replace(" ", "").replace("/", ""),
                                technology.replace(" ", "").replace("/", ""),
                            ),
                            pipeline_step,
                            SoftwaresTable(query_set),
                        ]
                    )
            ## if there is software for the pipeline step
            if len(vect_pipeline_step) > 0:
                all_tables.append(
                    [
                        technology.replace(" ", "").replace("/", ""),
                        technology,
                        vect_pipeline_step,
                    ]
                )

        context["all_softwares"] = all_tables
        context["nav_settings"] = True
        context["qc_settings"] = True  ## True for global softwares
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        return context


class SettingsView(LoginRequiredMixin, ListView):
    """
    Home page
    """

    #     model = Software
    #     context_object_name = 'software'
    template_name = "settings/settings.html"

    def get_queryset(self):
        """overwrite queryset to not get all software itens available in Software table"""
        return []

    def get_context_data(self, **kwargs):

        context = super(SettingsView, self).get_context_data(**kwargs)

        ### test all defaults first, if exist in database
        default_software = DefaultSoftware()
        default_software.test_all_defaults(
            self.request.user
        )  ## the user can have defaults yet

        all_tables = []  ## order by Technology, PipelineStep, table
        ## [ [unique_id, Technology, [ [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], ...],
        ##    [unique_id, Technology, [ [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], ...], etc
        ## Technology goes to NAV-container, PipelineStep goes to NAV-container, then table
        ## Mix parameters with software
        ### IMPORTANT, must have technology__name, because old versions don't
        for technology in ConstantsSettings.vect_technology:  ## run over all technology
            vect_pipeline_step = []
            for pipeline_step in ConstantsSettings.vect_pipeline_names:
                query_set = Software.objects.filter(
                    owner=self.request.user,
                    type_of_use=Software.TYPE_OF_USE_global,
                    type_of_software__in=[
                        Software.TYPE_SOFTWARE,
                        Software.TYPE_INSAFLU_PARAMETER,
                    ],
                    technology__name=technology,
                    pipeline_step__name=pipeline_step,
                    is_obsolete=False,
                )

                ### if there are software
                if query_set.count() > 0:
                    vect_pipeline_step.append(
                        [
                            "{}_{}".format(
                                pipeline_step.replace(" ", "").replace("/", ""),
                                technology.replace(" ", "").replace("/", ""),
                            ),
                            pipeline_step,
                            SoftwaresTable(query_set),
                        ]
                    )
            ## if there is software for the pipeline step
            if len(vect_pipeline_step) > 0:
                all_tables.append(
                    [
                        technology.replace(" ", "").replace("/", ""),
                        technology,
                        vect_pipeline_step,
                    ]
                )

        context["all_softwares"] = all_tables
        context["nav_settings"] = True
        context["main_settings"] = True  ## True for global softwares
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        return context


class UpdateParametersView(LoginRequiredMixin, UpdateView):
    model = Software
    form_class = SoftwareForm
    success_url = reverse_lazy("settings-index")
    template_name = "settings/software_update.html"

    ## Other solution to get the reference
    ## https://pypi.python.org/pypi?%3aaction=display&name=django-contrib-requestprovider&version=1.0.1
    def get_form_kwargs(self):
        """
        Set the request to pass in the form
        """
        kw = super(UpdateParametersView, self).get_form_kwargs()
        kw["request"] = self.request  # the trick!
        return kw

    def get_success_url(self):
        """
        get source_pk from update project, need to pass it in context
        """
        context = super(UpdateParametersView, self).get_context_data(**self.kwargs)
        type_of_use = context["software"].type_of_use

        if type_of_use == Software.TYPE_OF_USE_televir_global:
            return reverse_lazy("pathogenID_pipeline", args=(0,))

        return reverse_lazy("settings-index")

    def get_context_data(self, **kwargs):
        context = super(UpdateParametersView, self).get_context_data(**kwargs)

        context["error_cant_see"] = self.request.user != context["software"].owner
        context["type_of_use"] = context["software"].type_of_use
        context["nav_settings"] = True
        context["nav_modal"] = True  ## short the size of modal window
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        return context

    def form_valid(self, form):
        """
        form update
        """

        ## save it...
        with transaction.atomic():
            software = form.save(commit=False)
            paramers = Parameter.objects.filter(software=software)

            b_change = False
            for parameter in paramers:
                if not parameter.can_change:
                    continue
                if parameter.get_unique_id() in form.cleaned_data:
                    value_from_form = "{}".format(
                        form.cleaned_data[parameter.get_unique_id()]
                    )
                    if value_from_form != parameter.parameter:
                        b_change = True
                        parameter.parameter = value_from_form
                        parameter.save()

            if b_change:
                messages.success(
                    self.request,
                    "{} '".format("Software" if software.is_software() else "INSaFLU")
                    + software.name
                    + "' parameters was updated successfully",
                    fail_silently=True,
                )
            else:
                messages.success(
                    self.request,
                    "No parameters to update for {} '".format(
                        "Software" if software.is_software() else "INSaFLU"
                    )
                    + software.name
                    + "' parameters.",
                    fail_silently=True,
                )
        return super(UpdateParametersView, self).form_valid(form)

    ## static method, not need for now.
    form_valid_message = ""  ## need to have this


class UpdateParametersTelevirProjView(LoginRequiredMixin, UpdateView):
    model = Software
    form_class = SoftwareForm
    template_name = "settings/software_update.html"

    ## Other solution to get the reference
    ## https://pypi.python.org/pypi?%3aaction=display&name=django-contrib-requestprovider&version=1.0.1
    def get_form_kwargs(self):
        """
        Set the request to pass in the form
        """
        kw = super(UpdateParametersTelevirProjView, self).get_form_kwargs()
        kw["request"] = self.request  # the trick!
        kw["pk_televir_project"] = self.kwargs.get("pk_televir_project")
        return kw

    def get_success_url(self):
        """
        get source_pk from update project, need to pass it in context
        """
        project_pk = self.kwargs.get("pk_televir_project")
        return reverse_lazy("pathogenID_pipeline", kwargs={"level": project_pk})

    def get_context_data(self, **kwargs):
        context = super(UpdateParametersTelevirProjView, self).get_context_data(
            **kwargs
        )

        context["error_cant_see"] = self.request.user != context["software"].owner
        context["pk_televir_project"] = self.kwargs.get("pk_televir_project")
        context["nav_project"] = True
        context["nav_modal"] = True  ## short the size of modal window
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        return context

    def form_valid(self, form):
        """
        form update
        """
        ## save it...
        with transaction.atomic():
            software = form.save(commit=False)

            project_id = self.kwargs.get("pk_televir_project")
            project = None
            if not project_id is None:
                try:
                    project = Televir_Project.objects.get(pk=project_id)
                except Televir_Project.DoesNotExist:
                    messages.error(
                        self.request,
                        "Software '" + software.name + "' parameters was not updated",
                    )
                    return super(UpdateParametersTelevirProjView, self).form_valid(form)

            paramers = Parameter.objects.filter(
                software=software, televir_project=project
            )
            b_change = False
            for parameter in paramers:
                if not parameter.can_change:
                    continue
                if parameter.get_unique_id() in form.cleaned_data:
                    value_from_form = "{}".format(
                        form.cleaned_data[parameter.get_unique_id()]
                    )
                    if value_from_form != parameter.parameter:
                        b_change = True
                        parameter.parameter = value_from_form
                        parameter.save()

            if b_change:
                messages.success(
                    self.request,
                    "{} '".format("Software" if software.is_software() else "INSaFLU")
                    + software.name
                    + "' parameters were successfully updated for televir project '"
                    + project.name
                    + "'.",
                    fail_silently=True,
                )
            else:
                messages.success(
                    self.request,
                    "No parameters to update for {} '".format(
                        "Software" if software.is_software() else "INSaFLU"
                    )
                    + software.name
                    + "' for televir project '"
                    + project.name
                    + "'.",
                    fail_silently=True,
                )
        return super(UpdateParametersTelevirProjView, self).form_valid(form)

    ## static method, not need for now.
    form_valid_message = ""  ## need to have this


class UpdateParametersProjView(LoginRequiredMixin, UpdateView):
    model = Software
    form_class = SoftwareForm
    template_name = "settings/software_update.html"

    ## Other solution to get the reference
    ## https://pypi.python.org/pypi?%3aaction=display&name=django-contrib-requestprovider&version=1.0.1
    def get_form_kwargs(self):
        """
        Set the request to pass in the form
        """
        kw = super(UpdateParametersProjView, self).get_form_kwargs()
        kw["request"] = self.request  # the trick!
        kw["pk_project"] = self.kwargs.get("pk_proj")
        return kw

    def get_success_url(self):
        """
        get source_pk from update project, need to pass it in context
        """
        project_pk = self.kwargs.get("pk_proj")
        return reverse_lazy("project-settings", kwargs={"pk": project_pk})

    def get_context_data(self, **kwargs):
        context = super(UpdateParametersProjView, self).get_context_data(**kwargs)

        context["error_cant_see"] = self.request.user != context["software"].owner
        context["pk_project"] = self.kwargs.get("pk_proj")
        context["nav_project"] = True
        context["nav_modal"] = True  ## short the size of modal window
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        return context

    def form_valid(self, form):
        """
        form update
        """
        ## save it...
        with transaction.atomic():
            software = form.save(commit=False)

            project_id = self.kwargs.get("pk_proj")
            project = None
            if not project_id is None:
                try:
                    project = Project.objects.get(pk=project_id)
                except Project.DoesNotExist:
                    messages.error(
                        self.request,
                        "Software '" + software.name + "' parameters was not updated",
                    )
                    return super(UpdateParametersProjView, self).form_valid(form)

            paramers = Parameter.objects.filter(software=software, project=project)
            b_change = False
            for parameter in paramers:
                if not parameter.can_change:
                    continue
                if parameter.get_unique_id() in form.cleaned_data:
                    value_from_form = "{}".format(
                        form.cleaned_data[parameter.get_unique_id()]
                    )
                    if value_from_form != parameter.parameter:
                        b_change = True
                        parameter.parameter = value_from_form
                        parameter.save()

            if b_change:
                messages.success(
                    self.request,
                    "{} '".format("Software" if software.is_software() else "INSaFLU")
                    + software.name
                    + "' parameters were successfully updated for project '"
                    + project.name
                    + "'.",
                    fail_silently=True,
                )
            else:
                messages.success(
                    self.request,
                    "No parameters to update for {} '".format(
                        "Software" if software.is_software() else "INSaFLU"
                    )
                    + software.name
                    + "' for project '"
                    + project.name
                    + "'.",
                    fail_silently=True,
                )
        return super(UpdateParametersProjView, self).form_valid(form)

    ## static method, not need for now.
    form_valid_message = ""  ## need to have this


class UpdateParametersDatasetView(LoginRequiredMixin, UpdateView):
    model = Software
    form_class = SoftwareForm
    template_name = "settings/software_update.html"

    ## Other solution to get the reference
    ## https://pypi.python.org/pypi?%3aaction=display&name=django-contrib-requestprovider&version=1.0.1
    def get_form_kwargs(self):
        """
        Set the request to pass in the form
        """
        kw = super(UpdateParametersDatasetView, self).get_form_kwargs()
        kw["request"] = self.request  # the trick!
        kw["pk_dataset"] = self.kwargs.get("pk_dataset")
        return kw

    def get_success_url(self):
        """
        get source_pk from update project, need to pass it in context
        """
        dataset_pk = self.kwargs.get("pk_dataset")
        return reverse_lazy("dataset-settings", kwargs={"pk": dataset_pk})

    def get_context_data(self, **kwargs):
        context = super(UpdateParametersDatasetView, self).get_context_data(**kwargs)

        context["error_cant_see"] = self.request.user != context["software"].owner
        context["pk_dataset"] = self.kwargs.get("pk_dataset")
        context["nav_dataset"] = True
        context["nav_modal"] = True  ## short the size of modal window
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        return context

    def form_valid(self, form):
        """
        form update
        """
        ## save it...
        with transaction.atomic():
            software = form.save(commit=False)

            dataset_id = self.kwargs.get("pk_dataset")
            dataset = None
            if not dataset_id is None:
                try:
                    dataset = Dataset.objects.get(pk=dataset_id)
                except Dataset.DoesNotExist:
                    messages.error(
                        self.request,
                        "Software '" + software.name + "' parameters was not updated",
                    )
                    return super(UpdateParametersDatasetView, self).form_valid(form)

            paramers = Parameter.objects.filter(software=software, dataset=dataset)
            b_change = False
            for parameter in paramers:
                if not parameter.can_change:
                    continue
                if parameter.get_unique_id() in form.cleaned_data:
                    value_from_form = "{}".format(
                        form.cleaned_data[parameter.get_unique_id()]
                    )
                    if value_from_form != parameter.parameter:
                        b_change = True
                        parameter.parameter = value_from_form
                        parameter.save()

            if b_change:
                messages.success(
                    self.request,
                    "{} '".format("Software" if software.is_software() else "INSaFLU")
                    + software.name
                    + "' parameters were successfully updated for dataset '"
                    + dataset.name
                    + "'.",
                    fail_silently=True,
                )
            else:
                messages.success(
                    self.request,
                    "No parameters to update for {} '".format(
                        "Software" if software.is_software() else "INSaFLU"
                    )
                    + software.name
                    + "' for dataset '"
                    + dataset.name
                    + "'.",
                    fail_silently=True,
                )
        return super(UpdateParametersDatasetView, self).form_valid(form)

    ## static method, not need for now.
    form_valid_message = ""  ## need to have this


class UpdateParametersProjSampleView(LoginRequiredMixin, UpdateView):
    model = Software
    form_class = SoftwareForm
    template_name = "settings/software_update.html"

    ## Other solution to get the reference
    ## https://pypi.python.org/pypi?%3aaction=display&name=django-contrib-requestprovider&version=1.0.1
    def get_form_kwargs(self):
        """
        Set the request to pass in the form
        """
        kw = super(UpdateParametersProjSampleView, self).get_form_kwargs()
        kw["request"] = self.request  # the trick!
        kw["pk_project_sample"] = self.kwargs.get("pk_proj_sample")
        return kw

    def get_success_url(self):
        """
        get source_pk from update project, need to pass it in context
        """
        project_sample_pk = self.kwargs.get("pk_proj_sample")
        return reverse_lazy("sample-project-settings", kwargs={"pk": project_sample_pk})

    def get_context_data(self, **kwargs):
        context = super(UpdateParametersProjSampleView, self).get_context_data(**kwargs)

        context["error_cant_see"] = self.request.user != context["software"].owner
        context["nav_project"] = True
        context["pk_proj_sample"] = self.kwargs.get("pk_proj_sample")
        context["sample_project_settings"] = True
        context["nav_modal"] = True  ## short the size of modal window
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        return context

    def form_valid(self, form):
        """
        form update
        """

        ## save it...
        with transaction.atomic():
            software = form.save(commit=False)

            project_sample_id = self.kwargs.get("pk_proj_sample")
            project_sample = None
            if not project_sample_id is None:
                try:
                    project_sample = ProjectSample.objects.get(pk=project_sample_id)
                except ProjectSample.DoesNotExist:
                    messages.error(
                        self.request,
                        "Software '" + software.name + "' parameters was not updated",
                    )
                    return super(UpdateParametersProjSampleView, self).form_valid(form)

            paramers = Parameter.objects.filter(
                software=software, project_sample=project_sample
            )
            b_change_value = False
            for parameter in paramers:
                if not parameter.can_change:
                    continue
                if parameter.get_unique_id() in form.cleaned_data:
                    value_to_change = "{}".format(
                        form.cleaned_data[parameter.get_unique_id()]
                    )
                    if parameter.parameter != value_to_change:
                        b_change_value = True
                        parameter.parameter = value_to_change
                        parameter.save()

        ### re-run this sample
        if b_change_value:
            ### re-run data
            metaKeyAndValue = MetaKeyAndValue()
            manageDatabase = ManageDatabase()
            process_SGE = ProcessSGE()

            ### change flag to nor finished
            project_sample.is_finished = False
            project_sample.save()

            ### get the user
            user = project_sample.project.owner

            ### create a task to perform the analysis of snippy and freebayes
            try:
                (job_name_wait, job_name) = user.profile.get_name_sge_seq(
                    Profile.SGE_PROCESS_projects, Profile.SGE_GLOBAL
                )
                if project_sample.is_sample_illumina():
                    taskID = process_SGE.set_second_stage_snippy(
                        project_sample, user, job_name, [job_name_wait]
                    )
                else:
                    taskID = process_SGE.set_second_stage_medaka(
                        project_sample, user, job_name, [job_name_wait]
                    )

                ### set project sample queue ID
                manageDatabase.set_project_sample_metakey(
                    project_sample,
                    user,
                    metaKeyAndValue.get_meta_key_queue_by_project_sample_id(
                        project_sample.id
                    ),
                    MetaKeyAndValue.META_VALUE_Queue,
                    taskID,
                )

                ### need to collect global files again
                taskID = process_SGE.set_collect_global_files(
                    project_sample.project, user
                )
                manageDatabase.set_project_metakey(
                    project_sample.project,
                    user,
                    metaKeyAndValue.get_meta_key(
                        MetaKeyAndValue.META_KEY_Queue_TaskID_Project,
                        project_sample.project.id,
                    ),
                    MetaKeyAndValue.META_VALUE_Queue,
                    taskID,
                )
            except:
                pass

            messages.success(
                self.request,
                "{} '".format("Software" if software.is_software() else "INSaFLU")
                + software.name
                + "' parameters were successfully updated "
                + "for project '"
                + project_sample.project.name
                + "' and sample '"
                + project_sample.sample.name
                + "'.",
                fail_silently=True,
            )
        else:
            messages.success(
                self.request,
                "No parameters to update for {} '".format(
                    "Software" if software.is_software() else "INSaFLU"
                )
                + software.name
                + " for project '"
                + project_sample.project.name
                + "' and sample '"
                + project_sample.sample.name
                + "'.",
                fail_silently=True,
            )
        return super(UpdateParametersProjSampleView, self).form_valid(form)

    ## static method, not need for now.
    form_valid_message = ""  ## need to have this


class UpdateParametersSampleView(LoginRequiredMixin, UpdateView):
    model = Software
    form_class = SoftwareForm
    template_name = "settings/software_update.html"

    ## Other solution to get the reference
    ## https://pypi.python.org/pypi?%3aaction=display&name=django-contrib-requestprovider&version=1.0.1
    def get_form_kwargs(self):
        """
        Set the request to pass in the form
        """
        kw = super(UpdateParametersSampleView, self).get_form_kwargs()
        kw["request"] = self.request  # the trick!
        kw["pk_sample"] = self.kwargs.get("pk_sample")
        return kw

    def get_success_url(self):
        """
        get source_pk from update sample, need to pass it in context
        """
        sample_pk = self.kwargs.get("pk_sample")
        return reverse_lazy("sample-settings", kwargs={"pk": sample_pk})

    def get_context_data(self, **kwargs):
        context = super(UpdateParametersSampleView, self).get_context_data(**kwargs)

        context["error_cant_see"] = self.request.user != context["software"].owner
        context["nav_sample"] = True
        context["pk_sample"] = self.kwargs.get("pk_sample")
        context["sample_settings"] = True
        context["nav_modal"] = True  ## short the size of modal window
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        return context

    def form_valid(self, form):
        """
        form update
        """
        ## save it...
        with transaction.atomic():
            software = form.save(commit=False)

            sample_id = self.kwargs.get("pk_sample")
            sample = None
            if not sample_id is None:
                try:
                    sample = Sample.objects.get(pk=sample_id)
                except Sample.DoesNotExist:
                    messages.error(
                        self.request,
                        "Software '" + software.name + "' parameters were not updated",
                    )
                    return super(UpdateParametersSampleView, self).form_valid(form)

                ### can not do anything because the sample is running
                if sample.is_sample_in_the_queue:
                    messages.error(
                        self.request,
                        "Sample '{}' is in the queue to process. Software '{}' parameters were not updated".format(
                            sample.name, software.name
                        ),
                    )
                    return super(UpdateParametersSampleView, self).form_valid(form)

            paramers = Parameter.objects.filter(software=software, sample=sample)
            b_change_value = False
            for parameter in paramers:
                if not parameter.can_change:
                    continue
                if parameter.get_unique_id() in form.cleaned_data:
                    value_to_change = "{}".format(
                        form.cleaned_data[parameter.get_unique_id()]
                    )
                    if parameter.parameter != value_to_change:
                        b_change_value = True
                        parameter.parameter = value_to_change
                        parameter.save()

        ### re-run this sample
        if b_change_value:
            ### re-run data
            manageDatabase = ManageDatabase()
            process_SGE = ProcessSGE()

            ### create a task to perform the analysis of NanoFilt
            try:
                (job_name_wait, job_name) = sample.owner.profile.get_name_sge_seq(
                    Profile.SGE_PROCESS_clean_sample, Profile.SGE_SAMPLE
                )
                if sample.is_type_fastq_gz_sequencing():
                    taskID = process_SGE.set_run_trimmomatic_species(
                        sample, sample.owner, job_name
                    )
                else:
                    taskID = process_SGE.set_run_clean_minion(
                        sample, sample.owner, job_name
                    )

                ### set sample queue ID
                manageDatabase.set_sample_metakey(
                    sample,
                    sample.owner,
                    MetaKeyAndValue.META_KEY_Queue_TaskID,
                    MetaKeyAndValue.META_VALUE_Queue,
                    taskID,
                )
            except:
                sample.is_sample_in_the_queue = False
                sample.save()
                pass

            ## refresh sample list for this user
            if not job_name is None:
                process_SGE.set_create_sample_list_by_user(sample.owner, [job_name])

            messages.success(
                self.request,
                "{} '".format("Software" if software.is_software() else "INSaFLU")
                + software.name
                + "' parameters were successfully updated "
                + "for sample '"
                + sample.name
                + "'.",
                fail_silently=True,
            )
        else:
            messages.success(
                self.request,
                "No parameters to update for {} '".format(
                    "Software" if software.is_software() else "INSaFLU"
                )
                + software.name
                + " for sample '"
                + sample.name
                + "'.",
                fail_silently=True,
            )
        return super(UpdateParametersSampleView, self).form_valid(form)

    ## static method, not need for now.
    form_valid_message = ""  ## need to have this
