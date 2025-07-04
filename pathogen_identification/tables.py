import os
from typing import DefaultDict

import django_tables2 as tables
from crequest.middleware import CrequestMiddleware
from django.conf import settings
from django.urls import reverse
from django.utils.safestring import mark_safe
from django.utils.translation import ugettext_lazy as _

from constants.constants import Constants
from managing_files.manage_database import ManageDatabase
from managing_files.models import ProcessControler
from managing_files.models import ProjectSample as InsafluProjectSample
from pathogen_identification.constants_settings import ConstantsSettings as CS
from pathogen_identification.models import (
    ContigClassification,
    FinalReport,
    ParameterSet,
    PIProject_Sample,
    Projects,
    RawReference,
    RawReferenceCompoundModel,
    ReadClassification,
    ReferenceContigs,
    RunAssembly,
    RunMain,
    SampleQC,
    TeleFluProject,
    TelevirRunQC,
)
from pathogen_identification.utilities.reference_utils import (
    check_file_reference_submitted,
    check_reference_exists,
)
from pathogen_identification.utilities.televir_parameters import TelevirParameters
from pathogen_identification.utilities.utilities_general import (
    get_project_dir,
    get_project_dir_no_media_root,
    infer_run_media_dir,
)
from pathogen_identification.utilities.utilities_views import (
    RawReferenceCompound,
    RunMainWrapper,
)
from settings.constants_settings import ConstantsSettings as SettingsCS
from settings.models import Parameter, Software


class ProjectTable(tables.Table):
    #   Renders a normal value as an internal hyperlink to another page.
    #   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
    description = tables.Column(verbose_name="Description", orderable=False)
    settings = tables.Column(empty_values=(), orderable=False)

    samples = tables.Column("#Samples", orderable=False, empty_values=())
    last_change_date = tables.Column("Last Change date", empty_values=())
    creation_date = tables.Column("Creation date", empty_values=())
    results = tables.Column("Project Samples", orderable=False, empty_values=())
    technology = tables.Column(
        verbose_name="Technology", orderable=False, empty_values=()
    )
    finished_processes = tables.Column("Finished", orderable=False, empty_values=())
    running_processes = tables.Column("Running", orderable=False, empty_values=())
    queued_processes = tables.Column("Queued", orderable=False, empty_values=())

    class Meta:
        model = Projects

        fields = (
            "name",
            "results",
            "samples",
            "last_change_date",
            "creation_date",
            "description",
            "technology",
            "running_processes",
        )
        attrs = {"class": "table-striped table-bordered"}
        empty_text = "There are no Projects to show..."

        sequence = (
            "name",
            "results",
            "settings",
            "samples",
            "description",
            "technology",
            "running_processes",
            "queued_processes",
            "finished_processes",
        )

    def render_technology(self, record):
        """
        return a reference name
        """
        return record.technology

    def render_running_processes(self, record):
        """
        return number of running processes in this project"""

        running = 0
        parameter_sets = ParameterSet.objects.filter(
            project=record, sample__sample__is_deleted=False
        )
        for parameter_set in parameter_sets:
            if parameter_set.status == ParameterSet.STATUS_RUNNING:
                running += 1

        return running

    def render_queued_processes(self, record):
        """
        return number of queued processes in this project"""

        queued = 0
        parameter_sets = ParameterSet.objects.filter(
            project=record, sample__sample__is_deleted=False
        )
        for parameter_set in parameter_sets:
            if parameter_set.status == ParameterSet.STATUS_QUEUED:
                queued += 1

        mapping_runs = RunMain.objects.filter(
            project=record,
            run_type=RunMain.RUN_TYPE_MAP_REQUEST,
            status=RunMain.STATUS_PREP,
            parameter_set__status=ParameterSet.STATUS_PROXIED,
        ).count()

        queued += mapping_runs

        return queued

    def render_finished_processes(self, record):
        """
        return number of finished processes in this project"""

        finished = 0
        parameter_sets = ParameterSet.objects.filter(
            project=record, sample__sample__is_deleted=False
        )
        for parameter_set in parameter_sets:
            if parameter_set.status == ParameterSet.STATUS_FINISHED:
                finished += 1

        return finished

    def render_settings(self, record):
        color = ""
        project_settings_exist = Parameter.objects.filter(
            televir_project__pk=record.pk
        ).exists()

        if project_settings_exist:
            color = 'style="color: purple;"'

        parameters = (
            "<a href="
            + reverse("pathogenID_pipeline", kwargs={"level": record.pk})
            + ' data-toggle="tooltip" title="Manage settings">'
            + f'<span ><i class="padding-button-table fa fa-pencil padding-button-table" {color}></i></span></a>'
        )

        if project_settings_exist:
            parameters = parameters + (
                '<a href="#id_reset_modal" id="id_reset_parameters_modal" data-toggle="modal" data-toggle="tooltip" title="Reset"'
                + ' ref_name="'
                + record.name
                + '" pk="'
                + str(record.pk)
                + '"><i class="fa fa-power-off" style="color: orange;" ></i></span> </a>'
            )

        return mark_safe(parameters)

    def render_results(self, record):
        """
        return a reference name
        """
        results = (
            "<a href="
            + reverse("PIproject_samples", args=[record.pk])
            + ' data-toggle="tooltip" title="See Results">'
            + "Samples</a>"
        )

        return mark_safe(results)

    def render_name(self, record):

        current_request = CrequestMiddleware.get_request()
        user = current_request.user

        ## there's nothing to show
        count = ParameterSet.objects.filter(project__id=record.id).count()
        project_sample = record.name

        project_sample = (
            "<a href="
            + reverse("PIproject_samples", args=[record.pk])
            + ' data-toggle="tooltip" title="See Results">'
            + "{}</a>".format(record.name)
        )

        if user.username == Constants.USER_ANONYMOUS:
            return mark_safe(project_sample)
        if user.username == record.owner.username:
            return mark_safe(
                '<a href="#id_remove_modal" id="id_remove_reference_modal" data-toggle="modal" data-toggle="tooltip" title="Delete"'
                + ' ref_name="'
                + record.name
                + '" pk="'
                + str(record.pk)
                + '"><i class="fa fa-trash"></i></span> </a>'
                + project_sample
            )
        return project_sample

    def render_samples(self, record):
        """
        return a reference name
        """
        add_remove = ""

        nsamples = PIProject_Sample.objects.filter(
            project__id=record.id, is_deleted=False
        ).count()

        return mark_safe(
            "{}".format(nsamples)
            + " <a href="
            + reverse("add-sample-PIproject", args=[record.pk])
            + ' data-toggle="tooltip" title="Add samples" ><i class="fa fa-plus-square"></i> Add</a>'  # 		return mark_safe(tip_info + " ({}/{}/{}) ".format(n_processed, n_processing, n_error) + '<a href=# id="id_add_sample_message"' +\
            + add_remove
        )

    def render_creation_date(self, **kwargs):
        record = kwargs.pop("record")
        return record.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)

    def render_last_change_date(self, **kwargs):
        record = kwargs.pop("record")
        if record.last_change_date is None:
            return "Not set yet"
        return record.last_change_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)


class ProjectTableMetagenomics(ProjectTable):
    merge_explify = tables.Column("Actions", orderable=False, empty_values=())

    class Meta:
        model = Projects

        fields = (
            "name",
            "results",
            "samples",
            "merge_explify",
            "last_change_date",
            "creation_date",
            "description",
            "technology",
            "running_processes",
        )
        attrs = {"class": "table-striped table-bordered"}
        empty_text = "There are no Projects to show..."

        sequence = (
            "name",
            "results",
            "merge_explify",
            "settings",
            "samples",
            "description",
            "technology",
            "running_processes",
            "queued_processes",
            "finished_processes",
        )

    def render_merge_explify(self, record: Projects):
        """
        return merge tables modal button
        """

        process_controler = ProcessControler()

        if ProcessControler.objects.filter(
            owner__id=record.owner.pk,
            name=process_controler.get_name_televir_project_merge_explify(
                project_pk=record.pk,
            ),
            is_running=True,
        ).exists():
            return mark_safe(
                '<a href="#" '
                + 'data-toggle="tooltip" '
                + 'title="Running"'
                + '><i class="fa fa-spinner fa-spin"></i></span> </a>'
            )

        eye_blue = '><i class="fa fa-eye"></i></span> </a>'
        eye_purple = '><i class="fa fa-eye" style="color: purple;"></i></span> </a>'
        eye_show = eye_blue

        deploy_explify = (
            "<a "
            + 'href="#id_merge_televir_explify_modal" data-toggle="modal" data-toggle="tooltip" '
            + 'id="merge_explify_modal" '
            + 'title="Merge Explify"'
            + f"project_id={record.pk} "
            + f"ref_name={record.name} "
        )

        project_dir_structure = get_project_dir_no_media_root(record)
        project_dir = get_project_dir(record)
        merge_explify_file = os.path.join(
            project_dir, CS.EXPLIFY_MERGE_SUFFIX + f".{record.pk}.tsv"
        )
        merge_explify_file_provide = os.path.join(
            "/media/",
            project_dir_structure,
            CS.EXPLIFY_MERGE_SUFFIX + f".{record.pk}.tsv",
        )

        found_explify_result = os.path.isfile(merge_explify_file)
        download_button = ""

        if found_explify_result:
            # display icon and download on click

            download_button = (
                '<a rel="nofollow" href="'
                + merge_explify_file_provide
                + '" download="'
                + CS.EXPLIFY_MERGE_SUFFIX
                + f".{record.pk}.tsv"
                + '" '
                + 'data-toggle="tooltip" '
                + 'title="Download"'
                + '><i class="fa fa-download"></i></span> </a>'
            )

        else:
            #
            if ProcessControler.objects.filter(
                owner__id=record.owner.pk,
                name=process_controler.get_name_televir_project_merge_explify(
                    project_pk=record.pk,
                ),
            ).exists():
                eye_show = eye_purple

        deploy_explify = deploy_explify + eye_show + download_button

        return mark_safe(deploy_explify)


class SampleTableOne(tables.Table):
    color_runs = "#f2f2f2"
    cell_attrs = {
        "th": {
            "style": "text-align: center;",
        },
        "td": {"style": "text-align: center;"},
    }
    set_control = tables.Column(
        "Control", orderable=False, empty_values=(), attrs=cell_attrs
    )

    name = tables.Column(
        verbose_name="Sample Name",
        attrs={
            "th": {
                "style": "text-align: left;",
            },
            "td": {
                "style": "text-align: left;",
                "class": "sample-name",
            },
        },
    )
    report = tables.Column(
        verbose_name="Sample Report", orderable=False, empty_values=(), attrs=cell_attrs
    )
    runs = tables.Column(
        verbose_name="Workflows", orderable=False, empty_values=(), attrs=cell_attrs
    )
    deploy = tables.Column(
        verbose_name="Run", orderable=False, empty_values=(), attrs=cell_attrs
    )

    sorting = tables.Column(
        "Sorting", orderable=False, empty_values=(), attrs=cell_attrs
    )

    select_ref = tables.CheckBoxColumn(
        verbose_name="Select Samples",
        accessor="pk",
        orderable=False,
        attrs={
            "th": {
                "style": "background-color: #dce4f0; text-align: center;",
            },
            "td": {"style": "background-color: #dce4f0; text-align: center;"},
            "th__input": {"id": "checkBoxAll"},
        },
    )

    ref_management = tables.Column(
        "References",
        orderable=False,
        empty_values=(),
        attrs={
            "td": {"style": "text-align: center;"},
            "th": {"style": "text-align: center;"},
        },
    )

    combinations = tables.Column(
        verbose_name="Combinations",
        orderable=False,
        empty_values=(),
        attrs={
            "td": {"style": ""},
            "th": {"style": "background-color: #eaf5ff; text-align: center;"},
        },
    )
    mapping_runs = tables.Column(
        "Mapping Runs",
        orderable=False,
        empty_values=(),
        attrs={
            "th": {"style": "background-color: #eaf5ff; text-align: center;"},
        },
    )

    running_processes = tables.Column(
        "Running",
        orderable=False,
        empty_values=(),
        attrs={
            "th": {"style": "background-color: #eaf5ff; text-align: center;"},
        },
    )
    queued_processes = tables.Column(
        "Queued",
        orderable=False,
        empty_values=(),
        attrs={
            "th": {"style": "background-color: #eaf5ff; text-align: center;"},
        },
    )

    class Meta:
        model = PIProject_Sample

        attrs = {
            "class": "paleblue",
        }
        fields = (
            "set_control",
            "name",
            "report",
            "runs",
            "deploy",
            "sorting",
            "ref_management",
            "select_ref",
            "combinations",
            "mapping_runs",
            "running_processes",
            "queued_processes",
        )

    def render_set_control(self, record: PIProject_Sample):
        """
        return a reference name
        """
        is_sample_control = record.is_control

        if is_sample_control:
            return mark_safe(
                '<a href="#id_set_control_modal" id="id_set_control" data-toggle="modal" title="Remove"'
                + ' ref_name="'
                + record.name
                + '" pk="'
                + str(record.pk)
                + '"><i class="fa fa-circle"></i></span> </a>'
            )
        else:
            return mark_safe(
                '<a href="#id_set_control_modal" id="id_set_control" data-toggle="modal" title="Set as control"'
                + ' ref_name="'
                + record.name
                + '" pk="'
                + str(record.pk)
                + '"><i class="fa fa-circle-o"></i></span> </a>'
            )

    def render_report(self, record):
        current_request = CrequestMiddleware.get_request()
        user = current_request.user

        record_name = (
            '<a href="'
            + reverse(
                "televir_sample_compound_report", args=[record.project.pk, record.pk]
            )
            + '">'
            + " <fa class='fa fa-code-fork'></fa>"
            + " Combined Report"
            + "</a>"
        )
        if user.username == Constants.USER_ANONYMOUS:
            return mark_safe("report")
        if user.username == record.project.owner.username:
            return mark_safe(record_name)

    def render_runs(self, record):
        current_request = CrequestMiddleware.get_request()
        user = current_request.user

        record_name = (
            '<a href="'
            + reverse("sample_main", args=[record.project.pk, record.pk])
            + '">'
            + " <fa class='fa fa-reorder'></fa>"
            + " Workflow Panel"
            + "</a>"
        )
        if user.username == Constants.USER_ANONYMOUS:
            return mark_safe("report")
        if user.username == record.project.owner.username:
            return mark_safe(record_name)

    def render_deploy(self, record: PIProject_Sample):

        current_request = CrequestMiddleware.get_request()
        user = current_request.user

        active_runs = ParameterSet.objects.filter(
            sample=record,
            status__in=[ParameterSet.STATUS_RUNNING, ParameterSet.STATUS_QUEUED],
        )

        record_name = '<a><i class="fa fa-bug"></i></span> </a>'

        TELEVIR_DEPLOY_URL = "submit_televir_project_sample"
        if CS.DEPLOYMENT_DEFAULT == CS.DEPLOYMENT_TYPE_PIPELINE:
            TELEVIR_DEPLOY_URL = "submit_televir_runs_project_sample"

        if user.username != record.project.owner.username:
            return mark_safe(record_name)

        record_name = (
            '<a href="#" id="deploypi_sample_btn" class="sample-deploy" data-toggle="modal" data-toggle="tooltip" title="Run Televir Classic Workflow"'
            + ' ref_name="'
            + record.name
            + '"sample_id="'
            + str(record.pk)
            + '" deploy-url="'
            + reverse(
                TELEVIR_DEPLOY_URL,
            )
            + '"'
            + '"><i class="fa fa-flask"></i></span> </a>'
        )

        if (
            ParameterSet.objects.filter(
                sample=record,
                status=ParameterSet.STATUS_FINISHED,
            ).count()
            > 0
            and CS.METAGENOMICS
        ):

            ## encase following butons in a tooltip
            metagen_buttons = " <span class='tooltip-wrap' data-toggle='tooltip' style='display: inline-block; visibility: visible;' >"

            metagen_buttons = metagen_buttons + "</span>"

            record_name += metagen_buttons

        if active_runs.count() > 0:
            record_name += (
                '<a href="#id_kill_modal" id="id_kill_reference_modal" data-toggle="modal" data-toggle="tooltip" title="Cancel"'
                + ' ref_name="'
                + record.name
                + '" pk="'
                + str(record.pk)
                + '"><i class="fa fa-power-off"></i></span> </a>'
            )

        return mark_safe(record_name)

    def render_name(self, record: PIProject_Sample):
        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()
        user = current_request.user

        ### get the link for sample, to expand data
        sample_name = record.sample.name
        sample_name = (
            '<a href="'
            + reverse(
                "televir_sample_compound_report", args=[record.project.pk, record.pk]
            )
            + '">'
            + record.name
            + "</a>"
        )

        # sample_name = "<a>" + record.name + "</a>"

        if user.username == Constants.USER_ANONYMOUS:
            return mark_safe(sample_name)
        if user.username == record.project.owner.username:
            return mark_safe(
                '<a href="#id_remove_modal" id="id_remove_reference_modal" data-toggle="modal"'
                + ' ref_name="'
                + record.name
                + '" pk="'
                + str(record.pk)
                + '"><i class="fa fa-trash"></i></span> </a>'
                + sample_name
            )

        return mark_safe(sample_name)

    def render_sorting(self, record):
        current_request = CrequestMiddleware.get_request()
        user = current_request.user

        if user.username == Constants.USER_ANONYMOUS:
            return mark_safe("report")

        ### check if sorting
        process_controler = ProcessControler()

        process_running = ProcessControler.objects.filter(
            owner__id=user.pk,
            name=process_controler.get_name_televir_project_sample_sort(
                sample_pk=record.pk
            ),
            is_finished=False,
            is_error=False,
        )

        process_finished = ProcessControler.objects.filter(
            owner__id=user.pk,
            name=process_controler.get_name_televir_project_sample_sort(
                sample_pk=record.pk
            ),
            is_finished=True,
            is_error=False,
        )

        if process_running.exists() and not process_finished.exists():
            request_sorting = "<i class='fa fa-spinner fa-spin' title='Sorting'></i>"
            return mark_safe(request_sorting)

        ### check if sorted

        sample_runs = RunMain.objects.filter(sample=record)

        ## return empty square if no report
        if sample_runs.count() == 0:
            return mark_safe('<i class="fa fa-square-o" title="Empty"></i>')
        ## check sorted

        report_layout_params = TelevirParameters.get_report_layout_params(
            project_pk=record.project.pk
        )
        media_dir = None

        for run in sample_runs:
            try:
                media_dir = infer_run_media_dir(run)
                media_dir = os.path.dirname(media_dir)
                break
            except:
                continue

        if media_dir is None:
            return mark_safe('<i class="fa fa-square-o" title="Empty"></i>')

        analysis_df_path = os.path.join(
            media_dir,
            "overlap_analysis_{}.tsv".format(
                report_layout_params.shared_proportion_threshold
            ),
        )

        sorted = os.path.isfile(analysis_df_path)

        ## sorted icon, green if sorted, red if not
        if sorted:
            sorted_icon_assess = (
                ' <i class="fa fa-check" style="color: green;" title="Sorted"></i>'
            )

            return mark_safe(sorted_icon_assess)

        sorted_icon_assess = (
            ' <i class="fa fa-times" style="color: red;" title="un-sorted"></i>'
        )

        request_sorting = (
            ' <a href="#" id="sort_sample_btn" class="sort-sample" data-toggle="modal" data-toggle="tooltip" title="Sort"'
            + ' sample_id="'
            + str(record.pk)
            + '"'
            + ' sort-url="'
            + reverse("sort_sample_reports")
            + '"'
            + '><i class="fa fa-sort"></i></span> </a>'
        )
        return mark_safe(sorted_icon_assess + request_sorting)

    def render_select_ref(self, value, record: PIProject_Sample):
        return mark_safe(
            '<input class="select_sample-checkbox" name="select_ref" id="{}_{}" sample_id={} type="checkbox" value="{}"/>'.format(
                Constants.CHECK_BOX, record.id, record.id, value
            )
        )

    def render_ref_management(self, record: PIProject_Sample):
        nreferences = (
            RawReferenceCompoundModel.objects.filter(sample=record)
            .distinct("accid")
            .count()
        )
        references_management_button = (
            '<a href="'
            + reverse(
                "sample_references_management",
                args=[record.pk],
            )
            + '"'
            + 'title="Sample Reference Management">'
            + '<i class="fa fa-database"></i></span> </a>'
            + f" ({nreferences})"
        )

        return mark_safe(references_management_button)

    report = tables.LinkColumn(
        "sample_main", text="Report", args=[tables.A("project__pk"), tables.A("pk")]
    )

    def render_combinations(self, record: PIProject_Sample):
        return RunMain.objects.filter(
            sample__name=record.name,
            project=record.project,
            parameter_set__status=ParameterSet.STATUS_FINISHED,
            run_type=RunMain.RUN_TYPE_PIPELINE,
        ).count()

    def render_running_processes(self, record: PIProject_Sample):
        """
        return number of running processes in this project"""

        running = ParameterSet.objects.filter(
            sample=record, project=record.project, status=ParameterSet.STATUS_RUNNING
        ).count()

        return running

    def render_queued_processes(self, record: PIProject_Sample):
        """
        return number of running processes in this project"""

        queued = ParameterSet.objects.filter(
            sample=record, project=record.project, status=ParameterSet.STATUS_QUEUED
        ).count()

        mapping_runs = RunMain.objects.filter(
            sample=record,
            run_type__in=[
                RunMain.RUN_TYPE_MAP_REQUEST,
                RunMain.RUN_TYPE_COMBINED_MAPPING,
                RunMain.RUN_TYPE_MAP_REQUEST,
                RunMain.RUN_TYPE_PANEL_MAPPING,
            ],
            parameter_set__status__in=[
                ParameterSet.STATUS_PROXIED,
                ParameterSet.STATUS_QUEUED,
                ParameterSet.STATUS_NOT_STARTED,
            ],
        ).count()

        queued += mapping_runs

        return queued

    def render_mapping_runs(self, record):
        """
        row with two buttons, one for settings, one for deploy
        """

        return RunMain.objects.filter(
            sample__name=record.name,
            project=record.project,
            run_type__in=[
                RunMain.RUN_TYPE_MAP_REQUEST,
                RunMain.RUN_TYPE_COMBINED_MAPPING,
                RunMain.RUN_TYPE_MAP_REQUEST,
                RunMain.RUN_TYPE_PANEL_MAPPING,
            ],
            status__in=[
                RunMain.STATUS_DEFAULT,
                RunMain.STATUS_RUNNING,
                RunMain.STATUS_FINISHED,
            ],
        ).count()


from pathogen_identification.models import ReferenceSourceFile, ReferenceSourceFileMap


class ReferenceSourceFileTable(tables.Table):

    filename = tables.Column(verbose_name="Name", empty_values=())
    description = tables.Column(verbose_name="Description", empty_values=())
    owner = tables.Column(verbose_name="Owner", empty_values=())
    references = tables.Column(verbose_name="References", empty_values=())
    creation_date = tables.Column(verbose_name="Creation Date")

    class Meta:
        attrs = {"class": "paleblue"}

    def render_filename(self, record: ReferenceSourceFile):

        if record.owner is None:
            return record.file

        trash_button = (
            '<a href="#id_remove_modal" class="remove_file" id="id_remove_file_modal" data-toggle="modal" data-toggle="tooltip" title="Delete"'
            + '" pk="'
            + str(record.pk)
            + '" file-name="'
            + str(record.file)
            + '" remove-single-value-url="'
            + reverse("delete_reference_file")
            + '"><i class="fa fa-trash"></i></span> </a>  '
            + record.file
        )

        return mark_safe(f"{trash_button}")

    def render_description(self, record: ReferenceSourceFile):
        return record.description

    def render_owner(self, record: ReferenceSourceFile):
        if record.owner is None:
            return "system"

        return record.owner.username

    def render_references(self, record: ReferenceSourceFile):

        current_refs = ReferenceSourceFileMap.objects.filter(
            reference_source_file=record
        ).count()

        process_controler = ProcessControler()

        upload_processes = ProcessControler.objects.filter(
            owner=record.owner,
            name__icontains=process_controler.get_name_televir_file_upload(record.pk),
            is_finished=False,
        )

        if upload_processes.exists():
            return mark_safe(f"{current_refs} <i class='fa fa-spinner fa-spin'></i>")

        return current_refs

    def render_creation_date(self, record: ReferenceSourceFile):
        if record.creation_date is None:
            return "Not set"

        return record.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)


class TelevirReferencesTable(tables.Table):

    description = tables.Column(verbose_name="Description")
    accid = tables.Column(verbose_name="Accession ID")
    taxid = tables.Column(verbose_name="TaxID")
    source = tables.Column(verbose_name="Files", empty_values=(), orderable=False)
    create_teleflu_reference = tables.Column(
        verbose_name="Create INSaFLU Reference",
        orderable=False,
        empty_values=(),
        attrs={
            "th": {"style": "text-align: center; background-color: #dce4f0;"},
            "td": {"style": "text-align: center;"},
        },
    )

    def __init__(self, records, user_id=None):
        super(TelevirReferencesTable, self).__init__(records)
        self.user_id = user_id

    class Meta:
        # attrs = {"class": "paleblue"}

        fields = (
            "description",
            "accid",
            "taxid",
            "source",
        )

    def render_description(self, record: ReferenceSourceFileMap):
        return record.reference_source.description

    def render_accid(self, record: ReferenceSourceFileMap):
        return record.reference_source.accid

    def render_taxid(self, record: ReferenceSourceFileMap):
        return record.reference_source.taxid

    def render_source(self, record: ReferenceSourceFileMap):
        records_same_accid = ReferenceSourceFileMap.objects.filter(
            reference_source__accid=record.reference_source.accid,
            reference_source_file___owner__in=[
                None,
                record.reference_source_file.owner,
            ],
        ).distinct("reference_source_file")
        files_flat = [
            record.reference_source_file.file for record in records_same_accid
        ]
        return ", ".join(files_flat)

    def render_create_teleflu_reference(self, record: ReferenceSourceFileMap):

        user = CrequestMiddleware.get_request().user

        if check_reference_exists(record.reference_source.accid, user.pk):
            return ""

        if check_file_reference_submitted(ref_id=record.id, user_id=self.user_id):

            return '<i class="fa fa-spinner fa-spin"></i>'

        return mark_safe(
            '<a href="#create_teleflu_reference" '
            + 'id="add_teleflu_reference" '
            + 'class="add_teleflu_reference" '
            + "ref_id="
            + str(record.id)
            + " "
            + 'ref_accid="'
            + record.reference_source.accid
            + '" '
            + 'data-toggle="modal" '
            + 'title="Create INSaFLU Reference" '
            + 'ref-single-value-url="'
            + reverse("create_teleflu_reference")
            + '" '
            + '><i class="fa fa-plus"></i></span> </a>'
        )


class AddedReferenceTable(tables.Table):
    description = tables.Column(verbose_name="Description")
    accid = tables.Column(verbose_name="Accession ID")
    taxid = tables.Column(verbose_name="TaxID")

    class Meta:
        attrs = {"class": "paleblue"}

    def render_description(self, record: RawReference):
        description = record.description
        if description is None:
            description = "Not Available"

        return mark_safe(
            '<a href="#id_remove_modal" id="id_remove_reference_modal" data-toggle="modal" data-toggle="tooltip" title="Delete"'
            + '" pk="'
            + str(record.pk)
            + '" ref_name="'
            + record.accid
            + '"><i class="fa fa-trash"></i></span> </a>'
            + record.description
        )

    def render_accid(self, record: RawReference):
        return record.accid

    def render_taxid(self, record: RawReference):
        return record.taxid


class ReferenceSourceTable(tables.Table):
    select_ref = tables.CheckBoxColumn(
        accessor="pk",
        attrs={
            "th__input": {
                "id": "select-all-checkbox-modal",
                "class": "select-all-checkbox",
                "type": "checkbox",
            }
        },
        orderable=False,
    )
    description = tables.Column(verbose_name="Description")
    accid = tables.Column(verbose_name="Accession id")
    taxid = tables.Column(verbose_name="Taxid")

    class Meta:
        attrs = {"class": "paleblue"}

    def render_select_ref(self, value, record: ReferenceSourceFileMap):
        return mark_safe(
            '<input class="reference-checkbox"  name="select_source_ref" id="{}_{}" ref_id={}  type="checkbox" value="{}"/>'.format(
                Constants.CHECK_BOX,
                record.pk,
                record.pk,
                value,
            )
        )

    def render_description(self, record: ReferenceSourceFileMap):
        return record.description

    def render_accid(self, record: ReferenceSourceFileMap):
        return record.accid

    def render_taxid(self, record: ReferenceSourceFileMap):
        return record.taxid


class TeleFluReferenceTable(tables.Table):
    select_ref = tables.CheckBoxColumn(
        verbose_name=("Select One"),
        accessor="pk",
        orderable=False,
        attrs={"th__input": {"name": "teleflu_select_ref", "style": "display: none;"}},
    )
    description = tables.Column(verbose_name="Description", orderable=False)
    accid = tables.Column(verbose_name="Accession id", orderable=False)
    taxid = tables.Column(verbose_name="Taxid", orderable=False)
    e_rank = tables.Column(
        verbose_name="E Rank",
        orderable=False,
        empty_values=(),
        order_by=("-ensemble_ranking",),
    )

    class Meta:
        attrs = {"class": "paleblue"}

    def render_select_ref(self, value, record: RawReferenceCompound):
        return mark_safe(
            '<input class="teleflu_reference-checkbox"  name="teleflu_select_ref" id="{}_{}" ref_id={}  type="checkbox" value="{}"/>'.format(
                Constants.CHECK_BOX,
                record.selected_mapped_pk,
                record.selected_mapped_pk,
                value,
            )
        )

    def render_description(self, record: RawReferenceCompound):
        return record.description

    def render_accid(self, record: RawReferenceCompound):
        return record.accid

    def render_taxid(self, record: RawReferenceCompound):
        return record.taxid

    def render_e_rank(self, record: RawReferenceCompound):
        if record.ensemble_ranking is None:
            return "N/A"
        return round(record.ensemble_ranking, 3)

    def render_global_ranking(self, record: RawReferenceCompound):
        if record.global_ranking is None:
            return "Not available"
        return record.global_ranking

    def render_ensemble_ranking(self, record: RawReferenceCompound):
        if record.ensemble_ranking is None:
            return "Not available"
        return record.ensemble_ranking


class TeleFluInsaFLuProjectTable(tables.Table):
    header_attrs = {
        "th": {"style": "text-align: center; background-color: #dce4f0;"},
        "td": {"style": "text-align: center;"},
    }
    name = tables.Column(
        verbose_name="Insaflu Project", orderable=False, attrs=header_attrs
    )
    reference = tables.Column(
        verbose_name="Reference", orderable=False, attrs=header_attrs
    )
    results = tables.Column(
        verbose_name="Results", orderable=False, empty_values=(), attrs=header_attrs
    )
    samples = tables.Column(
        verbose_name="Samples", orderable=False, empty_values=(), attrs=header_attrs
    )
    last_change_date = tables.Column(
        "Last Change date", empty_values=(), attrs=header_attrs
    )
    parameters = tables.Column(
        "Parameters", orderable=False, empty_values=(), attrs=header_attrs
    )

    class Meta:
        ## CHANGE TH STYLE IN HEADER ONLY
        attrs = {
            "class": "paleblue",
        }

        sequence = (
            "samples",
            "reference",
            "results",
            "last_change_date",
            "name",
            "parameters",
        )

    def render_name(self, record: TeleFluProject):
        insaflu_project = record.insaflu_project

        if insaflu_project is None:
            return "no associated project"

        ## there's nothing to show
        count = InsafluProjectSample.objects.filter(
            project__id=insaflu_project.id,
            is_deleted=False,
            is_error=False,
            is_finished=True,
        ).count()

        project_name = '<i class="fa fa-eye-slash" aria-hidden="true"></i>'
        if count > 0:
            project_name = (
                "<a href="
                + reverse("show-sample-project-results", args=[insaflu_project.pk])
                + ' data-toggle="tooltip" title="See Results">'
                + "</>"
            )

        project_name = mark_safe(project_name + " " + record.name)

        return project_name

    def render_results(self, record: TeleFluProject):
        insaflu_project = record.insaflu_project
        if insaflu_project is None:
            return "no associated project"

        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()

        ## igv results button
        igv_results_exist = False

        if os.path.exists(record.project_igv_report_media):
            igv_results_exist = True

        igv_results = ""
        if igv_results_exist:
            ## open link to html in new tab
            igv_results = (
                '<a rel="nofollow" target="_blank" href="'
                + reverse("insaflu_project_igv", args=[record.pk])
                + '" title="IGV Report">'
                + '<i class="fa fa-eye"></i></span> </a>'
            )
        else:
            igv_results = '<i class="fa fa-eye-slash" aria-hidden="true"></i>'

        final_display = igv_results

        final_display = mark_safe(final_display)

        return final_display

    def render_description(self, record: TeleFluProject):
        return record.description

    def render_reference(self, record: TeleFluProject):
        return record.raw_reference.description

    def render_samples(self, record: TeleFluProject):
        """
        return a reference name
        """
        insaflu_project = record.insaflu_project
        if insaflu_project is None:
            return "no associated project"
        add_remove = ""

        n_processed = InsafluProjectSample.objects.filter(
            project__id=insaflu_project.id,
            is_deleted=False,
            is_error=False,
            is_finished=True,
        ).count()
        n_error = InsafluProjectSample.objects.filter(
            project__id=insaflu_project.id,
            is_deleted=False,
            is_error=True,
            is_finished=False,
        ).count()
        n_processing = InsafluProjectSample.objects.filter(
            project__id=insaflu_project.id,
            is_deleted=False,
            is_error=False,
            is_finished=False,
        ).count()
        tip_info = '<span ><i class="tip fa fa-info-circle" title="Processed: {}\nWaiting: {}\nError: {}"></i></span>'.format(
            n_processed, n_processing, n_error
        )
        return mark_safe(
            tip_info
            + " ({}/{}/{}) ".format(n_processed, n_processing, n_error)
            + "<a href="
            + reverse("add-sample-project", args=[insaflu_project.pk, 1])
            + ' data-toggle="tooltip" title="Add samples" ><i class="fa fa-plus-square"></i> Add</a>'  #         return mark_safe(tip_info + " ({}/{}/{}) ".format(n_processed, n_processing, n_error) + '<a href=# id="id_add_sample_message"' +\
            + add_remove
        )

    def render_last_change_date(self, record: TeleFluProject):
        return record.last_change_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)

    def render_parameters(self, record: TeleFluProject):
        """
        icon with link to extra info
        """
        insaflu_project = record.insaflu_project
        if insaflu_project is None:
            return "no associated project"
        ## there's nothing to show
        count_not_finished = InsafluProjectSample.objects.filter(
            project__id=insaflu_project.id,
            is_deleted=False,
            is_error=False,
            is_finished=False,
        ).count()

        sz_project_sample = ""

        if count_not_finished > 0:
            sz_project_sample = _("{} processing ".format(count_not_finished))
        else:  ## can change settings
            sz_project_sample += (
                "<a href="
                + reverse("project-settings", args=[insaflu_project.pk])
                + ' data-toggle="tooltip" title="Software settings">'
                + '<span ><i class="padding-button-table fa fa-magic padding-button-table"></i></span></a>'
            )

        return mark_safe(sz_project_sample)


class CompoundReferenceTable(tables.Table):

    select_ref = tables.CheckBoxColumn(
        accessor="pk",
        orderable=False,
        attrs={
            "th": {"style": "text-align: center;"},
            "td": {"style": "text-align: center;"},
        },
    )

    description = tables.Column(verbose_name="Description")

    accid = tables.Column(verbose_name="Accession id")
    taxid = tables.Column(
        verbose_name="Taxid",
        attrs={
            "td": {"style": "text-align: center;"},
            "th": {"style": "text-align: center;"},
        },
    )
    # runs is a integer column that is rendered and is orderable in reverse order
    runs = tables.Column(verbose_name="Runs", order_by="-run_count")
    # mapped column is a link column to the report
    mapped = tables.Column(
        verbose_name="Best Mapping", orderable=False, empty_values=()
    )

    class Meta:
        attrs = {"class": "paleblue"}

    def render_select_ref(self, value, record: RawReferenceCompound):
        return mark_safe(
            '<input class="class-ref-checkbox" ref_id={} name="select_ref" id="{}_{}" type="checkbox" value="{}"/>'.format(
                record.selected_mapped_pk,
                Constants.CHECK_BOX,
                record.selected_mapped_pk,
                value,
            )
        )

    def render_description(self, record: RawReferenceCompound):
        return record.description

    def render_create_teleflu_reference(self, record: RawReferenceCompound):
        return mark_safe(
            '<a href="#" '
            + 'id="add_teleflu_reference" '
            + "ref_id="
            + str(record.selected_mapped_pk)
            + " "
            + 'data-toggle="tooltip" '
            + 'title="Create INSaFLU Reference" '
            + 'url="'
            + reverse("create_teleflu_reference")
            + '" '
            + '><i class="fa fa-plus"></i></span> </a>'
        )

    def render_accid(self, record: RawReferenceCompound):
        return record.accid

    def render_taxid(self, record: RawReferenceCompound):
        return record.taxid

    def render_runs(self, record: RawReferenceCompound):
        return record.runs_str

    def render_mapped(self, record: RawReferenceCompound):
        return record.mapped_html


class CompoundReferenceScore(CompoundReferenceTable):
    # score = tables.Column(
    #    verbose_name="Score",
    #    orderable=True,
    #    empty_values=(),
    #    order_by=("-standard_score",),
    # )
    global_ranking = tables.Column(
        verbose_name="Global Ranking",
        orderable=True,
        empty_values=(),
        order_by=("global_ranking",),
    )
    ensemble_ranking = tables.Column(
        verbose_name="Ensemble Ranking",
        orderable=True,
        empty_values=(),
        order_by=("ensemble_ranking",),
    )

    class Meta:
        attrs = {"class": "paleblue"}
        # order_by = ("-standard_score",)
        sequence = (
            "select_ref",
            "description",
            "accid",
            "taxid",
            "runs",
            "global_ranking",
            "ensemble_ranking",
            "mapped",
        )

    def render_score(self, record: RawReferenceCompound):
        return round(record.standard_score, 3)

    def render_global_ranking(self, record: RawReferenceCompound):
        if record.global_ranking is None:
            return "Not available"
        return record.global_ranking

    def render_ensemble_ranking(self, record: RawReferenceCompound):
        if record.ensemble_ranking is None:
            return "Not available"
        return record.ensemble_ranking


class CompoundRefereceScoreWithScreening(CompoundReferenceScore):
    screenig = tables.Column(
        verbose_name="MM Ranking",
        orderable=True,
        empty_values=(),
        order_by=("-screening_count",),
        attrs={
            "th": {
                "class": "screen-header-style",  # Assign class name for header
            },
            "td": {
                "class": "screen-cell-style",  # Assign class name for cells
            },
        },
    )

    def render_screenig(self, record: RawReferenceCompoundModel):

        return record.screening_count


class RawReferenceTable_Basic(tables.Table):
    taxid = tables.Column(verbose_name="Taxid")
    accid = tables.Column(verbose_name="Accession id")
    description = tables.Column(verbose_name="Description")
    status = tables.Column(
        verbose_name="Status",
    )

    class Meta:
        model = RawReference
        attrs = {"class": "paleblue"}
        fields = (
            "taxid",
            "accid",
            "description",
            "status",
        )

        sequence = (
            "taxid",
            "accid",
            "description",
            "status",
        )

    def render_status(self, record):
        if record.status == RawReference.STATUS_MAPPING:
            return "Running"
        elif record.status == RawReference.STATUS_MAPPED:
            taxids_in_report = FinalReport.objects.filter(run=record.run).values_list(
                "taxid", flat=True
            )
            if record.taxid in taxids_in_report:
                return "Mapped"
            else:
                return "Mapped (0 reads)"

        elif record.status == RawReference.STATUS_FAIL:
            return "Fail"

        elif record.status == RawReference.STATUS_UNMAPPED:
            button = (
                " <a "
                + 'href="#" '
                + 'id="remap_reference" '
                + f"ref_id={record.pk} "
                + f"project_id={record.run.project.pk} "
                + '"><i class="fa fa-eye"></i></span> </a>'
            )
            return mark_safe("Unmapped" + button)


class RawReferenceTable(RawReferenceTable_Basic):

    classification_source = tables.Column(
        verbose_name="Classification Source",
    )
    counts = tables.Column(verbose_name="Counts")

    class Meta:
        model = RawReference
        attrs = {"class": "paleblue"}
        fields = (
            "taxid",
            "accid",
            "description",
            "counts",
            "status",
        )

        sequence = (
            "taxid",
            "accid",
            "description",
            "counts",
            "classification_source",
            "status",
        )

    def render_classification_source(self, record):

        if record.classification_source == "1":
            return "reads"

        if record.classification_source == "2":
            return "contigs"

        if record.classification_source == "3":
            return "reads / contigs"

        else:
            print("Unknown classification source")

            return record.classification_source


class RawReferenceTableNoRemapping(RawReferenceTable):
    class Meta:
        model = RawReference
        attrs = {"class": "paleblue"}
        fields = (
            "taxid",
            "accid",
            "description",
            "counts",
            "status",
        )

        sequence = (
            "taxid",
            "accid",
            "description",
            "counts",
            "classification_source",
            "status",
        )

    def render_status(self, record: RawReference):
        references_management_button = (
            '<a href="'
            + reverse(
                "sample_references_management",
                args=[record.run.sample.pk],
            )
            + '"'
            + 'title="Sample Reference Management">'
            + '<i class="fa fa-database"></i></span> </a>'
        )

        return mark_safe(references_management_button)


class SampleQCTable(tables.Table):
    class Meta:
        model = SampleQC
        attrs = {
            "class": "paleblue",
        }
        fields = (
            "encoding",
            "input_reads",
            "processed_reads",
            "sequence_length",
            "percent_gc",
        )


class ContigTable(tables.Table):
    contig = tables.Column(verbose_name="Contig")
    depth = tables.Column(verbose_name="Depth")
    depthr = tables.Column(verbose_name="Depth Covered")
    coverage = tables.Column(verbose_name="Coverage")

    class Meta:
        model = ReferenceContigs
        fields = ("contig", "depth", "depthr", "coverage")

    def render_contig(self, record):
        return record.contig

    def render_depth(self, record):
        return round(float(record.depth), 3)

    def render_depthr(self, record):
        return round(float(record.depthr), 3)

    def render_coverage(self, record):
        return round(float(record.coverage), 3)


class RunMappingTable(tables.Table):
    name = tables.Column(verbose_name="Run")
    progress_display = tables.Column(verbose_name="Report", orderable=False)
    nmapped = tables.Column(
        verbose_name="Mapped",
        orderable=False,
        empty_values=(),
        attrs={
            "th": {"title": "Success / Mapped / Total References"},
            "td": {"title": "Success / Mapped / Total References"},
        },
    )
    created = tables.Column(verbose_name="Created", orderable=False, empty_values=())
    success = tables.Column(verbose_name="Confirmed", orderable=False, empty_values=())

    extra_filtering = tables.Column(
        verbose_name="Extra filtering",
        orderable=False,
        empty_values=(),
        attrs={
            "td": {"style": "border-left: 1px solid #ddd; text-align: center;"},
            "th": {
                "style": "border-left: 1px solid #ddd; background-color: #dce4f0; text-align: center;"
            },
        },
    )

    enrichment = tables.Column(
        verbose_name="Enrichment",
        orderable=False,
        empty_values=(),
        attrs={
            "td": {"style": "text-align: center;"},
            "th": {"style": "background-color: #dce4f0; text-align: center;"},
        },
    )
    host_depletion = tables.Column(
        verbose_name="Depletion",
        orderable=False,
        empty_values=(),
        attrs={
            "td": {"style": "text-align: center;"},
            "th": {"style": "background-color: #dce4f0; text-align: center;"},
        },
    )

    mapping = tables.Column(
        verbose_name="Mapping",
        orderable=False,
        empty_values=(),
        attrs={
            "td": {"style": "border-right: 1px solid #ddd; text-align: center;"},
            "th": {
                "style": "border-right: 1px solid #ddd;background-color: #dce4f0; text-align: center;"
            },
        },
    )

    runtime = tables.Column(verbose_name="Runtime", orderable=False, empty_values=())

    class Meta:

        attrs = {
            "class": "paleblue",
        }

        sequence = (
            "name",
            "progress_display",
            "extra_filtering",
            "enrichment",
            "host_depletion",
            "mapping",
            "nmapped",
            "created",
            "success",
            "runtime",
        )

    def render_name(self, record: RunMainWrapper):
        prefix = ""

        if record.run_type == RunMain.RUN_TYPE_MAP_REQUEST:
            prefix = "Request - "
        elif record.run_type == RunMain.RUN_TYPE_COMBINED_MAPPING:
            prefix = "Combined - "

        elif record.run_type == RunMain.RUN_TYPE_PANEL_MAPPING:
            prefix = "Panel - "
            if record.panel is not None:
                prefix += record.panel.name + " - "

        return f"{prefix}{record.parameter_set.leaf.index}"

    def render_enrichment(self, record: RunMainWrapper):

        method_name = record.get_pipeline_software(
            SettingsCS.PIPELINE_NAME_viral_enrichment
        )

        return mark_safe(method_name)

    def render_host_depletion(self, record: RunMainWrapper):

        method_name = record.get_pipeline_software(
            SettingsCS.PIPELINE_NAME_host_depletion
        )

        return mark_safe(method_name)

    def render_nmapped(self, record: RunMainWrapper):
        # record = record_wrapped.record
        refs_all = RawReference.objects.filter(run=record.record).count()
        refs_mapped = RawReference.objects.filter(
            run=record.record, status=RawReference.STATUS_MAPPED
        ).count()
        success_mapped = (
            FinalReport.objects.filter(run=record.record).distinct("accid").count()
        )
        string_mapped = f"{success_mapped} / {refs_mapped} / {refs_all}"
        return mark_safe(string_mapped)

    def render_created(self, record: RunMain):
        if record.created_in is None:
            return "N/A"
        return record.created_in.strftime(settings.DATETIME_FORMAT_FOR_TABLE)

    def render_success(self, record: RunMainWrapper):
        success = False
        final_reports = FinalReport.objects.filter(run=record.record).count()
        for final_report in FinalReport.objects.filter(run=record.record):
            if final_report.classification_success != "none":
                success = True
        success = final_reports > 0
        if success:
            return mark_safe('<i class="fa fa-check"></i>')
        else:
            return mark_safe('<i class="fa fa-times"></i>')

    def render_extra_filtering(self, record: RunMainWrapper):

        method_name = record.get_pipeline_software(SettingsCS.PIPELINE_NAME_extra_qc)

        return mark_safe(method_name)

    def render_mapping(self, record: RunMainWrapper):

        method_name = record.get_pipeline_software(
            SettingsCS.PIPELINE_NAME_request_mapping
        )
        return mark_safe(method_name)

    def render_runtime(self, record: RunMainWrapper):

        return mark_safe(record.runtime)


class RunMainTable(tables.Table):

    name = tables.Column(verbose_name="Run")
    report = tables.Column(verbose_name="Report", orderable=False, empty_values=())
    success = tables.Column(verbose_name="Confirmed", orderable=False, empty_values=())

    extra_filtering = tables.Column(
        verbose_name="Extra filtering",
        orderable=False,
        empty_values=(),
        attrs={
            "td": {"style": "border-left: 1px solid #ddd; text-align: center;"},
            "th": {
                "style": "border-left: 1px solid #ddd; background-color: #dce4f0; text-align: center;"
            },
        },
    )

    enrichment = tables.Column(
        verbose_name="Enrichment",
        orderable=False,
        empty_values=(),
        attrs={
            "td": {"style": "text-align: center;"},
            "th": {"style": "background-color: #dce4f0; text-align: center;"},
        },
    )
    host_depletion = tables.Column(
        verbose_name="Depletion",
        orderable=False,
        empty_values=(),
        attrs={
            "td": {"style": "text-align: center;"},
            "th": {"style": "background-color: #dce4f0; text-align: center;"},
        },
    )

    assembly_method = tables.Column(
        verbose_name="Assembly Method",
        orderable=False,
        empty_values=(),
        attrs={
            "td": {"style": "text-align: center;"},
            "th": {"style": "background-color: #dce4f0; text-align: center;"},
        },
    )

    read_classification = tables.Column(
        verbose_name="Read Classification",
        orderable=False,
        empty_values=(),
        attrs={
            "td": {"style": "text-align: center;"},
            "th": {"style": "background-color: #dce4f0; text-align: center;"},
        },
    )

    contig_classification = tables.Column(
        verbose_name="Contig Classification",
        orderable=False,
        empty_values=(),
        attrs={
            "td": {"style": "text-align: center;"},
            "th": {"style": "background-color: #dce4f0; text-align: center;"},
        },
    )

    remapping = tables.Column(
        verbose_name="Remapping",
        orderable=False,
        empty_values=(),
        attrs={
            "td": {"style": "border-right: 1px solid #ddd; text-align: center;"},
            "th": {
                "style": "border-right: 1px solid #ddd; background-color: #dce4f0; text-align: center;"
            },
        },
    )

    runtime = tables.Column(verbose_name="Runtime", orderable=False, empty_values=())

    class Meta:
        # model = RunMain
        attrs = {
            "class": "paleblue",
        }

        sequence = (
            "name",
            "report",
            "extra_filtering",
            "enrichment",
            "host_depletion",
            "assembly_method",
            "contig_classification",
            "read_classification",
            "remapping",
            "success",
            "runtime",
        )

    def render_report(self, record: RunMainWrapper):

        if record.user.username == Constants.USER_ANONYMOUS:
            return mark_safe("report")

        run_log = record.run_progess_tracker()

        return mark_safe(run_log)

    def render_success(self, record: RunMainWrapper):
        success = False
        final_reports = FinalReport.objects.filter(run=record.record).count()
        for final_report in FinalReport.objects.filter(run=record.record):
            if final_report.classification_success != "none":
                success = True
        success = final_reports > 0
        if success:
            return mark_safe('<i class="fa fa-check"></i>')
        else:
            return mark_safe('<i class="fa fa-times"></i>')

    def render_runtime(self, record: RunMainWrapper):

        return mark_safe(record.runtime)

    def render_name(self, record: RunMainWrapper):
        prefix = ""

        if record.run_type == RunMain.RUN_TYPE_MAP_REQUEST:
            prefix = "Request - "
        elif record.run_type == RunMain.RUN_TYPE_COMBINED_MAPPING:
            prefix = "Combined - "

        elif record.run_type == RunMain.RUN_TYPE_PANEL_MAPPING:
            prefix = "Panel - "
            if record.panel is not None:
                prefix += record.panel.name + " - "

        return f"{prefix}{record.parameter_set.leaf.index}"

    def render_enrichment(self, record: RunMainWrapper):

        method_name = record.get_pipeline_software(
            SettingsCS.PIPELINE_NAME_viral_enrichment
        )

        return mark_safe(method_name)

    def render_host_depletion(self, record: RunMainWrapper):

        method_name = record.get_pipeline_software(
            SettingsCS.PIPELINE_NAME_host_depletion
        )

        return mark_safe(method_name)

    def render_extra_filtering(self, record: RunMainWrapper):

        method_name = record.get_pipeline_software(SettingsCS.PIPELINE_NAME_extra_qc)

        return mark_safe(method_name)

    def render_assembly_method(self, record: RunMainWrapper):

        method_name = record.get_pipeline_software(SettingsCS.PIPELINE_NAME_assembly)

        return mark_safe(method_name)

    def render_contig_classification(self, record: RunMainWrapper):

        method_name = record.get_pipeline_software(
            SettingsCS.PIPELINE_NAME_contig_classification
        )

        return mark_safe(method_name)

    def render_read_classification(self, record: RunMainWrapper):

        method_name = record.get_pipeline_software(
            SettingsCS.PIPELINE_NAME_read_classification
        )

        return mark_safe(method_name)

    def render_remapping(self, record: RunMainWrapper):

        method_name = record.get_pipeline_software(SettingsCS.PIPELINE_NAME_remapping)

        return mark_safe(method_name)
