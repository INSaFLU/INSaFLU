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
from pathogen_identification.constants_settings import ConstantsSettings as CS
from pathogen_identification.models import (
    ContigClassification,
    FinalReport,
    ParameterSet,
    PIProject_Sample,
    Projects,
    RawReference,
    ReadClassification,
    ReferenceContigs,
    RunAssembly,
    RunDetail,
    RunMain,
    SampleQC,
    TelevirRunQC,
)
from pathogen_identification.utilities.televir_parameters import TelevirParameters
from pathogen_identification.utilities.utilities_general import (
    get_project_dir,
    get_project_dir_no_media_root,
)
from pathogen_identification.utilities.utilities_views import (
    RawReferenceCompound,
    ReportSorter,
    check_sample_software_exists,
    duplicate_metagenomics_software,
)
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
        from crequest.middleware import CrequestMiddleware

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
    set_control = tables.Column("Control", orderable=False, empty_values=())

    name = tables.Column(verbose_name="Sample Name")
    report = tables.Column(
        verbose_name="Sample Report", orderable=False, empty_values=()
    )
    runs = tables.Column(verbose_name="Workflows", orderable=False, empty_values=())
    deploy = tables.Column(
        verbose_name="Run",
        orderable=False,
        empty_values=(),
        attrs={
            "th": {"style": "text-align: center;"},
            "td": {"style": "text-align: center;"},
        },
    )

    sorting = tables.Column(
        "Sorting",
        orderable=False,
        empty_values=(),
        attrs={
            "th": {
                "style": "text-align: center;",
            },
            "td": {"style": "text-align: center;"},
        },
    )

    select_ref = tables.CheckBoxColumn(
        accessor="pk",
        orderable=False,
        attrs={
            "th": {
                "style": "background-color: #dce4f0; border-left: 5px solid #ddd; text-align: center;",
            },
            "td": {"style": "text-align: center;"},
            "th__input": {"id": "checkBoxAll"},
        },
    )

    ref_management = tables.Column(
        "References",
        orderable=False,
        empty_values=(),
        attrs={
            "td": {"style": "border-left: 5px solid #ddd; text-align: center;"},
            "th": {"style": "background-color: #dce4f0; text-align: center;"},
        },
    )

    combinations = tables.Column(
        verbose_name="Combinations",
        orderable=False,
        empty_values=(),
        attrs={
            "td": {"style": "border-left: 5px solid #ddd;"},
            "th": {"style": "border-left: 5px solid #ddd; background-color: #eaf5ff;"},
        },
    )
    mapping_runs = tables.Column(
        "Mapping Runs",
        orderable=False,
        empty_values=(),
        attrs={
            "th": {"style": "background-color: #eaf5ff;"},
        },
    )

    running_processes = tables.Column(
        "Running",
        orderable=False,
        empty_values=(),
        attrs={
            "th": {"style": "background-color: #eaf5ff;"},
        },
    )
    queued_processes = tables.Column(
        "Queued",
        orderable=False,
        empty_values=(),
        attrs={
            "th": {"style": "background-color: #eaf5ff;"},
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

    def render_set_control(self, record):
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
            + " workflow panel"
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
            '<a href="#" id="deploypi_sample_btn" class="kill-button" data-toggle="modal" data-toggle="tooltip" title="Run Televir"'
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
            # if check_sample_software_exists(record) is False:
            #    duplicate_metagenomics_software(record.project, record)
            ## add light gray background using span

            color = ""

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

        final_report = FinalReport.objects.filter(sample=record).order_by("-coverage")

        ## return empty square if no report
        if final_report.count() == 0:
            return mark_safe('<i class="fa fa-square-o" title="Empty"></i>')
        ## check sorted

        report_layout_params = TelevirParameters.get_report_layout_params(
            project_pk=record.project.pk
        )
        report_sorter = ReportSorter(final_report, report_layout_params)
        sorted = report_sorter.check_analysis_exists()

        ## sorted icon, green if sorted, red if not
        sorted_icon = ""
        if sorted:
            sorted_icon = (
                ' <i class="fa fa-check" style="color: green;" title="Sorted"></i>'
            )
            return mark_safe(sorted_icon)
        else:
            sorted_icon = (
                ' <i class="fa fa-times" style="color: red;" title="un-sorted"></i>'
            )
            request_sorting = (
                ' <a href="#" id="sort_sample_btn" class="kill-button" data-toggle="modal" data-toggle="tooltip" title="Sort"'
                + ' sample_id="'
                + str(record.pk)
                + '"'
                + ' sort-url="'
                + reverse("sort_sample_reports")
                + '"'
                + '><i class="fa fa-sort"></i></span> </a>'
            )
            return mark_safe(sorted_icon + request_sorting)

    def render_select_ref(self, value, record: PIProject_Sample):
        return mark_safe(
            '<input class="select_sample-checkbox" name="select_ref" id="{}_{}" sample_id={} type="checkbox" value="{}"/>'.format(
                Constants.CHECK_BOX, record.id, record.id, value
            )
        )

    def render_ref_management(self, record: PIProject_Sample):
        references_management_button = (
            '<a href="'
            + reverse(
                "sample_references_management",
                args=[record.pk],
            )
            + '"'
            + 'title="Sample Reference Management">'
            + '<i class="fa fa-database"></i></span> </a>'
        )

        return mark_safe(references_management_button)

    report = tables.LinkColumn(
        "sample_main", text="Report", args=[tables.A("project__pk"), tables.A("pk")]
    )

    def render_combinations(self, record: PIProject_Sample):
        return RunMain.objects.filter(
            sample__name=record.name,
            project=record.project,
            parameter_set__status__in=[
                ParameterSet.STATUS_FINISHED,
            ],
            run_type=RunMain.RUN_TYPE_PIPELINE,
        ).count()

    def render_running_processes(self, record: PIProject_Sample):
        """
        return number of running processes in this project"""

        running = 0
        parameter_sets = ParameterSet.objects.filter(
            sample=record, project=record.project
        )
        for parameter_set in parameter_sets:
            if parameter_set.status == ParameterSet.STATUS_RUNNING:
                running += 1

        return running

    def render_queued_processes(self, record: PIProject_Sample):
        """
        return number of running processes in this project"""

        queued = 0
        parameter_sets = ParameterSet.objects.filter(
            sample=record, project=record.project
        )
        for parameter_set in parameter_sets:
            if parameter_set.status == ParameterSet.STATUS_QUEUED:
                queued += 1

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
                RunMain.RUN_TYPE_PANEL_MAPPING,
            ],
            status__in=[
                RunMain.STATUS_DEFAULT,
                RunMain.STATUS_RUNNING,
                RunMain.STATUS_FINISHED,
            ],
        ).count()


class SampleTableTwo(tables.Table):
    sorting = tables.Column("Sorting", orderable=False, empty_values=())
    ref_management = tables.Column("References", orderable=False, empty_values=())

    class Meta:
        model = PIProject_Sample

        attrs = {"class": "paleblue"}
        fields = (
            "sorting",
            "ref_management",
        )


class SampleTableThree(tables.Table):
    class Meta:
        model = PIProject_Sample

        attrs = {"class": "paleblue"}
        fields = (
            "combinations",
            "mapping_runs",
            "running_processes",
            "queued_processes",
        )


class AddedReferenceTable(tables.Table):
    description = tables.Column(verbose_name="Description")
    accid = tables.Column(verbose_name="Accession id")
    taxid = tables.Column(verbose_name="Taxid")

    class Meta:
        attrs = {"class": "paleblue"}

    def render_description(self, record: RawReferenceCompound):
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

    def render_accid(self, record: RawReferenceCompound):
        return record.accid

    def render_taxid(self, record: RawReferenceCompound):
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

    def render_select_ref(self, value, record: RawReferenceCompound):
        return mark_safe(
            '<input class="reference-checkbox"  name="select_source_ref" id="{}_{}" ref_id={}  type="checkbox" value="{}"/>'.format(
                Constants.CHECK_BOX, record.id, record.id, value
            )
        )

    def render_description(self, record: RawReferenceCompound):
        return record.description

    def render_accid(self, record: RawReferenceCompound):
        return record.accid

    def render_taxid(self, record: RawReferenceCompound):
        return record.taxid


class TeleFluReferenceTable(tables.Table):
    select_ref = tables.CheckBoxColumn(
        verbose_name=("Select One"),
        accessor="pk",
        orderable=False,
        attrs={
            "th__input": {
                "name": "teleflu_select_ref",
            }
        },
    )
    description = tables.Column(verbose_name="Description", orderable=False)
    accid = tables.Column(verbose_name="Accession id", orderable=False)
    taxid = tables.Column(verbose_name="Taxid", orderable=False)
    standard_score = tables.Column(
        verbose_name="Max Score",
        orderable=False,
        empty_values=(),
        order_by=("-standard_score",),
    )

    class Meta:
        attrs = {"class": "paleblue"}

    def render_select_ref(self, value, record: RawReferenceCompound):
        return mark_safe(
            '<input class="teleflu_reference-checkbox"  name="teleflu_select_ref" id="{}_{}" ref_id={}  type="checkbox" value="{}"/>'.format(
                Constants.CHECK_BOX, record.id, record.id, value
            )
        )

    def render_description(self, record: RawReferenceCompound):
        return record.description

    def render_accid(self, record: RawReferenceCompound):
        return record.accid

    def render_taxid(self, record: RawReferenceCompound):
        return record.taxid

    def render_standard_score(self, record: RawReferenceCompound):
        return round(record.standard_score, 3)


from managing_files.models import ProjectSample as InsafluProjectSample
from pathogen_identification.models import TeleFluProject


class TeleFluProjectTable(tables.Table):
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
                + record.project_igv_report_media.replace("/insaflu_web/INSaFLU", "")
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
        count = InsafluProjectSample.objects.filter(
            project__id=insaflu_project.id,
            is_deleted=False,
            is_error=False,
            is_finished=True,
        ).count()
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
        accessor="pk", attrs={"th__input": {"id": "checkBoxAll"}}, orderable=False
    )
    description = tables.Column(verbose_name="Description")
    # create_teleflu_reference = tables.Column(
    #    verbose_name="Create Reference",
    #    orderable=False,
    #    empty_values=(),
    #    attrs={
    #        "th": {"style": "text-align: center;"},
    #        "td": {"style": "text-align: center;"},
    #    },
    # )

    accid = tables.Column(verbose_name="Accession id")
    taxid = tables.Column(
        verbose_name="Taxid",
        attrs={
            "td": {"style": "text-align: center;"},
            "th": {"style": "text-align: center;"},
        },
    )
    # runs is a integer column that is rendered and is orderable in reverse order
    runs = tables.Column(verbose_name="Runs", order_by=("-run_count",))
    # mapped column is a link column to the report
    mapped = tables.Column(
        verbose_name="Best Mapping", orderable=False, empty_values=()
    )

    class Meta:
        attrs = {"class": "paleblue"}

    def render_select_ref(self, value, record: RawReferenceCompound):
        return mark_safe(
            '<input class="class-ref-checkbox" ref_id={} name="select_ref" id="{}_{}" type="checkbox" value="{}"/>'.format(
                record.id, Constants.CHECK_BOX, record.id, value
            )
        )

    def render_description(self, record: RawReferenceCompound):
        return record.description

    def render_create_teleflu_reference(self, record: RawReferenceCompound):
        return mark_safe(
            '<a href="#" '
            + 'id="add_teleflu_reference" '
            + "ref_id="
            + str(record.id)
            + " "
            + 'data-toggle="tooltip" '
            + 'title="Create Reference" '
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
    score = tables.Column(
        verbose_name="Score",
        orderable=True,
        empty_values=(),
        order_by=("-standard_score",),
    )

    class Meta:
        attrs = {"class": "paleblue"}
        # order_by = ("-standard_score",)

    def render_score(self, record: RawReferenceCompound):
        return round(record.standard_score, 3)


class CompoundRefereceScoreWithScreening(CompoundReferenceScore):
    screenig = tables.Column(
        verbose_name="MM Ranking",
        orderable=True,
        empty_values=(),
        order_by=("-screening",),
    )

    def render_screenig(self, record: RawReferenceCompound):
        screen_value = 0

        reference = RawReference.objects.filter(pk=record.pk).first()

        screens = RawReference.objects.filter(
            run__sample=reference.run.sample,
            run__run_type=RunMain.RUN_TYPE_SCREENING,
            accid=record.accid,
        )

        if screens.exists():
            screen_value = screens.first().counts

        return screen_value


class RawReferenceTable(tables.Table):
    taxid = tables.Column(verbose_name="Taxid")
    accid = tables.Column(verbose_name="Taxid representativde Accession id")
    description = tables.Column(verbose_name="Taxid representative Description")
    classification_source = tables.Column(verbose_name="Classification Source")
    counts = tables.Column(verbose_name="Counts")
    status = tables.Column(verbose_name="Status")

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
    report = tables.Column(verbose_name="Report", orderable=False, empty_values=())
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
    success = tables.Column(verbose_name="Success", orderable=False, empty_values=())

    extra_filtering = tables.Column(
        verbose_name="Extra filtering", orderable=False, empty_values=()
    )

    enrichment = tables.Column(
        verbose_name="Enrichment", orderable=False, empty_values=()
    )
    host_depletion = tables.Column(
        verbose_name="Depletion", orderable=False, empty_values=()
    )

    remapping = tables.Column(
        verbose_name="Remapping", orderable=False, empty_values=()
    )

    runtime = tables.Column(verbose_name="Runtime", orderable=False, empty_values=())

    class Meta:
        model = RunMain
        attrs = {
            "class": "paleblue",
            "title": "Run Mapping",
            "description": "Run Mapping",
        }

        fields = (
            "name",
            "report",
            "extra_filtering",
            "enrichment",
            "host_depletion",
        )

        sequence = (
            "name",
            "report",
            "nmapped",
            "created",
            "extra_filtering",
            "enrichment",
            "host_depletion",
            "remapping",
            "success",
            "runtime",
        )

        odering = ("created",)

    def get_software_extended_name(self, software_name):
        name_extended = software_name

        try:
            software = Software.objects.filter(name__iexact=software_name).first()
            name_extended = software.name_extended
        except:
            name_extended = software_name

        if name_extended is None:
            return software_name

        # remove parenthesis
        if "(" in name_extended:
            name_extended = name_extended.split("(")[0]

        # if "-" in name_extended:
        #    name_extended = name_extended.split("-")[0]

        return name_extended

    def render_name(self, record: RunMain):
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

    def render_nmapped(self, record: RunMain):
        refs_all = RawReference.objects.filter(run=record).count()
        refs_mapped = RawReference.objects.filter(
            run=record, status=RawReference.STATUS_MAPPED
        ).count()
        success_mapped = (
            FinalReport.objects.filter(run=record).distinct("taxid").count()
        )
        string_mapped = f"{success_mapped} / {refs_mapped} / {refs_all}"
        return mark_safe(string_mapped)

    def render_created(self, record: RunMain):
        if record.created_in is None:
            return "N/A"
        return record.created_in.strftime(settings.DATETIME_FORMAT_FOR_TABLE)

    def render_success(self, record):
        success = False
        final_reports = FinalReport.objects.filter(run=record).count()
        for final_report in FinalReport.objects.filter(run=record):
            if final_report.classification_success != "none":
                success = True
        success = final_reports > 0
        if success:
            return mark_safe('<i class="fa fa-check"></i>')
        else:
            return mark_safe('<i class="fa fa-times"></i>')

    def render_enrichment(self, record: RunMain):
        method = record.enrichment
        method_name = self.get_software_extended_name(method)

        return mark_safe(method_name)

    def render_host_depletion(self, record: RunMain):
        method = record.host_depletion
        method_name = self.get_software_extended_name(method)

        return mark_safe(method_name)

    def render_extra_filtering(self, record):
        try:
            run_qc = TelevirRunQC.objects.get(run=record)
            method = run_qc.method
            method_name = self.get_software_extended_name(method)
            return mark_safe(method_name)
        except:
            return mark_safe("None")

    def render_remapping(self, record: RunMain):
        method = record.remap
        method_name = self.get_software_extended_name(method)

        return mark_safe(method_name)

    def render_report(self, record: RunMain):
        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()
        user = current_request.user

        finished_preprocessing = record.report != "initial"
        finished_assembly = RunAssembly.objects.filter(run=record).count() > 0
        finished_classification = (
            ContigClassification.objects.filter(run=record).exists()
            and ReadClassification.objects.filter(run=record).exists()
        )

        finished_processing = (
            record.parameter_set.status == ParameterSet.STATUS_FINISHED
        )
        # finished_processing = FinalReport.objects.filter(run=record).count() > 0
        finished_remapping = record.report == "finished"
        report_link = (
            '<a href="'
            + reverse(
                "sample_detail",
                args=[record.project.pk, record.sample.pk, record.pk],
            )
            + '">'
            + "<i class='fa fa-bar-chart'></i>"
            + "</a>"
        )

        if finished_processing or finished_remapping:
            if user.username == Constants.USER_ANONYMOUS:
                return mark_safe("report")
            if user.username == record.project.owner.username:
                return mark_safe(report_link)

        else:
            runlog = " <a " + 'href="#" >'
            if finished_preprocessing:
                runlog += '<i class="fa fa-check"'
                runlog += 'title="Preprocessing finished"></i>'
            else:
                runlog += '<i class="fa fa-cog"'
                runlog += 'title="Preprocessing running."></i>'

            runlog += "</a>"

            ###

            runlog += " <a " + 'href="#" >'

            if finished_assembly:
                runlog += '<i class="fa fa-check"'
                runlog += 'title="Assembly finished"></i>'
            else:
                runlog += '<i class="fa fa-cog"'
                if finished_preprocessing:
                    runlog += 'title="Assembly running."></i>'
                else:
                    runlog += 'title="Assembly." style="color: gray;"></i>'
            runlog += "</a>"

            ###

            runlog += " <a " + 'href="#" >'

            if finished_classification:
                runlog += '<i class="fa fa-check"'
                runlog += 'title="Classification finished"></i>'
            else:
                runlog += '<i class="fa fa-cog"'
                if finished_assembly:
                    runlog += 'title="Classification running."></i>'
                else:
                    runlog += 'title="Classification." style="color: gray;"></i>'
            runlog += "</a>"

            runlog += " <a " + 'href="#" >'

            runlog += '<i class="fa fa-cog"'
            if finished_classification:
                runlog += 'title="Mapping to references."></i>'
            else:
                runlog += 'title="Validation mapping" style="color: gray;"></i>'
            runlog += "</a>"

            return mark_safe(runlog)


class RunMainTable(RunMappingTable):
    name = tables.Column(verbose_name="Run")
    report = tables.Column(verbose_name="Report", orderable=False, empty_values=())
    success = tables.Column(verbose_name="Success", orderable=False, empty_values=())

    extra_filtering = tables.Column(
        verbose_name="Extra filtering", orderable=False, empty_values=()
    )

    enrichment = tables.Column(
        verbose_name="Enrichment", orderable=False, empty_values=()
    )
    host_depletion = tables.Column(
        verbose_name="Depletion", orderable=False, empty_values=()
    )

    remapping = tables.Column(
        verbose_name="Remapping", orderable=False, empty_values=()
    )

    runtime = tables.Column(verbose_name="Runtime", orderable=False, empty_values=())

    class Meta:
        model = RunMain
        attrs = {
            "class": "paleblue",
        }

        fields = (
            "name",
            "report",
            "extra_filtering",
            "enrichment",
            "host_depletion",
            "assembly_method",
            "read_classification",
            "contig_classification",
        )

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

    def render_assembly_method(self, record: RunMain):
        method = record.assembly_method
        method_name = self.get_software_extended_name(method)

        return mark_safe(method_name)

    def render_contig_classification(self, record: RunMain):
        method = record.contig_classification
        method_name = self.get_software_extended_name(method)

        return mark_safe(method_name)

    def render_extra_filtering(self, record):
        try:
            run_qc = TelevirRunQC.objects.get(run=record)
            method = run_qc.method
            method_name = self.get_software_extended_name(method)
            return mark_safe(method_name)
        except:
            return mark_safe("None")

    def render_remapping(self, record: RunMain):
        method = record.remap
        method_name = self.get_software_extended_name(method)

        return mark_safe(method_name)
