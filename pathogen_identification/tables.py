import os
from typing import DefaultDict

import django_tables2 as tables
from constants.constants import Constants
from django.conf import settings
from django.urls import reverse
from django.utils.safestring import mark_safe
from managing_files.manage_database import ManageDatabase
from settings.models import Technology

from pathogen_identification.models import (
    FinalReport,
    ParameterSet,
    PIProject_Sample,
    Projects,
    RawReference,
    ReferenceContigs,
    RunMain,
    SampleQC,
)


class ProjectTable(tables.Table):
    #   Renders a normal value as an internal hyperlink to another page.
    #   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
    description = tables.Column(verbose_name="Description", orderable=False)
    samples = tables.Column("#Samples (P/W/E)", orderable=False, empty_values=())
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

        parameters = (
            "<a href="
            + reverse("pathogenID_pipeline", kwargs={"level": record.pk})
            + ' data-toggle="tooltip" title="Manage settings">'
            + '<span ><i class="padding-button-table fa fa-magic padding-button-table"></i></span></a>'
        )

        return mark_safe(parameters + " " + results)

    def render_name(self, record):
        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()
        user = current_request.user

        ## there's nothing to show
        count = PIProject_Sample.objects.filter(
            project__id=record.id, is_deleted=False, is_error=False, is_finished=True
        ).count()
        project_sample = record.name
        if count > 0:
            project_sample = (
                "<a href="
                + reverse("show-sample-project-results", args=[record.pk])
                + ' data-toggle="tooltip" title="See Results">'
                + "{}</a>".format(record.name)
            )
        else:
            project_sample = record.name

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
        # if (ProjectSample.objects.filter(project__id=record.id, is_deleted=False).count() > 0):
        # 	TODO
        # 	add_remove = ' <a href=' + reverse('remove-sample-project', args=[record.pk]) + '><span ><i class="fa fa-trash"></i></span> Remove</a>'
        # 	add_remove = ' <a href="#"><span ><i class="fa fa-trash"></i></span> Remove</a>'

        n_processed = PIProject_Sample.objects.filter(
            project__id=record.id, is_deleted=False, is_error=False, is_finished=True
        ).count()
        n_error = PIProject_Sample.objects.filter(
            project__id=record.id, is_deleted=False, is_error=True, is_finished=False
        ).count()
        n_processing = PIProject_Sample.objects.filter(
            project__id=record.id, is_deleted=False, is_error=False, is_finished=False
        ).count()
        tip_info = '<span ><i class="tip fa fa-info-circle" title="Processed: {}\nWaiting: {}\nError: {}"></i></span>'.format(
            n_processed, n_processing, n_error
        )
        return mark_safe(
            tip_info
            + " ({}/{}/{}) ".format(n_processed, n_processing, n_error)
            + "<a href="
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


class SampleTable(tables.Table):
    name = tables.Column(verbose_name="Sample Name")

    input = tables.Column(verbose_name="Input", orderable=False, empty_values=())
    combinations = tables.Column(
        verbose_name="Combinations", orderable=False, empty_values=()
    )
    report = tables.Column(verbose_name="Runs", orderable=False, empty_values=())
    running_processes = tables.Column("Running", orderable=False, empty_values=())
    queued_processes = tables.Column("Queued", orderable=False, empty_values=())

    class Meta:
        model = PIProject_Sample

        attrs = {"class": "paleblue"}
        fields = (
            "name",
            "report",
            "input",
            "combinations",
            "running_processes",
        )

    def render_running_processes(self, record):
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

    def render_queued_processes(self, record):
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

    def render_combinations(self, record):
        return RunMain.objects.filter(
            sample__name=record.name, project=record.project
        ).count()

    def render_report(self, record):
        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()
        user = current_request.user

        record_name = (
            '<a href="'
            + reverse("sample_main", args=[record.project.name, record.name])
            + '">'
            + "Run Panel"
            + "</a>"
        )
        if user.username == Constants.USER_ANONYMOUS:
            return mark_safe("report")
        if user.username == record.project.owner.username:
            return mark_safe(record_name)

    def render_name(self, record):
        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()
        user = current_request.user

        ### get the link for sample, to expand data
        sample_name = record.sample.name
        sample_name = (
            '<a href="'
            + reverse("sample_main", args=[record.project.name, record.name])
            + '">'
            + record.name
            + "</a>"
        )

        # sample_name = "<a>" + record.name + "</a>"

        if user.username == Constants.USER_ANONYMOUS:
            return mark_safe(sample_name)
        if (
            user.username == record.project.owner.username
        ):  ## it can't be in any active project
            ## it can have project samples not deleted but in projects deleted
            # for project_samples in record.project_samples.all().filter(
            #    is_deleted=False, is_error=False
            # ):
            #    if (
            #        not project_samples.is_deleted
            #        and not project_samples.project.is_deleted
            #    ):
            #        return mark_safe(sample_name)

            return mark_safe(
                '<a href="#id_remove_modal" id="id_remove_reference_modal" data-toggle="modal"'
                + ' ref_name="'
                + record.name
                + '" pk="'
                + str(record.pk)
                + '"><i class="fa fa-trash"></i></span> </a>'
                + sample_name
            )

        print(sample_name)
        return mark_safe(sample_name)

    report = tables.LinkColumn(
        "sample_main", text="Report", args=[tables.A("project__name"), tables.A("name")]
    )


class RawReferenceTable(tables.Table):
    taxid = tables.Column(verbose_name="Taxid")
    accid = tables.Column(verbose_name="Accid")
    description = tables.Column(verbose_name="Description")
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
                + '"><i class="fa fa-eye"></i></span> </a>'
            )
            return mark_safe("Unmapped" + button)


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


class RunMainTable(tables.Table):

    name = tables.Column(verbose_name="Run Name")
    report = tables.Column(verbose_name="Report", orderable=False, empty_values=())
    success = tables.Column(verbose_name="Success", orderable=False, empty_values=())
    enrichment = tables.Column(
        verbose_name="Enrichment", orderable=False, empty_values=()
    )
    host_depletion = tables.Column(
        verbose_name="Depletion", orderable=False, empty_values=()
    )
    assembly_method = tables.Column(
        verbose_name="Assembly", orderable=False, empty_values=()
    )
    read_classification = tables.Column(
        verbose_name="Read Classification", orderable=False, empty_values=()
    )
    contig_classification = tables.Column(
        verbose_name="Contig Classification", orderable=False, empty_values=()
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
            "enrichment",
            "host_depletion",
            "assembly_method",
            "read_classification",
            "contig_classification",
        )

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

    def render_report(self, record):
        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()
        user = current_request.user

        record_name = (
            '<a href="'
            + reverse(
                "sample_detail",
                args=[record.project.name, record.sample.name, record.name],
            )
            + '">'
            + "<i class='fa fa-bar-chart'></i>"
            + "</a>"
        )
        if user.username == Constants.USER_ANONYMOUS:
            return mark_safe("report")
        if user.username == record.project.owner.username:
            return mark_safe(record_name)

    # report = tables.LinkColumn(
    #    "sample_detail",
    #    text="<i class='fa fa-bar-chart'></i>",
    #    args=[tables.A("project.name"), tables.A("sample.name"), tables.A("name")],
    # )
