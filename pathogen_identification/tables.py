from typing import DefaultDict

import django_tables2 as tables
from constants.constants import Constants
from django.conf import settings
from django.urls import reverse
from django.utils.safestring import mark_safe

from pathogen_identification.models import (
    PIProject_Sample,
    Projects,
    ReferenceContigs,
    RunMain,
    SampleQC,
)


class ProjectTable(tables.Table):
    #   Renders a normal value as an internal hyperlink to another page.
    #   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
    description = tables.Column(verbose_name="Description")
    samples = tables.Column("#Samples (P/W/E)", orderable=False, empty_values=())
    last_change_date = tables.Column("Last Change date", empty_values=())
    creation_date = tables.Column("Creation date", empty_values=())

    class Meta:
        model = Projects

        fields = (
            "name",
            "last_change_date",
            "creation_date",
            "samples",
            "results",
        )
        attrs = {"class": "table-striped table-bordered"}
        empty_text = "There are no Projects to show..."

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
    class Meta:
        model = PIProject_Sample
        # attrs = {
        #    "class": "semantic-ui-table",
        # }
        # template_name = "django_tables2/bootstrap4.html"
        attrs = {"class": "paleblue"}
        fields = (
            "name",
            "combinations",
            "report",
            "input",
            "technology",
            "type",
        )

    def render_combinations(self, record):
        return RunMain.objects.filter(
            sample__name=record.name, project__name=record.project
        ).count()

    report = tables.LinkColumn(
        "sample_main", text="Report", args=[tables.A("project__name"), tables.A("name")]
    )


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
    class Meta:
        model = ReferenceContigs
        attrs = {
            "class": "paleblue",
        }
        fields = ("contig", "depth", "depthr", "coverage")


class RunMainTable(tables.Table):
    class Meta:
        model = RunMain
        attrs = {
            "class": "paleblue",
        }
        fields = (
            "name",
            "enrichment",
            "host_depletion",
            "assembly_method",
            "read_classification",
            "contig_classification",
            "finished",
            "runtime",
        )

    report = tables.LinkColumn(
        "sample_detail",
        text="Details",
        args=[tables.A("project"), tables.A("sample"), tables.A("name")],
    )

    def render_runtime(self, record):
        return float(record.runtime.split()[0])
