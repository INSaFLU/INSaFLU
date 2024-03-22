import codecs
import datetime
import os
from typing import Any, List

import networkx as nx
import numpy as np
import pandas as pd
from django import forms
from django.contrib.auth.models import User
from django.core.validators import RegexValidator
from django.db import models
from django.utils.safestring import mark_safe
from django.utils.translation import gettext_lazy as _

from managing_files.models import Sample
from pathogen_identification.constants_settings import ConstantsSettings as PICS
from pathogen_identification.data_classes import IntermediateFiles

# Create your models here.

# Create your models here.

no_space_validator = RegexValidator(
    r" ",
    _("No spaces allowed"),
    code="invalid_username",
    inverse_match=True,
)


class Projects(models.Model):
    owner = models.ForeignKey(User, on_delete=models.CASCADE)

    name = models.CharField(
        max_length=200,
        db_index=True,
        blank=False,
        verbose_name="Project name",
        default="nameless_project",
        validators=[no_space_validator],
    )
    description = models.TextField(default="", null=True, blank=True)
    technology = models.CharField(
        max_length=200,
        db_index=True,
        blank=True,
        null=True,
        verbose_name="Technology",
    )

    creation_date = models.DateTimeField(
        "uploaded date", db_index=True, auto_now_add=True
    )
    last_change_date = models.DateTimeField("Last change date", blank=True, null=True)
    is_deleted = models.BooleanField(default=False)
    number_passed_sequences = models.SmallIntegerField(
        default=-1
    )  ### has the number of sequences that passed the filters in this project

    ### if is deleted in file system
    is_deleted_in_file_system = models.BooleanField(
        default=False
    )  ## if this file was removed in file system
    date_deleted = models.DateTimeField(
        blank=True, null=True, verbose_name="Date attached"
    )  ## this date has the time of deleted by web page

    results = models.CharField(
        max_length=200, db_index=True, blank=True, null=True, verbose_name="Results"
    )

    running_processes = models.IntegerField(default=0)

    class Meta:
        ordering = ["name", "-creation_date"]

    def __str__(self):
        return self.name + " " + self.description

    def check_delete_schedule(self):
        time_since_last_change = datetime.datetime.now() - self.last_change_date

        if time_since_last_change > datetime.timedelta(weeks=24):
            parametersets = ParameterSet.objects.filter(project=self)
            for parameterset in parametersets:
                parameterset.delete_run_data()

    @property
    def samples_source(self):
        """
        return pk of INSaFLU samples that are part of this project"""
        project_samples = PIProject_Sample.objects.filter(project=self)
        samples = [project_sample.sample.pk for project_sample in project_samples]
        return samples


class SoftwareTree(models.Model):
    """"""

    model = models.IntegerField(default=0)
    version = models.IntegerField(default=0)
    date_created = models.DateTimeField(auto_now_add=True, blank=True, null=True)

    global_index = models.IntegerField(default=0)
    technology = models.CharField(
        max_length=100,
        name="technology",
        blank=True,
        null=True,
    )  # encoding

    owner = models.ForeignKey(User, on_delete=models.CASCADE, blank=True, null=True)
    project = models.ForeignKey(
        Projects, on_delete=models.CASCADE, blank=True, null=True
    )

    class Meta:
        ordering = ["global_index"]

    def get_current_version(self):
        return self.version


class SoftwareTreeNode(models.Model):
    INTERNAL_node = 0
    LEAF_node = 1
    id = models.AutoField(primary_key=True)

    software_tree = models.ForeignKey(SoftwareTree, on_delete=models.CASCADE)
    index = models.SmallIntegerField(default=-1)
    name = models.CharField(
        max_length=200,
        db_index=True,
        blank=True,
        null=True,
        verbose_name="Software name",
    )

    value = models.CharField(
        max_length=200,
        db_index=True,
        blank=True,
        null=True,
    )

    parent = models.ForeignKey(
        "self", on_delete=models.CASCADE, blank=True, null=True, related_name="children"
    )
    node_type = models.CharField(
        max_length=200,
        db_index=True,
        blank=True,
        null=True,
    )
    node_place = models.SmallIntegerField(
        default=INTERNAL_node
    )  ### if it is a software, a parameter or a parameter value

    class Meta:
        ordering = ["name"]

    def get_descendants(self, include_self: bool = True):
        """return all descendants of this node"""

        nodes = SoftwareTreeNode.objects.filter(
            software_tree=self.software_tree
        ).values_list("id", flat=True)
        edges = [
            (node.parent, node.id)
            for node in SoftwareTreeNode.objects.filter(
                software_tree=self.software_tree
            )
        ]

        graph = nx.DiGraph()
        graph.add_nodes_from(nodes)
        graph.add_edges_from(edges)

        descendants = nx.descendants(graph, self.id)

        if include_self:
            descendants = descendants | {self.id}

        return SoftwareTreeNode.objects.filter(id__in=descendants)


class PIProject_Sample(models.Model):
    """
    Main sample information. Connects to the RunMain and QC models.
    """

    sample = models.ForeignKey(Sample, on_delete=models.CASCADE)

    project = models.ForeignKey(
        Projects,
        related_name="project_samples",
        null=True,
        blank=True,
        on_delete=models.CASCADE,
    )

    name = models.CharField(
        max_length=200, db_index=True, blank=True, null=True
    )  # Create your models here. # Name of the sample
    name_extended = models.CharField(
        max_length=200, db_index=True, blank=True, null=True
    )  ## extra name to show in the settings HTML table

    type = models.CharField(
        max_length=10, blank=True, null=True
    )  # sample type: SE or PE
    combinations = models.IntegerField(blank=True, default=0)  # number of combinations
    input = models.TextField(blank=True, null=True)  #
    technology = models.CharField(
        max_length=100,
        name="technology",
        blank=True,
        null=True,
    )  # encoding

    report = models.CharField(max_length=1000, blank=True, null=True)  # report file

    creation_date = models.DateTimeField("uploaded date", auto_now_add=True)
    is_finished = models.BooleanField(default=False)
    is_deleted = models.BooleanField(default=False)
    is_error = models.BooleanField(default=False)  ## if some problem occurs
    is_mask_consensus_sequences = models.BooleanField(
        default=False
    )  ### True if the consensus is masked with SoftwareNames.SOFTWARE_MSA_MASKER
    alert_first_level = models.IntegerField(
        default=0
    )  ## has the number of alerts for high errors
    alert_second_level = models.IntegerField(
        default=0
    )  ## has the number of alerts for low errors
    seq_name_all_consensus = models.CharField(
        blank=True, null=True, max_length=200
    )  ## name of the sample when saved in AllConsensus.fasta file

    ### if is deleted in file system
    is_deleted_in_file_system = models.BooleanField(
        default=False
    )  ## if this file was removed in file system
    date_deleted = models.DateTimeField(
        blank=True, null=True, verbose_name="Date attached"
    )  ## this date has the time of deleted by web page

    running_processes = models.IntegerField(
        default=0
    )  ## has the number of running processes

    is_control = models.BooleanField(
        default=False
    )  ## if this sample is a control sample

    panel_pk_string = models.CharField(
        max_length=150, blank=True, null=True, default=""
    )  ## panel pk string

    class Meta:
        ordering = ["project__id", "-creation_date"]

    @property
    def panels_pks(self) -> List[int]:
        if not self.panel_pk_string:
            return []
        panels = self.panel_pk_string.split(",")
        return [int(panel) for panel in panels]

    @property
    def panels_added(self):
        panels = self.panels_pks
        return ReferencePanel.objects.filter(
            pk__in=panels, panel_type=ReferencePanel.PANEL_TYPE_MAIN
        )

    def remove_panel(self, panel_pk: int):
        panels = self.panels_pks
        panels.remove(panel_pk)
        self.panel_pk_string = ",".join([str(panel) for panel in panels])
        self.save()

    def add_panel(self, panel_pk: int):
        panels = self.panels_pks
        if panel_pk in panels:
            return
        panels.append(panel_pk)
        self.panel_pk_string = ",".join([str(panel) for panel in panels])
        self.save()

    def __str__(self):
        return self.sample.name

    def get_taxid_list(self):
        taxid_list = (
            FinalReport.objects.filter(sample=self)
            .distinct("taxid")
            .values_list("taxid", flat=True)
        )

        return taxid_list

    def get_media_dir(self):
        return os.path.join(
            PICS.media_directory,
            PICS.televir_subdirectory,
            str(self.project.owner.pk),
            str(self.project.pk),
            str(self.sample.pk),
        )


class ParameterSet(models.Model):
    STATUS_NOT_STARTED = 0
    STATUS_RUNNING = 1
    STATUS_FINISHED = 2
    STATUS_ERROR = 3
    STATUS_QUEUED = 4
    STATUS_KILLED = 5
    STATUS_PROXIED = 6
    STATUS_CHOICES = (
        (STATUS_NOT_STARTED, "Not started"),
        (STATUS_RUNNING, "Running"),
        (STATUS_FINISHED, "Finished"),
        (STATUS_ERROR, "Error"),
        (STATUS_QUEUED, "Queued"),
    )

    sample = models.ForeignKey(PIProject_Sample, on_delete=models.PROTECT)
    project = models.ForeignKey(Projects, on_delete=models.CASCADE, null=True)
    leaf = models.ForeignKey(SoftwareTreeNode, on_delete=models.PROTECT, null=True)
    status = models.IntegerField(choices=STATUS_CHOICES, default=STATUS_NOT_STARTED)

    def get_leaf_descendants(self):
        """return software tree nodes that are descendants of the leaf node"""

        if self.leaf is None:
            return []

        descendants = self.leaf.get_descendants(include_self=True)

        return descendants

    def register_subprocess(self):
        self.status = self.STATUS_RUNNING
        self.save()

    def register_finished(self):
        self.status = self.STATUS_FINISHED
        self.save()

    def register_error(self):
        self.status = self.STATUS_ERROR
        self.save()

    def delete_run_data(self):
        if self.status in [self.STATUS_FINISHED, self.STATUS_ERROR]:
            self.status = self.STATUS_NOT_STARTED
            self.save()

            runs = RunMain.objects.filter(parameter_set=self)

            for run in runs:
                if run.status == RunMain.STATUS_FINISHED:
                    continue
                if run.status == RunMain.STATUS_RUNNING:
                    run.status = RunMain.STATUS_KILLED
                    run.save()
                run.delete_data()

    def __str__(self):
        return self.sample.name


class Submitted(models.Model):
    parameter_set = models.ForeignKey(
        ParameterSet, on_delete=models.CASCADE, related_name="submitted_fastq_input"
    )
    date_submitted = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.parameter_set.sample.name + " " + self.parameter_set.leaf.index


class Processed(models.Model):
    parameter_set = models.ForeignKey(ParameterSet, on_delete=models.CASCADE)
    date_processed = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.parameter_set.sample.name + " " + str(self.parameter_set.leaf.index)


class SampleQC(models.Model):
    """Results of sample quality control."""

    ParameterSet = models.ForeignKey(
        ParameterSet, on_delete=models.CASCADE, related_name="sample_qc", default=None
    )

    sample = models.ForeignKey(
        PIProject_Sample, blank=True, null=True, on_delete=models.CASCADE
    )  ## sample
    software = models.CharField(max_length=100, blank=True, null=True)  # software used

    qc_type = models.CharField(
        max_length=100, name="qc_type", blank=True, null=True
    )  # qc type
    encoding = models.CharField(
        max_length=100, name="encoding", blank=True, null=True
    )  # encoding
    input_reads = models.CharField(
        max_length=100, name="input_reads", blank=True, null=True
    )  # input reads

    processed_reads = models.CharField(
        max_length=100, name="processed_reads", blank=True, null=True
    )  # processed reads

    percent_passed = models.FloatField(
        name="percent_passed", blank=True, null=True
    )  # percent passed

    sequence_length = models.CharField(
        max_length=100, blank=True, null=True
    )  # Read length distribution after processing
    percent_gc = models.FloatField(
        name="percent_gc", blank=True, null=True
    )  # percent GC content after filtering.
    # report = tables.LinkColumn("report", text="Report", args=["pk"])
    input_fastqc_report = models.FileField(
        upload_to="input_fastqc_report", blank=True, null=True
    )  # input fastqc report

    processed_fastqc_report = models.FileField(
        upload_to="processed_fastqc_report", blank=True, null=True
    )  # processed fastqc report

    def render_input_fastqc(self):
        html_path = self.input_fastqc_report.path

        if os.path.exists(html_path):
            html_string = codecs.open(html_path, "r").read()
            return mark_safe(html_string)
        return None

    def render_processed_fastqc(self):
        html_path = self.processed_fastqc_report.path
        if os.path.exists(html_path):
            html_string = codecs.open(html_path, "r").read()
            return mark_safe(html_string)
        return None

    class Meta:
        ordering = [
            "sample",
        ]

    def __str__(self):
        return self.sample.name


class QC_REPORT(models.Model):
    RAW = "input"
    PROCESSED = "processed"

    # sample_qc= models.ForeignKey(SampleQC, on_delete=models.CASCADE)

    sample = models.ForeignKey(
        PIProject_Sample, blank=True, null=True, on_delete=models.CASCADE
    )  ## sample

    report_source = models.CharField(max_length=200, blank=True, null=True)  # qc type

    QC_report = models.CharField(max_length=250, blank=True, null=True)  # qc type

    class Meta:
        ordering = [
            "sample",
        ]

    def __str__(self):
        return self.QC_report


class ReferencePanel(models.Model):

    PANEL_TYPE_MAIN = 0
    PANEL_TYPE_COPIED = 1

    panel_type = models.IntegerField(default=PANEL_TYPE_MAIN)
    parent = models.ForeignKey("self", on_delete=models.CASCADE, blank=True, null=True)

    owner = models.ForeignKey(User, on_delete=models.CASCADE)
    name = models.CharField(
        max_length=200,
        db_index=True,
        blank=False,
        verbose_name="Panel name",
        default="nameless_panel",
        validators=[no_space_validator],
    )

    icon = models.CharField(
        max_length=50,
        db_index=True,
        blank=True,
        null=True,
        default="",
        verbose_name="Icon",
    )
    description = models.TextField(default="", null=True, blank=True)
    creation_date = models.DateTimeField(auto_now_add=True, blank=True, null=True)
    is_deleted = models.BooleanField(default=False)
    project_sample = models.ForeignKey(
        PIProject_Sample, on_delete=models.CASCADE, blank=True, null=True
    )

    @property
    def references_count(self):
        return RawReference.objects.filter(panel=self).count()


class RunIndex(models.Model):
    project = models.ForeignKey(
        Projects,
        null=True,
        blank=True,
        on_delete=models.CASCADE,
    )
    sample = models.ForeignKey(Sample, blank=True, null=True, on_delete=models.CASCADE)
    name = models.CharField(
        max_length=100, db_index=True, blank=True, null=True
    )  # Create your models here.


class RunMain(models.Model):
    RUN_TYPE_PIPELINE = 0
    RUN_TYPE_STORAGE = 1
    RUN_TYPE_MAP_REQUEST = 2
    RUN_TYPE_SCREENING = 3
    RUN_TYPE_COMBINED_MAPPING = 4
    RUN_TYPE_PANEL_MAPPING = 5

    STATUS_DEFAULT = 0
    STATUS_PREP = 1
    STATUS_ERROR = 2
    STATUS_RUNNING = 3
    STATUS_FINISHED = 4
    STATUS_KILLED = 5

    run_type = models.IntegerField(default=RUN_TYPE_PIPELINE)
    status = models.IntegerField(default=STATUS_DEFAULT)
    panel = models.ForeignKey(
        ReferencePanel, on_delete=models.CASCADE, blank=True, null=True
    )
    created_in = models.DateTimeField(auto_now_add=True, blank=True, null=True)

    parameter_set = models.ForeignKey(
        ParameterSet, on_delete=models.CASCADE, related_name="run_main", default=None
    )

    project = models.ForeignKey(
        Projects,
        null=True,
        blank=True,
        on_delete=models.CASCADE,
    )

    suprun = models.CharField(max_length=100, blank=True, null=True)
    name = models.CharField(
        max_length=100, db_index=True, blank=True, null=True
    )  # Create your models here.

    sample = models.ForeignKey(
        PIProject_Sample, blank=True, null=True, on_delete=models.CASCADE
    )

    data_deleted = models.BooleanField(default=False)
    last_modified = models.DateTimeField(default=None, null=True, blank=True)

    params_file_path = models.CharField(max_length=250, blank=True, null=True)

    processed_reads_r1 = models.CharField(
        max_length=1000, blank=True, null=True
    )  # processed reads
    processed_reads_r2 = models.CharField(
        max_length=1000, blank=True, null=True
    )  # processed reads

    enrichment = models.CharField(
        max_length=20, blank=True, null=True
    )  # enrichment method if any
    enrichment_performed = models.BooleanField(blank=True)  # enrichment performed
    enrichment_args = models.CharField(
        max_length=150, blank=True, null=True
    )  # enrichment args

    enrichment_db = models.CharField(
        max_length=200, blank=True, null=True
    )  # enrichment db if any

    host_depletion = models.CharField(
        max_length=20, blank=True, null=True
    )  # host depletion method if any
    host_depletion_performed = models.BooleanField(
        blank=True
    )  # host depletion performed

    host_depletion_args = models.CharField(
        max_length=150, blank=True, null=True
    )  # enrichment args

    host_depletion_db = models.CharField(
        max_length=200, blank=True, null=True
    )  # enrichment db if any

    reads_after_processing = models.CharField(
        max_length=100, blank=True, null=True
    )  # reads after processing

    reads_proc_percent = models.CharField(
        max_length=100, blank=True, null=True
    )  # percent of reads after processing

    assembly_performed = models.CharField(
        max_length=10, blank=True, null=True
    )  # assembly method if any
    assembly_method = models.CharField(
        max_length=50, blank=True, null=True
    )  # assembly method if any

    assembly_max = models.CharField(
        max_length=100, blank=True, null=True
    )  # max length of contig.
    read_classification = models.CharField(
        max_length=50, blank=True, null=True
    )  # read classification method if any

    contig_classification = models.CharField(max_length=50, blank=True, null=True)

    remap = models.CharField(
        max_length=50, blank=True, null=True
    )  # remap method if any
    remap_args = models.CharField(max_length=150, blank=True, null=True)

    finished = models.CharField(max_length=10, blank=True, null=True)  # SE or PE
    runtime = models.CharField(max_length=100, blank=True, null=True)

    report = models.CharField(max_length=1000, blank=True, null=True)

    static_dir = models.CharField(max_length=250, blank=True, null=True)

    class Meta:
        ordering = [
            "name",
        ]

    def __str__(self):
        return self.name

    def sorted_reports_get(self):
        return FinalReport.objects.filter(run=self).order_by("-coverage")

    def intermediate_reports_get(self) -> IntermediateFiles:
        contig_classification = ContigClassification.objects.get(run=self)
        read_classification = ReadClassification.objects.get(run=self)
        run_remap = RunRemapMain.objects.get(run=self)

        intermediate_reports = IntermediateFiles(
            read_classification_report=read_classification.read_classification_report,
            contig_classification_report=contig_classification.contig_classification_report,
            remap_main_report=run_remap.merged_log,
            database_matches=run_remap.remap_plan,
        )

        return intermediate_reports

    def get_final_reports_df(self) -> pd.DataFrame:
        final_reports = FinalReport.objects.filter(run=self).exclude(coverage=0)

        columns_to_keep = [
            "taxid",
            "accid",
            "description",
            "ref_db",
            "reference_length",
            "reference_contig_str",
            "coverage",
            "windows_covered",
            "depth",
            "depthR",
            "mapped_reads",
            "ref_proportion",
            "mapped_proportion",
            "ngaps",
            "mapping_success",
            "classification_success",
        ]

        if len(final_reports) == 0:
            return pd.DataFrame(columns=columns_to_keep)

        final_reports_df = pd.DataFrame(list(final_reports.values()))
        final_reports_df = final_reports_df.drop(
            columns=[
                "id",
                "run_id",
            ]
        )

        final_reports_df = final_reports_df[columns_to_keep]
        final_reports_df["ref_db"] = final_reports_df["ref_db"].apply(
            lambda x: x.split("/")[-1]
        )

        return final_reports_df

    def delete_data(self):
        try:
            if os.path.isfile(self.processed_reads_r1):
                os.remove(self.processed_reads_r1)
                self.processed_reads_r1 = None

            if os.path.isfile(self.processed_reads_r2):
                os.remove(self.processed_reads_r2)
                self.processed_reads_r2 = None

            run_assembly = RunAssembly.objects.get(run=self)
            if run_assembly:
                run_assembly.delete_data()

            mapped_references = ReferenceMap_Main.objects.filter(run=self)

            for mapped_reference in mapped_references:
                mapped_reference.delete_data()

            self.data_deleted = True

            self.save()

        except Exception as e:
            print(e)


class RunReadsRegister(models.Model):
    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)

    qc_reads_r1 = models.CharField(max_length=1000, blank=True, null=True)  # qc reads
    qc_reads_r2 = models.CharField(max_length=1000, blank=True, null=True)  # qc reads

    enriched_reads_r1 = models.CharField(
        max_length=1000, blank=True, null=True
    )  # enriched reads
    enriched_reads_r2 = models.CharField(
        max_length=1000, blank=True, null=True
    )  # enriched reads

    depleted_reads_r1 = models.CharField(
        max_length=1000, blank=True, null=True
    )  # depleted reads
    depleted_reads_r2 = models.CharField(
        max_length=1000, blank=True, null=True
    )  # depleted reads


class TelevirRunQC(models.Model):
    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)
    performed = models.BooleanField(default=False)
    method = models.CharField(max_length=50, blank=True, null=True)
    args = models.CharField(max_length=50, blank=True, null=True)
    input_reads = models.CharField(max_length=50, blank=True, null=True)
    output_reads = models.CharField(max_length=50, blank=True, null=True)
    output_reads_percent = models.CharField(max_length=50, blank=True, null=True)

    class Meta:
        ordering = [
            "run",
        ]


class RunDetail(models.Model):
    name = models.CharField(
        max_length=100, db_index=True, blank=True, null=True
    )  # Create your models here.

    sample = models.ForeignKey(
        PIProject_Sample, blank=True, null=True, on_delete=models.CASCADE
    )

    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)

    max_depth = models.FloatField(blank=True, null=True)
    max_depthR = models.FloatField(blank=True, null=True)
    max_gaps = models.IntegerField(blank=True, null=True)
    max_prop = models.FloatField(blank=True, null=True)
    max_mapped = models.IntegerField(blank=True, null=True)

    input = models.CharField(max_length=300, blank=True, null=True)
    enriched_reads = models.IntegerField(blank=True, null=True, default=0)
    enriched_reads_percent = models.FloatField(blank=True, null=True, default=0)
    depleted_reads = models.IntegerField(blank=True, null=True, default=0)
    depleted_reads_percent = models.FloatField(blank=True, null=True, default=0)

    processed = models.CharField(max_length=300, blank=True, null=True)
    processed_percent = models.FloatField(blank=True, null=True)

    sift_preproc = models.BooleanField(blank=True)
    sift_remap = models.BooleanField(blank=True)
    sift_removed_pprc = models.CharField(max_length=300, blank=True, null=True)
    processing_final = models.CharField(max_length=300, blank=True, null=True)
    processing_final_percent = models.FloatField(blank=True, null=True)
    merged = models.BooleanField(blank=True)
    merged_number = models.IntegerField(blank=True, null=True)
    merged_files = models.CharField(max_length=1000, blank=True, null=True)

    class Meta:
        ordering = [
            "name",
        ]

    def __str__(self):
        return self.name


class RunAssembly(models.Model):
    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)
    sample = models.ForeignKey(
        PIProject_Sample, blank=True, null=True, on_delete=models.CASCADE
    )

    performed = models.BooleanField(default=False)

    assembly_contigs = models.CharField(max_length=1000, blank=True, null=True)

    method = models.CharField(
        max_length=20, blank=True, null=True
    )  # assembly method if any

    args = models.CharField(max_length=50, blank=True, null=True)  # assembly args

    contig_number = models.IntegerField(blank=True, null=True)

    contig_max = models.CharField(
        max_length=100, blank=True, null=True
    )  # max length of contig.
    contig_min = models.CharField(
        max_length=100, blank=True, null=True
    )  # min length of contig.
    contig_mean = models.CharField(
        max_length=100, blank=True, null=True
    )  # mean length of contig.

    contig_trim = models.CharField(
        max_length=100, blank=True, null=True
    )  # contig min len filter

    class Meta:
        ordering = [
            "method",
        ]

    def __str__(self):
        return self.method

    def delete_data(self):
        if os.path.isfile(self.assembly_contigs):
            os.remove(self.assembly_contigs)
            self.assembly_contigs = None

        self.save()


class ReadClassification(models.Model):
    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)
    sample = models.ForeignKey(
        PIProject_Sample, blank=True, null=True, on_delete=models.CASCADE
    )

    performed = models.BooleanField(default=False)

    method = models.CharField(
        max_length=20, blank=True, null=True
    )  # read classification method if any

    args = models.CharField(
        max_length=250, blank=True, null=True
    )  # read classification args
    db = models.CharField(
        max_length=500, blank=True, null=True
    )  # read classification db if any

    read_classification_report = models.CharField(
        max_length=1000, blank=True, null=True
    )  # read classification report

    classification_number = models.IntegerField(blank=True, null=True)
    classification_minhit = models.IntegerField(blank=True, null=True)

    success = models.BooleanField(default=False)

    class Meta:
        ordering = [
            "method",
        ]

    def __str__(self):
        return self.method


class ContigClassification(models.Model):
    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)
    sample = models.ForeignKey(
        PIProject_Sample, blank=True, null=True, on_delete=models.CASCADE
    )

    performed = models.BooleanField(default=False)

    method = models.CharField(
        max_length=250, blank=True, null=True
    )  # contig classification method if any

    args = models.CharField(
        max_length=250, blank=True, null=True
    )  # read classification args
    db = models.CharField(
        max_length=500, blank=True, null=True
    )  # read classification db if any

    contig_classification_report = models.CharField(
        max_length=1000, blank=True, null=True
    )  # contig classification report

    classification_number = models.IntegerField(blank=True, null=True)
    classification_minhit = models.IntegerField(blank=True, null=True)

    success = models.BooleanField(default=False)

    class Meta:
        ordering = [
            "method",
        ]

    def __str__(self):
        return self.method


class RawReference(models.Model):
    STATUS_MAPPED = 0
    STATUS_UNMAPPED = 1
    STATUS_MAPPING = 2
    STATUS_FAIL = 3

    STATUS_CHOICES = (
        (STATUS_MAPPED, "Mapped"),
        (STATUS_UNMAPPED, "Unmapped"),
        (STATUS_MAPPING, "Mapping"),
    )

    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)
    status = models.IntegerField(choices=STATUS_CHOICES, default=STATUS_UNMAPPED)

    taxid = models.CharField(max_length=100, blank=True, null=True)
    accid = models.CharField(max_length=100, blank=True, null=True)
    description = models.CharField(max_length=150, blank=True, null=True)
    counts = models.CharField(max_length=100, blank=True, null=True)
    classification_source = models.CharField(max_length=15, blank=True, null=True)
    panel = models.ForeignKey(
        ReferencePanel, on_delete=models.CASCADE, blank=True, null=True
    )

    @property
    def classification_source_str(self):
        if self.classification_source == "1":
            return "reads"

        if self.classification_source == "2":
            return "contigs"

        if self.classification_source == "3":
            return "reads / contigs"

        return "unknown"

    @property
    def read_counts(self):
        if self.classification_source == "1":
            return self.counts

        if self.classification_source == "2":
            return "0"

        if self.counts is None:
            return "0"

        if self.classification_source == "3":
            return self.counts.split("/")[0]

        return None

    @property
    def contig_counts(self):
        if self.classification_source == "1":
            return "0"

        if self.classification_source == "2":
            return self.counts

        if self.counts is None:
            return "0"

        if self.classification_source == "3":
            return self.counts.split("/")[1]

        return None

    @property
    def counts_int_array(self):
        return np.array([self.read_counts, self.contig_counts], dtype=np.int64)

    def update_raw_reference_status_mapped(self):
        self.status = self.STATUS_MAPPED
        self.save()

    def update_raw_reference_status_fail(self):
        self.status = self.STATUS_FAIL
        self.save()


from managing_files.models import Project as InsaFluProject
from managing_files.models import Reference as InsaFluReference


class MetaReference(models.Model):
    description = models.CharField(max_length=200, blank=True, null=True)
    created = models.DateTimeField(auto_now_add=True, blank=True, null=True)
    project = models.ForeignKey(
        Projects, on_delete=models.CASCADE, blank=True, null=True
    )
    file_path = models.CharField(max_length=300, blank=True, null=True)

    def __str__(self):
        return self.description

    @property
    def references(self) -> List[RawReference]:
        ref_maps = RawReferenceMap.objects.filter(reference=self)
        return [ref_map.raw_reference for ref_map in ref_maps]

    @property
    def taxids(self):
        return [ref.taxid for ref in self.references]

    @property
    def taxids_str(self):
        return "_".join(self.taxids)

    @property
    def accids(self):
        return [ref.accid for ref in self.references]

    @property
    def accids_str(self):
        return "_".join(self.accids)

    @property
    def descriptions(self):
        return [ref.description for ref in self.references]

    @property
    def description_first(self):
        return self.descriptions[0]

    @property
    def metaid(self):
        mapped_refs = self.references
        mapped_refs_ids = [ref.id for ref in mapped_refs]
        mapped_refids_sorted = sorted(mapped_refs_ids)
        return "_".join([str(refid) for refid in mapped_refids_sorted])


class RawReferenceMap(models.Model):
    reference = models.ForeignKey(
        MetaReference,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name="reference_map",
    )
    raw_reference = models.ForeignKey(
        RawReference,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name="raw_reference",
    )


class TeleFluProject(models.Model):
    televir_project = models.ForeignKey(Projects, on_delete=models.CASCADE)
    insaflu_project = models.ForeignKey(
        InsaFluProject,
        on_delete=models.CASCADE,
        blank=True,
        null=True,
        related_name="insaflu_project",
    )
    name = models.CharField(
        max_length=200,
        db_index=True,
        blank=False,
        verbose_name="Project name",
        default="nameless_project",
        validators=[no_space_validator],
    )
    description = models.TextField(default="", null=True, blank=True)
    raw_reference = models.ForeignKey(
        MetaReference,
        on_delete=models.CASCADE,
    )
    reference = models.ForeignKey(
        InsaFluReference, on_delete=models.CASCADE, blank=True, null=True
    )

    last_change_date = models.DateTimeField("Last change date", blank=True, null=True)
    is_deleted = models.BooleanField(default=False)

    @property
    def televir_references(self):
        if self.raw_reference is None:
            return None

        return self.raw_reference.references

    @property
    def owner(self):
        return self.televir_project.owner

    @property
    def technology(self):
        return self.televir_project.technology

    @property
    def project_media_directory(self):
        return os.path.join(
            PICS.media_directory,
            PICS.televir_subdirectory,
            str(self.televir_project.owner.pk),
            str(self.televir_project.pk),
            "teleflu",
            str(self.pk),
        )

    @property
    def project_static_directory(self):
        return os.path.join(
            PICS.static_directory,
            PICS.televir_subdirectory,
            str(self.televir_project.owner.pk),
            str(self.televir_project.pk),
            "teleflu",
            str(self.pk),
        )

    @property
    def project_teleflu_directory(self):
        return os.path.join(self.project_media_directory, "teleflu")

    @property
    def project_vcf_directory(self):
        return os.path.join(self.project_media_directory, "vcf")

    @property
    def project_vcf(self):
        return os.path.join(self.project_vcf_directory, "teleflu.vcf")

    @property
    def project_igv_report_media(self):
        return os.path.join(self.project_media_directory, "igv_report.html")

    @property
    def project_igv_report(self):
        return os.path.join(self.project_static_directory, "igv_report.html")

    @property
    def nsamples(self):

        return TeleFluSample.objects.filter(teleflu_project=self).count()


class TeleFluSample(models.Model):
    televir_sample = models.ForeignKey(PIProject_Sample, on_delete=models.CASCADE)
    teleflu_project = models.ForeignKey(
        TeleFluProject, on_delete=models.CASCADE, blank=True, null=True
    )


class TelefluMapping(models.Model):
    leaf = models.ForeignKey(
        SoftwareTreeNode, on_delete=models.CASCADE, blank=True, null=True
    )
    teleflu_project = models.ForeignKey(
        TeleFluProject, on_delete=models.CASCADE, blank=True, null=True
    )

    stacked_samples = models.ManyToManyField(TeleFluSample)

    @property
    def mapping_directory(self):
        return os.path.join(
            self.teleflu_project.project_teleflu_directory, str(self.leaf.index)
        )

    @property
    def stacked_samples_televir(self):
        return self.stacked_samples.values_list("televir_sample", flat=True)

    @property
    def mapping_vcf(self):
        return os.path.join(self.mapping_directory, "teleflu_stacked.vcf")

    @property
    def vcf_media_path(self):
        return self.mapping_vcf.replace(PICS.media_directory, "/media")

    @property
    def mapping_igv_report(self):
        return os.path.join(self.mapping_directory, "igv_report.html")

    @property
    def mapped_samples(self):

        accids = self.teleflu_project.raw_reference.accids
        samples = TeleFluSample.objects.filter(
            teleflu_project=self.teleflu_project
        ).values_list("televir_sample", flat=True)

        refs = RawReference.objects.filter(
            run__parameter_set__sample__in=samples,
            accid__in=accids,
            run__parameter_set__leaf=self.leaf,
            run__parameter_set__status=ParameterSet.STATUS_FINISHED,
        )

        sample_pks = list(set([ref.run.parameter_set.sample.pk for ref in refs]))

        return PIProject_Sample.objects.filter(pk__in=sample_pks)


class TelefluMappedSample(models.Model):

    teleflu_sample = models.ForeignKey(
        TeleFluSample, on_delete=models.CASCADE, blank=True, null=True
    )
    teleflu_mapping = models.ForeignKey(
        TelefluMapping, on_delete=models.CASCADE, blank=True, null=True
    )


class ReferenceTaxid(models.Model):
    taxid = models.CharField(max_length=100, blank=True, null=True)

    def __str__(self):
        return self.taxid


class ReferenceSourceFile(models.Model):
    file = models.CharField(max_length=100, blank=True, null=True)

    def __str__(self):
        return self.file


class ReferenceSource(models.Model):
    taxid = models.ForeignKey(
        ReferenceTaxid, blank=True, null=True, on_delete=models.CASCADE
    )
    accid = models.CharField(max_length=100, blank=True, null=True)
    description = models.CharField(max_length=300, blank=True, null=True)

    def __str__(self):
        return self.accid


class ReferenceSourceFileMap(models.Model):
    STATUS_GOOD = 0
    STATUS_CORRUPT = 1

    reference_source = models.ForeignKey(
        ReferenceSource, blank=True, null=True, on_delete=models.CASCADE
    )

    reference_source_file = models.ForeignKey(
        ReferenceSourceFile, blank=True, null=True, on_delete=models.CASCADE
    )

    status = models.IntegerField(default=STATUS_GOOD)

    @property
    def accid(self):
        return self.reference_source.accid

    @property
    def description(self):
        return self.reference_source.description

    @property
    def taxid(self):
        return self.reference_source.taxid


class ReferenceSourceMap(models.Model):
    STATUS_MAPPED = 0
    STATUS_UNMAPPED = 1
    STATUS_MAPPING = 2
    STATUS_FAIL = 3

    STATUS_CHOICES = (
        (STATUS_MAPPED, "Mapped"),
        (STATUS_UNMAPPED, "Unmapped"),
        (STATUS_MAPPING, "Mapping"),
    )

    reference_source = models.ForeignKey(
        ReferenceSource, blank=True, null=True, on_delete=models.CASCADE
    )

    raw_reference = models.ForeignKey(
        RawReference, blank=True, null=True, on_delete=models.CASCADE
    )

    status = models.IntegerField(choices=STATUS_CHOICES, default=STATUS_UNMAPPED)

    @property
    def taxid(self):
        return self.reference_source.taxid

    @property
    def accid(self):
        return self.reference_source.accid

    @property
    def description(self):
        return self.reference_source.description

    @property
    def counts(self):
        return self.raw_reference.counts

    @property
    def classification_source(self):
        return self.raw_reference.classification_source

    @property
    def classification_source_str(self):
        return self.raw_reference.classification_source_str

    @property
    def read_counts(self):
        return self.raw_reference.read_counts

    @property
    def contig_counts(self):
        return self.raw_reference.contig_counts

    @property
    def counts_int_array(self):
        return self.raw_reference.counts_int_array

    def update_reference_source_map_status_mapped(self):
        self.status = self.STATUS_MAPPED
        self.save()

    def update_reference_source_map_status_fail(self):
        self.status = self.STATUS_FAIL
        self.save()


class RunRemapMain(models.Model):
    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)
    sample = models.ForeignKey(
        PIProject_Sample, blank=True, null=True, on_delete=models.CASCADE
    )
    merged_log = models.CharField(max_length=350, blank=True, null=True)
    performed = models.BooleanField(default=False)

    method = models.CharField(
        max_length=350, blank=True, null=True
    )  # remap method if any

    found_total = models.IntegerField(blank=True, null=True)

    coverage_minimum = models.IntegerField(blank=True, null=True)
    coverage_maximum = models.IntegerField(blank=True, null=True)

    success = models.IntegerField(blank=True, null=True)
    remap_plan = models.CharField(max_length=1000, blank=True, null=True)

    class Meta:
        ordering = [
            "method",
        ]

    def __str__(self):
        return self.method


class ReferenceMap_Main(models.Model):
    reference = models.CharField(
        max_length=100, db_index=True, blank=True, null=True
    )  # Create your models here.

    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)

    sample = models.ForeignKey(
        PIProject_Sample, blank=True, null=True, on_delete=models.CASCADE
    )
    taxid = models.CharField(max_length=20, blank=True, null=True)
    reference_contig_str = models.CharField(max_length=100, blank=True, null=True)

    report = models.CharField(max_length=1000, blank=True, null=True)

    plotly_dotplot = models.TextField(blank=True, null=True)

    bam_file_path = models.CharField(max_length=1000, blank=True, null=True)
    bai_file_path = models.CharField(max_length=1000, blank=True, null=True)
    fasta_file_path = models.CharField(max_length=1000, blank=True, null=True)
    fai_file_path = models.CharField(max_length=1000, blank=True, null=True)
    mapped_subset_r1 = models.CharField(max_length=1000, blank=True, null=True)
    mapped_subset_r2 = models.CharField(max_length=1000, blank=True, null=True)
    mapped_subset_r1_fasta = models.CharField(max_length=1000, blank=True, null=True)
    mapped_subset_r2_fasta = models.CharField(max_length=1000, blank=True, null=True)
    vcf = models.CharField(max_length=1000, blank=True, null=True)

    class Meta:
        ordering = [
            "reference",
        ]

    def __str__(self):
        return self.reference

    def delete_data(self):
        if self.bam_file_path:
            if os.path.isfile(self.bam_file_path):
                os.remove(self.bam_file_path)
                self.bam_file_path = None

        if self.bai_file_path:
            if os.path.isfile(self.bai_file_path):
                os.remove(self.bai_file_path)
                self.bai_file_path = None

        if self.fasta_file_path:
            if os.path.isfile(self.fasta_file_path):
                os.remove(self.fasta_file_path)
                self.fasta_file_path = None

        if self.fai_file_path:
            if os.path.isfile(self.fai_file_path):
                os.remove(self.fai_file_path)
                self.fai_file_path = None

        if self.mapped_subset_r1:
            if os.path.isfile(self.mapped_subset_r1):
                os.remove(self.mapped_subset_r1)
                self.mapped_subset_r1 = None

        if self.mapped_subset_r2:
            if os.path.isfile(self.mapped_subset_r2):
                os.remove(self.mapped_subset_r2)
                self.mapped_subset_r2 = None

        if self.mapped_subset_r1_fasta:
            if os.path.isfile(self.mapped_subset_r1_fasta):
                os.remove(self.mapped_subset_r1_fasta)
                self.mapped_subset_r1_fasta = None

        if self.mapped_subset_r2_fasta:
            if os.path.isfile(self.mapped_subset_r2_fasta):
                os.remove(self.mapped_subset_r2_fasta)
                self.mapped_subset_r2_fasta = None

        if self.vcf:
            if os.path.isfile(self.vcf):
                os.remove(self.vcf)
                self.vcf = None

        self.save()


class FinalReport(models.Model):
    CONTROL_FLAG_NONE = 0
    CONTROL_FLAG_SOURCE = 1
    CONTROL_FLAG_PRESENT = 2
    CONTROL_FLAG_WARNING = 3

    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)
    sample = models.ForeignKey(
        PIProject_Sample, blank=True, null=True, on_delete=models.CASCADE
    )
    reference = models.CharField(
        max_length=100, db_index=True, blank=True, null=True
    )  # Create your models here.

    unique_id = models.CharField(max_length=20, blank=True, null=True)

    reference_length = models.IntegerField(blank=True, null=True)

    taxid = models.CharField(max_length=20, blank=True, null=True)
    simple_id = models.CharField(max_length=20, blank=True, null=True)
    description = models.CharField(max_length=450, blank=True, null=True)

    ref_db = models.CharField(max_length=400, blank=True, null=True)
    reference_contig_str = models.CharField(max_length=100, blank=True, null=True)

    accid = models.CharField(max_length=20, blank=True, null=True)
    coverage = models.FloatField(blank=True, null=True)
    windows_covered = models.CharField(max_length=20, blank=True, null=True)
    depth = models.FloatField(blank=True, null=True)
    depthR = models.FloatField(blank=True, null=True)
    mapped_reads = models.IntegerField(blank=True, null=True)
    ref_proportion = models.FloatField(blank=True, null=True)
    mapped_proportion = models.FloatField(blank=True, null=True)
    error_rate = models.FloatField(blank=True, null=True)
    quality_avg = models.FloatField(blank=True, null=True)
    ngaps = models.IntegerField(blank=True, null=True)
    mapping_success = models.CharField(max_length=20, blank=True, null=True)
    classification_success = models.CharField(max_length=20, blank=True, null=True)

    refa_dotplot = models.TextField(blank=True, null=True)
    refa_dotplot_exists = models.BooleanField(default=False)
    covplot = models.TextField(blank=True, null=True)
    covplot_exists = models.BooleanField(default=False)
    bam_path = models.CharField(max_length=1000, blank=True, null=True)
    bai_path = models.CharField(max_length=1000, blank=True, null=True)
    reference_path = models.CharField(max_length=1000, blank=True, null=True)
    reference_index_path = models.CharField(max_length=1000, blank=True, null=True)
    reference_assembly_paf = models.CharField(max_length=1000, blank=True, null=True)
    mapped_scaffolds_path = models.CharField(max_length=1000, blank=True, null=True)
    mapped_scaffolds_index_path = models.CharField(
        max_length=1000, blank=True, null=True
    )

    control_flag = models.IntegerField(default=CONTROL_FLAG_NONE)

    @property
    def in_control(self):
        return self.control_flag in [self.CONTROL_FLAG_PRESENT]


class ReferenceContigs(models.Model):
    reference = models.ForeignKey(
        ReferenceMap_Main, blank=True, null=True, on_delete=models.CASCADE
    )
    run = models.ForeignKey(RunMain, blank=True, null=True, on_delete=models.CASCADE)
    contig = models.CharField(max_length=100, blank=True, null=True)
    # length = models.CharField(max_length=100, blank=True, null=True)
    depth = models.CharField(max_length=100, blank=True, null=True)
    depthr = models.CharField(max_length=100, blank=True, null=True)
    coverage = models.CharField(max_length=100, blank=True, null=True)

    class Meta:
        ordering = [
            "reference",
        ]

    def __str__(self):
        return self.contig
