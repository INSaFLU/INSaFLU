import codecs
import os

from django.contrib.auth.models import User
from django.db import models
from django.utils.safestring import mark_safe
from managing_files.models import Sample
from django.core.validators import RegexValidator
from django.utils.translation import gettext_lazy as _
from django import forms

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


class SoftwareTree(models.Model):
    """"""

    version = models.IntegerField(default=0)
    date_created = models.DateTimeField(auto_now_add=True, blank=True, null=True)

    global_index = models.IntegerField(default=0)
    technology = models.CharField(
        max_length=100,
        name="technology",
        blank=True,
        null=True,
    )  # encoding

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


# class SoftwareTree_Path(models.Model):
#    software_tree = models.ForeignKey(SoftwareTree, on_delete=models.CASCADE)
#    software_tree_node = models.ForeignKey(SoftwareTreeNode, on_delete=models.CASCADE)
#    project = models.ForeignKey(Projects, on_delete=models.CASCADE)
#
#    class Meta:
#        ordering = ["software_tree_node"]
#


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

    class Meta:
        ordering = ["project__id", "-creation_date"]

    def __str__(self):
        return self.sample.name


class ParameterSet(models.Model):

    STATUS_NOT_STARTED = 0
    STATUS_RUNNING = 1
    STATUS_FINISHED = 2
    STATUS_ERROR = 3
    STATUS_QUEUED = 4
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

    def register_subprocess(self):
        self.status = self.STATUS_RUNNING
        self.save()

    def register_finished(self):
        self.status = self.STATUS_FINISHED
        self.save()

    def register_error(self):
        self.status = self.STATUS_ERROR
        self.save()

    def __str__(self):
        return self.sample.name + " " + str(self.leaf.index)


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
        max_length=50, blank=True, null=True
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
        max_length=50, blank=True, null=True
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
    remap_args = models.CharField(max_length=50, blank=True, null=True)

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
    description = models.CharField(max_length=100, blank=True, null=True)
    counts = models.CharField(max_length=100, blank=True, null=True)
    classification_source = models.CharField(max_length=15, blank=True, null=True)


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


class FinalReport(models.Model):

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
