import json
import os
from datetime import datetime
from operator import itemgetter

from django.conf import settings
# from django.db.models import Manager as GeoManager
from django.contrib.auth.models import User
# Create your models here.
from django.contrib.gis.db.models import GeoManager  # #  change to django  2.x
from django.contrib.gis.db.models import PointField
from django.db import models
from django.utils.safestring import mark_safe
from django.utils.translation import ugettext_lazy as _

from constants.constants import Constants, FileExtensions, FileType, TypePath
from constants.constants_mixed_infection import ConstantsMixedInfection
from constants.software_names import SoftwareNames
from fluwebvirus.formatChecker import ContentTypeRestrictedFileField
from manage_virus.constants_virus import ConstantsVirus
from manage_virus.models import IdentifyVirus
from settings.constants_settings import ConstantsSettings


def reference_directory_path(instance, filename):
    # file will be uploaded to MEDIA_ROOT/<filename>
    return "uploads/generic_data/user_{0}/{1}".format(instance.owner.id, filename)


def user_directory_path(instance, filename):
    # file will be uploaded to MEDIA_ROOT/<filename>
    return "uploads/generic_data/user_{0}/{1}".format(instance.owner.id, filename)


class SeasonReference(models.Model):
    """
    Each sample needs a dataset
    """

    name = models.CharField(max_length=100, db_index=True, blank=True, null=True)
    owner = models.ForeignKey(
        User,
        related_name="season_reference",
        blank=True,
        null=True,
        on_delete=models.CASCADE,
    )
    creation_date = models.DateTimeField(
        auto_now_add=True, verbose_name="Date creation"
    )

    def __str__(self):
        return self.name

    class Meta:
        ordering = [
            "name",
        ]


class MetaKey(models.Model):
    """
    Has meta tags to put values, for example, quality in the files, or samples
    """

    name = models.CharField(max_length=200, db_index=True, blank=True, null=True)

    def __str__(self):
        return self.name

    class Meta:
        ordering = [
            "name",
        ]


class Reference(models.Model):
    constants = Constants()

    ### species
    SPECIES_SARS_COV_2 = "SARS_COV_2"
    SPECIES_MPXV = "MPXV"
    SPECIES_INFLUENZA = "INFLUENZA"
    SPECIES_RSV = "RSV"
    SPECIES_NOT_SET = "NOT_SET"
    SPECIES_INFLUENZA_segment_four = "4"  ## Name of segment 4

    name = models.CharField(
        max_length=200, db_index=True, verbose_name="Reference name"
    )
    display_name = models.CharField(
        max_length=200, db_index=True, default="", verbose_name="Display name"
    )
    isolate_name = models.CharField(
        max_length=200, default="", verbose_name="Isolate Name"
    )
    creation_date = models.DateTimeField(
        auto_now_add=True, verbose_name="Uploaded Date"
    )

    ## Size 100K
    reference_fasta = ContentTypeRestrictedFileField(
        upload_to=reference_directory_path,
        content_types=["application/octet-stream"],
        max_upload_size=settings.MAX_REF_FASTA_FILE,
        blank=True,
        null=True,
        max_length=500,
    )
    reference_fasta_name = models.CharField(
        max_length=200, default="", verbose_name="Fasta file"
    )
    hash_reference_fasta = models.CharField(max_length=50, blank=True, null=True)

    ## Size 200K
    ## application/x-gameboy-rom because of 'gb' extension file of gbk
    reference_genbank = ContentTypeRestrictedFileField(
        upload_to=reference_directory_path,
        content_types=[
            "application/octet-stream",
            "application/x-gameboy-rom",
            "text/plain",
        ],
        max_upload_size=settings.MAX_REF_GENBANK_FILE,
        blank=True,
        null=True,
        max_length=500,
    )
    reference_genbank_name = models.CharField(
        max_length=200, default="", verbose_name="Genbank file"
    )
    hash_reference_genbank = models.CharField(max_length=50, blank=True, null=True)

    owner = models.ForeignKey(
        User, related_name="reference", blank=True, null=True, on_delete=models.CASCADE
    )
    is_obsolete = models.BooleanField(default=False, verbose_name="Obsolete")
    is_deleted = models.BooleanField(default=False, verbose_name="Deleted")
    number_of_locus = models.IntegerField(default=0, verbose_name="#Locus")

    season = models.ManyToManyField(SeasonReference)  ## can have the season
    description = models.CharField(
        max_length=500, default="", blank=True, null=True, verbose_name="Description"
    )

    ### if is deleted in file system
    is_deleted_in_file_system = models.BooleanField(
        default=False
    )  ## if this file was removed in file system
    date_deleted = models.DateTimeField(
        blank=True, null=True, verbose_name="Date attached"
    )  ## this date has the time of deleted by web page

    ### specie_tag, Has tag name of the specie;
    ### possible values SPECIES_SARS_COV_2, SPECIES_MPXV, etc...
    specie_tag = models.CharField(max_length=20, default="")

    def __str__(self):
        return self.name

    def get_reference_gbk(self, type_path):
        """
        get a path, type_path, from MEDIA_URL or MEDIA_ROOT
        """
        path_to_find = self.reference_genbank.name
        if type_path == TypePath.MEDIA_ROOT:
            if not path_to_find.startswith("/"):
                path_to_find = os.path.join(
                    getattr(settings, "MEDIA_ROOT", None), path_to_find
                )
        else:
            path_to_find = os.path.join(
                getattr(settings, "MEDIA_URL", None), path_to_find
            )
        return path_to_find

    def get_reference_fasta(self, type_path):
        """
        get a path, type_path, from MEDIA_URL or MEDIA_ROOT
        """
        path_to_find = self.reference_fasta.name
        if type_path == TypePath.MEDIA_ROOT:
            if not path_to_find.startswith("/"):
                path_to_find = os.path.join(
                    getattr(settings, "MEDIA_ROOT", None), path_to_find
                )
        else:
            path_to_find = os.path.join(
                getattr(settings, "MEDIA_URL", None), path_to_find
            )
        return path_to_find

    def get_reference_fasta_index(self, type_path):
        """
        get a fasta fai path, type_path, from MEDIA_URL or MEDIA_ROOT
        """
        return self.get_reference_fasta(type_path) + FileExtensions.FILE_FAI

    def get_reference_bed(self, type_path):
        """
        get a fasta bed path, type_path, from MEDIA_URL or MEDIA_ROOT
        """
        file_name = self.get_reference_gbk(type_path)
        return file_name[: file_name.rfind(".")] + FileExtensions.FILE_BED

    def get_reference_bed_index(self, type_path):
        """
        get a fasta bed.idx path, type_path, from MEDIA_URL or MEDIA_ROOT
        """
        return self.get_reference_bed(type_path) + FileExtensions.FILE_IDX

    def get_reference_fasta_web(self):
        """
        return web link for reference
        """
        out_file = self.get_reference_fasta(TypePath.MEDIA_ROOT)
        if os.path.exists(out_file):
            return mark_safe(
                '<a href="{}" download="{}"> {}</a>'.format(
                    self.get_reference_fasta(TypePath.MEDIA_URL),
                    os.path.basename(self.get_reference_fasta(TypePath.MEDIA_ROOT)),
                    self.constants.short_name(
                        self.reference_fasta_name, Constants.SHORT_NAME_LENGTH
                    ),
                )
            )
        return _("File not available.")

    def get_reference_gb_web(self):
        """
        return web link for reference
        """
        out_file = self.get_reference_fasta(TypePath.MEDIA_ROOT)
        if os.path.exists(out_file):
            return mark_safe(
                '<a href="{}" download="{}"> {}</a>'.format(
                    self.get_reference_gbk(TypePath.MEDIA_URL),
                    os.path.basename(self.get_reference_gbk(TypePath.MEDIA_ROOT)),
                    self.constants.short_name(
                        self.reference_genbank_name, Constants.SHORT_NAME_LENGTH
                    ),
                )
            )
        return _("File not available.")

    def get_gff3(self, type_path):
        """
        get GFF3 obtain form genbank
        :param type_path from MEDIA_URL or MEDIA_ROOT
        """
        path_to_find = self.reference_genbank.name
        path_to_find = (
            path_to_find[: path_to_find.rfind(".")] + FileExtensions.FILE_GFF3
        )
        if type_path == TypePath.MEDIA_ROOT:
            if not path_to_find.startswith("/"):
                path_to_find = os.path.join(
                    getattr(settings, "MEDIA_ROOT", None), path_to_find
                )
        else:
            path_to_find = os.path.join(
                getattr(settings, "MEDIA_URL", None), path_to_find
            )
        return path_to_find

    def get_gff3_with_gene_annotation(self, type_path):
        """
        get GFF3 obtain form genbank
        :param type_path from MEDIA_URL or MEDIA_ROOT
        """
        path_to_find = self.reference_genbank.name
        path_to_find = (
            path_to_find[: path_to_find.rfind(".")]
            + ".gene_annotation"
            + FileExtensions.FILE_GFF3
        )
        if type_path == TypePath.MEDIA_ROOT:
            if not path_to_find.startswith("/"):
                path_to_find = os.path.join(
                    getattr(settings, "MEDIA_ROOT", None), path_to_find
                )
        else:
            path_to_find = os.path.join(
                getattr(settings, "MEDIA_URL", None), path_to_find
            )
        return path_to_find

    def get_gff3_comulative_positions(self, type_path):
        """
        get GFF3 obtain form genbank
        :param type_path from MEDIA_URL or MEDIA_ROOT
        """
        path_to_find = self.reference_genbank.name
        path_to_find = (
            path_to_find[: path_to_find.rfind(".")]
            + ".comulative_positions"
            + FileExtensions.FILE_GFF3
        )
        if type_path == TypePath.MEDIA_ROOT:
            if not path_to_find.startswith("/"):
                path_to_find = os.path.join(
                    getattr(settings, "MEDIA_ROOT", None), path_to_find
                )
        else:
            path_to_find = os.path.join(
                getattr(settings, "MEDIA_URL", None), path_to_find
            )
        return path_to_find

    class Meta:
        verbose_name = "Reference"
        verbose_name_plural = "References"
        ordering = [
            "-creation_date",
        ]
        indexes = [
            models.Index(fields=["name"], name="name_idx"),
        ]


class MetaKeyReference(models.Model):
    """
    Relation ManyToMany in
    """

    meta_tag = models.ForeignKey(
        MetaKey, related_name="meta_key_reference", on_delete=models.CASCADE
    )
    reference = models.ForeignKey(
        Reference, related_name="meta_key_reference", on_delete=models.CASCADE
    )
    owner = models.ForeignKey(
        User, related_name="meta_key_reference", on_delete=models.CASCADE
    )
    creation_date = models.DateTimeField(
        "uploaded date", db_index=True, auto_now_add=True
    )
    value = models.CharField(default=Constants.META_KEY_VALUE_NOT_NEED, max_length=200)
    description = models.TextField(default="")

    class Meta:
        ordering = ["reference__id", "-creation_date"]

    def __str__(self):
        return self.meta_tag.name + " " + self.value + " " + self.description


class TagName(models.Model):
    """
    Has the tags to ticked the samples
    """

    name = models.CharField(max_length=100, db_index=True, blank=True, null=True)
    owner = models.ForeignKey(
        User, related_name="tag_name", blank=True, null=True, on_delete=models.CASCADE
    )
    creation_date = models.DateTimeField(
        auto_now_add=True, verbose_name="Uploaded Date"
    )
    is_meta_data = models.BooleanField(
        default=False
    )  ## if this tag belongs to meta data or not.

    def __str__(self):
        return self.name

    class Meta:
        ordering = [
            "name",
        ]


class DataSet(models.Model):
    """
    Each sample needs a dataset
    """

    name = models.CharField(max_length=100, db_index=True, blank=True, null=True)
    owner = models.ForeignKey(
        User, related_name="data_set", blank=True, null=True, on_delete=models.CASCADE
    )
    creation_date = models.DateTimeField(
        auto_now_add=True, verbose_name="Date creation"
    )

    def __str__(self):
        return self.name

    class Meta:
        ordering = [
            "creation_date",
            "name",
        ]


class VaccineStatus(models.Model):
    """
    Each sample needs a dataset
    """

    name = models.CharField(max_length=100, db_index=True, blank=True, null=True)
    owner = models.ForeignKey(
        User,
        related_name="vaccine_status",
        blank=True,
        null=True,
        on_delete=models.CASCADE,
    )
    creation_date = models.DateTimeField(
        auto_now_add=True, verbose_name="Date creation"
    )

    def __str__(self):
        return self.name

    class Meta:
        ordering = [
            "creation_date",
            "name",
        ]


class MixedInfectionsTag(models.Model):
    """
    Used to tag mixed infections
    """

    name = models.CharField(max_length=50, db_index=True, blank=True, null=True)

    def __str__(self):
        return self.name

    class Meta:
        ordering = [
            "name",
        ]


class MixedInfections(models.Model):
    """
    Used to identify mixed infections
    """

    tag = models.ForeignKey(
        MixedInfectionsTag,
        related_name="mixed_infections",
        blank=True,
        null=True,
        on_delete=models.CASCADE,
    )
    average_value = models.FloatField(default=0.0)
    description = models.TextField(default="")
    creation_date = models.DateTimeField("uploaded date", auto_now_add=True)
    last_change_date = models.DateTimeField("uploaded date", blank=True, null=True)
    has_master_vector = models.BooleanField(
        default=False
    )  ## if it has the master vector, has the vector to compare to all others

    ## and is not used in projectSample,
    ## It can change across time
    ## to trace the change of tag is set a metaValue in ProjectSample
    class Meta:
        ordering = [
            "tag",
        ]


class Sample(models.Model):
    """
    Sample, each sample has one or two files...
    """

    OUT_FILE_ABRICATE = "abricate.txt"
    DRAFT_CONTIGS = "_draft_contigs.fasta"
    SEGMENTS_CONTIGS = "_segments_to_contigs.tsv"
    SEGMENTS_READS = (
        "_segments_to_reads.tsv"  ### if the abricate is produced from reads, instead
    )

    ### consensus file
    OUT_CONSENSUS_FILE = "consensus.fasta"

    ### type of fastq.gz
    TYPE_OF_FASTQ_illumina = 0
    TYPE_OF_FASTQ_minion = 1
    TYPE_OF_FASTQ_not_defined = -1

    TYPE_SUBTYPE_LENGTH = 150

    ## to remove in future
    objects = models.Manager()  ## need to check this

    name = models.CharField(
        max_length=200, db_index=True, blank=True, null=True, verbose_name="Sample Name"
    )  ## This Id should match the prefix of the reads files (i.e. prefix_R1_001.fastq.gz /
    ##    prefix_R2_001.fastq.gz),
    date_of_onset = models.DateField("date of onset", blank=True, null=True)
    date_of_collection = models.DateField("date of collection", blank=True, null=True)
    date_of_receipt_lab = models.DateField("date of receipt lab", blank=True, null=True)
    week = models.IntegerField(blank=True, null=True)
    day = models.IntegerField(
        blank=True, null=True
    )  ## from "Date of onset” or “Date of collection”
    month = models.IntegerField(blank=True, null=True)
    year = models.IntegerField(blank=True, null=True)
    vaccine_status = models.ForeignKey(
        VaccineStatus,
        related_name="sample",
        blank=True,
        null=True,
        on_delete=models.CASCADE,
    )
    creation_date = models.DateTimeField(
        auto_now_add=True, verbose_name="Uploaded Date"
    )
    is_deleted = models.BooleanField(default=False)
    is_obsolete = models.BooleanField(default=False)
    owner = models.ForeignKey(
        User, related_name="sample", blank=True, null=True, on_delete=models.CASCADE
    )
    data_set = models.ForeignKey(
        DataSet, related_name="sample", blank=True, null=True, on_delete=models.CASCADE
    )
    geo_local = PointField(null=True, blank=True, srid=4326)
    ## 4326 which means latitude and longitude
    geo_manager = GeoManager()

    ### Type/Subtype Virus
    identify_virus = models.ManyToManyField(IdentifyVirus)
    type_subtype = models.CharField(
        verbose_name="Classification", max_length=150, blank=True, null=True
    )  ## has the type/subtype collected
    number_alerts = models.IntegerField(
        verbose_name="Alerts", default=0, blank=True, null=True
    )  ## has the number of alerts
    mixed_infections_tag = models.ForeignKey(
        MixedInfectionsTag,
        verbose_name="Putative Mixed Infection",
        related_name="sample",
        blank=True,
        null=True,
        on_delete=models.CASCADE,
    )
    ### has the tag of Yes/No mixed infection

    ## many to many relation
    tag_names = models.ManyToManyField(TagName, through="TagNames")

    ### files
    is_valid_1 = models.BooleanField(default=False)
    file_name_1 = models.CharField(max_length=200, blank=True, null=True)
    candidate_file_name_1 = models.CharField(max_length=200, blank=True, null=True)

    ## 30M
    path_name_1 = ContentTypeRestrictedFileField(
        upload_to=user_directory_path,
        blank=True,
        null=True,
        content_types=[
            "application/octet-stream",
            "application/gzip",
            "application/x-gzip",
        ],
        max_upload_size=(
            settings.MAX_FASTQ_FILE_WITH_DOWNSIZE
            if settings.DOWN_SIZE_FASTQ_FILES
            else settings.MAX_FASTQ_FILE_UPLOAD
        ),
        max_length=500,
    )
    is_valid_2 = models.BooleanField(default=False)
    file_name_2 = models.CharField(max_length=200, blank=True, null=True)
    candidate_file_name_2 = models.CharField(max_length=200, blank=True, null=True)
    path_name_2 = ContentTypeRestrictedFileField(
        upload_to=user_directory_path,
        blank=True,
        null=True,
        content_types=[
            "application/octet-stream",
            "application/gzip",
            "application/x-gzip",
        ],
        max_upload_size=(
            settings.MAX_FASTQ_FILE_WITH_DOWNSIZE
            if settings.DOWN_SIZE_FASTQ_FILES
            else settings.MAX_FASTQ_FILE_UPLOAD
        ),
        max_length=500,
    )

    ### type of fastq.gz, Default Illumina,
    type_of_fastq = models.SmallIntegerField(
        verbose_name="File type", default=TYPE_OF_FASTQ_not_defined
    )  ## has the type of the fastq files

    ## has files, the user can upload the files after
    has_files = models.BooleanField(default=False)

    ###    has the flag indicating that the sample can be processed by projects
    is_ready_for_projects = models.BooleanField(default=False)
    ###    has the flag indicating that the sample has end of processing, if False it is ready for process
    is_sample_in_the_queue = models.BooleanField(default=False)

    ### if is deleted in file system
    is_deleted_in_file_system = models.BooleanField(
        default=False
    )  ## if this file was removed in file system
    date_deleted = models.DateTimeField(
        blank=True, null=True, verbose_name="Date deleted"
    )  ## this date has the time of deleted by web page

    ### deleted processed fastq, but not delete the sample (trimmomatic/Rabbit/etc)
    is_deleted_processed_fastq = models.BooleanField(
        default=False
    )  ## if this files were removed in file system
    date_deleted_processed_fastq = models.DateTimeField(
        blank=True, null=True, verbose_name="Date fastq deleted"
    )  ## this date has the time of fastq deleted

    def __str__(self):
        return self.name

    class Meta:
        ordering = [
            "-creation_date",
        ]

    def get_file_names(self):
        sz_return = "" if self.file_name_1 == None else self.file_name_1
        sz_return += "/" if len(sz_return) > 0 and self.file_name_2 is not None else ""
        sz_return += "" if self.file_name_2 == None else self.file_name_2
        return sz_return

    def exist_file_2(self):
        if (
            self.path_name_2 is None
            or self.path_name_2.name is None
            or len(self.path_name_2.name) == 0
        ):
            return False
        return True

    def get_abricate_output(self, type_path, b_gzip_file=False):
        """
        type_path = [MEDIA_ROOT, MEDIA_URL]
        Return the file name of the abricate output base on fastq File input
        path it's a FileField instance, or a string
        Can be zipped
        """
        return os.path.join(
            self.__get_path__(type_path, True),
            self.OUT_FILE_ABRICATE + (FileExtensions.FILE_GZ if b_gzip_file else ""),
        )

    def get_draft_contigs_output(self, type_path):
        """
        type_path = [MEDIA_ROOT, MEDIA_URL]
        Return the file name of the abricate output base on fastq File input
        path it's a FileField instance, or a string
        """
        return os.path.join(
            self.__get_path__(type_path, True),
            "{}{}".format(self.name.replace(" ", "_"), self.DRAFT_CONTIGS),
        )

    def get_draft_contigs_abricate_output(self, type_path):
        """
        type_path = [MEDIA_ROOT, MEDIA_URL]
        Return the file name of the abricate output base on fastq File input
        path it's a FileField instance, or a string
        """
        return os.path.join(
            self.__get_path__(type_path, True),
            "{}{}".format(self.name.replace(" ", "_"), self.SEGMENTS_CONTIGS),
        )

    def get_draft_reads_abricate_output(self, type_path):
        """
        If the obreciate result came from Reads instead contigs. This appends in ONT reads. It is not produced contigs
        type_path = [MEDIA_ROOT, MEDIA_URL]
        Return the file name of the abricate output base on fastq File input
        path it's a FileField instance, or a string
        """
        return os.path.join(
            self.__get_path__(type_path, True),
            "{}{}".format(self.name.replace(" ", "_"), self.SEGMENTS_READS),
        )

    def get_trimmomatic_file(self, type_path, b_first_file):
        """
        get the trimmomatic files, it's going to be use for all processing
        type_path = [MEDIA_ROOT, MEDIA_URL]
        """
        if not b_first_file and not self.exist_file_2():
            return None
        path_out = os.path.join(
            self.__get_path__(type_path, b_first_file),
            Constants.DIR_PROCESSED_PROCESSED,
        )
        if type_path == TypePath.MEDIA_ROOT:
            os.makedirs(path_out, mode=0o755, exist_ok=True)
        return os.path.join(
            path_out, self.name + ("_1P.fastq.gz" if b_first_file else "_2P.fastq.gz")
        )

    def get_nanofilt_file(self, type_path):
        """
        get the trimmomatic files, it's going to be use for all processing
        type_path = [MEDIA_ROOT, MEDIA_URL]
        """
        b_first_file = True
        path_out = os.path.join(
            self.__get_path__(type_path, b_first_file),
            Constants.DIR_PROCESSED_PROCESSED,
        )
        if type_path == TypePath.MEDIA_ROOT:
            os.makedirs(path_out, mode=0o755, exist_ok=True)
        return os.path.join(path_out, self.name + ".fastq.gz")

    def get_fastq_trimmomatic(self, type_path, b_first_file):
        """
        return fastQC output first step
        """
        if not b_first_file and not self.exist_file_2():
            return None

        path_out = os.path.join(
            self.__get_path__(type_path, b_first_file),
            Constants.DIR_PROCESSED_PROCESSED,
        )
        if type_path == TypePath.MEDIA_ROOT:
            os.makedirs(path_out, mode=0o755, exist_ok=True)
        return os.path.join(
            path_out,
            self.name + ("_1P_fastqc.html" if b_first_file else "_2P_fastqc.html"),
        )

    def get_rabbitQC_nanofilt(self, type_path):
        """
        return rabbitQC output first step
        """
        b_first_file = True
        path_out = os.path.join(
            self.__get_path__(type_path, b_first_file),
            Constants.DIR_PROCESSED_PROCESSED,
        )
        if type_path == TypePath.MEDIA_ROOT:
            os.makedirs(path_out, mode=0o755, exist_ok=True)
        return os.path.join(path_out, self.name + "_rabbitQC.html")

    def get_fastq(self, type_path, b_first_file):
        """
        return fastq output first step, from MEDIA_URL or MEDIA_ROOT
        """
        if not b_first_file and not self.exist_file_2():
            return None
        return os.path.join(
            self.__get_path__(type_path, b_first_file),
            self.file_name_1 if b_first_file else self.file_name_2,
        )

    def get_fastq_available(self, type_path, b_first_file=True):
        """
        gets the fastq available, if not trimmomatic/nanofilt ran, return fastq
        try first trimmomatic/nanofilt, then return fastq
        """
        file_name = self.get_trimmomatic_file(type_path, b_first_file)
        if (file_name is not None) and os.path.exists(file_name):
            return file_name
        file_name = self.get_nanofilt_file(type_path)
        if (file_name is not None) and os.path.exists(file_name):
            return file_name
        return self.get_fastq(type_path, b_first_file)

    def is_original_fastq_removed(self):
        """
        Test if original fastq were removed already
        """
        return not os.path.exists(self.get_fastq(TypePath.MEDIA_ROOT, True))

    def is_processed_fastq_deleted(self):
        return self.is_deleted_processed_fastq

    def get_fastqc_output(self, type_path, b_first_file):
        """
        return fastqc output first quality control
        can be generic, also for rabbitQC
        """
        if not b_first_file and not self.exist_file_2():
            return None
        return os.path.join(
            self.__get_path__(type_path, b_first_file),
            (
                self.file_name_1.replace(
                    FileExtensions.FILE_FASTQ_GZ, "_fastqc.html"
                ).replace(FileExtensions.FILE_FQ_GZ, "_fastqc.html")
                if b_first_file
                else self.file_name_2.replace(
                    FileExtensions.FILE_FASTQ_GZ, "_fastqc.html"
                ).replace(FileExtensions.FILE_FQ_GZ, "_fastqc.html")
            ),
        )

    def get_rabbitQC_output(self, type_path):
        """
        return fastq output second step
        can be generic, also for rabbitQC
        """
        b_first_file = True
        return os.path.join(
            self.__get_path__(type_path, b_first_file),
            self.file_name_1.replace(
                FileExtensions.FILE_FASTQ_GZ, "_rabbitQC.html"
            ).replace(FileExtensions.FILE_FQ_GZ, "_rabbitQC.html"),
        )

    def __get_path__(self, type_path, b_first_file):
        """
        get a path, from MEDIA_URL or MEDIA_ROOT
        """
        if b_first_file:
            path_to_find = os.path.dirname(self.path_name_1.name)
            if type_path == TypePath.MEDIA_ROOT:
                if not path_to_find.startswith("/"):
                    path_to_find = os.path.join(
                        getattr(settings, "MEDIA_ROOT", None), path_to_find
                    )
            else:
                path_to_find = os.path.join(
                    getattr(settings, "MEDIA_URL", None), path_to_find
                )
        else:
            path_to_find = os.path.dirname(self.path_name_2.name)
            if type_path == TypePath.MEDIA_ROOT:
                if not path_to_find.startswith("/"):
                    path_to_find = os.path.join(
                        getattr(settings, "MEDIA_ROOT", None), path_to_find
                    )
            else:
                path_to_find = os.path.join(
                    getattr(settings, "MEDIA_URL", None), path_to_find
                )
        return path_to_find

    def get_type_sub_type(self):
        """
        return a string with Type/Subtpye
        """
        return IdentifyVirus().classify(self.identify_virus.all())

    def get_mixed_infection(self):
        """
        mixed infection based on the table static/mixed_infections/mixed_infections.xls
        return tuble (tag_mixed_infection, alert, message)
        tag_mixed_infection: ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO or ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES
        alert: positive number, or zero
        message: None, empty or the string with the error message
        """
        return ConstantsMixedInfection().get_mixed_infection(self.identify_virus.all())

    def get_is_ready_for_projects(self):
        """
        need to be true to be ready for projects
        """
        return (
            self.is_ready_for_projects and not self.is_obsolete and not self.is_deleted
        )

    def get_vect_tag_names_and_value(self):
        """
        return [[header1, value1], [header2, value2], [header3, value3], ...]
        """
        vect_out = []
        query_set = TagNames.objects.filter(sample=self)
        for tag_names in query_set:
            vect_out.append([tag_names.tag_name.name, tag_names.value])
        return sorted(vect_out, key=itemgetter(1))

    def get_tag_names(self):
        """
        get the tag names grouped by a number
        """
        query_set = TagNames.objects.filter(sample=self)
        if query_set.count() == 0:
            return None
        return query_set

    def is_type_fastq_gz_sequencing(self, type_of_fastq=TYPE_OF_FASTQ_illumina):
        """
        :param type_of_fastq, can be TYPE_OF_FASTQ_illumina, TYPE_OF_FASTQ_minion and so on
        """
        return self.type_of_fastq == type_of_fastq

    def set_type_of_fastq_sequencing(self, type_of_fastq):
        """
        :parm type_of_fastq Constants.FORMAT_FASTQ_illumina or Constants.FORMAT_FASTQ_ont
        """
        if type_of_fastq == Constants.FORMAT_FASTQ_illumina:
            self.type_of_fastq = Sample.TYPE_OF_FASTQ_illumina
        if type_of_fastq == Constants.FORMAT_FASTQ_ont:
            self.type_of_fastq = Sample.TYPE_OF_FASTQ_minion

    def get_type_technology(self):
        if self.type_of_fastq == Sample.TYPE_OF_FASTQ_illumina:
            return ConstantsSettings.TECHNOLOGY_illumina
        if self.type_of_fastq == Sample.TYPE_OF_FASTQ_minion:
            return ConstantsSettings.TECHNOLOGY_minion
        return "Not defined"


class TagNames(models.Model):
    value = models.CharField(max_length=150)
    tag_name = models.ForeignKey(TagName, on_delete=models.CASCADE)
    sample = models.ForeignKey(Sample, on_delete=models.CASCADE)

    def __str__(self):
        return self.value


class MetaKeySample(models.Model):
    """
    Relation ManyToMany in
    """

    meta_tag = models.ForeignKey(
        MetaKey, related_name="meta_key_sample", on_delete=models.CASCADE
    )
    sample = models.ForeignKey(
        Sample, related_name="meta_key_sample", on_delete=models.CASCADE
    )
    owner = models.ForeignKey(
        User, related_name="meta_key_sample", on_delete=models.CASCADE
    )
    creation_date = models.DateTimeField(
        "uploaded date", db_index=True, auto_now_add=True
    )
    value = models.CharField(default=Constants.META_KEY_VALUE_NOT_NEED, max_length=200)
    description = models.TextField(default="")

    class Meta:
        ordering = ["sample__id", "-creation_date"]

    def __str__(self):
        return self.meta_tag.name + " " + self.value + " " + self.description


class Project(models.Model):
    """
    Project for run a pipeline in the data,
    It's possible to run several times and add data to the project
    """

    ### has the path for main result
    PATH_MAIN_RESULT = "main_result"

    PROJECT_FILE_NAME_MAFFT = "Alignment_nt_All.fasta"
    PROJECT_FILE_NAME_FASTA = "All_nt.fasta"
    PROJECT_FILE_NAME_FASTTREE = "Tree_ML_All.nwk"
    PROJECT_FILE_NAME_FASTTREE_tree = "Tree_ML_All.tree"
    PROJECT_FILE_NAME_nex = "Alignment_nt_All.nex"
    PROJECT_FILE_NAME_COVERAGE = "coverage.tsv"
    PROJECT_FILE_NAME_TOTAL_VARIATIONS = "proportions_iSNVs_graph.tsv"
    PROJECT_FILE_NAME_TAB_VARIATIONS_SNIPPY = "validated_variants.tsv"
    PROJECT_FILE_NAME_TAB_VARIATIONS_FREEBAYES = "validated_minor_iSNVs.tsv"  ## remove del and ins and everything bigger than >50
    PROJECT_FILE_NAME_TAB_VARIATIONS_FREEBAYES_with_snps_indels = "validated_minor_inc_indels.tsv"  ## with snps, del and ins and everything bigger than >50
    ##                    MIGUEL
    ##                    "Minor intra-host variants (inc. indels):"
    ## freebayes_variants_file_snp_indel
    PERCENTAGE_validated_minor_variants = 51  ## only pass <= 50
    PROJECT_FILE_NAME_SAMPLE_RESULT_TSV = "Sample_list.tsv"  ### first column ID instead of 'sample name' to be compatible with Phandango e Microreact
    PROJECT_FILE_NAME_SAMPLE_RESULT_CSV = "Sample_list.csv"  ### first column ID instead of 'sample name' to be compatible with Phandango e Microreact
    PROJECT_FILE_NAME_SAMPLE_RESULT_SETTINGS_TSV = "Sample_list_settings.tsv"  ### first column ID instead of 'sample name' to be compatible with Phandango e Microreact
    PROJECT_FILE_NAME_SAMPLE_RESULT_SETTINGS_CSV = "Sample_list_settings.csv"  ### first column ID instead of 'sample name' to be compatible with Phandango e Microreact
    PROJECT_FILE_NAME_SAMPLE_RESULT_CSV_simple = "Sample_list_simple.csv"  ### first column must be ID because of manging_files.ajax_views.show_phylo_canvas
    ### this file is only used for to show the manging_files.views.ShowSampleProjectsView
    PROJECT_FILE_NAME_SAMPLE_RESULT_json = "Sample_list_simple.json"  ### first column ID instead of 'sample name' to be compatible with Phandango e Microreact, to download to
    PROJECT_FILE_NAME_SAMPLE_RESULT_all_consensus = (
        "AllConsensus.fasta"  ### all consensus sequences for a project sample
    )
    PROJECT_FILE_NAME_SAMPLE_mask_all_consensus = (
        "mask_all_consensus"  ### masking all consensus, defined by user
    )

    PROJECT_FILE_NAME_Pangolin_lineage = (
        "PangolinLineage.csv"  ### has the result of pangolin lineage
    )

    PROJECT_FILE_NAME_Aln2pheno_report_COG_UK = (
        "aln2pheno_final_report_COG_UK.tsv"  ### has results of aln2pheno
    )
    PROJECT_FILE_NAME_Aln2pheno_flagged_COG_UK = (
        "aln2pheno_flagged_mutation_report_COG_UK.tsv"  ### has results of aln2pheno
    )
    PROJECT_FILE_NAME_Aln2pheno_report_pokay = (
        "aln2pheno_final_report_pokay.tsv"  ### has results of aln2pheno
    )
    PROJECT_FILE_NAME_Aln2pheno_flagged_pokay = (
        "aln2pheno_flagged_mutation_report_pokay.tsv"  ### has results of aln2pheno
    )
    PROJECT_FILE_NAME_Aln2pheno_report_carabelli = "aln2pheno_final_report_EpitopeResidues_Carabelli_2023.tsv"  ### has results of aln2pheno
    PROJECT_FILE_NAME_Aln2pheno_flagged_carabelli = "aln2pheno_flagged_mutation_report_EpitopeResidues_Carabelli_2023.tsv"  ### has results of aln2pheno
    PROJECT_FILE_NAME_Aln2pheno_zip = "aln2pheno.zip"  ### has results of aln2pheno

    PROJECT_FILE_NAME_all_files_zipped = "AllFiles.zip"  ### Several files zipped

    ## put the type file here to clean if there isn't enough sequences to create the trees and alignments
    vect_clean_file = [
        PROJECT_FILE_NAME_MAFFT,
        PROJECT_FILE_NAME_FASTTREE,
        PROJECT_FILE_NAME_FASTTREE_tree,
        PROJECT_FILE_NAME_nex,
        PROJECT_FILE_NAME_FASTA,
    ]

    vect_exclude_clean_file_from_proteins = [PROJECT_FILE_NAME_FASTA]

    ## obsolete
    PROJECT_FILE_NAME_GRAPH_MINO_VAR_HTML = "graph_minor_var.html"
    PROJECT_FILE_NAME_GRAPH_MINO_VAR_PNG = "graph_minor_var.png"
    ## end obsolete

    ### this is only to join with other names
    PROJECT_FILE_NAME_FASTTREE_element = "Tree_ML"
    PROJECT_FILE_NAME_MAFFT_element_nt = "Alignment_nt"
    PROJECT_FILE_NAME_MAFFT_element_aa = "Alignment_aa"
    PROJECT_FILE_NAME_FASTA_element = "Sequences_nt"

    name = models.CharField(
        max_length=200,
        db_index=True,
        blank=True,
        null=True,
        verbose_name="Project name",
    )
    owner = models.ForeignKey(
        User, related_name="project", blank=True, null=True, on_delete=models.CASCADE
    )
    reference = models.ForeignKey(
        Reference,
        related_name="project",
        blank=True,
        null=True,
        on_delete=models.CASCADE,
    )
    creation_date = models.DateTimeField("Uploaded date", auto_now_add=True)
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

    def __str__(self):
        return self.name

    class Meta:
        ordering = [
            "-creation_date",
        ]

    def get_global_file_by_element(self, type_path, element, file_name):
        """
        type_path: constants.TypePath -> MEDIA_ROOT, MEDIA_URL
        element: element name or None
        file_name: Project.PROJECT_FILE_NAME_MAFFT, ....
        """
        if self.PROJECT_FILE_NAME_MAFFT == file_name:
            return os.path.join(
                self.__get_global_path__(type_path, element),
                "{}_{}.fasta".format(self.PROJECT_FILE_NAME_MAFFT_element_nt, element),
            )
        if self.PROJECT_FILE_NAME_FASTA == file_name:
            return os.path.join(
                self.__get_global_path__(type_path, element),
                "{}_{}.fasta".format(self.PROJECT_FILE_NAME_FASTA_element, element),
            )
        if self.PROJECT_FILE_NAME_FASTTREE == file_name:
            return os.path.join(
                self.__get_global_path__(type_path, element),
                "{}_{}.nwk".format(self.PROJECT_FILE_NAME_FASTTREE_element, element),
            )
        if self.PROJECT_FILE_NAME_FASTTREE_tree == file_name:
            return os.path.join(
                self.__get_global_path__(type_path, element),
                "{}_{}.tree".format(self.PROJECT_FILE_NAME_FASTTREE_element, element),
            )
        if self.PROJECT_FILE_NAME_nex == file_name:
            return os.path.join(
                self.__get_global_path__(type_path, element),
                "{}_{}.nex".format(self.PROJECT_FILE_NAME_MAFFT_element_nt, element),
            )
        return None

    def get_global_file_by_element_and_cds(self, type_path, element, CDS, file_name):
        """
        get file names for proteins
        type_path: constants.TypePath -> MEDIA_ROOT, MEDIA_URL
        element: element name or None
        file_name: Project.PROJECT_FILE_NAME_MAFFT, ....
        """
        CDS = self._clean_name(CDS)
        if self.PROJECT_FILE_NAME_MAFFT == file_name:  ## protein alignement
            return os.path.join(
                self.__get_global_path__(type_path, element),
                "{}_{}_{}.fasta".format(
                    self.PROJECT_FILE_NAME_MAFFT_element_aa, element, CDS
                ),
            )
        if self.PROJECT_FILE_NAME_FASTTREE == file_name:
            return os.path.join(
                self.__get_global_path__(type_path, element),
                "{}_{}_{}.nwk".format(
                    self.PROJECT_FILE_NAME_FASTTREE_element, element, CDS
                ),
            )
        if self.PROJECT_FILE_NAME_FASTTREE_tree == file_name:
            return os.path.join(
                self.__get_global_path__(type_path, element),
                "{}_{}_{}.tree".format(
                    self.PROJECT_FILE_NAME_FASTTREE_element, element, CDS
                ),
            )
        if self.PROJECT_FILE_NAME_nex == file_name:
            return os.path.join(
                self.__get_global_path__(type_path, element),
                "{}_{}_{}.nex".format(
                    self.PROJECT_FILE_NAME_MAFFT_element_aa, element, CDS
                ),
            )
        return None

    def _clean_name(
        self,
        name_to_clean,
        dict_to_clean={
            " ": "_",
            "(": "",
            ")": "",
            "$": "",
            "#": "",
            "&": "",
            "/": "",
            "\\": "",
            "-": "_",
        },
    ):
        """
        clean a name based on dictionary, dict_to_clean = { ' ' : '_', '(' : '' , ')' : '' }
        """
        for key in dict_to_clean:
            name_to_clean = name_to_clean.replace(key, dict_to_clean[key])
        return name_to_clean

    def get_clean_project_name(self):
        return self._clean_name(self.name)

    def get_global_file_by_project(self, type_path, file_name):
        """
        type_path: constants.TypePath -> MEDIA_ROOT, MEDIA_URL
        file_name: Project.PROJECT_FILE_NAME_MAFFT, ....
        """
        return os.path.join(self.__get_global_path__(type_path, None), file_name)

    def get_global_file_by_project_web(self, file_name):
        out_file = self.get_global_file_by_project(TypePath.MEDIA_ROOT, file_name)
        if os.path.exists(out_file):
            return mark_safe(
                '<a href="{}" download> {}</a>'.format(
                    self.get_global_file_by_project(TypePath.MEDIA_URL, file_name),
                    file_name,
                )
            )
        return _("File not available yet.")

    def __get_user_result_global_directory_path__(self, element):
        # file will be uploaded to MEDIA_ROOT/<filename>
        if not element is None and len(element) > 0:
            return (
                Constants.DIR_PROCESSED_FILES_PROJECT
                + "/user_{0}/project_{1}/{2}/{3}".format(
                    self.owner.id, self.pk, self.PATH_MAIN_RESULT, element
                )
            )
        return (
            Constants.DIR_PROCESSED_FILES_PROJECT
            + "/user_{0}/project_{1}/{2}".format(
                self.owner.id, self.pk, self.PATH_MAIN_RESULT
            )
        )

    def __get_global_path__(self, type_path, element):
        """
        get a path, from MEDIA_URL or MEDIA_ROOT
        """
        path_to_find = self.__get_user_result_global_directory_path__(element)
        if type_path == TypePath.MEDIA_ROOT:
            if not path_to_find.startswith("/"):
                path_to_find = os.path.join(
                    getattr(settings, "MEDIA_ROOT", None), path_to_find
                )
        ### URL
        else:
            path_to_find = os.path.join(
                getattr(settings, "MEDIA_URL", None), path_to_find
            )
        return path_to_find


class MetaKeyProject(models.Model):
    """
    Relation ManyToMany in
    """

    meta_tag = models.ForeignKey(
        MetaKey, related_name="meta_key_project", on_delete=models.CASCADE
    )
    project = models.ForeignKey(
        Project, related_name="meta_key_project", on_delete=models.CASCADE
    )
    owner = models.ForeignKey(
        User, related_name="meta_key_project", on_delete=models.CASCADE
    )
    creation_date = models.DateTimeField(
        "uploaded date", db_index=True, auto_now_add=True
    )
    value = models.CharField(default=Constants.META_KEY_VALUE_NOT_NEED, max_length=200)
    description = models.TextField(default="")

    class Meta:
        ordering = ["project__id", "-creation_date"]

    def __str__(self):
        return self.meta_tag.name + " " + self.value + " " + self.description


class CountVariations(models.Model):
    """
    has the number of variations for a sample
    """

    total = models.PositiveIntegerField(default=0)
    var_less_50 = models.PositiveIntegerField(default=0)
    var_bigger_50_90 = models.PositiveIntegerField(default=0)
    var_bigger_90 = models.PositiveIntegerField(default=0)

    def __str__(self):
        return "Total: {} Less 50:{} 50<Value<90 :{}  Bigger >90:{}".format(
            self.total, self.var_less_50, self.var_bigger_50_90, self.var_bigger_90
        )


class ProjectSample(models.Model):
    constants = Constants()

    OUT_FILE_ABRICATE = "abricate.txt"
    PATH_MAIN_RESULT = "main_result"
    PREFIX_FILE_COVERAGE = "coverage"
    FILE_CONSENSUS_FILE = "Consensus_"
    FILE_SNIPPY_TAB = "validated_variants_sample_"
    FILE_FREEBAYES_TAB = "validated_minor_iSNVs_sample_"
    FILE_FREEBAYES_TAB_with_indels = "validated_minor_inc_indels_sample_"

    project = models.ForeignKey(
        Project,
        related_name="project_samples",
        blank=True,
        null=True,
        on_delete=models.CASCADE,
    )
    sample = models.ForeignKey(
        Sample,
        related_name="project_samples",
        blank=True,
        null=True,
        on_delete=models.CASCADE,
    )
    mixed_infections = models.ForeignKey(
        MixedInfections,
        related_name="project_samples",
        blank=True,
        null=True,
        on_delete=models.CASCADE,
    )
    count_variations = models.ForeignKey(
        CountVariations,
        related_name="project_samples",
        blank=True,
        null=True,
        on_delete=models.CASCADE,
    )
    classification = models.CharField(
        verbose_name="Classification", max_length=150, blank=True, null=True
    )
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

    class Meta:
        ordering = ["project__id", "-creation_date"]

    def __str__(self):
        return self.project.name

    def get_global_file_by_element(
        self, type_path, prefix_file_name, sequence_name, extension
    ):
        """
        type_path: constants.TypePath -> MEDIA_ROOT, MEDIA_URL
        prefix_file_name: ProjectSample.PREFIX_FILE_COVERAGE
        sequence_name: sequence name
        extension: FileExtensions.FILE_PNG
        """
        return os.path.join(
            self.__get_path__(type_path, ProjectSample.PATH_MAIN_RESULT),
            "{}_{}{}".format(prefix_file_name, sequence_name, extension),
        )

    def get_file_output(self, type_path, file_type, software):
        """
        return file path output by software
        type_path: constants.TypePath -> MEDIA_ROOT, MEDIA_URL
        file_type: constants.FileType -> FILE_BAM, FILE_BAM_BAI, FILE_CONSENSUS_FA, ...
        software: SoftwareNames.SOFTWARE_FREEBAYES_name, SoftwareNames.SOFTWARE_SNIPPY_name
        """
        constants = Constants()
        return os.path.join(
            self.__get_path__(
                type_path, software.lower() if not software is None else None
            ),
            constants.get_extensions_by_file_type(self.sample.name, file_type),
        )

    def get_file_output_human(
        self, type_path, file_type, software, b_second_choice=False
    ):
        """
        return file path output by software, but with human name
        type_path: constants.TypePath -> MEDIA_ROOT, MEDIA_URL
        file_type: constants.FileType -> FILE_BAM, FILE_BAM_BAI, FILE_CONSENSUS_FA, ...
        software: SoftwareNames.SOFTWARE_FREEBAYES_name, SoftwareNames.SOFTWARE_SNIPPY_name
        """
        constants = Constants()

        return os.path.join(
            self.__get_path__(
                type_path, software.lower() if software != None else None
            ),
            constants.get_extensions_by_file_type(
                self.__get_human_name_file__(software, file_type, b_second_choice),
                file_type,
            ),
        )

    def __get_human_name_file__(self, software, file_type, b_second_choice=False):
        """
        get human file name
        """
        if (
            software == SoftwareNames.SOFTWARE_SNIPPY_name
            or software == SoftwareNames.SOFTWARE_Medaka_name
        ):
            if file_type == FileType.FILE_TAB:
                return "{}{}".format(ProjectSample.FILE_SNIPPY_TAB, self.sample.name)
        if software == SoftwareNames.SOFTWARE_FREEBAYES_name and not b_second_choice:
            if file_type == FileType.FILE_TAB:
                return "{}{}".format(ProjectSample.FILE_FREEBAYES_TAB, self.sample.name)
        if software == SoftwareNames.SOFTWARE_FREEBAYES_name and b_second_choice:
            if file_type == FileType.FILE_TAB:
                return "{}{}".format(
                    ProjectSample.FILE_FREEBAYES_TAB_with_indels, self.sample.name
                )
        return self.sample.name

    def __get_user_result_directory_path__(self, software):
        # file will be uploaded to MEDIA_ROOT/<filename>
        if software is None or not software:
            return "projects/result/user_{0}/project_{1}/sample_{2}".format(
                self.project.owner.id, self.project.id, self.sample.id
            )
        return "projects/result/user_{0}/project_{1}/sample_{2}/{3}".format(
            self.project.owner.id, self.project.id, self.sample.id, software
        )

    def __get_path__(self, type_path, software):
        """
        get a path, from MEDIA_URL or MEDIA_ROOT
        """
        path_to_find = self.__get_user_result_directory_path__(software)
        if type_path == TypePath.MEDIA_ROOT:
            if not path_to_find.startswith("/"):
                path_to_find = os.path.join(
                    getattr(settings, "MEDIA_ROOT", None), path_to_find
                )
        else:
            path_to_find = os.path.join(
                getattr(settings, "MEDIA_URL", None), path_to_find
            )
        return path_to_find

    def get_is_ready_to_proccess(self):
        """
        test if is ready to process
        """
        return (
            self.is_finished
            and not self.is_deleted
            and not self.is_error
            and self.sample.get_is_ready_for_projects()
        )

    def get_consensus_file(self, type_path):
        """
        get clean consensus file name
        """
        return os.path.join(
            self.__get_path__(type_path, ProjectSample.PATH_MAIN_RESULT),
            "{}{}{}".format(
                ProjectSample.FILE_CONSENSUS_FILE,
                self.sample.name,
                FileExtensions.FILE_FASTA,
            ),
        )

    def get_backup_consensus_file(self):
        """
        get clean consensus file name
        """
        return os.path.join(
            self.__get_path__(TypePath.MEDIA_ROOT, ProjectSample.PATH_MAIN_RESULT),
            "{}{}_backup{}".format(
                ProjectSample.FILE_CONSENSUS_FILE,
                self.sample.name,
                FileExtensions.FILE_FASTA,
            ),
        )

    def get_consensus_file_web(self, user_reject_file=False):
        """
        get consensus file web
        """
        out_file = self.get_consensus_file(TypePath.MEDIA_ROOT)
        if os.path.exists(out_file):
            if user_reject_file:
                return _("Rejected.")
            return mark_safe(
                '<a href="{}" download> {}</a>'.format(
                    self.get_consensus_file(TypePath.MEDIA_URL),
                    self.constants.short_name(os.path.basename(out_file), 15),
                )
            )
        return _("Not available.")

    def get_file_web(self, file_type, software, b_second_choice=False):
        """
        get file web from different softwares
        type_path: constants.TypePath -> MEDIA_ROOT, MEDIA_URL
        file_type: constants.FileType -> FILE_BAM, FILE_BAM_BAI, FILE_CONSENSUS_FA, ...
        software: SoftwareNames.SOFTWARE_FREEBAYES_name, SoftwareNames.SOFTWARE_SNIPPY_name
        :param b_second_choice some software can have second choices
        """
        out_file = self.get_file_output_human(
            TypePath.MEDIA_ROOT, file_type, software, b_second_choice
        )
        if os.path.exists(out_file):
            return mark_safe(
                '<a href="{}" download> {}</a>'.format(
                    self.get_file_output_human(
                        TypePath.MEDIA_URL, file_type, software, b_second_choice
                    ),
                    self.constants.short_name(os.path.basename(out_file), 20),
                )
            )
        return _("Not available.")

    def is_sample_illumina(self):
        """
        test if the sample is illumina or other
        """
        return self.sample.is_type_fastq_gz_sequencing()

    def is_sample_ont(self):
        """
        test if the sample is illumina or other
        """
        return self.sample.is_type_fastq_gz_sequencing(Sample.TYPE_OF_FASTQ_minion)

    def get_type_technology(self):
        """
        return type of technology
        """
        return self.sample.get_type_technology()

    def get_abricate_output(self, type_path, b_gzip_file=False):
        """
        type_path = [MEDIA_ROOT, MEDIA_URL]
        Return the file name of the abricate output base on fastq File input
        path it's a FileField instance, or a string
        Can be zipped
        """
        return os.path.join(
            self.__get_path__(type_path, True),
            self.OUT_FILE_ABRICATE + (FileExtensions.FILE_GZ if b_gzip_file else ""),
        )


class MetaKeyProjectSample(models.Model):
    """
    Relation ManyToMany in
    """

    meta_tag = models.ForeignKey(
        MetaKey, related_name="meta_key_project_sample", on_delete=models.CASCADE
    )
    project_sample = models.ForeignKey(
        ProjectSample, related_name="meta_key_project_sample", on_delete=models.CASCADE
    )
    owner = models.ForeignKey(
        User, related_name="meta_key_project_sample", on_delete=models.CASCADE
    )
    creation_date = models.DateTimeField(
        "uploaded date", db_index=True, auto_now_add=True
    )
    value = models.CharField(default=Constants.META_KEY_VALUE_NOT_NEED, max_length=200)
    description = models.TextField(default="")

    class Meta:
        ordering = ["project_sample__id", "-creation_date"]

    def __str__(self):
        return self.meta_tag.name + " " + self.value + " " + self.description


class Statistics(models.Model):
    """
    Has several percentils about the table CountVariations
    """

    tag = models.ForeignKey(
        TagName, related_name="statistics", on_delete=models.CASCADE
    )
    value = models.FloatField(default=0.0)

    def __str__(self):
        return "Tag: {} Value: {}".format(self.tag.name, self.value)


class Software(models.Model):
    """
    Has all the softwares
    It can be used for all users
    """

    BREAK_TAG = "$$$"
    name = models.CharField(max_length=100, blank=True, null=True)
    path_to_run = models.CharField(max_length=300)
    version = models.CharField(max_length=200)
    version_long = models.TextField(
        default=""
    )  ## has the version in json mode, in an dictonary

    ### last update of the software
    last_update = models.DateTimeField(
        auto_now_add=True, null=True, blank=True, verbose_name="Last update date"
    )

    class Meta:
        ordering = ["name"]

    def is_updated_today(self):
        return (
            not self.last_update is None
            and self.last_update.date() == datetime.now().date()
            and len(self.version_long) > 0
        )

    def set_last_update_today(self):
        self.last_update = datetime.now()

    def set_version(self, version):
        self.version = version

    def get_version(self):
        return self.version

    def set_version_long(self, dict_version):
        """
        set several versions in the version field, separated by dollar
        """
        if len(dict_version) > 0:
            self.version_long = json.dumps(dict_version)

    def get_version_long(self):
        """
        return dict ["soft name": "version", "soft1 name": "version"]
        """
        dt_result_version = {}
        if len(self.version) > 0:
            try:
                dt_result_version = json.loads(self.version_long)
            except Exception as e:
                pass
        return dt_result_version

    def is_same_version_long(self, dt_version):
        """test version as a all, all string"""
        if self.version_long == json.dumps(dt_version):
            return True
        return False


class UploadFiles(models.Model):
    """
    this class has the files that the user can upload, has he want,
    then the system make the relations with the samples
    """

    is_valid = models.BooleanField(
        default=False
    )  ## true if everything is OK with the file, without errors
    is_processed = models.BooleanField(
        default=False
    )  ## if samples file -> True when all files is attributed
    ## if fastq.gz -> True when the file is attributed
    is_deleted = models.BooleanField(default=False)  ## if this file is removed
    number_errors = models.IntegerField(
        default=0
    )  ## if has errors don't do anything, need to remove and upload again.
    number_files_processed = models.IntegerField(
        default=0
    )  ## samples_list, has the number of files already processed
    ## in fastq files, if this file is associated to a sample or not
    number_files_to_process = models.IntegerField(
        default=0
    )  ## samples_list, has the number of files to process. At the end this number must be equal to number_files_processed

    type_file = models.ForeignKey(
        MetaKey,
        related_name="upload_files",
        blank=True,
        null=True,
        on_delete=models.CASCADE,
    )  ## has the type of file
    ## constants.TYPE_FILE.TYPE_FILE_fastq_gz
    ## constants.TYPE_FILE.TYPE_FILE_sample_file
    ## constants.TYPE_FILE.TYPE_FILE_sample_file_metadata
    file_name = models.CharField(
        max_length=300, blank=True, null=True
    )  ## in fastq file, must have the same name in samples_list file
    creation_date = models.DateTimeField(
        auto_now_add=True, verbose_name="Uploaded Date"
    )
    attached_date = models.DateTimeField(
        blank=True, null=True, verbose_name="Date attached"
    )  ## only used in fastq.gz files

    ### if is deleted in file system
    is_deleted_in_file_system = models.BooleanField(
        default=False
    )  ## if this file was removed in file system
    date_deleted = models.DateTimeField(
        blank=True, null=True, verbose_name="Date attached"
    )  ## this date has the time of deleted by web page

    owner = models.ForeignKey(
        User, related_name="upload_files", on_delete=models.CASCADE
    )

    ### need to create a random name for this file
    path_name = ContentTypeRestrictedFileField(
        upload_to=user_directory_path,
        blank=True,
        null=True,
        content_types=[
            "application/octet-stream",
            "application/gzip",
            "application/x-gzip",
            "text/csv",
            "text/txt",
            "text/tsv",
            "application/vnd.ms-excel",
            "text/tab-separated-values",
            "text/plain",
        ],
        max_upload_size=(
            settings.MAX_FASTQ_FILE_WITH_DOWNSIZE
            if settings.DOWN_SIZE_FASTQ_FILES
            else settings.MAX_FASTQ_FILE_UPLOAD
        ),
        max_length=500,
    )

    samples = models.ManyToManyField(
        Sample
    )  ## if fastq file has the sample where it belongs
    ## if samples_file has all the relations with samples. Must be all created, files attributed, or deleted
    ##   to add other samples file
    upload_file = models.ForeignKey(
        "self", blank=True, null=True, on_delete=models.CASCADE
    )  ## in fastq file has the sample list where it belongs
    description = models.TextField(
        default=""
    )  ## has a json result.ProcessResults instance with errors or successes

    class Meta:
        ordering = ["-creation_date"]

    def __str__(self):
        return self.file_name

    def get_path_to_file(self, type_path):
        """
        get a path, type_path, from MEDIA_URL or MEDIA_ROOT
        """
        path_to_find = self.path_name.name
        if type_path == TypePath.MEDIA_ROOT:
            if not path_to_find.startswith("/"):
                path_to_find = os.path.join(
                    getattr(settings, "MEDIA_ROOT", None), path_to_find
                )
        else:
            path_to_find = os.path.join(
                getattr(settings, "MEDIA_URL", None), path_to_find
            )
        return path_to_find


class ProcessControler(models.Model):
    """
    this class has the info about the process
    """

    PREFIX_SAMPLE = "sample_"
    PREFIX_PROJECT_SAMPLE = "project_sample_"
    PREFIX_PROJECT = "project_"
    PREFIX_DATASET = "dataset_"
    PREFIX_TELEVIR_PROJECT = "televir_project_"
    PREFIX_TELEVIR_PROJECT_SAMPLE = "televir_project_sample"
    PREFIX_TELEVIR_REFERENCE_MAP = "televir_reference_map_"
    PREFIX_UPLOAD_FILES = "upload_files_"
    PREFIX_LINK_FILES_USER = "link_files_user_"
    PREFIX_COLLECT_ALL_SAMPLES_USER = "collect_all_samples_user_"
    PREFIX_COLLECT_ALL_PROJECTS_USER = "collect_all_projects_user_"

    ### flags of status
    FLAG_FINISHED = "flag_finished"
    FLAG_RUNNING = "flag_running"
    FLAG_ERROR = "flag_error"
    FLAG_KILLED = "flag_killed"

    owner = models.ForeignKey(
        User, related_name="process_controles", on_delete=models.CASCADE
    )
    creation_date = models.DateTimeField(auto_now_add=True, verbose_name="Submit job")
    close_date = models.DateTimeField(blank=True, null=True, verbose_name="Close date")
    is_finished = models.BooleanField(default=False)
    is_running = models.BooleanField(default=False)  ## when is running
    is_error = models.BooleanField(default=False)
    name = models.CharField(
        max_length=50, db_index=True, blank=True, null=True
    )  ## has the name of the process
    ## could be: sample_<pk>, project_sample_<pk>, project_<pk>
    name_sge_id = models.CharField(
        max_length=20, db_index=True, blank=True, null=True
    )  ## has the SGE id about this process

    class Meta:
        ordering = ["creation_date"]

    ### get names for samples, project and project samples
    def get_name_sample(self, sample):
        return "{}{}".format(ProcessControler.PREFIX_SAMPLE, sample.pk)

    def get_name_project_sample(self, project_sample):
        return "{}{}".format(ProcessControler.PREFIX_PROJECT_SAMPLE, project_sample.pk)

    def get_name_project(self, project):
        return "{}{}".format(ProcessControler.PREFIX_PROJECT, project.pk)

    def get_name_dataset(self, dataset):
        return "{}{}".format(ProcessControler.PREFIX_DATASET, dataset.pk)

    def get_name_upload_files(self, upload_files):
        return "{}{}".format(ProcessControler.PREFIX_UPLOAD_FILES, upload_files.pk)

    def get_name_link_files_user(self, user):
        return "{}{}".format(ProcessControler.PREFIX_LINK_FILES_USER, user.pk)

    def get_name_collect_all_samples_user(self, user):
        return "{}{}".format(ProcessControler.PREFIX_COLLECT_ALL_SAMPLES_USER, user.pk)

    def get_name_collect_all_projects_user(self, user):
        return "{}{}".format(ProcessControler.PREFIX_COLLECT_ALL_PROJECTS_USER, user.pk)

    def get_name_televir_project(self, project_pk):
        return "{}{}".format(ProcessControler.PREFIX_TELEVIR_PROJECT, project_pk)

    def get_name_televir_project_sample(self, project_pk, sample_pk):
        return "{}{}_sample_{}".format(
            ProcessControler.PREFIX_TELEVIR_PROJECT, project_pk, sample_pk
        )

    def get_name_televir_project_sample_metagenomics_run(self, sample_pk, leaf_pk):
        return "{}_combined_metagen_{}_{}".format(
            ProcessControler.PREFIX_TELEVIR_PROJECT, sample_pk, leaf_pk
        )
    
    def get_name_televir_project_sample_panel_map(self, sample_pk, leaf_pk):
        return "{}_televir_panel_map_{}_{}".format(
            ProcessControler.PREFIX_TELEVIR_PROJECT, sample_pk, leaf_pk
        )

    def get_name_televir_project_sample_metagenomics(self, sample_pk):
        return "{}_combined_metagen_{}".format(
            ProcessControler.PREFIX_TELEVIR_PROJECT, sample_pk
        )

    def get_name_televir_project_sample_sort(self, sample_pk):
        return "{}_report_sort_{}".format(
            ProcessControler.PREFIX_TELEVIR_PROJECT, sample_pk
        )

    def get_name_raw_televir_teleflu_ref_create(self, ref_id):
        return "{}_teleflu_ref_{}".format(
            ProcessControler.PREFIX_TELEVIR_PROJECT, ref_id
        )
    
    def get_name_file_televir_teleflu_ref_create(self, ref_id):
        return "{}_teleflu_file_ref_{}".format(
            ProcessControler.PREFIX_TELEVIR_PROJECT, ref_id
        )

    def get_name_televir_file_upload(self, file_id):
        return "televir_file_ref_{}".format(file_id)

    def get_name_televir_teleflu_igv_stack(self, teleflu_mapping_id):
        return "{}_teleflu_mapping_stack_{}".format(
            ProcessControler.PREFIX_TELEVIR_PROJECT, teleflu_mapping_id
        )

    def get_name_televir_teleflu_project_create(self, project_id):
        return "{}_teleflu_project_create_{}".format(
            ProcessControler.PREFIX_TELEVIR_PROJECT, project_id
        )

    def get_name_televir_teleflu_reference_create(self, project_id):
        return "{}_teleflu_project_reference_create_{}".format(
            ProcessControler.PREFIX_TELEVIR_PROJECT, project_id
        )

    def get_name_televir_teleflu_project_process(self, project_id):
        return "{}_teleflu_project_process_{}".format(
            ProcessControler.PREFIX_TELEVIR_PROJECT, project_id
        )

    def get_name_televir_project_merge_explify(self, project_pk):
        return "{}_report_merge_explify_{}".format(
            ProcessControler.PREFIX_TELEVIR_PROJECT, project_pk
        )

    def get_name_add_references_to_sample(self, sample_pk):
        return "{}_add_references_to_sample_{}".format(
            ProcessControler.PREFIX_TELEVIR_PROJECT, sample_pk
        )
    
    def get_name_update_televir_project(self, project_id):
        return "{}_update_project_{}".format(
            ProcessControler.PREFIX_TELEVIR_PROJECT, project_id
        )

    def get_name_televir_project_merge_explify_external(self, user_pk):
        return "{}_report_merge_explify_EXT_{}".format(
            ProcessControler.PREFIX_TELEVIR_PROJECT, user_pk
        )

    def get_name_televir_run(self, project_pk, sample_pk, leaf_pk) -> str:
        return "{}{}_sample_{}_leaf_{}".format(
            ProcessControler.PREFIX_TELEVIR_PROJECT, project_pk, sample_pk, leaf_pk
        )

    def get_name_televir_map(self, reference_pk):
        return "{}{}".format(
            ProcessControler.PREFIX_TELEVIR_REFERENCE_MAP, reference_pk
        )

    def __str__(self):
        return "PK:{} name:{}  is_finished:{}  is_running:{}  is_error:{}".format(
            self.pk, self.name, self.is_finished, self.is_running, self.is_error
        )
