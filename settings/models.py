from datasets.models import Dataset
from django.contrib.auth.models import User
from django.db import models
from managing_files.models import Project, ProjectSample, Sample
from pathogen_identification.models import Projects as TelevirProject

# Create your models here.


class Technology(models.Model):

    name = models.CharField(max_length=100, db_index=True, blank=True, null=True)
    name_extended = models.CharField(
        max_length=100, db_index=True, blank=True, null=True
    )  ## extra name to show in the settings HTML table

    class Meta:
        ordering = [
            "name",
        ]

    def __str__(self):
        return self.name


class PipelineStep(models.Model):
    """
    Pipeline step
    Each software must belong to a Pipeline...
    """

    name = models.CharField(max_length=100, db_index=True, blank=True, null=True)
    name_extended = models.CharField(
        max_length=100, db_index=True, blank=True, null=True
    )  ## extra name to show in the settings HTML table

    ###  small description of this pipeline step
    help_text = models.TextField(default="")  ## show in a modal dialog

    class Meta:
        ordering = [
            "name",
        ]

    def __str__(self):
        return self.name


class Software(models.Model):
    """
    Each user has it software parameters
    """

    ### where is going to be used
    TYPE_OF_USE_global = 0  ### Used in the samples
    TYPE_OF_USE_project = 1  ### Used in a particular project
    TYPE_OF_USE_project_sample = 2  ### Used in a particular project sample
    TYPE_OF_USE_sample = 3  ### Used in a particular sample
    TYPE_OF_USE_qc = 4  ### Used for  quality control
    TYPE_OF_USE_televir_global = 5  ### used for pathogen_identification.
    TYPE_OF_USE_televir_project = 6  ### Used for  pathogen_identification_projects.
    TYPE_OF_USE_dataset = 7  ### Used in a particular dataset
    ### if it is a software parameter or a general parameter (INSaFLU parameter)
    TYPE_SOFTWARE = 0  ### normal software
    TYPE_INSAFLU_PARAMETER = 1  ### it is a general parameter (INSaFLU parameter)

    name = models.CharField(max_length=100, db_index=True, blank=True, null=True)
    name_extended = models.CharField(
        max_length=100, db_index=True, blank=True, null=True
    )  ## extra name to show in the settings HTML table
    version = models.CharField(max_length=100, db_index=True, blank=True, null=True)
    version_parameters = models.SmallIntegerField(
        default=0
    )  ### every time you add/remove parameters that can be changed by a user, add a unit here.
    type_of_use = models.SmallIntegerField(
        default=TYPE_OF_USE_global
    )  ### where is possible to define
    type_of_software = models.SmallIntegerField(
        default=TYPE_SOFTWARE
    )  ### it is a software or a general parameter
    is_obsolete = models.BooleanField(default=False)  ## set to True if it's obsolete
    owner = models.ForeignKey(
        User,
        related_name="software_settings",
        blank=True,
        null=True,
        on_delete=models.PROTECT,
    )
    technology = models.ForeignKey(
        Technology,
        related_name="software_settings",
        blank=True,
        null=True,
        on_delete=models.PROTECT,
    )

    ### if this software is to run (also if the software can be ON/OFF in pipelines)
    can_be_on_off_in_pipeline = models.BooleanField(
        default=False
    )  ## set to True if can be ON/OFF in pipeline, otherwise always ON
    ### this flag is only used on global settings. For sample, project and project_sample is on parameters...
    is_to_run = models.BooleanField(
        default=True
    )  ## set to True if it is going to run, for example Trimmomatic can run or not

    ###  small description of software
    help_text = models.TextField(default="")

    ###  which part of pipeline is going to run
    pipeline_step = models.ForeignKey(
        PipelineStep,
        related_name="software_settings",
        blank=True,
        null=True,
        on_delete=models.PROTECT,
    )

    class Meta:
        ordering = [
            "name",
        ]

    def __str__(self):
        return self.name

    def is_used_in_global(self):
        return self.type_of_use == Software.TYPE_OF_USE_global

    def is_used_in_project(self):
        return self.type_of_use == Software.TYPE_OF_USE_project

    def is_used_in_sample(self):
        return self.type_of_use == Software.TYPE_OF_USE_sample

    def is_used_in_project_sample(self):
        return self.type_of_use == Software.TYPE_OF_USE_project_sample

    def is_software(self):
        """
        :return True if it is software type
        """
        return self.type_of_software == Software.TYPE_SOFTWARE


class Parameter(models.Model):
    """
    Can have more than one, for example

    SOFTWARE_TRIMMOMATIC_PARAMETERS = "SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33"
    """

    PARAMETER_int = 0
    PARAMETER_float = 1
    PARAMETER_char = 2
    PARAMETER_null = 3  ### its like only to show, not editable, Example: (TOPHRED33)
    PARAMETER_char_list = 4  ### combo list
    PARAMETER_check_box = 5  ### check box
    PARAMETER_none = 6  ### This case is when some "software/procedure" doesn't have parameters at all
    ### It is only has one parameter. Example: "Generate consensus"
    ### "Generate consensus" -> it is used for set ON/OFF consensus in the AllConsensus File

    name = models.CharField(max_length=100, db_index=True, blank=True, null=True)
    parameter = models.CharField(max_length=150, db_index=True, blank=True, null=True)
    type_data = models.SmallIntegerField()
    software = models.ForeignKey(
        Software, related_name="parameter", on_delete=models.PROTECT
    )

    ### this allow to have software parameters in projects
    project = models.ForeignKey(
        Project,
        related_name="parameter",
        on_delete=models.PROTECT,
        blank=True,
        null=True,
    )
    televir_project = models.ForeignKey(
        TelevirProject,
        related_name="parameter",
        on_delete=models.PROTECT,
        blank=True,
        null=True,
    )
    ### this allow to have software parameters in project sample
    project_sample = models.ForeignKey(
        ProjectSample,
        related_name="parameter",
        on_delete=models.PROTECT,
        blank=True,
        null=True,
    )
    ### this allow software to use in sample
    sample = models.ForeignKey(
        Sample,
        related_name="parameter",
        on_delete=models.PROTECT,
        blank=True,
        null=True,
    )
    ### this allow to have software parameters in datasets
    dataset = models.ForeignKey(
        Dataset,
        related_name="parameter",
        on_delete=models.PROTECT,
        blank=True,
        null=True,
    )

    union_char = models.CharField(
        max_length=10, default=""
    )  ### ":" in the case LEADING:3
    can_change = models.BooleanField(
        default=True
    )  ### TOPHRED33 can change, and it is not to show in change dialog
    sequence_out = models.SmallIntegerField()
    range_available = models.CharField(
        max_length=50, default=""
    )  ## description of the range

    range_max = models.CharField(
        max_length=50, default=""
    )  ## only used in int and float fields
    range_min = models.CharField(
        max_length=50, default=""
    )  ## only used in int and float fields
    description = models.CharField(
        max_length=500, default=""
    )  ## description of this size
    not_set_value = models.CharField(
        max_length=50, db_index=True, blank=True, null=True
    )  ## don't define value, in execution string, if this value

    ## if the software it run or not for the sample, project or project_sample
    ## The first parameter is mandatory, the others don't be take into account
    is_to_run = models.BooleanField(
        default=True
    )  ## set to True if it is going to run, for example Trimmomatic can run or not

    def __str__(self):
        return self.name

    class Meta:
        ordering = [
            "sequence_out",
        ]

    def is_integer(self):
        return self.type_data == Parameter.PARAMETER_int

    def is_float(self):
        return self.type_data == Parameter.PARAMETER_float

    def is_char(self):
        return self.type_data == Parameter.PARAMETER_char

    def is_null(self):
        return self.type_data == Parameter.PARAMETER_null

    def is_char_list(self):
        return self.type_data == Parameter.PARAMETER_char_list

    def is_check_box(self):
        return self.type_data == Parameter.PARAMETER_check_box

    def get_unique_id(self):
        """
        get unique ID for forms
        """
        return "{}_{}".format(self.name, self.sequence_out)
