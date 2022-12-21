"""
Created on 04/05/2020

@author: mmp
"""
import os

from constants.software_names import SoftwareNames
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Button, ButtonHolder, Div, Fieldset, Layout, Submit
from datasets.models import Dataset
from django import forms
from django.urls import reverse
from django.utils.html import escape
from django.utils.translation import ugettext_lazy as _
from managing_files.models import Project, ProjectSample
from pathogen_identification.utilities.utilities_pipeline import (
    Utility_Pipeline_Manager,
)
from pathogen_identification.models import Projects as TelevirProject
from utils.utils import Utils

from settings.default_parameters import DefaultParameters
from settings.models import Parameter, Sample, Software


## https://kuanyui.github.io/2015/04/13/django-crispy-inline-form-layout-with-bootstrap/
class SoftwareForm(forms.ModelForm):
    """
    Reference form, name, isolate_name and others
    """

    error_css_class = "error"
    utils = Utils()

    class Meta:
        model = Software
        fields = ()

    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop("request")
        ## try to get project or project sample, for settings update
        pk_project = kwargs.get("pk_project")
        pk_televir_project = kwargs.get("pk_televir_project")
        project = None
        televir_project = None

        ###
        self.televir_utiltity = Utility_Pipeline_Manager()
        self.televir_utiltity.get_software_list()
        self.televir_utiltity.get_software_db_dict()
        ###
        if not pk_project is None:
            kwargs.pop("pk_project")
            project = Project.objects.get(pk=pk_project)

        if not pk_televir_project is None:
            kwargs.pop("pk_televir_project")
            televir_project = TelevirProject.objects.get(pk=pk_televir_project)

        project_sample = None
        pk_project_sample = kwargs.get("pk_project_sample")
        if not pk_project_sample is None:
            kwargs.pop("pk_project_sample")
            project_sample = ProjectSample.objects.get(pk=pk_project_sample)

        sample = None
        pk_sample = kwargs.get("pk_sample")
        if not pk_sample is None:
            kwargs.pop("pk_sample")
            sample = Sample.objects.get(pk=pk_sample)

        dataset = None
        pk_dataset = kwargs.get("pk_dataset")
        if not pk_dataset is None:
            kwargs.pop("pk_dataset")
            dataset = Dataset.objects.get(pk=pk_dataset)

        ## end
        super(SoftwareForm, self).__init__(*args, **kwargs)

        ### return the parameters that is possible to change
        paramers = Parameter.objects.filter(
            software=self.instance,
            project=project,
            project_sample=project_sample,
            sample=sample,
            dataset=dataset,
            televir_project=televir_project,
        )
        dt_fields = {}
        vect_divs = []
        for parameter in paramers:
            if not parameter.can_change or parameter.is_null():
                dt_fields[parameter.get_unique_id()] = forms.CharField(
                    disabled=True,
                    required=False,
                )
                dt_fields[parameter.get_unique_id()].help_text = escape(
                    parameter.description
                )
            elif parameter.is_integer():
                dt_fields[parameter.get_unique_id()] = forms.IntegerField(
                    max_value=int(parameter.range_max),
                    min_value=int(parameter.range_min),
                    required=True,
                )
                help_text = parameter.description + " Range: {}.".format(
                    parameter.range_available
                )
                if not parameter.not_set_value is None:
                    help_text += (
                        " If value equal to {} this parameter is excluded.".format(
                            parameter.not_set_value
                        )
                    )
                dt_fields[parameter.get_unique_id()].help_text = escape(help_text)

            elif parameter.is_float():
                dt_fields[parameter.get_unique_id()] = forms.FloatField(
                    max_value=float(parameter.range_max),
                    min_value=float(parameter.range_min),
                    required=True,
                )
                help_text = parameter.description + " Range: {}.".format(
                    parameter.range_available
                )
                if not parameter.not_set_value is None:
                    help_text += (
                        " If value equal to {} this parameter is excluded.".format(
                            parameter.not_set_value
                        )
                    )
                dt_fields[parameter.get_unique_id()].help_text = escape(help_text)
            ### this is use for Medaka and Trimmomatic
            elif parameter.is_char_list():

                if parameter.software.name == SoftwareNames.SOFTWARE_NEXTSTRAIN_name:
                    list_data = [
                        data_ for data_ in SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_DESC
                    ]
                elif (
                    parameter.software.name
                    == SoftwareNames.SOFTWARE_Medaka_name_consensus
                ):
                    list_data = [
                        [data_, data_] for data_ in self.utils.get_all_medaka_models()
                    ]
                elif (
                    parameter.name == SoftwareNames.SOFTWARE_TRIMMOMATIC_illuminaclip
                    and parameter.software.name
                    == SoftwareNames.SOFTWARE_TRIMMOMATIC_name
                ):
                    list_data = [
                        [data_, data_]
                        for data_ in SoftwareNames.SOFTWARE_TRIMMOMATIC_addapter_vect_available
                    ]
                elif (
                    parameter.name == DefaultParameters.MASK_CLEAN_HUMAN_READS
                    and parameter.software.name
                    == SoftwareNames.SOFTWARE_CLEAN_HUMAN_READS_name
                ):
                    list_data = [
                        [data_, data_]
                        for data_ in SoftwareNames.SOFTWARE_CLEAN_HUMAN_READS_vect_available
                    ]
                elif (
                    parameter.name == "--db"
                    and parameter.software.pipeline_step.name
                    in self.televir_utiltity.steps_db_dependant
                ):
                    list_data = [
                        [data_, os.path.basename(data_)]
                        for data_ in self.televir_utiltity.get_from_software_db_dict(
                            parameter.software.name.lower(), []
                        )
                    ]
                else:
                    list_data = [[parameter.parameter, parameter.parameter]]

                dt_fields[parameter.get_unique_id()] = forms.ChoiceField(
                    choices=list_data
                )
                dt_fields[parameter.get_unique_id()].empty_label = None

                help_text = parameter.description
                dt_fields[parameter.get_unique_id()].help_text = escape(help_text)
            else:
                dt_fields[parameter.get_unique_id()] = forms.CharField(required=True)
                help_text = parameter.description
                if not parameter.not_set_value is None:
                    help_text += (
                        " If value equal to {} this parameter is excluded.".format(
                            parameter.not_set_value
                        )
                    )
                dt_fields[parameter.get_unique_id()].help_text = escape(help_text)

            dt_fields[parameter.get_unique_id()].label = parameter.name
            if not parameter.can_change and len(parameter.parameter) == 0:
                dt_fields[parameter.get_unique_id()].initial = parameter.name
            else:
                dt_fields[parameter.get_unique_id()].initial = parameter.parameter
            vect_divs.append(Div(parameter.get_unique_id(), css_class="col-sm-4"))

        ### set all fields
        self.fields.update(dt_fields)

        vect_rows_divs = []
        for _ in range(0, len(vect_divs), 3):
            if _ + 2 < len(vect_divs):
                vect_rows_divs.append(
                    Div(
                        vect_divs[_],
                        vect_divs[_ + 1],
                        vect_divs[_ + 2],
                        css_class="row",
                    )
                )
            elif _ + 1 < len(vect_divs):
                vect_rows_divs.append(
                    Div(vect_divs[_], vect_divs[_ + 1], css_class="row")
                )
            else:
                vect_rows_divs.append(Div(vect_divs[_], css_class="row"))

        self.helper = FormHelper()
        self.helper.form_method = "POST"

        #### message in form
        form_message = "Update {} parameters for -{}-".format(
            "software" if self.instance.is_software() else "INSaFLU", self.instance.name
        )
        if not project is None:
            form_message = "Update {} parameters for -{}- project:'{}'".format(
                "software" if self.instance.is_software() else "INSaFLU",
                self.instance.name,
                project.name,
            )
        if not dataset is None:
            form_message = "Update {} parameters for -{}- dataset:'{}'".format(
                "software" if self.instance.is_software() else "INSaFLU",
                self.instance.name,
                dataset.name,
            )
        if not project_sample is None:
            form_message = (
                "Update {} parameters for -{}- project:'{}' for sample:'{}'.".format(
                    "software" if self.instance.is_software() else "INSaFLU",
                    self.instance.name,
                    project_sample.project.name,
                    project_sample.sample.name,
                )
            )
        if not sample is None:
            form_message = "Update {} parameters for -{}- sample:'{}'.".format(
                "software" if self.instance.is_software() else "INSaFLU",
                self.instance.name,
                sample.name,
            )

        if len(vect_rows_divs) == 1:
            self.helper.layout = Layout(
                Fieldset(form_message, vect_rows_divs[0], css_class="article-content"),
                ButtonHolder(
                    Submit("save", "Save", css_class="btn-primary"),
                    Button(
                        "cancel",
                        "Cancel",
                        css_class="btn-secondary",
                        onclick='window.location.href="{}"'.format(
                            self._get_reverse(
                                project,
                                project_sample,
                                sample,
                                dataset,
                                self.instance,
                                televir_project,
                            )
                        ),
                    ),
                ),
            )
        elif len(vect_rows_divs) == 2:
            self.helper.layout = Layout(
                Fieldset(
                    form_message,
                    vect_rows_divs[0],
                    vect_rows_divs[1],
                    css_class="article-content",
                ),
                ButtonHolder(
                    Submit("save", "Save", css_class="btn-primary"),
                    Button(
                        "cancel",
                        "Cancel",
                        css_class="btn-secondary",
                        onclick='window.location.href="{}"'.format(
                            self._get_reverse(
                                project,
                                project_sample,
                                sample,
                                dataset,
                                self.instance,
                                televir_project,
                            )
                        ),
                    ),
                ),
            )
        elif len(vect_rows_divs) == 3:
            self.helper.layout = Layout(
                Fieldset(
                    form_message,
                    vect_rows_divs[0],
                    vect_rows_divs[1],
                    vect_rows_divs[2],
                    css_class="article-content",
                ),
                ButtonHolder(
                    Submit("save", "Save", css_class="btn-primary"),
                    Button(
                        "cancel",
                        "Cancel",
                        css_class="btn-secondary",
                        onclick='window.location.href="{}"'.format(
                            self._get_reverse(
                                project,
                                project_sample,
                                sample,
                                dataset,
                                self.instance,
                                televir_project,
                            )
                        ),
                    ),
                ),
            )
        elif len(vect_rows_divs) == 4:
            self.helper.layout = Layout(
                Fieldset(
                    form_message,
                    vect_rows_divs[0],
                    vect_rows_divs[1],
                    vect_rows_divs[2],
                    vect_rows_divs[3],
                    css_class="article-content",
                ),
                ButtonHolder(
                    Submit("save", "Save", css_class="btn-primary"),
                    Button(
                        "cancel",
                        "Cancel",
                        css_class="btn-secondary",
                        onclick='window.location.href="{}"'.format(
                            self._get_reverse(
                                project,
                                project_sample,
                                sample,
                                sample.dataset,
                                self.instance,
                                televir_project,
                            )
                        ),
                    ),
                ),
            )

    def _get_reverse(
        self, project, project_sample, sample, dataset, instance, televir_project
    ):
        """ """
        if instance.type_of_use == 5:
            return reverse("pathogenID_pipeline", args=[0])
        if not televir_project is None:
            return reverse("pathogenID_pipeline", args=[televir_project.pk])
        if not project is None:
            return reverse("project-settings", args=[project.pk])
        if not dataset is None:
            return reverse("dataset-settings", args=[dataset.pk])
        if not project_sample is None:
            return reverse("sample-project-settings", args=[project_sample.pk])
        if not sample is None:
            return reverse("sample-settings", args=[sample.pk])
        return reverse("settings")

    def clean(self):
        """
        Clean to check some values
        """

        cleaned_data = super(SoftwareForm, self).clean()
        #### THIS is only for NanoFilt
        if "--maxlength_5" in self.cleaned_data and "-l_2" in self.cleaned_data:
            max_length = self.cleaned_data.get("--maxlength_5")
            min_length = self.cleaned_data.get("-l_2")
            if max_length > 0 and (max_length - min_length) < 1:
                self.add_error(
                    "--maxlength_5",
                    _("Error: Max value must be bigger than Minimum value."),
                )
                self.add_error(
                    "-l_2", _("Error: Minimum value must be smaller than Max value.")
                )
                return cleaned_data
        if "--maxlength_5" in self.cleaned_data:
            max_length = self.cleaned_data.get("--maxlength_5")
            if max_length > 0 and max_length < DefaultParameters.NANOFILT_MINIMUN_MAX:
                self.add_error(
                    "--maxlength_5",
                    _(
                        "Error: Max value must be bigger than {} value.".format(
                            DefaultParameters.NANOFILT_MINIMUN_MAX
                        )
                    ),
                )
        #### END for NanoFilt

        return cleaned_data
