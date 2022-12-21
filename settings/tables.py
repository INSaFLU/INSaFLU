"""
Created on 03/05/2020

@author: mmp
"""
import logging
import os

import django_tables2 as tables
from constants.constants import Constants, TypePath
from constants.meta_key_and_values import MetaKeyAndValue
from constants.software_names import SoftwareNames
from django.urls import reverse
from django.utils.safestring import mark_safe
from extend_user.models import Profile
from managing_files.manage_database import ManageDatabase
from managing_files.models import ProjectSample
from utils.result import DecodeObjects

from settings.constants_settings import ConstantsSettings
from settings.default_parameters import DefaultParameters
from settings.default_software import DefaultSoftware
from settings.default_software_project_sample import DefaultProjectSoftware
from settings.models import Parameter, Software


class CheckBoxColumnWithName(tables.CheckBoxColumn):
    @property
    def header(self):
        default = {"type": "checkbox"}
        return self.verbose_name


class SoftwaresTable(tables.Table):

    #   Renders a normal value as an internal hyperlink to another page.
    #   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
    select_to_run = CheckBoxColumnWithName(
        verbose_name=("To Run"), accessor="pk", orderable=False
    )
    software = tables.Column("Software", orderable=False, empty_values=())
    version = tables.Column("Version", orderable=False, empty_values=())
    parameters = tables.Column("Parameters", orderable=False, empty_values=())
    options = tables.Column("Options", orderable=False, empty_values=())
    technology = tables.Column("Technology", orderable=False, empty_values=())
    constants = Constants()
    software_names = SoftwareNames()

    def __init__(
        self,
        query_set,
        project=None,
        project_sample=None,
        sample=None,
        televir_project=None,
        b_enable_options=True,
        dataset=None,
    ):
        tables.Table.__init__(self, query_set)
        self.project = project
        self.dataset = dataset
        self.project_sample = project_sample
        self.sample = sample
        self.b_enable_options = b_enable_options
        self.televir_project = televir_project

        self.count_project_sample = 0
        ### get number of samples inside of this project, if project exist
        if not project is None:
            self.count_project_sample = ProjectSample.objects.filter(
                project=project, is_deleted=False
            ).count()

    class Meta:
        model = Software()
        fields = (
            "select_to_run",
            "software",
            "technology",
            "version",
            "parameters",
            "options",
        )
        attrs = {"class": "table-striped table-bordered"}
        empty_text = "There are no Softwares to show..."

    def render_software(self, record):
        return record.name if record.name_extended is None else record.name_extended

    def render_select_to_run(self, value, record):
        ## test if its to run and get IDs from others
        is_to_run, sz_ids = self._is_to_run(record)
        ### When in sample you can not turn ON|OFF the software
        b_enable_options = self.b_enable_options
        if not self.sample is None:
            b_enable_options = True

        ### check if can ON/OFF SOFTWARE_GENERATE_CONSENSUS_name (consensus make exist)
        if (
            not self.project_sample is None
            and record.name == SoftwareNames.SOFTWARE_GENERATE_CONSENSUS_name
        ):
            b_enable_options = os.path.exists(
                self.project_sample.get_consensus_file(TypePath.MEDIA_ROOT)
            )

        ## need to remove # in href, otherwise still active
        sz_href = (
            '<a href="{}id_turn_software_on_off" {} '.format(
                "#" if record.can_be_on_off_in_pipeline and b_enable_options else "",
                'id="id_show_turn_on_off_modal"'
                if record.can_be_on_off_in_pipeline and b_enable_options
                else "",
            )
            + 'data-toggle="modal" type_of_use_id="{}" {} software_id="{}" {}'.format(
                record.type_of_use,
                f'televir_project_id="{self.televir_project.id}"'
                if self.televir_project
                else "",
                record.id,
                sz_ids,
            )
            + '><input name="select_to_run" id="{}_{}" type="checkbox" value="{}" {} {}/> </a>'.format(
                Constants.CHECK_BOX,
                record.id,
                record.id,
                "checked" if is_to_run else "",
                ""
                if record.can_be_on_off_in_pipeline and b_enable_options
                else "disabled",
            )
        )
        return mark_safe(sz_href)

    def render_parameters(self, **kwargs):
        """
        render parameters for the software
        """
        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()
        user = current_request.user

        record = kwargs.pop("record")
        technology_name = (
            ConstantsSettings.TECHNOLOGY_illumina
            if record.technology is None
            else record.technology.name
        )
        print(record.name)
        if (
            self.project is None
            and self.project_sample is None
            and self.sample is None
            and self.dataset is None
            and self.televir_project is None
        ):
            default_software = DefaultSoftware()
            return default_software.get_parameters(record.name, user, technology_name)
        elif self.dataset is not None:

            # default_software_projects = DefaultProjectSoftware()
            # logger.debug("Dataset parameters for {}".format(self.dataset))
            # TODO need to make this work...
            # parameters = default_software_projects.get_parameters(software_name=record.name, user=user,
            #    type_of_use=None, project=None, project_sample=None, sample=None,
            #    technology_name=ConstantsSettings.TECHNOLOGY_generic, dataset=self.dataset)
            # logger.debug("Dataset parameters:{}".format(parameters))
            # return parameters
            parameters_list = Parameter.objects.filter(dataset__name=self.dataset)

            if len(list(parameters_list)) == 1:
                return list(parameters_list)[0].parameter
            else:
                return "None"
        elif not self.project is None:
            default_software_projects = DefaultProjectSoftware()
            parameters = default_software_projects.get_parameters(
                record.name,
                user,
                Software.TYPE_OF_USE_project,
                self.project,
                None,
                None,
                technology_name,
            )
            if parameters == DefaultParameters.MASK_DONT_care:
                if record.name == SoftwareNames.SOFTWARE_MASK_CONSENSUS_BY_SITE_name:
                    return self.get_info_about_mask_consensus_by_sites()
            return parameters
        elif not self.project_sample is None:
            default_software_projects = DefaultProjectSoftware()
            parameters = default_software_projects.get_parameters(
                record.name,
                user,
                Software.TYPE_OF_USE_project_sample,
                None,
                self.project_sample,
                None,
                technology_name,
            )
            if parameters == DefaultParameters.MASK_DONT_care:
                if record.name == SoftwareNames.SOFTWARE_MASK_CONSENSUS_BY_SITE_name:
                    return self.get_info_about_mask_consensus_by_sites()
            return parameters
        elif not self.sample is None:
            default_software_projects = DefaultProjectSoftware()
            return default_software_projects.get_parameters(
                record.name,
                user,
                Software.TYPE_OF_USE_sample,
                None,
                None,
                self.sample,
                technology_name,
            )
        elif not self.televir_project is None:
            print(self.televir_project)
            default_software_projects = DefaultProjectSoftware()
            return default_software_projects.get_parameters(
                record.name,
                user,
                Software.TYPE_OF_USE_televir_project,
                None,
                None,
                None,
                technology_name,
                televir_project=self.televir_project,
            )
        return ""

    def render_technology(self, record):
        """return technology names"""
        return (
            ConstantsSettings.TECHNOLOGY_illumina
            if record.technology is None
            else record.technology.name
        )

    def render_options(self, record):

        ### if project
        ## Edit
        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()
        try:
            profile = Profile.objects.get(user=current_request.user)
            if profile.only_view_project:
                return "Options not available"
        except Profile.DoesNotExist:
            pass

        ### define tooltip message
        tooltip_reset = "Set default parameters"
        if record.name == SoftwareNames.SOFTWARE_MASK_CONSENSUS_BY_SITE_name:
            tooltip_reset = "Clean all positions"

        ## start define links
        b_enable_options = self.b_enable_options
        if not self.televir_project is None:  ## for televir projects
            if b_enable_options:  ## If b_enable_options is False it's false to All
                default_software_projects = DefaultProjectSoftware()
                b_enable_options = (
                    default_software_projects.can_change_values_for_this_software(
                        record,
                        project=None,
                        project_sample=None,
                        sample=None,
                        dataset=self.dataset,
                    )
                )

            str_links = (
                "<a href="
                + reverse(
                    "software-televir-project-update",
                    args=[record.pk, self.televir_project.pk],
                )
                + ' data-toggle="tooltip" title="Edit parameters" '
                + "{}".format(
                    "" if b_enable_options else "onclick='return false;' disable"
                )
                + '><span><i class="fa fa-2x fa-pencil padding-button-table '
                + "{}".format("" if b_enable_options else "disable_fa_icon")
                + '"></i></span></a>'
            )

            str_links += (
                '<a href="{}"'.format(
                    "#id_set_default_modal" if b_enable_options else ""
                )
                + ' id="id_default_parameter" data-toggle="modal" data-toggle="tooltip" title="{}"'.format(
                    tooltip_reset
                )
                + "{}".format(
                    "" if b_enable_options else 'onclick="return false;" disable'
                )
                + ' ref_name="'
                + record.name
                + '" pk="'
                + str(record.pk)
                + '" pk_televir_project="'
                + str(self.televir_project.pk)
                + '" type_software="{}'.format(
                    "software" if record.is_software() else "INSaFLU"
                )
                + '" televir_project="'
                + str(self.televir_project.name)
                + '"><span ><i class="fa fa-2x fa-power-off padding-button-table '
                + "{}".format("" if b_enable_options else "disable_fa_icon")
                + '"></i></span></a>'
            )
        elif not self.project is None:  ## for projects
            ### turn on/off buttons
            if (
                record.name == SoftwareNames.SOFTWARE_MASK_CONSENSUS_BY_SITE_name
            ):  ### allways true for this software
                b_enable_options = True
            elif b_enable_options:  ## If b_enable_options is False it's false to All
                default_software_projects = DefaultProjectSoftware()
                b_enable_options = (
                    default_software_projects.can_change_values_for_this_software(
                        record, self.project, None, None
                    )
                )

            if record.name == SoftwareNames.SOFTWARE_MASK_CONSENSUS_BY_SITE_name:
                (
                    has_data,
                    toolpit_masking,
                ) = self.get_tooltip_about_mask_consensus_by_sites()
                str_links = (
                    '<a href="#id_set_positions_to_mask_regions" id="showMaskModal" data-toggle="modal" project_id="{}" '.format(
                        self.project.id
                    )
                    + 'reference_name="{}" id_image=icon_mask_consensus_{} project_name="{}" >'.format(
                        self.project.reference.name, self.project.pk, self.project.name
                    )
                    + '<span><i id=icon_mask_consensus_{} class="padding-button-table {} fa fa-2x fa-pencil padding-button-table tip" title="{}"></i></span></a>'.format(
                        self.project.pk,
                        "warning_fa_icon" if has_data else "",
                        toolpit_masking,
                    )
                )
            else:
                str_links = (
                    "<a href="
                    + reverse(
                        "software-project-update", args=[record.pk, self.project.pk]
                    )
                    + ' data-toggle="tooltip" title="Edit parameters" '
                    + "{}".format(
                        "" if b_enable_options else "onclick='return false;' disable"
                    )
                    + '><span><i class="fa fa-2x fa-pencil padding-button-table '
                    + "{}".format("" if b_enable_options else "disable_fa_icon")
                    + '"></i></span></a>'
                )

            str_links += (
                '<a href="{}"'.format(
                    "#id_set_default_modal" if b_enable_options else ""
                )
                + ' id="id_default_parameter" data-toggle="modal" data-toggle="tooltip" title="{}"'.format(
                    tooltip_reset
                )
                + "{}".format(
                    "" if b_enable_options else 'onclick="return false;" disable'
                )
                + ' ref_name="'
                + record.name
                + '" pk="'
                + str(record.pk)
                + '" pk_proj="'
                + str(self.project.pk)
                + '" type_software="{}'.format(
                    "software" if record.is_software() else "INSaFLU"
                )
                + '" proj_name="'
                + str(self.project.name)
                + '"><span ><i class="fa fa-2x fa-power-off padding-button-table '
                + "{}".format("" if b_enable_options else "disable_fa_icon")
                + '"></i></span></a>'
            )
        elif not self.dataset is None:  ## for datasets
            ### turn on/off buttons
            if b_enable_options:  ## If b_enable_options is False it's false to All
                default_software_projects = DefaultProjectSoftware()
                b_enable_options = (
                    default_software_projects.can_change_values_for_this_software(
                        record,
                        project=None,
                        project_sample=None,
                        sample=None,
                        dataset=self.dataset,
                    )
                )

            str_links = (
                "<a href="
                + reverse("software-dataset-update", args=[record.pk, self.dataset.pk])
                + ' data-toggle="tooltip" title="Edit parameters" '
                + "{}".format(
                    "" if b_enable_options else "onclick='return false;' disable"
                )
                + '><span><i class="fa fa-2x fa-pencil padding-button-table '
                + "{}".format("" if b_enable_options else "disable_fa_icon")
                + '"></i></span></a>'
            )

            str_links += (
                '<a href="{}"'.format(
                    "#id_set_default_modal" if b_enable_options else ""
                )
                + ' id="id_default_parameter" data-toggle="modal" data-toggle="tooltip" title="{}"'.format(
                    tooltip_reset
                )
                + "{}".format(
                    "" if b_enable_options else 'onclick="return false;" disable'
                )
                + ' ref_name="'
                + record.name
                + '" pk="'
                + str(record.pk)
                + '" pk_dataset="'
                + str(self.dataset.pk)
                + '" type_software="{}'.format(
                    "software" if record.is_software() else "INSaFLU"
                )
                + '" dataset_name="'
                + str(self.dataset.name)
                + '"><span ><i class="fa fa-2x fa-power-off padding-button-table '
                + "{}".format("" if b_enable_options else "disable_fa_icon")
                + '"></i></span></a>'
            )
        elif not self.project_sample is None:  ## for project samples
            if b_enable_options:  ## If b_enable_options is False it's false to All
                default_software_projects = DefaultProjectSoftware()
                b_enable_options = (
                    default_software_projects.can_change_values_for_this_software(
                        record, None, self.project_sample, None
                    )
                )

            if record.name == SoftwareNames.SOFTWARE_MASK_CONSENSUS_BY_SITE_name:
                (
                    has_data,
                    toolpit_masking,
                ) = self.get_tooltip_about_mask_consensus_by_sites()
                str_links = (
                    '<a href="#id_set_positions_to_mask_regions" id="showMaskModal" data-toggle="modal" project_sample_id="{}" '.format(
                        self.project_sample.id
                    )
                    + 'reference_name="{}" id_image=icon_mask_consensus_{} project_name="{}" >'.format(
                        self.project_sample.project.reference.name,
                        self.project_sample.pk,
                        self.project_sample.project.name,
                    )
                    + '<span><i id=icon_mask_consensus_{} class="padding-button-table {} fa fa-2x fa-pencil padding-button-table tip" title="{}"></i></span></a>'.format(
                        self.project_sample.pk,
                        "warning_fa_icon" if has_data else "",
                        toolpit_masking,
                    )
                )
                b_enable_options = (
                    True  ## because for MaskConsensusBySite is always active
                )
            else:
                str_links = (
                    "<a href="
                    + reverse(
                        "software-project-sample-update",
                        args=[record.pk, self.project_sample.pk],
                    )
                    + ' data-toggle="tooltip" title="Edit parameters" '
                    + "{}".format(
                        "" if b_enable_options else "onclick='return false;' disable"
                    )
                    + '><span><i class="fa fa-2x fa-pencil padding-button-table '
                    + "{}".format("" if b_enable_options else "disable_fa_icon")
                    + '"></i></span></a>'
                )

            str_links += (
                '<a href="{}"'.format(
                    "#id_set_default_modal" if b_enable_options else ""
                )
                + ' id="id_default_parameter" data-toggle="modal" data-toggle="tooltip" title="{}"'.format(
                    tooltip_reset
                )
                + "{}".format(
                    "" if b_enable_options else "onclick='return false;' disable"
                )
                + ' ref_name="'
                + record.name
                + '" pk="'
                + str(record.pk)
                + '" pk_proj_sample="'
                + str(self.project_sample.pk)
                + '" type_software="{}'.format(
                    "software" if record.is_software() else "INSaFLU"
                )
                + '" proj_name="'
                + str(self.project_sample.project.name)
                + '"><span ><i class="fa fa-2x fa-power-off padding-button-table '
                + "{}".format("" if b_enable_options else "disable_fa_icon")
                + '"></i></span></a>'
            )
        elif not self.sample is None:  ## for samples
            is_to_run, sz_ids = self._is_to_run(record)
            if (
                b_enable_options and is_to_run
            ):  ## If b_enable_options is False it's false to All
                default_software_projects = DefaultProjectSoftware()
                b_enable_options = (
                    default_software_projects.can_change_values_for_this_software(
                        record, None, None, self.sample
                    )
                )
            else:
                b_enable_options = False
            str_links = (
                "<a href="
                + reverse("software-sample-update", args=[record.pk, self.sample.pk])
                + ' data-toggle="tooltip" title="Edit parameters" '
                + "{}".format(
                    "" if b_enable_options else "onclick='return false;' disable"
                )
                + '><span><i class="fa fa-2x fa-pencil padding-button-table '
                + "{}".format("" if b_enable_options else "disable_fa_icon")
                + '"></i></span></a>'
            )
            str_links += (
                '<a href="{}"'.format(
                    "#id_set_default_modal" if b_enable_options else ""
                )
                + ' id="id_default_parameter" data-toggle="modal" data-toggle="tooltip" title="{}"'.format(
                    tooltip_reset
                )
                + "{}".format(
                    "" if b_enable_options else "onclick='return false;' disable"
                )
                + ' ref_name="'
                + record.name
                + '" pk="'
                + str(record.pk)
                + '" pk_sample="'
                + str(self.sample.pk)
                + '" type_software="{}'.format(
                    "software" if record.is_software() else "INSaFLU"
                )
                + '" proj_name="'
                + str(self.sample.name)
                + '"><span ><i class="fa fa-2x fa-power-off padding-button-table '
                + "{}".format("" if b_enable_options else "disable_fa_icon")
                + '"></i></span></a>'
            )
        else:  ## for all
            ### test if all parameters are ON/OFF
            if b_enable_options:  ## If b_enable_options is False it's false to All
                default_software_projects = DefaultProjectSoftware()
                b_enable_options = (
                    default_software_projects.can_change_values_for_this_software(
                        record, None, None, None
                    )
                )
            str_links = (
                "<a href="
                + reverse("software-update", args=[record.pk])
                + ' data-toggle="tooltip" title="Edit parameters" '
                + "{}".format(
                    "" if b_enable_options else 'onclick="return false;" disable'
                )
                + '><span><i class="fa fa-2x fa-pencil padding-button-table '
                + "{}".format("" if b_enable_options else "disable_fa_icon")
                + '"></i></span></a>'
            )
            ## set default values
            str_links += (
                '<a href="{}"'.format(
                    "#id_set_default_modal" if b_enable_options else ""
                )
                + ' id="id_default_parameter" data-toggle="modal" data-toggle="tooltip" title="{}"'.format(
                    tooltip_reset
                )
                + ' ref_name="'
                + record.name
                + '"'
                "{}".format("" if b_enable_options else "onclick='return false;'")
                + ' type_software="{}'.format(
                    "software" if record.is_software() else "INSaFLU"
                )
                + '" pk="'
                + str(record.pk)
                + '"><span ><i class="fa fa-2x fa-power-off padding-button-table '
                + "{}".format("" if b_enable_options else "disable_fa_icon")
                + '"></i></span></a>'
            )
        return mark_safe(str_links)

    def has_change_parameters(self, record):
        """test if has some parameters that can be changed"""
        for parameter in Parameter.objects.filter(
            software=record,
            project=self.project,
            dataset=self.dataset,
            project_sample=self.project_sample,
            sample=self.sample,
        ):
            if parameter.can_change:
                return True
        return False

    def _is_to_run(self, record):
        """test if a software is to run and return the ids"""
        sz_ids = ""
        if not self.project is None:
            sz_ids += 'project_id="{}"'.format(self.project.id)
        if not self.dataset is None:
            sz_ids += ' dataset_id="{}"'.format(self.dataset.id)
        if not self.project_sample is None:
            sz_ids += ' project_sample_id="{}"'.format(self.project_sample.id)
        if not self.sample is None:
            sz_ids += ' sample_id="{}"'.format(self.sample.id)

        is_to_run = record.is_to_run
        if len(sz_ids) > 0:
            parameters = Parameter.objects.filter(
                software=record,
                project=self.project,
                dataset=self.dataset,
                project_sample=self.project_sample,
                sample=self.sample,
                televir_project=self.televir_project,
            )
            if len(parameters) > 0:
                is_to_run = parameters[0].is_to_run
        return is_to_run, sz_ids

    def get_info_about_mask_consensus_by_sites(self):
        """get info"""

        manageDatabase = ManageDatabase()
        decode_masking_consensus = DecodeObjects()
        meta_value = None
        if not self.project is None:
            meta_value = manageDatabase.get_project_metakey_last(
                self.project,
                MetaKeyAndValue.META_KEY_Masking_consensus,
                MetaKeyAndValue.META_VALUE_Success,
            )
        if not self.project_sample is None:
            meta_value = manageDatabase.get_project_sample_metakey_last(
                self.project_sample,
                MetaKeyAndValue.META_KEY_Masking_consensus,
                MetaKeyAndValue.META_VALUE_Success,
            )
            if meta_value is None:
                meta_value = manageDatabase.get_project_metakey_last(
                    self.project_sample.project,
                    MetaKeyAndValue.META_KEY_Masking_consensus,
                    MetaKeyAndValue.META_VALUE_Success,
                )
        if not meta_value is None:
            masking_consensus_original = decode_masking_consensus.decode_result(
                meta_value.description
            )
            if masking_consensus_original.has_masking_data():
                return "Has positions masked"
        return "Values not set, yet"

    def get_tooltip_about_mask_consensus_by_sites(self):
        """get info"""

        manageDatabase = ManageDatabase()
        decode_masking_consensus = DecodeObjects()
        meta_value = None
        if not self.project is None:
            meta_value = manageDatabase.get_project_metakey_last(
                self.project,
                MetaKeyAndValue.META_KEY_Masking_consensus,
                MetaKeyAndValue.META_VALUE_Success,
            )
        if not self.project_sample is None:
            meta_value = manageDatabase.get_project_sample_metakey_last(
                self.project_sample,
                MetaKeyAndValue.META_KEY_Masking_consensus,
                MetaKeyAndValue.META_VALUE_Success,
            )
            if meta_value is None:
                meta_value = manageDatabase.get_project_metakey_last(
                    self.project_sample.project,
                    MetaKeyAndValue.META_KEY_Masking_consensus,
                    MetaKeyAndValue.META_VALUE_Success,
                )
        if not meta_value is None:
            masking_consensus_original = decode_masking_consensus.decode_result(
                meta_value.description
            )
            if masking_consensus_original.has_masking_data():
                return (
                    True,
                    masking_consensus_original.get_message_mask_to_show_in_web_site(),
                )
        return (False, "Edit parameters")


class INSaFLUParametersTable(tables.Table):

    #   Renders a normal value as an internal hyperlink to another page.
    #   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
    select_to_run = CheckBoxColumnWithName(
        verbose_name=("To Run"), accessor="pk", orderable=False
    )
    software = tables.Column("Software", orderable=False, empty_values=())
    parameters = tables.Column("Parameters", orderable=False, empty_values=())
    options = tables.Column("Options", orderable=False, empty_values=())
    technology = tables.Column("Technology", orderable=False, empty_values=())
    constants = Constants()

    def __init__(
        self, query_set, project=None, project_sample=None, b_enable_options=True
    ):
        tables.Table.__init__(self, query_set)
        self.project = project
        self.project_sample = project_sample
        self.b_enable_options = b_enable_options

        self.count_project_sample = 0
        ### get number of samples inside of this project, if project exist
        if not project is None:
            self.count_project_sample = ProjectSample.objects.filter(
                project=project, is_deleted=False
            ).count()

    class Meta:
        model = Software()
        fields = ("select_to_run", "software", "technology", "parameters", "options")
        attrs = {"class": "table-striped table-bordered"}
        empty_text = "There are no INSaFLU parameters to show..."

    def render_software(self, record):
        return record.name if record.name_extended is None else record.name_extended

    def render_select_to_run(self, value, record):

        sz_ids = ""
        if not self.project is None:
            sz_ids += 'project_id="{}"'.format(self.project)
        # if (not self.project is None): sz_ids += 'project_id="{}"'.format(self.project)
        if not self.project_sample is None:
            sz_ids += ' project_sample_id="{}"'.format(self.project_sample)

        return mark_safe(
            '<input name="select_to_run" id="{}_{}" type="checkbox" value="{}" {} {} {}/>'.format(
                Constants.CHECK_BOX,
                record.id,
                record.id,
                "checked" if record.is_to_run else "",
                "" if record.can_be_on_off_in_pipeline else "disabled",
                sz_ids,
            )
        )

    def render_technology(self, record):
        """return technology names"""
        return (
            ConstantsSettings.TECHNOLOGY_illumina
            if record.technology is None
            else record.technology.name
        )

    def render_parameters(self, **kwargs):
        """
        render parameters for the software
        """
        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()
        user = current_request.user

        record = kwargs.pop("record")
        technology_name = (
            ConstantsSettings.TECHNOLOGY_illumina
            if record.technology is None
            else record.technology.name
        )
        if self.project is None and self.project_sample is None:
            default_software = DefaultSoftware()
            return default_software.get_parameters(record.name, user, technology_name)
        elif self.project_sample is None:
            default_software_projects = DefaultProjectSoftware()
            return default_software_projects.get_parameters(
                record.name,
                user,
                Software.TYPE_OF_USE_project,
                self.project,
                None,
                None,
                technology_name,
            )
        elif self.project is None:
            default_software_projects = DefaultProjectSoftware()
            return default_software_projects.get_parameters(
                record.name,
                user,
                Software.TYPE_OF_USE_project_sample,
                None,
                self.project_sample,
                None,
                technology_name,
            )
        return ""

    def render_options(self, record):

        ### if project
        ## Edit
        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()
        try:
            profile = Profile.objects.get(user=current_request.user)
            if profile.only_view_project:
                return "Options not available"
        except Profile.DoesNotExist:
            pass

        if not self.project is None:
            str_links = (
                "<a href="
                + reverse("software-project-update", args=[record.pk, self.project.pk])
                + ' data-toggle="tooltip" title="Edit parameters" '
                + "{}".format(
                    "" if self.b_enable_options else "onclick='return false;' disable"
                )
                + '><span><i class="fa fa-2x fa-pencil padding-button-table '
                + "{}".format("" if self.b_enable_options else "disable_fa_icon")
                + '"></i></span></a>'
            )
            str_links += (
                '<a href="{}"'.format(
                    "#id_set_default_modal" if self.b_enable_options else ""
                )
                + ' id="id_default_parameter" data-toggle="modal" data-toggle="tooltip" title="Set default parameters"'
                + "{}".format(
                    "" if self.b_enable_options else "onclick='return false;' disable"
                )
                + ' ref_name="'
                + record.name
                + '" pk="'
                + str(record.pk)
                + '" pk_proj="'
                + str(self.project.pk)
                + '" type_software="{}'.format(
                    "software" if record.is_software() else "INSaFLU"
                )
                + '" proj_name="'
                + str(self.project.name)
                + '"><span ><i class="fa fa-2x fa-power-off padding-button-table '
                + "{}".format("" if self.b_enable_options else "disable_fa_icon")
                + '"></i></span></a>'
            )
        elif not self.project_sample is None:
            str_links = (
                "<a href="
                + reverse(
                    "software-project-sample-update",
                    args=[record.pk, self.project_sample.pk],
                )
                + ' data-toggle="tooltip" title="Edit parameters" '
                + "{}".format(
                    "" if self.b_enable_options else "onclick='return false;' disable"
                )
                + '><span><i class="fa fa-2x fa-pencil padding-button-table '
                + "{}".format("" if self.b_enable_options else "disable_fa_icon")
                + '"></i></span></a>'
            )

            str_links += (
                '<a href="{}"'.format(
                    "#id_set_default_modal" if self.b_enable_options else ""
                )
                + ' id="id_default_parameter" data-toggle="modal" data-toggle="tooltip" title="Set default parameters"'
                + "{}".format(
                    "" if self.b_enable_options else "onclick='return false;' disable"
                )
                + ' ref_name="'
                + record.name
                + '" pk="'
                + str(record.pk)
                + '" pk_proj_sample="'
                + str(self.project_sample.pk)
                + '" type_software="{}'.format(
                    "software" if record.is_software() else "INSaFLU"
                )
                + '" proj_name="'
                + str(self.project_sample.project.name)
                + '"><span ><i class="fa fa-2x fa-power-off padding-button-table '
                + "{}".format("" if self.b_enable_options else "disable_fa_icon")
                + '"></i></span></a>'
            )
        else:
            str_links = (
                "<a href="
                + reverse("software-update", args=[record.pk])
                + ' data-toggle="tooltip" title="Edit parameters" '
                + '><span><i class="fa fa-2x fa-pencil padding-button-table '
                + '"></i></span></a>'
            )
            ## Remove
            str_links += (
                '<a href="{}"'.format(
                    "#id_set_default_modal" if self.b_enable_options else ""
                )
                + ' id="id_default_parameter" data-toggle="modal" data-toggle="tooltip" title="Set default parameters"'
                + ' ref_name="'
                + record.name
                + '" type_software="{}'.format(
                    "software" if record.is_software() else "INSaFLU"
                )
                + '" pk="'
                + str(record.pk)
                + '"><span ><i class="fa fa-2x fa-power-off padding-button-table '
                + '"></i></span></a>'
            )
        return mark_safe(str_links)
