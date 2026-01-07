from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union

from constants.software_names import SoftwareNames
from pathogen_identification.constants_settings import ConstantsSettings as PI_CS
from pathogen_identification.models import Projects, RunMain
from pathogen_identification.utilities.clade_objects import Clade
from pathogen_identification.utilities.mapping_flags import MappingFlagBuild
from settings.constants_settings import ConstantsSettings as CS
from settings.models import Parameter, Software


class WrongParameters(Exception):
    pass


@dataclass
class RemapParams:
    max_taxids: int
    max_accids: int
    min_coverage: int
    manual_references_include: bool = False


@dataclass
class PrinseqParams:
    entropy_threshold: float
    dust_threshold: float
    is_to_run: bool


@dataclass
class LayoutParams:
    read_overlap_threshold: float
    shared_proportion_threshold: float
    flag_str: str
    flag_build: MappingFlagBuild

    def __init__(self, read_overlap_threshold, shared_proportion_threshold, flag_str):
        self.read_overlap_threshold = read_overlap_threshold
        self.shared_proportion_threshold = shared_proportion_threshold
        self.flag_str = flag_str

    def __post_init__(self):
        self.update_flag_build(self.flag_str)

    def get_flag_build(self):
        """
        Get flag build"""
        return self.flag_build

    def update_flag_build(self, flag_str):
        """
        Update flag build based on flag_str
        """

        flag_build_list = [
            x for x in MappingFlagBuild.__subclasses__() if x.build_name == flag_str
        ]

        if len(flag_build_list) == 0:
            raise WrongParameters(
                f"Wrong flag build name {flag_str}. Available options are: {MappingFlagBuild.__subclasses__()}"
            )

        else:
            self.flag_build = flag_build_list[0]

    def __str__(self):
        return f"read_overlap_threshold: {self.read_overlap_threshold}, shared_proportion_threshold: {self.shared_proportion_threshold}, flag_str: {self.flag_str}"


class TelevirParameters:
    @staticmethod
    def technology_mincov(project: Projects):
        if project.technology in [CS.TECHNOLOGY_illumina, CS.TECHNOLOGY_illumina_old]:
            return PI_CS.CONSTANTS_ILLUMINA["minimum_coverage_threshold"]
        elif project.technology == CS.TECHNOLOGY_minion:
            return PI_CS.CONSTANTS_ONT["minimum_coverage_threshold"]
        else:
            raise Exception(f"Unknown technology {project.technology}")

    @staticmethod
    def retrieve_project_software(software_name: str, project: Projects):
        """
        Retrieve software parameters for a project
        """
        try:

            try:
                software = Software.objects.filter(
                    name=software_name,
                    owner=project.owner,
                    technology__name=project.technology,
                    type_of_use=Software.TYPE_OF_USE_televir_project_settings,
                    parameter__televir_project=project,
                ).distinct()[0]

            except IndexError:
                software = Software.objects.get(
                    name=software_name,
                    owner=project.owner,
                    technology__name=project.technology,
                    type_of_use=Software.TYPE_OF_USE_televir_settings,
                )

        except Software.DoesNotExist as exc:
            return [], None
            # raise Exception(f"Remap software not found for user {username}") from exc

        software_params = Parameter.objects.filter(
            software=software, televir_project=project
        )
        if software_params.count() == 0:
            software_params = Parameter.objects.filter(
                software=software, televir_project=None
            )
        if software_params.count() == 0:
            raise WrongParameters(
                f"Software parameters not found for user {project.owner.username} and project {project.name}"
            )

        return software_params, software

    @staticmethod
    def get_remap_software(project_pk: int):
        """
        Get remap software
        """

        project = Projects.objects.get(pk=project_pk)

        remap_params, _ = TelevirParameters.retrieve_project_software(
            SoftwareNames.SOFTWARE_REMAP_PARAMS_name, project
        )

        max_taxids = 0
        max_accids = 0
        manual_references_include = False

        for param in remap_params:
            if param.name == SoftwareNames.SOFTWARE_REMAP_PARAMS_max_taxids:
                max_taxids = int(param.parameter)
            elif param.name == SoftwareNames.SOFTWARE_REMAP_PARAMS_max_accids:
                max_accids = int(param.parameter)
            elif param.name == SoftwareNames.SOFTWARE_REMAP_PARAMS_include_manual:
                manual_references_include = param.parameter == "ON"

        remap = RemapParams(
            max_taxids=max_taxids,
            max_accids=max_accids,
            min_coverage=TelevirParameters.technology_mincov(project),
            manual_references_include=manual_references_include,
        )

        return remap

    # @staticmethod
    # def get_prinseq_software(username: str, project_name) -> PrinseqParams:
    #    """
    #    Get prinseq software
    #    """
    #    prinseq_params, prinseq_software = TelevirParameters.retrieve_project_software(
    #        SoftwareNames.SOFTWARE_PRINSEQ_name, username, project_name
    #    )
    #    entropy_threshold = 0
    #    dust_threshold = 0
    #    for param in prinseq_params:
    #        if param.name == SoftwareNames.SOFTWARE_PRINSEQ_lc_entropy:
    #            entropy_threshold = float(param.parameter)
    #        elif param.name == SoftwareNames.SOFTWARE_PRINSEQ_lc_dust:
    #            dust_threshold = float(param.parameter)
    #    prinseq = PrinseqParams(
    #        entropy_threshold=entropy_threshold,
    #        dust_threshold=dust_threshold,
    #        is_to_run=prinseq_software.is_to_run,
    #    )
    #    return prinseq

    @staticmethod
    def layout_config_get(run_params: List[Parameter]) -> LayoutParams:
        """
        Get layout parameters
        """

        report_layout_params = LayoutParams(
            PI_CS.clade_private_proportion,
            PI_CS.clade_shared_proportion_threshold,
            "viruses",
        )

        for param in run_params:
            if param.name == SoftwareNames.SOFTWARE_televir_report_layout_flag_name:
                report_layout_params.update_flag_build(param.parameter)
            elif (
                param.name
                == SoftwareNames.SOFTWARE_televir_report_layout_threshold_name
            ):
                report_layout_params.shared_proportion_threshold = float(
                    param.parameter
                )

        return report_layout_params

    @staticmethod
    def get_flag_build(project_pk) -> MappingFlagBuild:
        """
        Get flag build
        """

        # run_main = RunMain.objects.get(pk=run_pk)
        project = Projects.objects.get(pk=project_pk)
        username = project.owner.username

        (
            flag_build_params,
            _,
        ) = TelevirParameters.retrieve_project_software(
            SoftwareNames.SOFTWARE_televir_report_layout_name, project
        )

        report_layout_params = TelevirParameters.layout_config_get(flag_build_params)

        return report_layout_params.flag_build

    @staticmethod
    def get_report_layout_params(
        run_pk: Optional[int] = None, project_pk: Optional[int] = None
    ) -> LayoutParams:
        """
        Get overlap threshold
        """
        if run_pk:
            run_main = RunMain.objects.get(pk=run_pk)
            project = run_main.project
        elif project_pk:
            try:
                project = Projects.objects.get(pk=project_pk)
            except Projects.DoesNotExist:
                raise Exception("Project not found")
        else:
            raise Exception("No run or project pk")

        (
            layout_params,
            _,
        ) = TelevirParameters.retrieve_project_software(
            SoftwareNames.SOFTWARE_televir_report_layout_name, project
        )

        report_layout_config = TelevirParameters.layout_config_get(layout_params)

        return report_layout_config
