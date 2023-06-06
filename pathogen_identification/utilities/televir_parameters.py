from dataclasses import dataclass

from constants.software_names import SoftwareNames
from pathogen_identification.constants_settings import ConstantsSettings as CS
from pathogen_identification.models import Projects, RunMain
from pathogen_identification.utilities.mapping_flags import MappingFlagBuild
from settings.models import Parameter, Software


@dataclass
class RemapParams:
    max_taxids: int
    max_accids: int


@dataclass
class PrinseqParams:
    entropy_threshold: float
    dust_threshold: float
    is_to_run: bool


class TelevirParameters:
    @staticmethod
    def retrieve_project_software(software_name, username, project_name):
        """
        Retrieve software parameters for a project
        """
        try:
            project = Projects.objects.get(name=project_name)

            try:
                software = Software.objects.filter(
                    name=software_name,
                    owner__username=username,
                    technology__name=project.technology,
                    type_of_use=Software.TYPE_OF_USE_televir_project_settings,
                    parameter__televir_project=project,
                ).distinct()[0]

            except IndexError:
                software = Software.objects.get(
                    name=software_name,
                    owner__username=username,
                    technology__name=project.technology,
                    type_of_use=Software.TYPE_OF_USE_televir_settings,
                )

        except Software.DoesNotExist:
            raise Exception(f"Remap software not found for user {username}")

        software_params = Parameter.objects.filter(
            software=software, televir_project__name=project_name
        )
        if software_params.count() == 0:
            software_params = Parameter.objects.filter(
                software=software, televir_project__name=None
            )
        if software_params.count() == 0:
            raise Exception(
                f"Software parameters not found for user {username} and project {project_name}"
            )

        return software_params

    @staticmethod
    def get_remap_software(username: str, project_name):
        """
        Get remap software
        """

        remap_params = TelevirParameters.retrieve_project_software(
            SoftwareNames.SOFTWARE_REMAP_PARAMS_name, username, project_name
        )

        max_taxids = 0
        max_accids = 0
        for param in remap_params:
            if param.name == SoftwareNames.SOFTWARE_REMAP_PARAMS_max_taxids:
                max_taxids = int(param.parameter)
            elif param.name == SoftwareNames.SOFTWARE_REMAP_PARAMS_max_accids:
                max_accids = int(param.parameter)

        remap = RemapParams(max_taxids=max_taxids, max_accids=max_accids)

        return remap

    @staticmethod
    def get_prinseq_software(username: str, project_name) -> PrinseqParams:
        """
        Get prinseq software
        """

        prinseq_params = TelevirParameters.retrieve_project_software(
            SoftwareNames.SOFTWARE_PRINSEQ_name, username, project_name
        )

        entropy_threshold = 0
        dust_threshold = 0
        for param in prinseq_params:
            if param.name == SoftwareNames.SOFTWARE_PRINSEQ_lc_entropy:
                entropy_threshold = float(param.parameter)
            elif param.name == SoftwareNames.SOFTWARE_PRINSEQ_lc_dust:
                dust_threshold = float(param.parameter)

        prinseq = PrinseqParams(
            entropy_threshold=entropy_threshold,
            dust_threshold=dust_threshold,
            is_to_run=prinseq.is_to_run,
        )

        return prinseq

    @staticmethod
    def get_flag_build(run_pk) -> MappingFlagBuild:
        """
        Get flag build
        """

        run_main = RunMain.objects.get(pk=run_pk)
        project = run_main.project
        username = project.owner.username

        flag_build_params = TelevirParameters.retrieve_project_software(
            SoftwareNames.SOFTWARE_televir_report_layout_name,
            username,
            project.name,
        )

        flag_build_str = flag_build_params[0].parameter
        flag_build = [
            x
            for x in MappingFlagBuild.__subclasses__()
            if x.build_name == flag_build_str
        ]
        if len(flag_build) == 0:
            raise Exception(
                f"Flag build software parameters not found for user {username} and project {project_name}"
            )

        return flag_build[0]

    @staticmethod
    def get_read_overlap_threshold():
        """
        Get overlap threshold
        """

        return CS.READ_OVERLAP_THRESHOLD
