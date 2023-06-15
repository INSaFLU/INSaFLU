from dataclasses import dataclass

from constants.software_names import SoftwareNames
from pathogen_identification.models import Projects
from settings.models import Parameter, Software
from settings.constants_settings import ConstantsSettings as CS
from pathogen_identification.constants_settings import ConstantsSettings as PI_CS


@dataclass
class RemapParams:
    max_taxids: int
    max_accids: int
    min_quality: int
    max_mismatch: float
    min_coverage: int


@dataclass
class PrinseqParams:
    entropy_threshold: float
    dust_threshold: float
    is_to_run: bool


class TelevirParameters:
    def technology_mincov(project: Projects):
        if project.technology in [CS.TECHNOLOGY_illumina, CS.TECHNOLOGY_illumina_old]:
            return PI_CS.CONSTANTS_ILLUMINA["minimum_coverage_threshold"]
        elif project.technology == CS.TECHNOLOGY_minion:
            return PI_CS.CONSTANTS_ONT["minimum_coverage_threshold"]
        else:
            raise Exception(f"Unknown technology {project.technology}")

    @staticmethod
    def get_remap_software(username: str, project_name):
        """
        Get remap software
        """
        try:
            project = Projects.objects.get(name=project_name)

            try:
                remap = Software.objects.filter(
                    name=SoftwareNames.SOFTWARE_REMAP_PARAMS_name,
                    owner__username=username,
                    technology__name=project.technology,
                    type_of_use=Software.TYPE_OF_USE_televir_project_settings,
                    parameter__televir_project=project,
                ).distinct()[0]

            except IndexError:
                remap = Software.objects.get(
                    name=SoftwareNames.SOFTWARE_REMAP_PARAMS_name,
                    owner__username=username,
                    technology__name=project.technology,
                    type_of_use=Software.TYPE_OF_USE_televir_settings,
                )

        except Software.DoesNotExist:
            raise Exception(f"Remap software not found for user {username}")

        remap_params = Parameter.objects.filter(
            software=remap, televir_project__name=project_name
        )
        if remap_params.count() == 0:
            remap_params = Parameter.objects.filter(
                software=remap, televir_project__name=None
            )
        if remap_params.count() == 0:
            raise Exception(
                f"Remap software parameters not found for user {username} and project {project_name}"
            )
        max_taxids = 0
        max_accids = 0
        min_quality = 0
        max_mismatch = 0
        for param in remap_params:
            if param.name == SoftwareNames.SOFTWARE_REMAP_PARAMS_max_taxids:
                max_taxids = int(param.parameter)
            elif param.name == SoftwareNames.SOFTWARE_REMAP_PARAMS_max_accids:
                max_accids = int(param.parameter)
            elif param.name == SoftwareNames.SOFTWARE_REMAP_PARAMS_min_quality:
                min_quality = int(param.parameter)
            elif param.name == SoftwareNames.SOFTWARE_REMAP_PARAMS_max_mismatch:
                max_mismatch = float(param.parameter)

        remap = RemapParams(
            max_taxids=max_taxids,
            max_accids=max_accids,
            min_quality=min_quality,
            max_mismatch=max_mismatch,
            min_coverage=TelevirParameters.technology_mincov(project),
        )

        return remap

    @staticmethod
    def get_prinseq_software(username: str, project_name) -> PrinseqParams:
        """
        Get prinseq software
        """
        try:
            project = Projects.objects.get(name=project_name)

            try:
                prinseq = Software.objects.filter(
                    name=SoftwareNames.SOFTWARE_PRINSEQ_name,
                    owner__username=username,
                    technology__name=project.technology,
                    type_of_use=Software.TYPE_OF_USE_televir_project_settings,
                    parameter__televir_project=project,
                ).distinct()[0]

            except IndexError:
                prinseq = Software.objects.get(
                    name=SoftwareNames.SOFTWARE_PRINSEQ_name,
                    owner__username=username,
                    technology__name=project.technology,
                    type_of_use=Software.TYPE_OF_USE_televir_settings,
                )

        except Software.DoesNotExist:
            raise Exception(f"Prinseq software not found for user {username}")

        prinseq_params = Parameter.objects.filter(
            software=prinseq, televir_project__name=project_name
        )
        if prinseq_params.count() == 0:
            prinseq_params = Parameter.objects.filter(
                software=prinseq, televir_project__name=None
            )
        if prinseq_params.count() == 0:
            raise Exception(
                f"Prinseq software parameters not found for user {username} and project {project_name}"
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
