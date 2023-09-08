"""
Created on Sept 7, 2023

@author: js
"""

from abc import ABC
from typing import Dict

from fluwebvirus.settings import TelevirSetup


class TelevirConstantsBase(ABC):
    project_directory: str
    docker_app_directory: str
    docker_install_directory: str
    environments_directory: str
    conda_directory: str
    ref_db_directory: str
    ref_fasta_directory: str
    metadata_directory: str
    scripts_directory: str


class Televir_Directory_Constants_Docker(TelevirConstantsBase):
    """
    directory constants. To be changed on local installation without docker.
    """

    project_directory = "/tmp/televir/projects/"
    docker_app_directory = "/televir/mngs_benchmark/"
    docker_install_directory = "/televir/mngs_benchmark/mngs_environments/"
    environments_directory = "/televir/mngs_benchmark/mngs_environments/"
    conda_directory = "/opt/conda/"
    ref_db_directory = "/televir/mngs_benchmark/ref_db/"
    ref_fasta_directory = "/televir/mngs_benchmark/ref_fasta/"
    metadata_directory = "/televir/mngs_benchmark/metadata/"
    scripts_directory = "/insaflu_web/TELEVIR/deployment_scripts/scripts/"


class Televir_Directory_Constants_PreProduction(TelevirConstantsBase):
    """
    directory constants. To be changed on local installation without docker.
    """

    project_directory = "/usr/local/web_site/INSaFLU/media/televir_projects/"
    docker_app_directory = "/usr/local/web_site/televir/mngs_benchmark/"
    docker_install_directory = "/usr/local/web_site/televir/mngs_environments/"
    environments_directory = "/usr/local/web_site/televir/mngs_environments/"
    conda_directory = "/usr/local/software/insaflu/miniconda/"
    ref_db_directory = "/usr/local/web_site/televir/ref_db/"
    ref_fasta_directory = "/usr/local/web_site/televir/ref_fasta/"
    metadata_directory = "/usr/local/web_site/televir/metadata/"
    scripts_directory = "/usr/local/web_site/televir/deployment_scripts/scripts/"


class Televir_Directory_Constants_Production(TelevirConstantsBase):
    """
    directory constants. To be changed on local installation without docker.
    """

    project_directory = "/usr/local/web_site/media/televir_projects/"
    docker_app_directory = "/usr/local/web_site/televir/mngs_benchmark/"
    docker_install_directory = "/usr/local/web_site/televir/mngs_environments/"
    environments_directory = "/usr/local/web_site/televir/mngs_environments/"
    conda_directory = "/usr/local/software/insaflu/miniconda/"
    ref_db_directory = "/usr/local/web_site/televir/ref_db/"
    ref_fasta_directory = "/usr/local/web_site/televir/ref_fasta/"
    metadata_directory = "/usr/local/web_site/televir/metadata/"
    scripts_directory = "/usr/local/web_site/televir/deployment_scripts/scripts/"


def get_televir_directory_constants() -> TelevirConstantsBase:
    """
    get the directory constants
    """

    setup_map: Dict[int, TelevirConstantsBase] = {
        TelevirSetup.SETUP_DEVELOP: Televir_Directory_Constants_Docker(),
        TelevirSetup.SETUP_PREPRODUCTION: Televir_Directory_Constants_Docker(),
        TelevirSetup.SETUP_PRODUCTION: Televir_Directory_Constants_Production(),
    }

    return setup_map[TelevirSetup.CURRENT_SETUP]


Televir_Directory_Constants = get_televir_directory_constants()
