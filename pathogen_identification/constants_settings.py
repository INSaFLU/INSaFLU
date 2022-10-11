"""
Ceated on 06/05/2022
@author: joao santos
"""

from fluwebvirus.settings import MEDIA_ROOT, STATICFILES_DIRS
from settings.constants_settings import ConstantsSettings as CS


class ConstantsSettings:

    project_directory = "/tmp/televir/projects/"
    media_directory = MEDIA_ROOT
    static_directory = STATICFILES_DIRS[0]
    televir_subdirectory = "televir_projects"
    PAGINATE_NUMBER = 10
    docker_app_directory = "/televir/mngs_benchmark/"
    docker_install_directory = "/televir/mngs_benchmark/mngs_environments/"
    USER_TREE_INDEX = 0

    ###################################

    PIPELINE_STEP_ORDER = [
        CS.PIPELINE_NAME_viral_enrichment,
        CS.PIPELINE_NAME_assembly,
        CS.PIPELINE_NAME_contig_classification,
        CS.PIPELINE_NAME_read_classification,
        CS.PIPELINE_NAME_remapping,
    ]

    PIPELINE_STEPS_DB_DEPENDENT = [
        CS.PIPELINE_NAME_viral_enrichment,
        CS.PIPELINE_NAME_read_classification,
        CS.PIPELINE_NAME_contig_classification,
    ]

    ################################### Project_file_structure

    DIRS = {
        CS.PIPELINE_NAME_read_quality_analysis: "reads/clean/",
        "reads_depleted_dir": "reads/hd_filtered/",
        "reads_enriched_dir": "reads/enriched/",
        CS.PIPELINE_NAME_host_depletion: "host_depletion/",
        CS.PIPELINE_NAME_viral_enrichment: "enrichment/",
        CS.PIPELINE_NAME_assembly: "assembly/",
        CS.PIPELINE_NAME_contig_classification: "classification/assembly/",
        CS.PIPELINE_NAME_read_classification: "classification/reads/",
        CS.PIPELINE_NAME_remapping: "remap/",
        "log_dir": "logs/",
        "OUTD": "output/",
    }

    ################################## MODULES CONTROL

    ACTIONS = {
        "CLEAN": False,
        "QCONTROL": False,
        "ENRICH": True,
        "DEPLETE": False,
        "ASSEMBLE": True,
        "CLASSIFY": True,
        "REMAP": True,
        "PHAGE_DEPL": True,
        "VIRSORT": False,
        "SIFT": True,
    }

    ################################## TECHNOLOGY CONSTANTS

    CONSTANTS_ILLUMINA = {
        "minimum_coverage_threshold": 2,
        "max_output_number": 15,
        "taxid_limit": 12,
        "sift_query": "phage",
        "assembly_contig_min_length": 300,
    }

    CONSTANTS_ONT = {
        "minimum_coverage_threshold": 1,
        "max_output_number": 15,
        "taxid_limit": 12,
        "sift_query": "phage",
        "assembly_contig_min_length": 500,
    }
