"""
Ceated on 06/05/2022
@author: joao santos
"""

from fluwebvirus.settings import STATICFILES_DIRS


class ConstantsSettings:

    project_directory = "/mnt/sdc/TELEVIR/projects/"
    static_directory = STATICFILES_DIRS[0]
    static_directory_deployment = "projects"
    PAGINATE_NUMBER = 10
    docker_app_directory = "/televir/mngs_benchmark/"
    docker_install_directory = "/televir/mngs_benchmark/mngs_environments/"
    USER_TREE_INDEX = 0

    ################################### Project_file_structure

    DIRS = {
        "PREPROCESS": "reads/clean/",
        "reads_depleted_dir": "reads/hd_filtered/",
        "reads_enriched_dir": "reads/enriched/",
        "DEPLETION": "host_depletion/",
        "ENRICHMENT": "enrichment/",
        "ASSEMBLY": "assembly/",
        "CONTIG_CLASSIFICATION": "classification/assembly/",
        "READ_CLASSIFICATION": "classification/reads/",
        "REMAPPING": "remap/",
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
