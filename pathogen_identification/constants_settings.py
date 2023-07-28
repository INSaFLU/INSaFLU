"""
Ceated on 06/05/2022
@author: joao santos
"""

import os

import networkx as nx

from fluwebvirus.settings import MEDIA_ROOT, STATIC_ROOT, STATICFILES_DIRS
from pathogen_identification.utilities.mapping_flags import MapFlagViruses
from settings.constants_settings import ConstantsSettings as CS


class ConstantsSettings:
    """
    Constants for pathogen identification
    """

    media_directory = MEDIA_ROOT
    static_directory = STATIC_ROOT
    televir_subdirectory = "televir_projects"
    run_files_zipped = "run.zip"
    PAGINATE_NUMBER = 10

    USER_TREE_INDEX = 0

    ################################### Pipeline steps

    PIPELINE_STEPS_DB_DEPENDENT = [
        CS.PIPELINE_NAME_viral_enrichment,
        CS.PIPELINE_NAME_host_depletion,
        CS.PIPELINE_NAME_read_classification,
        CS.PIPELINE_NAME_contig_classification,
    ] 

    ################################### pipeline_deployment_type

    DEPLOYMENT_TYPE_TREE = "tree"
    DEPLOYMENT_TYPE_PIPELINE = "pipeline"

    DEPLOYMENT_DEFAULT = DEPLOYMENT_TYPE_PIPELINE

    ################################### Threads

    DEPLOYMENT_THREADS = 5
    MAPPING_THREADS = 4

    ################################### MAX time

    TIMEOUT = 60 * 60 * 5  # 5 hours
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

    ################################## FLAG BUILDS
    FLAG_BUILD_DEFAULT = MapFlagViruses

    ################################## Description filters

    DESCRIPTION_FILTERS = ["phage"]

    ################################## TAXONOMY

    READ_OVERLAP_THRESHOLD = 0.9
    SHARED_READS_THRESHOLD = 0.5

    ################################## TECHNOLOGY CONSTANTS

    CONSTANTS_ILLUMINA = {
        "minimum_coverage_threshold": 2,
        "max_output_number": 1,
        "taxid_limit": 15,
        "sift_query": "phage",
        "assembly_contig_min_length": 300,
    }

    CONSTANTS_ONT = {
        "minimum_coverage_threshold": 1,
        "max_output_number": 1,
        "taxid_limit": 15,
        "sift_query": "phage",
        "assembly_contig_min_length": 500,
    }

    ################################## SOFTWARE
