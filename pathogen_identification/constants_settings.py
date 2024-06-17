"""
Ceated on 06/05/2022
@author: joao santos
"""

import os
from typing import Dict, List

import networkx as nx

from fluwebvirus.settings import MEDIA_ROOT, STATIC_ROOT, STATICFILES_DIRS
from pathogen_identification.utilities.mapping_flags import (
    MapFlagProbes,
    MapFlagViruses,
)
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
    TELEVIR_REFERENCE_PAGINATE_NUMBER = 20

    USER_TREE_INDEX = 0

    ################################### Pipeline steps

    TEST_SOFTWARE = True

    ################################### Pipeline steps

    METAGENOMICS = True
    METAGENOMICS_file_limit = 1000000

    ################################### Pipeline_deployment_type

    DEPLOYMENT_TYPE_TREE = "tree"
    DEPLOYMENT_TYPE_PIPELINE = "pipeline"

    DEPLOYMENT_DEFAULT = DEPLOYMENT_TYPE_TREE

    ################################### Pipeline model

    PIPELINE_MODEL = 1

    ################################### Pipeline steps

    PIPELINE_STEPS_DB_DEPENDENT = [
        CS.PIPELINE_NAME_viral_enrichment,
        CS.PIPELINE_NAME_host_depletion,
        CS.PIPELINE_NAME_read_classification,
        CS.PIPELINE_NAME_contig_classification,
    ]

    #################################### Pipeline Groups

    PIPELINE_STEPS_WORKFLOWS = [
        CS.PIPELINE_NAME_extra_qc,
        CS.PIPELINE_NAME_host_depletion,
        CS.PIPELINE_NAME_viral_enrichment,
        CS.PIPELINE_NAME_assembly,
        CS.PIPELINE_NAME_contig_classification,
        CS.PIPELINE_NAME_read_classification,
        CS.PIPELINE_NAME_remapping,
        CS.PIPELINE_NAME_remap_filtering,
    ]

    PIPELINE_STEPS_MAPPINGS = [
        CS.PIPELINE_NAME_metagenomics_screening,
        CS.PIPELINE_NAME_request_mapping,
        CS.PIPELINE_NAME_map_filtering,
    ]

    PIPELINE_STEPS_GLOBAL = [
        CS.PIPELINE_NAME_reporting,
    ]

    ################################### Pipeline steps aggregate

    PIPELINE_STEPS_AGGREGATE = [
        CS.PIPELINE_NAME_remap_filtering,
        CS.PIPELINE_NAME_map_filtering,
    ]

    #################################### Tree sort default Parameters

    clade_private_proportion = 0.5
    clade_shared_proportion_threshold = 0.2

    ################################### Max Size

    MAX_LENGTH_SEQUENCE_TOTAL_REFERENCE_FASTA = 30 * 10e6

    ################################### Threads

    DEPLOYMENT_THREADS = 5
    MAPPING_THREADS = 4

    ################################### MAX time

    TIMEOUT = 60 * 60 * 5  # 5 hours
    ################################### Project_file_structure

    DIRS = {
        CS.PIPELINE_NAME_read_quality_analysis: "reads/clean/",
        CS.PIPELINE_NAME_extra_qc: "reads/clean/",
        "reads_depleted_dir": "reads/hd_filtered/",
        "reads_enriched_dir": "reads/enriched/",
        CS.PIPELINE_NAME_host_depletion: "host_depletion/",
        CS.PIPELINE_NAME_viral_enrichment: "enrichment/",
        CS.PIPELINE_NAME_assembly: "assembly/",
        CS.PIPELINE_NAME_contig_classification: "classification/assembly/",
        CS.PIPELINE_NAME_read_classification: "classification/reads/",
        CS.PIPELINE_NAME_metagenomics_screening: "metagenomics/",
        CS.PIPELINE_NAME_remapping: "remap/",
        CS.PIPELINE_NAME_request_mapping: "remap/",
        CS.PIPELINE_NAME_map_filtering: "remap/",
        CS.PIPELINE_NAME_metagenomics_screening: "metagenomics/",
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

    FLAGS_AVAILABLE = [
        MapFlagViruses,
        MapFlagProbes,
        # MapFlagBacteria,
    ]

    ################################## Description filters

    DESCRIPTION_FILTERS = ["phage"]

    ################################## TAXONOMY

    READ_OVERLAP_THRESHOLD = 0.9
    SHARED_READS_THRESHOLD = 0.5

    ################################## TECHNOLOGY CONSTANTS

    CONSTANTS_ILLUMINA = {
        "minimum_coverage_threshold": 1,
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

    ################################## constants

    PAIR_END = "PE"
    SINGLE_END = "SE"

    ################################## PATHOGEN IDENTIFICATION REPORTS

    SORT_GROUP_DISPLAY_DEFAULT_THRESHOLD_SHARED = 0.95

    ################################## ACTIONS DETAILS

    EXPLIFY_MERGE_SUFFIX = "merged_explify_project"

    @property
    def vect_pipeline_names_default(self) -> List[str]:
        vect_pipeline_names = CS.vect_pipeline_names

        return [
            pipeline_name
            for pipeline_name in vect_pipeline_names
            if pipeline_name
            not in [
                CS.PIPELINE_NAME_metagenomics_screening,
                CS.PIPELINE_NAME_request_mapping,
                CS.PIPELINE_NAME_metagenomics_settings,
            ]
        ]

    @property
    def vect_pipeline_names_condensed(self) -> Dict[str, List[str]]:
        constant_settings = CS()
        vect_pipeline_names = CS.vect_pipeline_names

        if self.METAGENOMICS is False:
            vect_pipeline_names = self.vect_pipeline_names_default

        pipeline_steps_dict = {
            pipeline_step: constant_settings.pipeline_step_to_pipeline_name(
                pipeline_step
            )
            for pipeline_step in vect_pipeline_names
        }

        pipeline_names_dict = CS.reverse_set_dict(pipeline_steps_dict)

        return pipeline_names_dict

    @property
    def vect_pipeline_groups(self) -> Dict[str, Dict[str, List[str]]]:
        constant_settings = CS()
        vect_pipeline_names = CS.vect_pipeline_names

        pipeline_steps_dict = {
            pipeline_step: constant_settings.pipeline_step_to_pipeline_name(
                pipeline_step
            )
            for pipeline_step in vect_pipeline_names
        }

        pipeline_groups_dict = {
            "Workflows": {
                p: n
                for p, n in pipeline_steps_dict.items()
                if p in self.PIPELINE_STEPS_WORKFLOWS
            },
            "Validation": {
                p: n
                for p, n in pipeline_steps_dict.items()
                if p in self.PIPELINE_STEPS_MAPPINGS
            },
            "Global": {
                p: n
                for p, n in pipeline_steps_dict.items()
                if p in self.PIPELINE_STEPS_GLOBAL
            },
        }

        pipeline_groups_dict = {
            k: constant_settings.reverse_set_dict(v)
            for k, v in pipeline_groups_dict.items()
        }

        return pipeline_groups_dict

    @property
    def vect_pipeline_televir_metagenomics_condensed(self) -> Dict[str, List[str]]:
        constant_settings = CS()

        pipeline_steps_dict = {
            pipeline_step: constant_settings.pipeline_step_to_pipeline_name(
                pipeline_step
            )
            for pipeline_step in constant_settings.vect_pipeline_televir_metagenomics_for_parameters
        }

        pipeline_names_dict = constant_settings.reverse_set_dict(pipeline_steps_dict)

        return pipeline_names_dict
