"""
Created on Dec 14, 2017

@author: mmp
"""

from re import S
from typing import Dict, List


class ConstantsSettings(object):
    """
    classdocs
    """

    ## setup
    SETUP_DEVELOP = 0
    SETUP_PREPRODUCTION = 1
    SETUP_PRODUCTION = 2

    CURRENT_SETUP = SETUP_DEVELOP
    ## constants

    PIPELINE_NAME_read_quality_analysis = "Read quality analysis and improvement"
    PIPELINE_NAME_type_and_subtype_analysis = "Classification"
    PIPELINE_NAME_variant_detection = "Mutation detection and consensus generation"
    PIPELINE_NAME_coverage_analysis = "Coverage analysis"
    PIPELINE_NAME_alignment = "Alignment/Phylogeny"
    PIPELINE_NAME_intra_host_minor_variant_detection = (
        "Intra-host minor variant detection"
    )
    PIPELINE_NAME_extra_qc = "Extra QC"
    PIPELINE_NAME_viral_enrichment = "Viral enrichment"
    PIPELINE_NAME_enrichment = "Enrichment"
    PIPELINE_NAME_host_depletion = "Host depletion"
    PIPELINE_NAME_contig_classification = "Contig classification"
    PIPELINE_NAME_read_classification = "Read classification"
    PIPELINE_NAME_metagenomics_screening = "Screening"
    PIPELINE_NAME_request_mapping = "Request Mapping"
    PIPELINE_NAME_assembly = "Assembly"
    PIPELINE_NAME_remapping = "Remapping"
    PIPELINE_NAME_remap_filtering = "Remap filtering"
    PIPELINE_NAME_map_filtering = "Map filtering"
    PIPELINE_NAME_reporting = "Reporting"
    PIPELINE_NAME_metagenomics_settings = "Metagenomics settings"

    ## values to upload to database
    vect_pipeline_names = [
        PIPELINE_NAME_read_quality_analysis,
        PIPELINE_NAME_type_and_subtype_analysis,
        PIPELINE_NAME_variant_detection,
        PIPELINE_NAME_coverage_analysis,
        PIPELINE_NAME_alignment,
        PIPELINE_NAME_intra_host_minor_variant_detection,
        PIPELINE_NAME_extra_qc,
        PIPELINE_NAME_viral_enrichment,
        PIPELINE_NAME_host_depletion,
        PIPELINE_NAME_assembly,
        PIPELINE_NAME_contig_classification,
        PIPELINE_NAME_read_classification,
        PIPELINE_NAME_metagenomics_screening,
        PIPELINE_NAME_metagenomics_settings,
        PIPELINE_NAME_request_mapping,
        PIPELINE_NAME_remapping,
        PIPELINE_NAME_remap_filtering,
        PIPELINE_NAME_map_filtering,
        PIPELINE_NAME_reporting,
    ]

    vect_pipeline_televir_classic = [
        PIPELINE_NAME_extra_qc,
        PIPELINE_NAME_viral_enrichment,
        PIPELINE_NAME_host_depletion,
        PIPELINE_NAME_assembly,
        PIPELINE_NAME_contig_classification,
        PIPELINE_NAME_read_classification,
        PIPELINE_NAME_remapping,
        PIPELINE_NAME_remap_filtering,
    ]

    vect_pipeline_televir_metagenomics = [
        PIPELINE_NAME_extra_qc,
        PIPELINE_NAME_viral_enrichment,
        PIPELINE_NAME_host_depletion,
        PIPELINE_NAME_metagenomics_settings,
        PIPELINE_NAME_request_mapping,
        PIPELINE_NAME_remap_filtering,
        PIPELINE_NAME_reporting,
    ]

    vect_pipeline_televir_metagenomics_for_parameters = (
        vect_pipeline_televir_metagenomics + [PIPELINE_NAME_metagenomics_screening]
    )

    vect_pipeline_televir_screening = [
        PIPELINE_NAME_extra_qc,
        PIPELINE_NAME_viral_enrichment,
        PIPELINE_NAME_host_depletion,
        PIPELINE_NAME_metagenomics_screening,
        PIPELINE_NAME_metagenomics_settings,
    ]

    vect_pipeline_televir_mapping_only = [
        PIPELINE_NAME_extra_qc,
        PIPELINE_NAME_viral_enrichment,
        PIPELINE_NAME_host_depletion,
        PIPELINE_NAME_request_mapping,
        PIPELINE_NAME_map_filtering,
        PIPELINE_NAME_reporting,
    ]

    vect_pipeline_televir_request_mapping = [
        PIPELINE_NAME_extra_qc,
        PIPELINE_NAME_viral_enrichment,
        PIPELINE_NAME_host_depletion,
        PIPELINE_NAME_request_mapping,
        PIPELINE_NAME_map_filtering,
        PIPELINE_NAME_reporting,
    ]

    ###############################
    ### technology available
    TECHNOLOGY_illumina_old = "Illumina"
    TECHNOLOGY_illumina = "Illumina/IonTorrent"
    TECHNOLOGY_minion = "ONT"
    TECHNOLOGY_generic_old = "Generic"
    TECHNOLOGY_generic = "Extra"
    TECHNOLOGY_Undefined = "Undefined"

    ### technology names
    vect_technology = [
        TECHNOLOGY_generic,  ## must be first to appear on first in NavTab on Settings
        TECHNOLOGY_illumina,
        TECHNOLOGY_minion,
        TECHNOLOGY_Undefined,
    ]

    @staticmethod
    def reverse_set_dict(dict: Dict[str, str]):
        new_dict: Dict[str, list] = {}
        for key, value in dict.items():
            if value not in new_dict:
                new_dict[value] = []

            new_dict[value].append(key)
        return new_dict

    @property
    def vect_pipeline_names_condensed(self) -> Dict[str, List[str]]:
        pipeline_steps_dict = {
            pipeline_step: self.pipeline_step_to_pipeline_name(pipeline_step)
            for pipeline_step in self.vect_pipeline_names
        }

        pipeline_names_dict = self.reverse_set_dict(pipeline_steps_dict)

        return pipeline_names_dict

    @property
    def vect_pipeline_televir_metagenomics_condensed(self) -> Dict[str, List[str]]:
        pipeline_steps_dict = {
            pipeline_step: self.pipeline_step_to_pipeline_name(pipeline_step)
            for pipeline_step in self.vect_pipeline_televir_metagenomics_for_parameters
        }

        pipeline_names_dict = self.reverse_set_dict(pipeline_steps_dict)

        return pipeline_names_dict

    def pipeline_step_to_pipeline_name(self, pipeline_step: str) -> str:
        """
        Translate pipeline step names - use to combine steps."""
        if pipeline_step == self.PIPELINE_NAME_metagenomics_settings:
            return self.PIPELINE_NAME_metagenomics_screening
        return pipeline_step
