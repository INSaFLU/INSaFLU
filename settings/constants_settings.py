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
    PIPELINE_NAME_host_depletion = "Host depletion"
    PIPELINE_NAME_contig_classification = "Contig classification"
    PIPELINE_NAME_read_classification = "Read classification"
    PIPELINE_NAME_assembly = "Assembly"
    PIPELINE_NAME_remapping = "Remapping"
    PIPELINE_NAME_remap_filtering = "Remap filtering"
    PIPELINE_NAME_reporting = "Reporting"

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
        PIPELINE_NAME_remapping,
        PIPELINE_NAME_remap_filtering,
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
    
    @property
    def vect_pipeline_names_condensed(self) -> Dict[str, List[str]]:

        pipeline_steps_dict= {
            pipeline_step: self.pipeline_step_to_pipeline_name(pipeline_step)
            for pipeline_step in self.vect_pipeline_names
        }

        def reverse_set_dict(dict: Dict[str, str]):
            new_dict: Dict[str, list] = {}
            for key, value in dict.items():
                if value not in new_dict:
                    new_dict[value]= []

                new_dict[value].append(key)
            return new_dict
        
        pipeline_names_dict= reverse_set_dict(pipeline_steps_dict)

        return pipeline_names_dict

    def pipeline_step_to_pipeline_name(self, pipeline_step: str) -> str:
        """
        Translate pipeline step names - use to combine steps."""
        if pipeline_step == self.PIPELINE_NAME_remap_filtering:
            return self.PIPELINE_NAME_remapping

        return pipeline_step

    ###### Relation between software and technology
    def get_list_software_names_by_technology(self, technology_name):
        """
        :param technology_name TECHNOLOGY_illumina, TECHNOLOGY_minion
        """
        return self.DICT_SOFTWARE_RELATION.get(technology_name, [])

    def is_software_in_technology(self, technology_name, software_name):
        """
        :param technology_name TECHNOLOGY_illumina, TECHNOLOGY_minion
        """
        return software_name in self.get_list_software_names_by_technology(
            technology_name
        )
