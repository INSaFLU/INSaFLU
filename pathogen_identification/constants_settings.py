"""
Ceated on 06/05/2022
@author: joao santos
"""

from fluwebvirus.settings import MEDIA_ROOT, STATICFILES_DIRS
from settings.constants_settings import ConstantsSettings as CS


class Pipeline_Makeup:
    """
    Pipeline steps
    """

    MAKEUP_FULL = 0
    MAKEUP_NO_ASSEMBLY = 1
    MAKEUP_NO_ENRICHMENT = 2
    MAKEUP_NO_DEPLETION = 3
    MAKEUP_NO_READS = 4
    MAKEUP_READS_ENRICHMENT = 5
    MAKEUP_READS_DEPLETION = 6
    MAKEUP_READS_ONLY = 7
    MAKEUP_CONTIGS_ENRICHMENT = 8
    MAKEUP_CONTIGS_DEPLETION = 9
    MAKEUP_CONTIGS_ONLY = 10

    FULL = [
        CS.PIPELINE_NAME_viral_enrichment,
        CS.PIPELINE_NAME_host_depletion,
        CS.PIPELINE_NAME_assembly,
        CS.PIPELINE_NAME_contig_classification,
        CS.PIPELINE_NAME_read_classification,
        CS.PIPELINE_NAME_remapping,
    ]

    NO_ENRICHMENT = [
        CS.PIPELINE_NAME_host_depletion,
        CS.PIPELINE_NAME_assembly,
        CS.PIPELINE_NAME_contig_classification,
        CS.PIPELINE_NAME_read_classification,
        CS.PIPELINE_NAME_remapping,
    ]

    NO_DEPLETION = [
        CS.PIPELINE_NAME_viral_enrichment,
        CS.PIPELINE_NAME_assembly,
        CS.PIPELINE_NAME_contig_classification,
        CS.PIPELINE_NAME_read_classification,
        CS.PIPELINE_NAME_remapping,
    ]

    NO_ASSEMBLY = [
        CS.PIPELINE_NAME_viral_enrichment,
        CS.PIPELINE_NAME_host_depletion,
        CS.PIPELINE_NAME_read_classification,
        CS.PIPELINE_NAME_remapping,
    ]

    NO_READS = [
        CS.PIPELINE_NAME_viral_enrichment,
        CS.PIPELINE_NAME_host_depletion,
        CS.PIPELINE_NAME_assembly,
        CS.PIPELINE_NAME_contig_classification,
        CS.PIPELINE_NAME_remapping,
    ]

    READS_ONLY_ENRICHMENT = [
        CS.PIPELINE_NAME_viral_enrichment,
        CS.PIPELINE_NAME_read_classification,
        CS.PIPELINE_NAME_remapping,
    ]

    READS_ONLY_DEPLETION = [
        CS.PIPELINE_NAME_host_depletion,
        CS.PIPELINE_NAME_read_classification,
        CS.PIPELINE_NAME_remapping,
    ]

    READS_ONLY = [
        CS.PIPELINE_NAME_read_classification,
        CS.PIPELINE_NAME_remapping,
    ]

    CONTIGS_ONLY_ENRICHMENT = [
        CS.PIPELINE_NAME_viral_enrichment,
        CS.PIPELINE_NAME_assembly,
        CS.PIPELINE_NAME_contig_classification,
        CS.PIPELINE_NAME_remapping,
    ]

    CONTIGS_ONLY_DEPLETION = [
        CS.PIPELINE_NAME_viral_enrichment,
        CS.PIPELINE_NAME_assembly,
        CS.PIPELINE_NAME_contig_classification,
        CS.PIPELINE_NAME_remapping,
    ]

    CONTIGS_ONLY = [
        CS.PIPELINE_NAME_assembly,
        CS.PIPELINE_NAME_contig_classification,
        CS.PIPELINE_NAME_remapping,
    ]

    MAKEUP = {
        MAKEUP_FULL: FULL,
        MAKEUP_NO_ASSEMBLY: NO_ASSEMBLY,
        MAKEUP_NO_READS: NO_READS,
        MAKEUP_NO_ENRICHMENT: NO_ENRICHMENT,
        MAKEUP_NO_DEPLETION: NO_DEPLETION,
        MAKEUP_READS_ENRICHMENT: READS_ONLY_ENRICHMENT,
        MAKEUP_READS_DEPLETION: READS_ONLY_DEPLETION,
        MAKEUP_READS_ONLY: READS_ONLY,
        MAKEUP_CONTIGS_ENRICHMENT: CONTIGS_ONLY_ENRICHMENT,
        MAKEUP_CONTIGS_DEPLETION: CONTIGS_ONLY_DEPLETION,
        MAKEUP_CONTIGS_ONLY: CONTIGS_ONLY,
    }

    def __init__(self):
        pass

    @staticmethod
    def get_makeup(makeup: int):
        return Pipeline_Makeup.MAKEUP[makeup]

    @staticmethod
    def get_makeup_name(makeup: int):
        return Pipeline_Makeup.MAKEUP[makeup][0]

    @staticmethod
    def get_makeup_list():
        return list(Pipeline_Makeup.MAKEUP.keys())

    @staticmethod
    def get_makeup_list_names():
        return list(Pipeline_Makeup.MAKEUP.values())

    def match_makeup_name_from_list(self, makeup_list: list):

        for makeup, mlist in Pipeline_Makeup.MAKEUP.items():

            if set(makeup_list) == set(mlist):
                return makeup
        return None

    def makeup_available(self, makeup: int):
        return makeup in Pipeline_Makeup.MAKEUP


class ConstantsSettings:
    """
    Constants for pathogen identification
    """

    project_directory = "/tmp/televir/projects/"
    media_directory = MEDIA_ROOT
    static_directory = STATICFILES_DIRS[0]
    televir_subdirectory = "televir_projects"
    PAGINATE_NUMBER = 10
    docker_app_directory = "/televir/mngs_benchmark/"
    docker_install_directory = "/televir/mngs_benchmark/mngs_environments/"
    USER_TREE_INDEX = 0

    ###################################

    PIPELINE_STEPS_DB_DEPENDENT = [
        CS.PIPELINE_NAME_viral_enrichment,
        CS.PIPELINE_NAME_host_depletion,
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
