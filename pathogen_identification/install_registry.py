from settings.constants_settings import ConstantsSettings as CS


class Params_Illumina:
    """
    Parameters for the Illumina pipeline
    """

    DATA_TYPE = "illumina"

    ################################## CONSTANT PARAMS

    CONSTANTS = {
        "minimum_coverage_threshold": 2,
        "max_output_number": 15,
        "taxid_limit": 12,
        "sift_query": "phage",
        "assembly_contig_min_length": 300,
    }

    ################################## SOFTWARE


class Params_Nanopore:
    """
    Class containing all the parameters for the nanopore pipeline
    """

    DATA_TYPE = "nanopore"

    CONSTANTS = {
        "minimum_coverage_threshold": 1,
        "max_output_number": 15,
        "taxid_limit": 12,
        "sift_query": "phage",
        "assembly_contig_min_length": 500,
    }
