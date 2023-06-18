"""
Ceated on 06/05/2022
@author: joao santos
"""

import os

import networkx as nx

from fluwebvirus.settings import MEDIA_ROOT, STATIC_ROOT, STATICFILES_DIRS
from pathogen_identification.models import Projects as PIprojects
from pathogen_identification.utilities.mapping_flags import MapFlagViruses
from settings.constants_settings import ConstantsSettings as CS
from settings.models import Software


class Pipeline_Makeup:
    """
    Pipeline steps
    """

    ROOT = "root"
    ASSEMBLY_SPECIAL_STEP = "ASSEMBLY_SPECIAL"
    VIRAL_ENRICHMENT_SPECIAL_STEP = "VIRAL_ENRICHMENT"

    dependencies_graph_edges = {
        CS.PIPELINE_NAME_viral_enrichment: [ROOT],
        VIRAL_ENRICHMENT_SPECIAL_STEP: [ROOT],
        CS.PIPELINE_NAME_host_depletion: [
            ROOT,
            CS.PIPELINE_NAME_viral_enrichment,
        ],
        CS.PIPELINE_NAME_read_classification: [
            ROOT,
            VIRAL_ENRICHMENT_SPECIAL_STEP,
            CS.PIPELINE_NAME_host_depletion,
        ],
        CS.PIPELINE_NAME_assembly: [
            ROOT,
            CS.PIPELINE_NAME_read_classification,
            CS.PIPELINE_NAME_host_depletion,
            VIRAL_ENRICHMENT_SPECIAL_STEP,
        ],
        ASSEMBLY_SPECIAL_STEP: [
            CS.PIPELINE_NAME_read_classification,
        ],
        CS.PIPELINE_NAME_contig_classification: [CS.PIPELINE_NAME_assembly],
        CS.PIPELINE_NAME_remapping: [
            CS.PIPELINE_NAME_contig_classification,
            CS.PIPELINE_NAME_read_classification,
            ASSEMBLY_SPECIAL_STEP,
        ],
    }

    dependencies_graph_root = CS.PIPELINE_NAME_remapping
    dependencies_graph_sink = ROOT

    def generate_dependencies_graph(self):
        """
        Generates a graph of dependencies between pipeline steps
        """
        G = nx.DiGraph()
        for (
            pipeline_step,
            dependencies,
        ) in self.dependencies_graph_edges.items():
            for dependency in dependencies:
                G.add_edge(pipeline_step, dependency)
        return G

    def process_path(self, dpath: list):
        """
        Processes the path to remove the root node
        """
        dpath = [
            x.replace(self.ASSEMBLY_SPECIAL_STEP, CS.PIPELINE_NAME_assembly).replace(
                self.VIRAL_ENRICHMENT_SPECIAL_STEP, CS.PIPELINE_NAME_viral_enrichment
            )
            for x in dpath
            if x != self.dependencies_graph_sink
        ]

        return dpath[::-1]

    def get_denpendencies_paths_dict(self):
        """
        Returns a dictionary with the dependencies between pipeline steps
        """
        G = self.generate_dependencies_graph()
        paths = nx.all_simple_paths(
            G,
            self.dependencies_graph_root,
            self.dependencies_graph_sink,
        )

        paths = {x: self.process_path(path) for x, path in enumerate(paths)}
        return paths

    def __init__(self):
        self.MAKEUP = self.get_denpendencies_paths_dict()

    def get_makeup(self, makeup: int):
        return self.MAKEUP.get(makeup, None)

    def get_makeup_name(self, makeup: int):
        return self.MAKEUP[makeup][0]

    def get_makeup_list(
        self,
    ):
        return list(self.MAKEUP.keys())

    def get_makeup_list_names(
        self,
    ):
        return list(self.MAKEUP.values())

    def match_makeup_name_from_list(self, makeup_list: list):
        for makeup, mlist in self.MAKEUP.items():
            if set(makeup_list) == set(mlist):
                return makeup
        return None

    def makeup_available(self, makeup: int):
        return makeup in self.MAKEUP

    def get_software_pipeline_list_including(
        self, software: Software, televir_project: PIprojects = None
    ):
        type_of_use = Software.TYPE_OF_USE_televir_global
        if televir_project:
            type_of_use = Software.TYPE_OF_USE_televir_project

        pipeline_steps_project = Software.objects.filter(
            type_of_use=type_of_use,
            technology=software.technology,
            parameter__televir_project=televir_project,
            is_to_run=True,
            owner=software.owner,
        ).values_list("pipeline_step__name", flat=True)
        return list(pipeline_steps_project)

    def get_software_pipeline_list_excluding(
        self, software: Software, televir_project: PIprojects = None
    ):
        type_of_use = Software.TYPE_OF_USE_televir_global
        if televir_project:
            type_of_use = Software.TYPE_OF_USE_televir_project

        pipeline_steps_project = (
            Software.objects.filter(
                type_of_use=type_of_use,
                technology=software.technology,
                parameter__televir_project=televir_project,
                is_to_run=True,
                owner=software.owner,
            )
            .exclude(pk=software.pk)
            .values_list("pipeline_step__name", flat=True)
        )

        return list(pipeline_steps_project)

    def get_pipeline_makeup_result_of_operation(
        self, software, turn_off=True, televir_project: PIprojects = None
    ):
        pipeline_steps_project = []

        if turn_off:
            pipeline_steps_project = self.get_software_pipeline_list_excluding(
                software, televir_project=televir_project
            )

        else:
            pipeline_steps_project = self.get_software_pipeline_list_including(
                software, televir_project=televir_project
            )

        return pipeline_steps_project


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
