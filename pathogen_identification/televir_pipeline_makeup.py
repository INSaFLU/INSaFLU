import networkx as nx

from collections import defaultdict
from typing import List, Optional

from pathogen_identification.models import Projects

from settings.constants_settings import ConstantsSettings as CS
tree = lambda: defaultdict(tree)

from pathogen_identification.models import (
    PIProject_Sample,
)
from settings.models import Software

def excluded_steps_decorator(function):
    """
    create excluded steps given project"""

    def wrapped(
        self,
        software: Software,
        televir_project: Optional[Projects] = None,
        project_sample: Optional[PIProject_Sample] = None,
    ):
        exclude_steps = [CS.PIPELINE_NAME_reporting]

        if project_sample is None:
            exclude_steps.append(CS.PIPELINE_NAME_metagenomics_screening)

        return function(
            self,
            software,
            televir_project,
            project_sample,
            exclude_steps,
        )

    return wrapped


def make_tree(lst):
    d = tree()
    for x in lst:
        curr = d
        for item in x:
            curr = curr[item]
    return d


def differences_tuple_list(lista, listb):
    """
    Return the differences between two lists
    """
    list_a = [tuple([str(x) for x in y]) for y in lista]
    list_a = set(list_a)

    list_b = [tuple([str(x) for x in y]) for y in listb]
    list_b = set(list_b)
    return list(list_a.symmetric_difference(list_b))


#################
# TREE UTILITIES


class PipelineTreeBase:
    ROOT = "root"
    READ_CLASSIFICATION_SPECIAL_STEP = "ASSEMBLY_SPECIAL"
    VIRAL_ENRICHMENT_SPECIAL_STEP = "VIRAL_ENRICHMENT"
    MAP_FILTERING_SPECIAL_STEP = "MAP_FILTERING"
    SINK = "sink"
    dependencies_graph_root = SINK
    dependencies_graph_sink = ROOT


class Pipeline_Graph_Metagenomics(PipelineTreeBase):
    def __init__(self):
        self.dependencies_graph_edges_metagenomics = {
            CS.PIPELINE_NAME_extra_qc: [self.ROOT],
            CS.PIPELINE_NAME_viral_enrichment: [self.ROOT, CS.PIPELINE_NAME_extra_qc],
            self.VIRAL_ENRICHMENT_SPECIAL_STEP: [self.ROOT, CS.PIPELINE_NAME_extra_qc],
            CS.PIPELINE_NAME_host_depletion: [
                self.ROOT,
                CS.PIPELINE_NAME_extra_qc,
                CS.PIPELINE_NAME_viral_enrichment,
            ],
            CS.PIPELINE_NAME_assembly: [
                self.ROOT,
                CS.PIPELINE_NAME_extra_qc,
                CS.PIPELINE_NAME_host_depletion,
                self.VIRAL_ENRICHMENT_SPECIAL_STEP,
            ],

            CS.PIPELINE_NAME_contig_classification: [CS.PIPELINE_NAME_assembly],

            CS.PIPELINE_NAME_read_classification: [
                self.ROOT,
                CS.PIPELINE_NAME_contig_classification,
                CS.PIPELINE_NAME_extra_qc,
                self.VIRAL_ENRICHMENT_SPECIAL_STEP,
                CS.PIPELINE_NAME_host_depletion,
            ],
            #self.READ_CLASSIFICATION_SPECIAL_STEP: [
            #    CS.PIPELINE_NAME_read_classification,
            #],
            CS.PIPELINE_NAME_remap_filtering: [
                CS.PIPELINE_NAME_contig_classification,
                CS.PIPELINE_NAME_read_classification,
                #self.READ_CLASSIFICATION_SPECIAL_STEP,
                CS.PIPELINE_NAME_host_depletion,
            ],
            CS.PIPELINE_NAME_remapping: [
                CS.PIPELINE_NAME_remap_filtering,
                CS.PIPELINE_NAME_contig_classification,
                CS.PIPELINE_NAME_read_classification,
                #self.READ_CLASSIFICATION_SPECIAL_STEP,
                CS.PIPELINE_NAME_host_depletion,
                self.VIRAL_ENRICHMENT_SPECIAL_STEP,
            ],
            CS.PIPELINE_NAME_map_filtering: [
                self.ROOT,
                CS.PIPELINE_NAME_extra_qc,
                #self.READ_CLASSIFICATION_SPECIAL_STEP,
                CS.PIPELINE_NAME_host_depletion,
                self.VIRAL_ENRICHMENT_SPECIAL_STEP,
            ],
            self.MAP_FILTERING_SPECIAL_STEP: [
                self.ROOT,
                CS.PIPELINE_NAME_extra_qc,
                CS.PIPELINE_NAME_host_depletion,
                self.VIRAL_ENRICHMENT_SPECIAL_STEP,
            ],
            CS.PIPELINE_NAME_request_mapping: [
                self.ROOT,
                CS.PIPELINE_NAME_extra_qc,
                # self.ASSEMBLY_SPECIAL_STEP,
                self.MAP_FILTERING_SPECIAL_STEP,
                CS.PIPELINE_NAME_host_depletion,
                self.VIRAL_ENRICHMENT_SPECIAL_STEP,
            ],
            CS.PIPELINE_NAME_metagenomics_screening: [
                self.ROOT,
                CS.PIPELINE_NAME_extra_qc,
                self.MAP_FILTERING_SPECIAL_STEP,
                # self.ASSEMBLY_SPECIAL_STEP,
                CS.PIPELINE_NAME_host_depletion,
                self.VIRAL_ENRICHMENT_SPECIAL_STEP,
            ],
            self.SINK: [
                CS.PIPELINE_NAME_request_mapping,
                CS.PIPELINE_NAME_metagenomics_screening,
                CS.PIPELINE_NAME_remapping,
                CS.PIPELINE_NAME_contig_classification,
                #self.READ_CLASSIFICATION_SPECIAL_STEP,
                CS.PIPELINE_NAME_read_classification,
            ],
        }


class Pipeline_Makeup(PipelineTreeBase):
    def __init__(self):
        super().__init__()

        self.dependencies_graph_edges = (
            Pipeline_Graph_Metagenomics().dependencies_graph_edges_metagenomics
        )

        self.MAKEUP = self.get_dependencies_paths_dict()

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

    def process_path(self, dpath: List[str]):
        """
        Processes the path to remove the root node
        """
        dpath = [
            x
            .replace(
                self.VIRAL_ENRICHMENT_SPECIAL_STEP, CS.PIPELINE_NAME_viral_enrichment
            )
            .replace(self.MAP_FILTERING_SPECIAL_STEP, CS.PIPELINE_NAME_map_filtering)
            for x in dpath
            if x not in [self.ROOT, self.SINK]
        ]

        return dpath[::-1]

    def get_dependencies_paths_dict(self):
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

    def get_makeup(self, makeup: int) -> list:
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

    @property
    def get_pipeline_names(self):
        return list(self.dependencies_graph_edges.keys())

    def match_makeup_name_from_list(
        self, makeup_list: list, ignore: List[str] = []
    ) -> Optional[int]:
        makeup_safe = [x for x in makeup_list if x in self.get_pipeline_names]
        if ignore:
            makeup_safe = [x for x in makeup_safe if x not in ignore]

        for makeup, mlist in self.MAKEUP.items():
            if set(makeup_safe) == set(mlist):
                return makeup
        return None

    def check_makeuplist_has_classification(self, makeup_list: list) -> bool:
        classification_steps = [
            CS.PIPELINE_NAME_contig_classification,
            CS.PIPELINE_NAME_read_classification,
        ]
        return any([x in makeup_list for x in classification_steps])

    def match_makeup_name_from_list_classification(
        self, makeup_list: list
    ) -> Optional[int]:
        ignore = [
            CS.PIPELINE_NAME_metagenomics_screening,
            CS.PIPELINE_NAME_request_mapping,
            CS.PIPELINE_NAME_map_filtering,
            CS.PIPELINE_NAME_remap_filtering,
        ]

        makeup_return = self.match_makeup_name_from_list(makeup_list, ignore=ignore)

        if makeup_return is None:
            return None

        if not self.check_makeuplist_has_classification(makeup_list):
            return None

        return makeup_return

    def makeup_available(self, makeup: int) -> bool:
        return makeup in self.MAKEUP

    @excluded_steps_decorator
    def get_software_pipeline_list_including(
        self,
        software: Software,
        televir_project: Optional[Projects] = None,
        project_sample: Optional[PIProject_Sample] = None,
        exclude_steps: List[str] = [],
    ):
        use_types = Software.TELEVIR_GLOBAL_TYPES

        if televir_project:
            use_types = Software.TELEVIR_PROJECT_TYPES

        pipeline_steps_project = (
            Software.objects.filter(
                type_of_use__in=use_types,
                technology=software.technology,
                parameter__televir_project=televir_project,
                parameter__televir_project_sample=project_sample,
                is_to_run=True,
                owner=software.owner,
            )
            .exclude(pipeline_step__name__in=exclude_steps)
            .values_list("pipeline_step__name", flat=True)
        )

        pipeline_steps_project = list(pipeline_steps_project)

        if software.pipeline_step.name not in exclude_steps:
            pipeline_steps_project.append(software.pipeline_step.name)

        return pipeline_steps_project

    @excluded_steps_decorator
    def get_software_pipeline_list_excluding(
        self,
        software: Software,
        televir_project: Optional[Projects] = None,
        project_sample: Optional[PIProject_Sample] = None,
        exclude_steps: List[str] = [],
    ):
        use_types = Software.TELEVIR_GLOBAL_TYPES
        if televir_project:
            use_types = Software.TELEVIR_PROJECT_TYPES

        pipeline_steps_project = (
            Software.objects.filter(
                type_of_use__in=use_types,
                technology=software.technology,
                parameter__televir_project=televir_project,
                parameter__televir_project_sample=project_sample,
                is_to_run=True,
                owner=software.owner,
            )
            .exclude(pk=software.pk)
            .exclude(pipeline_step__name__in=exclude_steps)
            .values_list("pipeline_step__name", flat=True)
        )

        return list(pipeline_steps_project)

    def get_pipeline_makeup_result_of_operation(
        self,
        software,
        turn_off=True,
        televir_project: Optional[Projects] = None,
        project_sample: Optional[PIProject_Sample] = None,
    ):
        pipeline_steps_project = []

        if turn_off:
            pipeline_steps_project = self.get_software_pipeline_list_excluding(
                software, televir_project=televir_project, project_sample=project_sample
            )

        else:
            pipeline_steps_project = self.get_software_pipeline_list_including(
                software, televir_project=televir_project, project_sample=project_sample
            )

        return pipeline_steps_project

