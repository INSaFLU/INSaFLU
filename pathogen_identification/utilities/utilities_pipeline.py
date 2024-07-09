import itertools as it
import logging
import os
from collections import defaultdict
from typing import Dict, List, Optional, Tuple, Union

import networkx as nx
import numpy as np
import pandas as pd
from django.contrib.auth.models import User
from django.db.models import Q, QuerySet

from constants.constants import Televir_Directory_Constants as Televir_Directories
from constants.constants import Televir_Metadata_Constants as Televir_Metadata
from pathogen_identification.constants_settings import ConstantsSettings
from pathogen_identification.host_library import Host
from pathogen_identification.models import (
    ParameterSet,
    PIProject_Sample,
    Projects,
    RawReference,
    RawReferenceCompoundModel,
    RunMain,
    SoftwareTree,
    SoftwareTreeNode,
)
from pathogen_identification.utilities.utilities_general import merge_classes
from pathogen_identification.utilities.utilities_televir_dbs import Utility_Repository
from pathogen_identification.utilities.utilities_views import RawReferenceCompound
from settings.constants_settings import ConstantsSettings as CS
from settings.models import Parameter, PipelineStep, Software, Technology
from utils.lock_atomic_transaction import LockedAtomicTransaction

tree = lambda: defaultdict(tree)


def exclued_steps_decorator(function):
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
    ASSEMBLY_SPECIAL_STEP = "ASSEMBLY_SPECIAL"
    VIRAL_ENRICHMENT_SPECIAL_STEP = "VIRAL_ENRICHMENT"
    SINK = "sink"
    dependencies_graph_root = SINK
    dependencies_graph_sink = ROOT


class Pipeline_Graph(PipelineTreeBase):
    """
    Pipeline steps
    """

    def __init__(self):
        self.dependencies_graph_edges = {
            CS.PIPELINE_NAME_extra_qc: [self.ROOT],
            CS.PIPELINE_NAME_viral_enrichment: [self.ROOT, CS.PIPELINE_NAME_extra_qc],
            self.VIRAL_ENRICHMENT_SPECIAL_STEP: [self.ROOT, CS.PIPELINE_NAME_extra_qc],
            CS.PIPELINE_NAME_host_depletion: [
                self.ROOT,
                CS.PIPELINE_NAME_extra_qc,
                CS.PIPELINE_NAME_viral_enrichment,
            ],
            CS.PIPELINE_NAME_read_classification: [
                self.ROOT,
                CS.PIPELINE_NAME_extra_qc,
                self.VIRAL_ENRICHMENT_SPECIAL_STEP,
                CS.PIPELINE_NAME_host_depletion,
            ],
            CS.PIPELINE_NAME_assembly: [
                self.ROOT,
                CS.PIPELINE_NAME_extra_qc,
                CS.PIPELINE_NAME_read_classification,
                CS.PIPELINE_NAME_host_depletion,
                self.VIRAL_ENRICHMENT_SPECIAL_STEP,
            ],
            self.ASSEMBLY_SPECIAL_STEP: [
                CS.PIPELINE_NAME_read_classification,
            ],
            CS.PIPELINE_NAME_contig_classification: [CS.PIPELINE_NAME_assembly],
            CS.PIPELINE_NAME_remap_filtering: [
                CS.PIPELINE_NAME_contig_classification,
                CS.PIPELINE_NAME_read_classification,
                self.ASSEMBLY_SPECIAL_STEP,
            ],
            CS.PIPELINE_NAME_remapping: [
                CS.PIPELINE_NAME_remap_filtering,
                CS.PIPELINE_NAME_contig_classification,
                CS.PIPELINE_NAME_read_classification,
                self.ASSEMBLY_SPECIAL_STEP,
            ],
            self.SINK: [CS.PIPELINE_NAME_remapping],
        }


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
            CS.PIPELINE_NAME_read_classification: [
                self.ROOT,
                CS.PIPELINE_NAME_extra_qc,
                self.VIRAL_ENRICHMENT_SPECIAL_STEP,
                CS.PIPELINE_NAME_host_depletion,
            ],
            CS.PIPELINE_NAME_assembly: [
                self.ROOT,
                CS.PIPELINE_NAME_extra_qc,
                CS.PIPELINE_NAME_read_classification,
                CS.PIPELINE_NAME_host_depletion,
                self.VIRAL_ENRICHMENT_SPECIAL_STEP,
            ],
            self.ASSEMBLY_SPECIAL_STEP: [
                CS.PIPELINE_NAME_read_classification,
            ],
            CS.PIPELINE_NAME_contig_classification: [CS.PIPELINE_NAME_assembly],
            CS.PIPELINE_NAME_remap_filtering: [
                CS.PIPELINE_NAME_contig_classification,
                CS.PIPELINE_NAME_read_classification,
                self.ASSEMBLY_SPECIAL_STEP,
                CS.PIPELINE_NAME_host_depletion,
            ],
            CS.PIPELINE_NAME_remapping: [
                CS.PIPELINE_NAME_remap_filtering,
                CS.PIPELINE_NAME_contig_classification,
                CS.PIPELINE_NAME_read_classification,
                self.ASSEMBLY_SPECIAL_STEP,
                CS.PIPELINE_NAME_host_depletion,
                self.VIRAL_ENRICHMENT_SPECIAL_STEP,
            ],
            CS.PIPELINE_NAME_map_filtering: [
                self.ROOT,
                CS.PIPELINE_NAME_extra_qc,
                self.ASSEMBLY_SPECIAL_STEP,
                CS.PIPELINE_NAME_host_depletion,
                self.VIRAL_ENRICHMENT_SPECIAL_STEP,
            ],
            CS.PIPELINE_NAME_request_mapping: [
                self.ROOT,
                CS.PIPELINE_NAME_extra_qc,
                CS.PIPELINE_NAME_map_filtering,
                self.ASSEMBLY_SPECIAL_STEP,
                CS.PIPELINE_NAME_host_depletion,
                self.VIRAL_ENRICHMENT_SPECIAL_STEP,
            ],
            CS.PIPELINE_NAME_metagenomics_screening: [
                self.ROOT,
                CS.PIPELINE_NAME_extra_qc,
                CS.PIPELINE_NAME_map_filtering,
                self.ASSEMBLY_SPECIAL_STEP,
                CS.PIPELINE_NAME_host_depletion,
                self.VIRAL_ENRICHMENT_SPECIAL_STEP,
            ],
            self.SINK: [
                CS.PIPELINE_NAME_request_mapping,
                CS.PIPELINE_NAME_metagenomics_screening,
                CS.PIPELINE_NAME_remapping,
                CS.PIPELINE_NAME_contig_classification,
                self.ASSEMBLY_SPECIAL_STEP,
                CS.PIPELINE_NAME_read_classification,
            ],
        }


class Pipeline_Makeup(PipelineTreeBase):
    def __init__(self):
        super().__init__()

        if ConstantsSettings.METAGENOMICS:
            self.dependencies_graph_edges = (
                Pipeline_Graph_Metagenomics().dependencies_graph_edges_metagenomics
            )

        else:
            self.dependencies_graph_edges = Pipeline_Graph().dependencies_graph_edges

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
            x.replace(self.ASSEMBLY_SPECIAL_STEP, CS.PIPELINE_NAME_assembly).replace(
                self.VIRAL_ENRICHMENT_SPECIAL_STEP, CS.PIPELINE_NAME_viral_enrichment
            )
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

    @exclued_steps_decorator
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

    @exclued_steps_decorator
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


class PipelineTree:

    technology: str
    nodes: list
    edges: dict
    leaves: list
    makeup: int
    graph: nx.DiGraph
    node_index: pd.DataFrame

    def __init__(
        self,
        technology: str,
        nodes: list,
        edges: dict,
        leaves: list,
        makeup: int,
        sorted=True,
        software_tree_pk: int = 0,
    ):
        self.technology = technology

        if sorted:
            self.node_index = pd.DataFrame([[x] for x in nodes], columns=["node"])
            self.nodes = nodes
        else:
            self.node_index = pd.DataFrame([[x[1]] for x in nodes], columns=["node"])
            self.node_index.index = [x[0] for x in nodes]
            self.nodes = self.node_index.node.tolist()

        #
        self.edges = edges
        self.leaves = leaves
        self.edge_dict = [(x[0], x[1]) for x in self.edges]
        self.makeup = makeup
        self.software_tree_pk = software_tree_pk
        self.logger = logging.getLogger(f"{__name__}_{self.technology}_{self.makeup}")
        self.dag_dict = {
            z: [
                self.edge_dict[x][1]
                for x in range(len(self.edge_dict))
                if self.edge_dict[x][0] == z
            ]
            for z in self.node_index.index
        }

    def __eq__(self, other):
        diff_nodes = differences_tuple_list(self.nodes, other.nodes)
        diff_edges = differences_tuple_list(self.edges, other.edges)

        if len(diff_nodes) == 0 and len(diff_edges) == 0:
            return True
        else:
            return False

    def get_parents_dict(self):
        """
        Get a dictionary of parents for each node
        """
        self.generate_graph()
        parents_dict = {}

        for i, row in self.node_index.iterrows():
            parents = list(self.graph.predecessors(i))
            if len(parents) == 0:
                parents = None
            else:
                parents = parents[0]
            parents_dict[i] = parents
        return parents_dict

    def nested_data(self):
        nested_layers = [self.leaves]
        parents_dict = self.get_parents_dict()

        root = False
        while not root:
            new_layer = []
            for node in nested_layers[-1]:
                if node in parents_dict.keys():
                    if parents_dict[node] == None:
                        root = True
                    else:
                        new_layer.append(parents_dict[node])

                else:
                    root = True
            if len(new_layer) > 0:
                nested_layers.append(new_layer)

        nested_layers.reverse()
        return nested_layers

    def node_from_index(self, nix):
        return self.node_index.loc[nix].node

    def get_path_explicit(self, path: list) -> list:
        """return nodes names for nodes index list"""

        return [(x, self.node_from_index(x)) for x in path]

    def generate_graph(self):
        """
        Generate a graph of pipeline
        """

        self.graph = nx.DiGraph()

        self.graph.add_edges_from(self.edge_dict)
        self.graph.add_nodes_from(self.node_index.index.tolist())

    def get_all_graph_paths(self, sample: Optional[PIProject_Sample] = None) -> dict:
        """
        Get all possible paths in the pipeline
        """

        self.generate_graph()
        all_paths = list(nx.all_simple_paths(self.graph, 0, self.leaves))
        all_paths_explicit = [self.get_path_explicit(path) for path in all_paths]
        path_dict = {
            all_paths_explicit[x][-1][0]: self.df_from_path(path, sample=sample)
            for x, path in enumerate(all_paths)
        }

        return path_dict

    def get_all_graph_paths_explicit(self) -> dict:
        """
        Get all possible paths in the pipeline
        explicit -> return nodes names for nodes index list
        """

        self.generate_graph()
        all_paths = list(nx.all_simple_paths(self.graph, 0, self.leaves))
        all_paths = [self.get_path_explicit(path) for path in all_paths]
        path_dict = {path[-1][0]: path for x, path in enumerate(all_paths)}

        return path_dict

    def paths_to_node(self, node: int):
        """
        Get path to node
        """
        self.generate_graph()
        paths = nx.all_simple_paths(self.graph, 0, node)
        paths = [self.get_path_explicit(path) for path in paths]

        return paths

    def get_specific_leaf_paths_explicit(self, leaves: List[int]) -> dict:
        """
        Get all possible paths in the pipeline
        explicit -> return nodes names for nodes index list
        """

        all_paths = self.get_all_graph_paths_explicit()
        leaf_paths = {leaf: all_paths[leaf] for leaf in leaves}

        return leaf_paths

    def check_if_leaf_steps_exist_list(
        self, paths: List[list], sample: PIProject_Sample
    ) -> list:

        leaves = []
        for path in paths:
            leaves.extend(self.check_if_leaf_step_exists(path, sample))

        leaves = list(set(leaves))
        return leaves

    def check_if_leaf_step_exists(self, path: list, sample: PIProject_Sample) -> list:
        """
        get leafs for a given path, path does not need to be complete
        """
        if len(path) <= 1:
            return []
        parameter_set_utils = Parameter_DB_Utility()
        utils_pipeline_manager = Utility_Pipeline_Manager()

        ps = ParameterSet.objects.filter(sample=sample).exclude(leaf=None)

        software_tree_pks = set([x.leaf.software_tree.pk for x in ps])
        software_trees = SoftwareTree.objects.filter(
            pk__in=software_tree_pks
        ).distinct()

        leaves_collected = []

        for software_tree in software_trees:
            pipeline_tree = parameter_set_utils.convert_softwaretree_to_pipeline_tree(
                software_tree
            )
            try:
                node_match = utils_pipeline_manager.match_path_to_tree_find_cutoff(
                    path, pipeline_tree
                )

                if node_match[0] == 0:
                    continue

                leaves_for_matched_node = pipeline_tree.leaves_from_node_using_graph(
                    node_match[0]
                )

                index_nodes = []

                for x in leaves_for_matched_node:

                    try:
                        index_nodes.append(
                            SoftwareTreeNode.objects.get(
                                software_tree=software_tree, index=x
                            )
                        )
                    except SoftwareTreeNode.DoesNotExist:
                        pass
                    except SoftwareTreeNode.MultipleObjectsReturned:
                        pass

                parameter_sets = ParameterSet.objects.filter(
                    leaf__in=index_nodes,
                    status=ParameterSet.STATUS_FINISHED,
                    sample=sample,
                ).distinct()

                leaves_collected.extend([x.pk for x in parameter_sets])

            except Exception as e:
                print("Exception", e)
                raise e

        leaves_collected = list(set(leaves_collected))
        return leaves_collected

    def df_from_path(
        self, path_extensive: list, sample: Optional[PIProject_Sample] = None
    ) -> pd.DataFrame:
        """
        Generate a dataframe from a path
        """
        path_extensive = self.get_path_explicit(path_extensive)

        path = [x[1] for x in path_extensive]
        df = []
        path = [x for x in path if x[0] != "root"]
        current_module = None
        for ix, node in enumerate(path):
            node_type = node[2]
            node_leaves = []
            if sample is not None:

                node_leaves = self.check_if_leaf_step_exists(
                    path_extensive[: ix + 1], sample=sample
                )

            if ix == 0 and node_type != "module":
                self.logger.info("First node must be a module")
                return pd.DataFrame()

            if node_type == "module":
                if current_module:
                    for param, value in current_module["params"].items():
                        df.append(
                            [
                                current_module.get("module"),
                                current_module.get("software"),
                                param,
                                value,
                                list(set(current_module.get("leaves"))),
                            ]
                        )

                current_module = {
                    "module": node[0],
                    "software": node[1],
                    "params": {
                        f"{node[1].upper()}_ARGS": "",
                    },
                    "leaves": node_leaves,
                }
                continue

            current_module["params"][node[0]] = node[1]
            current_module["leaves"].extend(node_leaves)

        if current_module:
            for param, value in current_module["params"].items():
                df.append(
                    [
                        current_module.get("module"),
                        current_module.get("software"),
                        param,
                        value,
                        list(set(current_module.get("leaves"))),
                    ]
                )
        df = pd.DataFrame(
            df, columns=["module", "software", "parameter", "value", "leaves"]
        )
        return df

    def leaves_from_node(self, node, leaves=[]):
        """ """

        try:
            if len(self.dag_dict[node]) == 0:
                return [node]

            for n in self.dag_dict[node]:
                leaves.extend(self.leaves_from_node(n))
        except KeyError:
            self.logger.error(f"Node {node} is not in the DAG")

        return leaves

    def leaves_from_node_compress_recursive(self, node):
        """ """
        leaves = []
        if len(self.compress_dag_dict[node]) == 0:
            return [node]

        for n in self.compress_dag_dict[node]:
            leaves.extend(self.leaves_from_node_compress_recursive(n))

        return leaves

    def leaves_from_node_using_graph(self, node):
        """
        Get all leaves from a node"""
        self.generate_graph()

        leaves = nx.descendants(self.graph, node)

        leaves = [x for x in leaves if self.graph.out_degree(x) == 0]

        return leaves

    def reduced_tree(self, leaves_list: list) -> Tuple[dict, pd.DataFrame]:
        """trims paths not leading to provided leaves"""
        root_node = ("root", None, None)

        for i, n in enumerate(leaves_list):
            if n not in self.leaves:
                self.logger.info(f"Node {n} is not a leaf")

                leaves_list.remove(n)

        if len(leaves_list) == 0:
            return self.dag_dict, pd.DataFrame()

        paths = self.get_all_graph_paths_explicit()
        compressed_paths = {z: paths[z] for z in leaves_list}

        new_nodes = it.chain(*[x for x in compressed_paths.values()])

        new_nodes = sorted(new_nodes, key=lambda x: x[0])
        new_nodes_index = [x[0] for x in new_nodes]
        new_nodes = [x[1] for x in new_nodes if x[1] != root_node]
        new_nodes_no_duplicates_same_order = list(dict.fromkeys(new_nodes))
        new_nodes = new_nodes_no_duplicates_same_order

        new_node_index = self.node_index[self.node_index.index.isin(new_nodes_index)]
        new_dag_dict = {
            parent: [
                child
                for child in self.dag_dict[parent]
                if child in new_node_index.index
            ]
            for parent in new_node_index.index
        }

        return new_dag_dict, new_node_index

    def simplify_tree(self, links, root, party: list, nodes_compress=[], edge_keep=[]):
        """ """
        party.append(root)

        if root in self.leaves:
            nodes_compress.append([party[0], tuple(party)])
            return

        subix = {}

        if len(links) > 1:
            nnode = tuple(party)
            ori = nnode[0]
            nodes_compress.append([ori, nnode])

        for i, g in enumerate(links):
            if len(links) != 1:
                edge_keep.append([ori, g])
                party = []

            self.simplify_tree(
                self.dag_dict[g],
                g,
                party,
                nodes_compress=nodes_compress,
                edge_keep=edge_keep,
            )

        return nodes_compress, edge_keep

    def compress_tree(self):
        """ """
        nodes_compress, edges_compress = self.simplify_tree(
            self.dag_dict[0], 0, [], nodes_compress=[], edge_keep=[]
        )
        #
        self.nodes_compress = nodes_compress
        self.edge_compress = edges_compress

        self.compress_dag_dict = {
            x[0]: [y[1] for y in self.edge_compress if y[0] == x[0]]
            for x in self.nodes_compress
        }

    def split_modules(self):
        if self.nodes_compress is None:
            self.compress_tree()

        nodes_compress = self.nodes_compress
        edge_compress = self.edge_compress

        new_nodes = []
        for node in nodes_compress:
            internal_nodes = []
            internal_edges = []

            internal_splits = [0]
            module_name = self.node_index.loc[node[0]].node

            for ix, internal_node in enumerate(node[1]):
                internal_name = self.node_index.loc[internal_node].node
                is_module = internal_name[2] == "module"

                if is_module and ix != 0:
                    internal_splits.append(ix)

            if len(internal_splits) > 1:
                internal_splits.append(len(node[1]))
                node_connections = [x for x in edge_compress if x[0] == node[0]]
                node_children = [x[1] for x in node_connections]
                for ix, split in enumerate(internal_splits[:-1]):
                    internal_nodes.append(
                        [node[1][split], (node[1][split : internal_splits[ix + 1]])]
                    )

                for ix, split in enumerate(internal_nodes[:-1]):
                    internal_edges.append([split[0], internal_nodes[ix + 1][0]])

                edge_compress = [x for x in edge_compress if x not in node_connections]
                new_nodes.extend(internal_nodes)
                edge_compress.extend(internal_edges)
                edge_compress.extend(
                    [[new_nodes[-1][0], child] for child in node_children]
                )

            else:
                new_nodes.append(node)

        ### change node name for leaves if in branch
        def check_leaf_in_list(list):
            for leaf in self.leaves:
                if leaf in list:
                    return leaf

            return None

        nodes_with_leaves = {}
        for node in new_nodes:
            if node[0] in self.leaves:
                continue
            leaf_found = check_leaf_in_list(node[1])
            if leaf_found is not None:
                nodes_with_leaves[node[0]] = leaf_found

        for node, leaf in nodes_with_leaves.items():
            new_nodes.append([leaf, ()])
            edge_compress.append([node, leaf])

        self.nodes_compress = new_nodes
        self.edge_compress = edge_compress

        self.compress_dag_dict = {
            z: [
                self.edge_compress[x][1]
                for x in range(len(self.edge_compress))
                if self.edge_compress[x][0] == z
            ]
            for z in list(set([x[0] for x in self.nodes_compress]))
        }

    def same_module_children(self, node, party, branches=[]):
        """ """
        children = self.compress_dag_dict[node]

        if len(children) == 0:
            return branches

        for child in children:
            child_name = self.node_index.loc[child].node
            new_party = party.copy()

            if child_name[2] != "module":
                new_party.append(child)
                self.same_module_children(child, new_party, branches=branches)

            else:
                if len(new_party) > 0:
                    branches.append({"branch": tuple(new_party), "exit": (node, child)})
                new_party = []

        return branches

    def get_module_tree(self):
        if self.nodes_compress is None:
            self.compress_tree()
            self.split_modules()

        original_nodes = pd.DataFrame(self.nodes_compress, columns=["node", "branch"])
        nodes_df = original_nodes.copy()
        original_edge_df = pd.DataFrame(self.edge_compress, columns=["parent", "child"])
        edge_df = original_edge_df.copy()

        def edit_branches(
            node,
            branches,
            nodes_df_small,
            edge_df_small,
            parent_node,
            original_nodes_df,
        ):
            """ """
            new_edges = []
            new_nodes = []
            for branch_meta in branches:
                branch = branch_meta["branch"]
                exit_edge = branch_meta["exit"]
                if len(branch) == 1:
                    continue

                recovered_branch = []
                for node in branch:
                    internal_nodes = original_nodes_df[
                        original_nodes_df.node == node
                    ].branch.values[0]
                    recovered_branch.extend(internal_nodes)

                recovered_branch = tuple(set(recovered_branch))
                branch = sorted(recovered_branch)
                nodes_df_small = nodes_df_small[~nodes_df_small.node.isin(branch)]
                edge_df_small = edge_df_small[
                    ~edge_df_small.parent.isin([x for x in branch if x != exit_edge[0]])
                ]
                edge_df_small = edge_df_small[~edge_df_small.child.isin(branch)]

                new_edges.append([parent_node, exit_edge[0]])
                new_nodes.append([exit_edge[0], tuple(branch)])

            return new_nodes, new_edges, nodes_df_small, edge_df_small

        def merge_new_branches(new_nodes, new_edges, nodes_df, edge_df):
            new_nodes = pd.DataFrame(new_nodes, columns=["node", "branch"])
            new_edges = pd.DataFrame(new_edges, columns=["parent", "child"])

            nodes_df = pd.concat([nodes_df, new_nodes], ignore_index=True)
            edge_df = pd.concat([edge_df, new_edges], ignore_index=True)

            return nodes_df, edge_df

        for node in self.nodes_compress:
            node_name = self.node_index.loc[node[0]].node

            if not node_name[1] is None:
                if not node_name[2] == "module":
                    continue

            parent_node = original_edge_df[
                original_edge_df.child == node[0]
            ].parent.values
            child_node = original_edge_df[
                original_edge_df.parent == node[0]
            ].child.values

            if len(parent_node) == 0:
                parent_node = [0]

            parent_node = parent_node[0]

            same_module_branches = self.same_module_children(
                node[0], [node[0]], branches=[]
            )

            if node_name == ("root", None, None):
                same_module_branches = [
                    {
                        "branch": [[node[0], self.compress_dag_dict[node[0]][0]]],
                        "exit": (0, self.compress_dag_dict[node[0]][0]),
                    }
                ]

            new_nodes, new_edges, nodes_df, edge_df = edit_branches(
                node[0],
                same_module_branches,
                nodes_df,
                edge_df,
                parent_node,
                original_nodes,
            )

            nodes_df, edge_df = merge_new_branches(
                new_nodes, new_edges, nodes_df, edge_df
            )

        self.nodes_compress = (
            nodes_df.drop_duplicates(subset=["node"]).to_numpy().tolist()
        )

        self.edge_compress = edge_df.drop_duplicates().to_numpy().tolist()

        self.compress_dag_dict = {
            z: [
                self.edge_compress[x][1]
                for x in range(len(self.edge_compress))
                if self.edge_compress[x][0] == z
            ]
            for z in list(set([x[0] for x in self.nodes_compress]))
        }


class Utility_Pipeline_Manager:
    """
    Takes a combined table and generates a pipeline tree.
    Combined table is a table with the following columns:
    - sample_name
    - sample_id
    - pipeline_step
    - software_name
    - parameter
    - value

    Uitility_Pipeline_Manager Comunicates with utility repository to complement the information
    with software specific installed databases. Creates a pipeline tree from the combined information.
    """

    software_name_list: list
    existing_pipeline_order: list
    combined_table: pd.DataFrame
    software_dbs_dict: dict
    technology: str
    new_variable: int
    pipeline_order: list
    pipeline_makeup: int

    def __init__(self):
        self.utility_repository = Utility_Repository(
            db_path=Televir_Directories.docker_app_directory,
            install_type="docker",
        )

        self.steps_db_dependant = ConstantsSettings.PIPELINE_STEPS_DB_DEPENDENT
        self.binaries = Televir_Metadata.BINARIES

        self.logger = logging.getLogger(__name__)
        if self.logger.hasHandlers():
            self.logger.handlers.clear()
        self.logger.setLevel(logging.ERROR)
        self.logger.addHandler(logging.StreamHandler())
        self.host_dbs = {}

    def input(self, combined_table: pd.DataFrame, technology="ONT"):
        """
        Input a combined table

        Args:

            combined_table (pd.DataFrame): Combined table with the following columns:
            - sample_name
            - sample_id
            - pipeline_step
            - software_name
            - parameter
            - value

            technology (str): Technology of the pipeline. Default is "ONT"
        """
        self.technology = technology

        combined_table = self.process_combined_table(combined_table)
        pipe_makeup_manager = Pipeline_Makeup()

        pipelines_available = combined_table.pipeline_step.unique().tolist()

        self.pipeline_makeup = pipe_makeup_manager.match_makeup_name_from_list(
            pipelines_available
        )

        if self.pipeline_makeup is None:
            self.logger.info("No pipeline makeup found")
            return False

        self.pipeline_order = pipe_makeup_manager.get_makeup(self.pipeline_makeup)

        self.existing_pipeline_order = [
            x for x in self.pipeline_order if x in pipelines_available
        ]

        self.software_name_list = combined_table.software_name.unique().tolist()

        self.combined_table = combined_table

        return True

    def process_combined_table(self, combined_table):
        """
        Process the combined table to make it compatible with the pipeline tree
        """
        combined_table = combined_table[
            combined_table.technology.str.contains(self.technology)
        ]

        def round_parameter_float(row):
            """
            Round float parameters to 2 decimals
            """
            new_param = row.parameter

            if row.type_data == Parameter.PARAMETER_float:
                if row.parameter and row.can_change:
                    new_param = round(float(row.parameter), 3)
                    new_param = str(new_param)

            return new_param

        combined_table["parameter"] = combined_table.apply(
            round_parameter_float, axis=1
        )
        return combined_table

    def generate_default_software_tree(self) -> PipelineTree:
        """
        Generate a default software tree
        """

        self.get_software_db_dict()

        self.generate_software_parameter_dict()

        return self.create_pipe_tree()

    def compare_software_trees(self, new_tree: PipelineTree):
        """
        Compare two software trees and return the differences
        """
        old_tree = self.generate_default_software_tree()

        return new_tree == old_tree

    def check_software_is_installed(self, software_name: str) -> bool:
        """
        Check if a software is installed
        """
        software_lower = software_name.lower()
        if software_lower in self.binaries["software"].keys():
            bin_path = os.path.join(
                Televir_Directories.docker_install_directory,
                self.binaries["software"][software_lower],
                "bin",
                software_lower,
            )
            return os.path.isfile(bin_path)
        else:
            for pipeline in [
                CS.PIPELINE_NAME_remapping,
                CS.PIPELINE_NAME_read_quality_analysis,
                CS.PIPELINE_NAME_extra_qc,
                CS.PIPELINE_NAME_assembly,
            ]:
                if os.path.exists(
                    os.path.join(
                        Televir_Directories.docker_install_directory,
                        self.binaries[pipeline]["default"],
                        "bin",
                        software_lower,
                    )
                ):
                    return True

        return False

    def check_software_db_available(self, software_name: str) -> bool:
        """
        Check if a software is installed
        """

        return self.utility_repository.check_exists(
            "software", "name", software_name.lower()
        )

    def set_software_list(self, software_list):
        self.software_name_list = software_list

    def get_software_list(self):
        self.software_name_list = Software.objects.filter(
            type_of_use__in=Software.TELEVIR_GLOBAL_TYPES,
        ).values_list("name", flat=True)

    def get_software_db_dict(self):
        software_list = self.utility_repository.get_list_unique_field(
            "software", "name"
        )

        self.software_dbs_dict = {
            software.lower(): self.get_software_dbs_if_exist(software)
            .path.unique()
            .tolist()
            for software in software_list
        }

    def get_host_dbs(self):
        software_list = self.utility_repository.get_list_unique_field(
            "software", "name"
        )
        hosts_dbs_dict = {
            software.lower(): self.get_software_dbs_if_exist(
                software,  # filters=[("tag", "host")]
            )
            for software in software_list
        }

        hosts_dbs_dict = {k: v for k, v in hosts_dbs_dict.items() if len(v) > 0}

        def recover_host(database: str) -> Host:
            for host in Host.__subclasses__():
                if database.startswith(host().host_name):
                    return host()

            return None

        def get_name_filename(row: pd.Series) -> pd.Series:
            host = recover_host(row.database)
            if host is None:
                row["host_name"] = np.nan
                row["host_filename"] = row.database
                row["file_str"] = f"{row.database}"

            else:
                row["host_name"] = host.host_name
                filename_simple = host.remote_filename
                row["host_filename"] = filename_simple
                row["file_str"] = f"{host.host_name} - {filename_simple}"

            return row

        for software in hosts_dbs_dict.keys():
            hosts_dbs_dict[software] = hosts_dbs_dict[software].apply(
                get_name_filename, axis=1
            )
            hosts_dbs_dict[software] = hosts_dbs_dict[software].dropna(
                subset=["host_name"], axis=0
            )

        self.host_dbs = hosts_dbs_dict

    def get_from_software_db_dict(self, software_name: str, empty=[]):
        possibilities = [software_name, software_name.lower()]
        if "_" in software_name:
            element = software_name.split("_")[0]

            possibilities.append(element)
            possibilities.append(element.lower())

        for possibility in possibilities:
            if possibility in self.software_dbs_dict.keys():
                return self.software_dbs_dict[possibility]

        return empty

    def get_from_host_db(self, software_name: str, empty=[]):
        possibilities = [software_name, software_name.lower()]
        if "_" in software_name:
            element = software_name.split("_")[0]

            possibilities.append(element)
            possibilities.append(element.lower())

        for possibility in possibilities:
            if possibility in self.host_dbs.keys():
                host_df = self.host_dbs[possibility]
                return list(
                    host_df[["path", "file_str"]].itertuples(index=False, name=None)
                )

        return empty

    def get_software_dbs_if_exist(
        self, software_name: str, filters: List[tuple] = []
    ) -> pd.DataFrame:
        fields = self.utility_repository.select_explicit_statement(
            "software", "name", software_name.lower(), filters=filters
        )

        try:
            fields = pd.read_sql(fields, self.utility_repository.engine)
            fields = fields.drop_duplicates(subset=["database"])

            return fields
        except Exception as e:
            self.logger.error(
                f"failed to fail to pandas read_sql {self.utility_repository.engine} software table for {software_name}. Error: {e}"
            )
            return pd.DataFrame(
                columns=["name", "path", "database", "installed", "env_path"]
            )

    def generate_argument_combinations(
        self, pipeline_software_dt: pd.DataFrame
    ) -> list:
        """
        Generate a list of argument combinations
        """
        pipeline_software_dict = {
            parameter_name: g.parameter.to_list()
            for parameter_name, g in pipeline_software_dt.groupby("parameter_name")
        }

        argument_combinations = []
        for flag, arg in pipeline_software_dict.items():
            flag_list = [" ".join([flag, x]) for x in arg]
            argument_combinations.append(flag_list)

        argument_combinations = list(it.product(*argument_combinations))
        argument_combinations = [" ".join(x) for x in argument_combinations]
        argument_combinations = list(set(argument_combinations))
        return argument_combinations

    def generate_software_parameter_dict(self) -> dict:
        """
        Generate a dictionary of software and parameters
        """
        step_dict = {c: g for c, g in self.combined_table.groupby("pipeline_step")}

        step_software_dict = {
            step: g.software_name.unique().tolist() for step, g in step_dict.items()
        }
        step_software_parameter_dict = {
            step: {
                software.lower(): {
                    f"{software.upper()}_ARGS": self.generate_argument_combinations(g),
                }
                for software, g in g.groupby("software_name")
            }
            for step, g in step_dict.items()
        }
        self.params_lookup = step_software_parameter_dict
        self.pipeline_software = {
            x: [f.lower() for f in y] for x, y in step_software_dict.items()
        }
        #

    def fill_dict(self, ix, branch_left) -> dict:
        """
        Generate a tree of pipeline
        """
        if ix == (len(self.existing_pipeline_order)):
            return branch_left

        if branch_left:
            return {
                node: self.fill_dict(ix, branch) for node, branch in branch_left.items()
            }

        else:
            current = self.existing_pipeline_order[ix]

        def fill_dict_chain_software_module(current) -> dict:
            soft_dict = {}
            suffix = ""
            soft = None

            param_combs = []
            if self.pipeline_software[current] == []:
                return {}

            for soft in self.pipeline_software[current]:
                if soft in self.params_lookup[current].keys():
                    # param_names = []

                    params_dict = self.params_lookup[current][soft]

                    for i, g in params_dict.items():
                        if not g:
                            continue

                        # param_names.append(i + suffix)
                        param_combs.append([(i + suffix, x, "param") for x in g])

                else:
                    param_combs.append([(f"{soft.upper()}_ARGS", "None", "param")])
                    # soft_dict[soft] = {(f"{soft.upper()}_ARGS", "None", "param"): {}}

            param_combs = list(it.product(*param_combs))
            param_tree = make_tree(param_combs)

            if soft is not None:
                soft_dict[soft] = param_tree

            return soft_dict

        def fill_dict_single_software_module(current) -> dict:
            soft_dict = {}
            suffix = ""

            for soft in self.pipeline_software[current]:
                if soft in self.params_lookup[current].keys():
                    param_names = []
                    param_combs = []
                    params_dict = self.params_lookup[current][soft]

                    for i, g in params_dict.items():
                        if not g:
                            continue

                        param_names.append(i + suffix)
                        param_combs.append([(i + suffix, x, "param") for x in g])

                    param_combs = list(it.product(*param_combs))

                    param_tree = make_tree(param_combs)

                    soft_dict[soft] = param_tree
                else:
                    soft_dict[soft] = {(f"{soft.upper()}_ARGS", "None", "param"): {}}

            return soft_dict

        if current in ConstantsSettings.PIPELINE_STEPS_AGGREGATE:
            soft_dict = fill_dict_chain_software_module(current)
        else:
            soft_dict = fill_dict_single_software_module(current)

        return {
            (current, soft, "module"): self.fill_dict(ix + 1, g)
            for soft, g in soft_dict.items()
        }

    def tree_index(self, tree, root, node_index=[], edge_dict=[], leaves=[]):
        """ """
        if len(tree) == 0:
            leaves.append(root)
            return

        subix = {}
        for i, g in tree.items():
            ix = len(node_index)
            node_index.append([ix, i])
            edge_dict.append([root, ix])
            subix[i] = ix

        #
        td = [
            self.tree_index(
                g, subix[i], node_index=node_index, edge_dict=edge_dict, leaves=leaves
            )
            for i, g in tree.items()
        ]

        return node_index, edge_dict, leaves

    def create_pipe_tree(self) -> PipelineTree:
        """ """
        self.pipeline_tree = self.fill_dict(0, {})

        node_index, edge_dict, leaves = self.tree_index(
            self.pipeline_tree,
            0,
            node_index=[(0, ("root", None, None))],
            edge_dict=[],
            leaves=[],
        )

        self.node_index = node_index

        return PipelineTree(
            technology=self.technology,
            nodes=[x[1] for x in node_index],
            edges=edge_dict,
            leaves=leaves,
            makeup=self.pipeline_makeup,
        )

    def generate_explicit_edge_dict(self, pipeline_tree: PipelineTree) -> dict:
        """ """
        nodes_dict = {(i, x): [] for i, x in enumerate(pipeline_tree.nodes)}

        for edge in pipeline_tree.edges:
            parent = (edge[0], pipeline_tree.nodes[edge[0]])
            child = (edge[1], pipeline_tree.nodes[edge[1]])
            nodes_dict[parent].append(child)

        nodes_dict = {
            x: pd.DataFrame([[z] for z in y], columns=["child"]).set_index("child")
            for x, y in nodes_dict.items()
        }

        return nodes_dict

    def node_index_dict(self, pipe_tree: PipelineTree) -> dict:
        """ """
        return {(i, x): i for i, x in enumerate(pipe_tree.nodes)}

    def match_path_to_tree_safe(self, explicit_path: list, pipe_tree: PipelineTree):
        """"""
        try:
            matched_path = self.match_path_to_tree(explicit_path, pipe_tree)
        except Exception as e:
            print(f"Path {explicit_path} not found in pipeline tree.")
            print("Exception:")
            print(e)
            return None

        return matched_path

    def match_path_to_tree(self, explicit_path: list, pipe_tree: PipelineTree):
        """"""

        self.logger.info("Matching path to tree")

        self.logger.info("Generating node index dict")
        nodes_index_dict = self.node_index_dict(pipe_tree)
        self.logger.info("Generating explicit edge dict")
        explicit_edge_dict = self.generate_explicit_edge_dict(pipe_tree)

        parent = explicit_path[0]
        parent_main = (0, ("root", None, None))
        child_main = None

        def match_nodes(node, node_list):
            for nd in node_list:
                if node[1] == nd[1]:
                    return nd

            return node

        self.logger.info("Initialize matching nodes")
        self.logger.info(f"Parent: {parent}")
        self.logger.info(f"Parent main: {parent_main}")
        self.logger.info(f"Child main: {child_main}")
        self.logger.info("Matching nodes iterating through explicit path")
        self.logger.info(f"leaves {pipe_tree.leaves}")

        for child in explicit_path[1:]:
            self.logger.info("--------------------")
            self.logger.info(f"Parent: {parent}")
            self.logger.info(f"Parent main: {parent_main}")
            self.logger.info(f"Child: {child}")

            try:
                child_main = match_nodes(
                    child, explicit_edge_dict[parent_main].index.tolist()
                )

            except KeyError:
                self.logger.info(f"{parent_main} not in parent tree edge dictionary.")
                return None

            self.logger.info(f"Child main: {child_main}")

            try:
                nodes_index_dict[child_main]
            except KeyError:
                self.logger.info(f"{child_main} node not in tree nodes")
                return None

            if nodes_index_dict[child_main] in pipe_tree.leaves:
                return nodes_index_dict[child_main]

            if child_main not in explicit_edge_dict[parent_main].index:
                self.logger.info(f"Child {child} not in parent {parent}")
                return None

            parent = child
            parent_main = child_main

    def match_path_to_tree_find_cutoff(
        self, explicit_path: list, pipe_tree: PipelineTree
    ) -> Tuple[int, tuple]:
        """"""

        self.logger.info("Matching path to tree")

        self.logger.info("Generating node index dict")
        nodes_index_dict = self.node_index_dict(pipe_tree)
        self.logger.info("Generating explicit edge dict")
        explicit_edge_dict = self.generate_explicit_edge_dict(pipe_tree)
        parent = explicit_path[0]
        parent_main = (0, ("root", None, None))
        child_main = None

        def match_nodes(node, node_list):
            for nd in node_list:
                if node[1] == nd[1]:
                    return nd

            return node

        self.logger.info("Initialize matching nodes")
        self.logger.info(f"Parent: {parent}")
        self.logger.info(f"Parent main: {parent_main}")
        self.logger.info(f"Child main: {child_main}")
        self.logger.info("Matching nodes iterating through explicit path")
        self.logger.info(f"leaves {pipe_tree.leaves}")

        for child in explicit_path[1:]:

            self.logger.info("--------------------")
            self.logger.info(f"Parent: {parent}")
            self.logger.info(f"Parent main: {parent_main}")
            self.logger.info(f"Child: {child}")

            try:
                child_main = match_nodes(
                    child, explicit_edge_dict[parent_main].index.tolist()
                )

            except KeyError:
                self.logger.info(f"{parent_main} not in parent tree edge dictionary.")
                return parent_main

            self.logger.info(f"Child main: {child_main}")

            try:
                nodes_index_dict[child_main]
            except KeyError:
                self.logger.info(f"{child_main} node not in tree nodes")
                return parent_main

            if nodes_index_dict[child_main] in pipe_tree.leaves:
                return child_main

            if child_main not in explicit_edge_dict[parent_main].index:
                self.logger.info(f"Child {child} not in parent {parent}")
                return parent_main

            parent = child
            parent_main = child_main

        return parent_main

    @staticmethod
    def pipe_tree_reconstruct(
        nodes_index_dict: dict, explicit_edge_dict: dict, technology: str, makeup: int
    ):
        """
        Reconstruct pipeline tree from nodes index dict and explicit edge dict
        """
        nodes = sorted(nodes_index_dict.keys(), key=lambda x: x[0])
        nodes = [x[1] for x in nodes]
        edges = []

        explicit_edge_dict = {
            x: y for x, y in explicit_edge_dict.items() if x in nodes_index_dict.keys()
        }

        for node, child_index in explicit_edge_dict.items():
            for child in child_index.index:
                edges.append([nodes_index_dict[node], nodes_index_dict[child]])

        leaves = []
        for node in nodes_index_dict.keys():
            if node in explicit_edge_dict.keys():
                if len(explicit_edge_dict[node]) == 0:
                    leaves.append(nodes_index_dict[node])

        return PipelineTree(
            technology=technology,
            nodes=nodes,
            edges=edges,
            leaves=leaves,
            makeup=makeup,
            sorted=True,
        )

    def match_path_to_tree_extend(
        self, explicit_path: list, pipe_tree: PipelineTree
    ) -> PipelineTree:
        """
        Match explicit path to pipeline tree and extend the tree if necessary
        """

        self.logger.info("Matching path to tree")
        tree_nodes = pipe_tree.nodes.copy()
        self.logger.info("Generating node index dict")
        nodes_index_dict = self.node_index_dict(pipe_tree)
        self.logger.info("Generating explicit edge dict")
        explicit_edge_dict = self.generate_explicit_edge_dict(pipe_tree)

        (
            nodes_index_dict_ext,
            explicit_edge_dict_ext,
            tree_nodes_ext,
        ) = self.extend_tree_dicts(
            explicit_path, nodes_index_dict, explicit_edge_dict, tree_nodes
        )

        extended_tree = self.pipe_tree_reconstruct(
            nodes_index_dict_ext,
            explicit_edge_dict_ext,
            pipe_tree.technology,
            pipe_tree.makeup,
        )

        return extended_tree

    def extend_tree_dicts(
        self,
        explicit_path: list,
        nodes_index_dict: dict,
        explicit_edge_dict: dict,
        tree_nodes: list,
    ):
        parent = explicit_path[0]
        parent_main = (0, ("root", None, None))
        child_main = None

        def match_nodes(node, node_list):
            for nd in node_list:
                if node[1] == nd[1]:
                    return nd

            raise KeyError

        def add_node(node, nodes: list):
            """Add node to nodes list and return index"""
            # if node[1] in nodes:
            #    return nodes.index(node[1])
            # else:
            nodes.append(node[1])
            return len(nodes) - 1

        def update_nodes_index(new_node, df: pd.DataFrame):
            """Add new node to nodes index dict"""
            df = df.append(
                pd.DataFrame([[new_node]], columns=["child"]).set_index("child")
            )
            return df

        def update_explicit_edge_dict(
            explicit_edge_dict: Dict[tuple, pd.DataFrame],
            parent_main: tuple,
            child_main: tuple,
        ):
            """Add new node to nodes index dict"""

            if parent_main in explicit_edge_dict.keys():
                explicit_edge_dict[parent_main] = update_nodes_index(
                    child_main, explicit_edge_dict[parent_main]
                )

            else:
                explicit_edge_dict[parent_main] = pd.DataFrame(
                    [[child_main]], columns=["child"]
                ).set_index("child")

            return explicit_edge_dict

        self.logger.info("Initialize matching nodes")
        self.logger.info(f"Parent: {parent}")
        self.logger.info(f"Parent main: {parent_main}")
        self.logger.info(f"Child main: {child_main}")
        self.logger.info("Matching nodes iterating through explicit path")

        if parent_main not in explicit_edge_dict.keys():
            explicit_edge_dict[parent_main] = pd.DataFrame(columns=["child"]).set_index(
                "child"
            )
            nodes_index_dict[parent_main] = 0

        for child in explicit_path[1:]:
            self.logger.info("--------------------")
            self.logger.info(f"Parent: {parent}")
            self.logger.info(f"Parent main: {parent_main}")
            self.logger.info(f"Child: {child}")

            try:
                child_main = match_nodes(
                    child, explicit_edge_dict[parent_main].index.tolist()
                )

            except KeyError:
                child_main = (add_node(child, tree_nodes), child[1])
                nodes_index_dict[child_main] = child_main[0]

                explicit_edge_dict = update_explicit_edge_dict(
                    explicit_edge_dict, parent_main, child_main
                )

            self.logger.info(f"Child main: {child_main}")

            try:
                nodes_index_dict[child_main]
            except KeyError:
                self.logger.info(f"{child_main} node not in tree nodes")
                # return None

            if child_main not in explicit_edge_dict[parent_main].index:
                self.logger.info(f"Child {child} not in parent {parent}")
                # return None

            parent = child
            parent_main = child_main

        if child_main not in explicit_edge_dict.keys():
            explicit_edge_dict[child_main] = pd.DataFrame(columns=["child"]).set_index(
                "child"
            )

        return nodes_index_dict, explicit_edge_dict, tree_nodes

    def update_tree_paths(self, pipe_tree: PipelineTree, explicit_paths: list):
        for path in explicit_paths:
            pipe_tree = self.match_path_to_tree_extend(path, pipe_tree)

    def compress_software_tree(self, software_tree: PipelineTree):
        software_tree.compress_tree()

        software_tree.split_modules()

        software_tree.get_module_tree()

        return software_tree


class Parameter_DB_Utility:
    """Comunicates with dango database.
    Functions:
        - create combined table of software and parameters.
        - get and store pipeline trees in tables SoftwareTree and SoftwareTree_node
    """

    def __init__(self):
        self.logger = logging.getLogger(__name__)
        if self.logger.hasHandlers():
            self.logger.handlers.clear()
        self.logger.setLevel(logging.ERROR)
        self.logger.addHandler(logging.StreamHandler())
        self.televir_constants = ConstantsSettings()

    def get_technologies_available(self):
        """
        Get technologies available for a user
        """

        technologies_available = (
            Software.objects.filter(type_of_use=Software.TYPE_OF_USE_televir_global)
            .values_list("technology__name", flat=True)
            .distinct()
        )

        technologies_available = list(set(technologies_available))

        return technologies_available

    def get_software_list(self):
        """
        Get software list available for a user
        """
        software_list = (
            Software.objects.filter(type_of_use=Software.TYPE_OF_USE_televir_global)
            .values_list("name", flat=True)
            .distinct()
        )
        software_list = list(set(software_list))
        return software_list

    @staticmethod
    def expand_parameters_table(combined_table, software_db_dict={}):
        def fix_row(row):
            if not row.parameter:
                return [""]

            if not row.can_change:
                return [row.parameter]

            if row.parameter_name == "--db" and software_db_dict:
                software_name = row.software_name
                possibilities = [software_name, software_name.lower()]
                if "_" in software_name:
                    possibilities.append(software_name.split("_")[0])

                for p in possibilities:
                    if p in software_db_dict.keys():
                        new_range = software_db_dict[p]
                        return new_range

            if not row.range_available:
                return [row.parameter]

            else:
                row_type = row.type_data
                range_min = row.range_min
                range_max = row.range_max
                new_range = row.parameter

                if row_type == Parameter.PARAMETER_int:
                    range_step = 1
                    range_min = int(range_min)
                    range_max = int(range_max)
                    range_step = int(range_step)
                    new_range = [
                        str(x) for x in range(range_min, range_max, range_step)
                    ]
                elif row_type == Parameter.PARAMETER_float:
                    range_step = 0.1
                    range_min = float(range_min)
                    range_max = float(range_max)
                    range_step = float(range_step)
                    new_range = [
                        str(round(x, 2))
                        for x in np.arange(range_min, range_max, range_step)
                    ]

                return new_range

        combined_table["parameter"] = combined_table.apply(fix_row, axis=1)
        combined_table = combined_table.reset_index(drop=True)

        combined_table = combined_table.explode("parameter")

        return combined_table

    def get_user_active_software_tables(self, technology: str, owner: User):
        """
        Get software tables for a user
        """

        software_available = Software.objects.filter(
            type_of_use__in=Software.TELEVIR_GLOBAL_TYPES,
            technology__name=technology,
            is_to_run=True,
            owner=owner,
        )

        parameters_available = Parameter.objects.filter(
            software__in=software_available,
            televir_project=None,
            televir_project_sample=None,
        )

        software_table = pd.DataFrame(software_available.values())

        parameters_table = pd.DataFrame(parameters_available.values())

        return software_table, parameters_table

    def get_software_tables(
        self,
        technology: str,
        user: User,
        project: Optional[Projects] = None,
        sample: Optional[PIProject_Sample] = None,
        metagenomics: bool = False,
        mapping_only: bool = False,
        screening: bool = False,
        request_mapping: bool = False,
    ):
        """
        Get software tables for a user
        """

        if metagenomics:
            steps = CS.vect_pipeline_televir_metagenomics
        elif mapping_only:
            steps = CS.vect_pipeline_televir_mapping_only
        elif screening:
            steps = CS.vect_pipeline_televir_screening
        elif request_mapping:
            steps = CS.vect_pipeline_televir_request_mapping
        else:
            steps = CS.vect_pipeline_televir_classic

        software_available = Software.objects.filter(
            technology__name=technology,
            pipeline_step__name__in=steps,
            is_to_run=True,
            owner=user,
        ).distinct()

        if not project and not sample:
            software_available = software_available.filter(
                type_of_use__in=Software.TELEVIR_GLOBAL_TYPES
            )

        elif sample is not None:
            software_available = software_available.filter(
                parameter__televir_project_sample=sample
            )

        elif project is not None:
            software_available = software_available.filter(
                parameter__televir_project=project,
                type_of_use__in=Software.TELEVIR_PROJECT_TYPES,
            )

        parameters_available = Parameter.objects.filter(
            software__in=software_available,
        ).distinct()

        software_table = pd.DataFrame(software_available.values())

        parameters_table = pd.DataFrame(parameters_available.values())

        return software_table, parameters_table

    def get_software_tables_project(self, owner: User, project: Projects):
        """
        Get software tables for a user
        """

        software_available = Software.objects.filter(
            owner=owner,
            type_of_use__in=Software.TELEVIR_PROJECT_TYPES,
            technology__name=project.technology,
            parameter__televir_project_sample=None,
            pipeline_step__name__in=self.televir_constants.vect_pipeline_names_default,
            is_to_run=True,
        )

        parameters_available = Parameter.objects.filter(
            software__in=software_available,
            televir_project=project,
            televir_project_sample=None,
            is_to_run=True,
        )

        software_table = pd.DataFrame(software_available.values())

        parameters_table = pd.DataFrame(parameters_available.values())

        return software_table, parameters_table

    def get_software_tables_project_sample_metagenomics(
        self, owner: User, project_sample: PIProject_Sample
    ):
        """
        Get software tables for a user
        """

        software_available = Software.objects.filter(
            owner=owner,
            type_of_use__in=Software.TELEVIR_PROJECT_TYPES,
            technology__name=project_sample.project.technology,
            parameter__televir_project_sample=project_sample,
            pipeline_step__name__in=CS.vect_pipeline_televir_metagenomics,
            is_to_run=True,
        )

        parameters_available = Parameter.objects.filter(
            software__in=software_available,
            televir_project=project_sample.project,
            televir_project_sample=project_sample,
            is_to_run=True,
        )

        software_table = pd.DataFrame(software_available.values())

        parameters_table = pd.DataFrame(parameters_available.values())

        return software_table, parameters_table

    def merge_software_tables(
        self, software_table: pd.DataFrame, parameters_table: pd.DataFrame
    ):
        """"""

        combined_table = pd.merge(
            software_table, parameters_table, left_on="id", right_on="software_id"
        )

        combined_table = combined_table.rename(
            columns={
                "id_x": "software_id",
                "id_y": "parameter_id",
                "name": "software_name",
                "name_x": "software_name",
                "name_y": "parameter_name",
                "is_to_run_x": "software_is_to_run",
                "is_to_run_y": "parameter_is_to_run",
            }
        )

        combined_table = combined_table[
            combined_table.type_of_use.isin(
                Software.TELEVIR_GLOBAL_TYPES + Software.TELEVIR_PROJECT_TYPES
            )
        ]

        combined_table["pipeline_step"] = combined_table["pipeline_step_id"].apply(
            lambda x: PipelineStep.objects.get(id=int(x)).name
        )

        combined_table["technology"] = combined_table["technology_id"].apply(
            lambda x: Technology.objects.get(id=int(x)).name
        )

        combined_table = combined_table.reset_index(drop=True)
        software_names = combined_table["software_name"].values
        can_change = combined_table["can_change"].values

        ## remove duplicate columns
        #
        combined_table = combined_table.loc[
            :, ~combined_table.T.duplicated(keep="last")
        ]

        combined_table["software_name"] = software_names
        combined_table["can_change"] = can_change

        return combined_table

    def generate_merged_table_safe(
        self,
        owner: User,
        technology: str,
        project: Optional[Projects] = None,
        sample: Optional[PIProject_Sample] = None,
        metagenomics: bool = False,
        mapping_only: bool = False,
        screening: bool = False,
        request_mapping: bool = False,
    ) -> pd.DataFrame:
        """
        Generate a software tree for a technology and a tree makeup"""

        software_table, parameters_table = self.get_software_tables(
            technology,
            owner,
            project=project,
            sample=sample,
            metagenomics=metagenomics,
            mapping_only=mapping_only,
            screening=screening,
            request_mapping=request_mapping,
        )

        if parameters_table.shape[0] == 0 or software_table.shape[0] == 0:
            if sample is not None:
                technology = sample.project.technology
            elif project is not None:
                technology = project.technology

            software_table, parameters_table = self.get_software_tables(
                technology,
                owner,
                project=project,
                sample=None,
                metagenomics=metagenomics,
                mapping_only=mapping_only,
                screening=screening,
                request_mapping=request_mapping,
            )

        if parameters_table.shape[0] == 0 or software_table.shape[0] == 0:
            if sample is not None:
                technology = sample.project.technology
            elif project is not None:
                technology = project.technology

            software_table, parameters_table = self.get_software_tables(
                technology,
                owner,
                project=None,
                sample=None,
                metagenomics=metagenomics,
                mapping_only=mapping_only,
                screening=screening,
                request_mapping=request_mapping,
            )

        if parameters_table.shape[0] == 0 or software_table.shape[0] == 0:
            return pd.DataFrame(
                columns=[
                    "software_id",
                    "parameter_id",
                    "technology",
                    "can_change",
                    "pipeline_step",
                    "software_name",
                ]
            )

        merged_table = self.merge_software_tables(software_table, parameters_table)

        return merged_table

    @staticmethod
    def convert_softwaretree_to_pipeline_tree(
        software_tree: SoftwareTree,
    ) -> PipelineTree:
        tree_nodes = SoftwareTreeNode.objects.filter(software_tree=software_tree)

        edges = []
        nodes = []
        leaves = []
        for node in tree_nodes:
            if node.parent:
                edges.append((node.parent.index, node.index))
            nodes.append((node.index, (node.name, node.value, node.node_type)))
            if node.node_place == 1:
                leaves.append(node.index)

        return PipelineTree(
            technology=software_tree.technology,
            nodes=[x[1] for x in sorted(nodes)],
            edges=edges,
            leaves=leaves,
            makeup=software_tree.global_index,
        )

    def retrace_from_leaf(self, leaf: SoftwareTreeNode) -> pd.DataFrame:
        """ """

        software_tree = leaf.software_tree
        parent = leaf.parent
        path = [(leaf.index, leaf.name, leaf.value, leaf.node_type)]
        while parent is not None:
            path.append((parent.index, parent.name, parent.value, parent.node_type))
            parent = parent.parent

        path = path[::-1]

        path_df = pd.DataFrame(path, columns=["index", "name", "value", "node_type"])

        return path_df

    def check_parameter_set_contains_module(
        self, leaf: SoftwareTreeNode, module: str
    ) -> bool:
        """
        retrace parameter set (leaf) settings,
        check if module is present in the path
        """

        path_df = self.retrace_from_leaf(leaf)
        if module in path_df.name.tolist():
            return True
        else:
            return False

    def check_ParameterSet_exists(
        self, sample: PIProject_Sample, leaf: SoftwareTreeNode, project: Projects
    ):
        self.logger.info("Checking if ParameterSet exists")
        return ParameterSet.objects.filter(
            sample=sample, leaf=leaf, project=project
        ).exists()

    def check_ParameterSet_available(
        self, sample: PIProject_Sample, leaf: SoftwareTreeNode, project: Projects
    ):
        if not self.check_ParameterSet_exists(sample, leaf, project):
            return True

        parameter_set = ParameterSet.objects.get(
            sample=sample, leaf=leaf, project=project
        )

        if parameter_set.status in [
            ParameterSet.STATUS_QUEUED,
            ParameterSet.STATUS_FINISHED,
            ParameterSet.STATUS_RUNNING,
        ]:
            return False

        return True

    def parameterset_update_status(
        self,
        sample: PIProject_Sample,
        leaf: SoftwareTreeNode,
        project: Projects,
        status: int,
    ):
        if not self.check_ParameterSet_exists(sample, leaf, project):
            return False

        parameter_set = ParameterSet.objects.get(
            sample=sample, leaf=leaf, project=project
        )

        parameter_set.status = status
        parameter_set.save()

        return True

    def check_ParameterSet_killed(
        self, sample: PIProject_Sample, leaf: SoftwareTreeNode, project: Projects
    ) -> bool:
        if not self.check_ParameterSet_exists(sample, leaf, project):
            return False

        parameter_set_killed = ParameterSet.objects.filter(
            sample=sample,
            leaf=leaf,
            project=project,
            status__in=[ParameterSet.STATUS_KILLED],
        ).exists()

        return parameter_set_killed

    def check_ParameterSet_available_to_run(
        self, sample: PIProject_Sample, leaf: SoftwareTreeNode, project: Projects
    ):
        if not self.check_ParameterSet_exists(sample, leaf, project):
            return True

        parameter_set = ParameterSet.objects.get(
            sample=sample, leaf=leaf, project=project
        )

        if parameter_set.status in [
            ParameterSet.STATUS_FINISHED,
            ParameterSet.STATUS_RUNNING,
        ]:
            return False

        return True

    def create_parameter_set(
        self, sample: PIProject_Sample, leaf: SoftwareTreeNode, project: Projects
    ):
        """
        Create a ParameterSet for a sample and leaf
        """
        self.logger.info("Creating ParameterSet")
        parameter_set = ParameterSet.objects.create(
            sample=sample, leaf=leaf, project=project
        )

    def check_ParameterSet_processed(
        self, sample: PIProject_Sample, leaf: SoftwareTreeNode, project: Projects
    ):
        """
        Check if ParameterSet is finished or running or queued.
        """

        if not self.check_ParameterSet_exists(sample, leaf, project):
            return False

        parameter_set = ParameterSet.objects.get(
            sample=sample, leaf=leaf, project=project
        )

        if parameter_set.status in [
            ParameterSet.STATUS_FINISHED,
            ParameterSet.STATUS_RUNNING,
            ParameterSet.STATUS_QUEUED,
        ]:
            return True
        else:
            return False

    def set_parameterset_to_queue(
        self, sample: PIProject_Sample, leaf: SoftwareTreeNode, project: Projects
    ):
        """
        Set ParameterSet to queue if it exists and is not finished or running.
        """

        try:
            parameter_set = ParameterSet.objects.get(
                sample=sample, leaf=leaf, project=project
            )

            if parameter_set.status not in [
                ParameterSet.STATUS_FINISHED,
                ParameterSet.STATUS_RUNNING,
            ]:
                parameter_set.status = ParameterSet.STATUS_QUEUED
                parameter_set.save()

        except ParameterSet.DoesNotExist:
            pass


class Utils_Manager:
    """Combines Utility classes to create a manager for the pipeline."""

    utilities_repository: Utility_Repository
    parameter_util: Parameter_DB_Utility
    utility_manager: Utility_Pipeline_Manager

    def __init__(self):
        ###
        self.parameter_util = Parameter_DB_Utility()

        self.utility_repository = Utility_Repository(
            db_path=Televir_Directories.docker_app_directory,
            install_type="docker",
        )

        self.utility_technologies = self.parameter_util.get_technologies_available()
        self.utility_manager = Utility_Pipeline_Manager()
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.ERROR)
        if self.logger.hasHandlers():
            self.logger.handlers.clear()

        self.logger.info("Utils_Manager initialized")

    def get_leaf_parameters(self, parameter_leaf: SoftwareTreeNode) -> pd.DataFrame:
        """ """
        pipeline_tree = self.parameter_util.convert_softwaretree_to_pipeline_tree(
            parameter_leaf.software_tree
        )

        if parameter_leaf.index not in pipeline_tree.leaves:
            raise Exception("Node is not a leaf")

        all_paths = pipeline_tree.get_all_graph_paths()

        return all_paths[parameter_leaf.index]

    def get_parameterset_leaves(
        self, parameterset: ParameterSet, pipeline_tree: PipelineTree
    ) -> list:
        """
        retrieve list of leaves for a parameterset, matched to a given pipeline tree explicitely (using full paths).
        """

        ps_pipeline_tree = self.parameter_util.convert_softwaretree_to_pipeline_tree(
            parameterset.leaf.software_tree
        )
        ps_leaves = ps_pipeline_tree.leaves_from_node(parameterset.leaf.index)

        ps_paths = ps_pipeline_tree.get_specific_leaf_paths_explicit(ps_leaves)

        new_matched_paths = {
            leaf: self.utility_manager.match_path_to_tree_safe(path, pipeline_tree)
            for leaf, path in ps_paths.items()
        }

        new_matched_paths = {
            k: v for k, v in new_matched_paths.items() if v is not None
        }

        return list(new_matched_paths.values())

    def collect_project_samples(self, project: Projects) -> dict:
        """
        Collect all samples from a project
        """
        samples = PIProject_Sample.objects.filter(project=project)
        submission_dict = {sample: [] for sample in samples if not sample.is_deleted}
        return submission_dict

    def sample_nodes_check_no_repeats(
        self, submission_dict: dict, available_path_nodes: dict, project: Projects
    ):
        utils = Utils_Manager()
        ### SUBMISSION
        runs_to_deploy = 0
        samples_available = []
        samples_leaf_dict = {sample: [] for sample in submission_dict.keys()}

        for sample in submission_dict.keys():
            for leaf, matched_path_node in available_path_nodes.items():
                exists = self.parameter_util.check_ParameterSet_exists(
                    sample=sample, leaf=matched_path_node, project=project
                )

                if exists:
                    available = utils.parameter_util.check_ParameterSet_available(
                        sample=sample, leaf=matched_path_node, project=project
                    )
                    if (
                        utils.parameter_util.check_ParameterSet_available(
                            sample=sample, leaf=matched_path_node, project=project
                        )
                        is False
                    ):
                        continue

                else:
                    self.parameter_util.create_parameter_set(
                        sample=sample, leaf=matched_path_node, project=project
                    )

                utils.parameter_util.set_parameterset_to_queue(
                    sample=sample, leaf=matched_path_node, project=project
                )
                runs_to_deploy += 1
                samples_available.append(sample)
                samples_leaf_dict[sample].append(matched_path_node)

        samples_leaf_dict = {x: g for x, g in samples_leaf_dict.items() if g}

        return samples_leaf_dict

    def sample_nodes_check_repeat_allowed(
        self, submission_dict: dict, available_path_nodes: dict, project: Projects
    ):
        utils = Utils_Manager()
        ### SUBMISSION
        runs_to_deploy = 0
        samples_available = []
        samples_leaf_dict = {sample: [] for sample in submission_dict.keys()}
        workflow_deployed_dict = {sample: {} for sample in submission_dict.keys()}

        for sample in submission_dict.keys():
            for leaf, matched_path_node in available_path_nodes.items():
                exists = self.parameter_util.check_ParameterSet_exists(
                    sample=sample, leaf=matched_path_node, project=project
                )

                if not exists:
                    self.parameter_util.create_parameter_set(
                        sample=sample, leaf=matched_path_node, project=project
                    )

                workflow_deployed = (
                    True
                    if utils.parameter_util.check_ParameterSet_available(
                        sample=sample, leaf=matched_path_node, project=project
                    )
                    is False
                    else False
                )

                utils.parameter_util.set_parameterset_to_queue(
                    sample=sample, leaf=matched_path_node, project=project
                )
                runs_to_deploy += 1
                samples_available.append(sample)
                samples_leaf_dict[sample].append(matched_path_node)
                workflow_deployed_dict[sample][matched_path_node] = workflow_deployed

        samples_leaf_dict = {x: g for x, g in samples_leaf_dict.items() if g}

        return samples_leaf_dict, workflow_deployed_dict

    def tree_subset(self, tree: PipelineTree, leaves: list) -> PipelineTree:
        """
        Return a subset of a tree
        """

        if len(leaves) == 0:
            return tree

        reduced_dag, reduced_node_index = tree.reduced_tree(leaves)

        reduced_tree = self.pipe_tree_from_dag_dict(
            reduced_dag, reduced_node_index, tree.technology, tree.makeup
        )

        reduced_tree.software_tree_pk = tree.software_tree_pk

        return reduced_tree

    def pipe_tree_from_dag_dict(
        self,
        dag_dict: dict,
        node_index: pd.DataFrame,
        technology: str,
        tree_makeup: int,
    ) -> PipelineTree:
        """
        Generate a pipeline tree from a dag dict
        """

        nodes = node_index.node.unique()

        nodes = []
        edge_list = []
        leaves = []
        for node in dag_dict:
            nodes.append(node)
            for child in dag_dict[node]:
                edge_list.append((node, child))

            if len(dag_dict[node]) == 0:
                leaves.append(node)

        return PipelineTree(
            nodes=node_index.reset_index().to_numpy().tolist(),
            edges=edge_list,
            leaves=leaves,
            technology=technology,
            makeup=tree_makeup,
            sorted=False,
        )

    ### Copied to softwareTreeUtils
    def check_pipeline_possible(self, combined_table: pd.DataFrame, tree_makeup: int):
        """
        Check if a pipeline is possible
        """

        pipeline_setup = Pipeline_Makeup()
        makeup_steps = pipeline_setup.get_makeup(tree_makeup)

        pipelines_available = combined_table.pipeline_step.unique().tolist()
        pipelines_available = [x for x in pipelines_available if x in makeup_steps]
        self.pipeline_makeup = pipeline_setup.match_makeup_name_from_list(
            pipelines_available
        )

        if not self.pipeline_makeup:
            return False

        return True

    def check_any_pipeline_possible(self, technology: str, user: User):
        """
        Check if a pipeline is possible
        """
        pipeline_setup = Pipeline_Makeup()

        combined_table = self.parameter_util.generate_merged_table_safe(
            user, technology
        )

        for makeup in pipeline_setup.get_makeup_list():
            if self.check_pipeline_possible(combined_table, makeup):
                return True

        return False

    def test_televir_pipelines_available(self, user_system: User):
        """
        Test if televir is available
        """

        software = Software.objects.filter(
            type_of_use=Software.TYPE_OF_USE_televir_global, owner=user_system
        )
        if software.count() == 0:
            return False

        for technology in self.utility_technologies:
            if self.check_any_pipeline_possible(technology, user_system):
                return True

        return False

    def module_tree(self, pipeline_tree: PipelineTree, leaves: List[int]):
        """
        Return a subset of a tree
        """

        reduced_tree = self.tree_subset(pipeline_tree, leaves)

        module_tree = self.utility_manager.compress_software_tree(reduced_tree)

        return module_tree


class SoftwareTreeUtils:
    def __init__(
        self,
        user: User,
        project: Optional[Projects] = None,
        sample: Optional[PIProject_Sample] = None,
    ):
        self.user = user
        self.project = project
        self.sample = sample
        if project:
            self.technology = project.technology
        else:
            self.technology = None
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.ERROR)

        self.utils_manager = Utils_Manager()
        self.parameter_util = Parameter_DB_Utility()
        self.utility_manager = Utility_Pipeline_Manager()

    ###############################################
    ###############################################  SOFTWARE TREE CONNECTIONS

    def set_project(self, project: Projects):
        self.project = project
        self.technology = project.technology

    def set_sample(self, sample: PIProject_Sample):
        self.sample = sample

    def query_software_tree(self, global_index: int) -> SoftwareTree:
        """
        Query software tree
        """

        try:
            software_tree = (
                SoftwareTree.objects.filter(
                    global_index=global_index,
                    technology=self.technology,
                    model=ConstantsSettings.PIPELINE_MODEL,
                    project=self.project,
                    owner=self.user,
                )
                .order_by("date_created")
                .last()
            )

        except SoftwareTree.DoesNotExist:
            software_tree = None
        return software_tree

    def check_default_software_tree_exists(
        self,
        global_index: int,
    ):
        try:
            software_tree = self.query_software_tree(
                global_index=global_index,
            )

            if software_tree:
                return True
            else:
                return False

        except SoftwareTree.DoesNotExist:
            return False

        except Exception as e:
            print(e)

    def get_software_tree_index(self, global_index: int) -> Optional[int]:
        """
        Get software tree index db
        """

        if self.check_default_software_tree_exists(global_index):
            software_tree = self.query_software_tree(
                global_index=global_index,
            )

            return software_tree.pk
        else:
            return None

    @staticmethod
    def software_pipeline_tree(software_tree: SoftwareTree) -> PipelineTree:
        tree_nodes = SoftwareTreeNode.objects.filter(software_tree=software_tree)

        edges = []
        nodes = []
        leaves = []
        for node in tree_nodes:
            if node.parent:
                edges.append((node.parent.index, node.index))
            nodes.append((node.index, (node.name, node.value, node.node_type)))
            if node.node_place == 1:
                leaves.append(node.index)

        return PipelineTree(
            technology=software_tree.technology,
            nodes=[x[1] for x in sorted(nodes)],
            edges=edges,
            leaves=leaves,
            makeup=software_tree.global_index,
        )

    def query_software_default_tree(
        self,
        global_index: int,
    ) -> PipelineTree:
        """
        Generate a default pipeline tree for a user
        """

        software_tree = self.query_software_tree(
            global_index=global_index,
        )

        if not software_tree:
            return None
        return self.software_pipeline_tree(software_tree)

    def update_software_tree(self, tree: PipelineTree):
        """
        Update SoftwareTree table
        """
        global_index = tree.makeup

        software_tree = self.query_software_tree(
            global_index=global_index,
        )

        if not software_tree:
            self.logger.info("Creating new software tree")
            software_tree = SoftwareTree(
                global_index=global_index,
                technology=tree.technology,
                version=0,
                model=ConstantsSettings.PIPELINE_MODEL,
                project=self.project,
                owner=self.user,
            )

            software_tree.save()

            tree.software_tree_pk = software_tree.pk

        self.update_SoftwareTree_nodes(software_tree, tree)

    def update_SoftwareTree_nodes(
        self, software_tree: SoftwareTree, tree: PipelineTree
    ):
        """
        Update the nodes of a software tree
        """

        parent_dict = tree.get_parents_dict()

        for index, row in tree.node_index.iterrows():
            index = int(index)
            node = row.node

            is_leaf = int(index in tree.leaves)
            name = node[0]
            value = node[1]
            node_type = node[2]

            try:
                tree_node = SoftwareTreeNode.objects.get(
                    software_tree=software_tree,
                    index=index,
                    name=name,
                    value=value,
                    node_type=node_type,
                )

            except SoftwareTreeNode.DoesNotExist:
                try:
                    with LockedAtomicTransaction(SoftwareTreeNode):
                        parent_node_index = parent_dict.get(index, None)
                        parent_node = None

                        if parent_node_index != None:
                            parent_node = tree.nodes[parent_node_index]
                            parent_name = parent_node[0]
                            parent_value = parent_node[1]
                            parent_type = parent_node[2]
                            parent_node = SoftwareTreeNode.objects.filter(
                                software_tree=software_tree,
                                index=parent_dict[index],
                                name=parent_name,
                                value=parent_value,
                                node_type=parent_type,
                            )
                            if parent_node.exists():
                                parent_node = parent_node.first()

                        tree_node = SoftwareTreeNode(
                            software_tree=software_tree,
                            index=index,
                            name=node[0],
                            value=node[1],
                            node_type=node[2],
                            parent=parent_node,
                            node_place=is_leaf,
                        )
                        tree_node.save()

                except Exception as e:
                    print(e)
                    print(index, name, value)

            except Exception as e:
                print(e)

    def generate_software_tree_safe(
        self,
        project: Projects,
        sample: Optional[PIProject_Sample] = None,
        metagenomics: bool = False,
        mapping_only: bool = False,
        screening: bool = False,
        request_mapping: bool = False,
    ) -> PipelineTree:
        """
        Generate a software tree for a technology and a tree makeup
        """

        merged_table = self.parameter_util.generate_merged_table_safe(
            project.owner,
            project.technology,
            project=project,
            sample=sample,
            metagenomics=metagenomics,
            mapping_only=mapping_only,
            screening=screening,
            request_mapping=request_mapping,
        )

        if merged_table.shape[0] == 0:
            return PipelineTree(
                technology=project.technology,
                nodes=[],
                edges={},
                leaves=[],
                makeup=-1,
            )

        return self.generate_tree_from_combined_table(merged_table)

    def generate_project_tree(self) -> PipelineTree:
        """
        Generate a project tree
        """

        return self.generate_software_tree_safe(self.project)

    def generate_tree_from_combined_table(
        self, combined_table: pd.DataFrame
    ) -> PipelineTree:
        utility_drone = Utility_Pipeline_Manager()
        input_success = utility_drone.input(combined_table, technology=self.technology)

        if not input_success:
            return PipelineTree(
                technology=self.technology,
                nodes=[],
                edges={},
                leaves=[],
                makeup=-1,
            )

        self.logger.info("Generating project tree")

        pipeline_tree = utility_drone.generate_default_software_tree()

        pipeline_tree = self.prep_tree_for_extend(pipeline_tree)

        return pipeline_tree

    def check_pipeline_possible(self, combined_table: pd.DataFrame, tree_makeup: int):
        """
        Check if a pipeline is possible
        """

        pipeline_setup = Pipeline_Makeup()
        makeup_steps = pipeline_setup.get_makeup(tree_makeup)

        pipelines_available = combined_table.pipeline_step.unique().tolist()
        pipelines_available = [x for x in pipelines_available if x in makeup_steps]
        self.pipeline_makeup = pipeline_setup.match_makeup_name_from_list(
            pipelines_available
        )

        if not self.pipeline_makeup:
            return False

        return True

    def check_any_pipeline_possible(self, technology: str, user: User):
        """
        Check if a pipeline is possible
        """
        pipeline_setup = Pipeline_Makeup()

        combined_table = self.parameter_util.generate_merged_table_safe(
            user, technology
        )

        for makeup in pipeline_setup.get_makeup_list():
            if self.check_pipeline_possible(combined_table, makeup):
                return True

        return False

    def test_televir_pipelines_available(self, user_system: User):
        """
        Test if televir is available
        """

        for technology in self.parameter_util.get_technologies_available():
            if self.check_any_pipeline_possible(technology, user_system):
                return True

        return False

    def get_project_pathnodes(self) -> dict:
        """
        Get all pathnodes for a project
        """
        # local_tree = self.generate_project_tree()
        if self.project is None:
            return {}
        local_tree = self.generate_software_tree_safe(self.project)

        return self.get_available_pathnodes(local_tree)

    def get_sample_pathnodes(
        self,
        metagenomics: bool = False,
        mapping_only: bool = False,
        screening: bool = False,
    ) -> dict:
        """
        Get all pathnodes for a project
        """

        local_tree = self.generate_software_tree_safe(
            self.project,
            self.sample,
            metagenomics=metagenomics,
            mapping_only=mapping_only,
            screening=screening,
        )

        if local_tree.makeup == -1:
            return {}

        return self.get_available_pathnodes(local_tree)

    def get_available_nodes_summary(
        self,
        metagenomics: bool = False,
        mapping_only: bool = False,
        screening: bool = False,
    ) -> dict:
        """return path as df for each leaf"""

        local_tree = self.generate_software_tree_safe(
            self.project,
            self.sample,
            metagenomics=metagenomics,
            mapping_only=mapping_only,
            screening=screening,
        )

        all_paths = local_tree.get_all_graph_paths()
        return all_paths

    def get_available_pathnodes(
        self, local_tree: PipelineTree
    ) -> Dict[int, SoftwareTreeNode]:
        """ """

        # pipeline_tree = utils.generate_software_tree(technology, tree_makeup)
        utils = Utils_Manager()

        local_paths = local_tree.get_all_graph_paths_explicit()
        pipeline_tree = self.generate_software_tree_extend(local_tree=local_tree)
        ### MANAGEMENT
        matched_paths = {
            leaf: utils.utility_manager.match_path_to_tree_safe(path, pipeline_tree)
            for leaf, path in local_paths.items()
        }

        available_paths = {
            leaf: path for leaf, path in matched_paths.items() if path is not None
        }

        available_path_nodes = {
            leaf: SoftwareTreeNode.objects.get(
                software_tree__pk=local_tree.software_tree_pk, index=leaf_index
            )
            for leaf, leaf_index in available_paths.items()
        }

        return available_path_nodes

    def check_runs_to_deploy_project(self) -> dict:
        """
        Check if there are runs to run. sets to queue if there are.
        """

        submission_dict = self.utils_manager.collect_project_samples(self.project)

        available_path_nodes = self.get_project_pathnodes()
        clean_samples_leaf_dict = self.utils_manager.sample_nodes_check_no_repeats(
            submission_dict, available_path_nodes, self.project
        )

        return clean_samples_leaf_dict

    def check_runs_to_submit_metagenomics_sample(
        self, sample: PIProject_Sample
    ) -> dict:
        """
        Check if there are runs to run. sets to queue if there are.
        """

        submission_dict = {sample: []}

        available_path_nodes = self.get_sample_pathnodes(
            metagenomics=True,
            screening=False,
            mapping_only=False,
        )

        clean_samples_leaf_dict, workflow_deployed_dict = (
            self.utils_manager.sample_nodes_check_repeat_allowed(
                submission_dict, available_path_nodes, self.project
            )
        )

        # samples_leaf_dict = {
        #    sample: available_path_nodes.values() for sample in submission_dict.keys()
        # }

        return clean_samples_leaf_dict

    def check_runs_to_submit_screening_sample(self, sample: PIProject_Sample) -> dict:
        """
        Check if there are runs to run. sets to queue if there are.
        """

        submission_dict = {sample: []}

        available_path_nodes = self.get_sample_pathnodes(
            metagenomics=False,
            screening=True,
            mapping_only=False,
        )

        clean_samples_leaf_dict, workflow_deployed_dict = (
            self.utils_manager.sample_nodes_check_repeat_allowed(
                submission_dict, available_path_nodes, self.project
            )
        )

        return clean_samples_leaf_dict

    def check_runs_to_submit_mapping_only(
        self, sample: PIProject_Sample
    ) -> Tuple[dict, dict]:
        """
        Check if there are runs to run. sets to queue if there are.
        """

        submission_dict = {sample: []}

        available_path_nodes = self.get_sample_pathnodes(
            metagenomics=False,
            screening=False,
            mapping_only=True,
        )

        clean_samples_leaf_dict, workflow_deployed_dict = (
            self.utils_manager.sample_nodes_check_repeat_allowed(
                submission_dict, available_path_nodes, self.project
            )
        )

        return clean_samples_leaf_dict, workflow_deployed_dict

    def check_runs_to_deploy_sample(self, sample: PIProject_Sample) -> dict:
        """
        Check if there are runs to run. sets to queue if there are.
        """

        submission_dict = {sample: []}

        available_path_nodes = self.get_sample_pathnodes(
            metagenomics=False,
            screening=False,
            mapping_only=False,
        )
        clean_samples_leaf_dict = self.utils_manager.sample_nodes_check_no_repeats(
            submission_dict, available_path_nodes, self.project
        )

        return clean_samples_leaf_dict

    def get_all_technology_pipelines(self, tree_makeup: int) -> dict:
        """
        Get all pipelines for a technology
        """

        pipeline_tree = self.generate_software_tree(tree_makeup)

        all_paths = pipeline_tree.get_all_graph_paths()
        return all_paths

    def generate_software_tree(self, tree_makeup: int):
        """
        Generate a software tree for a technology and a tree makeup
        """

        if self.check_default_software_tree_exists(global_index=tree_makeup):
            return self.query_software_default_tree(global_index=tree_makeup)
        else:
            raise Exception("No software tree for technology")

    def generate_software_tree_register(self, local_tree: PipelineTree):
        """
        Generate a software tree for a technology and a tree makeup
        """
        tree_makeup = local_tree.makeup
        technology = local_tree.technology

        if self.check_default_software_tree_exists(global_index=tree_makeup) is False:
            self.update_software_tree(local_tree)

        pipeline_tree = self.query_software_default_tree(global_index=tree_makeup)

        if len(pipeline_tree.nodes) == 0:
            self.update_software_tree(local_tree)

        return pipeline_tree

    def generate_software_tree_extend(self, local_tree: PipelineTree):
        """Generate Software Tree Register and extend with local paths"""
        local_paths = local_tree.get_all_graph_paths_explicit()
        pipeline_tree = self.generate_software_tree_register(local_tree)
        for leaf, path in local_paths.items():
            pipeline_tree = self.utility_manager.match_path_to_tree_extend(
                path, pipeline_tree
            )
        self.update_software_tree(pipeline_tree)
        pipeline_tree = self.prep_tree_for_extend(pipeline_tree)
        return pipeline_tree

    def prep_tree_for_extend(self, tree: PipelineTree):
        tree.software_tree_pk = self.get_software_tree_index(
            tree.makeup,
        )
        return tree


class RawReferenceUtils:
    def __init__(
        self,
        sample: Optional[PIProject_Sample] = None,
        project: Optional[Projects] = None,
    ):
        self.sample_registered = sample
        self.project_registered = project
        self.runs_found = 0
        self.list_tables: List[pd.DataFrame] = []
        self.merged_table: pd.DataFrame = pd.DataFrame(
            columns=["read_counts", "contig_counts", "taxid", "accid", "description"]
        )

    def references_table_from_query(
        self, references: Union[QuerySet, List[RawReference]]
    ) -> pd.DataFrame:
        table = []
        for ref in references:
            table.append(
                {
                    "taxid": ref.taxid,
                    "accid": ref.accid,
                    "description": ref.description,
                    "counts_str": ref.counts,
                    "read_counts": ref.read_counts,
                    "contig_counts": ref.contig_counts,
                }
            )

        if len(table) == 0:
            return pd.DataFrame(
                columns=[
                    "taxid",
                    "accid",
                    "description",
                    "counts_str",
                    "read_counts",
                    "contig_counts",
                ]
            )

        references_table = pd.DataFrame(table)

        references_table["read_counts"] = references_table["read_counts"].astype(float)
        references_table["contig_counts"] = references_table["contig_counts"].astype(
            float
        )

        # references_table = references_table[references_table["read_counts"] > 1]
        references_table = references_table[references_table["accid"] != "-"]
        references_table = references_table[references_table["accid"] != ""]
        references_table = references_table.sort_values(
            ["contig_counts", "read_counts"], ascending=[False, False]
        )

        references_table["sort_rank"] = range(1, references_table.shape[0] + 1)

        return references_table

    def run_references_standard_score_reads(
        self,
        references_table: pd.DataFrame,
    ) -> pd.DataFrame:

        if references_table.shape[0] == 0:
            return pd.DataFrame(
                columns=list(references_table.columns) + ["read_counts_standard_score"]
            )

        if max(references_table["read_counts"]) == 0:
            references_table["read_counts_standard_score"] = 1
            return references_table

        references_table["read_counts"] = references_table["read_counts"].astype(float)

        references_table["read_counts_standard_score"] = (
            references_table["read_counts"] - references_table["read_counts"].mean()
        ) / references_table["read_counts"].std()

        references_table["read_counts_standard_score"] = (
            references_table["read_counts_standard_score"]
            - references_table["read_counts_standard_score"].min()
        ) / (
            references_table["read_counts_standard_score"].max()
            - references_table["read_counts_standard_score"].min()
        )

        return references_table

    def run_references_standard_score_contigs(
        self,
        references_table: pd.DataFrame,
    ) -> pd.DataFrame:
        # references = RawReference.objects.filter(run=run)

        # references_table = references_table_from_query(references)

        if references_table.shape[0] == 0:
            return pd.DataFrame(
                columns=list(references_table.columns)
                + ["contig_counts_standard_score"]
            )

        references_table["contig_counts"] = references_table["contig_counts"].astype(
            float
        )
        if max(references_table["contig_counts"]) == 0:
            references_table["contig_counts_standard_score"] = 1
            return references_table

        references_table["contig_counts_standard_score"] = (
            references_table["contig_counts"] - references_table["contig_counts"].mean()
        ) / references_table["contig_counts"].std()

        references_table["contig_counts_standard_score"] = (
            references_table["contig_counts_standard_score"]
            - references_table["contig_counts_standard_score"].min()
        ) / (
            references_table["contig_counts_standard_score"].max()
            - references_table["contig_counts_standard_score"].min()
        )

        return references_table

    def merge_standard_scores(self, table: pd.DataFrame):
        if table.shape[0] == 0:
            return pd.DataFrame(columns=list(table.columns) + ["standard_score"])
        table["standard_score"] = table[
            "read_counts_standard_score"
        ]  # + table["contig_counts_standard_score"]
        return table

    def run_references_standard_scores(self, table):
        table = self.run_references_standard_score_reads(table)

        table = self.run_references_standard_score_contigs(table)
        table = self.merge_standard_scores(table)
        return table

    def merge_ref_tables_use_standard_score(
        self,
        list_tables: List[pd.DataFrame],
    ) -> pd.DataFrame:
        joint_tables = [
            self.run_references_standard_scores(table) for table in list_tables
        ]

        joint_tables = pd.concat(joint_tables)
        # group tables: average read_counts_standard_score, sum counts, read_counts, contig_counts
        if joint_tables.shape[0] == 0:
            return pd.DataFrame(columns=list(joint_tables.columns))

        joint_tables["standard_score"] = joint_tables["standard_score"].astype(float)
        joint_tables["contig_counts"] = joint_tables["contig_counts"].astype(float)
        joint_tables["read_counts"] = joint_tables["read_counts"].astype(float)
        joint_tables["contig_counts_standard_score"] = joint_tables[
            "contig_counts_standard_score"
        ].astype(float)
        joint_tables["sort_rank"] = joint_tables["sort_rank"].astype(float)

        joint_tables = joint_tables.groupby(["taxid"]).agg(
            {
                "taxid": "first",
                "accid": "first",
                "description": "first",
                "counts_str": "first",
                "standard_score": "mean",
                "contig_counts_standard_score": "mean",
                "read_counts": "sum",
                "contig_counts": "sum",
                "sort_rank": "mean",
            }
        )

        joint_tables = joint_tables.rename(columns={"sort_rank": "ensemble_ranking"})

        #############################################
        proxy_rclass = self.reference_table_renamed(
            joint_tables, {"read_counts": "counts"}
        ).reset_index(drop=True)
        proxy_aclass = self.reference_table_renamed(
            joint_tables, {"contig_counts": "counts"}
        ).reset_index(drop=True)

        targets, raw_targets = merge_classes(
            proxy_rclass, proxy_aclass, maxt=joint_tables.shape[0]
        )

        targets["global_ranking"] = range(1, targets.shape[0] + 1)

        def set_global_ranking_repeat_ranks(
            targets, rank_column="global_ranking", counts_column="counts"
        ):
            """
            Set the global ranking for repeated ranks
            """
            current_counts = None
            current_rank = 0
            for row in targets.iterrows():
                if current_counts != row[1][counts_column]:
                    current_rank += 1
                    current_counts = row[1][counts_column]

                targets.at[row[0], rank_column] = current_rank

            return targets

        targets = set_global_ranking_repeat_ranks(targets, rank_column="global_ranking")

        ####
        joint_tables = joint_tables.reset_index(drop=True)
        targets = targets.reset_index(drop=True)
        joint_tables["taxid"] = joint_tables["taxid"].astype(int)
        targets["taxid"] = targets["taxid"].astype(int)

        joint_tables = joint_tables.merge(
            targets[["taxid", "global_ranking"]], on=["taxid"], how="left"
        )
        ############################################# Reset the index
        joint_tables = joint_tables.reset_index(drop=True)

        joint_tables = joint_tables.sort_values("ensemble_ranking", ascending=True)

        joint_tables = joint_tables.reset_index(drop=True)

        return joint_tables

    def merge_ref_tables_use_ranking(
        self,
        list_tables: List[pd.DataFrame],
    ) -> pd.DataFrame:
        joint_tables = [
            self.run_references_standard_scores(table) for table in list_tables
        ]

        joint_tables = pd.concat(joint_tables)
        # group tables: average read_counts_standard_score, sum counts, read_counts, contig_counts
        if joint_tables.shape[0] == 0:
            return pd.DataFrame(columns=list(joint_tables.columns))

        joint_tables["standard_score"] = joint_tables["standard_score"].astype(float)
        joint_tables["contig_counts"] = joint_tables["contig_counts"].astype(float)
        joint_tables["read_counts"] = joint_tables["read_counts"].astype(float)
        joint_tables["contig_counts_standard_score"] = joint_tables[
            "contig_counts_standard_score"
        ].astype(float)

        joint_tables = joint_tables.groupby(["taxid", "accid", "description"]).agg(
            {
                "taxid": "first",
                "accid": "first",
                "description": "first",
                "counts_str": "first",
                "standard_score": "mean",
                "contig_counts_standard_score": "mean",
                "read_counts": "sum",
                "contig_counts": "sum",
                "sort_rank": "mean",
            }
        )

        # Define a function to calculate the final score
        def calculate_final_score(row):
            boost = 0
            if row["contig_counts"] > 0:
                boost = 1  # Define the boost value according to your needs
            return (
                row["standard_score"]
                + boost * row["contig_counts_standard_score"]
                + boost
            )

        # Apply the function to each row
        joint_tables["final_score"] = joint_tables.apply(calculate_final_score, axis=1)

        # Sort the table by the final score
        joint_tables = joint_tables.sort_values("final_score", ascending=False)

        # Reset the index
        joint_tables = joint_tables.reset_index(drop=True)

        joint_tables = joint_tables.sort_values(
            ["contig_counts", "standard_score"], ascending=[False, False]
        )
        joint_tables = joint_tables.reset_index(drop=True)

        return joint_tables

    def run_references_table(self, run: RunMain) -> pd.DataFrame:
        references = RawReference.objects.filter(run=run)

        references_table = self.references_table_from_query(references)

        return references_table

    def filter_runs(self):
        if self.sample_registered is None and self.project_registered is None:
            raise Exception("No sample or project registered")

        if self.sample_registered is None:
            sample_runs = RunMain.objects.filter(
                sample__project=self.project_registered
            )
        else:
            sample_runs = RunMain.objects.filter(sample=self.sample_registered)

        return sample_runs.exclude(
            run_type__in=[RunMain.RUN_TYPE_SCREENING, RunMain.RUN_TYPE_STORAGE]
        )

    def collect_references_all(self) -> QuerySet:
        sample_runs = self.filter_runs()

        references = RawReference.objects.filter(run__in=sample_runs)

        return references

    def sample_compound_refs_table(self) -> pd.DataFrame:

        compound_refs = RawReferenceCompoundModel.objects.filter(
            sample=self.sample_registered
        )

        compound_refs_table = pd.DataFrame(
            list(
                compound_refs.values(
                    "taxid",
                    "accid",
                    "description",
                    "ensemble_ranking",
                    "standard_score",
                )
            )
        )

        return compound_refs_table

    def sample_reference_tables(
        self, run_pks: Optional[List[int]] = None
    ) -> pd.DataFrame:
        sample_runs = self.filter_runs()
        if run_pks is not None:
            sample_runs = sample_runs.filter(pk__in=run_pks)
        self.runs_found = sample_runs.count()

        run_references_tables = [self.run_references_table(run) for run in sample_runs]

        # register tables
        self.list_tables.extend(run_references_tables)
        #
        if len(self.list_tables):
            run_references_tables = self.merge_ref_tables()
            # replace nan with 0
            run_references_tables = run_references_tables.fillna(0)
            self.merged_table = run_references_tables
        else:
            return pd.DataFrame(
                columns=[
                    "read_counts",
                    "contig_counts",
                    "taxid",
                    "accid",
                    "description",
                ]
            )

    def sample_reference_tables_filter(self, runs_filter: Optional[List[int]] = None):
        """
        Filter the sample reference tables to only include runs in the list
        """
        if runs_filter == []:
            runs_filter = None

        _ = self.sample_reference_tables(runs_filter)

    def compound_reference_update_standard_score(
        self, compound_ref: RawReferenceCompound
    ):
        """
        Update the standard score for a compound reference based on accid"""

        score = self.merged_table[self.merged_table.accid == compound_ref.accid]

        if score.shape[0] > 0:
            # compound_ref.standard_score = score.iloc[0]["standard_score"]
            compound_ref.standard_score = max(score["standard_score"])
            compound_ref.global_ranking = min(score["global_ranking"])
            compound_ref.ensemble_ranking = min(score["ensemble_ranking"])

    def update_scores_compound_references(
        self, compount_refs: List[RawReferenceCompound]
    ):
        """
        Update the standard score for a list of compound references based on accid"""
        for compound_ref in compount_refs:
            self.compound_reference_update_standard_score(compound_ref)

    def filter_reference_query_set(
        self, references: QuerySet, query_string: Optional[str] = ""
    ):
        """
        Filter a query set of references by a query string
        """
        if not query_string:
            return references.exclude(accid="-")

        references_select = references.filter(
            Q(description__icontains=query_string)
            | Q(accid__icontains=query_string)
            | Q(taxid__icontains=query_string)
        ).exclude(accid="-")

        return references_select

    def filter_reference_query_set_compound(
        self, references: QuerySet, query_string: Optional[str] = ""
    ):
        """
        Filter a query set of references by a query string
        """
        references_select = self.filter_reference_query_set(references, query_string)

        exclude_refs = []
        for ref in references_select:
            if RawReference.objects.filter(pk=ref.selected_mapped_pk).exists() is False:
                exclude_refs.append(ref.pk)

        references_select = references_select.exclude(pk__in=exclude_refs)

        return references_select

    def query_sample_compound_references(
        self, query_string: Optional[str] = None
    ) -> List[RawReferenceCompoundModel]:

        if self.sample_registered is not None:
            query_set = RawReferenceCompoundModel.objects.filter(
                sample=self.sample_registered
            ).order_by("ensemble_ranking")
        elif self.project_registered is not None:
            query_set = RawReferenceCompoundModel.objects.filter(
                sample__project=self.project_registered
            ).order_by("ensemble_ranking")
        else:
            query_set = RawReferenceCompoundModel.objects.none()

        return self.filter_reference_query_set_compound(query_set, query_string)

    def query_sample_references(self, query_string: Optional[str] = "") -> QuerySet:

        if self.sample_registered is not None:

            query_set = (
                RawReference.objects.filter(
                    run__sample__pk=self.sample_registered.pk,
                )
                .exclude(
                    run__run_type__in=[
                        RunMain.RUN_TYPE_STORAGE,
                        RunMain.RUN_TYPE_SCREENING,
                    ],
                    accid="-",
                )
                .distinct("accid")
            )
        elif self.project_registered is not None:
            query_set = (
                RawReference.objects.filter(
                    run__sample__project__pk=self.project_registered.pk,
                )
                .exclude(
                    run__run_type__in=[
                        RunMain.RUN_TYPE_STORAGE,
                        RunMain.RUN_TYPE_SCREENING,
                    ],
                    accid="-",
                )
                .distinct("accid")
            )
        else:
            query_set = RawReference.objects.none()

        return self.filter_reference_query_set(query_set, query_string)

    def register_compound_references(self, compound_refs: List[RawReferenceCompound]):
        """
        Register a list of compound references in the database
        """

        for compound_ref in compound_refs:
            self.register_compound_reference(compound_ref)

    def register_compound_reference(self, compound_ref: RawReferenceCompound):
        """
        Register a compound reference in the database
        """

        try:
            compound_ref_model = RawReferenceCompoundModel.objects.get(
                accid=compound_ref.accid, sample__id=compound_ref.sample_id
            )

            compound_ref_model.standard_score = compound_ref.standard_score
            compound_ref_model.global_ranking = compound_ref.global_ranking
            compound_ref_model.ensemble_ranking = compound_ref.ensemble_ranking
            compound_ref_model.manual_insert = compound_ref.manual_insert
            compound_ref_model.mapped_final_report = compound_ref.mapped_final_report
            compound_ref_model.mapped_raw_reference = compound_ref.mapped_raw_reference
            compound_ref_model.selected_mapped_pk = compound_ref.selected_mapped_pk
            compound_ref_model.run_count = compound_ref.run_count
            compound_ref_model.save()

        except RawReferenceCompoundModel.DoesNotExist:
            sample = PIProject_Sample.objects.get(pk=compound_ref.sample_id)

            description_short = compound_ref.description[:200]
            compound_ref_model = RawReferenceCompoundModel(
                taxid=compound_ref.taxid,
                description=description_short,
                accid=compound_ref.accid,
                sample=sample,
                standard_score=compound_ref.standard_score,
                global_ranking=compound_ref.global_ranking,
                ensemble_ranking=compound_ref.ensemble_ranking,
                manual_insert=compound_ref.manual_insert,
                mapped_final_report=compound_ref.mapped_final_report,
                mapped_raw_reference=compound_ref.mapped_raw_reference,
                selected_mapped_pk=compound_ref.selected_mapped_pk,
                run_count=compound_ref.run_count,
            )
            compound_ref_model.save()

            for run in compound_ref.runs:
                compound_ref_model.runs.add(run)

            for ref_pk in compound_ref.family:
                ref = RawReference.objects.get(pk=ref_pk)
                compound_ref_model.family.add(ref)

    def get_classification_runs(self):

        if self.sample_registered is not None:
            classification_runs = RunMain.objects.filter(
                sample=self.sample_registered, run_type=RunMain.RUN_TYPE_PIPELINE
            )

        elif self.project_registered is not None:

            classification_runs = RunMain.objects.filter(
                sample__project=self.project_registered,
                run_type=RunMain.RUN_TYPE_PIPELINE,
            )

        else:
            classification_runs = RunMain.objects.none()

        return classification_runs

    def create_compound(self, raw_references: List[RawReference]):

        raw_reference_compound = [
            RawReferenceCompound(raw_reference) for raw_reference in raw_references
        ]

        classification_runs = self.get_classification_runs()

        # pks of classification runs as integer list
        runs_pks = [run.pk for run in classification_runs]

        # if classification_runs.exists():
        self.sample_reference_tables_filter()

        self.update_scores_compound_references(raw_reference_compound)
        self.register_compound_references(raw_reference_compound)

    def retrieve_compound_references(
        self, query_string: Optional[str] = None
    ) -> QuerySet:
        """
        Retrieve compound references for a sample
        """
        compound_refs = self.query_sample_compound_references(query_string)

        if compound_refs.exists():

            return compound_refs

        compound_refs = self.create_compound_references(query_string=query_string)

        return compound_refs

    def create_compound_references(self, query_string: Optional[str] = None):
        """
        Create compound references for a sample_name, query_string):

        Returns a list of references that match the query string.
        :param query_string:
        :return:

        """
        try:
            references = self.query_sample_references(query_string)
        except Exception as e:
            print(e)

        if references.exists():
            self.create_compound(references)

            compound_refs = self.query_sample_compound_references(query_string)

        else:
            compound_refs = RawReferenceCompoundModel.objects.none()

        return compound_refs

    def reference_table_renamed(self, merged_table, rename_dict: dict):
        proxy_ref = merged_table.copy()

        for key, value in rename_dict.items():
            if key in proxy_ref.columns:
                if value in proxy_ref.columns:
                    # remove the column if it exists
                    proxy_ref = proxy_ref.drop(columns=[value])

                proxy_ref = proxy_ref.rename(columns={key: value})

        # proxy_ref = proxy_ref.rename(columns=rename_dict)

        proxy_ref["taxid"] = proxy_ref["taxid"].astype(int)
        proxy_ref["counts"] = proxy_ref["counts"].astype(float).astype(int)
        proxy_ref = proxy_ref[proxy_ref["counts"] > 0]
        proxy_ref = proxy_ref[proxy_ref["taxid"] > 0]
        proxy_ref = proxy_ref[proxy_ref["description"] != "-"]
        proxy_ref = proxy_ref[proxy_ref["accid"] != "-"]

        return proxy_ref

    def merge_ref_tables(self):
        run_references_tables = self.merge_ref_tables_use_standard_score(
            self.list_tables
        )

        run_references_tables = run_references_tables[run_references_tables.taxid != 0]

        return run_references_tables

    @staticmethod
    def simplify_by_description(df: pd.DataFrame):
        if "description" not in df.columns:
            return df

        df["description_first"] = df["description"].str.split(" ").str[0]

        df = df.sort_values("standard_score", ascending=False)
        df = df.drop_duplicates(subset=["description_first"], keep="first")

        df.drop(columns=["description_first"], inplace=True)

        return df

    # def collect_references_table_all(
    #    self,
    # ) -> pd.DataFrame:
    #    references = self.collect_references_all()
    #
    #    references_table = self.references_table_from_query(references)
    #    # references_table= sample_reference_tables()
    #    return references_table
