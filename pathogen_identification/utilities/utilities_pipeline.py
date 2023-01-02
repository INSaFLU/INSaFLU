import itertools as it
import logging
import os
from collections import defaultdict

import networkx as nx
import numpy as np
import pandas as pd
from django.contrib.auth.models import User
from django.db.models import Q
from pathogen_identification.constants_settings import (
    ConstantsSettings,
    Pipeline_Makeup,
)
from constants.constants import Televir_Metadata_Constants as Televir_Metadata
from pathogen_identification.models import (
    ParameterSet,
    PIProject_Sample,
    Projects,
    SoftwareTree,
    SoftwareTreeNode,
)
from pathogen_identification.utilities.utilities_televir_dbs import Utility_Repository
from settings.constants_settings import ConstantsSettings as CS
from settings.models import Parameter, PipelineStep, Software, Technology
from constants.constants import Televir_Directory_Constants as Televir_Directories

tree = lambda: defaultdict(tree)


def make_tree(lst):
    d = tree()
    for x in lst:
        curr = d
        for item in x:
            curr = curr[item]
    return d


class PipelineTree:

    technology: str
    nodes: list
    edges: dict
    leaves: list
    makeup: int

    def __init__(
        self, technology: str, nodes: list, edges: dict, leaves: list, makeup: int
    ):
        self.technology = technology
        self.nodes = nodes
        self.edges = edges
        self.leaves = leaves
        self.edge_dict = [(x[0], x[1]) for x in self.edges]
        self.makeup = makeup

        self.logger = logging.getLogger(__name__)
        self.logger.info(
            f"PipelineTree: {self.technology}, {self.nodes}, {self.edges}, {self.leaves}, {self.makeup}"
        )

    def get_parents_dict(self):
        """
        Get a dictionary of parents for each node
        """
        self.generate_graph()
        parents_dict = {}

        for node, name in enumerate(self.nodes):
            parents = list(self.graph.predecessors(node))
            if len(parents) == 0:
                parents = None
            else:
                parents = parents[0]
            parents_dict[node] = parents
        return parents_dict

    def node_from_index(self, nix):

        return self.nodes[nix]

    def get_path_explicit(self, path: list) -> list:
        """return nodes names for nodes index list"""

        return [(x, self.node_from_index(x)) for x in path]

    def generate_graph(self):
        """
        Generate a graph of pipeline
        """
        nodes_index = [i for i, x in enumerate(self.nodes)]

        self.graph = nx.DiGraph()
        self.graph.add_edges_from(self.edge_dict)
        self.graph.add_nodes_from(nodes_index)

    def get_all_graph_paths(self) -> dict:
        """
        Get all possible paths in the pipeline
        """

        self.generate_graph()
        all_paths = list(nx.all_simple_paths(self.graph, 0, self.leaves))
        all_paths_explicit = [self.get_path_explicit(path) for path in all_paths]

        path_dict = {
            all_paths_explicit[x][-1][0]: self.df_from_path(path)
            for x, path in enumerate(all_paths)
        }

        return path_dict

    def get_all_graph_paths_explicit(self) -> dict:
        """
        Get all possible paths in the pipeline
        """

        self.generate_graph()
        all_paths = list(nx.all_simple_paths(self.graph, 0, self.leaves))
        all_paths = [self.get_path_explicit(path) for path in all_paths]
        path_dict = {path[-1][0]: path for x, path in enumerate(all_paths)}

        return path_dict

    def df_from_path(self, path: list) -> pd.DataFrame:
        """
        Generate a dataframe from a path
        """
        path = self.get_path_explicit(path)
        path = [x[1] for x in path]
        df = []
        path = [x for x in path if x[0] != "root"]
        current_module = None
        for ix, node in enumerate(path):
            node_type = node[2]
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
                            ]
                        )

                current_module = {
                    "module": node[0],
                    "software": node[1],
                    "params": {
                        f"{node[1].upper()}_ARGS": "",
                    },
                }
                continue

            current_module["params"][node[0]] = node[1]

        if current_module:
            for param, value in current_module["params"].items():
                df.append(
                    [
                        current_module.get("module"),
                        current_module.get("software"),
                        param,
                        value,
                    ]
                )
        df = pd.DataFrame(df, columns=["module", "software", "parameter", "value"])

        return df


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
    with software specific installed databases. Creates a pipeline tree from the combined information."""

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
        self.pipeline_order = pipe_makeup_manager.get_makeup(self.pipeline_makeup)

        self.existing_pipeline_order = [
            x for x in self.pipeline_order if x in pipelines_available
        ]

        self.software_name_list = combined_table.software_name.unique().tolist()

        self.combined_table = combined_table

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

                    new_param = round(float(row.parameter), 2)
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

        old_tree_list = old_tree.nodes
        new_tree_list = new_tree.nodes

        old_tree_set = set(old_tree_list)
        new_tree_set = set(new_tree_list)

        diff = new_tree_set.symmetric_difference(old_tree_set)

        return diff

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
        self.logger.info(f"Checking software db available: {software_name}")

        return self.utility_repository.check_exists(
            "software", "name", software_name.lower()
        )

    def set_software_list(self, software_list):

        self.software_name_list = software_list

    def get_software_list(self):

        self.software_name_list = Software.objects.filter(
            type_of_use=Software.TYPE_OF_USE_televir_global
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

    def get_software_dbs_if_exist(self, software_name: str) -> pd.DataFrame:

        fields = self.utility_repository.select_explicit_statement(
            "software", "name", software_name.lower()
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
                    # f"{software.upper()}_DB": self.software_dbs_dict[software.lower()],
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

    def fill_dict(self, ix, branch_left) -> nx.DiGraph:
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
            return False

        print("matched_path: ", matched_path)

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

            if nodes_index_dict[child_main] in pipe_tree.leaves:
                return nodes_index_dict[child_main]
            try:
                nodes_index_dict[child_main]
            except KeyError:
                self.logger.info(f"{child_main} node not in tree nodes")
                return None

            if child_main not in explicit_edge_dict[parent_main].index:
                self.logger.info(f"Child {child} not in parent {parent}")
                return None

            parent = child
            parent_main = child_main


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
        # print(software_db_dict)
        # print("#####")

        def fix_row(row):

            if not row.parameter:
                return [""]
            if not row.can_change or not row.range_available:
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

                if row.parameter_name == "--db" and software_db_dict:
                    software_name = row.software_name
                    # print(software_name)
                    possibilities = [software_name, software_name.lower()]
                    if "_" in software_name:
                        possibilities.append(software_name.split("_")[0])

                    # print(possibilities)

                    for p in possibilities:
                        if p in software_db_dict:
                            new_range = software_db_dict[p]
                            break

                return new_range

        print(combined_table.shape)

        combined_table["parameter"] = combined_table.apply(fix_row, axis=1)
        combined_table = combined_table.reset_index(drop=True)

        combined_table = combined_table.explode("parameter")

        return combined_table

    def get_user_active_software_tables(self, technology: str, owner: User):
        """
        Get software tables for a user
        """

        software_available = Software.objects.filter(
            type_of_use=Software.TYPE_OF_USE_televir_global,
            technology__name=technology,
            is_to_run=True,
            owner=owner,
        )

        parameters_available = Parameter.objects.filter(
            software__in=software_available,
            televir_project=None,
        )

        software_table = pd.DataFrame(software_available.values())

        parameters_table = pd.DataFrame(parameters_available.values())

        return software_table, parameters_table

    def get_software_tables_global(self, technology: str):
        """
        Get software tables for a user
        """

        software_available = Software.objects.filter(
            type_of_use=Software.TYPE_OF_USE_televir_global,
            technology__name=technology,
        ).distinct()

        parameters_available = Parameter.objects.filter(
            software__in=software_available,
            televir_project=None,
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
            type_of_use=Software.TYPE_OF_USE_televir_project,
            technology__name=project.technology,
            is_to_run=True,
        )

        parameters_available = Parameter.objects.filter(
            software__in=software_available, televir_project=project, is_to_run=True
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
        ).rename(
            columns={
                "id_x": "software_id",
                "id_y": "parameter_id",
                "name_extended": "software_name",
                "name_y": "parameter_name",
                "is_to_run_x": "software_is_to_run",
                "is_to_run_y": "parameter_is_to_run",
            }
        )

        combined_table = combined_table[combined_table.type_of_use.isin([5, 6])]

        combined_table["pipeline_step"] = combined_table["pipeline_step_id"].apply(
            lambda x: PipelineStep.objects.get(id=int(x)).name
        )

        combined_table["technology"] = combined_table["technology_id"].apply(
            lambda x: Technology.objects.get(id=int(x)).name
        )

        combined_table = combined_table.reset_index(drop=True)

        combined_table = combined_table.loc[
            :, ~combined_table.T.duplicated(keep="last")
        ]

        return combined_table

    def generate_combined_parameters_table(self, technology: str):
        """"""
        software_table, parameters_table = self.get_software_tables_global(technology)

        return self.merge_software_tables(software_table, parameters_table)

    def generate_combined_parameters_table_project(
        self, owner: User, project: Projects
    ):
        """"""
        software_table, parameters_table = self.get_software_tables_project(
            owner, project
        )

        if parameters_table.shape[0] == 0:
            self.logger.info("No parameters for this project, using global")
            software_table, parameters_table = self.get_user_active_software_tables(
                project.technology, owner
            )

        merged_table = self.merge_software_tables(software_table, parameters_table)

        return merged_table

    def check_default_software_tree_exists(
        self, technology: Technology, global_index: int
    ):

        try:
            software_tree = (
                SoftwareTree.objects.filter(
                    global_index=global_index, technology=technology
                )
                .order_by("date_created")
                .last()
            )

            if software_tree:
                return True

        except SoftwareTree.DoesNotExist:
            return None

    def get_software_tree_index(self, technology: Technology, global_index: int):

        if self.check_default_software_tree_exists(technology, global_index):
            software_tree = (
                SoftwareTree.objects.filter(
                    global_index=global_index, technology=technology
                )
                .order_by("date_created")
                .last()
            )

            return software_tree.pk
        else:
            return None

    def check_ParameterSet_exists(
        self, sample: PIProject_Sample, leaf: SoftwareTreeNode, project: Projects
    ):
        self.logger.info("Checking if ParameterSet exists")
        try:
            parameter_set = ParameterSet.objects.get(
                sample=sample, leaf=leaf, project=project
            )
            return True

        except ParameterSet.DoesNotExist:
            return False

    def check_ParameterSet_available(
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

    def check_ParameterSet_processed(
        self, sample: PIProject_Sample, leaf: SoftwareTreeNode, project: Projects
    ):

        if not self.check_ParameterSet_exists(sample, leaf, project):
            return False

        parameter_set = ParameterSet.objects.get(
            sample=sample, leaf=leaf, project=project
        )

        if parameter_set.status == ParameterSet.STATUS_FINISHED:
            return True
        else:
            return False

    def get_software_tree_node_index(
        self, owner: User, technology: str, global_index: int, node_index: int
    ):

        software_tree_index = self.get_software_tree_index(technology, global_index)

        if software_tree_index:
            try:
                software_tree_node = SoftwareTreeNode.objects.get(
                    software_tree=software_tree_index, index=node_index
                )
                return software_tree_node.pk

            except SoftwareTreeNode.DoesNotExist:
                return None
        else:
            return None

    def query_software_default_tree(
        self, technology: str, global_index: int
    ) -> PipelineTree:
        """
        Generate a default software tree for a user
        """

        software_tree = (
            SoftwareTree.objects.filter(
                global_index=global_index, technology=technology
            )
            .order_by("date_created")
            .last()
        )

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
            technology=technology,
            nodes=[x[1] for x in sorted(nodes)],
            edges=edges,
            leaves=leaves,
            makeup=global_index,
        )

    def update_SoftwareTree_nodes(
        self, software_tree: SoftwareTree, tree: PipelineTree
    ):
        """
        Update the nodes of a software tree
        """
        parent_dict = tree.get_parents_dict()

        for index, node in enumerate(tree.nodes):
            is_leaf = int(index in tree.leaves)

            try:
                tree_node = SoftwareTreeNode.objects.get(
                    software_tree=software_tree, index=index
                )

            except SoftwareTreeNode.DoesNotExist:

                parent_node = parent_dict.get(index, None)

                if parent_node != None:

                    parent_node = SoftwareTreeNode.objects.get(
                        software_tree=software_tree, index=parent_dict[index]
                    )

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

    def update_software_tree(self, tree: PipelineTree):
        """
        Update SoftwareTree table
        """
        global_index = tree.makeup

        software_tree = (
            SoftwareTree.objects.filter(
                global_index=global_index, technology=tree.technology
            )
            .order_by("date_created")
            .last()
        )

        if software_tree:
            new_version = software_tree.version + 1

            self.logger.info("Creating new software tree")
            software_tree = SoftwareTree(
                global_index=global_index,
                technology=tree.technology,
                version=new_version,
            )

            software_tree.save()

            self.update_SoftwareTree_nodes(software_tree, tree)

        else:
            self.logger.info("Creating new software tree")
            software_tree = SoftwareTree(
                global_index=global_index,
                technology=tree.technology,
                version=0,
            )

            software_tree.save()

            self.update_SoftwareTree_nodes(software_tree, tree)


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
        self.logger.setLevel(logging.INFO)
        if self.logger.hasHandlers():
            self.logger.handlers.clear()

        self.logger.info("Utils_Manager initialized")

    def get_leaf_parameters(self, parameter_leaf: SoftwareTreeNode) -> pd.DataFrame:
        """ """
        pipeline_tree = self.generate_software_tree(
            parameter_leaf.software_tree.technology,
            parameter_leaf.software_tree.global_index,
        )

        if parameter_leaf.index not in pipeline_tree.leaves:
            raise Exception("Node is not a leaf")

        all_paths = pipeline_tree.get_all_graph_paths()

        return all_paths[parameter_leaf.index]

    def check_runs_to_deploy(self, user: User, project: Projects):
        """
        Check if there are runs to run
        """

        technology = project.technology
        samples = PIProject_Sample.objects.filter(project=project)
        local_tree = self.generate_project_tree(technology, project, user)

        self.logger.info("Checking runs to deploy")
        tree_makeup = local_tree.makeup

        pipeline_tree = self.generate_software_tree(technology, tree_makeup)
        self.logger.info("Pipeline tree generated")
        pipeline_tree_index = self.get_software_tree_index(technology, tree_makeup)
        self.logger.info("Pipeline tree index generated")
        local_paths = local_tree.get_all_graph_paths_explicit()
        sample = samples[0]

        runs_to_deploy = 0
        self.logger.info(
            "now going to start checking if existing paths correspond to branches in trees"
        )
        for sample in samples:

            for leaf, path in local_paths.items():

                try:
                    matched_path = self.utility_manager.match_path_to_tree(
                        path, pipeline_tree
                    )
                except Exception as e:
                    self.logger.info("Path not matched to tree")
                    self.logger.info(e)
                    continue

                self.logger.info("Matched path to tree")

                matched_path_node = SoftwareTreeNode.objects.get(
                    software_tree__pk=pipeline_tree_index, index=matched_path
                )

                exists = self.parameter_util.check_ParameterSet_exists(
                    sample=sample, leaf=matched_path_node, project=project
                )

                if exists:

                    if self.parameter_util.check_ParameterSet_processed(
                        sample=sample, leaf=leaf, project=project
                    ):
                        self.logger.info("parameter set processed")

                        continue

                else:
                    runs_to_deploy += 1
        if runs_to_deploy > 0:
            return True

        return False

    def get_all_technology_pipelines(self, technology: str, tree_makeup: int) -> dict:
        """
        Get all pipelines for a technology
        """

        pipeline_tree = self.generate_software_tree(technology, tree_makeup)

        all_paths = pipeline_tree.get_all_graph_paths()
        return all_paths

    def get_software_tree_index(self, technology: Technology, tree_makeup: int):
        """
        Get the software tree index from model
        """
        return self.parameter_util.get_software_tree_index(technology, tree_makeup)

    def generate_software_tree(self, technology, tree_makeup: int):
        """
        Generate a software tree for a technology and a tree makeup
        """

        if self.parameter_util.check_default_software_tree_exists(
            technology, global_index=tree_makeup
        ):
            return self.parameter_util.query_software_default_tree(
                technology, global_index=tree_makeup
            )

        else:
            raise Exception("No software tree for technology")

    def generate_software_base_tree(self, technology, tree_makeup: int):
        """
        Generate a software tree for a technology and a tree makeup
        """

        pipeline_setup = Pipeline_Makeup()
        makeup_steps = pipeline_setup.get_makeup(tree_makeup)

        combined_table = self.parameter_util.generate_combined_parameters_table(
            technology
        )

        combined_table = combined_table[combined_table.pipeline_step.isin(makeup_steps)]

        if len(combined_table) == 0 or "can_change" not in combined_table.columns:
            return PipelineTree(
                technology=technology,
                nodes=[],
                edges={},
                leaves=[],
                makeup=tree_makeup,
            )

        self.utility_manager.get_software_db_dict()

        full_table = self.parameter_util.expand_parameters_table(
            combined_table, software_db_dict=self.utility_manager.software_dbs_dict
        )

        self.utility_manager.input(full_table, technology=technology)

        pipeline_tree = self.utility_manager.generate_default_software_tree()

        if self.parameter_util.check_default_software_tree_exists(
            technology, global_index=tree_makeup
        ):
            existing_pipeline_tree = self.parameter_util.query_software_default_tree(
                technology, global_index=tree_makeup
            )

            tree_differences = self.utility_manager.compare_software_trees(
                existing_pipeline_tree
            )

            if len(tree_differences) > 0:
                self.parameter_util.update_software_tree(pipeline_tree)
        else:

            self.parameter_util.update_software_tree(pipeline_tree)

        return pipeline_tree

    def generate_project_tree(self, technology, project: Projects, owner: User):
        """
        Generate a project tree
        """

        combined_table = self.parameter_util.generate_combined_parameters_table_project(
            owner, project
        )

        utility_drone = Utility_Pipeline_Manager()
        utility_drone.input(combined_table, technology=technology)

        self.logger.info("Generating project tree")

        pipeline_tree = utility_drone.generate_default_software_tree()

        return pipeline_tree

    def generate_default_trees(self):
        """
        Generate default trees for all technologies and makeups
        """
        technology_trees = {}
        pipeline_makeup = Pipeline_Makeup()
        for technology in self.utility_technologies:
            for makeup in pipeline_makeup.get_makeup_list():

                technology_trees[technology] = self.generate_software_base_tree(
                    technology, makeup
                )
