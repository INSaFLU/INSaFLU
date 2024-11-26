import copy
import logging
import os
import re
import traceback  # for debugging
from typing import Callable, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from django.db.models import QuerySet

from constants.constants import Televir_Metadata_Constants as Televir_Metadata
from fluwebvirus.settings import STATIC_ROOT
from pathogen_identification.constants_settings import ConstantsSettings as PIConstants
from pathogen_identification.deployment_main import PathogenIdentificationDeploymentCore
from pathogen_identification.models import (
    FinalReport,
    ParameterSet,
    PIProject_Sample,
    Projects,
    RunMain,
    SoftwareTree,
    SoftwareTreeNode,
)
from pathogen_identification.modules.object_classes import Remap_Target
from pathogen_identification.modules.remap_class import Mapping_Instance
from pathogen_identification.modules.run_main import RunMainTree_class
from pathogen_identification.utilities.televir_parameters import TelevirParameters
from pathogen_identification.utilities.update_DBs_tree import (
    Update_Assembly,
    Update_Classification,
    Update_Remap,
    Update_RunMain_Initial,
    Update_RunMain_Secondary,
)
from pathogen_identification.utilities.utilities_pipeline import (
    Pipeline_Makeup,
    PipelineTree,
    Utils_Manager,
)
from pathogen_identification.utilities.utilities_views import ReportSorter
from settings.constants_settings import ConstantsSettings
from utils.utils import Utils


def logger_copy(x, memo):
    return x


copy._deepcopy_dispatch[logging.Logger] = logger_copy


class PathogenIdentification_TreeDeployment(PathogenIdentificationDeploymentCore):
    project: Projects
    sample: PIProject_Sample
    prefix: str
    rdir: str
    threads: int
    run_engine: RunMainTree_class
    params = dict
    run_params_db = pd.DataFrame()
    pk: int
    username: str
    prepped: bool
    sent: bool
    parameter_set: ParameterSet

    STATUS_ZERO = 0
    STATUS_PREPPED = 1
    STATUS_RUNNING = 2
    STATUS_SENT = 3

    assembly_udated: bool
    classification_updated: bool

    def __init__(
        self,
        sample: PIProject_Sample,
        deployment_root_dir: str = "/tmp/insaflu/insaflu_something",
        dir_branch: str = "deployment",
        prefix: str = "run",
        threads: int = 3,
    ) -> None:

        super().__init__(sample, deployment_root_dir, dir_branch, prefix, threads)

        self.sent = False
        self.assembly_udated = False
        self.classification_updated = False


class Tree_Node:
    module: str
    name: str
    node_index: int
    children: list
    parameters: pd.DataFrame
    software_tree_pk: int
    run_manager: PathogenIdentification_TreeDeployment
    tree_node: SoftwareTreeNode
    parameter_set: ParameterSet

    updated_classification: bool
    updated_assembly: bool

    def __init__(
        self,
        pipe_tree: PipelineTree,
        node_index: int,
        software_tree_pk: int,
        sample: PIProject_Sample,
    ):
        node_metadata = pipe_tree.node_index.loc[node_index].node

        self.module = node_metadata[0]
        self.node_index = node_index

        self.branch = pipe_tree.nodes_df.loc[node_index].branch
        self.children = pipe_tree.edge_df[
            pipe_tree.edge_df.parent == node_index
        ].child.tolist()

        self.software_tree_pk = software_tree_pk

        self.parameters = self.determine_params(pipe_tree, sample)
        self.leaves = pipe_tree.leaves_from_node_using_graph(node_index)

    def run_reference_overlap_analysis(self):
        # run = RunMain.objects.filter(parameter_set=self.parameter_set).first()
        final_report = FinalReport.objects.filter(
            sample=self.parameter_set.sample,  # run=run
        ).order_by("-coverage")
        #
        report_layout_params = TelevirParameters.get_report_layout_params(
            project_pk=self.parameter_set.project.pk
        )
        report_sorter = ReportSorter(
            self.parameter_set.sample, final_report, report_layout_params
        )
        report_sorter.sort_reports_save()

    def receive_run_manager(self, run_manager: PathogenIdentification_TreeDeployment):
        run_manager.prefix = f"run_leaf_{self.node_index}"
        run_manager.configure_constants()
        run_manager.import_params(self.parameters)

        run_manager.run_main_prep_check_first()

        self.run_manager = run_manager

    def _is_node_leaf(self):
        return len(self.children) == 0

    def generate_software_tree_node_entry(self, pipe_tree: PipelineTree):
        if not self._is_node_leaf():
            return

        node_metadata = pipe_tree.node_index.loc[self.node_index].node
        software_tree = SoftwareTree.objects.get(pk=self.software_tree_pk)

        try:
            tree_node = SoftwareTreeNode.objects.get(
                software_tree=software_tree,
                index=self.node_index,
                name=node_metadata[0],
                value=node_metadata[1],
                node_type=node_metadata[2],
            )
        except SoftwareTreeNode.DoesNotExist:
            tree_node = None

        return tree_node

    def setup_parameterset(
        self, project: Projects, sample: PIProject_Sample, node: SoftwareTreeNode
    ):
        utils_manager = Utils_Manager()

        try:
            parameter_set = ParameterSet.objects.get(
                project=project, sample=sample, leaf=node
            )

        except ParameterSet.DoesNotExist:
            parameter_set = ParameterSet()
            parameter_set.project = project
            parameter_set.sample = sample
            parameter_set.status = ParameterSet.STATUS_RUNNING
            parameter_set.leaf = node
            parameter_set.save()

        except ParameterSet.MultipleObjectsReturned:
            parameter_set = ParameterSet.objects.filter(
                project=project, sample=sample, leaf=node
            ).first()

        return parameter_set

    def register(self, project: Projects, sample: PIProject_Sample, tree: PipelineTree):
        tree_node = self.generate_software_tree_node_entry(tree)
        if tree_node is None:
            return False

        parameter_set = self.setup_parameterset(project, sample, tree_node)

        self.parameter_set = parameter_set

    def register_finished(
        self, project: Projects, sample: PIProject_Sample, tree: PipelineTree
    ):
        self.register(project, sample, tree)

        if self.parameter_set is None:
            return False

        self.parameter_set.status = ParameterSet.STATUS_FINISHED
        self.parameter_set.save()

        return True

    def register_failed(
        self, project: Projects, sample: PIProject_Sample, tree: PipelineTree
    ):
        self.register(project, sample, tree)

        if self.parameter_set is None:
            return False

        self.parameter_set.status = ParameterSet.STATUS_ERROR
        self.parameter_set.save()

        self.run_manager.delete_run()

        return True

    def register_running(
        self, project: Projects, sample: PIProject_Sample, tree: PipelineTree
    ):
        tree_node = self.generate_software_tree_node_entry(tree)

        if tree_node is None:
            return False

        parameter_set = self.setup_parameterset(project, sample, tree_node)
        if parameter_set is None:
            return False

        parameter_set.status = ParameterSet.STATUS_RUNNING
        parameter_set.save()

        self.tree_node = tree_node
        self.parameter_set = parameter_set

        return True

    def determine_params(self, pipe_tree: PipelineTree, sample: PIProject_Sample):
        arguments_list = []

        for node in self.branch:
            node_metadata = pipe_tree.node_index.loc[node].node

            path_to_node = pipe_tree.paths_to_node(node)
            ps_visited = pipe_tree.check_if_leaf_steps_exist_list(
                path_to_node, sample=sample
            )

            node_metadata = (
                node_metadata[0],
                node_metadata[1],
                node_metadata[2],
                ps_visited,
            )
            arguments_list.append(node_metadata)

            # ps_track.append(ps_visited)

        arguments_df = pd.DataFrame(
            arguments_list, columns=["parameter", "value", "flag", "leaves"]
        )

        arguments_df = arguments_df[arguments_df.parameter != "root"]

        if len(arguments_df) == 0:
            return pd.DataFrame()

        module_df = arguments_df[arguments_df.flag == "module"]
        module = module_df.parameter.values[0]
        software = module_df.value.values[0]
        parameters_df = arguments_df[arguments_df.flag == "param"]

        parameters_df["software"] = software
        parameters_df["module"] = module

        return parameters_df


class TreeNode_Iterator:
    def __init__(self):
        self.node_list = []

    def add_node(self, node: Tree_Node):
        self.node_list.append(node)

    def __getitem__(self, index) -> Tree_Node:
        return self.node_list[index]

    def __len__(self):
        return len(self.node_list)


from abc import ABC, abstractmethod
from functools import wraps


def check_planned(method):
    @wraps(method)
    def _impl(self, *method_args, **method_kwargs):
        if method_args[0].run_manager.run_engine.remap_prepped:
            return False
        method_output = method(self, *method_args, **method_kwargs)
        return method_output

    return _impl


class ClassificationMonitor(ABC):
    @abstractmethod
    def ready_to_merge(node: Tree_Node):
        pass

    @abstractmethod
    def classification_performed(node: Tree_Node):
        pass


class ClassificationMonitor_ContigOnly(ClassificationMonitor):
    @check_planned
    def ready_to_merge(self, node: Tree_Node):
        if node.run_manager.run_engine.contig_classification_performed:
            return True

        return False

    def classification_performed(self, node: Tree_Node):
        if node.run_manager.run_engine.contig_classification_performed:
            return True

        return False

    def __str__(self) -> str:
        return "ContigOnly"


class ClassificationMonitor_ContigAndReads(ClassificationMonitor):
    @check_planned
    def ready_to_merge(self, node: Tree_Node):
        if (
            node.run_manager.run_engine.contig_classification_performed
            and node.run_manager.run_engine.read_classification_performed
        ):
            return True

        return False

    def classification_performed(self, node: Tree_Node):
        if (
            node.run_manager.run_engine.contig_classification_performed
            and node.run_manager.run_engine.read_classification_performed
        ):
            return True

        return False

    def __str__(self) -> str:
        return "ContigAndReads"


class ClassificationMonitor_ReadsOnly(ClassificationMonitor):
    @check_planned
    def ready_to_merge(self, node: Tree_Node):
        if node.run_manager.run_engine.read_classification_performed:
            return True

        return False

    def classification_performed(self, node: Tree_Node):
        if node.run_manager.run_engine.read_classification_performed:
            return True

        return False

    def __str__(self) -> str:
        return "ReadsOnly"


class ClassificationMonitor_Factory:
    def __init__(self):
        self.tree_makeup_config = Pipeline_Makeup()

    def get_monitor(self, tree: PipelineTree):
        modules = self.tree_makeup_config.get_makeup(tree.makeup)

        if (
            ConstantsSettings.PIPELINE_NAME_contig_classification in modules
            and ConstantsSettings.PIPELINE_NAME_read_classification in modules
        ):
            return ClassificationMonitor_ContigAndReads()

        elif ConstantsSettings.PIPELINE_NAME_contig_classification in modules:
            return ClassificationMonitor_ContigOnly()

        elif ConstantsSettings.PIPELINE_NAME_read_classification in modules:
            return ClassificationMonitor_ReadsOnly()

        else:
            return None


class Tree_Progress:
    tree: PipelineTree
    current_nodes: List[Tree_Node]
    current_module: str

    sample: PIProject_Sample
    project: Projects
    updated_classification: bool

    def __init__(
        self,
        pipe_tree: PipelineTree,
        sample: PIProject_Sample,
        project: Projects,
    ):
        pipe_tree.nodes_df = pd.DataFrame(
            pipe_tree.nodes_compress, columns=["node", "branch"]
        ).set_index("node")

        pipe_tree.edge_df = pd.DataFrame(
            pipe_tree.edge_compress, columns=["parent", "child"]
        )

        self.logger = logging.getLogger(__name__)
        ## set logger level
        self.logger.setLevel(logging.INFO)

        self.tree = pipe_tree
        self.sample = sample
        self.project = project
        self.classification_monitor = ClassificationMonitor_Factory().get_monitor(
            pipe_tree
        )
        self.updated_classification = False

        self.initialize_nodes()
        self.determine_current_module_from_nodes()

    def setup_deployment_manager(self):
        utils = Utils()
        temp_dir = utils.get_temp_dir()

        prefix = f"{self.sample.sample.pk}_{self.sample.sample.name}"

        deployment_directory_structure = os.path.join(
            PIConstants.televir_subdirectory,
            f"{self.project.owner.pk}",
            f"{self.project.pk}",
            f"{self.sample.sample.pk}",
            prefix,
        )

        deployment_manager = PathogenIdentification_TreeDeployment(
            self.sample,
            # self.project,
            # self.project.owner.username,
            # self.project.technology,
            deployment_root_dir=temp_dir,
            dir_branch=deployment_directory_structure,
            threads=PIConstants.DEPLOYMENT_THREADS,
        )

        return deployment_manager

    def register_node_leaves(self, node: Tree_Node):
        if node.node_index in node.leaves:
            self.register_finished(node)

        for leaf in node.leaves:
            leaf_node = self.spawn_node_child(node, leaf)
            _ = self.register_node_safe(leaf_node)

    def update_node_leaves_dbs(self, node: Tree_Node):
        for leaf in node.leaves:
            leaf_node = self.spawn_node_child_prepped(node, leaf)
            self.register_node(leaf_node)
            self.update_node_dbs(leaf_node)

            node.run_manager.classification_updated = (
                leaf_node.run_manager.classification_updated
            )
            node.run_manager.assembly_udated = leaf_node.run_manager.assembly_udated

    def register_finished(self, node: Tree_Node):
        self.logger.info(f"Registering node {node.node_index} as finished")

        registration_success = node.register_finished(
            self.project, self.sample, self.tree
        )

        if registration_success:
            self.logger.info(f"Node {node.node_index} registered as finished")
        else:
            self.logger.info(f"Node {node.node_index} failed to register as finished")

        return registration_success

    def register_failed_children(self, node: Tree_Node):
        if node.node_index in node.leaves:
            # self.submit_node_run(node)
            registration_success = node.register_failed(
                self.project, self.sample, self.tree
            )

        for leaf in node.leaves:
            leaf_node = self.spawn_node_child(node, leaf)
            # self.submit_node_run(leaf_node)
            registration_success = leaf_node.register_failed(
                self.project, self.sample, self.tree
            )

    def initialize_nodes(self):
        origin_node = Tree_Node(
            self.tree,
            0,
            software_tree_pk=self.tree.software_tree_pk,
            sample=self.sample,
        )

        run_manager = self.setup_deployment_manager()

        origin_node.receive_run_manager(run_manager)

        self.register_node_leaves(origin_node)

        self.current_nodes = [origin_node]

    def get_current_module(self):
        return self.current_module

    def update_current_module(self, new_nodes: List[Tree_Node]):
        node_indices = [node.node_index for node in new_nodes]

        if len(new_nodes) == 0:
            self.current_module = "end"
        elif set(node_indices) == set(self.tree.leaves):
            self.current_module = "end"
        else:
            self.current_nodes = new_nodes
            self.determine_current_module_from_nodes()

    def determine_current_module_from_nodes(self):
        self.current_module = self.current_nodes[0].module

    def transfer_sample(self, node: Tree_Node, child: Tree_Node) -> Tree_Node:

        if node.node_index == 0:
            return child
        if "module" in node.run_manager.run_engine.config.keys():
            if "module" in child.run_manager.run_engine.config.keys():
                run_sample_copy = copy.deepcopy((node.run_manager.run_engine.sample))
                child.run_manager.run_engine.sample = run_sample_copy

        return child

    def spawn_node_child(self, node: Tree_Node, child: int) -> Tree_Node:
        new_node = Tree_Node(
            self.tree, child, node.software_tree_pk, sample=self.sample
        )

        run_manager_copy = copy.deepcopy(node.run_manager)
        run_manager_copy.classification_updated = (
            node.run_manager.classification_updated
        )
        run_manager_copy.assembly_udated = node.run_manager.assembly_udated
        new_node.receive_run_manager(run_manager_copy)

        return new_node

    def spawn_node_child_prepped(self, node: Tree_Node, child: int) -> Tree_Node:
        new_node = self.spawn_node_child(node, child)

        new_node.run_manager.update_config_prefix()
        new_node = self.transfer_sample(node, new_node)
        new_node.run_manager.update_engine()

        return new_node

    def register_node(self, node: Tree_Node):
        if node.run_manager.sent:
            return False

        registraction_success = node.register_running(
            self.project, self.sample, self.tree
        )

        return registraction_success

    def update_node_dbs(self, node: Tree_Node, step="initial"):
        try:
            db_updated = Update_RunMain_Initial(
                node.run_manager.run_engine, node.parameter_set
            )
            if not db_updated:
                return False

            if (
                node.run_manager.run_engine.qc_performed
                or node.run_manager.run_engine.enrichment_performed
                or node.run_manager.run_engine.depletion_performed
            ):
                node.run_manager.run_engine.export_sequences()
                db_updated = Update_RunMain_Secondary(
                    node.run_manager.run_engine, node.parameter_set
                )
                if not db_updated:
                    return False

            if (
                node.run_manager.run_engine.assembly_performed
                and node.run_manager.assembly_udated == False
            ):
                node.run_manager.run_engine.export_assembly()
                db_updated = Update_Assembly(
                    node.run_manager.run_engine, node.parameter_set
                )
                if not db_updated:
                    return False

                node.run_manager.assembly_udated = True

            if (
                self.classification_monitor.classification_performed(node)
                and node.run_manager.classification_updated == False
            ):
                node.run_manager.run_engine.plan_remap_prep_safe()
                node.run_manager.run_engine.export_intermediate_reports()
                node.run_manager.run_engine.generate_output_data_classes()
                db_updated = Update_Classification(
                    node.run_manager.run_engine, node.parameter_set
                )
                if not db_updated:
                    return False

                node.run_manager.classification_updated = True

            if node.run_manager.run_engine.remapping_performed:
                node.run_manager.run_engine.export_final_reports()
                node.run_manager.run_engine.Summarize()
                node.run_manager.run_engine.generate_output_data_classes()
                node.run_manager.run_engine.export_logdir()
                db_updated = Update_Remap(
                    node.run_manager.run_engine, node.parameter_set
                )
                if not db_updated:
                    return False

            return True
        except Exception as e:
            self.logger.error("Error updating node dbs, returning false.")
            self.logger.error(e)
            traceback.print_exc()
            return False

    def register_node_safe(self, node: Tree_Node):
        registration_success = self.register_node(node)

        return registration_success

    def merge_node_targets(self):
        """
        Merge targets from all nodes in the current node list

        :return: merged targets"""
        node_merged_targets = [
            n.run_manager.run_engine.merged_targets for n in self.current_nodes
        ]
        node_merged_targets = pd.concat(node_merged_targets, axis=0)
        node_merged_targets = node_merged_targets.drop_duplicates(subset=["taxid"])

        return node_merged_targets

    @staticmethod
    def get_remap_plans(nodes: List[Tree_Node]):
        """
        Get remap plans for all nodes in a list"""

        for n in nodes:
            # if n.run_manager.run_engine.remap_prepped is False:
            n.run_manager.run_engine.metadata_tool.reset()
            n.run_manager.run_engine.plan_remap_prep_safe()

        return nodes

    def get_node_node_targets(self, nodes_list: List[Tree_Node]) -> List[Remap_Target]:
        """
        Get merged targets from a list of nodes

        :param nodes_list: list of nodes"""
        accids_in_list = {}

        combined_list = []

        for node in nodes_list:
            for target in node.run_manager.run_engine.metadata_tool.remap_targets:
                if target.accid not in accids_in_list:
                    accids_in_list[target.accid] = True
                    combined_list.append(target)

        return combined_list

    @staticmethod
    def process_mapping_managerdf(df: pd.DataFrame):
        """
        Process mapping manager dataframe
        """
        if "accid" in df.columns:
            df = df.drop_duplicates(subset=["accid"])
        if "protein_id" in df.columns:
            df = df.drop_duplicates(subset=["protein_id"])
        if "accession_id" in df.columns:
            df = df.drop_duplicates(subset=["accession_id"])

        return df

    def group_nodes_by_sample_sources(self):
        nodes_by_sample_sources = {}
        for node in self.current_nodes:
            sample_source = node.run_manager.run_engine.sample.sources_list()
            if sample_source not in nodes_by_sample_sources.keys():
                nodes_by_sample_sources[sample_source] = []
            nodes_by_sample_sources[sample_source].append(node)
        return nodes_by_sample_sources

    def group_nodes_by_source_and_parameters(self) -> List[List[Tree_Node]]:
        source_paramaters_combinations = {}

        def check_node_in_combination_dict(node: Tree_Node):
            sample_source = node.run_manager.run_engine.sample.sources_list()
            parameters = node.parameters
            for source, register in source_paramaters_combinations.items():
                if set(sample_source) == set(register["source"]) and parameters.equals(
                    register["parameters"]
                ):
                    return source

            return None

        def combination_dict_new_entry(node, sample_source, parameters):
            new_entry_idx = len(source_paramaters_combinations)
            source_paramaters_combinations[new_entry_idx] = {
                "source": sample_source,
                "parameters": parameters,
                "nodes": [node],
            }

        def update_combination_dict(node: Tree_Node):
            sample_source = node.run_manager.run_engine.sample.sources_list()
            parameters = node.parameters
            source = check_node_in_combination_dict(node)

            if source is None:
                combination_dict_new_entry(node, sample_source, parameters)
            else:
                source_paramaters_combinations[source]["nodes"].append(node)

        for node in self.current_nodes:

            update_combination_dict(node)

        grouped_nodes = [
            register["nodes"]
            for source, register in source_paramaters_combinations.items()
        ]

        return grouped_nodes

    def group_nodes_by_module(self):
        nodes_by_module = {}
        for node in self.current_nodes:
            if node.module not in nodes_by_module:
                nodes_by_module[node.module] = []
            nodes_by_module[node.module].append(node)
        return nodes_by_module

    def process_subject(self, volonteer: Tree_Node, group_targets: List[Remap_Target]):
        original_targets = [
            copy.deepcopy(x)
            for x in volonteer.run_manager.run_engine.metadata_tool.remap_targets
        ]

        volonteer.run_manager.update_merged_targets(group_targets)

        run_success = self.run_node(volonteer)

        volonteer.run_manager.update_merged_targets(original_targets)

        mapped_instances_shared = (
            volonteer.run_manager.run_engine.remap_manager.mapped_instances
        )

        return mapped_instances_shared, run_success

    @staticmethod
    def transfer_read_classification_drone(node: Tree_Node, volonteer: Tree_Node):
        """
        Transfer read classification drone from volonteer to node"""
        node.run_manager.run_engine.read_classification_drone = (
            volonteer.run_manager.run_engine.read_classification_drone
        )
        node.run_manager.run_engine.read_classification_performed = True

        return node

    @staticmethod
    def transfer_contig_classification_drone(node: Tree_Node, volonteer: Tree_Node):
        """
        Transfer contig classification drone from volonteer to node"""
        node.run_manager.run_engine.contig_classification_drone = (
            volonteer.run_manager.run_engine.contig_classification_drone
        )
        node.run_manager.run_engine.contig_classification_performed = True

        return node

    def stacked_deployment_read_classification(
        self, nodes_by_sample_sources: List[List[Tree_Node]]
    ):
        self.loop_nodes_deploy(
            nodes_by_sample_sources, self.transfer_read_classification_drone
        )

    def stacked_deployment_contig_classification(
        self, nodes_by_sample_sources: List[List[Tree_Node]]
    ):
        self.loop_nodes_deploy(
            nodes_by_sample_sources, self.transfer_contig_classification_drone
        )

    def loop_nodes_deploy(
        self,
        nodes_by_sample_sources: List[List[Tree_Node]],
        transfer_function: Callable,
    ):
        new_nodes = []

        for nodes_subset in nodes_by_sample_sources:
            if len(nodes_subset) == 0:
                continue

            volonteer = nodes_subset[0]
            run_success = self.run_node(volonteer)

            if run_success:
                for node in nodes_subset:
                    node = transfer_function(node, volonteer)

                    new_nodes.append(node)
            else:
                for node in nodes_subset:
                    self.register_failed_children(node)

        if len(new_nodes) > 0:
            self.current_nodes = new_nodes

    def stacked_deployement_mapping(
        self,
        nodes_by_sample_sources: List[List[Tree_Node]],
    ):
        """
        deploy nodes remap by sample sources."""
        current_nodes = []

        for nodes in nodes_by_sample_sources:

            if len(nodes) == 0:
                continue

            nodes = self.get_remap_plans(nodes)

            group_targets = self.get_node_node_targets(nodes)

            self.logger.info(
                f"#### Total targets registered for remap: {len(group_targets)}"
            )

            volonteer = nodes[0]

            mapped_instances_shared, deployment_success = self.process_subject(
                volonteer, group_targets
            )

            nodes = self.update_mapped_instances(nodes, mapped_instances_shared)
            self.export_intermediate_reports_leaves(nodes)

            if deployment_success:
                current_nodes.extend(nodes)

            else:
                for node in nodes:
                    self.register_failed_children(node)

        if len(current_nodes) > 0:
            self.current_nodes = current_nodes

    @staticmethod
    def update_mapped_instances(
        nodes_to_update: List[Tree_Node],
        mapped_instances_shared: List[Mapping_Instance],
    ):
        new_nodes = []

        for node in nodes_to_update:
            node.run_manager.run_engine.update_mapped_instances(mapped_instances_shared)
            new_nodes.append(node)

        return new_nodes

    def export_intermediate_reports_leaves(self, nodes_to_update: List[Tree_Node]):
        for node in nodes_to_update:
            for leaf in node.leaves:
                leaf_node = self.spawn_node_child(node, leaf)
                leaf_node.run_manager.run_engine.export_intermediate_reports()

    def run_simplified_mapping(self):
        nodes_by_sample_sources = self.group_nodes_by_source_and_parameters()

        self.stacked_deployement_mapping(nodes_by_sample_sources)

    def run_simplified_classification_reads(self):
        nodes_by_sample_sources = self.group_nodes_by_source_and_parameters()

        self.stacked_deployment_read_classification(nodes_by_sample_sources)

    def run_simplified_classification_contigs(self):
        nodes_by_sample_sources = self.group_nodes_by_source_and_parameters()

        self.stacked_deployment_contig_classification(nodes_by_sample_sources)

    def update_tree_nodes(self):
        new_nodes = []
        for node in self.current_nodes:
            children = node.children
            for child in children:
                new_node = self.spawn_node_child_prepped(node, child)
                new_nodes.append(new_node)

        self.update_current_module(new_nodes)

    def run_node(self, node: Tree_Node):
        try:
            node.run_manager.run_main()
            return True
        except Exception as e:
            print("error")
            print(e)
            traceback.print_exc()

            return False

    def run_current_nodes(self):
        new_nodes = []

        for node in self.current_nodes:
            run_success = self.run_node(node)

            if run_success:
                new_nodes.append(node)

            else:
                self.register_node_leaves(node)
                self.register_failed_children(node)

        self.current_nodes = new_nodes

    def run_nodes_sequential(self):
        self.run_current_nodes()

    def run_nodes_remap(self):
        self.run_simplified_mapping()

    def register_leaves_finished(self):
        for node in self.current_nodes:
            for leaf in node.leaves:
                leaf_node = self.spawn_node_child(node, leaf)
                _ = leaf_node.register_running(self.project, self.sample, self.tree)

                _ = self.register_finished(leaf_node)

    def calculate_report_overlaps_runs(self):
        for node in self.current_nodes:
            for leaf in node.leaves:
                leaf_node = self.spawn_node_child(node, leaf)
                _ = leaf_node.register_running(self.project, self.sample, self.tree)
                leaf_node.run_reference_overlap_analysis()

    def do_nothing(self):
        pass

    def deploy_nodes(self):
        """
        Deploy nodes according to pipeline sstep.
        """

        map_actions = {
            ConstantsSettings.PIPELINE_NAME_extra_qc: self.run_nodes_sequential,
            ConstantsSettings.PIPELINE_NAME_read_classification: self.run_simplified_classification_reads,
            ConstantsSettings.PIPELINE_NAME_contig_classification: self.run_simplified_classification_contigs,
            ConstantsSettings.PIPELINE_NAME_viral_enrichment: self.run_nodes_sequential,
            ConstantsSettings.PIPELINE_NAME_host_depletion: self.run_nodes_sequential,
            ConstantsSettings.PIPELINE_NAME_assembly: self.run_nodes_sequential,
            ConstantsSettings.PIPELINE_NAME_remapping: self.run_nodes_remap,
            ConstantsSettings.PIPELINE_NAME_remap_filtering: self.run_nodes_sequential,
            ConstantsSettings.PIPELINE_NAME_map_filtering: self.run_nodes_sequential,
            ConstantsSettings.PIPELINE_NAME_metagenomics_screening: self.run_nodes_sequential,
        }

        if self.current_module in ["end"]:
            return

        if self.current_module == "root":
            self.register_node_leaves(self.current_nodes[0])
            self.update_tree_nodes()
            return

        self.logger.info(f"CURRENT MODULE, {self.current_module}")
        action = map_actions[self.current_module]

        action()

        for node in self.current_nodes:
            if (
                self.classification_monitor.ready_to_merge(node)
                and node.run_manager.classification_updated == False
            ):
                node.run_manager.run_engine.plan_remap_prep_safe()

            self.update_node_leaves_dbs(node)
            self.register_node_leaves(node)

        self.update_tree_nodes()

        return

    def run_current_nodes_batch_parallel(self, batch=2):
        """
        run nodes in parallel batches
        """
        import multiprocessing as mp

        node_batches = [
            self.current_nodes[i : i + batch]
            for i in range(0, len(self.current_nodes), batch)
        ]

        for node_batch in node_batches:
            nproc = len(node_batch)
            pool = mp.Pool(nproc)
            drones = [
                pool.apply_async(node.run_manager.run_main) for node in node_batch
            ]

            for drone in drones:
                drone.get()

            for node in node_batch:
                self.register_node_leaves(node)

            pool.close()
            pool.join()

    def cycle_process(self):
        current_module = self.get_current_module()

        while current_module != "end":
            self.deploy_nodes()

            current_module = self.get_current_module()

        self.register_leaves_finished()

        print("DONE")
        return

    def stacked_changes_log(self):
        """
        generate stacked dataframe of software by node by pipeline step.
        """
        current_module = self.get_current_module()
        modules_list = [ConstantsSettings.PIPELINE_NAME_read_quality_analysis]
        parent_dict = self.tree.get_parents_dict()

        software_dict = {leaf: [] for leaf in self.tree.leaves}
        software_dict = {parent_dict[leaf]: [] for leaf in self.tree.leaves}

        while current_module != "end":
            if self.current_module in [
                "root",
            ]:
                self.update_tree_nodes()
            nodes_by_sample_sources = self.group_nodes_by_source_and_parameters()
            for combination in nodes_by_sample_sources:
                for node in combination:
                    software = node.parameters["software"].values[0]
                    for leaf in node.leaves:
                        software_dict[leaf].append(software)

            self.update_tree_nodes()
            current_module = self.get_current_module()
            if current_module != "end":
                modules_list.append(current_module)

        stacked_df = [software_list for leaf, software_list in software_dict.items()]

        stacked_df = pd.DataFrame(stacked_df, columns=modules_list)
        stacked_df["leaves"] = self.tree.leaves

        return stacked_df


class TreeProgressGraph:
    progress_stackl_df_name = "{}_stacked_df.tsv"
    progress_graph_html_name = "{}_graph.html"

    def __init__(self, sample: PIProject_Sample):
        self.sample = sample
        self.media_dir = sample.media_dir
        self.project = sample.project
        self.stacked_df_path = os.path.join(
            self.media_dir, self.progress_stackl_df_name.format(sample.pk)
        )
        self.graph_html_path = os.path.join(
            self.media_dir, self.progress_graph_html_name.format(sample.pk)
        )

        self.pipeline_utils = Utils_Manager()

    def get_node_params_network(
        self, existing_parameter_sets: Union[QuerySet, List[ParameterSet]]
    ) -> pd.DataFrame:

        technologies = [ps.project.technology for ps in existing_parameter_sets]
        if len(set(technologies)) > 1:
            raise Exception("Multiple technologies found")

        parameter_makeups = [ps.leaf.software_tree.pk for ps in existing_parameter_sets]
        parameter_makeups = list(set(parameter_makeups))

        tree_list = [ps.leaf.software_tree for ps in existing_parameter_sets]
        trees_pk_list = [tree.pk for tree in tree_list]
        trees_pk_list = list(set(trees_pk_list))

        software_tree_dict = {
            tree_pk: SoftwareTree.objects.get(pk=tree_pk) for tree_pk in trees_pk_list
        }

        pipetrees_dict = {
            tree_pk: self.pipeline_utils.parameter_util.convert_softwaretree_to_pipeline_tree(
                tree
            )
            for tree_pk, tree in software_tree_dict.items()
        }

        pipeline_steps_to_r_colours = {
            "root2": "lightblue",
            ConstantsSettings.PIPELINE_NAME_read_quality_analysis: "cadetblue",
            ConstantsSettings.PIPELINE_NAME_extra_qc: "cadetblue",
            ConstantsSettings.PIPELINE_NAME_viral_enrichment: "darkgreen",
            ConstantsSettings.PIPELINE_NAME_host_depletion: "blueviolet",
            ConstantsSettings.PIPELINE_NAME_assembly: "brown",
            ConstantsSettings.PIPELINE_NAME_contig_classification: "darkorange",
            ConstantsSettings.PIPELINE_NAME_read_classification: "deeppink",
            ConstantsSettings.PIPELINE_NAME_metagenomics_screening: "dodgerblue",
            ConstantsSettings.PIPELINE_NAME_map_filtering: "dodgerblue",
            ConstantsSettings.PIPELINE_NAME_remapping: "khaki",
            ConstantsSettings.PIPELINE_NAME_remap_filtering: "darkslategray",
            ConstantsSettings.PIPELINE_NAME_metagenomics_screening: "dodgerblue",
            ConstantsSettings.PIPELINE_NAME_request_mapping: "khaki",
            ConstantsSettings.PIPELINE_NAME_map_filtering: "darkslategray",
            ConstantsSettings.PIPELINE_NAME_metagenomics_screening: "dodgerblue",
            "Combined analysis": "darkslategray",
            "leaves2": "lightblue",
        }

        def merge_names(row: pd.Series):
            parent = row["parent"]
            child = row["child"]

            parent = node_dict[parent]
            child = node_dict[child]
            if parent == "NA":
                row["parent"] = "NA"
            else:
                parent_software = row["software_parent"]
                if parent_software is None:
                    parent_software = "input"
                row["parent"] = f"{parent_software}_{parent}"
            row["child"] = f"{row['software_child']}_{child}"
            return row

        network_df = [["NA", "0", "root2", "input", "input", "lightblue"]]
        for tree_pk, tree in pipetrees_dict.items():
            #
            tree.compress_tree()
            tree.split_modules()
            tree.get_module_tree()

            # network_df = [["NA", "0", "root", "input", "input", "lightblue"]]
            for edge in tree.edge_compress:
                parent_actual_node = [
                    x for x in tree.nodes_compress if x[0] == edge[0]
                ][0][1][0]
                child_group = [x for x in tree.nodes_compress if x[0] == edge[1]][0][1]

                if len(child_group) == 0:
                    continue
                child_actual_node = child_group[0]

                child_metadata = tree.node_index.loc[child_actual_node].node

                if child_metadata[-1] != "module":
                    continue

                parent_metadata = tree.node_index.loc[parent_actual_node].node

                parent = edge[0]
                child = edge[1]
                if parent == 0:
                    parent = "0"
                else:
                    parent = f"{tree_pk}_{parent}"
                child = f"{tree_pk}_{child}"
                module = child_metadata[0]
                colour = pipeline_steps_to_r_colours[module]
                software_child = child_metadata[1]
                software_parent = parent_metadata[1]
                network_df.append(
                    [parent, child, module, software_parent, software_child, colour]
                )

                def find_leaf_in_group(group):
                    for node in group:
                        if node in tree.leaves:
                            return node
                    return None

                leaf_in_group = find_leaf_in_group(child_group)

                if leaf_in_group:
                    color = pipeline_steps_to_r_colours["leaves2"]
                    network_df.append(
                        [
                            child,
                            f"{tree_pk}_{leaf_in_group}",
                            "leaves2",
                            software_child,
                            f"workflow {leaf_in_group}",
                            color,
                        ]
                    )

            network_df = pd.DataFrame(
                network_df,
                columns=[
                    "parent",
                    "child",
                    "module",
                    "software_parent",
                    "software_child",
                    "colour",
                ],
            )

            unique_nodes = list(network_df["child"].values) + list(
                network_df["parent"].values
            )
            unique_nodes = list(set(unique_nodes))

            # replace node names with numbers paste to software

            node_dict = {node: i for i, node in enumerate(unique_nodes)}
            node_dict["NA"] = "NA"

            network_df = network_df.apply(merge_names, axis=1)

        return network_df

    def get_node_params(
        self, existing_parameter_sets: Union[QuerySet, List[ParameterSet]]
    ) -> pd.DataFrame:
        """
        setup the trees for the progress graph
        """
        pipeline_utils = Utils_Manager()

        technologies = [ps.project.technology for ps in existing_parameter_sets]
        if len(set(technologies)) > 1:
            raise Exception("Multiple technologies found")

        parameter_makeups = [ps.leaf.software_tree.pk for ps in existing_parameter_sets]
        parameter_makeups = list(set(parameter_makeups))

        tree_list = [ps.leaf.software_tree for ps in existing_parameter_sets]
        trees_pk_list = [tree.pk for tree in tree_list]
        trees_pk_list = list(set(trees_pk_list))

        software_tree_dict = {
            tree_pk: SoftwareTree.objects.get(pk=tree_pk) for tree_pk in trees_pk_list
        }

        # pipetrees_dict = {
        #    tree_pk: pipeline_utils.parameter_util.convert_softwaretree_to_pipeline_tree(
        #        tree
        #    )
        #    for tree_pk, tree in software_tree_dict.items()
        # }

        stacked_df_dict = {}

        for ps in existing_parameter_sets:
            node_leaves = ps.get_leaf_descendants()
            leaf_node = node_leaves[0]
            leaf_node_index = leaf_node.index
            node_params = self.pipeline_utils.get_leaf_parameters(leaf_node)
            node_params = node_params[["module", "software"]]
            node_params = node_params.set_index("module")
            node_params = node_params.T
            node_params["leaves"] = leaf_node_index

            stacked_df_dict[ps.pk] = node_params

        # concatenate the stacked dfs
        ## columns are not the same.
        column_order = [
            ConstantsSettings.PIPELINE_NAME_read_quality_analysis,
            ConstantsSettings.PIPELINE_NAME_extra_qc,
            ConstantsSettings.PIPELINE_NAME_viral_enrichment,
            ConstantsSettings.PIPELINE_NAME_host_depletion,
            ConstantsSettings.PIPELINE_NAME_assembly,
            ConstantsSettings.PIPELINE_NAME_contig_classification,
            ConstantsSettings.PIPELINE_NAME_read_classification,
            ConstantsSettings.PIPELINE_NAME_metagenomics_screening,
            ConstantsSettings.PIPELINE_NAME_remapping,
            "leaves",
        ]

        all_columns = set()
        for _, stacked_df in stacked_df_dict.items():
            all_columns.update(stacked_df.columns)

        all_columns = list(all_columns)
        all_columns = [column for column in column_order if column in all_columns]

        for _, stacked_df in stacked_df_dict.items():
            missing_columns = [
                column for column in all_columns if column not in stacked_df.columns
            ]
            for column in missing_columns:
                stacked_df[column] = np.nan

        stacked_df_dict = {
            makeup: stacked_df[all_columns]
            for makeup, stacked_df in stacked_df_dict.items()
        }

        stacked_df = pd.concat(stacked_df_dict.values(), axis=0)
        stacked_df.to_csv(self.stacked_df_path, sep="\t")

        return stacked_df

    def get_combined_progress_df(self):
        ## setup a deployment and record the progress
        existing_parameter_sets = ParameterSet.objects.filter(
            project=self.project,
            status__in=[
                ParameterSet.STATUS_FINISHED,
            ],
            run_main__run_type=RunMain.RUN_TYPE_PIPELINE,
            sample=self.sample,
        )

        current_status = {ps.pk: ps.status for ps in existing_parameter_sets}

        #
        if existing_parameter_sets.count() == 0:
            return pd.DataFrame()

        #
        test_df = self.get_node_params_network(existing_parameter_sets)
        #

        for ps in existing_parameter_sets:
            ps.status = current_status[ps.pk]
            ps.save()

        test_df.to_csv(self.stacked_df_path, sep="\t", index=False)

        return test_df

    def get_tree_progress_df(self, tree: PipelineTree):
        ## setup a deployment and record the progress

        deployment_tree = Tree_Progress(tree, self.sample, self.project)

        stacked_df = deployment_tree.stacked_changes_log()
        #

        return stacked_df

    @staticmethod
    def extract_graph_data(html_filepath) -> Optional[str]:
        """
        graph data is stored in line cotaining: data-for="""

        with open(html_filepath, "r") as f:
            lines = f.readlines()

        for line in lines:
            if "data-for=" in line:
                return line

        return None

    def generate_graph(self):
        stacked_df = self.get_combined_progress_df()
        if stacked_df.shape[0] == 0:
            if os.path.exists(self.graph_html_path):
                os.remove(self.graph_html_path)
            return

        Rgraph_cmd = [
            Televir_Metadata.BINARIES["ROOT"]
            + Televir_Metadata.BINARIES["software"]["collapsibleTree"]
            + "/bin/"
            + "Rscript",
            "--vanilla",
            os.path.join(STATIC_ROOT, "R", "pipeline_dendrograph_network.R"),
            self.stacked_df_path,
            self.graph_html_path,
            ",".join(stacked_df.columns).replace(" ", "."),
        ]

        _ = os.system(" ".join(Rgraph_cmd))

    @staticmethod
    def extract_graph_id(graph_data: str) -> Optional[str]:
        """
        graph data is stored in line cotaining: data-for="""

        if graph_data is None:
            return None

        graph_id = graph_data.split("data-for=")[1].split('"')[1]

        return graph_id

    @staticmethod
    def extract_graph_json(graph_data: str) -> Optional[str]:
        """
        graph data is stored in line cotaining: data-for="""

        if graph_data is None:
            return None

        graph_json = graph_data.split(">")[1].split("<")[0]

        return graph_json

    def get_graph_data(self) -> Tuple[Optional[str], Optional[str]]:
        if not os.path.exists(self.graph_html_path):
            return None, None

        graph_data = self.extract_graph_data(self.graph_html_path)

        if graph_data is None:
            return None, None

        graph_id = self.extract_graph_id(graph_data)
        graph_json = self.extract_graph_json(graph_data)
        if graph_json is not None:
            graph_json = re.sub(r"_\d+", "", graph_json)

        return graph_json, graph_id
