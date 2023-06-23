import copy
import logging
import os
import shutil
from copy import _copy_immutable, _deepcopy_dispatch
from typing import List

import pandas as pd

from constants.constants import Televir_Metadata_Constants as Televir_Metadata
from constants.constants import TypePath
from pathogen_identification.constants_settings import ConstantsSettings as PIConstants
from pathogen_identification.models import (
    ParameterSet,
    PIProject_Sample,
    Projects,
    SoftwareTree,
    SoftwareTreeNode,
    RunMain,
    FinalReport,
)
from pathogen_identification.utilities.televir_parameters import TelevirParameters
from pathogen_identification.utilities.utilities_views import ReportSorter

from pathogen_identification.modules.remap_class import Mapping_Instance
from pathogen_identification.modules.run_main import RunMainTree_class
from pathogen_identification.utilities.update_DBs_tree import (
    Update_Assembly,
    Update_Classification,
    Update_Remap,
    Update_RunMain_Initial,
    Update_RunMain_Secondary,
    get_run_parents,
)
from pathogen_identification.utilities.utilities_pipeline import (
    PipelineTree,
    Utils_Manager,
)
from settings.constants_settings import ConstantsSettings
from utils.utils import Utils
from fluwebvirus.settings import STATIC_ROOT
from typing import Optional, Tuple


def logger_copy(x, memo):
    return x


copy._deepcopy_dispatch[logging.Logger] = logger_copy


class PathogenIdentification_Deployment_Manager:
    project: str
    prefix: str
    rdir: str
    threads: int
    run_engine: RunMainTree_class
    params = dict
    run_params_db = pd.DataFrame()
    pk: int = 0
    username: str
    prepped: bool = False
    sent: bool = False
    parameter_set: ParameterSet

    STATUS_ZERO = 0
    STATUS_PREPPED = 1
    STATUS_RUNNING = 2
    STATUS_SENT = 3

    def __init__(
        self,
        sample: PIProject_Sample,  # sample name
        project: Projects,
        username: str = "admin",
        technology: str = "ONT",
        deployment_root_dir: str = "/tmp/insaflu/insaflu_something",
        dir_branch: str = "deployment",
        threads: int = 3,
    ) -> None:
        self.username = username
        self.project = project.name
        self.sample = sample.name

        self.deployment_root_dir = deployment_root_dir
        self.dir_branch = dir_branch
        self.dir = os.path.join(self.deployment_root_dir, dir_branch)

        self.technology = technology
        self.install_registry = Televir_Metadata

        self.threads = threads
        self.file_r1 = sample.sample.get_fastq_available(TypePath.MEDIA_ROOT, True)
        if sample.sample.exist_file_2():
            self.file_r2 = sample.sample.get_fastq_available(TypePath.MEDIA_ROOT, False)
        else:
            self.file_r2 = ""

    def input_read_project_path(self, filepath) -> str:
        """copy input reads to project directory and return new path"""

        if not os.path.isfile(filepath):
            return ""

        rname = os.path.basename(filepath)
        new_rpath = os.path.join(self.dir, "reads") + "/" + rname
        shutil.copy(filepath, new_rpath)
        return new_rpath

    def configure(self) -> bool:
        """generate config dictionary for run_main, and copy input reads to project directory."""
        self.get_constants()

        self.generate_config_file()
        self.prep_test_env()

        new_r1_path = self.input_read_project_path(self.file_r1)
        new_r2_path = self.input_read_project_path(self.file_r2)

        self.config["sample_name"] = self.sample
        self.config["r1"] = new_r1_path
        self.config["r2"] = new_r2_path
        self.config["type"] = ["SE", "PE"][int(os.path.isfile(self.config["r2"]))]

        return True

    def get_constants(self):
        """set constants for technology"""
        if self.technology == ConstantsSettings.TECHNOLOGY_illumina:
            self.constants = PIConstants.CONSTANTS_ILLUMINA
        if self.technology == ConstantsSettings.TECHNOLOGY_minion:
            self.constants = PIConstants.CONSTANTS_ONT

    def generate_config_file(self):
        self.config = {
            "project": self.project,
            "source": self.install_registry.SOURCE,
            "deployment_root_dir": self.deployment_root_dir,
            "sub_directory": self.dir_branch,
            "directories": {},
            "threads": self.threads,
            "prefix": self.prefix,
            "project_name": self.project,
            "metadata": {
                x: os.path.join(self.install_registry.METADATA["ROOT"], g)
                for x, g in self.install_registry.METADATA.items()
            },
            "technology": self.technology,
            "bin": self.install_registry.BINARIES,
            "actions": {},
        }

        for dr, g in PIConstants.DIRS.items():
            self.config["directories"][dr] = os.path.join(self.dir, g)

        for dr, g in PIConstants.ACTIONS.items():
            self.config["actions"][dr] = g

        self.config.update(self.constants)

    def update_config(self):
        self.config["prefix"] = self.prefix

    def prep_test_env(self):
        """
        from main directory bearing scripts, params.py and main.sh, create metagenome run directory

        :return:
        """
        os.makedirs(self.dir, exist_ok=True)
        os.makedirs(
            os.path.join(PIConstants.media_directory, self.dir_branch),
            exist_ok=True,
        )
        os.makedirs(
            os.path.join(PIConstants.static_directory, self.dir_branch),
            exist_ok=True,
        )

        for directory in self.config["directories"].values():
            os.makedirs(directory, exist_ok=True)

    def close(self):
        if os.path.exists(self.dir):
            shutil.rmtree(self.dir)

    def import_params(self, run_params_db: pd.DataFrame):
        self.run_params_db = run_params_db

    def run_main_prep(self):
        if self.prepped or self.run_params_db.empty:
            return

        self.run_engine = RunMainTree_class(
            self.config, self.run_params_db, self.username
        )
        self.run_engine.Prep_deploy()
        self.prepped = True

    def run_main(self):
        self.run_engine.Run_QC()
        self.run_engine.Run_PreProcess()
        self.run_engine.Sanitize_reads()
        self.run_engine.Run_Assembly()
        self.run_engine.Run_Contig_classification()
        self.run_engine.Run_Read_classification()
        self.run_engine.plan_remap_prep_safe()
        # self.run_engine.Run_Classification()
        self.run_engine.Run_Remapping()

    def update_engine(self):
        self.update_config()
        self.run_engine.Update(self.config, self.run_params_db)

    def update_merged_targets(self, merged_targets: pd.DataFrame):
        self.run_engine.update_merged_targets(merged_targets)

    def delete_run_media(self):
        """delete project media directory"""

        if self.prepped:
            if os.path.isdir(self.run_engine.media_dir):
                shutil.rmtree(self.run_engine.media_dir, ignore_errors=True)

    def delete_run_static(self):
        """delete project static directory"""

        if self.prepped:
            if os.path.isdir(self.run_engine.static_dir):
                shutil.rmtree(self.run_engine.static_dir, ignore_errors=True)

    def delete_run_record(self, parameter_set: ParameterSet):
        """delete project record in database"""

        if self.prepped:
            _, runmain, _ = get_run_parents(self.run_engine, parameter_set)

            if runmain is not None:
                runmain.delete()

    def delete_run(self, parameter_set: ParameterSet):
        """delete project record in database"""

        self.delete_run_media()
        self.delete_run_static()
        self.delete_run_record(parameter_set)


class Tree_Node:
    module: str
    name: str
    node_index: int
    children: list
    parameters: pd.DataFrame
    software_tree_pk: int
    run_manager: PathogenIdentification_Deployment_Manager
    tree_node: SoftwareTreeNode
    parameter_set: ParameterSet

    def __init__(self, pipe_tree: PipelineTree, node_index: int, software_tree_pk: int):
        node_metadata = pipe_tree.node_index.loc[node_index].node

        self.module = node_metadata[0]
        self.node_index = node_index

        self.branch = pipe_tree.nodes_df.loc[node_index].branch
        self.children = pipe_tree.edge_df[
            pipe_tree.edge_df.parent == node_index
        ].child.tolist()

        self.parameters = self.determine_params(pipe_tree)
        self.software_tree_pk = software_tree_pk
        self.leaves = pipe_tree.leaves_from_node_compress(node_index)

    def run_reference_overlap_analysis(self):
        run = RunMain.objects.get(parameter_set=self.parameter_set)
        final_report = FinalReport.objects.filter(
            sample=self.parameter_set.sample, run=run
        ).order_by("-coverage")
        #
        print(final_report)
        report_layout_params = TelevirParameters.get_report_layout_params(run_pk=run.pk)
        report_sorter = ReportSorter(final_report, report_layout_params)
        report_sorter.sort_reports()

    def receive_run_manager(
        self, run_manager: PathogenIdentification_Deployment_Manager
    ):
        run_manager.prefix = f"run_leaf_{self.node_index}"
        run_manager.configure()
        run_manager.import_params(self.parameters)

        run_manager.run_main_prep()

        self.run_manager = run_manager

    def _is_node_leaf(self):
        return len(self.children) == 0

    def generate_software_tree_node_entry(self, pipe_tree: PipelineTree):
        if not self._is_node_leaf():
            return

        node_metadata = pipe_tree.node_index.loc[self.node_index].node
        software_tree = SoftwareTree.objects.get(pk=self.software_tree_pk)

        tree_node = SoftwareTreeNode.objects.filter(
            software_tree=software_tree,
            index=self.node_index,
            name=node_metadata[0],
            value=node_metadata[1],
            node_type=node_metadata[2],
        )

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

        return parameter_set

    def register_finished(
        self, project: Projects, sample: PIProject_Sample, tree: PipelineTree
    ):
        tree_node = self.generate_software_tree_node_entry(tree)
        if tree_node is None:
            return False

        parameter_set = self.setup_parameterset(project, sample, tree_node)

        if parameter_set is None:
            return False

        print("Updating run status FINISHED")
        print(parameter_set)

        parameter_set.status = ParameterSet.STATUS_FINISHED
        parameter_set.save()

        return True

    def register_failed(
        self, project: Projects, sample: PIProject_Sample, tree: PipelineTree
    ):
        tree_node = self.generate_software_tree_node_entry(tree)
        if tree_node is None:
            return False

        parameter_set = self.setup_parameterset(project, sample, tree_node)

        if parameter_set is None:
            return False

        parameter_set.status = ParameterSet.STATUS_ERROR
        parameter_set.save()

        self.run_manager.delete_run(parameter_set)

        return True

    def register(self, project: Projects, sample: PIProject_Sample, tree: PipelineTree):
        tree_node = self.generate_software_tree_node_entry(tree)
        if tree_node is None:
            return False

        parameter_set = self.setup_parameterset(project, sample, tree_node)

        if parameter_set is None:
            return False

        print("SETTING STATUS TO RUNNING")
        print(parameter_set)
        parameter_set.status = ParameterSet.STATUS_RUNNING
        parameter_set.save()

        self.tree_node = tree_node
        self.parameter_set = parameter_set

        return True

    def determine_params(self, pipe_tree):
        arguments_list = []
        for node in self.branch:
            node_metadata = pipe_tree.node_index.loc[node].node
            arguments_list.append(node_metadata)

        arguments_df = pd.DataFrame(
            arguments_list, columns=["parameter", "value", "flag"]
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


class PathogenDeployment_Iterator:
    def __init__(self):
        self.deplyment_list = []

    def add_deployment(self, deployment: PathogenIdentification_Deployment_Manager):
        self.deplyment_list.append(deployment)

    def __getitem__(self, index):
        return self.deplyment_list[index]

    def __len__(self):
        return len(self.deplyment_list)


class TreeNode_Iterator:
    def __init__(self):
        self.node_list = []

    def add_node(self, node: Tree_Node):
        self.node_list.append(node)

    def __getitem__(self, index) -> Tree_Node:
        return self.node_list[index]

    def __len__(self):
        return len(self.node_list)


class Tree_Progress:
    tree: PipelineTree
    current_nodes: List[Tree_Node]
    current_module: str

    sample: PIProject_Sample
    project: Projects

    def __init__(
        self, pipe_tree: PipelineTree, sample: PIProject_Sample, project: Projects
    ):
        pipe_tree.nodes_df = pd.DataFrame(
            pipe_tree.nodes_compress, columns=["node", "branch"]
        ).set_index("node")
        pipe_tree.edge_df = pd.DataFrame(
            pipe_tree.edge_compress, columns=["parent", "child"]
        )

        self.logger = logging.getLogger(__name__)

        self.tree = pipe_tree
        self.sample = sample
        self.project = project

        self.initialize_nodes()
        self.determine_current_module()

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

        deployment_manager = PathogenIdentification_Deployment_Manager(
            self.sample,
            self.project,
            self.project.owner.username,
            self.project.technology,
            temp_dir,
            deployment_directory_structure,
            PIConstants.DEPLOYMENT_THREADS,
        )

        return deployment_manager

    def register_node_leaves(self, node: Tree_Node):
        print("registering node: ", node.node_index, "leaves: ", node.leaves)

        if node.node_index in node.leaves:
            print(node, node.leaves)
            # self.submit_node_run(node)
            self.register_finished(node)

        for leaf in node.leaves:
            leaf_node = self.spawn_node_child(node, leaf)
            self.submit_node_run(leaf_node)

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
            self.tree, 0, software_tree_pk=self.tree.software_tree_pk
        )

        run_manager = self.setup_deployment_manager()

        origin_node.receive_run_manager(run_manager)

        self.register_node_leaves(origin_node)

        self.current_nodes = [origin_node]

    def get_current_module(self):
        return self.current_module

    def determine_current_module(self):
        self.current_module = self.current_nodes[0].module

    def spawn_node_child(self, node: Tree_Node, child: int):
        new_node = Tree_Node(self.tree, child, node.software_tree_pk)
        # node.run_manager.run_engine.logger = None

        run_manager_copy = copy.deepcopy(node.run_manager)
        new_node.receive_run_manager(run_manager_copy)
        new_node.run_manager.update_engine()

        return new_node

    @staticmethod
    def disable_logger(logger: logging.Logger):
        logger.disabled = True

    @staticmethod
    def enable_logger(logger: logging.Logger):
        logger.disabled = False

    @staticmethod
    def copy_class_instance(instance):
        return copy.deepcopy(instance)

    def register_node(self, node: Tree_Node):
        print("registering node")

        if node.run_manager.sent:
            return False

        registraction_success = node.register(self.project, self.sample, self.tree)

        return registraction_success

    def update_node_dbs(self, node: Tree_Node, step="initial"):
        try:
            db_updated = Update_RunMain_Initial(
                node.run_manager.run_engine, node.parameter_set
            )
            if not db_updated:
                return False

            if (
                node.run_manager.run_engine.enrichment_performed
                or node.run_manager.run_engine.depletion_performed
            ):
                db_updated = Update_RunMain_Secondary(
                    node.run_manager.run_engine, node.parameter_set
                )
                if not db_updated:
                    return False

            if node.run_manager.run_engine.assembly_performed:
                print("##### UPDATING ASSEMBLY DBS ######")
                print(node.run_manager.run_engine.assembly_drone.assembly_method.name)
                db_updated = Update_Assembly(
                    node.run_manager.run_engine, node.parameter_set
                )
                if not db_updated:
                    return False

            if (
                node.run_manager.run_engine.read_classification_performed
                and node.run_manager.run_engine.contig_classification_performed
            ):
                print("##### UPDATING CLASSIFICATION DBS ######")
                print(step)
                db_updated = Update_Classification(
                    node.run_manager.run_engine, node.parameter_set, tag=step
                )
                if not db_updated:
                    return False

                self.updated_classification = True

            if node.run_manager.run_engine.remapping_performed:
                print("exporting remapping")
                node.run_manager.run_engine.export_sequences()
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
            return False

    def submit_node_run(self, node: Tree_Node):
        print("Submitting node run: " + str(node.node_index))

        registration_success = self.register_node(node)

        print("Registration success: " + str(registration_success))
        if not registration_success:
            return

        self.update_node_dbs(node, step=node.module)

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
        planned_nodes = []
        for n in nodes:
            n.run_manager.run_engine.plan_remap_prep()

        return nodes

    def merge_node_targets_list(self, targetdf_list: List[Tree_Node]):
        """
        Merge targets from a list of nodes

        :param targetdf_list: list of nodes
        """
        node_merged_targets = [
            n.run_manager.run_engine.merged_targets for n in targetdf_list
        ]

        node_merged_targets = pd.concat(node_merged_targets, axis=0).reset_index()

        return node_merged_targets

    def get_node_node_targets(self, nodes_list: List[Tree_Node]):
        """
        Get merged targets from a list of nodes

        :param nodes_list: list of nodes"""
        node_merged_targets = self.merge_node_targets_list(nodes_list)

        node_merged_targets = self.process_mapping_managerdf(node_merged_targets)

        return node_merged_targets

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
            sample_source = node.run_manager.run_engine.sample.sources_list()
            print("node: " + str(node.node_index))
            print("Sample source: " + str(sample_source))
            update_combination_dict(node)

        print("########### GROUPED NODES ###########")
        print(source_paramaters_combinations)

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

    def process_subject(self, volonteer: Tree_Node, group_targets: pd.DataFrame):
        original_targets = copy.deepcopy(
            volonteer.run_manager.run_engine.merged_targets
        )

        volonteer.run_manager.update_merged_targets(group_targets)

        run_success = self.run_node(volonteer)

        volonteer.run_manager.update_merged_targets(original_targets)

        mapped_instances_shared = (
            volonteer.run_manager.run_engine.remap_manager.mapped_instances
        )

        return mapped_instances_shared, run_success

    def stacked_deployment_classification(
        self, nodes_by_sample_sources: List[List[Tree_Node]]
    ):
        new_nodes = []
        print("########## STACKED DEPLOYMENT CLASSIFICATION ##########")
        print(nodes_by_sample_sources)
        for nodes_subset in nodes_by_sample_sources:
            if len(nodes_subset) == 0:
                continue

            volonteer = nodes_subset[0]
            print("volonteer: " + str(volonteer.node_index))
            run_success = self.run_node(volonteer)
            print("run_success: " + str(run_success))

            if run_success:
                for node in nodes_subset:
                    node.run_manager.run_engine.read_classification_drone = (
                        volonteer.run_manager.run_engine.read_classification_drone
                    )
                    node.run_manager.run_engine.read_classification_performed = True
                    new_nodes.append(node)
            else:
                for node in nodes_subset:
                    self.register_failed_children(node)

        if len(new_nodes) > 0:
            self.current_nodes = new_nodes

    def stacked_deployement(self, nodes_by_sample_sources: List[List[Tree_Node]]):
        """
        deploy nodes remap by sample sources."""
        current_nodes = []

        for nodes in nodes_by_sample_sources:
            if len(nodes) == 0:
                continue

            nodes = self.get_remap_plans(nodes)
            group_targets = self.get_node_node_targets(nodes)

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

        self.stacked_deployement(nodes_by_sample_sources)

    def run_simplified_classification(self):
        nodes_by_sample_sources = self.group_nodes_by_source_and_parameters()

        self.stacked_deployment_classification(nodes_by_sample_sources)

    def update_nodes(self):
        new_nodes = []
        for node in self.current_nodes:
            children = node.children
            for child in children:
                new_node = self.spawn_node_child(node, child)
                new_nodes.append(new_node)

        if len(new_nodes) == 0:
            self.current_module = "end"
        else:
            self.current_nodes = new_nodes
            self.determine_current_module()

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

    def run_node(self, node: Tree_Node):
        try:
            node.run_manager.run_main()
            return True
        except Exception as e:
            print("error")
            print(e)

            return False

    def run_nodes_sequential(self):
        self.run_current_nodes()
        self.register_current_nodes()
        self.update_nodes()

    def run_nodes_simply(self):
        self.run_simplified_mapping()
        self.register_current_nodes()
        self.update_nodes()

    def run_nodes_classification_reads(self):
        self.run_simplified_classification()
        print("RAN CLASSIFICATION")
        self.register_current_nodes()
        self.update_nodes()

    def run_nodes_classification_contigs(self):
        self.run_simplified_classification()
        print("RAN CLASSIFICATION")
        self.register_current_nodes()
        self.update_nodes()

    def register_current_nodes(self):
        for node in self.current_nodes:
            self.register_node_leaves(node)

    def register_leaves_finished(self):
        for node in self.current_nodes:
            for leaf in node.leaves:
                leaf_node = self.spawn_node_child(node, leaf)
                _ = leaf_node.register(self.project, self.sample, self.tree)

                success_register = self.register_finished(leaf_node)
                print("success_register")

    def calculate_report_overlaps(self):
        for node in self.current_nodes:
            for leaf in node.leaves:
                leaf_node = self.spawn_node_child(node, leaf)
                _ = leaf_node.register(self.project, self.sample, self.tree)
                leaf_node.run_reference_overlap_analysis()

    def deploy_nodes(self):
        if self.current_module == "end":
            return
        if self.current_module in [
            "root",
        ]:
            self.update_nodes()
            return

        if self.current_module == ConstantsSettings.PIPELINE_NAME_read_classification:
            self.run_nodes_classification_reads()
            return

        # if self.current_module == ConstantsSettings.PIPELINE_NAME_contig_classification:
        #    self.run_nodes_classification_reads()

        if self.current_module == ConstantsSettings.PIPELINE_NAME_remapping:
            self.run_nodes_simply()
        else:
            self.run_nodes_sequential()

        return

    def run_current_nodes_batch_parallel(self, batch=2):
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

        print(self.current_nodes)
        self.calculate_report_overlaps()
        self.register_leaves_finished()

        print("DONE")

    def stacked_changes_log(self):
        """
        generate stacked dataframe of software by node by pipeline step.
        """
        current_module = self.get_current_module()
        modules_list = [
            "QC",
        ]
        parent_dict = self.tree.get_parents_dict()

        software_dict = {leaf: [] for leaf in self.tree.leaves}
        software_dict = {parent_dict[leaf]: [] for leaf in self.tree.leaves}

        while current_module != "end":
            if self.current_module in [
                "root",
            ]:
                self.update_nodes()
            nodes_by_sample_sources = self.group_nodes_by_source_and_parameters()
            for combination in nodes_by_sample_sources:
                for node in combination:
                    software = node.parameters["software"].values[0]
                    for leaf in node.leaves:
                        software_dict[leaf].append(software)

            self.update_nodes()
            current_module = self.get_current_module()
            if current_module != "end":
                modules_list.append(current_module)

        stacked_df = [software_list for lead, software_list in software_dict.items()]

        stacked_df = pd.DataFrame(stacked_df, columns=modules_list)
        stacked_df["leaves"] = self.tree.leaves

        return stacked_df


class TreeProgressGraph:
    progress_stackl_df_name = "{}_stacked_df.tsv"
    progress_graph_html_name = "{}_graph.html"

    def __init__(self, sample: PIProject_Sample):
        self.sample = sample
        self.media_dir = sample.get_media_dir()
        self.project = sample.project
        self.stacked_df_path = os.path.join(
            self.media_dir, self.progress_stackl_df_name.format(sample.pk)
        )
        self.graph_html_path = os.path.join(
            self.media_dir, self.progress_graph_html_name.format(sample.pk)
        )

        self.pipeline_utils = Utils_Manager()

    def setup_tree(self):
        existing_parameter_sets = ParameterSet.objects.filter(
            project=self.project,
            status=ParameterSet.STATUS_RUNNING,
            sample=self.sample,
        )
        ## create a tree that contains all the parameter sets
        (
            combined_tree,
            additional_leaves,
        ) = self.pipeline_utils.get_parameterset_leaves_list(existing_parameter_sets)

        # CRUNCH TREE
        module_tree = self.pipeline_utils.module_tree(combined_tree, additional_leaves)
        return module_tree

    def get_progress_df(self):
        ## setup a deployment and record the progress
        module_tree = self.setup_tree()
        deployment_tree = Tree_Progress(module_tree, self.sample, self.project)

        stacked_df = deployment_tree.stacked_changes_log()
        stacked_df.to_csv(self.stacked_df_path, sep="\t")

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
        self.get_progress_df()

        Rgraph_cmd = [
            Televir_Metadata.BINARIES["software"]["collapsibleTree"] + "/bin/" + "Rscript",
            "--vanilla",
            os.path.join(STATIC_ROOT, "R", "pipeline_dendrograph.R"),
            self.stacked_df_path,
            self.graph_html_path,
        ]

        result = os.system(" ".join(Rgraph_cmd))

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
            self.generate_graph()

        graph_data = self.extract_graph_data(self.graph_html_path)

        if graph_data is None:
            return None, None

        graph_id = self.extract_graph_id(graph_data)
        graph_json = self.extract_graph_json(graph_data)

        return graph_json, graph_id
