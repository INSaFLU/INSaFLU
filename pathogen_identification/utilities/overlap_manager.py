import itertools
import itertools as it
import os
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from scipy.spatial.distance import pdist, squareform

from pathogen_identification.utilities.clade_objects import Clade, CladeFilter
from pathogen_identification.utilities.phylo_tree import PhyloTreeManager

## pairwise matrix by individual reads
from pathogen_identification.utilities.utilities_general import readname_from_fasta


def accid_from_metadata(metadata: pd.DataFrame, read_name: str) -> str:
    """
    Return accession id of read from metadata
    """
    filename = os.path.basename(read_name)
    try:
        return metadata[metadata["filename"] == filename]["accid"].values[0]
    except IndexError:
        return read_name


import logging


class ReadOverlapManager:
    distance_matrix_filename: str = "distance_matrix_{}.tsv"
    clade_statistics_filename: str = "clade_statistics_{}.tsv"
    accid_statistics_filename: str = "accid_statistics_{}.tsv"
    tree_plot_filename: str = "tree_{}.png"
    overlap_matrix_plot_filename: str = "overlap_matrix_{}.png"
    min_freq: float = 0.05
    max_reads: int = 100000

    def __init__(
        self,
        metadata_df: pd.DataFrame,
        reference_clade: Clade,
        media_dir: str,
        pid: str,
    ):
        self.metadata = metadata_df
        self.fasta_list = metadata_df["file"].tolist()
        self.clade_filter = CladeFilter(reference_clade=reference_clade)
        self.excluded_leaves = []
        self.media_dir = media_dir

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

        self.distance_matrix_path = os.path.join(
            self.media_dir, self.distance_matrix_filename.format(pid)
        )
        self.clade_statistics_path = os.path.join(
            self.media_dir, self.clade_statistics_filename.format(pid)
        )
        self.accid_statistics_path = os.path.join(
            self.media_dir, self.accid_statistics_filename.format(pid)
        )
        self.tree_plot_path = os.path.join(
            self.media_dir, self.tree_plot_filename.format(pid)
        )
        self.overlap_matrix_plot_path = os.path.join(
            self.media_dir, self.overlap_matrix_plot_filename.format(pid)
        )

        self.metadata["filename"] = self.metadata["file"].apply(
            lambda x: x.split("/")[-1]
        )

        self.metadata["filepath"] = self.metadata["file"]

        self.tree_manager = self.prep_tree_for_clade_analysis()
        if not os.path.exists(self.tree_plot_path):
            self.tree_manager.plot_tree(self.tree_plot_path)

        self.tree_plot_exists = os.path.exists(self.tree_plot_path)
        self.overlap_matrix_plot_exists = os.path.exists(self.overlap_matrix_plot_path)

    def all_accs_analyzed(self):
        if not os.path.exists(self.accid_statistics_path):
            return False

        accid_df = pd.read_csv(self.accid_statistics_path, sep="\t")

        for accid in self.metadata["accid"].unique():
            if accid not in accid_df["accid"].tolist():
                return False

        return True

    def parse_for_data(self):
        # self.read_names_dict: Dict[str, List[str]] = self.get_accid_readname_dict()
        self.read_profile_matrix: pd.DataFrame = self.generate_read_matrix()
        self.overlap_matrix: pd.DataFrame = self.pairwise_shared_count(
            self.read_profile_matrix
        )

        accid_df = self.metadata[["accid", "description"]]
        accid_df["read_count"] = accid_df["accid"].apply(
            lambda x: self.get_accession_total_counts(x)
        )
        # sort table by accid and then by read count
        accid_df = accid_df.sort_values(["accid", "read_count"], ascending=False)
        # drop duplicates of accid
        accid_df = accid_df.drop_duplicates(subset=["accid"])

        accid_df["proportion"] = accid_df["accid"].apply(
            lambda x: self.get_proportion_counts(x)
        )
        accid_df.to_csv(self.accid_statistics_path, sep="\t", index=False)

    def overlap_heatmap_plot(self) -> None:
        """
        Plot heatmap of read overlap between all pairs of lists
        """
        plt.figure(figsize=(15, 6))
        sns.heatmap(self.overlap_matrix, annot=True)
        plt.savefig(self.overlap_matrix_plot_path)

    def update_excluded_leaves(self, read_profile_matrix: pd.DataFrame):
        all_node_leaves = self.tree_manager.all_clades_leaves()

        for _, leaves in all_node_leaves.items():
            for leaf in leaves:
                if leaf not in read_profile_matrix.index:
                    self.excluded_leaves.append(leaf)

    def get_accid_readname_dict(self):
        """
        Return dictionary of read names and descriptions
        """
        readname_dict = {}
        for fasta_file in self.fasta_list:
            accid = accid_from_metadata(self.metadata, fasta_file)
            read_names = readname_from_fasta(fasta_file)

            readname_dict[accid] = read_names
        return readname_dict

    @staticmethod
    def readoverlap_2_files(lista, listb) -> list:
        """
        Return list of read names that are in both lists
        """
        return list(set(lista).intersection(set(listb)))

    def readoverlap_allpairs(
        self,
        read_lists: List[List[str]],
    ) -> Dict[Tuple[str, str], List[str]]:
        """
        Return dictionary of read overlap between all pairs of lists
        """
        read_overlap_dict = {}
        for i in range(len(read_lists)):
            for j in range(len(read_lists)):
                read_overlap = self.readoverlap_2_files(read_lists[i], read_lists[j])
                percent_i = len(read_overlap) / len(read_lists[i])
                percent_j = len(read_overlap) / len(read_lists[j])
                summary = f"{percent_i:.2f} - {percent_j:.2f}"

                read_overlap_dict[(i, j)] = percent_i

        return read_overlap_dict

    @staticmethod
    def all_reads_set(files_readnames: List[list]) -> list:
        all_reads = list(it.chain.from_iterable(files_readnames))
        all_reads = list(set(all_reads))

        return all_reads

    @staticmethod
    def read_profile_get(accid: str, readname_dict: dict, all_reads: list) -> list:
        """
        Return list of 1s and 0s for presence/absence of read in accid
        """
        acc_read_dict = {read: 1 for read in readname_dict[accid]}

        return [acc_read_dict.get(read, 0) for read in all_reads]

    def read_profile_dict_get(self, readname_dict: dict, all_reads: list) -> dict:
        """
        Return dictionary of read profiles for all accids
        """
        read_profile_dict = {}
        for accid in readname_dict.keys():
            read_profile_dict[accid] = self.read_profile_get(
                accid, readname_dict, all_reads
            )
        return read_profile_dict

    def readoverlap_allpairs_df(self):
        """
        Return dataframe of read overlap between all pairs of lists
        """
        proportions_matrix = self.read_overlap_proportions()

        read_overlap_as_pairs = [
            [
                proportions_matrix.index[i],
                proportions_matrix.columns[j],
                proportions_matrix.iloc[i, j],
            ]
            for i in range(len(proportions_matrix))
            for j in range(len(proportions_matrix.columns))
        ]

        pair_overlap_matrix = pd.DataFrame(
            data=read_overlap_as_pairs,
            columns=["accid_A", "accid_B", "overlap"],
        )

        return pair_overlap_matrix

    @staticmethod
    def read_profile_matrix_get(read_profile_dict: dict) -> pd.DataFrame:
        """
        Return dataframe of read profiles
        """
        return pd.DataFrame(read_profile_dict).T

    @staticmethod
    def pairwise_shared_reads(read_profile_matrix: pd.DataFrame) -> pd.DataFrame:
        """
        Return dataframe of pairwise shared reads
        """

        if read_profile_matrix.shape[1] == 0:
            return pd.DataFrame()

        return pd.DataFrame(
            squareform(pdist(read_profile_matrix, metric="jaccard")),
            columns=read_profile_matrix.index,
            index=read_profile_matrix.index,
        )

    @staticmethod
    def matrix_lower_triangle(matrix: pd.DataFrame) -> pd.DataFrame:
        """
        Return dataframe of lower triangle of matrix
        """
        return matrix.where(np.tril(np.ones(matrix.shape)).astype(bool))

    @staticmethod
    def matrix_to_phylotriangle(distance_matrix: pd.DataFrame):
        """
        Return tree from matrix
        """
        distmat = distance_matrix.values.tolist()
        distmat = [x[: i + 1] for i, x in enumerate(distmat)]
        distmat = DistanceMatrix(list(distance_matrix.index), distmat)

        return distmat

    def tree_from_distance_matrix(self, distance_matrix: pd.DataFrame):
        """
        Return tree from distance matrix
        """
        distmat = self.matrix_to_phylotriangle(distance_matrix)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(distmat)
        tree.rooted = True
        tree.ladderize()
        return tree

    ####################
    ## Construct tree ##
    ####################

    def filter_read_matrix(self, read_profile_matrix: pd.DataFrame) -> pd.DataFrame:
        """
        Filter read matrix, reads as columns, accids as rows
        """

        if read_profile_matrix.shape[1] == 0:
            self.logger.info("No reads with frequency > min_freq")
            return read_profile_matrix

        read_counts = read_profile_matrix.sum(axis=0)
        read_freqs = read_counts / read_profile_matrix.shape[0]
        # filter out reads that are only present in one accession
        read_profile_matrix_filtered = read_profile_matrix.loc[:, read_counts > 1]
        # filter reads with less than min_freq
        read_profile_matrix_filtered = read_profile_matrix_filtered.loc[
            :, read_freqs > self.min_freq
        ]

        if read_profile_matrix_filtered.shape[1] == 0:
            self.logger.info("No reads with frequency > min_freq")
            return read_profile_matrix

        if self.max_reads:
            if read_profile_matrix_filtered.shape[1] > self.max_reads:
                self.logger.info(
                    f"More than {self.max_reads} reads ({read_profile_matrix_filtered.shape[1]}) - sampling"
                )
                ## sample reads
                read_profile_matrix_filtered = read_profile_matrix_filtered.sample(
                    n=self.max_reads, axis=1
                )

        return read_profile_matrix_filtered

    def generate_read_matrix(self):
        """
        Generate read matrix
        """
        self.logger.info("generating read matrix")

        readname_dict = self.get_accid_readname_dict()
        all_reads = self.all_reads_set(list(readname_dict.values()))
        read_profile_dict = self.read_profile_dict_get(readname_dict, all_reads)
        read_profile_matrix = self.read_profile_matrix_get(read_profile_dict)
        read_profile_matrix = self.filter_read_matrix(read_profile_matrix)
        return read_profile_matrix

    def pairwise_shared_count(
        self,
        read_profile_matrix: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Return dataframe of pairwise shared read counts,
        use matrix multiplication to sum shared reads from binary matrix for each pair.
        """

        binary_matrix = np.array(read_profile_matrix)

        shared_reads = []
        for i in range(binary_matrix.shape[0]):
            row = binary_matrix[i].reshape(1, -1)
            prod0 = row @ binary_matrix.T

            shared_reads.append(prod0)

        if len(shared_reads) == 0:
            return pd.DataFrame(
                index=read_profile_matrix.index, columns=read_profile_matrix.index
            )
        shared_reads = np.concatenate(shared_reads, axis=0)

        shared_reads = pd.DataFrame(
            shared_reads,
            index=read_profile_matrix.index,
            columns=read_profile_matrix.index,
        )

        return shared_reads

    def read_overlap_proportions(self) -> pd.DataFrame:
        """
        Return dataframe of pairwise shared read proportions."""
        proxy_overlap = self.overlap_matrix

        def standardize_by_index_values(row: pd.Series) -> pd.Series:
            """ """
            # get value of first index
            row = row / row[row.name]
            return row

        proportions_matrix = proxy_overlap.apply(standardize_by_index_values, axis=1)
        ## fill upper triable with transposed lower triangle

        return proportions_matrix

    def reads_dict_from_matrix(self, read_profile_matrix: pd.DataFrame) -> dict:
        """
        Return dictionary of reads for each accession
        """
        reads_dict = {}
        for accid in read_profile_matrix.index:
            reads_dict[accid] = read_profile_matrix.columns[
                read_profile_matrix.loc[accid] == 1
            ].values.tolist()
        return reads_dict

    def get_accession_total_counts(self, accid: str):
        """
        Get total counts for accession
        """
        return self.read_profile_matrix.loc[accid].sum()

    def get_proportion_counts(self, accid: str):
        """
        Get proportion counts for accession
        """

        return (
            self.get_accession_total_counts(accid)
            / self.read_profile_matrix.sum().sum()
        )

    def check_all_accessions_in_distance_matrix(self, distance_matrix):
        for accid in self.metadata["accid"].values:
            if accid not in distance_matrix.columns:
                return False

        return True

    def generate_distance_matrix(self):
        """
        Generate distance matrix
        """
        if os.path.isfile(self.distance_matrix_path):
            distance_matrix = pd.read_csv(self.distance_matrix_path, index_col=0)
        else:
            self.parse_for_data()
            distance_matrix = self.pairwise_shared_reads(self.read_profile_matrix)

        if not self.check_all_accessions_in_distance_matrix(distance_matrix):
            self.parse_for_data()
            distance_matrix = self.pairwise_shared_reads(self.read_profile_matrix)

        try:  # Written only on job submisision. File not written on query.
            distance_matrix.to_csv(self.distance_matrix_path)
        except:
            pass

        return distance_matrix

    def generate_tree(self):
        """
        Generate tree
        """
        distance_matrix = self.generate_distance_matrix()

        tree = self.tree_from_distance_matrix(distance_matrix)
        return tree

    ####################
    ## private clades ##
    ####################

    def clade_shared_by_pair(self, leaves: list) -> Tuple[float, float, float]:
        group = self.read_profile_matrix.loc[leaves]

        group_pairwise_shared = self.pairwise_shared_count(group)

        # divide shared rows by group row sums
        group_pairwise_shared = group_pairwise_shared.div(group.sum(axis=1), axis=1)
        print(group_pairwise_shared.shape)

        # get lower triangle of shared
        lower_triangle = self.matrix_lower_triangle(group_pairwise_shared)
        # get upper triangle from matrix, without diagonal
        upper_triangle = group_pairwise_shared.where(
            np.triu(np.ones(group_pairwise_shared.shape), k=1).astype(bool)
        )

        # flatten both to list and append
        group_pairwise_shared = lower_triangle.values.flatten().tolist()
        group_pairwise_shared_copy = upper_triangle.values.flatten().tolist()

        group_pairwise_shared.extend(group_pairwise_shared_copy)
        # remove nan values
        group_pairwise_shared = [x for x in group_pairwise_shared if str(x) != "nan"]

        if "NC_021505.1" in leaves and "AP013070.1" in leaves:
            print("##########")
            print(leaves)
            print(group_pairwise_shared)

        min_shared = min(group_pairwise_shared)
        max_shared = max(group_pairwise_shared)
        std_shared = np.std(group_pairwise_shared)

        return min_shared, max_shared, std_shared

    def clade_shared_by_pair_old(self, leaves: list) -> pd.DataFrame:
        """
        return tuple of proportions of reads shared by each pair of leaves
        """
        overlap_df = self.readoverlap_allpairs_df()

        subset_clade = overlap_df.loc[
            (overlap_df.accid_A.isin(leaves)) & (overlap_df.accid_B.isin(leaves))
        ]

        subset_clade = subset_clade.loc[subset_clade.accid_A != subset_clade.accid_B]

        pair_proportions_df = []

        accids = subset_clade.accid_A.unique()
        combinations = itertools.combinations(accids, 2)

        for pair in combinations:
            pair_proportions = subset_clade.loc[
                (subset_clade.accid_A.isin(pair)) & (subset_clade.accid_B.isin(pair))
            ]

            pair_proportions_df.append(
                [pair[0], pair[1], pair_proportions.overlap.tolist()]
            )

        combinations = pd.DataFrame(
            pair_proportions_df, columns=["accid_A", "accid_B", "proportion_shared"]
        )

        combinations["proportion_max"] = combinations.proportion_shared.apply(
            lambda x: max(x)
        )
        combinations["proportion_min"] = combinations.proportion_shared.apply(
            lambda x: min(x)
        )
        combinations["proportion_std"] = combinations.proportion_shared.apply(
            lambda x: np.std(x)
        )
        return combinations

    def clade_total_counts(self, leaves: list) -> float:
        """
        return total counts of clade
        """
        return int(self.read_profile_matrix.loc[leaves].sum().sum())

    def clade_private_proportions_old(self, leaves: list) -> float:
        """ """
        group = self.read_profile_matrix.loc[leaves]
        group_sum = group.sum(axis=0)
        group_sum_as_bool = group_sum > 0
        group_sum_as_bool_list = group_sum_as_bool.tolist()

        group_reads = self.read_profile_matrix.iloc[:, group_sum_as_bool_list]

        group_reads_sum_all = group_reads.sum(axis=0)
        group_reads_sum_group = group_reads.loc[leaves].sum(axis=0)
        group_reads_sum_group = group_reads_sum_group.fillna(0)

        read_proportions = group_reads_sum_group / group_reads_sum_all
        read_proportions = read_proportions.fillna(0)

        # proportion of reads in group that are private
        proportion_private = read_proportions.sum() / len(read_proportions)

        return proportion_private

    def clade_private_proportions(self, leaves: list) -> float:
        """ """
        group = self.read_profile_matrix.loc[leaves]
        group_sum = group.sum(axis=0)
        group_sum_as_bool = group_sum > 0
        group_sum_as_bool_list = group_sum_as_bool.tolist()

        sum_all = self.read_profile_matrix.iloc[:, group_sum_as_bool_list].sum(axis=0)
        sum_group = self.read_profile_matrix.loc[leaves]
        sum_group = sum_group.iloc[:, group_sum_as_bool_list].sum(axis=0)

        private_reads = sum_group - sum_all

        private_reads = sum(private_reads == 0)

        print(
            leaves,
            private_reads,
            sum(group_sum_as_bool_list),
            private_reads / sum(group_sum_as_bool_list),
        )

        proportion_private = private_reads / sum(group_sum_as_bool_list)

        return proportion_private

    def clade_reads_matrix(self, filter_names=[], remove_leaves=True) -> pd.DataFrame:
        """
        Return dataframe reads per clade"""
        clade_read_matrix = []
        belonging = []
        for clade, leaves in self.all_clade_leaves_filtered.items():
            if len(leaves) == 0:
                continue

            # continue clade is leaf:
            if len(leaves) == 1 and remove_leaves:
                continue

            if filter_names:
                if clade.name not in filter_names:
                    continue

            reads_in_clade = self.read_profile_matrix.loc[leaves]
            reads_in_clade_sum = reads_in_clade.sum(axis=0)
            reads_in_clade_sum_as_bool = reads_in_clade_sum > 0
            reads_in_clade_sum_as_int_list = reads_in_clade_sum_as_bool.astype(
                int
            ).tolist()
            clade_read_matrix.append(reads_in_clade_sum_as_int_list)
            belonging.append(clade.name)

        clade_read_matrix = pd.DataFrame(
            clade_read_matrix,
            index=belonging,
            columns=self.read_profile_matrix.columns,
        )
        ## sort rows by row sum in descending order
        # clade_read_matrix = clade_read_matrix.loc[
        #    clade_read_matrix.sum(axis=1).sort_values(ascending=False).index
        # ]

        shared_clade_matrix = self.pairwise_shared_count(clade_read_matrix)

        ## divide rows of shared_clade_matrix by clade_read_matrix row sums
        shared_clade_matrix = shared_clade_matrix.div(
            clade_read_matrix.sum(axis=1), axis=0
        )

        # set diagonal to 1
        np.fill_diagonal(shared_clade_matrix.values, 0)

        return shared_clade_matrix

    def pairwise_clade_shared_reads(self, clades_filter=[]) -> pd.DataFrame:
        """
        Return dataframe of pairwise shared reads between all pairs of clades"""
        clade_read_matrix = self.clade_reads_matrix(
            filter_names=clades_filter, remove_leaves=False
        )

        return clade_read_matrix

    def plot_pairwise_shared_clade_reads(self, clades_filter=[]) -> None:
        """
        Plot heatmap of pairwise shared reads between all pairs of leaves
        """
        pairwise_shared_clade = self.pairwise_clade_shared_reads(
            clades_filter=clades_filter
        )
        plt.figure(figsize=(15, 6))
        sns.heatmap(pairwise_shared_clade, annot=True)
        # center figure
        plt.tight_layout()
        plt.savefig(self.overlap_matrix_plot_path)

    def node_statistics(self) -> dict:
        self.parse_for_data()
        self.update_excluded_leaves(self.read_profile_matrix)

        node_stats_dict = {}

        if self.read_profile_matrix.shape[1] == 0:
            self.logger.info("No reads with frequency > min_freq")
            return node_stats_dict

        for node, leaves in self.all_clade_leaves_filtered.items():
            if len(leaves) == 0:
                node_stats_dict[node] = Clade(
                    name=node,
                    leaves=leaves,
                    group_counts=0,
                    private_proportion=0,
                    shared_proportion_std=0,
                    shared_proportion_min=0,
                    shared_proportion_max=0,
                )

                continue

            proportion_private = self.clade_private_proportions(leaves)
            clade_counts = self.clade_total_counts(leaves)

            if len(leaves) == 1:
                node_stats_dict[node] = Clade(
                    name=node,
                    leaves=leaves,
                    group_counts=clade_counts,
                    private_proportion=proportion_private,
                    shared_proportion_std=0,
                    shared_proportion_min=0,
                    shared_proportion_max=0,
                )

                continue

            combinations = self.clade_shared_by_pair_old(leaves)

            node_stats_dict[node] = Clade(
                name=node,
                leaves=leaves,
                private_proportion=proportion_private,
                group_counts=clade_counts,
                shared_proportion_min=min(combinations.proportion_max),
                shared_proportion_max=max(combinations.proportion_max),
                shared_proportion_std=np.std(combinations.proportion_max),
            )

        return node_stats_dict

    def all_clades_summary(
        self, node_stats_dict: Dict[Phylo.BaseTree.Clade, Clade]
    ) -> pd.DataFrame:
        clade_summary = []
        for node, stats in node_stats_dict.items():
            clade_summary.append(
                [
                    node,
                    stats.private_proportion,
                    stats.group_counts,
                    stats.shared_proportion_std,
                    stats.shared_proportion_min,
                    stats.shared_proportion_max,
                ]
            )

        clade_summary = pd.DataFrame(
            clade_summary,
            columns=[
                "node",
                "private_proportion",
                "total_counts",
                "shared_proportion_std",
                "shared_proportion_min",
                "shared_proportion_max",
            ],
        )

        return clade_summary

    def clade_dict_from_summary(
        self,
        clade_summary: pd.DataFrame,
        tree,
        inner_node_leaf_dict: Dict[Phylo.BaseTree.Clade, list],
    ) -> Dict[Phylo.BaseTree.Clade, Clade]:
        """
        Create a dictionary of clades from the summary dataframe
        """
        node_stats_dict = {}

        # set the index to the node name
        clade_summary = clade_summary.set_index("node")

        for node, stats in clade_summary.iterrows():
            try:
                node_phylo_clade = [
                    x for x in inner_node_leaf_dict.keys() if x.name == node
                ][0]
            except IndexError:
                # self.logger.error("node not found in tree")
                continue

            node_stats_dict[node_phylo_clade] = Clade(
                name=str(node),
                leaves=inner_node_leaf_dict[node_phylo_clade],
                private_proportion=stats.private_proportion,
                group_counts=stats.total_counts,
                shared_proportion_std=stats.shared_proportion_std,
                shared_proportion_min=stats.shared_proportion_min,
                shared_proportion_max=stats.shared_proportion_max,
            )

        return node_stats_dict

    def prep_tree_for_clade_analysis(self) -> PhyloTreeManager:
        njtree = self.generate_tree()
        ### inner node to leaf dict
        tree_manager = PhyloTreeManager(njtree)
        # inner_node_leaf_dict = tree_manager.clades_get_leaves_clades()
        return tree_manager

    @property
    def all_clade_leaves_filtered(self) -> Dict[Phylo.BaseTree.Clade, list]:
        all_node_leaves = self.tree_manager.all_clades_leaves()

        all_node_leaves = {
            node: [leaf for leaf in leaves if leaf not in self.excluded_leaves]
            for node, leaves in all_node_leaves.items()
        }

        return all_node_leaves

    def get_node_statistics(self, force=False) -> Dict[Phylo.BaseTree.Clade, Clade]:
        if os.path.isfile(self.clade_statistics_path) and not force:
            clade_summary = pd.read_csv(self.clade_statistics_path)
            node_statistics_dict = self.clade_dict_from_summary(
                clade_summary=clade_summary,
                tree=self.tree_manager.tree,
                inner_node_leaf_dict=self.all_clade_leaves_filtered,
            )

        else:
            node_statistics_dict = self.node_statistics()

            clade_summary = self.all_clades_summary(
                node_stats_dict=node_statistics_dict
            )

            clade_summary.to_csv(self.clade_statistics_path, index=False)

        return node_statistics_dict

    def node_private_reads(self, inner_node_leaf_dict: dict) -> dict:
        """
        Calculate statistics per internal node to use for private clade analysis.
        """
        node_private_dict = {}
        for node, leaves in inner_node_leaf_dict.items():
            group = self.read_profile_matrix.loc[leaves]
            group_sum = group.sum(axis=0)
            group_sum_as_bool = group_sum > 0
            group_sum_as_bool_list = group_sum_as_bool.tolist()

            group_reads = self.read_profile_matrix.iloc[:, group_sum_as_bool_list]
            group_reads_sum_all = group_reads.sum(axis=0)
            group_reads_sum_group = group_reads.loc[leaves].sum(axis=0)
            group_reads_sum_group = group_reads_sum_group.fillna(0)

            read_proportions = group_reads_sum_group / group_reads_sum_all
            read_proportions = read_proportions.fillna(0)

            # proportion of reads in group that are private
            proportion_private = read_proportions.sum() / len(read_proportions)
            node_private_dict[node] = proportion_private

        return node_private_dict

    def filter_clades(self, clades_dict):
        """
        Return list of clades with private reads above threshold"""

        print("clades_dict")
        print(self.clade_filter.reference_clade)

        for clade, clade_obj in clades_dict.items():
            print("###########")
            print(clade)
            print(clade_obj)
            print(self.clade_filter.filter_clade(clade_obj))

        clades_filtered = [
            clade
            for clade, clade_obj in clades_dict.items()
            if self.clade_filter.filter_clade(clade_obj)
            and self.safe_clade_name(clade) != "None"
        ]
        return clades_filtered

    @staticmethod
    def safe_clade_name(clade: Phylo.BaseTree.Clade) -> str:
        """
        Return clade name if exists, otherwise return None
        """
        try:
            return clade.name
        except:
            return "None"

    def leaf_clades_to_pandas(
        self,
        leaf_clades: Dict[str, Phylo.BaseTree.Clade],
        statistics_dict: Dict[Phylo.BaseTree.Clade, Clade],
    ) -> pd.DataFrame:
        """
        Return dataframe of leaf clades
        """
        leaf_clades_dict = []
        for leaf, clade in leaf_clades.items():
            if clade is None:
                leaf_clades_dict.append(
                    (
                        leaf,
                        "None",
                        0,
                        0,
                        0,
                    )
                )
            else:
                try:
                    leaf_clades_dict.append(
                        (
                            leaf,
                            self.safe_clade_name(clade),
                            statistics_dict[clade].group_counts,
                            statistics_dict[clade].private_proportion,
                            statistics_dict[clade].shared_proportion_max,
                        )
                    )
                except KeyError:
                    # self.logger.error("clade not found")
                    leaf_clades_dict.append(
                        (
                            leaf,
                            "None",
                            0,
                            0,
                            0,
                        )
                    )

        ##

        leaf_clades_df = pd.DataFrame(
            leaf_clades_dict,
            columns=[
                "leaf",
                "clade",
                "total_counts",
                "private_proportion",
                "shared_proportion",
            ],
        )

        ##

        accids_df = pd.read_csv(self.accid_statistics_path, sep="\t")
        #
        leaf_clades_df = leaf_clades_df.merge(
            accids_df, right_on="accid", left_on="leaf"
        )

        def copy_leaf_to_clade_if_none(row):
            if row.clade == "None":
                row.clade = row.leaf
            return row

        leaf_clades_df = leaf_clades_df.apply(copy_leaf_to_clade_if_none, axis=1)

        def set_single_count_to_zero(row):
            if row.clade == "single":
                row.group_count = 0
            return row

        def set_No_clade_to_leaf(row):
            if row.clade == "None":
                row.clade = row.leaf
            return row

        leaf_clades_df.apply(set_single_count_to_zero, axis=1)
        leaf_clades_df.apply(set_No_clade_to_leaf, axis=1)

        leaf_clades_df.sort_values(
            by=["total_counts", "clade", "read_count"],
            ascending=[False, True, False],
            inplace=True,
        )

        leaf_clades_df.reset_index(drop=True, inplace=True)
        return leaf_clades_df

    def get_leaf_clades(self, force=False) -> pd.DataFrame:
        statistics_dict_all = self.get_node_statistics(force=force)

        selected_clades = self.filter_clades(statistics_dict_all)

        print("#### selected_clades")
        print(selected_clades)
        leaf_clades = self.tree_manager.leaf_clades_clean(selected_clades)

        clades = self.leaf_clades_to_pandas(leaf_clades, statistics_dict_all)

        return clades
