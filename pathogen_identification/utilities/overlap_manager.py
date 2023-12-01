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
    overlap_pca_plot_filename: str = "overlap_pca_{}.png"
    min_freq: float = 0
    max_reads: int = 500000

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
        self.pid = pid

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
        self.overlap_pca_plot_path = os.path.join(
            self.media_dir, self.overlap_pca_plot_filename.format(pid)
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

    def get_media_path_heatmap_clade(self, clade_str: str):
        clade_str_keep = clade_str.replace(" ", "_").lower()
        return os.path.join(
            self.media_dir, f"heatmap_clade_{clade_str_keep}.{self.pid}.png"
        )

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
        accid_df = self.metadata[["accid", "description"]]

        self.read_profile_matrix: pd.DataFrame = self.generate_read_matrix()
        self.total_read_counts = self.read_profile_matrix.sum(axis=0)
        self.read_profile_matrix_filtered: pd.DataFrame = self.filter_read_matrix(
            self.read_profile_matrix
        )

        print(self.read_profile_matrix.index)

        self.overlap_matrix: pd.DataFrame = self.pairwise_shared_count(
            self.read_profile_matrix_filtered
        )

        accid_df["read_count"] = accid_df.accid.apply(self.get_accession_total_counts)

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
        # read_profile_matrix_filtered = read_profile_matrix.loc[:, read_counts > 1]
        # filter reads with less than min_freq
        read_profile_matrix_filtered = read_profile_matrix.loc[
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
        ## create list of duplicated rows
        return read_profile_matrix

    def get_private_reads_no_duplicates(
        self,
    ):
        """ """

        accid_df = pd.read_csv(self.accid_statistics_path, sep="\t")
        # duplicate_groups = self.duplicate_groups_from_dataframe(
        #    self.read_profile_matrix
        # )
        similar_groups = self.very_similar_groups_from_dataframe()

        if "private_reads" not in accid_df.columns:
            accid_df["private_reads"] = 0
        print(similar_groups)
        print("duplicate groups")
        for duplicate_group in similar_groups:
            print("#")
            print(duplicate_group)
            # group_private_counts = self.get_accession_private_counts(duplicate_group)
            (
                group_private_counts,
                total_reads,
                proportion_private,
            ) = self.clade_private_proportions(list(duplicate_group))
            print(group_private_counts)
            print(accid_df.loc[accid_df.accid.isin(duplicate_group), "private_reads"])
            accid_df.loc[
                accid_df.accid.isin(duplicate_group), "private_reads"
            ] = group_private_counts
            print(accid_df.loc[accid_df.accid.isin(duplicate_group), "private_reads"])

        accid_df.to_csv(self.accid_statistics_path, sep="\t", index=False)

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
        return self.read_profile_matrix_filtered.loc[accid].sum()

    def get_proportion_counts(self, accid: str):
        """
        Get proportion counts for accession
        """

        return (
            self.get_accession_total_counts(accid)
            / self.read_profile_matrix_filtered.sum().sum()
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
            distance_matrix = self.pairwise_shared_reads(
                self.read_profile_matrix_filtered
            )

        if not self.check_all_accessions_in_distance_matrix(distance_matrix):
            self.parse_for_data()
            distance_matrix = self.pairwise_shared_reads(
                self.read_profile_matrix_filtered
            )

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
        group = self.read_profile_matrix_filtered.loc[leaves]

        group_pairwise_shared = self.pairwise_shared_count(group)

        # divide shared rows by group row sums
        group_pairwise_shared = group_pairwise_shared.div(group.sum(axis=1), axis=1)

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
        return int(self.read_profile_matrix_filtered.loc[leaves].sum().sum())

    def clade_private_proportions_old(self, leaves: list) -> float:
        """ """
        group = self.read_profile_matrix_filtered.loc[leaves]
        group_sum = group.sum(axis=0)
        group_sum_as_bool = group_sum > 0
        group_sum_as_bool_list = group_sum_as_bool.tolist()

        group_reads = self.read_profile_matrix_filtered.iloc[:, group_sum_as_bool_list]

        group_reads_sum_all = group_reads.sum(axis=0)
        group_reads_sum_group = group_reads.loc[leaves].sum(axis=0)
        group_reads_sum_group = group_reads_sum_group.fillna(0)

        read_proportions = group_reads_sum_group / group_reads_sum_all
        read_proportions = read_proportions.fillna(0)

        # proportion of reads in group that are private
        proportion_private = read_proportions.sum() / len(read_proportions)

        return proportion_private

    def clade_private_proportions(self, leaves: list) -> Tuple[float, float, float]:
        """ """
        group = self.read_profile_matrix.loc[leaves]
        group_sum = group.sum(axis=0)
        group_sum_as_bool = group_sum > 0
        group_sum_as_bool_list = group_sum_as_bool.tolist()

        sum_all = self.read_profile_matrix_filtered.iloc[:, group_sum_as_bool_list].sum(
            axis=0
        )
        sum_group = self.read_profile_matrix.loc[leaves]
        sum_group = sum_group.iloc[:, group_sum_as_bool_list].sum(axis=0)

        private_reads = sum_group == sum_all

        private_reads = sum(private_reads)
        print(private_reads)
        print(len(group_sum_as_bool_list))
        print(sum(group_sum_as_bool_list))

        if sum(group_sum_as_bool_list) == 0:
            private_reads = 0
            total_reads = 0
            proportion_private = 0
        else:
            total_reads = sum(group_sum_as_bool_list)
            proportion_private = private_reads / total_reads

        return private_reads, total_reads, proportion_private

    def get_accession_private_counts(self, duplicate_group: tuple) -> int:
        """
        Get private counts for accession
        """

        duplicate_group_counts = self.read_profile_matrix.loc[
            self.read_profile_matrix.index.isin(duplicate_group) == True
        ].sum(axis=0)

        duplicate_counts_as_bool = duplicate_group_counts > 0

        private_counts = duplicate_group_counts == self.total_read_counts
        private_counts = private_counts[duplicate_counts_as_bool]

        private_counts = sum(private_counts)

        return private_counts

    def between_clade_reads_matrix(
        self, filter_names=[], remove_leaves=True, sort_private=False
    ) -> pd.DataFrame:
        """
        Return dataframe reads per clade"""
        clade_read_matrix = []
        belonging = []
        private_sort = {}
        for clade, leaves in self.all_clade_leaves_filtered.items():
            if len(leaves) == 0:
                continue

            # continue clade is leaf:
            if len(leaves) == 1 and remove_leaves:
                continue

            if filter_names:
                if clade.name not in filter_names:
                    continue

            reads_in_clade = self.read_profile_matrix_filtered.loc[leaves]
            reads_in_clade_sum = reads_in_clade.sum(axis=0)
            reads_in_clade_sum_as_bool = reads_in_clade_sum > 0
            reads_in_clade_sum_as_int_list = reads_in_clade_sum_as_bool.astype(
                int
            ).tolist()
            clade_read_matrix.append(reads_in_clade_sum_as_int_list)
            belonging.append(clade.name)
            if sort_private:
                (
                    private_reads,
                    total_reads,
                    proportion_private,
                ) = self.clade_private_proportions(leaves)

                private_sort[clade.name] = proportion_private

        clade_read_matrix = pd.DataFrame(
            clade_read_matrix,
            index=belonging,
            columns=self.read_profile_matrix_filtered.columns,
        )

        if sort_private:
            private_sort = {
                k: v
                for k, v in sorted(
                    private_sort.items(), key=lambda item: item[1], reverse=True
                )
            }
            clade_read_matrix = clade_read_matrix.reindex(private_sort.keys())

        return clade_read_matrix

    def within_clade_reads_matrix(self, leaves: list) -> pd.DataFrame:
        """
        Return dataframe reads per clade"""
        reads_in_clade = self.read_profile_matrix_filtered.loc[leaves]

        return reads_in_clade

    def square_and_fill_diagonal(self, clade_read_matrix: pd.DataFrame) -> pd.DataFrame:
        shared_clade_matrix = self.pairwise_shared_count(clade_read_matrix)

        ## divide rows of shared_clade_matrix by clade_read_matrix row sums

        shared_clade_matrix = shared_clade_matrix.div(
            clade_read_matrix.sum(axis=1),
            axis=0,
        )

        # set diagonal to 0
        np.fill_diagonal(shared_clade_matrix.values, 0)

        # fill na with 0
        shared_clade_matrix = shared_clade_matrix.fillna(0)

        return shared_clade_matrix

    def between_clade_shared_reads(self, clades_filter=[]) -> pd.DataFrame:
        """
        Return dataframe of pairwise shared reads between all pairs of clades"""
        clade_read_matrix = self.between_clade_reads_matrix(
            filter_names=clades_filter, remove_leaves=False
        )
        shared_clade_matrix = self.square_and_fill_diagonal(clade_read_matrix)

        return shared_clade_matrix

    def within_clade_shared_reads(self, clade) -> pd.DataFrame:
        """ """
        clade_node = [
            x for x in self.all_clade_leaves_filtered.keys() if x.name == clade
        ][0]
        leaves = self.all_clade_leaves_filtered[clade_node]
        print("#########################")
        print(clade_node)
        print(leaves)
        clade_read_matrix = self.within_clade_reads_matrix(leaves)
        print(clade_read_matrix.sum(axis=1))
        shared_clade_matrix = self.square_and_fill_diagonal(clade_read_matrix)
        return shared_clade_matrix

    def duplicate_groups_from_dataframe(self, pairwise_shared_clade: pd.DataFrame):
        """
        find duplicate entries in pairwise shared clade
        """
        duplicate_groups = (
            pairwise_shared_clade.groupby(pairwise_shared_clade.columns.tolist())
            .apply(lambda x: tuple(x.index))
            .tolist()
        )
        # duplicate_groups = [x for x in duplicate_groups if len(x) > 1]
        return duplicate_groups

    def very_similar_groups_from_dataframe(self) -> List[tuple]:
        """
        find very similar entries entries in pairwise shared clade
        """

        shared_read_matrix = self.square_and_fill_diagonal(
            self.read_profile_matrix_filtered
        )
        threshold = 0.95

        clusters_assigment_dict = {}
        clusternum = 0

        for i in range(shared_read_matrix.shape[0]):
            for j in range(shared_read_matrix.shape[1]):
                if i == j:
                    continue
                if (
                    shared_read_matrix.iloc[i, j] >= threshold
                    and shared_read_matrix.iloc[j, i] >= threshold
                ):
                    assingments = [
                        clusters_assigment_dict.get(i, None),
                        clusters_assigment_dict.get(j, None),
                    ]
                    assingments = [x for x in assingments if x is not None]
                    if len(assingments) == 0:
                        clusters_assigment_dict[i] = clusternum
                        clusters_assigment_dict[j] = clusternum
                        clusternum += 1
                    elif len(assingments) == 1:
                        clusters_assigment_dict[i] = assingments[0]
                        clusters_assigment_dict[j] = assingments[0]

        clusters = {}
        for k, v in clusters_assigment_dict.items():
            clusters.setdefault(v, []).append(k)

        clusters = [v for k, v in clusters.items()]
        print(clusters)
        for i in range(shared_read_matrix.shape[0]):
            if i not in clusters_assigment_dict.keys():
                clusters.append([i])

        clusters = [x for x in clusters if len(x) > 1]
        clusters = [
            tuple([self.read_profile_matrix_filtered.index[y] for y in x])
            for x in clusters
        ]

        return clusters

    def plot_pairwise_shared_clade_reads(
        self, pairwise_shared_clade: pd.DataFrame, subplot=False, clade_str=""
    ) -> None:
        """
        Plot heatmap of pairwise shared reads between all pairs of leaves
        """

        plt.figure(figsize=(15, 6))
        sns.heatmap(pairwise_shared_clade, annot=True)
        # center figure
        plt.tight_layout()
        if subplot is False:
            plot_path = self.overlap_matrix_plot_path
        else:
            plot_path = self.get_media_path_heatmap_clade(clade_str)
        plt.savefig(plot_path, bbox_inches="tight")

    def plot_pca_full(self, accid_df: pd.DataFrame) -> None:
        """
        plot pca of all hits using reads as features
        for colors, assign colors to hits in clades with private proportion > 0.5.
        draw two plots: 1/2 pc and 2/3 pc
        """

        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler

        pca = PCA(n_components=3)
        pca.fit(self.read_profile_matrix_filtered)

        pca_df = pd.DataFrame(
            pca.transform(self.read_profile_matrix_filtered),
            index=self.read_profile_matrix_filtered.index,
            columns=["pc1", "pc2", "pc3"],
        )

        pca_df["accid"] = pca_df.index
        pca_df["clade"] = "None"

        def find_clade(string):
            if string in accid_df.accid.tolist():
                return accid_df.loc[accid_df.accid == string, "clade"].values[0]

            return "None"

        pca_df["clade"] = pca_df["accid"].apply(find_clade)
        pca_df["clade"] = pca_df["clade"].astype("category")

        fig, ax = plt.subplots(1, 2, figsize=(15, 6))
        # scatterplot with large points
        sns.scatterplot(
            data=pca_df,
            x="pc1",
            y="pc2",
            hue="clade",
            s=100,
            ax=ax[0],
        )
        sns.scatterplot(
            data=pca_df,
            x="pc2",
            y="pc3",
            hue="clade",
            s=100,
            ax=ax[1],
        )
        plt.tight_layout()
        plt.savefig(self.overlap_pca_plot_path)

    def node_statistics(self) -> dict:
        self.parse_for_data()
        self.update_excluded_leaves(self.read_profile_matrix_filtered)

        node_stats_dict = {}

        if self.read_profile_matrix_filtered.shape[1] == 0:
            self.logger.info("No reads with frequency > min_freq")
            return node_stats_dict

        for node, leaves in self.all_clade_leaves_filtered.items():
            if len(leaves) == 0:
                node_stats_dict[node] = Clade(
                    name=node,
                    leaves=leaves,
                    group_counts=0,
                    private_proportion=0,
                    total_proportion=0,
                    shared_proportion_std=0,
                    shared_proportion_min=0,
                    shared_proportion_max=0,
                    overlap_df=pd.DataFrame(),
                )

                continue

            (
                private_reads,
                total_reads,
                proportion_private,
            ) = self.clade_private_proportions(leaves)
            total_proportion = total_reads / self.read_profile_matrix_filtered.shape[1]

            clade_counts = self.clade_total_counts(leaves)

            if len(leaves) == 1:
                node_stats_dict[node] = Clade(
                    name=node,
                    leaves=leaves,
                    group_counts=clade_counts,
                    private_proportion=proportion_private,
                    total_proportion=total_proportion,
                    shared_proportion_std=0,
                    shared_proportion_min=0,
                    shared_proportion_max=0,
                    overlap_df=pd.DataFrame(),
                )

                continue

            combinations = self.clade_shared_by_pair_old(leaves)

            node_stats_dict[node] = Clade(
                name=node,
                leaves=leaves,
                private_proportion=proportion_private,
                total_proportion=total_proportion,
                group_counts=clade_counts,
                shared_proportion_min=min(combinations.proportion_max),
                shared_proportion_max=max(combinations.proportion_max),
                shared_proportion_std=np.std(combinations.proportion_max),
                overlap_df=combinations,
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
                    stats.total_proportion,
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
                "total_proportion",
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
                total_proportion=stats.private_proportion,
                group_counts=stats.total_counts,
                shared_proportion_std=stats.shared_proportion_std,
                shared_proportion_min=stats.shared_proportion_min,
                shared_proportion_max=stats.shared_proportion_max,
                overlap_df=pd.DataFrame(),
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
            group = self.read_profile_matrix_filtered.loc[leaves]
            group_sum = group.sum(axis=0)
            group_sum_as_bool = group_sum > 0
            group_sum_as_bool_list = group_sum_as_bool.tolist()

            group_reads = self.read_profile_matrix_filtered.iloc[
                :, group_sum_as_bool_list
            ]
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

        leaf_clades = self.tree_manager.leaf_clades_clean(selected_clades)

        clades = self.leaf_clades_to_pandas(leaf_clades, statistics_dict_all)

        return clades
