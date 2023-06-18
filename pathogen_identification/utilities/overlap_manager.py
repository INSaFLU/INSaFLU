from typing import Tuple
import itertools as it
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from scipy.spatial.distance import pdist, squareform
import pandas as pd
import os

## pairwise matrix by individual reads
from pathogen_identification.utilities.utilities_general import readname_from_fasta

from typing import List, Dict

from dataclasses import dataclass
import itertools
import numpy as np
from abc import ABC, abstractmethod
from pathogen_identification.utilities.clade_objects import Clade, CladeFilter


def accid_from_metadata(metadata: pd.DataFrame, read_name: str) -> str:
    """
    Return accession id of read from metadata
    """
    filename = os.path.basename(read_name)
    try:
        return metadata[metadata["filename"] == filename]["accid"].values[0]
    except IndexError:
        return read_name


class ReadOverlapManager:
    def __init__(
        self,
        fasta_list: List[str],
        metadata_df: pd.DataFrame,
        reference_clade: Clade,
    ):
        self.metadata = metadata_df
        self.fasta_list = fasta_list
        self.clade_filter = CladeFilter(reference_clade=reference_clade)

        self.metadata["filename"] = self.metadata["file"].apply(
            lambda x: x.split("/")[-1]
        )

        self.metadata["filepath"] = self.metadata["file"]

        self.read_profile_matrix: pd.DataFrame = self.generate_read_matrix()
        self.overlap_matrix: pd.DataFrame = self.readoverlap_allpairs_df()

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
        return [1 if read in readname_dict[accid] else 0 for read in all_reads]

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
        files_lists = self.metadata.file.values
        read_lists = [readname_from_fasta(f) for f in files_lists]

        read_overlap_dict = self.readoverlap_allpairs(read_lists)

        overlap_matrix = pd.DataFrame(
            data=[[x[0], x[1], g] for x, g in read_overlap_dict.items()],
            columns=["A", "B", "overlap"],
        )

        overlap_matrix["accid_A"] = overlap_matrix["A"].apply(
            lambda x: self.metadata["accid"].values[x]
        )
        overlap_matrix["accid_B"] = overlap_matrix["B"].apply(
            lambda x: self.metadata["accid"].values[x]
        )

        return overlap_matrix

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

    def generate_read_matrix(self):
        """
        Generate read matrix
        """
        files_readnames = [
            readname_from_fasta(fasta_file) for fasta_file in self.fasta_list
        ]
        readname_dict = self.get_accid_readname_dict()
        all_reads = self.all_reads_set(files_readnames)
        read_profile_dict = self.read_profile_dict_get(readname_dict, all_reads)
        read_profile_matrix = self.read_profile_matrix_get(read_profile_dict)
        return read_profile_matrix

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

    def generate_distance_matrix(self):
        """
        Generate distance matrix
        """
        distance_matrix = self.pairwise_shared_reads(self.read_profile_matrix)
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

    def clade_shared_by_pair(self, leaves: list) -> pd.DataFrame:
        """
        return tuple of proportions of reads shared by each pair of leaves
        """

        subset_clade = self.overlap_matrix.loc[
            (self.overlap_matrix.accid_A.isin(leaves))
            & (self.overlap_matrix.accid_B.isin(leaves))
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

    def clade_private_proportions(self, leaves: list) -> float:
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

    def node_statistics(self, inner_node_leaf_dict: dict) -> dict:
        node_stats_dict = {}
        for node, leaves in inner_node_leaf_dict.items():
            if len(leaves) == 0:
                node_stats_dict[node] = Clade(
                    name=node,
                    leaves=leaves,
                    private_proportion=0,
                    shared_proportion_std=0,
                    shared_proportion_min=0,
                    shared_proportion_max=0,
                )

                continue

            proportion_private = self.clade_private_proportions(leaves)

            if len(leaves) == 1:
                node_stats_dict[node] = Clade(
                    name=node,
                    leaves=leaves,
                    private_proportion=proportion_private,
                    shared_proportion_std=1,
                    shared_proportion_min=0,
                    shared_proportion_max=0,
                )

                continue

            pairwise_shared_clade = self.clade_shared_by_pair(leaves)

            node_stats_dict[node] = Clade(
                name=node,
                leaves=leaves,
                private_proportion=proportion_private,
                shared_proportion_min=min(pairwise_shared_clade.proportion_max),
                shared_proportion_max=max(pairwise_shared_clade.proportion_max),
                shared_proportion_std=np.std(pairwise_shared_clade.proportion_std),
            )

        return node_stats_dict

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
        clades_filtered = [
            clade
            for clade, clade_obj in clades_dict.items()
            if self.clade_filter.filter_clade(clade_obj)
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
        self, leaf_clades: Dict[str, Phylo.BaseTree.Clade]
    ) -> pd.DataFrame:
        """
        Return dataframe of leaf clades
        """

        leaf_clades_dict = [
            (leaf, self.safe_clade_name(clade)) for leaf, clade in leaf_clades.items()
        ]
        leaf_clades_df = pd.DataFrame(leaf_clades_dict, columns=["leaf", "clade"])

        leaf_clades_df["read_count"] = leaf_clades_df["leaf"].apply(
            lambda x: self.get_accession_total_counts(x)
        )
        leaf_clades_df["proportion"] = leaf_clades_df["leaf"].apply(
            lambda x: self.get_proportion_counts(x)
        )
        # sort by clade and then read count

        def group_count(clade):
            return leaf_clades_df[leaf_clades_df["clade"] == clade].read_count.sum()

        leaf_clades_df["group_count"] = leaf_clades_df["clade"].apply(group_count)

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
            by=["group_count", "clade", "read_count"],
            ascending=[False, True, False],
            inplace=True,
        )

        leaf_clades_df.reset_index(drop=True, inplace=True)
        return leaf_clades_df
