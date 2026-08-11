import itertools as it
import logging
import os
import re as _re
import shutil
import subprocess
import tempfile
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from scipy.spatial.distance import pdist, squareform
from scipy.stats import entropy as scipy_entropy

from pathogen_identification.utilities.clade_objects import Clade, CladeFilter
from pathogen_identification.utilities.phylo_tree import PhyloTreeManager
## pairwise matrix by individual reads

from collections import Counter

import numpy as np


def shannon_diversity(proportions: list[float]) -> float:
    """
    Calculate Shannon diversity index given a list of proportions.
    """
    return -sum(p * (p if p == 0 else np.log(p)) for p in proportions)


def shannon_diversity_from_counts(counts: list[int]) -> float:
    """
    Calculate Shannon diversity index given a list of counts.
    """
    total = sum(counts)
    if total == 0:
        return 0.0
    proportions = [count / total for count in counts]
    return shannon_diversity(proportions)

def shannon_diversity_from_list(taxa: list[str]) -> float:
    """
    Calculate Shannon diversity index given a list of taxa.
    """

    if not taxa:
        return 0.0
    counts = Counter(taxa)
    return shannon_diversity_from_counts(list(counts.values()))

def pairwise_shared_count(
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


def square_and_fill_diagonal(clade_read_matrix: pd.DataFrame) -> pd.DataFrame:
    shared_clade_matrix = pairwise_shared_count(clade_read_matrix)

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


def very_similar_groups_from_dataframe(
    read_profile_matrix_filtered: pd.DataFrame, threshold=0.95
) -> List[tuple]:
    """
    find very similar entries entries in pairwise shared clade
    """

    shared_read_matrix = square_and_fill_diagonal(read_profile_matrix_filtered)

    clusters_assigment_dict = {}
    clusternum = 0

    for i in range(shared_read_matrix.shape[0]):
        for j in range(shared_read_matrix.shape[1]):
            if i == j:
                continue
            if (
                shared_read_matrix.iloc[i, j] >= threshold
                or shared_read_matrix.iloc[j, i] >= threshold
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
    for i in range(shared_read_matrix.shape[0]):
        if i not in clusters_assigment_dict.keys():
            clusters.append([i])

    # clusters = [x for x in clusters if len(x) > 1]
    clusters = [
        tuple([read_profile_matrix_filtered.index[y] for y in x]) for x in clusters
    ]

    return clusters


def clade_private_proportions_old(read_profile_matrix_filtered, leaves: list) -> float:
    """ """
    group = read_profile_matrix_filtered.loc[leaves]
    group_sum = group.sum(axis=0)
    group_sum_as_bool = group_sum > 0
    group_sum_as_bool_list = group_sum_as_bool.tolist()

    group_reads = read_profile_matrix_filtered.iloc[:, group_sum_as_bool_list]

    group_reads_sum_all = group_reads.sum(axis=0)
    group_reads_sum_group = group_reads.loc[leaves].sum(axis=0)
    group_reads_sum_group = group_reads_sum_group.fillna(0)

    read_proportions = group_reads_sum_group / group_reads_sum_all
    read_proportions = read_proportions.fillna(0)

    # proportion of reads in group that are private
    proportion_private = read_proportions.sum() / len(read_proportions)

    return proportion_private


def clade_private_proportions(
    read_profile_matrix: pd.DataFrame, leaves: list
) -> Tuple[float, float, float]:
    """ """
    group = read_profile_matrix.loc[leaves]
    group_sum = group.sum(axis=0)
    group_sum_as_bool = group_sum > 0
    group_sum_as_bool_list = group_sum_as_bool.tolist()

    sum_all = read_profile_matrix.iloc[:, group_sum_as_bool_list].sum(axis=0)
    sum_group = read_profile_matrix.loc[leaves]
    sum_group = sum_group.iloc[:, group_sum_as_bool_list].sum(axis=0)

    private_reads = sum_group == sum_all

    private_reads = sum(private_reads)

    if sum(group_sum_as_bool_list) == 0:
        private_reads = 0
        total_reads = 0
        proportion_private = 0
    else:
        total_reads = sum(group_sum_as_bool_list)
        proportion_private = private_reads / total_reads

    return private_reads, total_reads, proportion_private


def pairwise_shared_reads_distance(read_profile_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Return dataframe of pairwise shared reads
    this is a symmetrical matrix with the max proportion of shared reads as their mutual distance.
    """
    #########
    pairwise_props = square_and_fill_diagonal(read_profile_matrix)

    for i in range(pairwise_props.shape[0]):
        pairwise_props.iloc[i, i] = 1
        for j in range(i + 1, pairwise_props.shape[1]):
            shared_ij = pairwise_props.iloc[i, j]
            shared_ji = pairwise_props.iloc[j, i]
            shared_pair = (shared_ij, shared_ji)
            pairwise_props.iloc[i, j] = max(shared_pair)
            pairwise_props.iloc[j, i] = max(shared_pair)

    pairwise_props = 1 - pairwise_props

    return pairwise_props


def pairwise_shared_reads(read_profile_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Return dataframe of pairwise shared reads
    this is a symmetrical matrix with the max proportion of shared reads as their mutual distance.
    """

    pairwise_props = square_and_fill_diagonal(read_profile_matrix)

    return pairwise_props


def pairwise_shared_reads_old(read_profile_matrix: pd.DataFrame) -> pd.DataFrame:
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


class MappingResultsParser:

    read_profile_matrix_filtered: pd.DataFrame
    overlap_matrix: pd.DataFrame
    accid_statistics_filename: str = "accid_statistics_{}.tsv"
    min_freq: float = 0
    max_reads: int

    def __init__(
        self,
        metadata_df: pd.DataFrame,
        media_dir: str,
        pid: str,
        max_reads: int = 500000,
    ):
        self.metadata = metadata_df
        self.media_dir = media_dir
        self.pid = pid
        self.max_reads = max_reads

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    @property
    def accid_statistics_path(self):
        return os.path.join(
            self.media_dir, self.accid_statistics_filename.format(self.pid)
        )

    @staticmethod
    def accid_from_metadata(metadata: pd.DataFrame, read_name: str) -> str:
        """
        Return accession id of read from metadata
        """
        filename = os.path.basename(read_name)
        try:
            return metadata[metadata["filename"] == filename]["accid"].values[0]
        except IndexError:
            return read_name

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

    def prep_accid_table(self):
        #########
        accid_df = self.metadata[["accid", "description"]]

        accid_df["read_count"] = accid_df.accid.apply(self.get_accession_total_counts)

        # sort table by accid and then by read count
        accid_df = accid_df.sort_values(["accid", "read_count"], ascending=False)

        # drop duplicates of accid
        accid_df = accid_df.drop_duplicates(subset=["accid"])

        # get proportion of reads
        accid_df["proportion"] = accid_df["accid"].apply(
            lambda x: self.get_proportion_counts(x)
        )
        accid_df.to_csv(self.accid_statistics_path, sep="\t", index=False)

        return accid_df

    def filter_read_matrix(self, read_profile_matrix: pd.DataFrame) -> pd.DataFrame:
        """
        Filter read matrix, reads as columns, accids as rows
        """

        if read_profile_matrix.shape[1] == 0:
            self.logger.info("No reads with frequency > min_freq")
            return read_profile_matrix

        read_counts = read_profile_matrix.sum(axis=0)
        read_freqs = read_counts / read_profile_matrix.shape[0]
        # filter reads by frequency
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

class ReadOverlapManager(MappingResultsParser):
    distance_matrix_filename: str = "distance_matrix_{}.tsv"
    shared_prop_matrix_filename: str = "shared_prop_matrix_{}.tsv"
    clade_shared_prop_matrix_filename: str = "clade_shared_prop_matrix_{}.tsv"
    clade_statistics_filename: str = "clade_statistics_{}.tsv"

    tree_plot_filename: str = "tree_{}.png"
    overlap_matrix_plot_filename: str = "overlap_matrix_{}.png"
    overlap_pca_plot_filename: str = "overlap_pca_{}.png"

    def __init__(
        self,
        metadata_df: pd.DataFrame,
        reference_clade: Clade,
        media_dir: str,
        pid: str,
        force_tree_rebuild: bool = False,
        max_reads: int = 500000,
        clustering_model_type: str = "Fixed",
    ):

        super().__init__(metadata_df, media_dir, pid, max_reads=max_reads)

        self.clade_filter = CladeFilter(reference_clade=reference_clade)
        self.excluded_leaves = []
        self.force_tree_rebuild = force_tree_rebuild
        self.clustering_model_type = clustering_model_type
        self._accid_statistics_df = None
        self.distance_mat = pd.DataFrame()

        self.metadata["filename"] = self.metadata["file"].apply(
            lambda x: x.split("/")[-1]
        )

        self.metadata["filepath"] = self.metadata["file"]

        try:
            self.tree_manager = self.prep_tree_for_clade_analysis()
        except Exception as e:
            print(e)

        self.tree_plot_exists = os.path.exists(self.tree_plot_path)
        self.tree_plot_path_render = os.path.join(
            "/media/", self.tree_plot_path.split("/media/")[-1]
        )
        self.overlap_matrix_plot_exists = os.path.exists(self.overlap_matrix_plot_path)
        self.analysis_exists = os.path.exists(self.media_dir)

    @property
    def distance_matrix_path(self):
        return os.path.join(
            self.media_dir, self.distance_matrix_filename.format(self.pid)
        )

    @property
    def shared_prop_matrix_path(self):
        return os.path.join(
            self.media_dir, self.shared_prop_matrix_filename.format(self.pid)
        )

    @property
    def clade_shared_prop_matrix_path(self):
        return os.path.join(
            self.media_dir, self.clade_shared_prop_matrix_filename.format(self.pid)
        )

    @property
    def clade_statistics_path(self):
        return os.path.join(
            self.media_dir, self.clade_statistics_filename.format(self.pid)
        )

    @property
    def tree_plot_path(self):
        return os.path.join(self.media_dir, self.tree_plot_filename.format(self.pid))

    @property
    def overlap_matrix_plot_path(self):
        return os.path.join(
            self.media_dir, self.overlap_matrix_plot_filename.format(self.pid)
        )

    @property
    def overlap_pca_plot_path(self):
        return os.path.join(
            self.media_dir, self.overlap_pca_plot_filename.format(self.pid)
        )

    def build_tree(self):

        self.generate_tree()
        self.generate_shared_proportion_matrix()
        self.generate_clade_shared_proportion_matrix()        

        if not os.path.exists(self.tree_plot_path):
            try:
                self.tree_manager.plot_tree(self.tree_plot_path)
            except Exception as e:
                print(e)

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

    @staticmethod
    def readoverlap_2_files(lista, listb) -> list:
        """
        Return list of read names that are in both lists
        """
        return list(set(lista).intersection(set(listb)))

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

    def generate_shared_proportion_matrix(self):
        """
        Generate shared proportion matrix
        """
        proportion_matrix = square_and_fill_diagonal(self.read_profile_matrix_filtered)
        proportion_matrix.to_csv(self.shared_prop_matrix_path)


    def generate_clade_shared_proportion_matrix(self):
        """
        Generate shared proportion matrix
        """
        clade_shared_proportion_matrix = self.between_clade_shared_reads()
        clade_shared_proportion_matrix.to_csv(self.clade_shared_prop_matrix_path)

    ####################
    ## Construct tree ###########################################################
    ####################

    def tree_from_distance_matrix(self, distance_matrix: pd.DataFrame):
        """
        Return tree from distance matrix
        """
        distmat = self.matrix_to_phylotriangle(distance_matrix)
        constructor = DistanceTreeConstructor()

        if distance_matrix.empty is True:
            return Phylo.BaseTree.Tree(rooted="True")

        if distance_matrix.shape[0] <= 1:
            tree = constructor.nj(distmat)
        else:
            tree = constructor.upgma(distmat)

        tree.rooted = False
        tree.ladderize()
        return tree

    def get_private_reads_no_duplicates(
        self,
    ):
        """ """

        accid_df = pd.read_csv(self.accid_statistics_path, sep="\t")

        similar_groups = very_similar_groups_from_dataframe(
            self.read_profile_matrix_filtered
        )

        if "private_reads" not in accid_df.columns:
            accid_df["private_reads"] = 0

        for duplicate_group in similar_groups:

            (
                group_private_counts,
                _,
                _,
            ) = clade_private_proportions(
                self.read_profile_matrix_filtered, list(duplicate_group)
            )

            accid_df.loc[accid_df.accid.isin(duplicate_group), "private_reads"] = (
                group_private_counts
            )

        accid_df.to_csv(self.accid_statistics_path, sep="\t", index=False)

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

    def check_all_accessions_in_distance_matrix(self, distance_matrix):
        for accid in self.metadata["accid"].values:
            if accid not in distance_matrix.columns:
                return False

        return True

    def _binary_path(self) -> str:
        from constants.constants import \
            Televir_Metadata_Constants as Deployment_Params
        return Deployment_Params().get_software_binary("mapping_to_matrix")

    def _symlink_bams_to_temp(self, temp_dir: str) -> Dict[str, str]:
        accid_to_src = {}
        
        for _, row in self.metadata.drop_duplicates('accid').iterrows():
            accid = row["accid"]
            bam = row["bam"]
            if not os.path.isfile(bam):
                self.logger.warning(f"BAM not found for {accid}: {bam}")
                continue
            link_path = os.path.join(temp_dir, f"{accid}.bam")
            os.symlink(bam, link_path)
            accid_to_src[accid] = bam
        return accid_to_src

    def _run_mapping_to_matrix(self, output_dir: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
        binary = self._binary_path()
        os.makedirs(output_dir, exist_ok=True)
        cmd = [
            binary,
            "--input-directory", output_dir,
            "--max-reads", str(self.max_reads),
            "--frequency-threshold", str(self.min_freq) if self.min_freq > 0 else "0.1",
            "--output-directory", output_dir,
            "--no-cluster-analysis",
        ]

        self.logger.info(f"Running mapping_to_matrix: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            self.logger.error(f"mapping_to_matrix failed:\n{result.stderr}")
            raise RuntimeError(f"mapping_to_matrix failed: {result.stderr}")
                
        distance_path = os.path.join(output_dir, "distance_matrix.tsv")
        presence_path = os.path.join(output_dir, "presence_absence_matrix.tsv")

        if not os.path.isfile(distance_path):
            distance_matrix = pd.DataFrame()
        else:
            distance_matrix = pd.read_csv(distance_path, sep="\t", index_col=0)

        # skip first line it is buffered and not a valid header
        if not os.path.isfile(presence_path):
            presence_matrix = pd.DataFrame()
        else:
            presence_matrix = pd.read_csv(presence_path, sep="\t", index_col=0, skiprows=1, header=None)

        if presence_matrix.shape[0] == 0 or presence_matrix.shape[1] == 0:
            presence_matrix = pd.DataFrame(
                index=distance_matrix.index, columns=distance_matrix.columns
            )
            return distance_matrix, presence_matrix
        
        # remove path and .bam extension from index
        distance_matrix.index = distance_matrix.index.astype(str).str.replace(r"\.bam$", "", regex=True)
        distance_matrix.index = distance_matrix.index.astype(str).str.replace(r".*/", "", regex=True)
        distance_matrix.columns = distance_matrix.columns.astype(str).str.replace(r"\.bam$", "", regex=True)
        distance_matrix.columns = distance_matrix.columns.astype(str).str.replace(r".*/", "", regex=True)

        presence_matrix.set_index(distance_matrix.index, inplace=True)

        return distance_matrix, presence_matrix

    def generate_distance_matrix(self, force=False):

        if os.path.isfile(self.distance_matrix_path) and not force:
            distance_matrix = pd.read_csv(self.distance_matrix_path, index_col=0)
        else:
            with tempfile.TemporaryDirectory() as tmpdir:
                self._symlink_bams_to_temp(tmpdir)
                distance_matrix, presence_matrix = self._run_mapping_to_matrix(tmpdir)
                self.read_profile_matrix_filtered = presence_matrix
                self.total_read_counts = presence_matrix.sum(axis=0)
                self.overlap_matrix = pairwise_shared_count(presence_matrix)
                
                self._accid_statistics_df = self.prep_accid_table()

        if not self.check_all_accessions_in_distance_matrix(distance_matrix):
            with tempfile.TemporaryDirectory() as tmpdir:
                self._symlink_bams_to_temp(tmpdir)
                distance_matrix, presence_matrix = self._run_mapping_to_matrix(tmpdir)
                self.read_profile_matrix_filtered = presence_matrix
                self.total_read_counts = presence_matrix.sum(axis=0)
                self.overlap_matrix = pairwise_shared_count(presence_matrix)
                self._accid_statistics_df = self.prep_accid_table()

        try:
            distance_matrix.to_csv(self.distance_matrix_path)
        except:
            pass

        return distance_matrix

    def  generate_tree(self):
        """
        This method is used to generate a tree structure using a distance matrix.

        The steps are as follows:
        1. Generate a distance matrix using the method 'generate_distance_matrix()'. This matrix represents the distance
           between different nodes.
        2. Use this distance matrix to generate the tree structure using the method 'tree_from_distance_matrix()'.
        Returns:
            tree: a tree structure generated from the distance matrix
        """
        # Generate the distance matrix
        self.distance_mat = self.generate_distance_matrix(force=True)

        # Generate the tree from the distance matrix
        tree = self.tree_from_distance_matrix(self.distance_mat)

        return tree
    

    def recalculate_all_min_pairwise_dist(self):
        if not os.path.exists(self.distance_matrix_path):
            return

        distance_matrix = self.generate_distance_matrix()
        proximity_matrix = 1 - distance_matrix

        wp_values = proximity_matrix.values
        index_arr = list(proximity_matrix.index)
        index_to_pos = {name: i for i, name in enumerate(index_arr)}

        def calc_node_stats(node):
            leaves_parted = self.tree_manager.get_leaves_parted(node)
            if len(leaves_parted) < 2:
                return {"Min_Pairwise_Dist": 0, "Min_Shared": 0, "Min_Dist": 0}

            left = leaves_parted[0]
            right = leaves_parted[1]
            
            left_pos = [index_to_pos[l.name] for l in left]
            right_pos = [index_to_pos[l.name] for l in right]
            all_pos = left_pos + right_pos

            sub = wp_values[np.ix_(left_pos, right_pos)]
            pair_dists = sub #np.maximum(sub, sub.T)
            min_dist = float(pair_dists.min()) if pair_dists.size > 0 else 0.0

            all_internal = set(left) | set(right)
            other_keys = [l for l in index_arr if l not in all_internal]
            if other_keys:
                other_pos = [index_to_pos[l] for l in other_keys]
                sub_ab = wp_values[np.ix_(all_pos, other_pos)]
                ext_kept = float(sub_ab.min())
            else:
                ext_kept = 0.0

            return {
                "Min_Dist": min_dist,
                "Min_Shared": ext_kept,
                "Min_Pairwise_Dist": min_dist if min_dist < 1 else 1.0 / max(min_dist, 1e-10),
            }
        
        self.node_stats = {}
        for node in self.tree_manager.tree.find_clades(order="postorder"):
            stats = calc_node_stats(node)
            self.node_stats[node] = stats

    ####################
    ## private clades ##
    #################### 
    def individual_shared_reads_summary(self, combinations_df: pd.DataFrame) -> pd.DataFrame:
        """
        Return dataframe of individual shared reads summary
        """

        individual_sets = []
        for leaf in combinations_df.accid_A.unique():
            pairs_set = combinations_df[
                (combinations_df.accid_A == leaf)
                | (combinations_df.accid_B == leaf)
            ]
            pairs_set = pairs_set.aggregate(
                {
                    "proportion_max": lambda x: max(x),
                    "proportion_min": lambda x: max(x),
                }
            )
            individual_sets.append(pairs_set)
        individual_sets = pd.DataFrame(individual_sets)
        individual_sets["proportion_max"] = individual_sets[
            "proportion_max"
        ].fillna(0)
        individual_sets["proportion_min"] = individual_sets[
            "proportion_min"
        ].fillna(0)
        individual_sets["proportion_std"] = individual_sets["proportion_max"].std()

        return individual_sets

    def clade_shared_individual_summary(self, leaves: list) -> pd.DataFrame:
        combinations = self.clade_shared_by_pair(leaves)

        ########################################################
        ### calculate max per sample shared reads per sample
        individual_sets = self.individual_shared_reads_summary(combinations)
        return individual_sets

    def clade_shared_by_pair(self, leaves: list) -> pd.DataFrame:

        group = self.read_profile_matrix_filtered.loc[leaves]
        group_pairwise_shared = pairwise_shared_count(group)

        group_pairwise_shared /= group.sum(axis=1)

        combinations = []

        for i in range(group_pairwise_shared.shape[0]):
            group_pairwise_shared.iloc[i, i] = 0
            for j in range(i + 1, group_pairwise_shared.shape[1]):
                shared_ij = group_pairwise_shared.iloc[i, j]
                shared_ji = group_pairwise_shared.iloc[j, i]
                shared_pair = (shared_ij, shared_ji)

                combinations.append(
                    [
                        group_pairwise_shared.index[i],
                        group_pairwise_shared.columns[j],
                        max(shared_pair),
                        min(shared_pair),
                        np.std(shared_pair),
                    ]
                )

        combinations = pd.DataFrame(
            combinations,
            columns=[
                "accid_A",
                "accid_B",
                "proportion_max",
                "proportion_min",
                "proportion_std",
            ],
        )

        return combinations

    def clade_total_counts(self, leaves: list) -> float:
        """
        return total counts of clade
        """
        return int(self.read_profile_matrix_filtered.loc[leaves].sum().sum())

    def get_accession_private_counts(self, duplicate_group: tuple) -> int:
        duplicate_group_counts = self.read_profile_matrix_filtered.loc[
            self.read_profile_matrix_filtered.index.isin(duplicate_group) == True
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
                    _,
                    _,
                    proportion_private,
                ) = clade_private_proportions(self.read_profile_matrix_filtered, leaves)

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

    def between_clade_shared_reads(self, clades_filter=[]) -> pd.DataFrame:
        """
        Return dataframe of pairwise shared reads between all pairs of clades"""
        clade_read_matrix = self.between_clade_reads_matrix(
            filter_names=clades_filter, remove_leaves=False
        )
        shared_clade_matrix = square_and_fill_diagonal(clade_read_matrix)

        return shared_clade_matrix

    def within_clade_shared_reads(self, clade) -> pd.DataFrame:
        """ """
        clade_node = [
            x for x in self.all_clade_leaves_filtered.keys() if x.name == clade
        ][0]
        leaves = self.all_clade_leaves_filtered[clade_node]

        clade_read_matrix = self.within_clade_reads_matrix(leaves)
        shared_clade_matrix = square_and_fill_diagonal(clade_read_matrix)
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

    def plot_pairwise_shared_clade_reads(
        self, pairwise_shared_clade: pd.DataFrame, subplot=False, clade_str=""
    ) -> None:
        """
        Plot heatmap of pairwise shared reads between all pairs of leaves
        """
        if pairwise_shared_clade.shape[0] < 2:
            return

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

        if self.read_profile_matrix_filtered.shape[1] <= 3:
            return

        if self.read_profile_matrix_filtered.shape[0] <= 3:
            return

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
                    private_counts=0,
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
            ) = clade_private_proportions(self.read_profile_matrix_filtered, leaves)
            total_proportion = total_reads / self.read_profile_matrix_filtered.shape[1]

            clade_counts = self.clade_total_counts(leaves)

            if len(leaves) == 1:
                node_stats_dict[node] = Clade(
                    name=node,
                    leaves=leaves,
                    group_counts=clade_counts,
                    private_counts=private_reads,
                    private_proportion=proportion_private,
                    total_proportion=total_proportion,
                    shared_proportion_std=0,
                    shared_proportion_min=0,
                    shared_proportion_max=0,
                    overlap_df=pd.DataFrame(),
                )

                continue
            
            combinations = self.clade_shared_by_pair(leaves)
            individual_sets = self.individual_shared_reads_summary(combinations)

            node_stats_dict[node] = Clade(
                name=node,
                leaves=leaves,
                private_proportion=proportion_private,
                total_proportion=total_proportion,
                group_counts=clade_counts,
                private_counts=private_reads,
                shared_proportion_min=min(individual_sets.proportion_max),
                shared_proportion_max=max(individual_sets.proportion_max),
                shared_proportion_std=np.std(individual_sets.proportion_max),
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
                    stats.private_counts,
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
                "private_counts",
                "shared_proportion_std",
                "shared_proportion_min",
                "shared_proportion_max",
            ],
        )

        return clade_summary

    def clade_dict_from_summary(
        self,
        clade_summary: pd.DataFrame,
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
            private_counts = 0
            if "private_counts" in clade_summary.columns:
                private_counts = int(stats.private_counts)

            node_stats_dict[node_phylo_clade] = Clade(
                name=str(node),
                leaves=inner_node_leaf_dict[node_phylo_clade],
                private_proportion=stats.private_proportion,
                total_proportion=stats.private_proportion,
                group_counts=stats.total_counts,
                private_counts=private_counts,
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
                            statistics_dict[clade].private_counts,
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
                            0,
                        )
                    )

        ## create dataframe
        leaf_clades_df = pd.DataFrame(
            leaf_clades_dict,
            columns=[
                "leaf",
                "clade",
                "total_counts",
                "private_counts",
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

    # ── ML-based clustering ──────────────────────────────────────────

    def _ml_client(self):
        from pathogen_identification.utilities.ml_api_client import MLAPIClient
        return MLAPIClient()

    def _get_accid_statistics(self) -> pd.DataFrame:
        if self._accid_statistics_df is None:
            self._accid_statistics_df = pd.read_csv(self.accid_statistics_path, sep="\t")
        return self._accid_statistics_df

    def _best_taxid_for_leaf(self, leaf: str) -> Optional[int]:
        row = self.metadata[self.metadata["accid"] == leaf]
        if row.empty:
            return None
        return int(row.iloc[0].get("taxid", None)) if "taxid" in row.columns else None

    def _best_taxid_for_leaves(self, leaves: list) -> Optional[int]:
        """Return best taxid across a set of leaves (by coverage from accid_statistics)."""
        try:
            accid_df = self._get_accid_statistics()
            subset = accid_df[accid_df["accid"].isin(leaves)]
            if subset.empty:
                return None
            subset = subset.sort_values("read_count", ascending=False)
            return int(subset.iloc[0]["taxid"]) if "accid" in subset.columns else None
        except Exception:
            return None

    def _node_features(
        self,
        clade: Phylo.BaseTree.Clade,
        leaves: List[str],
        accid_df: pd.DataFrame,
    ) -> Dict[str, float]:
        features: Dict[str, float] = {}

        features["n_leaves"] = float(len(leaves))

        shared_df = self.clade_shared_individual_summary(leaves)

        node_stats = self.node_stats.get(clade, {})
        features["Min_Shared"] = float(node_stats.get("Min_Shared", 0.0))
        features["Min_Dist"] = float(node_stats.get("Min_Dist", 0.0))

        if accid_df.empty:
            features["tax_diversity"] = 0.0
        else:

            features["tax_diversity"] = shannon_diversity_from_list(
                accid_df[
                    accid_df['accid'].isin(leaves)
                ]['order'].dropna().tolist()
            )

        return features

    def _predict_stop_traversal(
        self, features: Dict[str, float], model_type: str
    ) -> Tuple[Optional[bool], float]:
        try:
            client = self._ml_client() 

            result = client.predict_composition_stop_traversal(features, model=model_type)
            return bool(result["stop_traversal"]), float(result.get("probability", 0.0))
        except Exception:
            import traceback
            traceback.print_exc()

            self.logger.warning("ML API unavailable, falling back to fixed traversal")
            return None, 0.0

    def _traverse_with_fixed(
        self,
        clade: Phylo.BaseTree.Clade,
        accid_df: pd.DataFrame,
        results: List[Dict[str, Any]],
        min_shared_threshold: Optional[float] = None,
    ):
        if min_shared_threshold is None:
            min_shared_threshold = self.clade_filter.reference_clade.shared_proportion_min

        leaves = [l.name for l in self.tree_manager.get_node_leaves(clade) if l.name]

        if not leaves:
            return

        shared_df = self.clade_shared_by_pair(leaves)
        min_shared = float(shared_df["proportion_max"].min()) if not shared_df.empty else 0.0

        best_taxid = self._best_taxid_for_leaves(leaves)

        if min_shared >= min_shared_threshold or clade.is_terminal():
            results.append(
                {
                    "node": clade.name or f"clade_{len(results)}",
                    "n_leaves": len(leaves),
                    "leaves": leaves,
                    "best_taxid_match": best_taxid,
                }
            )
        else:
            for child in clade.clades:
                if child.is_terminal():
                    taxid = self._best_taxid_for_leaf(child.name)
                    results.append(
                        {
                            "node": child.name or f"leaf_{child.name}",
                            "n_leaves": 1,
                            "leaves": [child.name],
                            "best_taxid_match": taxid,
                        }
                    )
                else:
                    self._traverse_with_fixed(
                        child, accid_df, results, min_shared_threshold=min_shared_threshold
                    )

    def _traverse_with_prediction(
        self,
        clade: Phylo.BaseTree.Clade,
        accid_df: pd.DataFrame,
        model_type: str,
        results: Optional[List[Dict[str, Any]]] = None,
    ):
        if results is None:
            results = []

        if self.clustering_model_type == "Fixed":
            self._traverse_with_fixed(clade, accid_df, results)
        else:
            self._traverse_with_api_prediction(clade, accid_df, model_type, results)

    def _traverse_with_api_prediction(
        self,
        clade: Phylo.BaseTree.Clade,
        accid_df: pd.DataFrame,
        model_type: str,
        results: Optional[List[Dict[str, Any]]] = None,
    ):
        if results is None:
            results = []
        
        leaves = [l.name for l in self.tree_manager.get_node_leaves(clade) if l.name]

        if not leaves:
            return results

        features = self._node_features(clade, leaves, accid_df)
        stop, prob = self._predict_stop_traversal(features, model_type)

        if stop is None:
            raise RuntimeError("ML API unavailable and no fallback model provided.")

        best_taxid = self._best_taxid_for_leaves(leaves)

        if stop or clade.is_terminal():
            results.append(
                {
                    "node": clade.name or f"clade_{len(results)}",
                    "n_leaves": len(leaves),
                    "leaves": leaves,
                    "best_taxid_match": best_taxid,
                }
            )
        else:
            for child in clade.clades:
                if child.is_terminal():
                    taxid = self._best_taxid_for_leaf(child.name)
                    results.append(
                        {
                            "node": child.name or f"leaf_{child.name}",
                            "n_leaves": 1,
                            "leaves": [child.name],
                            "best_taxid_match": taxid,
                        }
                    )
                else:
                    self._traverse_with_api_prediction(child, accid_df, model_type, results)

        return results

    def predict_clades_composition(self, model_type: str) -> pd.DataFrame:
        accid_df = self._get_accid_statistics()
        from pathogen_identification.modules.metadata_handler import RunMetadataHandler
        accid_df['taxid'] = accid_df['accid'].apply(lambda x: RunMetadataHandler.get_accid_taxid(x) if pd.notnull(x) else None)
        accid_df['family'] = accid_df['taxid'].apply(lambda x: RunMetadataHandler._get_taxid_taxonomy(x, 'family') if pd.notnull(x) else None)
        accid_df['order'] = accid_df['taxid'].apply(lambda x: RunMetadataHandler._get_taxid_taxonomy(x, 'order') if pd.notnull(x) else None)
        cluster_results: List[Dict[str, Any]] = []
        root = self.tree_manager.tree.root
        self._traverse_with_prediction(root, accid_df, model_type, cluster_results)

        rows = []

        for cluster in cluster_results:
            for leaf in cluster["leaves"]:
                rows.append(
                    {
                        "leaf": leaf,
                        "clade": cluster["node"],
                        "total_counts": cluster.get("n_leaves", 0),
                        "private_counts": 0,
                        "private_proportion": 0.0,
                        "shared_proportion": 0.0,
                    }
                )

        leaf_clades_df = pd.DataFrame(
            rows,
            columns=["leaf", "clade", "total_counts", "private_counts", "private_proportion", "shared_proportion"],
        )

        try:
            accids_df = self._get_accid_statistics()
            leaf_clades_df = leaf_clades_df.merge(accids_df, right_on="accid", left_on="leaf", how="left")
        except Exception:
            pass

        if "clade" in leaf_clades_df.columns:
            leaf_clades_df["clade"] = leaf_clades_df["clade"].fillna(leaf_clades_df["leaf"])

        sort_cols = [c for c in ["total_counts", "clade", "read_count"] if c in leaf_clades_df.columns]
        if sort_cols:
            leaf_clades_df.sort_values(by=sort_cols, ascending=[False, True, False], inplace=True)
        leaf_clades_df.reset_index(drop=True, inplace=True) 

        return leaf_clades_df

    def get_leaf_clades(self, force=False) -> pd.DataFrame:
        #if self.clustering_model_type and self.clustering_model_type != "Fixed":
        self.recalculate_all_min_pairwise_dist()
        
        clades_df = self.predict_clades_composition(self.clustering_model_type)

        statistics_dict_all = self.get_node_statistics(force=force)

        return clades_df

