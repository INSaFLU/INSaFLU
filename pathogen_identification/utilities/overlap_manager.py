
import itertools as it
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from scipy.spatial.distance import pdist, squareform
import pandas as pd
import os
## pairwise matrix by individual reads
from pathogen_identification.utilities.utilities_general import readname_from_fasta

from typing import List

def accid_from_metadata(metadata: pd.DataFrame, read_name: str) -> str:
    """
    Return accession id of read from metadata
    """
    filename= os.path.basename(read_name)
    try:
        return metadata[metadata["filename"] == filename]["accid"].values[0]
    except IndexError:
        return read_name


class ReadOverlapManager:

    def __init__(self, fasta_list: List[str], metadata_df: pd.DataFrame, threshold: float= 1):

        self.metadata = metadata_df
        self.threshold = threshold
        self.fasta_list = fasta_list

        self.metadata["filename"]= self.metadata["file"].apply(lambda x: x.split("/")[-1])

        self.read_profile_matrix: pd.DataFrame = self.generate_read_matrix()

    def get_accid_readname_dict(self):
        """
        Return dictionary of read names and descriptions
        """
        readname_dict = {}
        for fasta_file in self.fasta_list:
            accid= accid_from_metadata(self.metadata, fasta_file)
            read_names = readname_from_fasta(fasta_file)
            
            readname_dict[accid] = read_names
        return readname_dict

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

    
    def read_profile_dict_get(self,readname_dict: dict, all_reads: list) -> dict:
        """
        Return dictionary of read profiles for all accids
        """
        read_profile_dict = {}
        for accid in readname_dict.keys():
            read_profile_dict[accid] = self.read_profile_get(accid, readname_dict, all_reads)
        return read_profile_dict

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
        return pd.DataFrame(squareform(pdist(read_profile_matrix, metric = "jaccard")), columns = read_profile_matrix.index, index = read_profile_matrix.index)

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
        distmat= distance_matrix.values.tolist()
        distmat= [x[:i+1] for i,x in enumerate(distmat)]
        distmat= DistanceMatrix(list(distance_matrix.index), distmat)

        return distmat


    def tree_from_distance_matrix(self, distance_matrix: pd.DataFrame):
        """
        Return tree from distance matrix
        """
        distmat= self.matrix_to_phylotriangle(distance_matrix)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(distmat)
        tree.rooted= True
        tree.ladderize()
        return tree



    ####################
    ## Construct tree ##
    ####################

    def generate_read_matrix(self):
        """
        Generate read matrix
        """
        files_readnames = [readname_from_fasta(fasta_file) for fasta_file in self.fasta_list]
        readname_dict= self.get_accid_readname_dict()
        all_reads= self.all_reads_set(files_readnames)
        read_profile_dict= self.read_profile_dict_get(readname_dict, all_reads)
        read_profile_matrix= self.read_profile_matrix_get(read_profile_dict)
        return read_profile_matrix
    
    def generate_distance_matrix(self):
        """
        Generate distance matrix
        """
        distance_matrix= self.pairwise_shared_reads(self.read_profile_matrix)
        return distance_matrix
    
    def generate_tree(self):
        """
        Generate tree
        """
        distance_matrix= self.generate_distance_matrix()
        tree= self.tree_from_distance_matrix(distance_matrix)
        return tree

    ####################
    ## private clades ##
    ####################

    
    def node_private_reads(self, inner_node_leaf_dict: dict) -> dict:

        node_private_dict= {}
        for node, leaves in inner_node_leaf_dict.items():
            group= self.read_profile_matrix.loc[leaves]
            group_sum= group.sum(axis=0)
            group_sum_as_bool= group_sum > 0
            group_sum_as_bool_list= group_sum_as_bool.tolist()

            group_reads= self.read_profile_matrix.iloc[:,group_sum_as_bool_list]
            group_reads_sum_all= group_reads.sum(axis=0)
            group_reads_sum_group= group_reads.loc[leaves].sum(axis=0)
            group_reads_sum_group= group_reads_sum_group.fillna(0)

            read_proportions= group_reads_sum_group / group_reads_sum_all
            read_proportions= read_proportions.fillna(0)

            # proportion of reads in group that are private
            proportion_private= read_proportions.sum() / len(read_proportions)
            node_private_dict[node]= proportion_private

        return node_private_dict

    
    def filter_clades_by_private_reads(self, private_read_dict: dict) -> list:
        """
        Return list of clades with private reads above threshold
        """
        return [node for node, proportion in private_read_dict.items() if proportion >= self.threshold]


    def leaf_clades_to_pandas(self, leaf_clades: Dict[str, Phylo.BaseTree.Clade]) -> pd.DataFrame:
        """
        Return dataframe of leaf clades
        """

        leaf_clades_dict= [(leaf, clade.name) for leaf, clade in leaf_clades.items()]
        leaf_clades_df= pd.DataFrame(leaf_clades_dict, columns=["leaf", "clade"])
        leaf_clades_df["read_count"]= leaf_clades_df["leaf"].apply(lambda x: self.get_accession_total_counts(x))
        # sort by clade and then read count
        
        def group_count(clade):
            return leaf_clades_df[leaf_clades_df["clade"]==clade].read_count.sum()

        leaf_clades_df["group_count"]= leaf_clades_df["clade"].apply(group_count)

        def set_single_count_to_zero(row):
            if row.clade == "single":
                row.group_count = 0
            return row
        
        leaf_clades_df.apply(set_single_count_to_zero, axis=1)

        leaf_clades_df.sort_values(by=["group_count","clade", "read_count"], ascending= [False, True, False], inplace=True)
        
        leaf_clades_df.reset_index(drop=True, inplace=True)
        return leaf_clades_df
