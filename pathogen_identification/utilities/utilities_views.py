

import pandas as pd
import os
from pathogen_identification.models import FinalReport, ReferenceMap_Main
from pathogen_identification.utilities.phylo_tree import PhyloTreeManager
from pathogen_identification.utilities.utilities_general import simplify_name
from typing import List
from pathogen_identification.utilities.overlap_manager import ReadOverlapManager
# import Django BaseManager
#from django.db.models import BaseManager


class ReportSorter:

    def __init__(self, reports: List[FinalReport], threshold: float):
        self.reports= reports
        self.threshold= threshold
        self.report_dict= {report.accid: report for report in reports}
        self.metadata_df= self.prep_metadata_df()

        self.fasta_files= self.metadata_df.file.tolist()
    
    def retrieved_mapped_subset(self, report: FinalReport):
        """
        Return subset of retrieved and mapped reads
        """
        try:
            simple_accid= simplify_name(report.accid)

            mapped_ref= ReferenceMap_Main.objects.get(reference= simple_accid, run= report.run)

            if not mapped_ref.mapped_subset_r1_fasta:
                return None
            
            if os.path.exists(mapped_ref.mapped_subset_r1_fasta):
                return mapped_ref.mapped_subset_r1_fasta

        except ReferenceMap_Main.DoesNotExist:
            return None


    def prep_metadata_df(self):
        """
        Return metadata dataframe
        columns: filename, description, accid
        """
        metadata_dict= []

        for report in self.reports:
            accid= report.accid
            description= report.description
            mapped_subset_r1= self.retrieved_mapped_subset(report)
            if not mapped_subset_r1:
                continue
            filename= mapped_subset_r1
            metadata_dict.append({"file": filename, "description": description, "accid": accid})
        
        metadata_df= pd.DataFrame(metadata_dict)
        

        if metadata_df.empty:
            metadata_df= pd.DataFrame(columns=["file", "description", "accid"])

        return metadata_df
    

    def read_overlap_analysis(self):
        """
        Return read overlap analysis as dataframe
        columns: leaf (accid), clade, read_count, group_count
        """
        overlap_manager= ReadOverlapManager(self.fasta_files, self.metadata_df, threshold= self.threshold)

        njtree= overlap_manager.generate_tree()

        ### inner node to leaf dict
        tree_manager= PhyloTreeManager(njtree)
        inner_node_leaf_dict= tree_manager.clades_get_leaves_clades()

        private_read_dict= overlap_manager.node_private_reads(inner_node_leaf_dict)
        private_clades = overlap_manager.filter_clades_by_private_reads(private_read_dict)
        leaf_clades= tree_manager.leaf_clades_clean(private_clades)
        clades= overlap_manager.leaf_clades_to_pandas(leaf_clades)
        
        return clades

    def sort_reports(self):
        """
        Return sorted reports
        """

        if self.metadata_df.empty:
            return [self.reports]

        overlap_analysis= self.read_overlap_analysis()
        overlap_groups= list(overlap_analysis.groupby(["group_count","clade"]))[::-1]
        sorted_reports= []

        for group in overlap_groups:
            group_df= group[1]
            group_accids= group_df.leaf.tolist()
            sorted_reports.append([self.report_dict[accid] for accid in group_accids])
        
        return sorted_reports







            

