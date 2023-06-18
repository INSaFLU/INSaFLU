import pandas as pd
import os
from pathogen_identification.models import FinalReport, ReferenceMap_Main
from pathogen_identification.utilities.phylo_tree import PhyloTreeManager
from pathogen_identification.utilities.utilities_general import (
    simplify_name,
    infer_run_media_dir,
)
from typing import List
from pathogen_identification.utilities.overlap_manager import ReadOverlapManager
from pathogen_identification.utilities.clade_objects import Clade

# import Django BaseManager
# from django.db.models import BaseManager


class ReportSorter:
    analysis_filename = "overlap_analysis.tsv"

    def __init__(
        self,
        reports: List[FinalReport],
        private_threshold: float,
        min_shared_threshold: float,
        force=False,
    ):
        self.reports = reports
        self.threshold = private_threshold
        self.reference_clade = self.generate_reference_clade(
            private_threshold, min_shared_threshold
        )
        self.report_dict = {report.accid: report for report in reports}
        self.metadata_df = self.prep_metadata_df()

        self.fasta_files = self.metadata_df.file.tolist()
        self.run = self.infer_run()
        self.run_media_dir = self.inferred_run_media_dir()
        self.analysis_df_path = os.path.join(self.run_media_dir, self.analysis_filename)

        self.force = force

    @staticmethod
    def generate_reference_clade(private_threshold: float, min_shared_threshold: float):
        """
        Return reference clade"""
        ref_clade = Clade(
            name="ref",
            leaves=[],
            private_proportion=private_threshold,
            shared_proportion_std=min_shared_threshold,
            shared_proportion_min=min_shared_threshold,
            shared_proportion_max=min_shared_threshold,
        )
        return ref_clade

    def infer_run(self):
        """
        Return run
        """
        if not self.reports:
            return None

        return self.reports[0].run

    def inferred_run_media_dir(self):
        """
        Return run media directory
        """
        if not self.run:
            return None

        return infer_run_media_dir(self.run)

    def retrieved_mapped_subset(self, report: FinalReport):
        """
        Return subset of retrieved and mapped reads
        """
        try:
            simple_accid = simplify_name(report.accid)

            mapped_ref = ReferenceMap_Main.objects.get(
                reference=simple_accid, run=report.run
            )

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
        metadata_dict = []

        for report in self.reports:
            accid = report.accid
            description = report.description
            mapped_subset_r1 = self.retrieved_mapped_subset(report)
            if not mapped_subset_r1:
                continue
            filename = mapped_subset_r1
            metadata_dict.append(
                {"file": filename, "description": description, "accid": accid}
            )

        metadata_df = pd.DataFrame(metadata_dict)

        if metadata_df.empty:
            metadata_df = pd.DataFrame(columns=["file", "description", "accid"])

        return metadata_df

    def read_overlap_analysis(self):
        """
        Return read overlap analysis as dataframe
        columns: leaf (accid), clade, read_count, group_count
        """
        overlap_manager = ReadOverlapManager(
            self.fasta_files, self.metadata_df, threshold=self.threshold
        )

        njtree = overlap_manager.generate_tree()

        ### inner node to leaf dict
        tree_manager = PhyloTreeManager(njtree)
        inner_node_leaf_dict = tree_manager.clades_get_leaves_clades()

        # private_read_dict = overlap_manager.node_private_reads(inner_node_leaf_dict)
        # private_clades = overlap_manager.filter_clades_by_private_reads(
        #    private_read_dict
        # )

        statistics_dict = overlap_manager.node_statistics(inner_node_leaf_dict)

        selected_clades = overlap_manager.filter_clades(statistics_dict)

        leaf_clades = tree_manager.leaf_clades_clean(selected_clades)
        clades = overlap_manager.leaf_clades_to_pandas(leaf_clades)

        return clades

    def check_all_accids_analyzed(self, df: pd.DataFrame):
        """
        Return True if all accids have been analyzed
        """
        for accid in self.report_dict:
            if accid not in df.leaf.tolist():
                return False

        return True

    def check_analyzed(self):
        """
        Return True if all accids have been analyzed
        """
        if not os.path.isfile(self.analysis_df_path):
            return False

        try:
            analysis_df = pd.read_csv(self.analysis_df_path, sep="\t")

            if not os.path.isfile(analysis_df):
                return True

            if self.check_all_accids_analyzed(analysis_df):
                return False

            return True

        except pd.errors.EmptyDataError:
            return False

    def sort_reports(self):
        """
        Return sorted reports
        """

        overlap_analysis = self.read_overlap_analysis()
        overlap_analysis.to_csv(self.analysis_df_path, sep="\t", index=False)

    def get_reports(self):
        """
        Return sorted reports
        """

        if self.metadata_df.empty:
            return [self.reports]

        if not self.check_all_accids_analyzed(self.analysis_df_path):
            return [self.reports]

        overlap_analysis = pd.read_csv(self.analysis_df_path, sep="\t")

        overlap_groups = list(overlap_analysis.groupby(["group_count", "clade"]))[::-1]

        sorted_reports = []

        for group in overlap_groups:
            group_df = group[1]
            group_accids = group_df.leaf.tolist()
            sorted_reports.append([self.report_dict[accid] for accid in group_accids])

        return sorted_reports
