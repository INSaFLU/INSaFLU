import os
from typing import List

import pandas as pd

from pathogen_identification.models import FinalReport, ReferenceMap_Main
from pathogen_identification.utilities.overlap_manager import ReadOverlapManager
from pathogen_identification.utilities.phylo_tree import PhyloTreeManager
from pathogen_identification.utilities.utilities_general import (
    infer_run_media_dir,
    simplify_name,
)


def set_control_reports(project_pk: int):
    """
    set control reports
    """

    try:
        project = Projects.objects.get(pk=project_pk)

        control_samples = PIProject_Sample.objects.filter(
            project=project, is_control=True
        )

        control_reports = FinalReport.objects.filter(sample__in=control_samples)

        control_report_taxids = control_reports.values_list("taxid", flat=True)
        control_report_taxids_set = set(control_report_taxids)
        print(control_report_taxids_set)

        other_reports = FinalReport.objects.filter(sample__project=project).exclude(
            sample__in=control_samples
        )

        for sample_report in other_reports:
            if sample_report.taxid in control_report_taxids_set:
                sample_report.control_flag = FinalReport.CONTROL_FLAG_PRESENT
            else:
                sample_report.control_flag = FinalReport.CONTROL_FLAG_NONE

            sample_report.save()

        for report in control_reports:
            report.control_flag = FinalReport.CONTROL_FLAG_NONE
            report.save()

    except Exception as e:
        print(e)
        pass


class ReportSorter:
    analysis_filename = "overlap_analysis.tsv"

    def __init__(self, reports: List[FinalReport], threshold: float, force=False):
        self.reports = reports
        self.threshold = threshold
        self.report_dict = {
            report.accid: report
            for report in reports
            if self.retrieved_mapped_subset(report)
        }
        self.metadata_df = self.prep_metadata_df()

        self.fasta_files = self.metadata_df.file.tolist()
        self.run = self.infer_run()
        self.run_media_dir = self.inferred_run_media_dir()
        self.analysis_df_path = os.path.join(self.run_media_dir, self.analysis_filename)

        self.force = force

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

        private_read_dict = overlap_manager.node_private_reads(inner_node_leaf_dict)
        private_clades = overlap_manager.filter_clades_by_private_reads(
            private_read_dict
        )
        leaf_clades = tree_manager.leaf_clades_clean(private_clades)
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

            if self.check_all_accids_analyzed(analysis_df):
                return True

            return False

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

        if not self.check_analyzed():
            return [self.reports]

        overlap_analysis = pd.read_csv(self.analysis_df_path, sep="\t")

        overlap_groups = list(overlap_analysis.groupby(["group_count", "clade"]))[::-1]

        sorted_reports = []

        for group in overlap_groups:
            group_df = group[1]
            group_accids = group_df.leaf.tolist()
            sorted_reports.append([self.report_dict[accid] for accid in group_accids])

        return sorted_reports


import os
from typing import List

import pandas as pd

from pathogen_identification.models import FinalReport, ReferenceMap_Main
from pathogen_identification.utilities.overlap_manager import ReadOverlapManager
from pathogen_identification.utilities.phylo_tree import PhyloTreeManager
from pathogen_identification.utilities.utilities_general import (
    infer_run_media_dir,
    simplify_name,
)

# import Django BaseManager
# from django.db.models import BaseManager


class ReportSorter:
    analysis_filename = "overlap_analysis.tsv"

    def __init__(self, reports: List[FinalReport], threshold: float, force=False):
        self.reports = reports
        self.threshold = threshold
        self.report_dict = {
            report.accid: report
            for report in reports
            if self.retrieved_mapped_subset(report)
        }
        self.metadata_df = self.prep_metadata_df()

        self.fasta_files = self.metadata_df.file.tolist()
        self.run = self.infer_run()
        self.run_media_dir = self.inferred_run_media_dir()
        self.analysis_df_path = os.path.join(self.run_media_dir, self.analysis_filename)

        self.force = force

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

        private_read_dict = overlap_manager.node_private_reads(inner_node_leaf_dict)
        private_clades = overlap_manager.filter_clades_by_private_reads(
            private_read_dict
        )
        leaf_clades = tree_manager.leaf_clades_clean(private_clades)
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

            if self.check_all_accids_analyzed(analysis_df):
                return True

            return False

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

        if not self.check_analyzed():
            return [self.reports]

        overlap_analysis = pd.read_csv(self.analysis_df_path, sep="\t")

        overlap_groups = list(overlap_analysis.groupby(["group_count", "clade"]))[::-1]

        sorted_reports = []

        for group in overlap_groups:
            group_df = group[1]
            group_accids = group_df.leaf.tolist()
            sorted_reports.append([self.report_dict[accid] for accid in group_accids])

        return sorted_reports
