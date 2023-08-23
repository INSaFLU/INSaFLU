import os
from typing import List

import pandas as pd

from pathogen_identification.models import (
    FinalReport,
    ReferenceMap_Main,
    Projects,
    PIProject_Sample,
    RunAssembly,
    RunMain,
)
from pathogen_identification.utilities.televir_parameters import LayoutParams

from pathogen_identification.utilities.overlap_manager import ReadOverlapManager
from pathogen_identification.utilities.phylo_tree import PhyloTreeManager
from pathogen_identification.utilities.utilities_general import (
    infer_run_media_dir,
    simplify_name,
)
from pathogen_identification.utilities.televir_parameters import TelevirParameters

from pathogen_identification.utilities.clade_objects import (
    Clade,
)
import datetime

def set_control_reports(project_pk: int):
    """
    set control reports
    """

    try:
        project = Projects.objects.get(pk=project_pk)

        control_samples = PIProject_Sample.objects.filter(
            project=project, is_control=True
        )

        control_reports = FinalReport.objects.filter(sample__in=control_samples).distinct("taxid")

        control_report_taxids = control_reports.values_list("taxid", flat=True)
        control_report_taxids_set = set(control_report_taxids)

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
            report.control_flag = FinalReport.CONTROL_FLAG_SOURCE
            report.save()

    except Exception as e:
        print(e)
        pass


def recover_assembly_contigs(run_main: RunMain, run_assembly: RunAssembly):
    """
    check contigs exist, if not, replace path with media path check again, if so, replace with media path.
    """
    ##
    assembly_contigs = run_assembly.assembly_contigs

    if not assembly_contigs:
        return

    assembly_contigs_exist = os.path.exists(assembly_contigs)

    if assembly_contigs_exist:
        return

    media_dir = infer_run_media_dir(run_main)

    if not media_dir:
        return

    if not assembly_contigs_exist:
        assembly_contigs = os.path.basename(assembly_contigs)
        assembly_contigs = os.path.join(media_dir, "assembly", assembly_contigs)
        assembly_contigs_exist = os.path.exists(assembly_contigs)
        if assembly_contigs_exist:
            run_assembly.assembly_contigs = assembly_contigs
            run_assembly.save()



class ReportSorter:
    analysis_filename = "overlap_analysis_{}.tsv"
    all_clade_filename = "all_clades_{}.tsv"

    def __init__(
        self,
        reports: List[FinalReport],
        report_layout_params: LayoutParams,
        force=False,
        level= 0,
    ):
        self.reports = reports
        self.reference_clade = self.generate_reference_clade(report_layout_params)
        self.report_dict = {
            report.accid: report
            for report in reports
            if self.retrieved_mapped_subset(report)
        }
        self.excluded_dict= {
            report.accid: report
            for report in reports
            if not self.retrieved_mapped_subset(report)
        }
        self.metadata_df = self.prep_metadata_df()

        self.fasta_files = self.metadata_df.file.tolist()
        self.level= level
        
        self.model= self.set_level(reports, level)
        self.run = self.infer_run()

        if self.model is not None:
            self.media_dir= self.infer_media_dir()
            #self.run_media_dir = self.inferred_run_media_dir()

            self.all_clades_df_path = os.path.join(
                self.media_dir, self.all_clade_filename
            )

            self.analysis_df_path = os.path.join(
                self.media_dir,
                self.analysis_filename.format(
                    report_layout_params.read_overlap_threshold
                ),
            )
            self.force = force

        else:
            self.media_dir = None
            self.analysis_df_path = None
            self.all_clades_df_path = None
            self.force = False

    def set_level(self, final_report_list: List[FinalReport], level):
        
        if len(final_report_list) == 0:
            return None
        final_report= final_report_list[0]
        if level == 0:
            return final_report.sample
        elif level == 1:
            return final_report.run
        
        return None

    def infer_media_dir(self):
        """
        Return media directory
        """
        if self.level == 0:
            return self.inferred_sample_media_dir()
        elif self.level == 1:
            return self.inferred_run_media_dir()
    
    @staticmethod
    def generate_reference_clade(layout_params: LayoutParams):
        """
        Return reference clade"""
        ref_clade = Clade(
            name="ref",
            leaves=[],
            private_proportion=layout_params.read_overlap_threshold,
            group_counts=0,
            shared_proportion_std=layout_params.shared_proportion_threshold,
            shared_proportion_min=layout_params.shared_proportion_threshold,
            shared_proportion_max=layout_params.shared_proportion_threshold,
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
        if not self.model:
            raise Exception("No run found")

        return infer_run_media_dir(self.run)
    
    def inferred_sample_media_dir(self):
        """
        Return run media directory
        """
        if not self.model:
            raise Exception("No model found")

        rundir= infer_run_media_dir(self.run)
        sample_dir= os.path.dirname(rundir)

        return sample_dir
        

    def retrieved_mapped_subset(self, report: FinalReport):
        """
        Return subset of retrieved and mapped reads
        """
        if report.depth == 0:
            return None

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
            self.metadata_df,
            self.reference_clade,
            self.media_dir,
            str(self.model.pk),
        )
        ### time operations
        print("generating tree")

        ### generate tree
        start= datetime.datetime.now()

        njtree = overlap_manager.generate_tree()
        end = datetime.datetime.now()
        
        # time in seconds
        time = (end-start).total_seconds()
        print("time to generate tree: ", time)

        ### inner node to leaf dict
        tree_manager = PhyloTreeManager(njtree)
        # inner_node_leaf_dict = tree_manager.clades_get_leaves_clades()
        all_node_leaves = tree_manager.all_clades_leaves()

        ### get statistics
        statistics_dict_all = overlap_manager.get_node_statistics(
            njtree, all_node_leaves
        )
        end= datetime.datetime.now()
        # time in seconds
        time = (end-start).total_seconds()
        #print("time to get statistics: ", time)

        selected_clades = overlap_manager.filter_clades(statistics_dict_all)

        leaf_clades = tree_manager.leaf_clades_clean(selected_clades)
        clades = overlap_manager.leaf_clades_to_pandas(leaf_clades, statistics_dict_all)

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

        overlap_manager = ReadOverlapManager(
            self.metadata_df,
            self.reference_clade,
            self.media_dir,
            str(self.model.pk),
        )

        if not os.path.exists(self.analysis_df_path):
            return False

        if not os.path.exists(overlap_manager.distance_matrix_path):
            return False

        if not os.path.exists(overlap_manager.clade_statistics_path):
            return False

        if not os.path.exists(overlap_manager.clade_statistics_path):
            return False

        if not overlap_manager.all_accs_analyzed():
            return False

        return True

    def sort_reports_save(self):
        """
        Return sorted reports
        """
        if self.model is not None:
            overlap_analysis = self.read_overlap_analysis()
            overlap_analysis.to_csv(self.analysis_df_path, sep="\t", index=False)

    def get_sorted_reports(self) -> List[List[FinalReport]]:
        overlap_analysis = self.read_overlap_analysis()

        overlap_groups = list(overlap_analysis.groupby(["total_counts", "clade"]))[::-1]

        sorted_reports = []

        for group in overlap_groups:
            group_df = group[1]

            group_accids = group_df.leaf.tolist()
            group_accids = [x for x in group_accids if x in self.report_dict]
            group_list = [self.report_dict[accid] for accid in group_accids]
            # sort by coverage
            group_list.sort(key=lambda x: x.coverage, reverse=True)
            if len(group_list):
                sorted_reports.append(group_list)

        # sort groups by max coverage among group

        def get_max_coverage(group: List[FinalReport]):
            return max(y.coverage for y in group)

        sorted_groups = sorted(sorted_reports, key=get_max_coverage, reverse=True)

        return sorted_groups

    def get_reports(self):
        """
        Return sorted reports
        """
        if self.model is None:
            return [self.reports]

        if self.metadata_df.empty:
            return [self.reports]

        if not self.check_analyzed():
            return [self.reports]

        return self.get_sorted_reports()
    
    def check_excluded_exist(self) -> bool:
        """return True if there are excluded reports"""

        if len(self.excluded_dict) > 0:
            return True

        return False

    def get_reports_empty(self) -> List[FinalReport]:
        """return reports in excluded report_dict"""

        return list(self.excluded_dict.values())

def calculate_reports_overlaps(sample: PIProject_Sample):
    final_reports = FinalReport.objects.filter(sample=sample)
    report_layout_params = TelevirParameters.get_report_layout_params(
        project_pk=sample.project.pk
    )
    report_sorter = ReportSorter(final_reports, report_layout_params)
    report_sorter.sort_reports_save()

from typing import List, Union
from django.db.models.query import QuerySet
from braces.views import FormValidMessageMixin, LoginRequiredMixin
from django.views import generic
from django.views.generic import ListView
from pathogen_identification.models import RunMain, RunDetail, FinalReport, ParameterSet


class FinalReportCompound(LoginRequiredMixin, generic.TemplateView):
    def __init__(self, report: FinalReport):
        """
        copy all attributes from report
        """

        for attr in dir(FinalReport):
            if not attr.startswith("__"):
                if attr == "objects":
                    continue
                try:
                    setattr(self, attr, getattr(report, attr))
                except Exception as e:
                    raise e

        self.found_in = self.get_identical_reports_ps(report)
        self.run_detail = self.get_report_rundetail(report)
        self.run_main = self.get_report_runmain(report)
        self.run_index = self.run_main.pk
        self.data_exists = self.check_data_exists(report)
        self.control_flag = report.control_flag


    def get_identical_reports_ps(self, report: FinalReport) -> list:
        reports_unique = FinalReport.objects.filter(
            run__project__pk=report.run.project.pk,
            sample__pk=report.sample.pk,
            accid=report.accid,
        )

        sets = [r.run.parameter_set.leaf.index for r in reports_unique]
        return ", ".join([str(s) for s in sets])

    def check_data_exists(self, report: FinalReport) -> bool:
        return report.run.data_deleted == False

    def get_report_rundetail(self, report: FinalReport) -> RunDetail:
        return RunDetail.objects.get(sample=report.sample, run=report.run)

    def get_report_runmain(self, report: FinalReport) -> RunMain:
        return report.run


def final_report_best_cov_by_accid(reports: QuerySet) -> QuerySet:
    """
    get the best coverage report for each accid
    """

    pk_to_keep = []
    for accid in set(reports.values_list("accid", flat=True)):
        best_report = reports.filter(accid=accid).order_by("-coverage").first()
        pk_to_keep.append(best_report.pk)

    return reports.filter(pk__in=pk_to_keep)
