import datetime
import logging
import os
from typing import Dict, List, Union

import pandas as pd
from braces.views import FormValidMessageMixin, LoginRequiredMixin
from django.db.models.query import QuerySet
from django.views import generic
from django.views.generic import ListView

from pathogen_identification.models import (
    FinalReport,
    PIProject_Sample,
    Projects,
    ReferenceMap_Main,
    RunAssembly,
    RunDetail,
    RunMain,
)
from pathogen_identification.utilities.clade_objects import Clade
from pathogen_identification.utilities.overlap_manager import ReadOverlapManager
from pathogen_identification.utilities.phylo_tree import PhyloTreeManager
from pathogen_identification.utilities.televir_parameters import (
    LayoutParams,
    TelevirParameters,
)
from pathogen_identification.utilities.utilities_general import (
    infer_run_media_dir,
    simplify_name,
)
from settings.constants_settings import ConstantsSettings
from settings.models import Parameter, Software


class EmptyRemapMain:
    run = None
    sample = None
    merged_log = None
    performed = False
    found_total = 0
    coverage_minimum = 0
    coverage_maximum = 0
    success = False
    coverage_mean = ""


class FinalReportWrapper:
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

        self.private_reads = 0

    def update_private_reads(self, private_reads: int):
        self.private_reads = private_reads


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
        self.private_reads = 0

    def update_private_reads(self, private_reads: int):
        self.private_reads = private_reads

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


class FinalReportGroup:
    name: str
    total_counts: str
    shared_proportion: str
    private_proportion: str
    group_list: List[FinalReportWrapper]

    def __init__(
        self,
        name: str,
        total_counts: int,
        shared_proportion: float,
        private_proportion: float,
        group_list: List[FinalReportWrapper],
        heatmap_path: str = "",
        heatmap_exists: bool = False,
    ):
        self.name = name
        self.total_counts = f"total counts {total_counts}"
        self.shared_proportion = f"shared proportion {shared_proportion:.2f}"
        self.private_proportion = f"private proportion {private_proportion:.2f}"
        self.group_list = group_list
        self.heatmap_path = heatmap_path
        self.heatmap_exists = heatmap_exists
        self.max_private_reads = 0

    def reports_have_private_reads(self) -> bool:
        for report in self.group_list:
            if report.private_reads > 0:
                return True
        return False

    def update_max_private_reads(self):
        for report in self.group_list:
            if report.private_reads > self.max_private_reads:
                self.max_private_reads = report.private_reads


def check_sample_software_exists(sample: PIProject_Sample) -> bool:
    """
    Return True if sample software exists
    """

    parameters_exist = Parameter.objects.filter(
        televir_project=sample.project,
        televir_project_sample=sample,
    ).exists()

    return parameters_exist


def check_project_params_exist(project: Projects) -> bool:
    """
    check if project parameters exist
    """

    query_set = Parameter.objects.filter(
        televir_project=project.pk, televir_project_sample=None
    )
    if query_set.count() == 0:
        return False
    return True


def duplicate_metagenomics_software(project: Projects, sample: PIProject_Sample):
    owner = project.owner
    project_exists = check_project_params_exist(project)
    if project_exists:
        project_call = project
        types_of_use = Software.TELEVIR_PROJECT_TYPES
    else:
        project_call = None
        types_of_use = Software.TELEVIR_GLOBAL_TYPES
    project_call = project if check_project_params_exist(project) else None
    query_set = Software.objects.filter(
        owner=owner,
        type_of_use__in=types_of_use,
        type_of_software__in=[
            Software.TYPE_SOFTWARE,
            Software.TYPE_INSAFLU_PARAMETER,
        ],
        is_obsolete=False,
        technology__name=project.technology,
        parameter__televir_project=project_call,
        parameter__televir_project_sample=None,
        pipeline_step__name__in=[
            ConstantsSettings.PIPELINE_NAME_metagenomics_combine,
            ConstantsSettings.PIPELINE_NAME_remapping,
            ConstantsSettings.PIPELINE_NAME_remap_filtering,
            ConstantsSettings.PIPELINE_NAME_reporting,
        ],
    )
    project = Projects.objects.get(pk=project.pk)
    for software in query_set:
        software_parameters = Parameter.objects.filter(
            software=software,
        )

        software.pk = None
        if software.type_of_use == Software.TYPE_OF_USE_televir_global:
            software.type_of_use = Software.TYPE_OF_USE_televir_project
        else:
            software.type_of_use = Software.TYPE_OF_USE_televir_project_settings

        try:
            Software.objects.get(
                name=software.name,
                type_of_use=software.type_of_use,
                parameter__televir_project=project,
                parameter__televir_project_sample=sample,
                pipeline_step=software.pipeline_step,
            )

        except Software.MultipleObjectsReturned:
            pass

        except Software.DoesNotExist:
            software.save()
            for parameter in software_parameters:
                parameter.pk = None
                parameter.software = software
                parameter.televir_project = project
                parameter.televir_project_sample = sample
                parameter.save()


def set_control_reports(project_pk: int):
    """
    set control reports
    """

    try:
        project = Projects.objects.get(pk=project_pk)

        control_samples = PIProject_Sample.objects.filter(
            project=project, is_control=True
        )

        control_reports = FinalReport.objects.filter(
            sample__in=control_samples
        ).distinct("taxid")

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
    report_dict: Dict[str, List[FinalReport]]
    excluded_dict: Dict[str, List[FinalReport]]
    error_rate_available: bool

    def __init__(
        self,
        reports: List[FinalReport],
        report_layout_params: LayoutParams,
        force=False,
        level=0,
    ):
        self.reports = reports
        self.max_error_rate = 0
        self.max_quality_avg = 1
        self.max_mapped_prop = 0
        self.max_coverage = 0
        self.max_windows_covered = 1
        self.error_rate_available = self.assess_error_rate_available()
        self.quality_avg_available = self.assess_quality_avg_available()
        self.assess_max_mapped_prop()
        self.assess_max_coverage()
        self.reference_clade = self.generate_reference_clade(report_layout_params)
        self.report_dict = {
            report.accid: report
            for report in reports
            if self.retrieved_mapped_subset(report)
        }
        self.excluded_dict = {
            report.accid: report
            for report in reports
            if not self.retrieved_mapped_subset(report)
        }
        self.metadata_df = self.prep_metadata_df()

        self.fasta_files = self.metadata_df.file.tolist()
        self.level = level

        self.model = self.set_level(reports, level)
        self.run = self.infer_run()

        self.logger = logging.getLogger(__name__)
        self.logger.info("ReportSorter: {}".format(self.run))
        # set logger level
        self.logger.setLevel(logging.INFO)

        if self.model is not None:
            self.media_dir = self.infer_media_dir()

            self.all_clades_df_path = os.path.join(
                self.media_dir, self.all_clade_filename
            )

            self.analysis_df_path = os.path.join(
                self.media_dir,
                self.analysis_filename.format(
                    report_layout_params.shared_proportion_threshold
                ),
            )
            self.force = force
            self.overlap_manager = ReadOverlapManager(
                self.metadata_df,
                self.reference_clade,
                self.media_dir,
                str(self.model.pk),
            )
            self.tree_plot_path = self.overlap_manager.tree_plot_path
            # remove everything before media dir
            self.tree_plot_exists = os.path.exists(self.tree_plot_path)
            self.tree_plot_path = "/media/" + self.tree_plot_path.split("media/")[-1]

            self.overlap_heatmap_exists = os.path.exists(
                self.overlap_manager.overlap_matrix_plot_path
            )
            self.overlap_heatmap_path = (
                "/media/"
                + self.overlap_manager.overlap_matrix_plot_path.split("media/")[-1]
            )
            self.overlap_pca_exists = os.path.exists(
                self.overlap_manager.overlap_pca_plot_path
            )
            self.overlap_pca_path = (
                "/media/"
                + self.overlap_manager.overlap_pca_plot_path.split("media/")[-1]
            )

        else:
            self.media_dir = None
            self.analysis_df_path = None
            self.all_clades_df_path = None
            self.force = False
            self.tree_plot_exists = False
            self.tree_plot_path = None
            self.overlap_heatmap_exists = False
            self.overlap_heatmap_path = None

    def update_max_error_rate(self, report: FinalReport):
        """
        update max error rate"""
        if report.error_rate is None:
            return

        if report.error_rate > self.max_error_rate:
            self.max_error_rate = report.error_rate

    def update_max_quality_avg(self, report: FinalReport):
        """
        update max quality avg"""
        if report.quality_avg is None:
            return

        if report.quality_avg > self.max_quality_avg:
            self.max_quality_avg = report.quality_avg

    def update_max_mapped_prop(self, report: FinalReport):
        """
        update max quality avg"""
        if report.mapped_proportion is None:
            return

        if report.mapped_proportion > self.max_mapped_prop:
            self.max_mapped_prop = report.mapped_proportion

    def update_max_coverage(self, report: FinalReport):
        """
        update max quality avg"""
        if report.coverage is None:
            return

        if report.coverage > self.max_coverage:
            self.max_coverage = report.coverage

    def update_max_windows_covered(self, report: FinalReport):
        """
        update max quality avg"""
        if report.windows_covered is None:
            return

        windows_covered = report.windows_covered
        if "/" in report.windows_covered:
            windows_covered = report.windows_covered.split("/")
            try:
                windows_covered = int(windows_covered[0]) / int(windows_covered[1])
            except Exception as e:
                print(e)
                windows_covered = 0
        else:
            windows_covered = int(windows_covered)

        if windows_covered > self.max_windows_covered:
            self.max_windows_covered = windows_covered

    def assess_max_mapped_prop(self):
        for report in self.reports:
            self.update_max_mapped_prop(report)

    def assess_max_coverage(self):
        for report in self.reports:
            self.update_max_coverage(report)

    def assess_error_rate_available(self):
        if not self.reports:
            return False
        for report in self.reports:
            if report.error_rate is None:
                return False
            self.update_max_error_rate(report)

        return True

    def assess_quality_avg_available(self):
        if not self.reports:
            return False

        for report in self.reports:
            if report.quality_avg is None:
                return False
            self.update_max_quality_avg(report)

        return True

    def set_level(self, final_report_list: List[FinalReport], level):
        if len(final_report_list) == 0:
            return None
        final_report = final_report_list[0]
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
            total_proportion=0,
            group_counts=0,
            shared_proportion_std=layout_params.shared_proportion_threshold,
            shared_proportion_min=layout_params.shared_proportion_threshold,
            shared_proportion_max=layout_params.shared_proportion_threshold,
            overlap_df=pd.DataFrame(),
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

        rundir = infer_run_media_dir(self.run)
        sample_dir = os.path.dirname(rundir)

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

    def read_overlap_analysis(self, force: bool = False):
        """
        Return read overlap analysis as dataframe
        columns: leaf (accid), clade, read_count, group_count
        """
        ### time operations
        self.logger.info("generating tree")

        ### generate tree
        # njtree = overlap_manager.generate_tree()

        ### inner node to leaf dict
        # tree_manager = PhyloTreeManager(njtree)
        # inner_node_leaf_dict = tree_manager.clades_get_leaves_clades()

        ### get statistics
        # statistics_dict_all = overlap_manager.get_node_statistics(force=force)
        # selected_clades = overlap_manager.filter_clades(statistics_dict_all)
        # leaf_clades = overlap_manager.tree_manager.leaf_clades_clean(selected_clades)
        # clades = overlap_manager.leaf_clades_to_pandas(leaf_clades, statistics_dict_all)
        clades = self.overlap_manager.get_leaf_clades(force=force)
        print(clades)
        self.update_report_excluded_dicts(self.overlap_manager)

        return clades

    def update_report_excluded_dicts(self, overlap_manager: ReadOverlapManager):
        new_report_dict = {}

        for accid, report in self.report_dict.items():
            if accid in overlap_manager.excluded_leaves:
                self.excluded_dict[accid] = report
            else:
                new_report_dict[accid] = report

        self.report_dict = new_report_dict

    def check_all_accids_analyzed(self, df: pd.DataFrame):
        """
        Return True if all accids have been analyzed
        """
        for accid in self.report_dict:
            if accid not in df.leaf.tolist():
                return False

        return True

    def check_parsed(self):
        overlap_manager = ReadOverlapManager(
            self.metadata_df,
            self.reference_clade,
            self.media_dir,
            str(self.model.pk),
        )

        if not os.path.exists(overlap_manager.distance_matrix_path):
            return False

        if not os.path.exists(overlap_manager.clade_statistics_path):
            return False

        if not os.path.exists(overlap_manager.clade_statistics_path):
            return False

        if not overlap_manager.all_accs_analyzed():
            return False

        return True

    def check_analysis_exists(self):
        """
        Return True if analysis exists
        """
        if not os.path.exists(self.analysis_df_path):
            return False

        return True

    def check_analyzed(self):
        """
        Return True if all accids have been analyzed
        """

        if not self.check_analysis_exists():
            return False

        if not self.check_parsed():
            return False

        return True

    def sort_reports_save(self, force=False):
        """
        Return sorted reports
        """
        if self.model is not None:
            overlap_analysis = self.read_overlap_analysis(force=True)
            overlap_analysis.to_csv(self.analysis_df_path, sep="\t", index=False)

            overlap_groups = list(overlap_analysis.groupby(["total_counts", "clade"]))[
                ::-1
            ]

            clades_to_keep = []

            for group in overlap_groups:
                group_df = group[1]

                group_accids = group_df.leaf.tolist()
                group_accids = [x for x in group_accids if x in self.report_dict]
                group_list = [self.report_dict[accid] for accid in group_accids]
                # sort by coverage
                group_list.sort(key=lambda x: x.coverage, reverse=True)
                name = group_df.clade.iloc[0]
                if len(group_list):
                    clades_to_keep.append(name)
                    if len(group_list) > 1:
                        pairwise_shared_within_clade = (
                            self.overlap_manager.within_clade_shared_reads(clade=name)
                        )
                        self.get_private_reads_no_duplicates(
                            pairwise_shared_within_clade
                        )

                        self.overlap_manager.plot_pairwise_shared_clade_reads(
                            pairwise_shared_within_clade, subplot=True, clade_str=name
                        )

            pairwise_shared_among_clade = (
                self.overlap_manager.between_clade_shared_reads(
                    clades_filter=clades_to_keep
                )
            )
            self.overlap_manager.plot_pairwise_shared_clade_reads(
                pairwise_shared_among_clade
            )

            self.overlap_manager.plot_pca_full()

    def get_private_reads_no_duplicates(
        self, pairwise_shared_among_clade: pd.DataFrame
    ):
        """"""
        accid_df = pd.read_csv(self.overlap_manager.accid_statistics_path, sep="\t")

        duplicate_groups = self.overlap_manager.duplicate_groups_from_dataframe(
            pairwise_shared_among_clade
        )
        print("duplicate groups")
        print(duplicate_groups)
        accid_df["private_reads"] = 0

        for duplicate_group in duplicate_groups:
            group_private_counts = self.overlap_manager.get_accession_private_counts(
                duplicate_group
            )

            accid_df.loc[
                accid_df.accid.isin(duplicate_group), "private_reads"
            ] = group_private_counts

        print(accid_df)
        accid_df.to_csv(
            self.overlap_manager.accid_statistics_path, sep="\t", index=False
        )

    def wrap_report(self, report: FinalReport) -> FinalReportWrapper:
        return FinalReportWrapper(report)

    def wrap_group_reports(self, report_group: FinalReportGroup) -> FinalReportGroup:
        report_group.group_list = [
            self.wrap_report(report) for report in report_group.group_list
        ]
        return report_group

    def wrap_group_list_reports(
        self, report_groups: List[FinalReportGroup]
    ) -> List[FinalReportGroup]:
        return [self.wrap_group_reports(report_group) for report_group in report_groups]

    def get_reports_private_reads(
        self, report_groups: List[FinalReportGroup]
    ) -> List[FinalReportGroup]:
        """
        Update reports with private reads
        """
        accid_df = pd.read_csv(self.overlap_manager.accid_statistics_path, sep="\t")
        if "private_reads" not in accid_df.columns:
            report_groups = self.wrap_group_list_reports(report_groups)
            return report_groups

        for report_group in report_groups:
            for report in report_group.group_list:
                if report.accid not in accid_df.accid.tolist():
                    continue
                report = self.wrap_report(report)
                accid_df_accid = accid_df[accid_df.accid == report.accid]
                private_reads = accid_df_accid.private_reads.iloc[0]
                report.update_private_reads(private_reads)

            report_group.update_max_private_reads()

        return report_groups

    def sort_group_by_private_reads(self, group: FinalReportGroup) -> FinalReportGroup:
        """
        sort group by private reads
        """
        group.group_list.sort(key=lambda x: x.private_reads, reverse=True)

        return group

    def sort_group_list_reports(
        self, report_groups: List[FinalReportGroup]
    ) -> List[FinalReportGroup]:
        """
        sort group list reports
        """
        return [
            self.sort_group_by_private_reads(report_group)
            for report_group in report_groups
        ]

    def get_sorted_reports(self) -> List[FinalReportGroup]:
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
            name = group_df.clade.iloc[0]
            group_heatmap = self.overlap_manager.get_media_path_heatmap_clade(name)
            group_heatmap_exists = os.path.exists(group_heatmap)

            if group_heatmap_exists:
                group_heatmap = "/media/" + group_heatmap.split("media/")[-1]
            if len(group_list):
                report_group = FinalReportGroup(
                    name=name,
                    total_counts=group_df.total_counts.iloc[0],
                    shared_proportion=group_df.shared_proportion.iloc[0],
                    private_proportion=group_df.private_proportion.iloc[0],
                    group_list=group_list,
                    heatmap_path=group_heatmap,
                    heatmap_exists=group_heatmap_exists,
                )
                sorted_reports.append(report_group)

        # sort groups by max coverage among group

        def get_private_proportion(group: FinalReportGroup):
            return group.private_proportion

        sorted_groups = sorted(sorted_reports, key=get_private_proportion, reverse=True)
        sorted_groups = self.get_reports_private_reads(sorted_groups)
        sorted_groups = self.sort_group_list_reports(sorted_groups)

        return sorted_groups

    def get_reports(self) -> List[FinalReportGroup]:
        """
        Return sorted reports
        """

        if self.model is None:
            return []

        if self.metadata_df.empty:
            return []

        if not self.check_analyzed():
            return [
                FinalReportGroup(
                    name="Full report, no overlap analysis",
                    total_counts=0,
                    shared_proportion=0,
                    private_proportion=0,
                    group_list=self.reports,
                )
            ]

        return self.get_sorted_reports()

    def get_reports_compound(self) -> List[FinalReportGroup]:
        reports = self.get_reports()

        if len(reports) == 0:
            return []

        for report_groups in reports:
            report_groups.group_list = [
                FinalReportCompound(report) for report in report_groups.group_list
            ]

        return reports

    def check_excluded_exist(self) -> bool:
        """return True if there are excluded reports"""

        if len(self.excluded_dict) > 0:
            return True

        return False

    def get_reports_empty(self) -> FinalReportGroup:
        """return reports in excluded report_dict"""

        reports = list(self.excluded_dict.values())

        report_group = FinalReportGroup(
            name="Excluded",
            total_counts=0,
            shared_proportion=0,
            private_proportion=0,
            group_list=reports,
        )

        report_group = self.wrap_group_reports(report_group)

        return report_group


def calculate_reports_overlaps(sample: PIProject_Sample, force=False):
    """
    calculate reports overlaps
    """
    final_reports = FinalReport.objects.filter(sample=sample)
    report_layout_params = TelevirParameters.get_report_layout_params(
        project_pk=sample.project.pk
    )
    report_sorter = ReportSorter(final_reports, report_layout_params)
    report_sorter.sort_reports_save()


def final_report_best_cov_by_accid(reports: QuerySet) -> QuerySet:
    """
    get the best coverage report for each accid
    """

    pk_to_keep = []
    for accid in set(reports.values_list("accid", flat=True)):
        best_report = reports.filter(accid=accid).order_by("-coverage").first()
        pk_to_keep.append(best_report.pk)

    return reports.filter(pk__in=pk_to_keep)
