import json
import logging
import os
from typing import Dict, List, Optional, Union

import pandas as pd
from django.db import transaction
from django.db.models import Q
from django.db.models.query import QuerySet
from django.urls import reverse
from django.utils.safestring import mark_safe

from constants.constants import Constants
from fluwebvirus.settings import STATIC_ROOT
from pathogen_identification.constants_settings import \
    ConstantsSettings as PIConstantsSettings
from pathogen_identification.models import (ContigClassification, FinalReport,
                                            ParameterSet, PIProject_Sample,
                                            Projects, RawReference,
                                            RawReferenceCompoundModel,
                                            ReadClassification,
                                            ReferenceMap_Main, ReferencePanel,
                                            ReferenceSourceFileMap,
                                            RunAssembly, RunDetail, RunMain,
                                            SoftwareTree, SoftwareTreeNode)
from pathogen_identification.utilities.clade_objects import Clade
from pathogen_identification.utilities.overlap_manager import \
    ReadOverlapManager
from pathogen_identification.utilities.televir_parameters import (
    LayoutParams, TelevirParameters)
from pathogen_identification.utilities.utilities_general import (
    infer_run_media_dir, merge_classes, simplify_name)
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager
from settings.constants_settings import ConstantsSettings
from settings.models import Parameter, Software


class SampleReferenceManager:
    def __init__(self, sample: PIProject_Sample):
        self.sample = sample
        self.software_tree: SoftwareTree = self.proxy_tree_prepare()
        self.software_tree_node_storage: Optional[SoftwareTreeNode] = (
            self.proxy_leaf_prepare()
        )
        self.prep_storage()

    def add_reference(self, reference: ReferenceSourceFileMap):
        ref_description = reference.description

        if ref_description is None:
            ref_description = reference.accid

        if len(ref_description) > 150:
            ref_description = ref_description[:150]

        new_reference = RawReference(
            run=self.storage_run,
            accid=reference.accid,
            taxid=reference.taxid,
            description=ref_description,
            status=RawReference.STATUS_UNMAPPED,
            counts=0,
            classification_source="none",
        )

        new_reference.save()

    def proxy_tree_prepare(self):
        try:
            software_tree = SoftwareTree.objects.get(
                model=-1,
                version=0,
                technology=self.sample.project.technology,
                owner=self.sample.project.owner,
                project=self.sample.project,
            )
        except SoftwareTree.DoesNotExist:
            software_tree = SoftwareTree.objects.create(
                model=-1,
                version=0,
                technology=self.sample.project.technology,
                owner=self.sample.project.owner,
                project=self.sample.project,
            )
            software_tree.save()

        return software_tree

    def copy_panel(self, panel: ReferencePanel) -> ReferencePanel:
        """
        copy panel
        """
        with transaction.atomic():
            original_pk = panel.pk
            panel.pk = None
            panel.panel_type = ReferencePanel.PANEL_TYPE_COPIED
            panel.save()
            try:
                panel.parent = ReferencePanel.objects.get(pk=original_pk)
                panel.save()

            except ReferencePanel.DoesNotExist:
                pass

        return panel

    def proxy_leaf_prepare(self):
        try:
            software_tree_node = SoftwareTreeNode.objects.get(
                software_tree=self.software_tree,
                name="storage",
            )

        except SoftwareTreeNode.DoesNotExist:
            software_tree_node = SoftwareTreeNode.objects.create(
                index=-1,
                software_tree=self.software_tree,
                name="storage",
                parent=None,
            )
            software_tree_node.save()

        except Exception as e:
            print(e)
            print("multiple objects returned")
            software_tree_node = None

        return software_tree_node

    def proxy_parameter_set_prepare(self):
        try:
            parameter_set_management = ParameterSet.objects.get(
                sample__project=self.sample.project,
                sample=self.sample,
                leaf=self.software_tree_node_storage,
                status=ParameterSet.STATUS_PROXIED,
            )

        except ParameterSet.DoesNotExist:
            parameter_set_management = ParameterSet.objects.create(
                sample=self.sample,
                leaf=self.software_tree_node_storage,
                status=ParameterSet.STATUS_PROXIED,
                project=self.sample.project,
            )
            parameter_set_management.save()

    @property
    def parameter_set_storage(self):
        return ParameterSet.objects.get(
            sample__project=self.sample.project,
            sample=self.sample,
            leaf=self.software_tree_node_storage,
            status=ParameterSet.STATUS_PROXIED,
        )

    def prep_storage(self):
        self.proxy_parameter_set_prepare()

        try:
            storage_run = RunMain.objects.get(
                name="storage",
                run_type=RunMain.RUN_TYPE_STORAGE,
                project=self.sample.project,
                sample=self.sample,
            )
        except RunMain.DoesNotExist:
            storage_run = RunMain.objects.create(
                name="storage",
                sample=self.sample,
                run_type=RunMain.RUN_TYPE_STORAGE,
                project=self.sample.project,
                parameter_set=self.parameter_set_storage,
                host_depletion_performed=False,
                enrichment_performed=False,
            )
            storage_run.save()

    @property
    def storage_run(self):
        return RunMain.objects.get(
            name="storage",
            run_type=RunMain.RUN_TYPE_STORAGE,
            project=self.sample.project,
            sample=self.sample,
        )

    def anchor_parameter_set_leaf(self, leaf: SoftwareTreeNode):
        """
        mapping run from leaf
        """

        try:
            parameter_set = ParameterSet.objects.get(
                sample=self.sample,
                leaf=leaf,
                project=self.sample.project,
            )

            if parameter_set.status == ParameterSet.STATUS_FINISHED:
                parameter_set.status = ParameterSet.STATUS_QUEUED
                parameter_set.save()

        except ParameterSet.DoesNotExist:
            parameter_set = ParameterSet.objects.create(
                sample=self.sample,
                leaf=leaf,
                status=ParameterSet.STATUS_PROXIED,
                project=self.sample.project,
            )
            parameter_set.save()
        except Exception as e:
            print(e)

        return parameter_set

    def create_mapping_run(
        self, leaf: SoftwareTreeNode, run_type: int, panel_pk: Optional[int] = None
    ) -> RunMain:
        parameter_set = self.anchor_parameter_set_leaf(leaf)

        if panel_pk is not None:
            panel = ReferencePanel.objects.get(pk=panel_pk)

        mapping_run = RunMain.objects.create(
            name=leaf.name,
            sample=self.sample,
            run_type=run_type,
            project=self.sample.project,
            parameter_set=parameter_set,
            host_depletion_performed=False,
            enrichment_performed=False,
            status=RunMain.STATUS_PREP,
            panel=panel if panel_pk is not None else None,
        )

        return mapping_run

    def mapping_run_from_leaf(self, leaf: SoftwareTreeNode) -> RunMain:
        """ """
        return self.create_mapping_run(leaf, RunMain.RUN_TYPE_COMBINED_MAPPING)

    def mapping_request_run_from_leaf(self, leaf: SoftwareTreeNode) -> RunMain:
        """ """
        return self.create_mapping_run(leaf, RunMain.RUN_TYPE_MAP_REQUEST)

    def mapping_request_panel_run_from_leaf(
        self, leaf: SoftwareTreeNode, panel_pk: int
    ) -> RunMain:
        """ """
        return self.create_mapping_run(
            leaf, RunMain.RUN_TYPE_PANEL_MAPPING, panel_pk=panel_pk
        )

    def screening_run_from_leaf(self, leaf: SoftwareTreeNode) -> RunMain:
        """ """
        return self.create_mapping_run(leaf, RunMain.RUN_TYPE_SCREENING)


class RunMainWrapper:

    def __init__(self, run: RunMain):

        self.name = f"run {run.parameter_set.leaf.index}"

        self.record = run
        self.user = run.project.owner
        self.project = run.project
        self.sample = run.sample
        self.parameter_set = run.parameter_set
        self.pk = run.pk

        if run.parameter_set is None:
            self.params_df = pd.DataFrame()

        else:
            utils_manager = Utils_Manager()
            self.params_df = utils_manager.get_leaf_parameters(run.parameter_set.leaf)

            self.params_df.drop_duplicates(["module", "software"], inplace=True)
            self.params_df.set_index("module", inplace=True)

        self.run_type = run.run_type
        self.panel = run.panel
        self.created_in = run.created_in
        self.runtime = run.runtime

    @staticmethod
    def capitalize_software(software_name: str) -> str:

        software_name_list = software_name.split(" ")
        software_name_list = [word.capitalize() for word in software_name_list]

        return " ".join(software_name_list)

    def get_pipeline_software(self, pipeline_name: str):

        software_name = "None"

        if pipeline_name in self.params_df.index:

            software_name = self.params_df.loc[pipeline_name, "software"]
            software_name = str(software_name)

        if "(" in software_name:
            software_name = software_name.split("(")[0]

        if "_" in software_name:
            software_name = software_name.split("_")[0]

        return self.capitalize_software(software_name)

    def progress_display(self) -> str:

        if self.user.username == Constants.USER_ANONYMOUS:
            return mark_safe("report")

        run_log = self.run_progess_tracker()

        return mark_safe(run_log)

    def run_progess_tracker(self) -> str:

        finished_preprocessing = self.record.report != "initial"
        finished_assembly = RunAssembly.objects.filter(run=self.record).count() > 0
        finished_classification = (
            ContigClassification.objects.filter(run=self.record).exists()
            and ReadClassification.objects.filter(run=self.record).exists()
        )

        finished_processing = (
            self.record.parameter_set.status == ParameterSet.STATUS_FINISHED
        )
        #

        finished_remapping = self.record.report == "finished"

        report_link = (
            '<a href="'
            + reverse(
                "sample_detail",
                args=[
                    self.record.project.pk,
                    self.record.sample.pk,
                    self.record.pk,
                ],
            )
            + '">'
            + "<i class='fa fa-bar-chart'></i>"
            + "</a>"
        )

        if finished_processing or finished_remapping:

            return report_link

        else:
            runlog = " <a " + 'href="#" >'
            if finished_preprocessing:
                runlog += '<i class="fa fa-check"'
                runlog += 'title="Preprocessing finished"></i>'
            else:
                runlog += '<i class="fa fa-cog"'
                runlog += 'title="Preprocessing running."></i>'

            runlog += "</a>"

            ###

            runlog += " <a " + 'href="#" >'

            if finished_assembly:
                runlog += '<i class="fa fa-check"'
                runlog += 'title="Assembly finished"></i>'
            else:
                runlog += '<i class="fa fa-cog"'
                if finished_preprocessing:
                    runlog += 'title="Assembly running."></i>'
                else:
                    runlog += 'title="Assembly." style="color: gray;"></i>'
            runlog += "</a>"

            ###

            runlog += " <a " + 'href="#" >'

            if finished_classification:
                runlog += '<i class="fa fa-check"'
                runlog += 'title="Classification finished"></i>'
            else:
                runlog += '<i class="fa fa-cog"'
                if finished_assembly:
                    runlog += 'title="Classification running."></i>'
                else:
                    runlog += 'title="Classification." style="color: gray;"></i>'
            runlog += "</a>"

            runlog += " <a " + 'href="#" >'

            runlog += '<i class="fa fa-cog"'
            if finished_classification:
                runlog += 'title="Mapping to references."></i>'
            else:
                runlog += 'title="Validation mapping" style="color: gray;"></i>'
            runlog += "</a>"

            return runlog

        return ""


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
    accid: str
    sample: PIProject_Sample
    coverage: float
    mapped_proportion: float

    media_fields = [
        "covplot",
        "refa_dotplot",
    ]

    # can take either FinalReport or EmptyRemapMain
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
                    if attr in self.media_fields:
                        setattr(self, attr, self.prep_for_static(getattr(report, attr)))
                except Exception as e:
                    raise e

        self.private_reads = 0
        self.control_flag = report.control_flag
        self.control_flag_str = report.control_flag_str
        self.first_in_group = False
        self.row_class_name = "secondary-row"
        self.display = "none"

    @staticmethod
    def prep_for_static(filepath: str) -> str:

        if STATIC_ROOT in filepath:
            return filepath.split(STATIC_ROOT)[-1]

        return filepath

    def update_private_reads(self, private_reads: int):
        self.private_reads = private_reads


class FinalReportCompound:
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
        self.control_flag_str = report.control_flag_str
        # self.control_flag_str = inform_control_flag(report, self.control_flag_str)
        self.private_reads = 0

        self.row_class_name = report.row_class_name
        self.first_in_group = report.first_in_group
        self.display = report.display

    def update_private_reads(self, private_reads: int):
        self.private_reads = private_reads

    def get_identical_reports_ps(self, report: FinalReport) -> str:
        references_found_in = RawReference.objects.filter(
            run__project__pk=report.run.project.pk,
            run__run_type=RunMain.RUN_TYPE_PIPELINE,
            run__sample__pk=report.sample.pk,
            taxid=report.taxid,
        )

        sets = set([r.run.parameter_set.leaf.index for r in references_found_in])

        if len(sets) == 0:
            return "M"

        return ", ".join([str(s) for s in sets])

    def check_data_exists(self, report: FinalReport) -> bool:
        if report.run is None:
            return False
        return report.run.data_deleted == False

    def get_report_rundetail(self, report: FinalReport) -> RunDetail:
        return RunDetail.objects.get(sample=report.sample, run=report.run)

    def get_report_runmain(self, report: FinalReport) -> Optional[RunMain]:
        return report.run


class FinalReportGroup:
    analysis_empty = False

    name: str
    total_counts: str
    private_counts: int
    shared_proportion: float
    private_proportion: float
    group_list: List[FinalReportWrapper]

    def __init__(
        self,
        name: str,
        total_counts: int,
        private_counts: int,
        shared_proportion: float,
        private_proportion: float,
        group_list: List[FinalReportWrapper],
        heatmap_path: str = "",
        heatmap_exists: bool = False,
        private_counts_exist: bool = True,
        analysis_empty=False,
    ):
        self.name = name
        self.total_counts = f"total counts {total_counts}"
        self.private_counts = private_counts
        self.shared_proportion = shared_proportion
        self.private_proportion = round(private_proportion, 2)
        self.group_list = group_list
        self.heatmap_path = heatmap_path
        self.heatmap_exists = heatmap_exists
        self.max_private_reads = 0
        self.max_coverage = 0
        self.private_counts_exist = private_counts_exist
        self.js_heatmap_ready = False
        self.js_heatmap_data = None
        self.analysis_empty = analysis_empty
        self.has_multiple = True if len(group_list) > 1 else False
        self.toggle = (
            "off"
            if shared_proportion
            > PIConstantsSettings.SORT_GROUP_DISPLAY_DEFAULT_THRESHOLD_SHARED
            else "on"
        )

        self.update_max_coverage()

    def reports_have_private_reads(self) -> bool:
        for report in self.group_list:
            if report.private_reads > 0:
                return True
        return False

    def update_max_coverage(self):
        for report in self.group_list:
            if report.coverage > self.max_coverage:
                self.max_coverage = report.coverage

    def update_max_private_reads(self):
        for report in self.group_list:
            if report.private_reads > self.max_private_reads:
                self.max_private_reads = report.private_reads

    def set_unsorted(self):

        for report in self.group_list:
            report.row_class_name = "unsorted-row"
            report.first_in_group = False
            report.display = "table-row"


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
        pipeline_step__name__in=ConstantsSettings.vect_pipeline_televir_metagenomics_for_parameters,
        technology__name=project.technology,
        parameter__televir_project=project_call,
        parameter__televir_project_sample=None,
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
    report_dict: Dict[str, FinalReport]
    excluded_dict: Dict[str, FinalReport]
    error_rate_available: bool

    def __init__(
        self,
        sample: PIProject_Sample,
        reports: List[FinalReport],
        report_layout_params: LayoutParams,
        force=False,
        level=0,
    ):
        self.reports: List[FinalReport] = reports
        self.sample = sample
        self.analysis_empty = False
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

        self.report_dict = {}
        self.excluded_dict = {}

        for report in reports:
            if report.accid is None:
                continue
            if self.retrieved_mapped_subset(report):
                self.report_dict[report.accid] = report
            else:
                self.excluded_dict[report.accid] = report

        self.level = level

        self.reports_availble = len(reports) > 0

        self.media_dir = sample.media_dir

        self.metadata_df = self.prep_metadata_df()
        self.fasta_files = self.metadata_df.file.tolist()

        self.overlap_manager = ReadOverlapManager(
            self.metadata_df,
            self.reference_clade,
            self.media_dir,
            str(self.sample.pk),
            force_tree_rebuild=force,
            max_reads=PIConstantsSettings.MAX_READS_INPUT,
        )

        self.all_clades_df_path = os.path.join(self.media_dir, self.all_clade_filename)

        self.analysis_df_path = os.path.join(
            self.media_dir,
            self.analysis_filename.format(
                report_layout_params.shared_proportion_threshold
            ),
        )

        self.tree_plot_path = self.overlap_manager.tree_plot_path
        #
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
            "/media/" + self.overlap_manager.overlap_pca_plot_path.split("media/")[-1]
        )

        #####################################
        #####################################
        self.logger = logging.getLogger(__name__)
        self.logger.info("ReportSorter: {}".format(self.media_dir))
        self.logger.setLevel(logging.DEBUG)

    def build_tree(self):

        self.overlap_manager.build_tree()

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
            private_counts=0,
            shared_proportion_std=layout_params.shared_proportion_threshold,
            shared_proportion_min=layout_params.shared_proportion_threshold,
            shared_proportion_max=layout_params.shared_proportion_threshold,
            overlap_df=pd.DataFrame(),
        )

        return ref_clade

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

        except Exception as e:
            print(e)
            return None

    def retrieved_mapped_bam(self, report: FinalReport):
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

            if not mapped_ref.bam_file_path:
                return None

            if os.path.exists(mapped_ref.bam_file_path):
                return mapped_ref.bam_file_path

        except ReferenceMap_Main.DoesNotExist:
            return None

        except Exception as e:
            print(e)
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
            mapped_bam = self.retrieved_mapped_bam(report)
            # f"samtools view -F 0x4 {self.read_map_sorted_bam} | cut -f 1 | sort | uniq > {self.mapped_reads_file}"
            if not mapped_subset_r1:
                continue
            filename = mapped_subset_r1
            metadata_dict.append(
                {
                    "file": filename,
                    "description": description,
                    "accid": accid,
                    "bam": mapped_bam,
                }
            )

        metadata_df = pd.DataFrame(metadata_dict)

        if metadata_df.empty:
            metadata_df = pd.DataFrame(columns=["file", "description", "accid", "bam"])

        return metadata_df

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
            str(self.sample.pk),
            max_reads=PIConstantsSettings.MAX_READS_INPUT,
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
        if self.analysis_df_path is None:
            return False

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

    def sort_group_by_private_reads(self, group: FinalReportGroup) -> FinalReportGroup:
        """
        sort group by private reads
        """
        group.group_list.sort(key=lambda x: x.private_reads, reverse=True)

        if len(group.group_list) == 0:
            return group

        group.group_list[0].first_in_group = True
        group.group_list[0].row_class_name = "primary-row"
        group.group_list[0].display = "table-row"

        return group

    def sort_group_by_coverage(self, group: FinalReportGroup) -> FinalReportGroup:
        """
        sort group by private reads
        """
        group.group_list.sort(key=lambda x: x.coverage, reverse=True)

        if len(group.group_list) == 0:
            return group

        group.group_list[0].first_in_group = True
        group.group_list[0].row_class_name = "primary-row"
        group.group_list[0].display = "table-row"

        return group

    def sort_group_list_reports(
        self, report_groups: List[FinalReportGroup]
    ) -> List[FinalReportGroup]:
        """
        sort group list reports
        """
        return [
            self.sort_group_by_coverage(report_group) for report_group in report_groups
        ]

    def read_overlap_analysis(self, force: bool = False):
        """
        Return read overlap analysis as dataframe
        columns: leaf (accid), clade, read_count, group_count
        """
        ### time operations
        self.logger.info("generating tree")

        clades = self.overlap_manager.get_leaf_clades(force=force)

        self.update_report_excluded_dicts(self.overlap_manager)

        return clades

    def sort_reports_save(self, force=False):
        """
        Return sorted reports
        """
        if self.reports_availble is None:
            return self.return_no_analysis()

        try:
            assert self.overlap_manager.tree_manager is not None
        except AttributeError as e:
            print(e)
            return self.return_no_analysis()

        overlap_analysis = self.read_overlap_analysis(force=True)
        self.overlap_manager.plot_pca_full(overlap_analysis)
        overlap_analysis.to_csv(self.analysis_df_path, sep="\t", index=False)

        self.overlap_manager.get_private_reads_no_duplicates()

        overlap_groups = list(overlap_analysis.groupby(["total_counts", "clade"]))[::-1]

        clades_to_keep = []

        ## plot pairwise shared reads
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

                    self.overlap_manager.plot_pairwise_shared_clade_reads(
                        pairwise_shared_within_clade, subplot=True, clade_str=name
                    )

        pairwise_shared_among_clade = self.overlap_manager.between_clade_shared_reads(
            clades_filter=clades_to_keep
        )

        self.overlap_manager.plot_pairwise_shared_clade_reads(
            pairwise_shared_among_clade
        )

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
            new_group_list = []
            for report_original in report_group.group_list:
                report = self.wrap_report(report_original)
                if report.accid not in accid_df.accid.tolist():
                    new_group_list.append(report)
                    continue

                accid_df_accid = accid_df[accid_df.accid == report.accid]
                private_reads = accid_df_accid.private_reads.iloc[0]
                report.update_private_reads(private_reads)
                new_group_list.append(report)
            report_group.group_list = new_group_list
            report_group.update_max_private_reads()

        return report_groups

    def get_sorted_reports(self) -> List[FinalReportGroup]:

        if os.path.exists(self.analysis_df_path):
            overlap_analysis = pd.read_csv(self.analysis_df_path, sep="\t")
        else:
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

                private_counts_present = "private_counts" in group_df.columns
                private_counts = 0
                if private_counts_present:
                    private_counts = group_df.private_counts.iloc[0]

                report_group = FinalReportGroup(
                    name=name,
                    total_counts=group_df.total_counts.iloc[0],
                    private_counts=private_counts,
                    shared_proportion=group_df.shared_proportion.iloc[0],
                    private_proportion=group_df.private_proportion.iloc[0],
                    group_list=group_list,
                    heatmap_path=group_heatmap,
                    heatmap_exists=group_heatmap_exists,
                    private_counts_exist=private_counts_present,
                )
                sorted_reports.append(report_group)

        # sort groups by max coverage among group

        if len(sorted_reports) == 0:
            return self.return_no_analysis()

        def get_private_proportion(group: FinalReportGroup):
            return group.private_proportion

        def get_group_max_coverage(group: FinalReportGroup):
            return group.max_coverage

        sorted_groups: List[FinalReportGroup] = sorted(
            sorted_reports, key=get_group_max_coverage, reverse=True
        )

        sorted_groups = self.get_reports_private_reads(sorted_groups)
        sorted_groups = self.sort_group_list_reports(sorted_groups)
        sorted_groups = self.prep_heatmap_data_several(sorted_groups)

        return sorted_groups

    def read_shared_matrix(self) -> pd.DataFrame:
        """
        read accession shared reads matrix
        """
        try:
            distance_matrix = pd.read_csv(
                self.overlap_manager.shared_prop_matrix_path, index_col=0
            )
            return distance_matrix
        except Exception as e:
            return pd.DataFrame()

    def read_clade_shared_matrix(self):
        """
        read clade shared matrix
        """
        try:
            clade_shared_matrix = pd.read_csv(
                self.overlap_manager.clade_shared_prop_matrix_path, index_col=0
            )

            # fill diagonal with 1s
            for i in range(clade_shared_matrix.shape[0]):
                clade_shared_matrix.iloc[i, i] = 1

            return clade_shared_matrix

        except Exception as e:
            return None

    def prep_heatmap_data_within_clade(
        self, report_group: FinalReportGroup, distance_matrix: pd.DataFrame
    ):
        """
        prepare heatmap data to be used to create javascript heatmap"""

        group_members = report_group.group_list
        group_members_accids = [member.accid for member in group_members]

        distance_matrix = distance_matrix[
            distance_matrix.index.isin(group_members_accids)
        ]
        distance_matrix = distance_matrix[group_members_accids]

        json_data = self.prep_heatmap_data(distance_matrix)

        return json_data

    def prep_heatmap_data(self, distance_matrix: pd.DataFrame):

        distance_matrix = distance_matrix.fillna(0)

        # Convert the DataFrame to JSON
        json_data = []

        # Calculate row sums
        row_sums = distance_matrix.sum(axis=1)

        # Sort rows by row sums in descending order
        sorted_row_index = row_sums.sort_values(ascending=False).index

        # Sort rows
        distance_matrix = distance_matrix.reindex(sorted_row_index)

        # Sort columns using the same order as rows
        distance_matrix = distance_matrix.reindex(sorted_row_index, axis=1)

        for ix, row in distance_matrix.iterrows():
            for col, value in row.items():
                json_data.append({"x": ix, "y": col, "value": value})

        json_data = json.dumps(json_data)

        return json_data

    def prep_heatmap_data_several(self, report_groups: List[FinalReportGroup]):
        """
        prepare heatmap data to be used to create javascript heatmap
        """
        distance_matrix = self.read_shared_matrix()

        if distance_matrix.empty is True:
            return report_groups

        for report_group in report_groups:
            json_data = self.prep_heatmap_data_within_clade(
                report_group, distance_matrix
            )
            report_group.js_heatmap_data = json_data
            report_group.js_heatmap_ready = True

        return report_groups

    def clade_heatmap_json(self, to_keep=Optional[List[str]]):

        distance_matrix = self.read_clade_shared_matrix()

        if distance_matrix is None:
            return None

        try:
            distance_matrix = distance_matrix.loc[to_keep, to_keep]
        except KeyError as e:
            print(e)
            return None

        return self.prep_heatmap_data(distance_matrix)

    def return_no_analysis(self) -> List[FinalReportGroup]:
        report_group = FinalReportGroup(
            name="Full report, no overlap analysis",
            total_counts=0,
            private_counts=0,
            shared_proportion=0,
            private_proportion=0,
            group_list=[self.wrap_report(report) for report in self.reports],
            analysis_empty=True,
        )
        self.analysis_empty = True
        report_group = self.wrap_group_reports(report_group)
        report_group.set_unsorted()
        return [report_group]

    def get_reports(self) -> List[FinalReportGroup]:
        """
        Return sorted reports
        """
        if self.reports_availble is None:
            return self.return_no_analysis()

        if self.metadata_df.empty:
            return self.return_no_analysis()

        if not self.check_analyzed():
            report_group = FinalReportGroup(
                name="Full report, no overlap analysis",
                total_counts=0,
                private_counts=0,
                shared_proportion=0,
                private_proportion=0,
                group_list=[self.wrap_report(report) for report in self.reports],
                analysis_empty=True,
            )

            self.analysis_empty = True
            report_group = self.wrap_group_reports(report_group)
            report_group.set_unsorted()
            return [report_group]

        return self.get_sorted_reports()

    def get_reports_compound(self) -> List[FinalReportGroup]:
        reports = self.get_reports()

        if len(reports) == 0:
            return []

        for report_group in reports:
            new_list = []
            for wapped_report in report_group.group_list:
                report_compound = FinalReportCompound(wapped_report)
                report_compound.update_private_reads(wapped_report.private_reads)
                new_list.append(report_compound)

            report_group.group_list = new_list

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
            private_counts=0,
            shared_proportion=0,
            private_proportion=0,
            group_list=reports,
        )

        report_group = self.wrap_group_reports(report_group)
        report_group.set_unsorted()

        return report_group


def calculate_reports_overlaps(sample: PIProject_Sample, force=False):
    """
    calculate reports overlaps
    """
    final_reports = FinalReport.objects.filter(sample=sample)
    report_layout_params = TelevirParameters.get_report_layout_params(
        project_pk=sample.project.pk
    )
    report_sorter = ReportSorter(
        sample, final_reports, report_layout_params, force=force
    )

    report_sorter.build_tree()

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


class ReferenceManager:
    def __init__(self):
        pass


class RawReferenceCompound:
    def __init__(self, raw_reference: RawReference):

        self.sample_id = raw_reference.run.sample.pk
        self.selected_mapped_pk = raw_reference.id
        self.taxid = raw_reference.taxid
        self.accid = raw_reference.accid
        self.description = raw_reference.description
        self.family = []
        self.runs = []
        self.manual_insert = False
        self.mapped_final_report: Optional[FinalReport] = None
        self.mapped_raw_reference: Optional[RawReference] = None
        self.standard_score = 0
        self.ensemble_ranking = None
        self.global_ranking = None

        if raw_reference.run.sample is not None:
            self.find_across_sample(raw_reference.run.sample)
            self.find_mapped(raw_reference.run.sample)

        self.determine_runs()

    def find_across_sample(self, sample: PIProject_Sample):
        """
        find across sample
        """
        across_sample = RawReference.objects.filter(
            taxid=self.taxid,
            accid=self.accid,
            description=self.description,
            run__sample=sample,
        ).values_list("pk", flat=True)

        for raw_reference in across_sample:
            self.family.append(raw_reference)

    def determine_runs(self):
        """
        determine runs
        """
        family_refs = (
            RawReference.objects.filter(
                pk__in=self.family, run__run_type=RunMain.RUN_TYPE_PIPELINE
            )
            .exclude(run=None)
            .distinct("run")
        )
        self.runs = [x.run for x in family_refs]

    def determine_manual_insert(self):
        """
        determine if manual insert
        """
        self.manual_insert = (
            RawReference.objects.filter(pk__in=self.family)
            .exclude(run__run_type=RunMain.RUN_TYPE_PIPELINE)
            .exists()
        )

    def find_mapped(self, sample: PIProject_Sample):
        """
        find mapped, get the first, sorted by coverage
        """

        self.mapped_final_report = (
            FinalReport.objects.filter(
                taxid=self.taxid,
                accid=self.accid,
                sample=sample,
            )
            .order_by("-coverage")
            .first()
        )

        self.mapped_raw_reference = (
            RawReference.objects.filter(
                taxid=self.taxid,
                accid=self.accid,
                run__sample=sample,
                status=RawReference.STATUS_MAPPED,
            )
            .exclude(
                run__run_type__in=[RunMain.RUN_TYPE_SCREENING, RunMain.RUN_TYPE_STORAGE]
            )
            .first()
        )

    @property
    def mapped_html(self):

        if self.mapped_final_report is None and self.mapped_raw_reference is None:
            return mark_safe(
                '<a><i class="fa fa-times" title="unmapped"></i> Unmapped</a>'
            )

        if self.mapped_final_report is not None:
            run = self.mapped_final_report.run

            return mark_safe(
                '<a href="'
                + reverse(
                    "sample_detail",
                    args=[run.sample.project.pk, run.sample.pk, run.pk],
                )
                + '" title="workflow link">'
                + '<i class="fa fa-check-circle"></i>'
                + " Mapped"
                + "</a>"
            )

        elif self.mapped_raw_reference is not None:
            run = self.mapped_raw_reference.run
            return mark_safe(
                '<a href="'
                + reverse(
                    "sample_detail",
                    args=[run.sample.project.pk, run.sample.pk, run.pk],
                )
                + '" title="workflow link">'
                + "<i class='fa fa-circle-o'></i>"
                + " Mapped, 0 reads"
                + "</a>"
            )

    @property
    def run_count(self):
        return len(self.runs)

    @property
    def runs_str(self):
        return ", ".join([str(r.parameter_set.leaf.index) for r in self.runs])


class RawReferenceUtils:
    def __init__(
        self,
        sample: Optional[PIProject_Sample] = None,
        project: Optional[Projects] = None,
    ):
        self.sample_registered = sample
        self.project_registered = project
        self.runs_found = 0
        self.list_tables: List[pd.DataFrame] = []
        self.merged_table: pd.DataFrame = pd.DataFrame(
            columns=["read_counts", "contig_counts", "taxid", "accid", "description"]
        )

    def references_table_from_query(
        self, references: Union[QuerySet, List[RawReference]]
    ) -> pd.DataFrame:
        table = []
        for ref in references:
            table.append(
                {
                    "taxid": ref.taxid,
                    "accid": ref.accid,
                    "description": ref.description,
                    "counts_str": ref.counts,
                    "read_counts": ref.read_counts,
                    "contig_counts": ref.contig_counts,
                }
            )

        if len(table) == 0:
            return pd.DataFrame(
                columns=[
                    "taxid",
                    "accid",
                    "description",
                    "counts_str",
                    "read_counts",
                    "contig_counts",
                ]
            )

        references_table = pd.DataFrame(table)

        references_table["read_counts"] = references_table["read_counts"].astype(float)
        references_table["contig_counts"] = references_table["contig_counts"].astype(
            float
        )

        # references_table = references_table[references_table["read_counts"] > 1]
        references_table = references_table[references_table["accid"] != "-"]
        references_table = references_table[references_table["accid"] != ""]
        references_table = references_table.sort_values(
            ["contig_counts", "read_counts"], ascending=[False, False]
        )

        references_table["sort_rank"] = range(1, references_table.shape[0] + 1)

        return references_table

    def run_references_standard_score_reads(
        self,
        references_table: pd.DataFrame,
    ) -> pd.DataFrame:

        if references_table.shape[0] == 0:
            return pd.DataFrame(
                columns=list(references_table.columns) + ["read_counts_standard_score"]
            )

        if max(references_table["read_counts"]) == 0:
            references_table["read_counts_standard_score"] = 1
            return references_table

        references_table["read_counts"] = references_table["read_counts"].astype(float)

        references_table["read_counts_standard_score"] = (
            references_table["read_counts"] - references_table["read_counts"].mean()
        ) / references_table["read_counts"].std()

        references_table["read_counts_standard_score"] = (
            references_table["read_counts_standard_score"]
            - references_table["read_counts_standard_score"].min()
        ) / (
            references_table["read_counts_standard_score"].max()
            - references_table["read_counts_standard_score"].min()
        )

        return references_table

    def run_references_standard_score_contigs(
        self,
        references_table: pd.DataFrame,
    ) -> pd.DataFrame:
        # references = RawReference.objects.filter(run=run)

        # references_table = references_table_from_query(references)

        if references_table.shape[0] == 0:
            return pd.DataFrame(
                columns=list(references_table.columns)
                + ["contig_counts_standard_score"]
            )

        references_table["contig_counts"] = references_table["contig_counts"].astype(
            float
        )
        if max(references_table["contig_counts"]) == 0:
            references_table["contig_counts_standard_score"] = 1
            return references_table

        references_table["contig_counts_standard_score"] = (
            references_table["contig_counts"] - references_table["contig_counts"].mean()
        ) / references_table["contig_counts"].std()

        references_table["contig_counts_standard_score"] = (
            references_table["contig_counts_standard_score"]
            - references_table["contig_counts_standard_score"].min()
        ) / (
            references_table["contig_counts_standard_score"].max()
            - references_table["contig_counts_standard_score"].min()
        )

        return references_table

    def merge_standard_scores(self, table: pd.DataFrame):
        if table.shape[0] == 0:
            return pd.DataFrame(columns=list(table.columns) + ["standard_score"])
        table["standard_score"] = table[
            "read_counts_standard_score"
        ]  # + table["contig_counts_standard_score"]
        return table

    def run_references_standard_scores(self, table):
        table = self.run_references_standard_score_reads(table)

        table = self.run_references_standard_score_contigs(table)
        table = self.merge_standard_scores(table)
        return table

    def merge_ref_tables_use_standard_score(
        self,
        list_tables: List[pd.DataFrame],
    ) -> pd.DataFrame:
        joint_tables = [
            self.run_references_standard_scores(table) for table in list_tables
        ]

        joint_tables = pd.concat(joint_tables)
        # group tables: average read_counts_standard_score, sum counts, read_counts, contig_counts
        if joint_tables.shape[0] == 0:
            return pd.DataFrame(columns=list(joint_tables.columns))

        joint_tables["standard_score"] = joint_tables["standard_score"].astype(float)
        joint_tables["contig_counts"] = joint_tables["contig_counts"].astype(float)
        joint_tables["read_counts"] = joint_tables["read_counts"].astype(float)
        joint_tables["contig_counts_standard_score"] = joint_tables[
            "contig_counts_standard_score"
        ].astype(float)
        joint_tables["sort_rank"] = joint_tables["sort_rank"].astype(float)

        joint_tables = joint_tables.groupby(["taxid"]).agg(
            {
                "taxid": "first",
                "accid": "first",
                "description": "first",
                "counts_str": "first",
                "standard_score": "mean",
                "contig_counts_standard_score": "mean",
                "read_counts": "sum",
                "contig_counts": "sum",
                "sort_rank": "mean",
            }
        )

        joint_tables = joint_tables.rename(columns={"sort_rank": "ensemble_ranking"})

        #############################################
        proxy_rclass = self.reference_table_renamed(
            joint_tables, {"read_counts": "counts"}
        ).reset_index(drop=True)
        proxy_aclass = self.reference_table_renamed(
            joint_tables, {"contig_counts": "counts"}
        ).reset_index(drop=True)

        targets, _ = merge_classes(
            proxy_rclass, proxy_aclass, maxt=joint_tables.shape[0]
        )

        targets["global_ranking"] = range(1, targets.shape[0] + 1)

        def set_global_ranking_repeat_ranks(
            targets, rank_column="global_ranking", counts_column="counts"
        ):
            """
            Set the global ranking for repeated ranks
            """
            current_counts = None
            current_rank = 0
            for row in targets.iterrows():
                if current_counts != row[1][counts_column]:
                    current_rank += 1
                    current_counts = row[1][counts_column]

                targets.at[row[0], rank_column] = current_rank

            return targets

        targets = set_global_ranking_repeat_ranks(targets, rank_column="global_ranking")

        ####
        joint_tables = joint_tables.reset_index(drop=True)
        targets = targets.reset_index(drop=True)
        joint_tables["taxid"] = joint_tables["taxid"].astype(int)
        targets["taxid"] = targets["taxid"].astype(int)

        joint_tables = joint_tables.merge(
            targets[["taxid", "global_ranking"]], on=["taxid"], how="left"
        )
        ############################################# Reset the index
        joint_tables = joint_tables.reset_index(drop=True)

        joint_tables = joint_tables.sort_values("ensemble_ranking", ascending=True)

        joint_tables = joint_tables.reset_index(drop=True)

        return joint_tables

    def merge_ref_tables_use_ranking(
        self,
        list_tables: List[pd.DataFrame],
    ) -> pd.DataFrame:
        # joint_tables = [
        #    self.run_references_standard_scores(table) for table in list_tables
        # ]

        joint_tables = pd.concat(list_tables)
        # group tables: average read_counts_standard_score, sum counts, read_counts, contig_counts
        if joint_tables.shape[0] == 0:
            return pd.DataFrame(columns=list(joint_tables.columns))

        joint_tables["contig_counts"] = joint_tables["contig_counts"].astype(float)
        joint_tables["read_counts"] = joint_tables["read_counts"].astype(float)

        joint_tables = joint_tables.groupby(["taxid", "accid", "description"]).agg(
            {
                "taxid": "first",
                "accid": "first",
                "description": "first",
                "counts_str": "first",
                "read_counts": "sum",
                "contig_counts": "sum",
                "sort_rank": "mean",
            }
        )
        joint_tables["standard_score"] = 0
        joint_tables = joint_tables.rename(columns={"sort_rank": "ensemble_ranking"})

        proxy_rclass = self.reference_table_renamed(
            joint_tables, {"read_counts": "counts"}
        ).reset_index(drop=True)
        proxy_aclass = self.reference_table_renamed(
            joint_tables, {"contig_counts": "counts"}
        ).reset_index(drop=True)

        targets, _ = merge_classes(
            proxy_rclass, proxy_aclass, maxt=joint_tables.shape[0]
        )

        targets["global_ranking"] = range(1, targets.shape[0] + 1)

        def set_global_ranking_repeat_ranks(
            targets, rank_column="global_ranking", counts_column="counts"
        ):
            """
            Set the global ranking for repeated ranks
            """
            current_counts = None
            current_rank = 0
            for row in targets.iterrows():
                if current_counts != row[1][counts_column]:
                    current_rank += 1
                    current_counts = row[1][counts_column]

                targets.at[row[0], rank_column] = current_rank

            return targets

        targets = set_global_ranking_repeat_ranks(targets, rank_column="global_ranking")

        ####
        joint_tables = joint_tables.reset_index(drop=True)
        targets = targets.reset_index(drop=True)
        joint_tables["taxid"] = joint_tables["taxid"].astype(int)
        targets["taxid"] = targets["taxid"].astype(int)

        joint_tables = joint_tables.merge(
            targets[["taxid", "global_ranking"]], on=["taxid"], how="left"
        )
        ############################################# Reset the index

        joint_tables = joint_tables.sort_values("ensemble_ranking", ascending=True)

        joint_tables = joint_tables.reset_index(drop=True)

        return joint_tables

    def run_references_table(self, run: RunMain) -> pd.DataFrame:
        references = RawReference.objects.filter(run=run)

        references_table = self.references_table_from_query(references)

        return references_table

    def filter_runs(self):
        if self.sample_registered is None and self.project_registered is None:
            raise Exception("No sample or project registered")

        if self.sample_registered is None:
            sample_runs = RunMain.objects.filter(
                sample__project=self.project_registered
            )
        else:
            sample_runs = RunMain.objects.filter(sample=self.sample_registered)

        return sample_runs.exclude(
            run_type__in=[RunMain.RUN_TYPE_SCREENING, RunMain.RUN_TYPE_STORAGE]
        )

    def collect_references_all(self) -> QuerySet:
        sample_runs = self.filter_runs()

        references = RawReference.objects.filter(run__in=sample_runs)

        return references

    def sample_compound_refs_table(self) -> pd.DataFrame:

        compound_refs = RawReferenceCompoundModel.objects.filter(
            sample=self.sample_registered
        )

        compound_refs_table = pd.DataFrame(
            list(
                compound_refs.values(
                    "taxid",
                    "accid",
                    "description",
                    "ensemble_ranking",
                    "standard_score",
                )
            )
        )

        return compound_refs_table

    def sample_reference_tables(
        self, run_pks: Optional[List[int]] = None
    ) -> pd.DataFrame:
        sample_runs = self.filter_runs()
        if run_pks is not None:
            sample_runs = sample_runs.filter(pk__in=run_pks)
        self.runs_found = sample_runs.count()

        run_references_tables = [self.run_references_table(run) for run in sample_runs]

        # register tables
        self.list_tables.extend(run_references_tables)
        #
        if len(self.list_tables):
            run_references_tables = self.merge_ref_tables()
            # replace nan with 0
            run_references_tables = run_references_tables.fillna(0)
            self.merged_table = run_references_tables
        else:
            return pd.DataFrame(
                columns=[
                    "read_counts",
                    "contig_counts",
                    "taxid",
                    "accid",
                    "description",
                ]
            )

    def sample_reference_tables_filter(self, runs_filter: Optional[List[int]] = None):
        """
        Filter the sample reference tables to only include runs in the list
        """
        if runs_filter == []:
            runs_filter = None

        _ = self.sample_reference_tables(runs_filter)

    def compound_reference_update_ranks(self, compound_ref: RawReferenceCompound):
        """
        Update the standard score for a compound reference based on accid"""

        score = self.merged_table[self.merged_table.accid == compound_ref.accid]

        if score.shape[0] > 0:
            # compound_ref.standard_score = score.iloc[0]["standard_score"]
            compound_ref.standard_score = max(score["standard_score"])
            compound_ref.global_ranking = min(score["global_ranking"])
            compound_ref.ensemble_ranking = min(score["ensemble_ranking"])

    def update_scores_compound_references(
        self, compount_refs: List[RawReferenceCompound]
    ):
        """
        Update the standard score for a list of compound references based on accid"""
        for compound_ref in compount_refs:
            self.compound_reference_update_ranks(compound_ref)

    def filter_reference_query_set(
        self, references: QuerySet, query_string: Optional[str] = ""
    ):
        """
        Filter a query set of references by a query string
        """
        if not query_string:
            return references.exclude(accid="-")

        references_select = references.filter(
            Q(description__icontains=query_string)
            | Q(accid__icontains=query_string)
            | Q(taxid__icontains=query_string)
        ).exclude(accid="-")

        return references_select

    def filter_reference_query_set_compound(
        self, references: QuerySet, query_string: Optional[str] = ""
    ):
        """
        Filter a query set of references by a query string
        """
        references_select = self.filter_reference_query_set(references, query_string)

        exclude_refs = []

        for ref in references_select:
            if RawReference.objects.filter(pk=ref.selected_mapped_pk).exists() is False:
                exclude_refs.append(ref.pk)

        references_select = references_select.exclude(pk__in=exclude_refs)

        return references_select

    def query_sample_compound_references_regressive(
        self, query_string: Optional[str] = None
    ) -> List[RawReferenceCompoundModel]:
        """
        Query compound references using sample, project or none, in that order. Filter by query string
        """

        if self.sample_registered is not None:
            query_set = RawReferenceCompoundModel.objects.filter(
                sample=self.sample_registered
            ).order_by("ensemble_ranking")
        elif self.project_registered is not None:
            query_set = RawReferenceCompoundModel.objects.filter(
                sample__project=self.project_registered
            ).order_by("ensemble_ranking")
        else:
            query_set = RawReferenceCompoundModel.objects.none()

        return self.filter_reference_query_set_compound(query_set, query_string)

    def query_sample_references(self, query_string: Optional[str] = "") -> QuerySet:

        if self.sample_registered is not None:

            query_set = (
                RawReference.objects.filter(
                    run__sample__pk=self.sample_registered.pk,
                )
                .exclude(
                    run__run_type__in=[
                        RunMain.RUN_TYPE_STORAGE,
                        RunMain.RUN_TYPE_SCREENING,
                    ],
                    accid="-",
                )
                .distinct("accid")
            )
        elif self.project_registered is not None:
            query_set = (
                RawReference.objects.filter(
                    run__sample__project__pk=self.project_registered.pk,
                )
                .exclude(
                    run__run_type__in=[
                        RunMain.RUN_TYPE_STORAGE,
                        RunMain.RUN_TYPE_SCREENING,
                    ],
                    accid="-",
                )
                .distinct("accid")
            )
        else:
            query_set = RawReference.objects.none()

        return self.filter_reference_query_set(query_set, query_string)

    def register_compound_references(self, compound_refs: List[RawReferenceCompound]):
        """
        Register a list of compound references in the database
        """

        for compound_ref in compound_refs:
            self.register_compound_reference(compound_ref)

    def register_compound_reference(self, compound_ref: RawReferenceCompound):
        """
        Register a compound reference in the database
        """

        try:
            compound_ref_model = RawReferenceCompoundModel.objects.get(
                accid=compound_ref.accid, sample__id=compound_ref.sample_id
            )

            compound_ref_model.standard_score = compound_ref.standard_score
            compound_ref_model.global_ranking = compound_ref.global_ranking
            compound_ref_model.ensemble_ranking = compound_ref.ensemble_ranking
            compound_ref_model.manual_insert = compound_ref.manual_insert
            compound_ref_model.mapped_final_report = compound_ref.mapped_final_report
            compound_ref_model.mapped_raw_reference = compound_ref.mapped_raw_reference
            compound_ref_model.selected_mapped_pk = compound_ref.selected_mapped_pk
            compound_ref_model.run_count = compound_ref.run_count
            compound_ref_model.save()

        except RawReferenceCompoundModel.DoesNotExist:
            sample = PIProject_Sample.objects.get(pk=compound_ref.sample_id)

            description_short = compound_ref.description[:200]
            compound_ref_model = RawReferenceCompoundModel(
                taxid=compound_ref.taxid,
                description=description_short,
                accid=compound_ref.accid,
                sample=sample,
                standard_score=compound_ref.standard_score,
                global_ranking=compound_ref.global_ranking,
                ensemble_ranking=compound_ref.ensemble_ranking,
                manual_insert=compound_ref.manual_insert,
                mapped_final_report=compound_ref.mapped_final_report,
                mapped_raw_reference=compound_ref.mapped_raw_reference,
                selected_mapped_pk=compound_ref.selected_mapped_pk,
                run_count=compound_ref.run_count,
            )
            compound_ref_model.save()

            # for run in compound_ref.runs:
            #    compound_ref_model.runs.add(run)
            compound_ref_model.runs.add(*compound_ref.runs)

            refs = RawReference.objects.filter(pk__in=compound_ref.family)

            ## add family
            compound_ref_model.family.add(*refs)

            # for ref_pk in compound_ref.family:
            #    ref = RawReference.objects.get(pk=ref_pk)
            #    compound_ref_model.family.add(ref)

    def get_classification_runs(self):

        if self.sample_registered is not None:
            classification_runs = RunMain.objects.filter(
                sample=self.sample_registered, run_type=RunMain.RUN_TYPE_PIPELINE
            )

        elif self.project_registered is not None:

            classification_runs = RunMain.objects.filter(
                sample__project=self.project_registered,
                run_type=RunMain.RUN_TYPE_PIPELINE,
            )

        else:
            classification_runs = RunMain.objects.none()

        return classification_runs

    def create_compound(self, raw_references: List[RawReference]):

        raw_reference_compound = [
            RawReferenceCompound(raw_reference) for raw_reference in raw_references
        ]

        self.sample_reference_tables_filter()

        self.update_scores_compound_references(raw_reference_compound)
        self.register_compound_references(raw_reference_compound)

    def retrieve_compound_references(
        self, query_string: Optional[str] = None
    ) -> QuerySet:
        """
        Retrieve compound references for a sample, query, or create them if they do not exist
        """
        compound_refs = self.query_sample_compound_references_regressive(query_string)

        if compound_refs.exists():

            return compound_refs

        compound_refs = self.create_compound_references(query_string=query_string)

        return compound_refs

    def create_compound_references(self, query_string: Optional[str] = None):
        """
        Create compound references for a sample_name, query_string):

        Returns a list of references that match the query string.
        :param query_string:
        :return:

        """
        try:
            references = self.query_sample_references(query_string)
        except Exception as e:
            print(e)

        if references.exists():
            self.create_compound(references)

            compound_refs = self.query_sample_compound_references_regressive(
                query_string
            )

        else:
            compound_refs = RawReferenceCompoundModel.objects.none()

        return compound_refs

    def reference_table_renamed(self, merged_table, rename_dict: dict):
        proxy_ref = merged_table.copy()

        for key, value in rename_dict.items():
            if key in proxy_ref.columns:
                if value in proxy_ref.columns:
                    # remove the column if it exists
                    proxy_ref = proxy_ref.drop(columns=[value])

                proxy_ref = proxy_ref.rename(columns={key: value})

        # proxy_ref = proxy_ref.rename(columns=rename_dict)

        proxy_ref["taxid"] = proxy_ref["taxid"].astype(int)
        proxy_ref["counts"] = proxy_ref["counts"].astype(float).astype(int)
        proxy_ref = proxy_ref[proxy_ref["counts"] > 0]
        proxy_ref = proxy_ref[proxy_ref["taxid"] > 0]
        proxy_ref = proxy_ref[proxy_ref["description"] != "-"]
        proxy_ref = proxy_ref[proxy_ref["accid"] != "-"]

        return proxy_ref

    def merge_ref_tables(self):
        run_references_tables = self.merge_ref_tables_use_ranking(self.list_tables)

        run_references_tables = run_references_tables[run_references_tables.taxid != 0]

        return run_references_tables

    @staticmethod
    def simplify_by_description(df: pd.DataFrame):
        if "description" not in df.columns:
            return df

        df["description_first"] = df["description"].str.split(" ").str[0]

        df = df.sort_values("standard_score", ascending=False)
        df = df.drop_duplicates(subset=["description_first"], keep="first")

        df.drop(columns=["description_first"], inplace=True)

        return df
