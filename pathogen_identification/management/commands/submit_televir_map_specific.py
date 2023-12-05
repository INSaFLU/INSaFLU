import argparse
import logging
import os
import shutil

import pandas as pd
from django.core.management.base import BaseCommand

from constants.constants import Televir_Metadata_Constants as Televir_Metadata
from managing_files.models import ProcessControler
from pathogen_identification.constants_settings import MEDIA_ROOT, ConstantsSettings
from pathogen_identification.install_registry import Params_Illumina, Params_Nanopore
from pathogen_identification.models import (
    FinalReport,
    Projects,
    RawReference,
    RunAssembly,
    RunMain,
    SoftwareTreeNode,
)
from pathogen_identification.modules.metadata_handler import Metadata_handler
from pathogen_identification.modules.object_classes import (
    Read_class,
    Sample_runClass,
    Software_detail,
    SoftwareDetailCompound,
    SoftwareRemap,
)
from pathogen_identification.modules.remap_class import (
    Mapping_Instance,
    Mapping_Manager,
)
from pathogen_identification.utilities.televir_parameters import TelevirParameters
from pathogen_identification.utilities.update_DBs import (
    Update_FinalReport,
    Update_ReferenceMap_Update,
)
from pathogen_identification.utilities.utilities_general import simplify_name_lower
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager
from pathogen_identification.utilities.utilities_views import (
    ReportSorter,
    TelevirParameters,
    recover_assembly_contigs,
)
from settings.constants_settings import ConstantsSettings as CS
from utils.process_SGE import ProcessSGE


class RunEngine:
    remap_manager: Mapping_Manager
    mapping_instance: Mapping_Instance
    metadata_tool: Metadata_handler

    ##  metadata
    sift_query: str
    max_remap: int
    taxid_limit: int

    ## directories.
    root: str

    input_reads_dir: str
    filtered_reads_dir: str
    depleted_reads_dir: str

    log_dir: str

    dir_classification: str = f"classification_reports"
    dir_plots: str = f"plots"
    igv_dir: str = f"igv"

    def __init__(
        self, config: dict, method_args: pd.DataFrame, project_name: str, username: str
    ):
        self.sample_name = config["sample_name"]
        self.type = config["type"]
        self.project_name = project_name
        self.username = username
        self.prefix = "none"
        self.config = config
        self.taxid = config["taxid"]
        self.accid = config["accid"]
        self.threads = config["threads"]
        self.house_cleaning = False
        self.clean = config["clean"]

        self.full_report = os.path.join(
            self.config["directories"][CS.PIPELINE_NAME_remapping], "full_report.csv"
        )

        self.logger_level_main = logging.INFO
        self.logger_level_detail = logging.ERROR
        self.logger = logging.getLogger("main {}".format(self.prefix))
        self.logger.setLevel(self.logger_level_main)

        logFormatter = logging.Formatter(
            fmt="{} {} %(levelname)s :%(message)s".format(
                config["sample_name"], self.prefix
            )
        )

        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        self.logger.addHandler(consoleHandler)
        self.logger.propagate = False

        ######## DIRECTORIES ########

        self.deployment_root_dir = config["deployment_root_dir"]
        self.substructure_dir = config["sub_directory"]
        self.deployment_dir = os.path.join(
            self.deployment_root_dir, self.substructure_dir
        )

        self.media_dir = os.path.join(
            ConstantsSettings.media_directory, self.substructure_dir
        )
        self.static_dir = os.path.join(
            ConstantsSettings.static_directory, self.substructure_dir
        )

        self.media_dir_logdir = os.path.join(
            self.media_dir,
            "logs",
        )
        #####################################

        self.r1 = Read_class(
            config["r1"],
            config["directories"][CS.PIPELINE_NAME_read_quality_analysis],
            config["directories"]["reads_enriched_dir"],
            config["directories"]["reads_depleted_dir"],
            bin=get_bindir_from_binaries(
                config["bin"], CS.PIPELINE_NAME_read_quality_analysis
            ),
        )

        self.r2 = Read_class(
            config["r2"],
            config["directories"][CS.PIPELINE_NAME_read_quality_analysis],
            config["directories"]["reads_enriched_dir"],
            config["directories"]["reads_depleted_dir"],
            bin=get_bindir_from_binaries(
                config["bin"], CS.PIPELINE_NAME_read_quality_analysis
            ),
        )

        self.sample = Sample_runClass(
            self.r1,
            self.r2,
            self.sample_name,
            self.project_name,
            self.username,
            self.config["technology"],
            self.type,
            0,
            ",".join(
                [os.path.basename(self.r1.current), os.path.basename(self.r2.current)]
            ),
            bin=get_bindir_from_binaries(
                config["bin"], CS.PIPELINE_NAME_read_quality_analysis
            ),
            threads=self.threads,
        )

        self.contigs = config["contig_file"]

        ### mapping parameters
        self.min_scaffold_length = config["assembly_contig_min_length"]
        self.minimum_coverage = int(config["minimum_coverage_threshold"])
        self.maximum_coverage = 1000000000

        ### metadata
        remap_params = TelevirParameters.get_remap_software(
            self.username, self.project_name
        )
        self.metadata_tool = Metadata_handler(
            self.username,
            self.config,
            sift_query=config["sift_query"],
            prefix=self.prefix,
            rundir=self.deployment_dir,
        )

        self.max_remap = remap_params.max_accids
        self.taxid_limit = remap_params.max_taxids
        self.remap_params = remap_params

        ### methods
        self.remapping_method = Software_detail(
            CS.PIPELINE_NAME_remapping,
            method_args,
            config,
            self.prefix,
        )

        self.remap_filtering_method = SoftwareDetailCompound(
            CS.PIPELINE_NAME_remap_filtering,
            method_args,
            config,
            self.prefix,
        )

        ###

        self.software_remap = SoftwareRemap(
            self.remapping_method,
            self.remap_filtering_method,
        )

        ###

        self.media_dir_classification = os.path.join(
            self.media_dir,
            self.dir_classification,
        )

        self.static_dir_plots = os.path.join(
            self.substructure_dir,
            self.dir_plots,
        )

        self.media_dir_igv = os.path.join(
            self.static_dir,
            self.igv_dir,
        )

        os.makedirs(
            self.media_dir_classification,
            exist_ok=True,
        )

        os.makedirs(
            os.path.join(ConstantsSettings.static_directory, self.static_dir_plots),
            exist_ok=True,
        )

        os.makedirs(
            self.media_dir_igv,
            exist_ok=True,
        )

        self.filtered_reads_dir = config["directories"][
            CS.PIPELINE_NAME_read_quality_analysis
        ]
        self.log_dir = config["directories"]["log_dir"]

    def export_sequences(self):
        self.sample.export_reads(self.media_dir)

    def generate_targets(self):
        result_df = pd.DataFrame(columns=["qseqid", "taxid"])
        if self.taxid:
            result_df = pd.DataFrame([{"taxid": self.taxid, "qseqid": ""}])
        elif self.accid:
            result_df = pd.DataFrame([{"qseqid": self.accid, "taxid": ""}])

        self.metadata_tool.match_and_select_targets(
            result_df,
            pd.DataFrame(columns=["qseqid", "taxid"]),
            self.max_remap,
            self.taxid_limit,
        )

    def deploy_REMAPPING(self):
        """ """
        self.remap_manager = Mapping_Manager(
            self.metadata_tool.remap_targets,
            self.sample.r1,
            self.sample.r2,
            self.software_remap,
            self.contigs,
            self.type,
            self.prefix,
            self.threads,
            get_bindir_from_binaries(self.config["bin"], CS.PIPELINE_NAME_remapping),
            self.logger_level_detail,
            self.house_cleaning,
            self.remap_params,
            logdir=self.config["directories"]["log_dir"],
        )

        self.logger.info(
            f"{self.prefix} remapping # targets: {len(self.metadata_tool.remap_targets)}"
        )

        print("moving to : ", self.static_dir_plots)
        print("moving to : ", self.media_dir_igv)

        self.remap_manager.run_mappings_move_clean(
            self.static_dir_plots, self.media_dir_igv
        )
        self.remap_manager.export_reference_fastas_if_failed(self.media_dir_igv)
        self.remap_manager.export_mapping_files(self.media_dir_igv)

        self.remap_manager.merge_mapping_reports()
        self.remap_manager.collect_final_report_summary_statistics()

    def run(self):
        self.deploy_REMAPPING()
        print("remap_manager.report")
        self.report = self.remap_manager.report
        self.export_final_reports()

    def export_final_reports(self):
        ### main report
        self.report.to_csv(
            self.full_report,
            index=False,
            sep="\t",
            header=True,
        )


def get_bindir_from_binaries(binaries, key, value: str = ""):
    if value == "":
        try:
            return os.path.join(binaries["ROOT"], binaries[key]["default"], "bin")
        except KeyError:
            return ""
    else:
        try:
            return os.path.join(binaries["ROOT"], binaries[key][value], "bin")
        except KeyError:
            return ""


class Input_Generator:
    method_args: pd.DataFrame

    def __init__(self, reference: RawReference, output_dir: str, threads: int = 4):
        self.utils = Utils_Manager()
        self.reference = reference
        self.install_registry = Televir_Metadata

        self.dir_branch = os.path.join(
            ConstantsSettings.televir_subdirectory,
            f"{reference.run.project.owner.pk}",
            f"{reference.run.project.pk}",
            f"{reference.run.sample.pk}",
            "remapping",
            f"{reference.pk}",
        )

        self.technology = reference.run.project.technology
        self.threads = threads
        self.prefix = reference.accid
        self.project = reference.run.project.name
        self.clean = False
        self.deployment_root_dir = os.path.join(
            output_dir, "remapping", f"ref_{reference.pk}"
        )
        self.dir = os.path.join(self.deployment_root_dir, self.dir_branch)

        os.makedirs(self.dir, exist_ok=True)

        self.r1_path = reference.run.sample.sample.path_name_1.path
        self.r2_path = (
            reference.run.sample.sample.path_name_2.path
            if reference.run.sample.sample.exist_file_2()
            else ""
        )
        self.contigs_path = self.find_run_contigs(reference.run)

        self.taxid = reference.taxid
        self.accid = reference.accid

        if self.technology == "ONT":
            self.params = Params_Nanopore
        else:
            self.params = Params_Illumina

    def find_run_contigs(self, run_main: RunMain) -> str:
        if not run_main:
            return ""

        try:
            run_assembly = RunAssembly.objects.get(run=run_main)
            recover_assembly_contigs(run_main, run_assembly)
            assembly_contigs = run_assembly.assembly_contigs
        except RunAssembly.DoesNotExist:
            run_assembly = None
            assembly_contigs = ""

        return assembly_contigs

    def input_read_project_path(self, filepath) -> str:
        if not os.path.isfile(filepath):
            return ""
        rname = os.path.basename(filepath)

        new_rpath = os.path.join(self.dir, "reads") + "/" + rname
        shutil.copy(filepath, new_rpath)
        return new_rpath

    def generate_method_args(self):
        parameter_set = self.reference.run.parameter_set

        pipeline_tree = self.utils.parameter_util.convert_softwaretree_to_pipeline_tree(
            parameter_set.leaf.software_tree
        )
        ps_leaves = self.utils.get_parameterset_leaves(parameter_set, pipeline_tree)
        parameter_leaf_index = ps_leaves[0]
        parameter_leaf = SoftwareTreeNode.objects.get(
            index=parameter_leaf_index, software_tree=parameter_set.leaf.software_tree
        )

        run_df = self.utils.get_leaf_parameters(parameter_leaf)

        self.method_args = run_df[run_df.module == CS.PIPELINE_NAME_remapping]

        if self.method_args.empty:
            raise ValueError(
                f"no remapping parameters found for {self.reference.accid} in leaf {parameter_leaf}"
            )

    def generate_config(self):
        self.config = {
            "sample_name": simplify_name_lower(
                os.path.basename(self.r1_path).replace(".fastq.gz", "")
            ),
            "source": self.install_registry.SOURCE,
            "technology": self.technology,
            "deployment_root_dir": self.deployment_root_dir,
            "sub_directory": self.dir_branch,
            "directories": {},
            "threads": self.threads,
            "prefix": self.prefix,
            "project_name": self.project,
            "metadata": {
                x: os.path.join(self.install_registry.METADATA["ROOT"], g)
                for x, g in self.install_registry.METADATA.items()
            },
            "bin": self.install_registry.BINARIES,
            "taxid": self.taxid,
            "accid": self.accid,
            "clean": self.clean,
        }

        for dr, g in ConstantsSettings.DIRS.items():
            self.config["directories"][dr] = os.path.join(self.dir, g)
            os.makedirs(self.config["directories"][dr], exist_ok=True)

        self.config["r1"] = self.input_read_project_path(self.r1_path)
        self.config["r2"] = self.input_read_project_path(self.r2_path)
        self.config["contig_file"] = self.contigs_path
        self.config["type"] = [
            ConstantsSettings.SINGLE_END,
            ConstantsSettings.PAIR_END,
        ][int(os.path.isfile(self.config["r2"]))]

        self.config.update(self.params.CONSTANTS)

    def update_raw_reference_status_mapped(self):
        self.reference.update_raw_reference_status_fail()

    def update_raw_reference_status_fail(self):
        self.reference.update_raw_reference_status_fail()

    def engine_report_modify_mapping_success(self, run_class: RunEngine):
        def render_classification_source(record: RawReference):
            return record.classification_source_str

        run_class.report["mapping_success"] = render_classification_source(
            self.reference
        )

    def update_final_report(self, run_class: RunEngine):
        run = self.reference.run
        sample = run.sample

        self.engine_report_modify_mapping_success(run_class)

        Update_FinalReport(run_class, run, sample)

        for ref_map in run_class.remap_manager.mapped_instances:
            Update_ReferenceMap_Update(ref_map, run, sample)

    def run_reference_overlap_analysis(self):
        run = self.reference.run
        sample = run.sample
        final_report = FinalReport.objects.filter(sample=sample, run=run).order_by(
            "-coverage"
        )
        #
        report_layout_params = TelevirParameters.get_report_layout_params(run_pk=run.pk)
        report_sorter = ReportSorter(final_report, report_layout_params)
        report_sorter.sort_reports_save()


class Command(BaseCommand):
    help = "deploy run"

    def add_arguments(self, parser):
        parser.add_argument(
            "--ref_id",
            type=int,
            help="user deploying the run (pk)",
        )

        parser.add_argument(
            "--project_id",
            type=int,
            help="project (pk)",
        )

        parser.add_argument(
            "-o",
            "--outdir",
            type=str,
            help="output directory",
        )

    def handle(self, *args, **options):
        ###
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        raw_reference_id = int(options["ref_id"])
        project_pk = int(options["project_id"])

        reference = RawReference.objects.get(pk=raw_reference_id)
        project = Projects.objects.get(pk=project_pk)
        user = reference.run.project.owner
        project_name = project.name

        ######## register map
        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_map(reference.pk),
            ProcessControler.FLAG_RUNNING,
        )

        ########

        input_generator = Input_Generator(
            reference, options["outdir"], threads=ConstantsSettings.MAPPING_THREADS
        )

        try:
            input_generator.generate_method_args()
            input_generator.generate_config()

            run_engine = RunEngine(
                input_generator.config,
                input_generator.method_args,
                project_name,
                user.username,
            )
            run_engine.generate_targets()
            run_engine.run()
            run_engine.export_sequences()
            input_generator.update_raw_reference_status_mapped()
            input_generator.update_final_report(run_engine)
            input_generator.run_reference_overlap_analysis()

            ######## register map sucess
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_map(reference.pk),
                ProcessControler.FLAG_FINISHED,
            )

            ########

        except Exception as e:
            print(e)
            input_generator.update_raw_reference_status_fail()

            ######## register map sucess
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_map(reference.pk),
                ProcessControler.FLAG_ERROR,
            )

            ########
