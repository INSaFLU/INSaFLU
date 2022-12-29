import logging
import os
import shutil
import time
from random import randint
from typing import Type

import numpy as np
import pandas as pd
from pathogen_identification.constants_settings import ConstantsSettings
from pathogen_identification.modules.assembly_class import Assembly_class
from pathogen_identification.modules.classification_class import Classifier
from pathogen_identification.modules.metadata_handler import Metadata_handler
from pathogen_identification.modules.object_classes import (
    Assembly_results,
    Contig_classification_results,
    Read_class,
    Read_classification_results,
    Remap_main,
    Run_detail_report,
    RunCMD,
    Sample_runClass,
    Software_detail,
)
from pathogen_identification.modules.preprocess_class import Preprocess
from pathogen_identification.modules.remap_class import Mapping_Manager
from settings.constants_settings import ConstantsSettings as CS


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


class RunDetail_main:
    suprun: str = None

    threads: int
    config: dict
    prefix: str
    project_name: str
    username: str
    ## input
    sample_name: str
    type: str
    r1_suffix: str
    r2_suffix: str

    r1: Read_class
    r2: Read_class

    sample: Sample_runClass
    ##  metadata
    metadata_tool: Metadata_handler
    sift_query: str
    max_remap: int
    taxid_limit: int

    ## actions
    quality_control: bool
    sift: bool
    depletion: bool
    enrichment: bool
    assembly: bool
    classification: bool
    remapping: bool
    house_cleaning: bool

    ## methods
    preprocess_method: Software_detail
    depletion_method: Software_detail
    enrichment_method: Software_detail
    assembly_method: Software_detail
    assembly_classification_method: Software_detail
    read_classification_method: Software_detail
    remapping_method: Software_detail
    remap_manager = Mapping_Manager

    ## directories.
    root: str

    input_reads_dir: str
    filtered_reads_dir: str
    depleted_reads_dir: str

    log_dir: str

    dir_classification: str = f"classification_reports"
    dir_plots: str = f"plots"
    igv_dir: str = f"igv"

    ## output content
    report: pd.DataFrame

    def __init__(self, config: dict, method_args: pd.DataFrame, username: str):

        self.project_name = config["project_name"]
        self.username = username
        self.prefix = config["prefix"]
        self.suprun = self.prefix

        self.method_args = method_args
        self.config = config
        self.cmd = RunCMD(
            get_bindir_from_binaries(
                config["bin"], CS.PIPELINE_NAME_read_quality_analysis
            )
        )
        self.threads = config["threads"]

        self.logger_level_main = logging.INFO
        self.logger_level_detail = logging.INFO
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
        self.runtime = 0
        self.start_time = time.perf_counter()
        self.exec_time = 0

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

        ######### INPUT
        self.sample_name = config["sample_name"]
        self.type = config["type"]
        self.run_detail_report = pd.DataFrame()
        self.aclass_summary = pd.DataFrame([[0]], columns=["counts"])
        self.rclass_summary = pd.DataFrame([[0]], columns=["counts"])
        self.merged_targets = pd.DataFrame(columns=["taxid"])
        self.raw_targets = pd.DataFrame(columns=["taxid"])
        self.remap_plan = pd.DataFrame()
        self.report = pd.DataFrame()

        self.r1 = Read_class(
            config["r1"],
            config["directories"][CS.PIPELINE_NAME_read_quality_analysis],
            config["directories"]["reads_enriched_dir"],
            config["directories"]["reads_depleted_dir"],
            bin=get_bindir_from_binaries(
                config["bin"], CS.PIPELINE_NAME_read_quality_analysis
            ),
        )

        self.r1.cmd = RunCMD(
            logdir=self.log_dir, bin=self.r1.cmd.bin, prefix="r1", task="housekeeping"
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

        self.r2.cmd = RunCMD(
            logdir=self.log_dir, bin=self.r2.cmd.bin, prefix="r2", task="housekeeping"
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

        ### mapping parameters
        self.min_scaffold_length = config["assembly_contig_min_length"]
        self.minimum_coverage = int(config["minimum_coverage_threshold"])
        self.maximum_coverage = 1000000000

        ### metadata
        self.metadata_tool = Metadata_handler(
            self.config, sift_query=config["sift_query"], prefix=self.prefix
        )

        self.max_remap = config["max_output_number"]
        self.taxid_limit = config["taxid_limit"]

        ### methods
        self.preprocess_method = Software_detail(
            CS.PIPELINE_NAME_read_quality_analysis,
            self.method_args,
            config,
            self.prefix,
        )
        self.assembly_method = Software_detail(
            CS.PIPELINE_NAME_assembly,
            method_args,
            config,
            self.prefix,
        )
        self.depletion_method = Software_detail(
            CS.PIPELINE_NAME_host_depletion,
            method_args,
            config,
            self.prefix,
        )

        self.enrichment_method = Software_detail(
            CS.PIPELINE_NAME_viral_enrichment,
            method_args,
            config,
            self.prefix,
        )

        self.contig_classification_method = Software_detail(
            CS.PIPELINE_NAME_contig_classification,
            method_args,
            config,
            self.prefix,
        )
        self.read_classification_method = Software_detail(
            CS.PIPELINE_NAME_read_classification,
            method_args,
            config,
            self.prefix,
        )

        self.remapping_method = Software_detail(
            CS.PIPELINE_NAME_remapping,
            method_args,
            config,
            self.prefix,
        )

        ### actions
        self.subsample = False
        self.quality_control = config["actions"]["QCONTROL"]
        self.sift = config["actions"]["SIFT"]
        self.depletion = bool(self.depletion_method.name != "None")
        self.depletion = bool(self.depletion_method.name != "None")
        self.enrichment = bool(self.enrichment_method.name != "None")
        self.assembly = bool(self.assembly_method.name != "None")
        self.classification = config["actions"]["CLASSIFY"]
        self.remapping = config["actions"]["REMAP"]
        self.house_cleaning = config["actions"]["CLEAN"]

        ### drones
        self.depletion_drone = Classifier(
            Software_detail("NONE", method_args, config, self.prefix),
            logging_level=self.logger_level_detail,
        )
        self.enrichment_drone = Classifier(
            Software_detail("NONE", method_args, config, self.prefix),
            logging_level=self.logger_level_detail,
        )

        ### output files
        self.params_file_path = os.path.join(
            self.media_dir_classification,
            f"{self.prefix}_params.csv",
        )
        self.remap_plan_path = os.path.join(
            self.media_dir_classification,
            f"{self.prefix}_remap_plan.csv",
        )
        self.full_report = os.path.join(
            self.media_dir_classification,
            f"{self.prefix}_full_report.tsv",
        )
        self.assembly_classification_summary = os.path.join(
            self.media_dir_classification,
            f"{self.prefix}_aclass_summary.tsv",
        )
        self.read_classification_summary = os.path.join(
            self.media_dir_classification,
            f"{self.prefix}_rclass_summary.tsv",
        )
        self.merged_classification_summary = os.path.join(
            self.media_dir_classification,
            f"{self.prefix}_mclass_summary.tsv",
        )

    def Update(self, config: dict, method_args: pd.DataFrame):

        self.method_args = pd.concat((self.method_args, method_args))
        # with open(config_json) as json_file:
        #    config = json.load(json_file)

        self.config = config
        self.prefix = config["prefix"]
        self.type = config["type"]
        self.logger = logging.getLogger("{}".format(self.prefix))
        self.logger.setLevel(self.logger_level_main)
        logFormatter = logging.Formatter(fmt="{} :%(message)s".format(self.prefix))
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        self.logger.addHandler(consoleHandler)
        self.logger.info(f"prefix: {self.prefix}")
        self.logger.info(f"type: {self.type}")
        self.start_time = time.perf_counter()

        # directories
        self.filtered_reads_dir = config["directories"][
            CS.PIPELINE_NAME_read_quality_analysis
        ]
        self.log_dir = config["directories"]["log_dir"]

        self.sample.r1.update(
            self.sample.r1,
            clean_dir=config["directories"][CS.PIPELINE_NAME_read_quality_analysis],
            enriched_dir=config["directories"]["reads_enriched_dir"],
            depleted_dir=config["directories"]["reads_depleted_dir"],
        )

        self.sample.r2.update(
            self.sample.r2,
            clean_dir=config["directories"][CS.PIPELINE_NAME_read_quality_analysis],
            enriched_dir=config["directories"]["reads_enriched_dir"],
            depleted_dir=config["directories"]["reads_depleted_dir"],
        )

        ### actions
        self.subsample = False
        self.quality_control = config["actions"]["QCONTROL"]
        self.sift = config["actions"]["SIFT"]
        self.depletion = config["actions"]["DEPLETE"]
        self.enrichment = config["actions"]["ENRICH"]
        self.assembly = config["actions"]["ASSEMBLY"]
        self.classification = config["actions"]["CLASSIFY"]
        self.remapping = config["actions"]["REMAPPING"]
        self.house_cleaning = config["actions"]["CLEANING"]

        ### methods

        self.preprocess_method = Software_detail(
            CS.PIPELINE_NAME_read_quality_analysis,
            self.method_args,
            config,
            self.prefix,
        )
        self.assembly_method = Software_detail(
            CS.PIPELINE_NAME_assembly,
            method_args,
            config,
            self.prefix,
        )

        self.depletion_method = Software_detail(
            CS.PIPELINE_NAME_host_depletion,
            method_args,
            config,
            self.prefix,
        )

        self.enrichment_method = Software_detail(
            CS.PIPELINE_NAME_viral_enrichment,
            method_args,
            config,
            self.prefix,
        )

        self.contig_classification_method = Software_detail(
            CS.PIPELINE_NAME_contig_classification,
            method_args,
            config,
            self.prefix,
        )
        self.read_classification_method = Software_detail(
            CS.PIPELINE_NAME_read_classification,
            method_args,
            config,
            self.prefix,
        )

        self.remapping_method = Software_detail(
            CS.PIPELINE_NAME_remapping,
            method_args,
            config,
            self.prefix,
        )

        ### output files
        self.params_file_path = os.path.join(
            self.media_dir_classification,
            f"{self.prefix}_params.csv",
        )
        self.remap_plan_path = os.path.join(
            self.media_dir_classification,
            f"{self.prefix}_remap_plan.csv",
        )
        self.full_report = os.path.join(
            self.media_dir_classification,
            f"{self.prefix}_full_report.tsv",
        )
        self.assembly_classification_summary = os.path.join(
            self.media_dir_classification,
            f"{self.prefix}_aclass_summary.tsv",
        )
        self.read_classification_summary = os.path.join(
            self.media_dir_classification,
            f"{self.prefix}_rclass_summary.tsv",
        )
        self.merged_classification_summary = os.path.join(
            self.media_dir_classification,
            f"{self.prefix}_mclass_summary.tsv",
        )

    def Update_exec_time(self):
        """
        Update the execution time of the pipeline.
        """
        self.exec_time = self.exec_time + time.perf_counter() - self.start_time


class Run_Deployment_Methods(RunDetail_main):
    def __init__(
        self, config_json: os.PathLike, method_args: pd.DataFrame, username: str
    ):
        super().__init__(config_json, method_args, username)
        self.mapped_instances = []

    def Prep_deploy(self, fake_run: bool = False):

        self.preprocess_drone = Preprocess(
            self.sample.r1.current,
            self.sample.r2.current,
            self.filtered_reads_dir,
            self.type,
            self.preprocess_method,
            self.sample.r1.clean,
            self.sample.r2.clean,
            self.threads,
            self.subsample,
            logging_level=self.logger_level_detail,
            log_dir=self.log_dir,
        )

        self.depletion_drone = Classifier(
            self.depletion_method,
            self.sample.r1.current,
            type=self.type,
            r2=self.sample.r2.current,
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(
                self.config["bin"], CS.PIPELINE_NAME_remapping
            ),
            logging_level=self.logger_level_detail,
            log_dir=self.log_dir,
        )

        self.enrichment_drone = Classifier(
            self.enrichment_method,
            self.sample.r1.current,
            type=self.type,
            r2=self.sample.r2.current,
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(
                self.config["bin"], CS.PIPELINE_NAME_remapping
            ),
            logging_level=self.logger_level_detail,
            log_dir=self.log_dir,
        )

        self.assembly_drone = Assembly_class(
            self.sample.r1.current,
            self.assembly_method,
            self.type,
            min_scaffold_length=self.min_scaffold_length,
            r2=self.sample.r2.current,
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(
                self.config["bin"], CS.PIPELINE_NAME_remapping
            ),
            logging_level=self.logger_level_detail,
            log_dir=self.log_dir,
        )

        self.contig_classification_drone = Classifier(
            self.contig_classification_method,
            self.assembly_drone.assembly_file_fasta_gz,
            r2="",
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(
                self.config["bin"], CS.PIPELINE_NAME_remapping
            ),
            logging_level=self.logger_level_detail,
            log_dir=self.log_dir,
        )

        self.read_classification_drone = Classifier(
            self.read_classification_method,
            self.sample.r1.current,
            type=self.type,
            r2=self.sample.r2.current,
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(
                self.config["bin"], CS.PIPELINE_NAME_remapping
            ),
            logging_level=self.logger_level_detail,  #
            log_dir=self.log_dir,
        )

        self.remap_manager = Mapping_Manager(
            [],
            self.sample.r1,
            self.sample.r2,
            self.remapping_method,
            self.assembly_drone.assembly_file_fasta_gz,
            self.type,
            self.prefix,
            self.threads,
            self.minimum_coverage,
            get_bindir_from_binaries(self.config["bin"], CS.PIPELINE_NAME_remapping),
            self.logger_level_detail,
            True,
            logdir=self.config["directories"]["log_dir"],
        )

    def deploy_QC(self, fake_run: bool = False):

        self.logger.info(f"r1 reads: {self.sample.r1.get_current_fastq_read_number()}")
        self.logger.info(f"r2 reads: {self.sample.r2.get_current_fastq_read_number()}")

        if fake_run:
            self.preprocess_drone.fake_run()
        else:
            self.preprocess_drone.run()

    def deploy_HD(self):
        self.depletion_drone = Classifier(
            self.depletion_method,
            self.sample.r1.current,
            type=self.type,
            r2=self.sample.r2.current,
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(
                self.config["bin"], CS.PIPELINE_NAME_remapping
            ),
            logging_level=self.logger_level_detail,
            log_dir=self.log_dir,
        )

        self.depletion_drone.run()

    def deploy_EN(self):
        self.enrichment_drone = Classifier(
            self.enrichment_method,
            self.sample.r1.current,
            type=self.type,
            r2=self.sample.r2.current,
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(
                self.config["bin"], CS.PIPELINE_NAME_remapping
            ),
            logging_level=self.logger_level_detail,
            log_dir=self.log_dir,
        )
        self.enrichment_drone.run()

    def deploy_ASSEMBLY(self, fake_run: bool = False):
        self.assembly_drone = Assembly_class(
            self.sample.r1.current,
            self.assembly_method,
            self.type,
            min_scaffold_length=self.min_scaffold_length,
            r2=self.sample.r2.current,
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(
                self.config["bin"], CS.PIPELINE_NAME_remapping
            ),
            logging_level=self.logger_level_detail,
            log_dir=self.log_dir,
        )
        if fake_run:
            self.assembly_drone.fake_run()
        else:
            self.assembly_drone.run()

    def deploy_CONTIG_CLASSIFICATION(self):
        self.contig_classification_drone = Classifier(
            self.contig_classification_method,
            self.assembly_drone.assembly_file_fasta_gz,
            r2="",
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(
                self.config["bin"], CS.PIPELINE_NAME_remapping
            ),
            logging_level=self.logger_level_detail,
            log_dir=self.log_dir,
        )

        self.contig_classification_drone.run()

    def deploy_READ_CLASSIFICATION(self):
        self.read_classification_drone = Classifier(
            self.read_classification_method,
            self.sample.r1.current,
            type=self.type,
            r2=self.sample.r2.current,
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(
                self.config["bin"], CS.PIPELINE_NAME_remapping
            ),
            logging_level=self.logger_level_detail,  #
            log_dir=self.log_dir,
        )
        self.read_classification_drone.run()

    def deploy_REMAPPING(self):

        self.remap_manager = Mapping_Manager(
            self.metadata_tool.remap_targets,
            self.sample.r1,
            self.sample.r2,
            self.remapping_method,
            self.assembly_drone.assembly_file_fasta_gz,
            self.type,
            self.prefix,
            self.threads,
            self.minimum_coverage,
            get_bindir_from_binaries(self.config["bin"], CS.PIPELINE_NAME_remapping),
            self.logger_level_detail,
            True,
            logdir=self.config["directories"]["log_dir"],
        )

        self.logger.info(
            f"{self.prefix} remapping # targets: {len(self.metadata_tool.remap_targets)}"
        )

        self.remap_manager.run_mappings_move_clean(
            self.static_dir_plots, self.media_dir_igv
        )

        self.remap_manager.merge_mapping_reports()
        self.remap_manager.collect_final_report_summary_statistics()


class RunMain_class(Run_Deployment_Methods):
    def __init__(
        self, config_json: os.PathLike, method_args: pd.DataFrame, username: str
    ):
        super().__init__(config_json, method_args, username)

    def Run_Full_Pipeline(self):

        self.Prep_deploy()
        self.Run_QC()
        self.Run_PreProcess()
        self.Sanitize_reads()
        self.Run_Assembly()
        self.Run_Classification()
        self.Run_Remapping()

    def Run_QC(self):

        self.logger.info("Starting Pipeline")

        self.logger.info(f"quality control: {self.quality_control}")
        self.logger.info(f"enrichment: {self.enrichment}")
        self.logger.info(f"depletion: {self.depletion}")
        self.logger.info(f"assembly: {self.assembly}")
        self.logger.info(f"classification: {self.classification}")
        self.logger.info(f"sift: {self.sift}")
        self.logger.info(f"remapping: {self.remapping}")

        if self.quality_control:
            self.deploy_QC()

            self.sample.r1.is_clean()
            self.sample.r2.is_clean()

            self.sample.qc_soft = self.preprocess_drone.preprocess_method.name
            self.sample.input_fastqc_report = self.preprocess_drone.input_qc_report
            self.sample.processed_fastqc_report = (
                self.preprocess_drone.processed_qc_report
            )

            self.sample.reads_after_processing = self.sample.current_total_read_number()
            self.sample.get_qc_data()
            self.sample.r1.clean_read_names()
            self.sample.r2.clean_read_names()

        else:
            self.deploy_QC(fake_run=True)

            shutil.copy(self.sample.r1.current, self.sample.r1.clean)

            if self.sample.r2.exists:
                shutil.copy(self.sample.r2.current, self.sample.r2.clean)

            self.sample.qc_soft = "none"
            self.sample.input_fastqc_report = self.preprocess_drone.input_qc_report
            self.sample.processed_fastqc_report = (
                self.preprocess_drone.processed_qc_report
            )

            self.sample.r1.is_clean()
            self.sample.r2.is_clean()
            self.sample.reads_after_processing = self.sample.current_total_read_number()
            self.sample.get_fake_qc_data()

            self.sample.r1.clean_read_names()
            self.sample.r2.clean_read_names()

        self.Update_exec_time()

    def Run_PreProcess(self):

        if self.enrichment:
            self.deploy_EN()

            self.logger.info(
                "r1 current reads: "
                + str(self.sample.r1.get_current_fastq_read_number())
            )

            self.sample.r1.enrich(self.enrichment_drone.classified_reads_list)
            self.sample.r2.enrich(self.enrichment_drone.classified_reads_list)

            self.logger.info(
                "r1 current reads after enrichment: "
                + str(self.sample.r1.get_current_fastq_read_number())
            )

        if self.depletion:
            self.deploy_HD()

            self.sample.r1.deplete(self.depletion_drone.classified_reads_list)
            self.sample.r2.deplete(self.depletion_drone.classified_reads_list)

        self.Update_exec_time()
        self.generate_output_data_classes()

    def Sanitize_reads(self):
        if self.enrichment or self.depletion or self.assembly:
            self.logger.info(
                "r1 current before trim: "
                + str(self.sample.r1.get_current_fastq_read_number())
            )
            self.sample.trimmomatic_sort()
            self.sample.remove_duplicates()
            self.logger.info(
                "r1 current after trim: "
                + str(self.sample.r1.get_current_fastq_read_number())
            )

            self.generate_output_data_classes()

    def Run_Assembly(self):

        if self.assembly:
            self.deploy_ASSEMBLY()
        else:
            self.deploy_ASSEMBLY(fake_run=True)

        self.Update_exec_time()
        self.generate_output_data_classes()

    def Run_Classification(self):

        if self.classification:
            self.deploy_READ_CLASSIFICATION()
            self.deploy_CONTIG_CLASSIFICATION()

            self.metadata_tool.match_and_select_targets(
                self.read_classification_drone.classification_report,
                self.contig_classification_drone.classification_report,
                self.max_remap,
                self.taxid_limit,
            )
            self.aclass_summary = self.metadata_tool.aclass
            self.rclass_summary = self.metadata_tool.rclass
            self.merged_targets = self.metadata_tool.merged_targets
            self.raw_targets = self.metadata_tool.raw_targets
            self.remap_plan = self.metadata_tool.remap_plan

            self.export_intermediate_reports()

        self.Update_exec_time()
        self.generate_output_data_classes()

    def Run_Remapping(self):
        if self.remapping:
            self.deploy_REMAPPING()
            self.report = self.remap_manager.report
            self.export_final_reports()

        self.Update_exec_time()

    #### SUMMARY FUNCTIONS ####

    def export_logdir(self):

        if os.path.exists(self.media_dir_logdir):
            shutil.rmtree(self.media_dir_logdir)

        shutil.copytree(
            self.log_dir,
            self.media_dir_logdir,
        )

    def export_final_reports(self):

        ### main report
        self.report.to_csv(
            self.full_report,
            index=False,
            sep="\t",
            header=True,
        )

    def save_df_check_exists(self, df: pd.DataFrame, path: str):
        if not os.path.exists(path):
            df.to_csv(path, index=False, sep="\t", header=True)

    def export_intermediate_reports(self):
        export_dict = {
            self.params_file_path: self.method_args,
            self.remap_plan_path: self.remap_plan,
            self.assembly_classification_summary: self.aclass_summary,
            self.read_classification_summary: self.rclass_summary,
            self.merged_classification_summary: self.merged_targets,
        }
        for output_df_path, df in export_dict.items():

            self.save_df_check_exists(df, output_df_path)

    def export_sequences(self):
        self.sample.export_reads(self.media_dir)
        self.assembly_drone.export_assembly(self.media_dir)

    def Summarize(self):

        self.logger.info(f"prefix: {self.prefix}")
        with open(os.path.join(self.log_dir, self.prefix + "_latest.fofn"), "w") as f:
            f.write(self.sample.r1.current + "\n")
            if self.type == "PE":
                f.write(self.sample.r2.current + "\n")

        with open(os.path.join(self.log_dir, "reads_latest.stats"), "w") as f:
            f.write(f"CLEAN\t{self.sample.r1.read_number_clean}\n")
            f.write(f"ENRICHED\t{self.sample.r1.read_number_enriched}\n")

    def generate_output_data_classes(self):
        ### transfer to sample class
        processed_reads = self.sample.reads_after_processing

        filtered_reads = (
            self.sample.r1.read_number_filtered + self.sample.r2.read_number_filtered
        )

        final_processing_reads = (
            self.sample.r1.current_fastq_read_number()
            + self.sample.r2.current_fastq_read_number()
        )

        filtered_reads_perc = (int(filtered_reads) / processed_reads) * 100
        final_processing_percent = (final_processing_reads / processed_reads) * 100

        ### transfer to assembly class / drone.
        print(self.aclass_summary)
        minhit_assembly = self.aclass_summary["counts"].min()
        if not minhit_assembly or not self.aclass_summary.shape[0]:
            minhit_assembly = 0

        minhit_reads = self.rclass_summary["counts"].min()
        if np.isnan(minhit_reads):
            minhit_reads = 0

        files = list(
            set(
                [
                    os.path.basename(t.reference.target.file)
                    for t in self.remap_manager.mapped_instances
                ]
            )
        )

        enriched_reads = len(self.enrichment_drone.classified_reads_list)
        depleted_reads = len(self.depletion_drone.classified_reads_list)

        if self.type == "PE":
            enriched_reads = enriched_reads * 2
            depleted_reads = depleted_reads * 2

        self.run_detail_report = Run_detail_report(
            self.remap_manager.max_depth,
            self.remap_manager.max_depthR,
            self.remap_manager.max_gaps,
            self.remap_manager.max_prop,
            self.remap_manager.max_mapped,
            f"{processed_reads:,}",
            enriched_reads,
            enriched_reads / processed_reads,
            depleted_reads,
            depleted_reads / processed_reads,
            f"{filtered_reads:,}",
            f"{filtered_reads_perc:.2f}",
            False,
            self.sift,
            f"{self.metadata_tool.sift_report.loc[0]['removed']:,}",
            f"{final_processing_reads:,}",
            round(final_processing_percent, 3),
            self.remapping,
            self.merged_targets.taxid.nunique(),
            ", ".join(files),
        )

        self.contig_classification_results = Contig_classification_results(
            True
            if self.contig_classification_drone.classifier_method.name != "None"
            else False,
            self.contig_classification_drone.classifier_method.name,
            self.contig_classification_drone.classifier_method.args,
            self.contig_classification_drone.classifier_method.db_name,
            self.aclass_summary.shape[0],
            minhit_assembly,
            self.aclass_summary.shape[0] > 0,
        )

        self.read_classification_results = Read_classification_results(
            True
            if self.read_classification_drone.classifier_method.name != "None"
            else False,
            self.read_classification_drone.classifier_method.name,
            self.read_classification_drone.classifier_method.args,
            self.read_classification_drone.classifier_method.db_name,
            self.rclass_summary.shape[0],
            minhit_reads,
            self.rclass_summary.shape[0] > 0,
        )

        self.assembly_report = Assembly_results(
            self.assembly,
            self.assembly_drone.assembly_method.name,
            self.assembly_drone.assembly_method.args,
            self.assembly_drone.assembly_number,
            f"{self.assembly_drone.assembly_min:,}",
            f"{int(self.assembly_drone.assembly_mean):,}",
            f"{self.assembly_drone.assembly_max:,}",
            f"{int(self.min_scaffold_length):,}",
        )

        self.remap_main = Remap_main(
            True,
            self.report.shape[0],
            self.remapping_method.name,
            len(self.mapped_instances),
            self.minimum_coverage,
            self.maximum_coverage,
        )
