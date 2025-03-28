import logging
import os
import shutil
import time
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from pathogen_identification.constants_settings import ConstantsSettings
from pathogen_identification.models import PIProject_Sample, Projects
from pathogen_identification.modules.assembly_class import Assembly_class
from pathogen_identification.modules.classification_class import Classifier
from pathogen_identification.modules.metadata_handler import RunMetadataHandler
from pathogen_identification.modules.object_classes import (
    Assembly_results,
    Contig_classification_results,
    Read_class,
    Read_classification_results,
    Remap_main,
    Remap_Target,
    Run_detail_report,
    RunCMD,
    RunQC_report,
    Sample_runClass,
    SoftwareDetail,
    SoftwareDetailCompound,
    SoftwareRemap,
    SoftwareUnit,
)
from pathogen_identification.modules.preprocess_class import Preprocess
from pathogen_identification.modules.remap_class import (
    Mapping_Instance,
    Mapping_Manager,
)
from pathogen_identification.utilities.televir_parameters import (
    RemapParams,
    TelevirParameters,
)
from pathogen_identification.utilities.utilities_views import RawReferenceUtils
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


@dataclass
class RunMainConfig:
    project_name: str
    sample_name: str
    r1: str
    r2: str
    type: str
    prefix: str
    source: str
    deployment_root_dir: str
    sub_directory: str
    directories: dict
    threads: int
    metadata: dict
    technology: str
    bin: dict
    actions: dict


class RunDetail_main:
    suprun: str = None

    RUN_TYPE_PIPELINE = 1
    RUN_TYPE_MAPPING_REQUEST = 2
    RUN_TYPE_SCREENING = 3
    RUN_TYPE_COMBINED_MAPPING = 4
    RUN_TYPE_PANEL_MAPPING = 5

    threads: int
    config: dict
    prefix: str
    project_name: str
    username: str
    run_type: int
    run_pk: Optional[int]
    ## input
    sample_name: str
    type: str
    r1_suffix: str
    r2_suffix: str

    r1: Read_class
    r2: Read_class

    sample: Sample_runClass
    sample_registered: PIProject_Sample
    ##  metadata
    metadata_tool: RunMetadataHandler
    sift_query: str
    max_remap: int
    taxid_limit: int

    ## actions
    modules: List[str]

    quality_control: bool
    sift: bool
    depletion: bool
    enrichment: bool
    assembly: bool
    classification: bool
    contig_classification: bool
    read_classification: bool
    metagenomics_classification: bool
    remapping: bool
    house_cleaning: bool
    depletion_report = None

    ## methods
    preprocess_method: SoftwareUnit
    depletion_method: SoftwareDetail
    enrichment_method: SoftwareDetail
    assembly_method: SoftwareDetail
    assembly_classification_method: SoftwareDetail
    read_classification_method: SoftwareDetail
    metagenomics_classification_method: SoftwareDetail
    remapping_method: SoftwareDetail
    remap_manager: Mapping_Manager
    remap_params: RemapParams
    software_remap: SoftwareRemap
    active_methods: Dict[str, SoftwareUnit]
    ## directories.
    root: str

    input_reads_dir: str
    filtered_reads_dir: str
    depleted_reads_dir: str

    log_dir: str

    dir_classification: str = "classification_reports"
    dir_plots: str = "plots"
    igv_dir: str = "igv"

    ## output content
    report: pd.DataFrame

    def set_logger(self):
        self.logger_level_main = logging.ERROR
        self.logger_level_detail = logging.ERROR

        self.logger = logging.getLogger("main {}".format(self.prefix))
        self.logger.setLevel(self.logger_level_main)

        logFormatter = logging.Formatter(
            fmt="{} %(levelname)s :%(message)s".format(self.prefix)
        )

        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        self.logger.addHandler(consoleHandler)
        self.logger.propagate = False

    def delete_logger(self):
        self.logger = None

    def set_preprocess_check(self, config: dict, method_args: pd.DataFrame):
        self.preprocess_method = SoftwareDetail(
            CS.PIPELINE_NAME_extra_qc,
            method_args,
            config,
            self.prefix,
        )

        self.check_preprocess_exists()

    def check_preprocess_exists(self):
        self.quality_control = self.preprocess_method.check_exists()

    def set_depletion_check(self, config: dict, method_args: pd.DataFrame):
        self.depletion_method = SoftwareDetail(
            CS.PIPELINE_NAME_host_depletion,
            method_args,
            config,
            self.prefix,
        )

        self.check_depletion_exists()

    def check_depletion_exists(self):
        self.depletion = self.depletion_method.check_exists()

    def set_enrichment_check(self, config: dict, method_args: pd.DataFrame):
        self.enrichment_method = SoftwareDetail(
            CS.PIPELINE_NAME_viral_enrichment,
            method_args,
            config,
            self.prefix,
        )

        self.enrichment = self.enrichment_method.check_exists()

    def check_enrichment_exists(self):
        self.enrichment = self.enrichment_method.check_exists()

    def set_assembly_check(self, config: dict, method_args: pd.DataFrame):
        self.assembly_method = SoftwareDetail(
            CS.PIPELINE_NAME_assembly,
            method_args,
            config,
            self.prefix,
        )

        self.assembly = self.assembly_method.check_exists()

    def check_assembly_exists(self):
        self.assembly = self.assembly_method.check_exists()

    def set_contig_classification_check(self, config: dict, method_args: pd.DataFrame):
        self.contig_classification_method = SoftwareDetail(
            CS.PIPELINE_NAME_contig_classification,
            method_args,
            config,
            self.prefix,
        )

        self.contig_classification = self.contig_classification_method.check_exists()

    def check_contig_classification_exists(self):
        self.contig_classification = self.contig_classification_method.check_exists()

    def set_read_classification_check(self, config: dict, method_args: pd.DataFrame):
        self.read_classification_method = SoftwareDetail(
            CS.PIPELINE_NAME_read_classification,
            method_args,
            config,
            self.prefix,
        )

        self.read_classification = self.read_classification_method.check_exists()

    def check_read_classification_exists(self):
        self.read_classification = self.read_classification_method.check_exists()

    def set_metagenomics_classification_check(
        self, config: dict, method_args: pd.DataFrame
    ):
        self.metagenomics_classification_method = SoftwareDetail(
            CS.PIPELINE_NAME_metagenomics_screening,
            method_args,
            config,
            self.prefix,
        )
        #

        self.check_metagenomics_classification_exists()

    def check_metagenomics_classification_exists(self):
        self.metagenomics_classification = (
            self.metagenomics_classification_method.check_exists()
        )

        if self.metagenomics_classification:
            self.run_type = self.RUN_TYPE_SCREENING

    def pick_remapping_method(
        self, config: dict, method_args: pd.DataFrame
    ) -> SoftwareDetail:
        if CS.PIPELINE_NAME_remapping in method_args["module"].unique():
            return SoftwareDetail(
                CS.PIPELINE_NAME_remapping,
                method_args,
                config,
                self.prefix,
            )
        elif CS.PIPELINE_NAME_request_mapping in method_args["module"].unique():
            return SoftwareDetail(
                CS.PIPELINE_NAME_request_mapping,
                method_args,
                config,
                self.prefix,
            )
        else:
            return SoftwareDetail(
                CS.PIPELINE_NAME_remapping,
                method_args,
                config,
                self.prefix,
            )

    def set_remapping_check(self, config: dict, method_args: pd.DataFrame):
        self.remapping_method = self.pick_remapping_method(config, method_args)
        self.check_remapping_exists()

    def check_remapping_exists(self):
        self.remapping = self.remapping_method.check_exists()

    def set_remapping_filtering_check(self, config: dict, method_args: pd.DataFrame):
        self.remap_filtering_method = SoftwareDetailCompound(
            [CS.PIPELINE_NAME_remap_filtering, CS.PIPELINE_NAME_map_filtering],
            method_args,
            config,
            self.prefix,
        )

        self.check_remap_filtering_exists()

    def check_remap_filtering_exists(self):
        self.remapping_filtering = self.remap_filtering_method.check_exists()

    def set_settings_dict(self):
        self.settings_dict = {
            CS.PIPELINE_NAME_extra_qc: self.set_preprocess_check,
            CS.PIPELINE_NAME_host_depletion: self.set_depletion_check,
            CS.PIPELINE_NAME_viral_enrichment: self.set_enrichment_check,
            CS.PIPELINE_NAME_assembly: self.set_assembly_check,
            CS.PIPELINE_NAME_contig_classification: self.set_contig_classification_check,
            CS.PIPELINE_NAME_read_classification: self.set_read_classification_check,
            CS.PIPELINE_NAME_metagenomics_screening: self.set_metagenomics_classification_check,
            CS.PIPELINE_NAME_remapping: self.set_remapping_check,
            CS.PIPELINE_NAME_remap_filtering: self.set_remapping_filtering_check,
        }

    def software_check_map(self):
        self.module_software_check_map = {
            CS.PIPELINE_NAME_extra_qc: self.quality_control,
            CS.PIPELINE_NAME_host_depletion: self.depletion,
            CS.PIPELINE_NAME_viral_enrichment: self.enrichment,
            CS.PIPELINE_NAME_assembly: self.assembly,
            CS.PIPELINE_NAME_contig_classification: self.contig_classification,
            CS.PIPELINE_NAME_read_classification: self.read_classification,
            CS.PIPELINE_NAME_metagenomics_screening: self.metagenomics_classification,
            CS.PIPELINE_NAME_remapping: self.remapping,
            CS.PIPELINE_NAME_remap_filtering: self.remapping_filtering,
        }

    def check_software_print(self):
        self.software_check_map()

        for module, software in self.module_software_check_map.items():
            print(f"{module} : {software}")

    def set_methods(self, config: dict, method_args: pd.DataFrame):
        self.set_settings_dict()

        for _, software_get_function in self.settings_dict.items():
            software_get_function(config, method_args)

        self.software_remap = SoftwareRemap(
            self.remapping_method,
            self.remap_filtering_method,
        )

    def update_methods(self, config: dict, method_args: pd.DataFrame):
        self.set_settings_dict()

        for _, software_get_function in self.settings_dict.items():
            software_get_function(config, method_args)

        ###
        self.software_remap.remap_software = self.remapping_method
        self.software_remap.remap_filters.args_df = self.remap_filtering_method.args_df
        self.software_remap.remap_filters.config = self.remap_filtering_method.config
        self.software_remap.remap_filters.prefix = self.remap_filtering_method.prefix
        self.software_remap.remap_filters.software_list.extend(
            self.remap_filtering_method.software_list
        )

    def __init__(self, config: dict, method_args: pd.DataFrame, project_pk: int):
        project = Projects.objects.get(pk=project_pk)
        self.project_name = project.name
        self.username = project.owner.username

        self.prefix = config["prefix"]
        self.suprun = self.prefix
        self.run_type = self.RUN_TYPE_PIPELINE
        self.project_pk = project_pk
        self.run_pk = None
        self.ps_pk = None

        self.method_args = method_args
        self.modules = list(self.method_args["module"].unique())
        self.config = config
        self.log_dir = config["directories"]["log_dir"]
        print("logdir", self.log_dir, "run type", self.run_type)

        self.cmd = RunCMD(
            get_bindir_from_binaries(
                config["bin"],
                CS.PIPELINE_NAME_read_quality_analysis,
            ),
            logdir=self.log_dir,
            prefix=self.prefix,
            task="MAIN",
        )
        self.threads = config["threads"]

        self.set_logger()
        self.runtime = 0
        self.start_time = time.perf_counter()
        self.exec_time = 0

        # activity log

        self.qc_performed: bool = False
        self.enrichment_performed: bool = False
        self.depletion_performed: bool = False
        self.assembly_performed: bool = False
        self.read_classification_performed: bool = False
        self.contig_classification_performed: bool = False
        self.read_metagenomics_classification_performed: bool = False
        self.remap_prepped: bool = False
        self.remapping_performed: bool = False

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
            self.static_dir,
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
            os.path.join(self.static_dir_plots),
            exist_ok=True,
        )

        os.makedirs(
            self.media_dir_igv,
            exist_ok=True,
        )

        self.filtered_reads_dir = config["directories"][
            CS.PIPELINE_NAME_read_quality_analysis
        ]

        ######### INPUT
        self.sample_registered = config["sample_registered"]
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
            prefix=self.prefix,
            ps_pk=self.ps_pk,
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
            prefix=self.prefix,
            ps_pk=self.ps_pk,
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

        remap_params = TelevirParameters.get_remap_software(self.project_pk)

        self.set_metadata_tool()

        self.max_remap = remap_params.max_accids
        self.taxid_limit = remap_params.max_taxids
        self.remap_params = remap_params

        ### set software methods and actions

        self.set_methods(config, method_args)
        self.check_software_print()

        ### set default actions
        self.subsample = False
        self.sift = config["actions"]["SIFT"]

        self.classification = self.contig_classification or self.read_classification
        # self.remapping = config["actions"]["REMAP"]

        self.house_cleaning = config["actions"]["CLEAN"]

        ### drones
        self.depletion_drone = Classifier(
            SoftwareDetail("NONE", method_args, config, self.prefix),
            logging_level=self.logger_level_detail,
            log_dir=self.log_dir,
            prefix="drone",
        )
        self.enrichment_drone = Classifier(
            SoftwareDetail("NONE", method_args, config, self.prefix),
            logging_level=self.logger_level_detail,
            log_dir=self.log_dir,
            prefix="drone",
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

    def set_metadata_tool(self):
        self.metadata_tool = RunMetadataHandler(
            self.username,
            self.config,
            sift_query=self.config["sift_query"],
            prefix=self.prefix,
            rundir=self.deployment_dir,
        )

    def update_reads(self):
        self.sample.r1.update(
            self.prefix,
            clean_dir=self.config["directories"][
                CS.PIPELINE_NAME_read_quality_analysis
            ],
            enriched_dir=self.config["directories"]["reads_enriched_dir"],
            depleted_dir=self.config["directories"]["reads_depleted_dir"],
        )

        self.sample.r2.update(
            self.prefix,
            clean_dir=self.config["directories"][
                CS.PIPELINE_NAME_read_quality_analysis
            ],
            enriched_dir=self.config["directories"]["reads_enriched_dir"],
            depleted_dir=self.config["directories"]["reads_depleted_dir"],
        )

    def Update(self, config: dict, method_args: pd.DataFrame):
        self.method_args = pd.concat((self.method_args, method_args))
        # with open(config_json) as json_file:
        #    config = json.load(json_file)

        self.config = config
        self.prefix = config["prefix"]
        # self.type = config["type"]
        self.logger = logging.getLogger("{}".format(self.prefix))
        self.logger.setLevel(self.logger_level_main)
        logFormatter = logging.Formatter(fmt="{} :%(message)s".format(self.prefix))
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        self.logger.addHandler(consoleHandler)
        self.logger.info(f"prefix: {self.prefix}")
        self.logger.info(f"type: {self.type}")

        # directories
        self.filtered_reads_dir = config["directories"][
            CS.PIPELINE_NAME_read_quality_analysis
        ]
        self.log_dir = config["directories"]["log_dir"]

        ### set software methods and actions

        self.update_methods(config, method_args)

        ### set default actions
        self.sift = config["actions"]["SIFT"]
        self.house_cleaning = config["actions"]["CLEAN"]

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

    def update_merged_targets(self, targets_list: List[Remap_Target]):
        """
        Update the merged classification summary file.
        """

        self.metadata_tool.remap_targets = targets_list

    def Update_exec_time(self):
        """
        Update the execution time of the pipeline.
        """
        self.exec_time = time.perf_counter() - self.start_time


class Run_Deployment_Methods(RunDetail_main):
    def __init__(self, config_dict: dict, method_args: pd.DataFrame, project_pk: int):
        super().__init__(config_dict, method_args, project_pk)
        self.mapped_instances = []

    def Prep_deploy(self, remap_prep=True):
        if self.qc_performed is False:
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

        if self.depletion_performed is False:
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

        if self.enrichment_performed is False:
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

        if self.assembly_performed is False:
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

        if self.contig_classification_performed is False:
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

        if self.read_classification_performed is False:
            print("preparing read classification")
            print(self.sample.r1.current)
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

        if remap_prep and self.remap_prepped is False:
            self.prep_REMAPPING()

    def deploy_QC(self, fake_run: bool = False):
        self.logger.info(f"r1 reads: {self.sample.r1.get_current_fastq_read_number()}")
        self.logger.info(f"r2 reads: {self.sample.r2.get_current_fastq_read_number()}")

        if fake_run:
            self.preprocess_drone.fake_run()
        else:
            self.update_reads()

            if self.preprocess_method.check_processed_exist():
                r1_proc, r2_proc = self.preprocess_method.retrieve_qc_reads()
                self.sample.r1.clean_exo = r1_proc
                if self.type == ConstantsSettings.PAIR_END:
                    self.sample.r2.clean_exo = r2_proc

            self.preprocess_drone.run()

    def deploy_HD(self):
        self.update_reads()
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
        self.update_reads()

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

    def deploy_METAGENOMICS_CLASSIFICATION_reads(self):
        self.metagenomics_classification_drone = Classifier(
            self.metagenomics_classification_method,
            self.sample.r1.current,
            r2="",
            prefix=self.prefix,
            threads=self.threads,
            bin=get_bindir_from_binaries(
                self.config["bin"], CS.PIPELINE_NAME_remapping
            ),
            logging_level=self.logger_level_detail,
            log_dir=self.log_dir,
        )
        self.metagenomics_classification_drone.run()

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

    def prep_REMAPPING(self):

        self.remap_manager = Mapping_Manager(
            self.metadata_tool.remap_targets,
            self.sample.r1,
            self.sample.r2,
            self.software_remap,
            self.assembly_drone.assembly_file_fasta_gz,
            self.type,
            self.prefix,
            self.threads,
            get_bindir_from_binaries(self.config["bin"], CS.PIPELINE_NAME_remapping),
            self.logger_level_detail,
            True,
            remap_params=self.remap_params,
            logdir=self.log_dir,
        )

        self.remap_prepped = True

    def deploy_REMAPPING(self):

        print(
            f"{self.prefix} remapping # targets: {len(self.metadata_tool.remap_targets)}"
        )

        self.remap_manager.run_mappings_move_clean(
            self.static_dir_plots, self.media_dir_igv
        )

        self.remap_manager.merge_mapping_reports()
        self.remap_manager.collect_final_report_summary_statistics()


class RunEngine_class(Run_Deployment_Methods):
    def __init__(
        self, config_json: os.PathLike, method_args: pd.DataFrame, project_pk: int
    ):
        super().__init__(config_json, method_args, project_pk)

        self.logger.info("Starting Pipeline")

        self.logger.info(f"quality control: {self.quality_control}")
        self.logger.info(f"enrichment: {self.enrichment}")
        self.logger.info(f"depletion: {self.depletion}")
        self.logger.info(f"assembly: {self.assembly}")
        self.logger.info(f"classification: {self.classification}")
        self.logger.info(f"sift: {self.sift}")
        self.logger.info(f"remapping: {self.remapping}")
        self.logger.info(f"current reads: {self.sample.r1.current}")

    def Run_Full_Pipeline(self):
        self.Prep_deploy()
        self.Run_QC()
        self.Run_PreProcess()
        self.Sanitize_reads()
        self.Run_Assembly()
        self.Run_Classification()
        self.Run_Remapping()

    def process_QC(self):
        self.sample.input_fastqc_report = self.preprocess_drone.input_qc_report
        self.sample.processed_fastqc_report = self.preprocess_drone.processed_qc_report

        self.sample.r1.is_clean()
        self.sample.r2.is_clean()
        self.sample.reads_after_processing = self.sample.current_total_read_number()
        self.sample.get_fake_qc_data()
        self.sample.r1.clean_read_names()
        self.sample.r2.clean_read_names()

    def Run_QC(self):
        if self.quality_control:
            print("Deploying QC")
            self.deploy_QC()

            self.sample.qc_soft = self.preprocess_drone.preprocess_method.name
            self.sample.input_fastqc_report = self.preprocess_drone.input_qc_report
            self.sample.processed_fastqc_report = (
                self.preprocess_drone.processed_qc_report
            )

            self.process_QC()

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

            self.process_QC()

        self.Update_exec_time()
        self.generate_output_data_classes()

    def Run_PreProcess(self):
        self.logger.info(
            "r1 current reads: " + str(self.sample.r1.get_current_fastq_read_number())
        )
        if self.enrichment:
            self.deploy_EN()

            self.sample.r1.enrich(self.enrichment_drone.classified_reads_list)
            self.sample.r2.enrich(self.enrichment_drone.classified_reads_list)

            self.logger.info(
                "r1 current reads after enrichment: "
                + str(self.sample.r1.get_current_fastq_read_number())
            )

        if self.depletion:
            self.deploy_HD()

            self.logger.info(
                f"depleted reads: {len(self.depletion_drone.classified_reads_list)}"
            )

            self.sample.r1.deplete(self.depletion_drone.classified_reads_list)
            self.sample.r2.deplete(self.depletion_drone.classified_reads_list)

            self.logger.info(
                "r1 current reads after depletion: "
                + str(self.sample.r1.get_current_fastq_read_number())
            )

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
                self.remap_params.max_accids,
                self.remap_params.max_taxids,
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

    def export_assembly(self):
        self.assembly_drone.export_assembly(self.media_dir)

    def Summarize(self):
        self.logger.info(f"prefix: {self.prefix}")
        with open(os.path.join(self.log_dir, self.prefix + "_latest.fofn"), "w") as f:
            f.write(self.sample.r1.current + "\n")
            if self.type == ConstantsSettings.PAIR_END:
                f.write(self.sample.r2.current + "\n")

        with open(os.path.join(self.log_dir, "reads_latest.stats"), "w") as f:
            f.write(f"CLEAN\t{self.sample.r1.read_number_clean}\n")
            f.write(f"ENRICHED\t{self.sample.r1.read_number_enriched}\n")

    def generate_output_data_classes(self):
        ### transfer to sample class
        input_reads = self.sample.reads_before_processing
        processed_reads = self.sample.reads_after_processing

        filtered_reads = (
            self.sample.r1.read_number_filtered + self.sample.r2.read_number_filtered
        )

        final_processing_reads = (
            self.sample.r1.current_fastq_read_number()
            + self.sample.r2.current_fastq_read_number()
        )

        filtered_reads_perc = (int(filtered_reads) / input_reads) * 100
        final_processing_percent = (final_processing_reads / input_reads) * 100

        ### transfer to assembly class / drone.
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

        # enriched_reads = len(self.enrichment_drone.classified_reads_list)
        # depleted_reads = len(self.depletion_drone.classified_reads_list)
        enriched_reads = (
            self.sample.r1.enriched_read_number + self.sample.r2.enriched_read_number
        )
        depleted_reads = (
            self.sample.r1.depleted_read_number + self.sample.r2.depleted_read_number
        )
        # if self.type == ConstantsSettings.PAIR_END:
        #    enriched_reads = enriched_reads * 2
        #    depleted_reads = depleted_reads * 2

        self.run_detail_report = Run_detail_report(
            self.remap_manager.max_depth,
            self.remap_manager.max_depthR,
            self.remap_manager.max_gaps,
            self.remap_manager.max_prop,
            self.remap_manager.max_mapped,
            f"{processed_reads:,}",
            enriched_reads,
            enriched_reads / input_reads,
            depleted_reads,
            depleted_reads / input_reads,
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

        self.qc_report = RunQC_report(
            performed=self.quality_control,
            method=self.preprocess_drone.preprocess_method.name,
            args=self.preprocess_drone.preprocess_method.args,
            input_reads=self.sample.reads_before_processing,
            output_reads=self.sample.reads_after_processing,
            output_reads_percent=self.sample.reads_after_processing
            / self.sample.reads_before_processing,
        )

        self.contig_classification_results = Contig_classification_results(
            (
                True
                if self.contig_classification_drone.classifier_method.name != "None"
                else False
            ),
            self.contig_classification_drone.classifier_method.name,
            self.contig_classification_drone.classifier_method.args,
            self.contig_classification_drone.classifier_method.db_name,
            self.aclass_summary.shape[0],
            minhit_assembly,
            self.aclass_summary.shape[0] > 0,
        )

        self.read_classification_results = Read_classification_results(
            (
                True
                if self.read_classification_drone.classifier_method.name != "None"
                else False
            ),
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


class RunMainTree_class(Run_Deployment_Methods):
    def __init__(self, config_dict: dict, method_args: pd.DataFrame, project_pk: int):
        super().__init__(config_dict, method_args, project_pk)

        self.logger.info("Starting Pipeline")

        self.logger.info(f"quality control: {self.quality_control}")
        self.logger.info(f"enrichment: {self.enrichment}")
        self.logger.info(f"depletion: {self.depletion}")
        self.logger.info(f"assembly: {self.assembly}")
        self.logger.info(f"classification: {self.classification}")
        self.logger.info(f"sift: {self.sift}")
        self.logger.info(f"remapping: {self.remapping}")
        self.logger.info(f"current reads: {self.sample.r1.current}")

    def Run_Full_Pipeline(self):
        self.Prep_deploy()
        self.Run_QC()
        self.Run_PreProcess()
        self.Sanitize_reads()
        self.Run_Assembly()
        self.Run_Classification()
        self.Run_Remapping()

    def process_QC(self):
        self.sample.input_fastqc_report = self.preprocess_drone.input_qc_report
        self.sample.processed_fastqc_report = self.preprocess_drone.processed_qc_report

        self.sample.r1.is_clean()
        self.sample.r2.is_clean()
        self.sample.reads_after_processing = self.sample.current_total_read_number()
        self.sample.get_fake_qc_data()
        self.sample.r1.clean_read_names()
        self.sample.r2.clean_read_names()

    def Run_QC(self):
        if self.quality_control and not self.qc_performed:
            self.deploy_QC()
            self.sample.qc_soft = self.preprocess_drone.preprocess_method.name

            self.process_QC()

            self.qc_performed = True

        elif (
            self.sample.r1.current_status == "raw"
            and self.sample.r2.current_status == "raw"
        ):
            self.deploy_QC(fake_run=True)

            shutil.copy(self.sample.r1.current, self.sample.r1.clean)

            if self.sample.r2.exists:
                shutil.copy(self.sample.r2.current, self.sample.r2.clean)

            self.sample.qc_soft = "none"

            self.process_QC()

            # self.qc_performed = True

        self.Update_exec_time()
        self.generate_output_data_classes()

    def Run_PreProcess(self):
        self.logger.info(
            "r1 current reads: " + str(self.sample.r1.get_current_fastq_read_number())
        )
        print("RUNNING PREPROCESS", self.enrichment)

        if self.enrichment:
            print("ENRICHMENT EXISTS")
            print(self.enrichment_method.check_enriched_exist())

            if self.enrichment_method.check_enriched_exist():

                r1_proc, r2_proc = self.enrichment_method.retrieve_enriched_reads()
                enriched_read_number = self.enrichment_method.get_enriched_read_number()
                if self.type == ConstantsSettings.PAIR_END:
                    enriched_read_number = enriched_read_number / 2
                    self.sample.r2.enriched_read_number = enriched_read_number
                self.sample.r1.enriched_read_number = enriched_read_number

                self.sample.r1.enriched_exo = r1_proc
                shutil.copy(r1_proc, self.sample.r1.enriched)

                if self.type == ConstantsSettings.PAIR_END:
                    self.sample.r2.enriched_exo = r2_proc
                    shutil.copy(r2_proc, self.sample.r2.enriched)

                self.sample.r1.is_enriched()
                self.sample.r2.is_enriched()

            else:
                self.deploy_EN()

                self.sample.r1.enrich(self.enrichment_drone.classified_reads_list)
                self.sample.r2.enrich(self.enrichment_drone.classified_reads_list)
                self.enrichment_drone.deployed = True

            self.logger.info(
                "r1 current reads after enrichment: "
                + str(self.sample.r1.get_current_fastq_read_number())
            )

            self.enrichment_performed = True

        if self.depletion:
            print(self.enrichment_method.check_depleted_exist())
            print("DEPLETION EXISTS")
            if self.depletion_method.check_depleted_exist():
                r1_proc, r2_proc = self.depletion_method.retrieve_depleted_reads()
                depleted_read_number = self.depletion_method.get_depleted_read_number()
                if self.type == ConstantsSettings.PAIR_END:
                    depleted_read_number = depleted_read_number / 2
                    self.sample.r2.depleted_read_number = depleted_read_number
                self.sample.r1.depleted_read_number = depleted_read_number

                self.sample.r1.depleted_exo = r1_proc
                shutil.copy(r1_proc, self.sample.r1.depleted)

                if self.type == ConstantsSettings.PAIR_END:
                    self.sample.r2.depleted_exo = r2_proc
                    shutil.copy(r2_proc, self.sample.r2.depleted)

                self.sample.r1.is_depleted()
                self.sample.r2.is_depleted()
                self.depletion_drone.deployed = True

            else:
                self.deploy_HD()

                ###########################
                ###########################

                from pathogen_identification.utilities.televir_bioinf import (
                    TelevirBioinf,
                )

                # televir_bioinf = TelevirBioinf()
                # alignment_file = self.depletion_drone.classifier.report_path
                # self.depletion_report = televir_bioinf.alignment_agg_by_target(
                #    alignment_file
                # )
                ###########################
                ###########################

                self.sample.r1.deplete(self.depletion_drone.classified_reads_list)
                self.sample.r2.deplete(self.depletion_drone.classified_reads_list)

            self.depletion_performed = True

        self.Update_exec_time()
        self.generate_output_data_classes()

    def Sanitize_reads(self):
        if self.enrichment or self.depletion:
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
        if self.assembly and self.assembly_performed is False:
            self.deploy_ASSEMBLY()
            self.assembly_performed = True

        self.Update_exec_time()
        self.generate_output_data_classes()

    def Run_Classification(self):
        if self.classification:
            self.deploy_READ_CLASSIFICATION()
            self.deploy_CONTIG_CLASSIFICATION()

            self.metadata_tool.match_and_select_targets(
                self.read_classification_drone.classification_report,
                self.contig_classification_drone.classification_report,
                self.remap_params.max_accids,
                self.remap_params.max_taxids,
            )
            self.import_from_remap_prep()

            self.read_classification_performed = True
            self.contig_classification_performed = True

            self.export_intermediate_reports()

        self.Update_exec_time()
        self.generate_output_data_classes()

    def Prep_Metagenomics_Classification(self):
        reference_utils = RawReferenceUtils(self.sample_registered)

        ### ############################################################# ###
        # reference_table = reference_utils.sample_reference_tables()

        self.metadata_tool.merge_sample_references_ensemble(
            self.sample_registered, max_remap=1
        )

        self.prep_REMAPPING()
        self.remap_manager.generate_remap_targets_fasta()

    def Run_Metagenomics_Classification(self):
        if self.metagenomics_classification:
            if self.remapping_method.dir == "None":
                self.software_remap.set_output_dir(
                    self.metagenomics_classification_method.dir
                )

            self.Prep_Metagenomics_Classification()

            self.metagenomics_classification_method.set_db(
                self.remap_manager.combined_fasta_gz_path
            )

            self.deploy_METAGENOMICS_CLASSIFICATION_reads()
            self.read_classification_drone = self.metagenomics_classification_drone
            self.read_classification_performed = True
            self.read_metagenomics_classification_performed = True

            self.plan_remap_prep_safe()

        self.Update_exec_time()
        self.generate_output_data_classes()

    def Run_Contig_classification(self):
        """
        This is a special case where we only want to run contig classification"""

        if self.contig_classification and not self.contig_classification_performed:
            self.deploy_CONTIG_CLASSIFICATION()
            self.contig_classification_performed = True

        self.Update_exec_time()
        self.generate_output_data_classes()

    def Run_Read_classification(self):
        """
        This is a special case where we only want to run read classification"""

        if self.read_classification and not self.read_classification_performed:
            self.deploy_READ_CLASSIFICATION()
            self.read_classification_performed = True

        self.Update_exec_time()
        self.generate_output_data_classes()

    def plan_remap_prep_safe(self):
        self.plan_remap_prep()
        # self.export_intermediate_reports()
        self.remap_prepped = True

    def plan_remap_prep(self):

        print("########### PLANNING REMAP PREP ###########")
        print("remap prep: ", self.remap_prepped)

        print("MAX TAXIDS: ", self.remap_params.max_taxids)
        print("MAX ACCIDS: ", self.remap_params.max_accids)

        self.metadata_tool.match_and_select_targets(
            self.read_classification_drone.classification_report,
            self.contig_classification_drone.classification_report,
            max_remap=self.remap_params.max_accids,
            taxid_limit=self.remap_params.max_taxids,
        )

        print("remap targets: ", self.metadata_tool.remap_targets)

        self.import_from_remap_prep()

    def import_from_remap_prep(self):
        self.aclass_summary = self.metadata_tool.aclass
        self.rclass_summary = self.metadata_tool.rclass
        self.merged_targets = self.metadata_tool.merged_targets

        print("RAW TARGETS: ", self.metadata_tool.raw_targets)
        self.raw_targets = self.metadata_tool.raw_targets
        self.remap_plan = self.metadata_tool.remap_plan

    def plan_combined_remapping(self):
        self.metadata_tool.merge_sample_references_ensemble(
            self.sample_registered,
            max_taxids=self.remap_params.max_taxids,
            max_remap=self.remap_params.max_accids,
        )

        self.import_from_remap_prep()
        self.generate_output_data_classes()
        self.prep_REMAPPING()
        self.remap_prepped = True

    def Run_Remapping(self, prep=True):
        print("remapping: ", self.remapping)
        print(self.remap_prepped)
        if not self.remap_prepped:
            return

        print("remapping: ", self.remapping)
        print(self.remap_prepped)
        if self.remapping is True and self.remapping_performed is False:
            if self.remap_prepped == False:
                self.plan_remap_prep_safe()

            print("merged targets: ", self.merged_targets)

            if prep:
                self.prep_REMAPPING()

            self.deploy_REMAPPING()

            self.remapping_performed = True

        self.Update_exec_time()
        self.generate_output_data_classes()

    #### SUMMARY FUNCTIONS ####

    def update_mapped_instances(self, mapped_instance: List[Mapping_Instance]):
        """Update the remap manager with the new mapped instances, register."""
        self.prep_REMAPPING()
        self.remap_manager.update_mapped_instances(mapped_instance)

        self.remapping_performed = True
        self.generate_output_data_classes()
        self.Update_exec_time()

    def export_logdir(self):
        if os.path.exists(self.media_dir_logdir):
            shutil.rmtree(self.media_dir_logdir)

        shutil.copytree(
            self.log_dir,
            self.media_dir_logdir,
        )

    def export_final_reports(self):
        # main report
        self.report.to_csv(
            self.full_report,
            index=False,
            sep="\t",
            header=True,
        )

    def save_df_check_exists(self, df: pd.DataFrame, path: str):
        dirname = os.path.dirname(path)
        print("############### SAVING ", path)
        print("path exists: ", os.path.exists(path))
        print("dirname exists: ", os.path.exists(dirname))
        print(df.shape)

        if not os.path.exists(dirname):
            os.makedirs(dirname, exist_ok=True)

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

    def export_assembly(self):
        self.assembly_drone.export_assembly(self.media_dir)

    def Summarize(self):
        self.logger.info(f"prefix: {self.prefix}")
        with open(os.path.join(self.log_dir, self.prefix + "_latest.fofn"), "w") as f:
            f.write(self.sample.r1.current + "\n")
            if self.type == ConstantsSettings.PAIR_END:
                f.write(self.sample.r2.current + "\n")

        with open(os.path.join(self.log_dir, "reads_latest.stats"), "w") as f:
            f.write(f"CLEAN\t{self.sample.r1.read_number_clean}\n")
            f.write(f"ENRICHED\t{self.sample.r1.read_number_enriched}\n")

    def generate_output_data_classes(self):
        # merge mapping results if exist.
        #
        self.remap_manager.merge_mapping_reports()
        self.remap_manager.collect_final_report_summary_statistics()
        self.report = self.remap_manager.report
        # transfer to sample class
        processed_reads = self.sample.reads_before_processing

        filtered_reads = (
            self.sample.r1.read_number_filtered + self.sample.r2.read_number_filtered
        )

        final_processing_reads = (
            self.sample.r1.current_fastq_read_number()
            + self.sample.r2.current_fastq_read_number()
        )

        if processed_reads == 0:
            filtered_reads_perc = 0
            final_processing_percent = 0
            processed_reads = 1
        else:
            filtered_reads_perc = (int(filtered_reads) / processed_reads) * 100
            final_processing_percent = (final_processing_reads / processed_reads) * 100

        # transfer to assembly class / drone.

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

        enriched_reads = (
            self.sample.r1.read_number_enriched + self.sample.r2.read_number_enriched
        )
        depleted_reads = (
            self.sample.r1.depleted_read_number + self.sample.r2.depleted_read_number
        )

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

        self.qc_report = RunQC_report(
            performed=self.qc_performed,
            method=self.preprocess_drone.preprocess_method.name,
            args=self.preprocess_drone.preprocess_method.args,
            input_reads=self.sample.reads_before_processing,
            output_reads=self.sample.reads_after_processing,
            output_reads_percent=self.sample.reads_after_processing
            / self.sample.reads_before_processing,
        )

        self.contig_classification_results = Contig_classification_results(
            (
                True
                if self.contig_classification_drone.classifier_method.name != "None"
                else False
            ),
            self.contig_classification_drone.classifier_method.name,
            self.contig_classification_drone.classifier_method.args,
            self.contig_classification_drone.classifier_method.db_name,
            self.aclass_summary.shape[0],
            minhit_assembly,
            self.aclass_summary.shape[0] > 0,
        )

        self.read_classification_results = Read_classification_results(
            (
                True
                if self.read_classification_drone.classifier_method.name != "None"
                else False
            ),
            self.read_classification_drone.classifier_method.name,
            self.read_classification_drone.classifier_method.args,
            self.read_classification_drone.classifier_method.db_name,
            self.rclass_summary.shape[0],
            minhit_reads,
            self.rclass_summary.shape[0] > 0,
        )

        self.assembly_report = Assembly_results(
            self.assembly_performed,
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
