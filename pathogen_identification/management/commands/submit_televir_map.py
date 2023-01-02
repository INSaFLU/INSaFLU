import argparse
import logging
import os
import shutil

import pandas as pd
from django.core.management.base import BaseCommand
from managing_files.models import ProcessControler
from pathogen_identification.constants_settings import MEDIA_ROOT, ConstantsSettings
from constants.constants import Televir_Metadata_Constants as Televir_Metadata
from pathogen_identification.install_registry import Params_Illumina, Params_Nanopore
from pathogen_identification.models import RawReference
from pathogen_identification.modules.metadata_handler import Metadata_handler
from pathogen_identification.modules.object_classes import (
    Read_class,
    Sample_runClass,
    Software_detail,
)
from pathogen_identification.modules.remap_class import (
    Mapping_Instance,
    Mapping_Manager,
)
from pathogen_identification.utilities.update_DBs import (
    Update_FinalReport,
    Update_ReferenceMap,
)
from pathogen_identification.utilities.utilities_general import simplify_name
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager
from settings.constants_settings import ConstantsSettings as CS
from utils.process_SGE import ProcessSGE


class RunMain:

    remap_manager: Mapping_Manager
    mapping_instance: Mapping_Instance
    metadata_tool: Metadata_handler

    ##  metadata
    sift_query: str
    max_remap: int
    taxid_limit: int

    def __init__(self, config: dict, method_args: pd.DataFrame):

        self.sample_name = config["sample_name"]
        self.type = config["type"]
        self.project_name = "none"
        self.username = "none"
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
        self.remapping_method = Software_detail(
            CS.PIPELINE_NAME_remapping,
            method_args,
            config,
            self.prefix,
        )

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
            self.remapping_method,
            "Dummy",
            self.type,
            self.prefix,
            self.threads,
            self.minimum_coverage,
            get_bindir_from_binaries(self.config["bin"], CS.PIPELINE_NAME_remapping),
            self.logger_level_detail,
            self.house_cleaning,
            logdir=self.config["directories"]["log_dir"],
        )

        self.logger.info(
            f"{self.prefix} remapping # targets: {len(self.metadata_tool.remap_targets)}"
        )

        self.remap_manager.run_mappings()
        self.remap_manager.merge_mapping_reports()
        self.remap_manager.collect_final_report_summary_statistics()

    def run(self):
        self.deploy_REMAPPING()
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

        self.technology = reference.run.project.technology
        self.threads = threads
        self.prefix = reference.accid
        self.project = reference.run.project.name
        self.clean = False
        self.deployment_root_dir = os.path.join(output_dir, self.technology)
        self.dir_branch = os.path.join(self.deployment_root_dir, self.prefix)
        self.dir = self.dir_branch

        os.makedirs(self.dir, exist_ok=True)

        self.r1_path = reference.run.sample.sample.path_name_1.path
        self.r2_path = (
            reference.run.sample.sample.path_name_2.path
            if reference.run.sample.sample.exist_file_2()
            else ""
        )

        self.taxid = reference.taxid
        self.accid = reference.accid

        if self.technology == "ONT":
            self.params = Params_Nanopore
        else:
            self.params = Params_Illumina

    def input_read_project_path(self, filepath):
        if not os.path.isfile(filepath):
            return ""
        rname = os.path.basename(filepath)
        new_rpath = os.path.join(self.dir, "reads") + "/" + rname
        shutil.copy(filepath, new_rpath)
        return new_rpath

    def generate_method_args(self):

        parameter_leaf = self.reference.run.parameter_set.leaf
        run_df = self.utils.get_leaf_parameters(parameter_leaf)

        self.method_args = run_df[run_df.module == CS.PIPELINE_NAME_remapping]

        if self.method_args.empty:
            raise ValueError(
                f"no remapping parameters found for {self.reference.accid} in leaf {parameter_leaf}"
            )

    def generate_config(self):

        self.config = {
            "sample_name": simplify_name(
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
            self.config["directories"][dr] = os.path.join(self.dir_branch, g)
            os.makedirs(self.config["directories"][dr], exist_ok=True)

        self.config["r1"] = self.input_read_project_path(self.r1_path)
        self.config["r2"] = self.input_read_project_path(self.r2_path)
        self.config["type"] = ["SE", "PE"][int(os.path.isfile(self.config["r2"]))]

        self.config.update(self.params.CONSTANTS)

    def update_raw_reference_status_mapped(self):
        self.reference.status = RawReference.STATUS_MAPPED
        self.reference.save()

    def update_raw_reference_status_fail(self):
        self.reference.status = RawReference.STATUS_FAIL
        self.reference.save()

    def update_final_report(self, run_class: RunMain):

        run = self.reference.run
        sample = run.sample

        Update_FinalReport(run_class, run, sample)

        for ref_map in run_class.remap_manager.mapped_instances:

            Update_ReferenceMap(ref_map, run, sample)


class Command(BaseCommand):
    help = "deploy run"

    def add_arguments(self, parser):
        parser.add_argument(
            "--ref_id",
            type=int,
            help="user deploying the run (pk)",
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

        reference = RawReference.objects.get(pk=raw_reference_id)
        user = reference.run.project.owner

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

            run_engine = RunMain(input_generator.config, input_generator.method_args)
            run_engine.generate_targets()
            run_engine.run()

            input_generator.update_raw_reference_status_mapped()
            input_generator.update_final_report(run_engine)

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
