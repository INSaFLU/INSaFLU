import argparse
import logging
import os
import shutil

import pandas as pd
from pathogen_identification.constants_settings import (MEDIA_ROOT,
                                                        ConstantsSettings)
from pathogen_identification.install_registry import Deployment_Params as DP
from pathogen_identification.install_registry import (Params_Illumina,
                                                      Params_Nanopore)
from pathogen_identification.modules.metadata_handler import Metadata_handler
from pathogen_identification.modules.object_classes import (Read_class,
                                                            Remap_Target,
                                                            RunCMD,
                                                            Sample_runClass,
                                                            Software_detail)
from pathogen_identification.modules.remap_class import (Mapping_Instance,
                                                         Mapping_Manager)
from pathogen_identification.utilities.utilities_general import simplify_name
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


class Input_Generator:
    def __init__(self):
        self.install_registry = DP
        args = self.parse_args()

        self.technology = args.technology
        self.threads = args.threads
        self.prefix = args.prefix
        self.project = args.prefix
        self.clean = False if args.keep else True
        self.deployment_root_dir = os.path.join(MEDIA_ROOT, self.technology)
        self.dir_branch = os.path.join(self.deployment_root_dir, self.prefix)
        self.dir = self.dir_branch

        os.makedirs(self.dir, exist_ok=True)
        print(self.dir)

        self.r1_path = args.r1
        self.r2_path = args.r2
        self.taxid = args.taxid
        self.accid = args.accid

        if args.technology == "ONT":
            self.params = Params_Nanopore
        else:
            self.params = Params_Illumina

    def parse_args(self):
        args = argparse.ArgumentParser()
        args.add_argument("--taxid", type=str, required=False, default="none")
        args.add_argument("--accid", type=str, required=False, default="none")
        args.add_argument("--r1", type=str, required=True)
        args.add_argument("--r2", type=str, required=False, default="")
        args.add_argument("--threads", type=int, required=False, default=1)

        args.add_argument("--prefix", type=str, required=False, default="test")
        args.add_argument("--technology", type=str, required=False, default="ONT")
        args.add_argument("--keep", required=False, default=True, action="store_false")

        return args.parse_args()

    def input_read_project_path(self, filepath):
        if not os.path.isfile(filepath):
            return ""
        rname = os.path.basename(filepath)
        new_rpath = os.path.join(self.dir, "reads") + "/" + rname
        shutil.copy(filepath, new_rpath)
        return new_rpath

    def generate_method_args(self):
        self.method_args = pd.DataFrame(
            {
                "software": "snippy",
                "module": CS.PIPELINE_NAME_remapping,
                "parameter": "SNIPPY_ARGS",
                "value": "--cpus 1 --outdir snippy --ref {reference} --R1 {r1} --R2 {r2}",
                "description": "Snippy arguments",
            },
            index=[0],
        )

        self.method_args = pd.DataFrame(
            {
                "software": "minimap2",
                "module": CS.PIPELINE_NAME_remapping,
                "parameter": "MINIMAP2_ARGS",
                "value": "",
            },
            index=[0],
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

        print(self.r1_path)
        self.config["r1"] = self.input_read_project_path(self.r1_path)
        self.config["r2"] = self.input_read_project_path(self.r2_path)
        self.config["type"] = ["SE", "PE"][int(os.path.isfile(self.config["r2"]))]

        self.config.update(self.params.CONSTANTS)


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
        self.house_cleaning = True
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

        print(result_df)

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

        if self.clean:
            self.remap_manager.clean_mapping_files()

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


if __name__ == "__main__":

    import argparse

    import pandas as pd

    input_generator = Input_Generator()

    input_generator.generate_method_args()
    input_generator.generate_config()

    run_engine = RunMain(input_generator.config, input_generator.method_args)
    run_engine.generate_targets()
    run_engine.run()

    print("done")
