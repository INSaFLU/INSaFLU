import itertools as it
import logging
import os
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from random import randint
from typing import Type

import matplotlib
import pandas as pd
from numpy import ERR_CALL
from pathogen_identification.utilities.utilities_general import fastqc_parse

matplotlib.use("Agg")
import gzip
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from constants.constants import Televir_Metadata_Constants as Televir_Metadata


class Temp_File:
    """
    Temporary file.
    """

    def __init__(self, temp_dir: str, prefix: str = "temp", suffix: str = ""):
        """
        Initialize.
        """

        self.temp_dir = temp_dir
        self.prefix = prefix
        self.suffix = suffix

        self.temp_file = os.path.join(
            self.temp_dir, f"{self.prefix}_{randint(1000000, 9999999)}{self.suffix}"
        )

    def __enter__(self):
        """
        Enter.
        """

        return self.temp_file

    def __exit__(self):
        """
        Exit.
        """

        if os.path.exists(self.temp_file):
            os.remove(self.temp_file)


class Operation_Temp_Files:
    """
    Operation on temporary files.
    """

    def __init__(self, temp_dir: str, prefix: str = "temp"):
        """
        Initialize.
        """

        self.temp_dir = temp_dir
        self.prefix = prefix
        seed = randint(1000000, 9999999)

        self.script = os.path.join(self.temp_dir, f"{self.prefix}_{seed}.sh")

        self.log = os.path.join(self.temp_dir, f"{self.prefix}_{seed}.log")

        self.flag = os.path.join(self.temp_dir, f"{self.prefix}_{seed}.flag")

    def __enter__(self):
        """
        Enter.
        """

        return self

    def __exit__(
        self,
        exc_type: Type[BaseException],
        exc_value: BaseException,
        traceback,
    ):
        """
        Exit.
        """

        if os.path.exists(self.script):
            os.remove(self.script)

        if os.path.exists(self.log):
            os.remove(self.log)

        if os.path.exists(self.flag):
            os.remove(self.flag)

    def write_bash_script(self, cmd: str):
        """
        Write bash script.
        """

        with open(self.script, "w") as f:
            f.write("#!/bin/bash")
            f.write("\n")
            f.write(cmd)
            f.write("\n")
            f.write("touch " + self.flag)

    def run_bash_script(self):
        """
        Run bash script.
        """

        start_time = time.perf_counter()

        proc_prep = subprocess.Popen(
            "bash " + self.script + " &> " + self.log,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        out, err = proc_prep.communicate()

        found_flag = False

        while not found_flag:
            time.sleep(1)
            found_flag = os.path.exists(self.flag)

        err = open(self.log).read()

        exec_time = time.perf_counter() - start_time

        return out, err, exec_time


class RunCMD:
    """
    Run command line commands.
    """

    def __init__(self, bin, logdir: str = "", prefix: str = "run", task: str = "NONE"):
        """
        Initialize.
        """
        if bin:
            if bin[-1] != "/":
                bin += "/"

        self.bin = bin
        self.logs = []
        self.task = task

        self.logger = logging.getLogger(f"{prefix}_{task}")
        self.logger.setLevel(logging.ERROR)
        self.logger.addHandler(logging.StreamHandler())
        self.logger.propagate = False

        self.logfile = os.path.join(logdir, f"{prefix}_{task}.log")
        self.logdir = logdir
        self.prefix = prefix

    def flag_error(self, subprocess_errorlog, cmd: str = ""):
        """
        Check if error in subprocess.
        """

        if subprocess_errorlog:
            try:
                subprocess_errorlog = subprocess_errorlog.decode("utf-8")
                if "Killed" in subprocess_errorlog:
                    print(f"Killed {cmd}")
                if "[error]" in subprocess_errorlog.lower():
                    return True

            except Exception as e:
                return False

        return False

    @staticmethod
    def process_cmd_log(cmd_out):
        """
        Process command log.
        """

        if not cmd_out:
            return ""

        if isinstance(cmd_out, bytes):
            try:  # python 3
                cmd_out = cmd_out.decode()
            except Exception as e:
                print("log error:", e)
                print("out:", cmd_out)
                cmd_out = ""

        if len(cmd_out) > 300:
            cmd_out = cmd_out[:70] + "..." + cmd_out[-70:].strip()

        return cmd_out

    def system_deploy(self, cmd: str):

        start_time = time.perf_counter()

        proc_prep = subprocess.Popen(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        out, err = proc_prep.communicate()

        exec_time = time.perf_counter() - start_time

        return out, err, exec_time

    def dispose_output_carefully(self, cmd, out, err, exec_time):
        """
        Dispose output carefully.
        """

        if self.flag_error(err, cmd):
            self.logger.error(f"errror in command: {self.bin}{cmd}")
            raise Exception(err.decode("utf-8"))

        self.output_disposal(cmd, err, out, exec_time, self.bin)

    def bash_cmd_string(self, cmd: str):
        return f"{cmd}"

    def bash_software_cmd_string(self, cmd: str):
        return f"{self.bin}{cmd}"

    def python_cmd_string(self, cmd: str):
        return f"python {self.bin}{cmd}"

    def java_cmd_string(self, cmd: str):
        java_bin = os.path.join(
            Televir_Metadata.BINARIES["ROOT"],
            Televir_Metadata.BINARIES["software"]["java"],
            "bin",
            "java",
        )

        return f"{java_bin} -cp {self.bin} {cmd}"

    def output_disposal(self, cmd: str, err: str, out: str, exec_time: float, bin: str):
        if self.logdir:
            with open(os.path.join(self.logdir, self.logfile), "a") as f:
                software = cmd.split(" ")[0]
                f.write(f"exec\t{software}\t{exec_time}\n")
                f.write(f"bin\t{bin}\n")
                f.write(f"{cmd}\n")

                out = self.process_cmd_log(out)
                err = self.process_cmd_log(err)

                f.write(f"out\t{out}\n")
                f.write(f"err\t{err}\n")

        else:
            self.logs.append(cmd)
            self.logs.append(out)

    def run(self, cmd):
        """
        Run software.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        self.logger.info(f"running: {self.bin}{cmd}")

        cmd_string = self.bash_software_cmd_string(cmd)
        out, err, exec_time = self.system_deploy(cmd_string)
        self.dispose_output_carefully(cmd, out, err, exec_time)

    def run_python(self, cmd):
        """
        Run software with python.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        self.logger.info(f"running command: python {self.bin}{cmd}")

        cmd_string = self.python_cmd_string(cmd)
        out, err, exec_time = self.system_deploy(cmd_string)
        self.dispose_output_carefully(cmd, out, err, exec_time)

    def run_java(self, cmd):
        """
        Run software with java.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        self.logger.info(f"running command: java {self.bin}{cmd}")

        cmd_string = self.java_cmd_string(cmd)
        out, err, exec_time = self.system_deploy(cmd_string)
        self.dispose_output_carefully(cmd, out, err, exec_time)

    def run_bash(self, cmd):
        """
        Run bash command.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        self.logger.info(f"running command: {cmd}")

        cmd_string = self.bash_cmd_string(cmd)
        out, err, exec_time = self.system_deploy(cmd_string)
        self.dispose_output_carefully(cmd, out, err, exec_time)

    def run_bash_return(self, cmd):
        """
        Run bash command and return stdout.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        self.logger.info(f"running command: {cmd}")

        cmd_string = self.bash_cmd_string(cmd)
        out, err, exec_time = self.system_deploy(cmd_string)
        self.dispose_output_carefully(cmd, out, err, exec_time)

        return out

    def run_script_software(self, cmd):
        """
        Run bash script.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        operation_files = Operation_Temp_Files(self.logdir)

        with operation_files as op_files:

            cmd_string = self.bash_software_cmd_string(cmd)

            op_files.write_bash_script(cmd_string)

            out, err, exec_time = op_files.run_bash_script()

            if self.flag_error(err):
                self.logger.error(f"errror in command: {self.bin}{cmd}")

        self.output_disposal(cmd, err, out, exec_time, "")

    def run_script(self, cmd):
        """
        Run bash script.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        operation_files = Operation_Temp_Files(self.logdir)

        with operation_files as op_files:

            cmd_string = self.bash_cmd_string(cmd)

            op_files.write_bash_script(cmd_string)

            out, err, exec_time = op_files.run_bash_script()

            if self.flag_error(err):
                self.logger.error(f"errror in command: {self.bin}{cmd}")

        self.output_disposal(cmd, err, out, exec_time, "")

    def run_script_return(self, cmd):
        """
        Run bash script.
        """

        if isinstance(cmd, list):
            cmd = " ".join(cmd)

        operation_files = Operation_Temp_Files(self.logdir)

        with operation_files as op_files:

            cmd_string = self.bash_cmd_string(cmd)

            op_files.write_bash_script(cmd_string)

            out, err, exec_time = op_files.run_bash_script()

            if self.flag_error(err):
                self.logger.error(f"errror in command: {self.bin}{cmd}")

        self.output_disposal(cmd, err, out, exec_time, "")

        return out.strip()


class Read_class:
    def __init__(
        self, filepath, clean_dir: str, enriched_dir: str, depleted_dir: str, bin: str
    ):
        """
        Initialize.

        Args:
            filepath: path to file.
            clean_dir: path to clean directory.
            enriched_dir: path to enriched directory.
            depleted_dir: path to depleted directory.
            bin: path to bin directory.

        """
        self.cmd = RunCMD(bin, prefix="read", task="housekeeping")

        self.exists = os.path.isfile(filepath)

        self.filepath = filepath
        self.current = filepath
        self.prefix = self.determine_read_name(filepath)
        self.clean = os.path.join(clean_dir, self.prefix + ".clean.fastq.gz")
        self.enriched = os.path.join(enriched_dir, self.prefix + ".enriched.fastq.gz")
        self.depleted = os.path.join(depleted_dir, self.prefix + ".depleted.fastq.gz")
        self.current_status = "raw"
        self.read_number_raw = self.get_current_fastq_read_number()
        self.read_number_clean = 0
        self.read_number_enriched = 0
        self.read_number_depleted = 0
        self.read_number_filtered = 0

    def update(self, clean_dir: str, enriched_dir: str, depleted_dir: str):
        self.clean = os.path.join(clean_dir, self.prefix + ".clean.fastq.gz")
        self.enriched = os.path.join(enriched_dir, self.prefix + ".enriched.fastq.gz")
        self.depleted = os.path.join(depleted_dir, self.prefix + ".depleted.fastq.gz")

    def get_read_names_fastq(self, filepath):
        """
        Get read names from fastq file.
        """

        read_names = []

        if self.exists:
            with gzip.open(filepath, "rt") as f:
                for line in f:
                    if line.startswith("@"):
                        read_names.append(line.split()[0][1:])

        return read_names

    def determine_read_name(self, filepath):

        if not self.exists:
            return "none"

        if "gz" not in filepath:
            self.cmd.run("bgzip -f %s" % filepath)

        filename = os.path.basename(filepath)

        filename = filename.replace(".gz", "")

        filename = filename.replace(".fasta", "").replace(".fa", "")

        return filename

    def read_filter_move(self, input: str, read_list: list, output: str = ""):

        if not self.exists:
            return

        temp_reads_keep = os.path.join(
            os.path.dirname(output), f"keep_temp_{randint(1,1999)}.lst"
        )
        with open(temp_reads_keep, "w") as f:
            f.write("\n".join(read_list))

        self.read_filter(input, output, temp_reads_keep)

        # os.remove(temp_reads_keep)

    def read_filter_inplace(self, input: str, read_list: str):

        if not self.exists:
            return

        output_dir = os.path.dirname(input)
        tempreads = os.path.join(output_dir, f"temp_{randint(1,1999)}.fq.gz")

        self.read_filter(input, tempreads, read_list)

        if os.path.isfile(tempreads) and os.path.getsize(tempreads):
            os.remove(input)
            os.rename(tempreads, input)

    def read_filter(self, input: str, output: str, read_list: str):
        """
        filter read file using exisiting lsit of reads.
        Args:

            input: path to input file.
            output: path to output file.
            read_list: path to file containing read list.
        """
        if not self.exists:
            return

        cmd = "seqtk subseq %s %s | gzip > %s" % (input, read_list, output)

        self.cmd.run(cmd)

    def enrich(self, read_list):
        """
        filter reads and set current status to enriched.
        """

        if len(read_list) > 0:
            self.read_filter_move(self.current, read_list, self.enriched)
            self.is_enriched()

    def deplete(self, read_list):
        """
        filter reads and aset current status to depleted.
        """

        current_reads = self.get_read_names_fastq(self.current)

        read_list_to_keep = list(set(current_reads) - set(read_list))

        if len(read_list) > 0:
            self.read_filter_move(self.current, read_list_to_keep, self.depleted)
            self.is_depleted()

    def is_clean(self):
        """
        Set current status to clean.
        """
        self.current = self.clean
        self.read_number_clean = self.get_current_fastq_read_number()
        self.current_status = "clean"
        self.filepath = os.path.dirname(self.current)

    def is_enriched(self):
        """
        Set current status to enriched.
        """
        self.current = self.enriched
        self.read_number_enriched = self.get_current_fastq_read_number()
        self.current_status = "enriched"
        self.filepath = os.path.dirname(self.current)
        self.read_number_filtered = self.read_number_enriched

    def is_depleted(self):
        """
        Set current status to depleted.
        """
        self.current = self.depleted
        self.read_number_depleted = self.get_current_fastq_read_number()
        self.current_status = "depleted"
        self.filepath = os.path.dirname(self.current)
        self.read_number_filtered = self.read_number_depleted

    def get_current_fastq_read_number(self):
        """
        Get number of reads in current fastq file."""

        if not self.exists:
            return 0

        cmd = "zcat %s | wc -l" % self.current
        rnumber = self.cmd.run_bash_return(cmd).decode("utf-8")
        return int(rnumber) // 4

    def clean_read_names(self):
        """
        Clean read names in current fastq file.
        """

        if not self.exists:
            return

        temp_fq = os.path.join(
            os.path.dirname(self.current), f"temp_clean_{randint(1,1000)}.fastq"
        )
        temp_fq_gz = temp_fq + ".gz"

        cmd_unzip = "gunzip -c %s > %s" % (self.current, temp_fq)

        cmd = [
            "sed",
            "-i",
            "'/^[@+]/ s/\/[12]$//g'",
            temp_fq,
        ]

        cmd_zip = "bgzip %s" % temp_fq

        self.cmd.run_bash(cmd_unzip)
        self.cmd.run_bash(cmd)
        self.cmd.run(cmd_zip)

        if os.path.isfile(temp_fq_gz) and os.path.getsize(temp_fq_gz) > 100:
            os.remove(self.current)
            os.rename(temp_fq_gz, self.current)

    def remove_duplicate_reads(self):
        """
        Remove duplicate reads from current fastq file.
        """

        if not self.exists:
            return

        temp = os.path.join(
            os.path.dirname(self.current), f"temp_rmdup_{randint(1,1000)}.fq.gz"
        )

        cmd = "seqkit rmdup -n %s -o %s" % (self.current, temp)

        self.cmd.run(cmd)

        if os.path.isfile(temp) and os.path.getsize(temp):
            os.remove(self.current)
            os.rename(temp, self.current)

    def current_fastq_read_number(self):
        """
        Get cached number of reads in current fastq file."""

        if self.current_status == "raw":
            return self.read_number_raw
        elif self.current_status == "clean":
            return self.read_number_clean
        elif self.current_status == "enriched":
            return self.read_number_enriched
        elif self.current_status == "depleted":
            return self.read_number_depleted
        else:
            return 0

    def export_reads(self, directory):
        """
        Export reads to directory.
        """
        if not self.exists:
            return

        if not os.path.isdir(directory):
            os.makedirs(directory)

        if os.path.exists(os.path.join(directory, os.path.basename(self.current))):
            os.remove(os.path.join(directory, os.path.basename(self.current)))

        shutil.move(self.current, directory)
        self.current = os.path.join(directory, os.path.basename(self.current))

    def __str__(self):
        return self.filepath


class Sample_runClass:

    r1: Read_class
    r2: Read_class
    report: str
    qcdata: dict
    reads_before_processing: int = 0
    reads_after_processing: int = 0
    sample_dir: str = ""
    qc_soft: str = ""
    input_fastqc_report: os.PathLike
    processed_fastqc_report: os.PathLike
    threads: int = 1
    user_name: str
    reads_subdirectory: str = "reads"

    def __init__(
        self,
        r1: Read_class,
        r2: Read_class,
        sample_name: str,
        project_name: str,
        user_name: str,
        technology: str,
        type: str,
        combinations: str,
        input: str,
        bin: str,
        threads: int = 1,
    ) -> None:
        self.cmd = RunCMD(bin)
        self.r1 = r1
        self.r2 = r2
        self.sample_name = sample_name
        self.project_name = project_name
        self.user_name = user_name
        self.technology = technology
        self.type = type
        self.combinations = combinations
        self.input = input
        self.threads = threads

        self.QCdir = os.path.dirname(self.r1.clean)

        self.reads_before_processing = self.r1.read_number_raw + self.r2.read_number_raw

    def current_total_read_number(self):
        print("current total")
        print(self.r1.get_current_fastq_read_number())
        print(self.r2.get_current_fastq_read_number())
        return (
            self.r1.get_current_fastq_read_number()
            + self.r2.get_current_fastq_read_number()
        )

    def QC_summary(self):
        """get QC summary
        runs fastqc parse on zip files in QCdir

        :param QCdir:
        :return: qc_summary
        """
        qc_summary = {
            "input": fastqc_parse(os.path.join(self.QCdir, "input_data.zip")),
            "processed": fastqc_parse(os.path.join(self.QCdir, "processed_data.zip")),
        }

        return qc_summary

    def get_qc_data(self, reads_dir: str = "reads/clean/"):
        """get qc data for sample_class. Update sample_class.qc_data.

        :param sample_class:
        :return: None
        """

        self.qcdata = self.QC_summary()

    def get_fake_qc_data(self):
        """get fake qc data for sample_class. Update sample_class.qc_data.

        :param sample_class:
        :return: None
        """

        fake_df = pd.DataFrame(
            [
                ["Total_sequences", 1000000],
                ["Sequence length", 5000],
                ["%GC", 50],
                ["Q20", 90],
            ],
            columns=["measure", "value"],
        ).set_index("measure")

        self.qcdata = {
            "input": fake_df,
            "processed": fake_df,
        }

    def remove_duplicates(self):
        self.r1.remove_duplicate_reads()
        self.r2.remove_duplicate_reads()

    def clean_unique(self):
        if self.type == "SE":
            self.clean_unique_SE()
        else:
            self.clean_unique_PE()

    def clean_unique_SE(self):
        WHERETO = os.path.dirname(self.r1.current)
        unique_reads = os.path.join(WHERETO, "unique_reads.lst")

        cmd = [
            "zgrep",
            "'^@'",
            self.r1.current,
            "|",
            "cut",
            "-f1",
            "-d",
            " ",
            "|",
            'grep -v "+|/*|/&"',
            "|",
            "sed",
            r"'s/\s.*$//'",
            "|",
            "sed",
            "'s/^@//g'",
            "|",
            "sort",
            "|",
            "uniq",
            ">",
            unique_reads,
        ]

        self.cmd.run_bash(cmd)
        if not os.path.exists(unique_reads) or os.path.getsize(unique_reads) == 0:
            print(
                f"No unique reads found in {self.r1.current}, skipping unique read cleaning"
            )
            return

        self.r1.read_filter_inplace(self.r1.current, unique_reads)

    def clean_unique_PE(self):

        WHERETO = os.path.dirname(self.r1.current)
        common_reads = os.path.join(WHERETO, "common_reads.lst")

        cmd_find_common = [
            f"seqkit common -n -i {self.r1.current} {self.r1.current} | paste - - - - | cut -f1 | sort | uniq | sed 's/^@//g' > {common_reads}"
        ]

        self.cmd.run_script(cmd_find_common)

        if os.path.getsize(common_reads) == 0:
            return

        self.r1.read_filter_inplace(self.r1.current, common_reads)
        self.r2.read_filter_inplace(self.r2.current, common_reads)

    def trimmomatic_sort(self):
        if self.type == "SE":
            self.trimmomatic_sort_SE()
        else:
            self.trimmomatic_sort_PE()

    def trimmomatic_sort_SE(self):

        tempdir = os.path.dirname(self.r1.current)
        tempfq = os.path.join(tempdir, f"temp{randint(1,1999)}")

        cmd_trimsort = [
            "trimmomatic",
            "SE",
            "-threads",
            f"{self.threads}",
            self.r1.current,
            tempfq,
            "MINLEN:20",
        ]

        self.cmd.run(cmd_trimsort)

        if tempfq in os.listdir(tempdir):
            if os.path.getsize(tempfq) > 100:
                os.remove(self.r1.current)
                os.rename(tempfq, self.r1.current)

    def trimmomatic_sort_PE(self):
        if self.type == "SE":
            return

        tempdir = os.path.dirname(self.r1.current)
        tempfq = os.path.join(tempdir, f"temp{randint(1,1999)}")

        cmd_trimsort = [
            f"trimmomatic PE -phred33 -threads {self.threads} {self.r1.current} {self.r2.current} -baseout {tempfq}.fastq.gz MINLEN:20"
        ]

        self.cmd.run(cmd_trimsort)

        if tempfq + "_1P.fastq.gz" in os.listdir(tempdir):
            os.remove(self.r1.current)
            os.remove(self.r2.current)
            os.rename(tempfq + "_1P.fastq.gz", self.r1.current)
            os.rename(tempfq + "_2P.fastq.gz", self.r2.current)

    def export_reads(self, reads_dir: str = "reads/clean/"):
        """export reads to reads_dir

        :param reads_dir:
        :return: None
        """
        if self.type == "SE":
            self.r1.export_reads(reads_dir)
        else:
            self.r1.export_reads(reads_dir)
            self.r2.export_reads(reads_dir)


class Software_detail:
    def __init__(self, module, args_df: pd.DataFrame, config: dict, prefix: str):
        """

        Args:
            module: name of module.
            args_df: dataframe containing module arguments.
            config: dictionary containing module configuration.
            prefix: prefix of module.
        """
        if module not in args_df.module.unique():
            self.module = "None"
            self.name = "None"
            self.args = "None"
            self.db = "None"
            self.db_name = "None"
            self.bin = "None"
            self.dir = "None"
            self.output_dir = "None"
        else:
            method_details = args_df[(args_df.module == module)]
            self.module = module
            self.name = method_details.software.values[0]

            try:
                args_string = method_details[
                    method_details.parameter.str.contains("ARGS")
                ].value.values[0]

                self.args, self.db = self.excise_db_from_args(args_string)
                self.db_name = self.db.split("/")[-1]

            except IndexError:
                self.args = ""
                self.db = ""
                self.db_name = ""

            try:
                db_name = method_details[
                    method_details.parameter.str.contains("DB")
                ].value.values[0]
                if ".gz" in db_name:
                    self.db = os.path.join(config["source"]["REF_FASTA"], db_name)
                    self.db_name = db_name
                else:
                    self.db = os.path.join(config["source"]["DBDIR_MAIN"], db_name)
                    self.db_name = db_name

            except IndexError:
                pass

            try:
                self.bin = os.path.join(
                    config["bin"]["ROOT"], config["bin"]["software"][self.name], "bin"
                )
            except KeyError:
                try:
                    self.bin = os.path.join(
                        config["bin"]["ROOT"], config["bin"][module]["default"], "bin"
                    )
                except KeyError:
                    self.bin = ""

            try:
                self.dir = config["directories"][module]
            except KeyError:
                self.dir = ""

            print(
                f"Module: {self.module}, Software: {self.name}, Args: {self.args}, DB: {self.db}, Bin: {self.bin}, Dir: {self.dir}"
            )

            self.output_dir = os.path.join(self.dir, f"{self.name}.{prefix}")

    def excise_db_from_args(self, args: str):
        """replace db name in args with db path

        :param args: string of args
        :return: string of args with db path
        """
        args_list = args.split(" ")
        db_idx = None
        db = ""

        for i in range(len(args_list)):
            if "--db" in args_list[i]:
                db_idx = i
                db = args_list[i + 1]
                break

        if db_idx is not None:
            args_list = args_list[:db_idx] + args_list[db_idx + 2 :]
            args_list = " ".join(args_list)

            return args_list, db

        else:
            args_list = " ".join(args_list)
            return args_list, db

    def __str__(self):
        return f"{self.module}:{self.name}:{self.args}:{self.db}:{self.bin}"


class Bedgraph:
    """Class to store and work with bedgraph files


    Methods:
    read_bedgraph: reads bedgraph wiith coverage column.
    get_coverage_array: returns numpy array of coverage column
    get_bins: generates unzipped bins from start & end positions.
    plot_coverage: barplot of coverage by window in bdgraph.
    """

    def __init__(self, bedgraph_file, max_bars=7000, nbins=300):
        self.max_bars = max_bars
        self.nbins = nbins
        self.bedgraph = self.read_bedgraph(bedgraph_file)
        self.reduce_number_bars()
        self.bar_to_histogram()

    def read_bedgraph(self, coverage_file) -> pd.DataFrame:
        coverage = pd.read_csv(coverage_file, sep="\t", header=None).rename(
            columns={0: "read_id", 1: "start", 2: "end", 3: "coverage"}
        )

        return coverage

    def get_bins(self) -> np.ndarray:
        """
        Get the bins.

        :param coverage_file: The coverage file.
        """
        self.bedgraph["width"] = self.bedgraph.end - self.bedgraph.start

    def get_coverage_array(self, coverage: pd.DataFrame) -> np.ndarray:
        """
        Get the coverage of the remapping.

        :param coverage_file: The coverage file.
        """
        coverage_values = np.array(coverage.coverage.to_list())

        return coverage_values

    def get_bar_coordinates(self):
        """
        Get the bar coordinates.
        """
        self.bedgraph["width"] = self.bedgraph.end - self.bedgraph.start

        self.bedgraph["x"] = self.bedgraph.start + self.bedgraph.end / 2
        self.bedgraph["y"] = self.bedgraph.coverage

        return self.bedgraph

    def reduce_number_bars(self):
        """
        Reduce the number of bars.
        """
        self.bedgraph = self.bedgraph[self.bedgraph.coverage > 0]

        if self.bedgraph.shape[0] > self.max_bars:

            self.bedgraph = self.bedgraph.sample(self.max_bars)

    def merge_bedgraph_rows(self):
        """
        Merge the rows of the bedgraph.
        """

        for ix in range(1, self.bedgraph.shape[0]):
            if self.bedgraph.iloc[ix - 1].end < (self.bedgraph.iloc[ix].start - 1):
                self.bedgraph.iloc[ix].end = self.bedgraph.iloc[ix].start - 1

    def bar_to_histogram(self):
        """
        Bar to histogram.
        """

        self.bedgraph["coord"] = (self.bedgraph.start + self.bedgraph.end) / 2
        self.coverage = [
            [self.bedgraph.iloc[x]["coord"]] * self.bedgraph.iloc[x]["coverage"]
            for x in range(self.bedgraph.shape[0])
        ]
        self.coverage = list(it.chain.from_iterable(self.coverage))

    def plot_coverage(self, output_file, borders=50, tlen=0):
        """
        Plot the coverage of the remapping.

        :param coverage_file: The coverage file. bedgraph produced with samtools.
        :param output_file: The output file.
        """

        fig, ax = plt.subplots(figsize=(11, 3))

        if len(self.coverage) <= 1:
            return

        start_time = time.perf_counter()

        ax.hist(
            self.coverage,
            bins=self.nbins,
            color="skyblue",
            edgecolor="none",
        )

        ax.set_xlabel("Reference")
        ax.set_ylabel("Coverage")

        ##
        xmax = self.bedgraph.end.max()
        if tlen:
            xmax = tlen
        ax.set_xlim(0 - borders, xmax + borders)
        ##

        fig.savefig(output_file, bbox_inches="tight")
        ax.cla()
        fig.clf()
        plt.close("all")


@dataclass
class Remap_Target:
    accid: str
    acc_simple: str
    taxid: str
    file: str
    run_prefix: str
    description: str = ""
    accid_in_file: list = field(default_factory=lambda: [])
    reads: bool = False
    contigs: bool = False

    def __post_init__(self):
        self.name = f"{self.run_prefix}_{self.acc_simple}_{self.taxid}_{os.path.splitext(os.path.basename(self.file))[0]}"


@dataclass(frozen=True)
class Run_detail_report:
    max_depth: int
    max_depthR: int
    max_gaps: int
    max_prop: int
    max_mapped: int
    input: str
    enriched_reads: int
    enriched_reads_percent: float
    depleted_reads: int
    depleted_reads_percent: float
    processed: str
    processed_percent: str
    sift_preproc: bool
    sift_remap: bool
    sift_removed_pprc: str
    processing_final: str
    processing_final_percent: str
    merged: bool
    merged_number: int
    merged_files: str


@dataclass(frozen=True)
class Contig_classification_results:
    performed: bool
    method: str
    args: str
    db: str
    classification_number: int
    classification_minhit: int
    success: bool


@dataclass(frozen=True)
class Read_classification_results:
    performed: bool
    method: str
    args: str
    db: str
    classification_number: int
    classification_minhit: int
    success: bool


@dataclass(frozen=True)
class Remap_main:
    performed: bool
    success: int
    method: str
    found_total: int
    coverage_min: int
    coverage_max: int


@dataclass(frozen=True)
class Assembly_results:
    performed: bool
    assembly_soft: str
    assembly_args: str
    assembly_number: int
    assembly_min: int
    assembly_mean: int
    assembly_max: int
    assembly_trim: int
