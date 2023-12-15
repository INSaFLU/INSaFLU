import itertools as it
import logging
import os
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from random import randint
from typing import List, Type

import matplotlib
import pandas as pd
from numpy import ERR_CALL

from pathogen_identification.constants_settings import ConstantsSettings
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

    def __init__(self, temp_dir: str, prefix: str = "temp", suffix: str = ".txt"):
        """
        Initialize.
        """

        self.temp_dir = temp_dir
        self.prefix = prefix
        self.suffix = suffix

        self.path = os.path.join(
            self.temp_dir, f"{self.prefix}_{randint(1000000, 9999999)}{self.suffix}"
        )
        self.file = os.path.basename(self.path)

    def __enter__(self):
        """
        Enter.
        """

        open(self.path, "w").close()
        return self.path

    def __exit__(
        self,
        exc_type: Type[BaseException],
        exc_value: BaseException,
        traceback,
    ):
        """
        Exit.
        """

        if os.path.exists(self.path):
            os.remove(self.path)


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
        print("Operation temp dir: " + self.temp_dir)

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
        time_delay = 0

        print(f"timeout: {ConstantsSettings.TIMEOUT} seconds")

        while not found_flag:
            time.sleep(1)
            found_flag = os.path.exists(self.flag)
            time_delay += 1

            if time_delay > ConstantsSettings.TIMEOUT:
                proc_prep.kill()
                err = "Timeout"
                out = "Timeout"
                exec_time = time.perf_counter() - start_time

                return out, err, exec_time

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
        self.prefix = prefix

        # set logger

        self.logger = logging.getLogger(f"{prefix}_{task}")
        self.logger.setLevel(logging.ERROR)

        # remove handlers
        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)

        self.logger.addHandler(
            logging.FileHandler(os.path.join(logdir, f"{prefix}_{task}.log"))
        )
        # set handler for stdout
        self.logger.propagate = False

        self.logfile = os.path.join(logdir, f"{prefix}_{task}.log")
        self.logdir = logdir

    def set_logger(self, logger):
        """
        Set logger.
        """

        self.logger = logger

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

        if isinstance(out, bytes):
            out = out.decode("utf-8", errors="ignore")

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
                f.write(f"####################\n")

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
        print(f"running: {self.bin}{cmd}")

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
        self,
        filepath,
        clean_dir: str,
        enriched_dir: str,
        depleted_dir: str,
        bin: str,
        prefix: str = "r0",
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
        print("Initializing Read_class")
        print("clean_dir", clean_dir)
        self.cmd = RunCMD(bin, prefix="read", task="housekeeping", logdir=clean_dir)

        self.exists = os.path.isfile(filepath)

        self.filepath = filepath
        self.current = filepath
        self.suffix = prefix
        self.prefix = self.determine_file_name(filepath)
        self.clean = os.path.join(clean_dir, self.prefix + ".clean.fastq.gz")
        self.enriched = os.path.join(enriched_dir, self.prefix + ".enriched.fastq.gz")
        self.depleted = os.path.join(depleted_dir, self.prefix + ".depleted.fastq.gz")
        self.current_status = "raw"
        self.read_number_raw = self.get_current_fastq_read_number()
        self.read_number_clean = 0
        self.read_number_enriched = 0
        self.enriched_read_number = 0
        self.read_number_depleted = 0
        self.depleted_read_number = 0
        self.read_number_filtered = 0
        self.history = [self.current]

    def create_link(self, file_path, new_path):
        if os.path.isfile(file_path):
            if os.path.isfile(new_path):
                os.remove(new_path)
            os.symlink(file_path, new_path)

    def update(self, new_suffix, clean_dir: str, enriched_dir: str, depleted_dir: str):
        self.prefix = self.prefix.replace(self.suffix, new_suffix)
        self.suffix = new_suffix
        new_clean = os.path.join(clean_dir, self.prefix + ".clean.fastq.gz")
        if os.path.isfile(self.clean):
            if new_clean != self.clean:
                self.create_link(self.clean, new_clean)
        self.clean = new_clean

        new_enriched = os.path.join(enriched_dir, self.prefix + ".enriched.fastq.gz")
        if os.path.isfile(self.enriched):
            if new_enriched != self.enriched:
                self.create_link(self.enriched, new_enriched)
        self.enriched = new_enriched

        new_depleted = os.path.join(depleted_dir, self.prefix + ".depleted.fastq.gz")
        if os.path.isfile(self.depleted):
            if new_depleted != self.depleted:
                self.create_link(self.depleted, new_depleted)
        self.depleted = new_depleted

        if self.current_status == "raw":
            self.current = self.filepath
        elif self.current_status == "clean":
            self.current = self.clean
        elif self.current_status == "enriched":
            self.current = self.enriched
        elif self.current_status == "depleted":
            self.current = self.depleted

    def get_read_names_fastq(self):
        """
        Get read names from fastq file.
        """
        filepath = self.current

        read_names = []
        counter = 0

        if self.exists:
            with gzip.open(filepath, "rt") as f:
                for line in f:
                    if counter == 0:
                        read_names.append(line.split()[0][1:])

                    counter += 1

                    if counter == 4:
                        counter = 0

        return read_names

    def determine_file_name(self, filepath):
        if not self.exists:
            return "none"

        if "gz" not in filepath:
            self.cmd.run("bgzip -f %s" % filepath)

        filename = os.path.basename(filepath)

        filename = filename.replace(".gz", "")

        filename = filename.replace(".fasta", "").replace(".fa", "")

        return filename

    def read_filter_move(self, read_list: list, output: str = ""):
        if not self.exists:
            return

        temp_file = Temp_File(os.path.dirname(output), suffix=".lst")

        with temp_file as tpf:
            with open(tpf, "w") as f:
                f.write("\n".join(read_list))

            self.read_filter(output, tpf)

    def read_filter_inplace(self, read_list: str):
        if not self.exists:
            return

        output_dir = os.path.dirname(self.current)

        tempreads = os.path.join(output_dir, f"temp_{randint(1,1999)}.fq.gz")

        self.read_filter(tempreads, read_list)

        if os.path.isfile(tempreads) and os.path.getsize(tempreads) > 100:
            os.remove(self.current)
            os.rename(tempreads, self.current)

    def read_filter(self, output: str, read_list: str):
        """
        filter read file using exisiting lsit of reads.
        Args:

            input: path to input file.
            output: path to output file.
            read_list: path to file containing read list.
        """
        if not self.exists:
            return

        cmd = "seqtk subseq %s %s | gzip > %s" % (self.current, read_list, output)

        self.cmd.run(cmd)

    def update_history(self):
        """
        Update history of read file.
        """
        self.history.append(self.current)

    def enrich(self, read_list):
        """
        filter reads and set current status to enriched.
        """

        reads_start = self.get_current_fastq_read_number()

        if len(read_list) > 0:
            self.read_filter_move(read_list, self.enriched)
            self.is_enriched()

        reads_end = self.get_current_fastq_read_number()

        self.enriched_read_number = reads_start - reads_end

    def deplete(self, read_list):
        """
        filter reads and aset current status to depleted.
        """
        reads_start = self.get_current_fastq_read_number()
        current_reads = self.get_read_names_fastq()

        read_list_to_keep = list(set(current_reads) - set(read_list))

        if len(read_list) > 0:
            self.read_filter_move(read_list_to_keep, self.depleted)
            self.remove_duplicate_reads()
            self.is_depleted()

        reads_end = self.get_current_fastq_read_number()

        self.depleted_read_number = reads_start - reads_end

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
        rnumber = self.cmd.run_bash_return(cmd)
        return int(rnumber) // 4

    def clean_fastq_headers_python(self, fastq, temp_fq):
        """
        Clean fastq header using python.
        """

        if not self.exists:
            return

        counter = 0  # counter to keep track of headers
        with open(temp_fq, "w") as f:
            for line in open(fastq):
                line = line.strip()
                if line.startswith("@") and counter == 0:
                    if line[-2:] == "/1" or line[-2:] == "/2":
                        line = line[:-2]

                f.write(line + "\n")

                counter += 1
                if counter == 4:
                    counter = 0

    def clean_read_names(self):
        """
        Clean read names in current fastq file.
        """

        if not self.exists:
            return

        temp_fq = os.path.join(
            os.path.dirname(self.current), f"temp_clean_{randint(1,1000)}.fastq"
        )
        final_temp = os.path.join(
            os.path.dirname(self.current), f"temp_second_{randint(1,1000)}.fastq"
        )

        temp_fq_gz = final_temp + ".gz"

        cmd_unzip = "gunzip -c %s > %s" % (self.current, temp_fq)
        cmd_zip = "bgzip %s" % final_temp

        self.cmd.run_bash(cmd_unzip)
        self.clean_fastq_headers_python(temp_fq, final_temp)
        self.cmd.run(cmd_zip)
        os.remove(temp_fq)

        if os.path.isfile(temp_fq_gz) and os.path.getsize(temp_fq_gz) > 100:
            os.remove(self.current)
            os.rename(temp_fq_gz, self.current)

    def clean_read_names_old(self):
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

        final_file = os.path.join(directory, os.path.basename(self.current))

        if os.path.exists(self.current):
            if os.path.exists(final_file):
                os.remove(final_file)

            shutil.move(self.current, directory)

        self.current = final_file

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
        self.cmd = RunCMD(
            bin, prefix="sample", task="housekeeping", logdir=os.path.dirname(r1.clean)
        )
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

    def sources_list(self):
        return tuple(
            [
                self.r1.history[-1],
                self.r2.history[-1],
            ]
        )

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
        if self.type == ConstantsSettings.SINGLE_END:
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

        self.r1.read_filter_inplace(unique_reads)

    def clean_unique_PE(self):
        WHERETO = os.path.dirname(self.r1.current)
        common_reads = os.path.join(WHERETO, "common_reads.lst")

        cmd_find_common = [
            f"seqkit common -n -i {self.r1.current} {self.r1.current} | paste - - - - | cut -f1 | sort | uniq | sed 's/^@//g' > {common_reads}"
        ]

        self.cmd.run_script(cmd_find_common)

        if os.path.getsize(common_reads) == 0:
            return

        self.r1.read_filter_inplace(common_reads)
        self.r2.read_filter_inplace(common_reads)

    def trimmomatic_sort(self):
        if self.type == ConstantsSettings.SINGLE_END:
            self.trimmomatic_sort_SE()
        else:
            self.trimmomatic_sort_PE()

    def trimmomatic_sort_SE(self):
        tempdir = os.path.dirname(self.r1.current)
        tempfq = os.path.join(tempdir, f"temp{randint(1,1999)}")

        cmd_trimsort = [
            "trimmomatic",
            ConstantsSettings.SINGLE_END,
            "-phred33",
            "-threads",
            f"{self.threads}",
            self.r1.current,
            tempfq,
            "MINLEN:20",
        ]

        self.cmd.run_script_software(cmd_trimsort)

        if tempfq in os.listdir(tempdir):
            if os.path.getsize(tempfq) > 100:
                bgzip_cmd = [
                    "bgzip",
                    tempfq,
                ]
                self.cmd.run_script_software(bgzip_cmd)
                os.remove(self.r1.current)
                os.rename(tempfq, self.r1.current)

    def trimmomatic_sort_PE(self):
        if self.type == ConstantsSettings.SINGLE_END:
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
        if self.type == ConstantsSettings.SINGLE_END:
            self.r1.export_reads(reads_dir)
        else:
            self.r1.export_reads(reads_dir)
            self.r2.export_reads(reads_dir)


class SoftwareUnit:
    name: str

    SOFTWARE_NOT_FOUND = "None"

    def __init__(
        self,
        module: str = SOFTWARE_NOT_FOUND,
        name: str = SOFTWARE_NOT_FOUND,
        args: str = SOFTWARE_NOT_FOUND,
        db: str = SOFTWARE_NOT_FOUND,
        db_name: str = SOFTWARE_NOT_FOUND,
        bin: str = SOFTWARE_NOT_FOUND,
        dir: str = SOFTWARE_NOT_FOUND,
        output_dir: str = SOFTWARE_NOT_FOUND,
    ):
        """ """
        self.module = module
        self.name = name
        self.args = args
        self.db = db
        self.db_name = db_name
        self.bin = bin
        self.dir = dir
        self.output_dir = output_dir

    def check_exists(self):
        return bool(self.name != Software_detail.SOFTWARE_NOT_FOUND)

    def get_bin(self, config: dict):
        try:
            self.bin = os.path.join(
                config["bin"]["ROOT"], config["bin"]["software"][self.name], "bin"
            )
        except KeyError:
            try:
                self.bin = os.path.join(
                    config["bin"]["ROOT"], config["bin"][self.module]["default"], "bin"
                )
            except KeyError:
                self.bin = ""


class Software_detail(SoftwareUnit):
    def __init__(self, module, args_df: pd.DataFrame, config: dict, prefix: str):
        """

        Args:
            module: name of module.
            args_df: dataframe containing module arguments.
            config: dictionary containing module configuration.
            prefix: prefix of module.
        """
        super().__init__()

        if module in args_df.module.unique():
            method_details = args_df[(args_df.module == module)]
            self.module = module
            self.name = method_details.software.values[0]

            self.extract_args(method_details)

            self.extract_db(method_details, config)

            self.get_bin(config)

            self.get_dir_from_config(config)

            print(
                f"Module: {self.module}, Software: {self.name}, Args: {self.args}, DB: {self.db}, Bin: {self.bin}, Dir: {self.dir}"
            )

            self.output_dir = os.path.join(self.dir, f"{self.name}.{prefix}")

    def get_dir_from_config(self, config: dict):
        try:
            self.dir = config["directories"][self.module]
        except KeyError:
            self.dir = ""

    def extract_db(self, method_details: pd.DataFrame, config: dict):
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

    def set_db(self, filepath, name=""):
        if name == "":
            name = os.path.basename(filepath)

        self.db = filepath
        self.db_name = name

    def extract_args(self, method_details: pd.DataFrame):
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


class SoftwareDetailCompound:
    def __init__(self, module, args_df: pd.DataFrame, config: dict, prefix: str):
        """

        Args:
            module: name of module.
            args_df: dataframe containing module arguments.
            config: dictionary containing module configuration.
            prefix: prefix of module.
        """
        self.module = module
        self.args_df = args_df
        self.config = config
        self.prefix = prefix

        self.software_list: List[Software_detail] = []
        if module in args_df.module.unique():
            self.fill_software_list()

    def fill_software_list(self):
        print(self.args_df)
        module_df = self.args_df[self.args_df.module == self.module]

        for _, software_df in module_df.groupby("software"):
            software = Software_detail(
                self.module, software_df, self.config, self.prefix
            )
            self.software_list.append(software)

    def check_exists(self):
        return any([x.check_exists() for x in self.software_list])


class SoftwareRemap:
    def __init__(
        self, remap_software: Software_detail, remap_filter: SoftwareDetailCompound
    ):
        self.remap_software = remap_software
        self.remap_filters = remap_filter

    @property
    def output_dir(self):
        return self.remap_software.dir

    def set_output_dir(self, output_dir):
        self.remap_software.dir = output_dir


@dataclass
class MappingStats:
    def __init__(
        self,
        error_rate: float,
        quality_avg: float,
    ):
        self.error_rate = error_rate
        self.quality_avg = quality_avg


class Bedgraph:
    """Class to store and work with bedgraph files


    Methods:
    read_bedgraph: reads bedgraph wiith coverage column.
    get_coverage_array: returns numpy array of coverage column
    get_bins: generates unzipped bins from start & end positions.
    plot_coverage: barplot of coverage by window in bdgraph.
    """

    def __init__(self, bedgraph_file, max_bars=1000, nbins=500):
        self.max_bars = max_bars
        self.nbins = nbins
        self.bedgraph = self.read_bedgraph(bedgraph_file)
        if self.bedgraph.end.max() <= nbins:
            self.nbins = int(self.bedgraph.end.max())

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

    def get_coverage_array(self, new_bed_coordinates: pd.DataFrame) -> np.ndarray:
        """
        Get the coverage of the remapping.

        :param coverage_file: The coverage file.
        """

        def new_coordinates(x):
            """ """
            average_bedgraph = self.bedgraph[
                (self.bedgraph.end >= x.start) & (self.bedgraph.start <= x.end)
            ]

            average_bedgraph.iloc[0]["start"] = x.start
            average_bedgraph.iloc[average_bedgraph.shape[0] - 1]["end"] = x.end

            if "width" not in average_bedgraph.columns:
                average_bedgraph["width"] = (
                    average_bedgraph.end - average_bedgraph.start
                )

            average_bedgraph["coverage"] = (
                average_bedgraph.coverage * average_bedgraph.width
            )

            coverage = average_bedgraph.coverage.sum() / (x.end - x.start)

            return coverage

        new_bed_coordinates["coverage"] = new_bed_coordinates.apply(
            new_coordinates, axis=1
        )
        new_bed_coordinates.fillna(0, inplace=True)

        return new_bed_coordinates

    @staticmethod
    def get_bar_coordinates(new_bed_coordinates: pd.DataFrame):
        """
        Get the bar coordinates.
        """
        new_bed_coordinates["width"] = (
            new_bed_coordinates.end - new_bed_coordinates.start
        )
        new_bed_coordinates["coord"] = (
            new_bed_coordinates.start + new_bed_coordinates.end
        ) / 2

        return new_bed_coordinates

    def standardize_bedgraph_mean(self):
        """ """
        chromosome = self.bedgraph.read_id.unique()[0]
        new_bed_range = [0, max(self.bedgraph.end)]
        new_bed_bins = np.linspace(new_bed_range[0], new_bed_range[1], self.nbins)
        new_bed_coordinates = [
            [chromosome, int(x), int(y)]
            for x, y in zip(new_bed_bins[:-1], new_bed_bins[1:])
        ]
        new_bed_coordinates = pd.DataFrame(
            new_bed_coordinates, columns=["read_id", "start", "end"]
        )
        new_bed_coordinates = self.get_coverage_array(new_bed_coordinates)
        new_bed_coordinates = self.get_bar_coordinates(new_bed_coordinates)

        return new_bed_coordinates

    def reduce_number_bars(self):
        """
        Reduce the number of bars.
        """
        new_bedgraph = self.bedgraph[self.bedgraph.coverage > 0]

        if new_bedgraph.shape[0] > self.max_bars:
            new_bedgraph = new_bedgraph.sample(self.max_bars)

        new_bedgraph = self.merge_bedgraph_rows(new_bedgraph)

    @staticmethod
    def merge_bedgraph_rows(bedgraph):
        """
        Merge the rows of the bedgraph.
        """

        for ix in range(1, bedgraph.shape[0]):
            if bedgraph.iloc[ix - 1].end < (bedgraph.iloc[ix].start - 1):
                bedgraph.iloc[ix].end = bedgraph.iloc[ix].start - 1

        return bedgraph

    def plot_coverage_bar(self, output_file, borders=50, tlen=0):
        """
        Plot the coverage of the remapping.

        :param coverage_file: The coverage file. bedgraph produced with samtools.
        :param output_file: The output file.
        """

        new_bedgraph = self.standardize_bedgraph_mean()

        fig, ax = plt.subplots(figsize=(11, 3))

        if len(new_bedgraph.shape) <= 1:
            return

        start_time = time.perf_counter()

        ax.bar(
            new_bedgraph.coord,
            new_bedgraph.coverage,
            width=new_bedgraph.width,
            color="skyblue",
            edgecolor="none",
        )

        ax.set_xlabel("Reference", fontsize=9)
        ax.set_ylabel(f"Coverage ({self.nbins} windows)", fontsize=9)

        ##
        xmax = new_bedgraph.end.max()
        if tlen:
            xmax = tlen
        ax.set_xlim(0 - borders, xmax + borders)
        ##

        plt.text(
            0.05,
            0.9,
            f"Bar width: {round((new_bedgraph.end.max() - new_bedgraph.start.min()) / self.nbins)} bp",
            horizontalalignment="left",
            verticalalignment="top",
            transform=ax.transAxes,
            in_layout=True,
            backgroundcolor="white",
            fontsize=9,
        )
        fig.savefig(output_file, bbox_inches="tight")

        ax.cla()
        fig.clf()
        plt.close("all")

    def plot_coverage(self, output_file, borders=50, tlen=0):
        self.plot_coverage_bar(output_file, borders, tlen)


@dataclass
class Remap_Target:
    accid: str
    acc_simple: str
    taxid: str
    file: str
    run_prefix: str
    description: str = ""
    accid_in_file: List[str] = field(default_factory=lambda: [])
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
class RunQC_report:
    performed: bool
    method: str
    args: str
    input_reads: int
    output_reads: int
    output_reads_percent: float


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
