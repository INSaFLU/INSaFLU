#!/bin/python

import logging
import os
import subprocess
import sys
from typing import Type

from pathogen_detection.object_classes import RunCMD


class Preprocess:
    def __init__(
        self,
        r1,
        r2,
        preprocess_dir,
        preprocess_type,
        preprocess_method,
        preprocess_name_fastq_gz,
        preprocess_name_r2_fastq_gz="",
        threads: int = 1,
        subsample: bool = False,
        subsample_reads: int = 0,
        logging_level: int = logging.ERROR,
    ):
        """

        Args:
            r1: path to r1 fastq
            r2: path to r2 fastq
            preprocess_dir: path to preprocess directory
            preprocess_type: PE or SE
            preprocess_method: preprocess method
            preprocess_name_fastq_gz: name of preprocessed fastq
            preprocess_name_r2_fastq_gz: name of preprocessed r2 fastq
            threads: number of threads
            subsample: subsample reads
            logging_level: logging level

        :param r1: str
        :param r2: str
        :param preprocess_dir: str
        :param preprocess_type: str
        :param preprocess_method: PreprocessMethod
        :param preprocess_name_fastq_gz: str
        :param preprocess_name_r2_fastq_gz: str
        :param threads: int
        :param subsample: bool
        :param logging_level: int

        """
        self.r1 = r1
        self.r2 = r2
        self.threads = str(threads)
        self.preprocess_dir = preprocess_dir
        self.bin = preprocess_method.bin
        self.args = preprocess_method.args
        self.subsample = subsample
        self.subsample_number = subsample_reads

        self.preprocess_type = preprocess_type
        self.preprocess_method = preprocess_method

        self.preprocess_name_fastq_gz = preprocess_name_fastq_gz
        self.preprocess_name_fastq = self.preprocess_name_fastq_gz.replace(".gz", "")
        self.preprocess_name = self.preprocess_name_fastq.replace(".fastq", "")
        self.preprocess_name_r2_fastq_gz = preprocess_name_r2_fastq_gz
        self.preprocess_name_r2_fastq = self.preprocess_name_r2_fastq_gz.replace(
            ".gz", ""
        )
        self.preprocess_name_r2 = self.preprocess_name_r2_fastq.replace(".fastq", "")
        self.cmd = RunCMD(self.bin)

        self.input_qc_report = os.path.join(self.preprocess_dir, "input_data.html")
        self.processed_qc_report = os.path.join(
            self.preprocess_dir, "processed_data.html"
        )

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging_level)
        self.logger.addHandler(logging.StreamHandler())
        self.logger.info("Preprocess class initialized")
        self.logger.info("Preprocess type: {}".format(self.preprocess_type))
        self.logger.info("Preprocess method: {}".format(self.preprocess_method.name))
        self.logger.info("args: {}".format(self.args))

    def check_gz_file_not_empty(self, file):
        """
        Check if gzipped file has lines starting with read or sequence indicators.
        """
        cmd = "gunzip -c {} | grep '^>\|^@' | wc -l".format(file)
        number_of_sequences = int(self.cmd.run_bash_return(cmd))
        print("Number of sequences: {}".format(number_of_sequences))
        if number_of_sequences > 0:
            return True
        else:
            return False

    def check_already_preprocessed(self):
        """
        Check if preprocessed files exist
        """
        if self.preprocess_type == "PE":
            return self.check_already_preprocessed_PE()
        elif self.preprocess_type == "SE":
            return self.check_already_preprocessed_SE()

    def check_already_preprocessed_PE(self):
        """
        Check if preprocessed files exist
        """
        if (
            self.r1_output_check()
            and self.r2_output_check()
            and self.check_fastqc_performed()
        ):
            return True
        else:
            return False

    def check_already_preprocessed_SE(self):
        """
        Check if preprocessed files exist
        """
        if self.r1_output_check() and self.check_fastqc_performed():
            return True
        else:
            return False

    def r1_output_check(self):
        """
        Check if filtered output exists
        """
        if os.path.exists(
            self.preprocess_name_fastq_gz
        ) and self.check_gz_file_not_empty(self.preprocess_name_fastq_gz):
            return True
        else:
            return False

    def r2_output_check(self):
        """
        Check if filtered output exists
        """
        if os.path.exists(
            self.preprocess_name_r2_fastq_gz
        ) and self.check_gz_file_not_empty(self.preprocess_name_r2_fastq_gz):
            return True
        else:
            return False

    def check_fastqc_performed(self):
        """
        Check if fastqc was performed
        """
        if os.path.exists(
            os.path.join(self.preprocess_dir, "input_data.html")
        ) and os.path.exists(os.path.join(self.preprocess_dir, "processed_data.html")):
            return True
        else:
            return False

    def run(self):

        if self.check_already_preprocessed():
            self.logger.info(
                "Preprocessed files already exist. Skipping preprocessing."
            )
            return

        if self.subsample:
            self.subsample_reads()

        self.fastqc_input()

        self.preprocess_QC()
        self.clean_read_names()

        self.fastqc_processed()

    def r2_output_check(self):
        """
        Check if filtered output exists
        """
        if os.path.exists(self.preprocess_name_r2_fastq_gz) and os.path.getsize(
            self.preprocess_name_r2_fastq_gz
        ):
            return True
        else:
            return False

    def preprocess_QC(self):
        """
        Configure preprocess
        """
        if self.preprocess_method.name == "trimmomatic":
            self.run_trimmomatic()
        elif self.preprocess_method.name == "nanofilt":
            self.run_nanofilt()
        else:
            raise ValueError(
                "preprocess method {} not supported".format(self.preprocess_method.name)
            )

    def fastqc_input(self, suffix="input_data"):
        """
        Fastqc
        """
        if self.preprocess_type == "PE":
            self.fastqc_PE()
        elif self.preprocess_type == "SE":
            self.fastqc_SE()
        else:
            raise ValueError("read type {} not supported".format(self.preprocess_type))

        subprocess.run(
            [
                "mv",
                os.path.join(self.preprocess_dir, "stdin_fastqc.html"),
                os.path.join(self.preprocess_dir, f"{suffix}.html"),
            ]
        )

        subprocess.run(
            [
                "mv",
                os.path.join(self.preprocess_dir, "stdin_fastqc.zip"),
                os.path.join(self.preprocess_dir, f"{suffix}.zip"),
            ]
        )

    def fastqc_PE(self):
        """
        Fastqc PE
        """
        fastq_cmd = [
            "zcat",
            self.r1,
            self.r2,
            "|",
            os.path.join(self.cmd.bin, "fastqc"),
            "stdin",
            "--outdir",
            self.preprocess_dir,
        ]

        self.cmd.run_bash(fastq_cmd)

    def fastqc_SE(self):
        """
        Fastqc SE
        """
        fastq_cmd = [
            "zcat",
            self.r1,
            "|",
            os.path.join(self.cmd.bin, "fastqc"),
            "stdin",
            "--outdir",
            self.preprocess_dir,
        ]

        self.cmd.run_bash(fastq_cmd)

    def fastqc_processed(self, suffix="processed_data"):
        """
        Fastqc
        """
        if self.preprocess_type == "PE":
            self.fastqc_processed_PE()
        elif self.preprocess_type == "SE":
            self.fastqc_processed_SE()
        else:
            raise ValueError("read type {} not supported".format(self.preprocess_type))

        subprocess.run(
            [
                "mv",
                os.path.join(self.preprocess_dir, "stdin_fastqc.html"),
                os.path.join(self.preprocess_dir, f"{suffix}.html"),
            ]
        )

        subprocess.run(
            [
                "mv",
                os.path.join(self.preprocess_dir, "stdin_fastqc.zip"),
                os.path.join(self.preprocess_dir, f"{suffix}.zip"),
            ]
        )

    def fastqc_processed_PE(self):
        """
        Fastqc PE
        """
        fastq_cmd = [
            "zcat",
            self.preprocess_name_fastq_gz,
            self.preprocess_name_r2_fastq_gz,
            "|",
            os.path.join(self.cmd.bin, "fastqc"),
            "stdin",
            "--outdir",
            self.preprocess_dir,
        ]

        self.cmd.run_bash(fastq_cmd)

    def fastqc_processed_SE(self):
        """
        Fastqc SE
        """
        fastq_cmd = [
            "zcat",
            self.preprocess_name_fastq_gz,
            "|",
            os.path.join(self.cmd.bin, "fastqc"),
            "stdin",
            "--outdir",
            self.preprocess_dir,
        ]

        self.cmd.run_bash(fastq_cmd)

    def run_trimmomatic(self):
        """
        Trimmomatic
        """
        if self.preprocess_type == "PE":
            self.trimmomatic_PE()
        elif self.preprocess_type == "SE":
            self.trimmomatic_SE()
        else:
            raise ValueError("read type {} not supported".format(self.preprocess_type))

    def trimmomatic_PE(self):
        """
        Trimmomatic PE
        """
        trimmomatic_cmd = [
            "trimmomatic",
            "PE",
            "-threads",
            self.threads,
            self.r1,
            self.r2,
            "-baseout",
            self.preprocess_name_fastq_gz,
            self.args,
        ]

        self.cmd.run(trimmomatic_cmd)

        subprocess.run(
            ["cp", self.preprocess_name + "_1P.fastq.gz", self.preprocess_name_fastq_gz]
        )
        subprocess.run(
            [
                "cp",
                self.preprocess_name + "_2P.fastq.gz",
                self.preprocess_name_r2_fastq_gz,
            ]
        )

    def trimmomatic_SE(self):
        """
        Trimmomatic SE
        """
        trimmomatic_cmd = [
            "trimmomatic",
            "SE",
            "-threads",
            self.threads,
            self.r1,
            self.preprocess_name_fastq,
            self.args,
        ]

        self.cmd.run(trimmomatic_cmd)

    def run_nanofilt(self):
        """
        Nanofilt
        """
        if self.preprocess_type == "PE":
            self.nanofilt_PE()
        elif self.preprocess_type == "SE":
            self.nanofilt_SE()
        else:
            raise ValueError("read type {} not supported".format(self.preprocess_type))

    def nanofilt_PE(self):
        """
        Nanofilt PE
        """
        nanofilt_cmd = [
            "NanoFilt",
            "-1",
            self.preprocess_name_fastq,
            "-2",
            self.preprocess_name_fastq,
            "-o",
            self.preprocess_name_fastq,
            "-p",
            "pe",
            self.args,
        ]

        self.cmd.run(nanofilt_cmd)

    def nanofilt_SE(self):
        """
        Nanofilt SE
        """
        temp_file_path = os.path.dirname(self.preprocess_name_fastq) + "/temp_r1.fastq"
        unzip_cmd = ["gunzip", "-c", self.r1, ">", temp_file_path]
        nanofilt_cmd = [
            "NanoFilt",
            self.args,
            "--readtype",
            "1D",
            temp_file_path,
            ">",
            self.preprocess_name_fastq,
        ]

        gzip_cmd = ["gzip", self.preprocess_name_fastq]

        self.cmd.run_bash(unzip_cmd)
        self.cmd.run(nanofilt_cmd)
        if os.path.isfile(self.preprocess_name_fastq_gz):
            os.remove(self.preprocess_name_fastq_gz)
        self.cmd.run_bash(gzip_cmd)

    def clean_read_names(self):
        """
        Clean read names
        """
        if self.preprocess_type == "PE":
            self.clean_read_names_single(self.preprocess_name_fastq)
            self.clean_read_names_single(self.preprocess_name_r2_fastq)

        elif self.preprocess_type == "SE":
            self.clean_read_names_single(self.preprocess_name_fastq)
        else:
            raise ValueError("read type {} not supported".format(self.preprocess_type))

    def clean_read_names_single(self, read_name):
        """
        Clean read names PE
        """

        self.cmd.run_bash(["gunzip", read_name + ".gz"])

        self.cmd.run_bash(
            [
                "sed",
                "-i",
                "'/^[@+]/ s/\/[12]$//g'",
                read_name,
            ]
        )
        self.cmd.run_bash(
            [
                "sed",
                "-i",
                "'/^[@+]/ s/\s.*$//g'",
                read_name,
            ]
        )

        self.cmd.run_bash(["sed", "-i", "'/^[@+]/ s/-/_/g'", read_name])

        self.cmd.run_bash(
            [
                "gzip",
                read_name,
            ]
        )

    def subsample_reads(self):
        if self.preprocess_type == "SE":
            self.subsample_SE()
        if self.preprocess_type == "PE":
            self.subsample_PE()
            self.clean_unpaired()

    def subsample_SE(self):
        """
        Subsample SE, use seqtk to subsample.
        """
        from random import randint

        seed = randint(1, 100000)
        tempdir = os.path.dirname(self.r1)
        tempname = f"temp_subsample_{seed}"

        temp1 = os.path.join(tempdir, tempname + "_r1.fq")

        seqtk_subsample_r1_cmd = [
            "seqtk",
            "sample",
            "-s",
            str(seed),
            self.r1,
            str(self.subsample_number),
            ">",
            temp1,
        ]

        compress_and_relocate_r1_cmd = ["bgzip", "-c", temp1, ">", self.r1]

        self.cmd.run(seqtk_subsample_r1_cmd)
        self.cmd.run(compress_and_relocate_r1_cmd)

        os.system(f"rm {temp1}")

    def subsample_PE(self):
        """
        subsample Paired reads. i) use seqkit to find common reads, ii) subsample that, iii) use seqtk to subsample.
        """
        if self.subsample_number == 0:
            return

        WHERETO = os.path.dirname(self.r1.current)
        common_reads = os.path.join(WHERETO, "common_reads.lst")

        cmd_find_common = [
            "seqkit",
            "common",
            "-n",
            "-i",
            self.r1.current,
            self.r1.current,
            "|",
            "paste",
            "- - - -",
            "|",
            "cut",
            "-f1",
            "|",
            "uniq",
            "|",
            "sed",
            "'s/^@//g'",
            "|",
            "shuf",
            "-n",
            str(self.subsample_number),
            ">",
            common_reads,
        ]

        self.cmd.run(cmd_find_common)
        if os.path.getsize(common_reads) == 0:
            self.logger.info("No common reads found")
            return

        self.r1.read_filter_inplace(self.r1.current, common_reads)
        self.r2.read_filter_inplace(self.r2.current, common_reads)

    def clean_unpaired(self):
        """
        use seqtk to find common reads. i) use seqkit to find common reads, ii) keep those.
        iii) run trimmomatic with minimal arguments to sort reads and remove unpaired reads.
        """

        read_dir = os.path.dirname(self.r1)
        common_reads = os.path.join(read_dir, "common_reads.tsv")
        temp_name = os.path.join(read_dir, "temp1")
        temp1 = temp_name + ".fq.gz"
        temp2 = os.path.join(read_dir, "temp2.fq.gz")

        find_common_seqtk_cmd = [
            "seqtk",
            "subseq",
            self.r1,
            self.r2,
            "|",
            "paste",
            "- - - -",
            "|",
            "cut",
            "-f1",
            "|",
            "uniq",
            "|",
            "sed",
            "'s/^@//g'",
            ">",
            common_reads,
        ]

        subseq_r1_cmd = [
            "seqtk",
            "subseq",
            self.r1,
            common_reads,
            "|",
            "gzip",
            ">",
            temp1,
        ]

        subseq_r2_cmd = [
            "seqtk",
            "subseq",
            self.r2,
            common_reads,
            "|",
            "gzip",
            ">",
            temp2,
        ]

        mv_r1 = [
            "mv",
            temp1,
            self.r1,
        ]

        mv_r2 = [
            "mv",
            temp2,
            self.r2,
        ]

        trimmomatic_cmd = [
            "trimmomatic",
            "PE",
            "-phred33",
            "-threads",
            self.threads,
            self.r1,
            self.r2,
            "-baseout",
            temp1,
            "MINLEN:20",
        ]

        mv_trimmomatic_output_r1_cmd = ["mv", temp1, self.r1]

        mv_trimmomatic_output_r2_cmd = ["mv", temp2, self.r2]

        self.cmd.run(find_common_seqtk_cmd)
        if os.path.getsize(common_reads) > 0:
            self.cmd.run(subseq_r2_cmd)
            self.cmd.run(subseq_r1_cmd)
            os.system(" ".join(mv_r1))
            os.system(" ".join(mv_r2))
            self.cmd.run(trimmomatic_cmd)
            self.cmd.run(" ".join(mv_trimmomatic_output_r1_cmd))
            self.cmd.run(" ".join(mv_trimmomatic_output_r2_cmd))

        else:
            raise ValueError("No common reads found")
