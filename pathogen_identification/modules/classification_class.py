import csv
import logging
import os
import re
import shutil
from abc import ABC, abstractmethod
from random import randint
from typing import Any, Type

import pandas as pd

from pathogen_identification.constants_settings import ConstantsSettings
from pathogen_identification.modules.object_classes import RunCMD, SoftwareDetail


def read_sam_file(
    file: str, sep: str = "\t", columns_to_keep: list = None
) -> pd.DataFrame:
    """
    read file line by line, keep only lines that do not start with @, return dataframe.

    """

    with open(file, "r") as f:
        line = f.readline()
        while line.startswith("@"):
            line = f.readline()

        if columns_to_keep:
            columns_to_keep = [int(c) for c in columns_to_keep]
            data = [line.split(sep)[c] for c in columns_to_keep]
        else:
            data = line.split(sep)
            columns_to_keep = list(range(len(data)))

        data = [data]
        for line in f:
            data.append([line.split(sep)[c] for c in columns_to_keep])

    return pd.DataFrame(data, columns=columns_to_keep)


def check_report_empty(file, comment="@"):
    if not os.path.exists(file):
        return True

    with open(file, "r") as f:
        lines = f.readlines()

    lines = [l for l in lines if not l.startswith(comment)]
    if len(lines) == 0:
        return True
    else:
        return False


class Classifier_init(ABC):
    """
    list of clasifiers that overwrite Classifier_init:
    - run_kaiju
    - run_CLARK
    - run_FastViromeExplorer

    """

    method_name = ""
    report_suffix = ""
    full_report_suffix = ""

    def __init__(
        self,
        db_path: str,
        query_path: str,
        out_path: str,
        args="",
        r2: str = "",
        prefix: str = "",
        bin: str = "",
        log_dir="",
    ):
        self.db_path = db_path
        self.query_path = query_path
        self.out_path = out_path
        self.r2 = r2
        self.prefix = prefix
        self.args = args
        self.cmd = RunCMD(bin, logdir=log_dir, prefix=prefix, task="classification")

        self.report_path = os.path.join(self.out_path, self.prefix + self.report_suffix)
        self.full_report_path = os.path.join(
            self.out_path, self.prefix + self.full_report_suffix
        )

    def filter_samfile_read_names(self, same=True, output_sam="", sep=",", idx=0):
        if not output_sam:
            output_sam = os.path.join(self.out_path, f"temp{randint(1,1999)}.sam")

        read_name_filter_regex = re.compile("^[A-Za-z0-9_-]*$")  # (r"@|=&$\t")

        with open(self.report_path, "r") as f:
            column_number = 0
            with open(output_sam, "w") as f2:
                for line in f:
                    if line.startswith(tuple(["@HD", "@SQ", "@PG", "@RG", "@CO"])):
                        f2.write(line)
                    line_columns = len(line.split(sep))

                    if column_number == 0:
                        column_number = line_columns

                    if line_columns > column_number:
                        continue

                    if read_name_filter_regex.search(
                        line.split(sep)[idx].replace(":", "_")
                    ):
                        f2.write(line)

        if same:
            os.remove(self.report_path)
            os.rename(output_sam, self.report_path)

    def filter_secondary_alignments(self):
        """
        Filter secondary alignments from bam file."""
        cmd = [
            "samtools",
            "view",
            "-F0x900",
            self.report_path,
            ">",
            self.report_path.replace(".sam", ".filtered.sam"),
        ]

        try:
            self.cmd.run(cmd)
            if os.path.isfile(self.report_path.replace(".sam", ".filtered.sam")):
                os.remove(self.report_path)
                shutil.copy(
                    self.report_path.replace(".sam", ".filtered.sam"), self.report_path
                )
        except:
            pass

    @abstractmethod
    def run_SE(self, threads: str = "3"):
        pass

    @abstractmethod
    def run_PE(self, threads: str = "3"):
        pass

    @abstractmethod
    def get_report(self) -> pd.DataFrame:
        pass

    @abstractmethod
    def get_report_simple(self) -> pd.DataFrame:
        pass


class run_kaiju(Classifier_init):
    method_name = "kaiju"
    report_suffix = ".kaiju.out"
    full_report_suffix = ".kaiju"
    classification_type = "viral"

    def __init__(
        self,
        db_path: str,
        query_path: str,
        out_path: str,
        args="",
        r2: str = "",
        prefix: str = "",
        bin: str = "",
        log_dir="",
    ):
        super().__init__(db_path, query_path, out_path, args, r2, prefix, bin, log_dir)
        self.report_path = os.path.join(self.out_path, self.prefix + self.report_suffix)

        try:
            dbdir = os.path.dirname(self.db_path)

            self.nodes = os.path.join(dbdir, "nodes.dmp")
            if not os.path.isfile(self.nodes):
                raise FileNotFoundError(
                    "nodes.dmp not found, tried looking in db directory."
                )
        except FileNotFoundError:
            raise FileNotFoundError("Kaiju nodes file not found")

    def run_SE(self, threads: int = 3):
        cmd = f"kaiju -t {self.nodes} -f {self.db_path} -i {self.query_path} -o {self.report_path} -z {threads} {self.args}"
        self.cmd.run(cmd)

    def run_PE(self, threads: int = 3):
        cmd = f"kaiju -t {self.nodes} -f {self.db_path} -i {self.query_path} -j {self.r2} -o {self.report_path} -z {threads} {self.args}"
        self.cmd.run(cmd)

    def get_report(self) -> pd.DataFrame:
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(self.report_path, sep="\t", header=None).rename(
            columns={
                0: "CU",
                1: "seqid",
                2: "taxid",
                3: "length",
                4: "LCAmap",
            }
        )

    def get_report_simple(self) -> pd.DataFrame:
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(
            self.report_path, sep="\t", header=None, usecols=[1, 2], comment="@"
        ).rename(columns={1: "qseqid", 2: "taxid"})


class run_FastViromeExplorer(Classifier_init):
    method_name = "FastViromeExplorer"
    report_suffix = ".sam"
    report_name = "FastViromeExplorer-reads-mapped-sorted.sam"

    full_report_suffix = ".tsv"
    full_report_name = "abundance.tsv"

    def __init__(
        self,
        db_path: str,
        query_path: str,
        out_path: str,
        args="",
        r2: str = "",
        prefix: str = "",
        bin: str = "",
        log_dir: str = "",
    ):
        super().__init__(db_path, query_path, out_path, args, r2, prefix, bin, log_dir)

        self.report_path = os.path.join(self.out_path, self.report_name)
        self.full_report_path = os.path.join(self.out_path, self.full_report_name)

    def run_SE(self, threads: int = 3):
        cmd = [
            "FastViromeExplorer",
            self.args,
            "-1",
            self.query_path,
            "-i",
            self.db_path,
            "-l",
            self.db_path.replace(".idx", "-list.txt"),
            "-o",
            self.out_path,
        ]

        self.cmd.run_java(cmd)

    def run_PE(self, threads: int = 3):
        cmd = [
            "FastViromeExplorer",
            self.args,
            "-1",
            self.query_path,
            "-2",
            self.r2,
            "-i",
            self.db_path,
            "-l",
            self.db_path.replace(".idx", "-list.txt"),
            "-o",
            self.out_path,
        ]

        self.cmd.run_java(cmd)

    def get_report(self) -> pd.DataFrame:
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(self.report_path, sep="\t", header=None).rename(
            columns={
                0: "qseqid",
                1: "flag",
                2: "acc",
                3: "pos",
                4: "mapq",
                5: "cigar",
                6: "rnext",
                7: "pnext",
                8: "tlen",
                9: "seq",
                10: "qual",
            }
        )

    def get_report_simple(self) -> pd.DataFrame:
        """
        read classifier output, return only query and reference sequence id columns.
        """

        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(
            self.report_path, sep="\t", header=None, usecols=[0, 2], comment="@"
        ).rename(columns={0: "qseqid", 2: "acc"})


class run_CLARK(Classifier_init):
    method_name = "CLARK"
    report_name = "results.csv"
    report_suffix = ".csv"

    full_report_name = "CLARK"
    full_report_suffix = ".tsv"

    def __init__(
        self,
        db_path: str,
        query_path: str,
        out_path: str,
        args="",
        r2: str = "",
        prefix: str = "",
        bin: str = "",
        log_dir: str = "",
    ):
        super().__init__(db_path, query_path, out_path, args, r2, prefix, bin, log_dir)

        self.report_path = os.path.join(self.out_path, self.report_name)
        self.cmd.bin = self.cmd.bin.replace("/bin", "")

    def run_SE(self, threads: int = 3):
        cmd = f"classify_metagenome.sh -O {self.query_path} -R {os.path.splitext(self.report_path)[0]} --gzipped {self.args} -n {threads}"
        self.cmd.run(cmd)

    def run_PE(self, threads: int = 3):
        cmd = f"classify_metagenome.sh -P {self.query_path} {self.r2} -R {os.path.splitext(self.report_path)[0]} --gzipped {self.args}  -n {threads}"
        self.cmd.run(cmd)

    def read_report(self):
        self.filter_samfile_read_names(sep=",")
        report = pd.read_csv(self.report_path, sep=",")
        report.columns = report.columns.str.strip()
        report = report.rename(
            columns={
                "Object_ID": "qseqid",
                "Length": "length",
                "Assignment": "taxid",
                "1st_assignment": "taxid",
            }
        )

        if "2nd_assignment" in report.columns:
            report["taxid"] = (
                report["taxid"].astype(str) + "," + report["2nd_assignment"].astype(str)
            )

            report = report.drop(columns=["2nd_assignment"])
            report["taxid"] = report["taxid"].apply(lambda x: x.split(","))
            report = report.explode("taxid").drop_duplicates(
                subset=["qseqid", "taxid"], keep="first"
            )
            report = report[report.taxid != "nan"]

        report = report.dropna().reset_index(drop=True)
        report["taxid"] = report["taxid"].astype(float).astype(int)
        return report

    def get_report(self) -> pd.DataFrame:
        return self.read_report()

    def get_report_simple(self) -> pd.DataFrame:
        """
        read classifier output, return only query and reference sequence id columns.
        """
        return self.read_report()


class run_blast(Classifier_init):
    method_name = "blast"
    report_suffix = ".blast_results.tsv"
    full_report_suffix = ".blast_full_results.tsv"

    def unzip_query(self, output_path: str = ""):
        """
        unzip query file.
        """
        cmd = [
            "gunzip",
            "-c",
            self.query_path,
            ">",
            output_path,
        ]

        self.cmd.run_bash(cmd)

    def run_SE(self, threads: int = 3, *args, **kwargs):
        """
        run single read file classification.
        """
        temp_read = os.path.join(self.out_path, self.prefix + ".temp_read.fasta")

        self.unzip_query(temp_read)
        db_dir = os.path.dirname(self.db_path)
        db_name = os.path.basename(self.db_path)
        cwd = os.getcwd()

        os.chdir(db_dir)

        cmd = [
            "blastn",
            "-db",
            db_name,
            "-query",
            temp_read,
            "-outfmt",
            "6",
            "-out",
            self.report_path,
            "-num_threads",
            str(threads),
            self.args,
        ]

        self.cmd.run(cmd)

        os.chdir(cwd)

    def run_PE(self, threads: int = 3, *args, **kwargs):
        self.run_PE(threads)

    def get_report(self) -> pd.DataFrame:
        """
        read classifier output. return pandas dataframe with standard query sequence id and accession column names.
        """

        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(self.report_path, sep="\t", header=None).rename(
            columns={
                0: "qseqid",
                1: "acc",
                2: "pident",
                3: "length",
                4: "mismatch",
                5: "gapopen",
                6: "qstart",
                7: "qend",
                8: "sstart",
                9: "send",
                10: "evalue",
                11: "bitscore",
            }
        )

    def get_report_simple(self) -> pd.DataFrame:
        """
        read classifier output, return only query and reference sequence id columns.
        """
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(
            self.report_path, sep="\t", header=None, usecols=[0, 1]
        ).rename(
            columns={
                0: "qseqid",
                1: "acc",
            }
        )


class run_blast_p(Classifier_init):
    method_name = "blast_p"
    report_suffix = ".blast_results.tsv"
    full_report_suffix = ".blast_results_full.tsv"

    def unzip_query(self, output_path: str = ""):
        """
        unzip query file.
        """
        cmd = [
            "gunzip",
            "-c",
            self.query_path,
            ">",
            output_path,
        ]

        self.cmd.run(cmd)

    def run_PE(self, threads: int = 3):
        pass

    def run_SE(self, threads: int = 3):
        """
        run single read file classification.
        """
        temp_read = os.path.join(self.out_path, self.prefix + ".temp_read.fasta")

        self.unzip_query(temp_read)
        cmd = [
            "blastp",
            "-db",
            self.db_path,
            "-query",
            temp_read,
            "-outfmt",
            "6",
            "-out",
            self.report_path,
            "-num_threads",
            str(threads),
            self.args,
        ]

        self.cmd.run(cmd)
        os.remove(temp_read)

    def get_report(self) -> pd.DataFrame:
        """
        read classifier output. return pandas dataframe with standard query sequence id and accession column names.
        """
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(self.report_path, sep="\t", header=None).rename(
            columns={
                0: "qseqid",
                1: "acc",
                2: "pident",
                3: "length",
                4: "mismatch",
                5: "gapopen",
                6: "qstart",
                7: "qend",
                8: "sstart",
                9: "send",
                10: "evalue",
                11: "bitscore",
            }
        )

    def get_report_simple(self) -> pd.DataFrame:
        """
        read classifier output, return only query and reference sequence id columns.
        """
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(
            self.report_path, sep="\t", header=None, usecols=[0, 1]
        ).rename(
            {
                0: "qseqid",
                1: "acc",
            }
        )


import json

from constants.constants import Televir_Metadata_Constants


class run_voyager(Classifier_init):
    method_name = "voyager"
    report_suffix = ".json"
    full_report_suffix = ".json"

    @staticmethod
    def extract_json_taxonomy(json_file: str):
        """
        Poorly structured json file, need to extract 'taxonomy' field."""
        lines_to_keep = ["{"]
        with open(json_file, "r") as f:
            line = f.readline()
            keep = False
            while line:
                if "taxonomy" in line:
                    keep = True
                if "]" in line:
                    keep = False
                if keep:
                    lines_to_keep.append(line.replace("\t", ""))
                line = f.readline()
        if keep:
            lines_to_keep.append("]}")
        else:
            lines_to_keep.append("}")
        lines_to_keep = "".join(lines_to_keep)
        data = json.loads(lines_to_keep)
        return data

    def parse_json_taxonomy(self, json_file: str) -> pd.DataFrame:
        """
        parse voyager json output file.
        """
        order = [
            "species",
            "subgenus",
            "genus",
            "family",
            "suborder",
            "order",
            "class",
            "phylum",
            "kingdom",
            "clade",
        ]
        data = self.extract_json_taxonomy(json_file)
        if "taxonomy" not in data:
            return pd.DataFrame(columns=["taxid", "description"])
        taxonomy = data["taxonomy"]
        output = []
        for taxa in taxonomy:
            for clade in order:
                if clade in taxa.keys():
                    output.append(taxa[clade])
                    break
        output = pd.DataFrame(output)
        output.columns = ["taxid", "description"]

        output["taxid"] = output["taxid"].astype(str)
        output["description"] = output["description"].astype(str)
        return output

    def run_SE(self, threads: int = 3):
        televir_dirs = Televir_Metadata_Constants()
        voyager_bindir = televir_dirs.get_software_bin_directory("voyager")

        voyager_dir = voyager_bindir.replace("/bin", "")
        print(self.db_path)
        print(self.args)
        cmd = [
            os.path.join(voyager_dir, "voyager-cli"),
            "-x",
            self.db_path,
            "--input",
            self.query_path,
            "--output",
            self.report_path,
        ]

        self.cmd.run_bash(cmd)

    def run_PE(self, threads: int = 3):
        televir_dirs = Televir_Metadata_Constants()
        voyager_bindir = televir_dirs.get_software_bin_directory("voyager")

        voyager_dir = voyager_bindir.replace("/bin", "")
        cmd = [
            os.path.join(voyager_dir, "voyager-cli"),
            "-x",
            self.db_path,
            "--input",
            self.query_path,
            "--output",
            self.report_path,
        ]

        self.cmd.run_bash(cmd)

    def get_report(self) -> pd.DataFrame:
        """
        read classifier output. return pandas dataframe with standard query sequence id and accession column names.
        """
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        report = self.parse_json_taxonomy(self.report_path)

        return report

    def get_report_simple(self) -> pd.DataFrame:
        """
        read classifier output, return only query and reference sequence id columns.
        """
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        report = self.parse_json_taxonomy(self.report_path)

        return report


class run_metaphlan(Classifier_init):
    method_name = "metaphlan"
    report_suffix = ".tsv"
    full_report_suffix = ".tsv"

    def __post_init__(self):
        """
        Post initialization method to set up the report path.
        """
        self.parse_args_db_edits()

    @staticmethod
    def parse_metaphlan_output(file: str) -> pd.DataFrame:
        """
        parse metaphlan output file.
        """
        data = []
        with open(file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                line = line.strip().split("\t")
                if len(line) < 2:
                    continue
                taxonomy = line[0].split("|")
                taxids = line[1].split("|")
                tax_index = -1
                if "SGD" in taxonomy[tax_index]:
                    tax_index = -2
                data.append([taxonomy[tax_index], taxids[tax_index], line[2]])
        if len(data) == 0:
            return pd.DataFrame(columns=["description", "taxid", "abundance"])
        data = pd.DataFrame(data)
        data.columns = ["description", "taxid", "abundance"]
        data["taxid"] = data["taxid"].astype(str)
        data["description"] = data["description"].astype(str)
        data["abundance"] = data["abundance"].astype(float)
        data = data[data["taxid"] != "0"]
        data = data[data["taxid"] != ""]
        return data

    @staticmethod
    def parse_metaphlan_output_as_tree(file: str) -> pd.DataFrame:
        nodes = {}
        edges = {}
        with open(file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                line = line.strip().split("\t")
                if len(line) < 2:
                    continue
                taxonomy = line[0].split("|")
                taxids = line[1].split("|")
                for i in range(len(taxonomy)):
                    if taxonomy[i] == "unclassified":
                        continue
                    if "SGD" in taxonomy[i]:
                        continue
                    if taxids[i] == "0":
                        continue
                    if taxids[i] not in nodes and taxids[i] != "":
                        nodes[taxids[i]] = (taxids[i], taxonomy[i], i)
                if len(taxids) < 2:
                    continue
                for i in range(len(taxids) - 1):
                    if taxonomy[i] == "unclassified":
                        continue
                    if "SGD" in taxonomy[i] or "SGD" in taxonomy[i + 1]:
                        continue
                    if taxids[i] == "0" or taxids[i + 1] == "0" or taxids[i] == "":
                        continue
                    if taxids[i + 1] == "":
                        continue
                    if taxids[i] not in edges:
                        edges[taxids[i]] = []
                    edges[taxids[i]].append((taxids[i + 1], taxonomy[i + 1], i + 1))
        if len(nodes) > 0:
            leaves = [nodes[taxid] for taxid in nodes.keys() if taxid not in edges]
            return pd.DataFrame(leaves, columns=["taxid", "description", "level"])
        return pd.DataFrame(columns=["taxid", "description", "level"])

    def parse_args_db_edits(self):

        args = self.args.split(" ")
        where_db_edits = -1
        if "--edits" not in self.args:
            return
        where_db_edits = self.args.split(" ").index("--edits")

        edits = args[where_db_edits + 1]
        edits = edits.split(";")
        args = args[:where_db_edits] + args[where_db_edits + 2 :]
        for edit in edits:
            args.append(edit)

        # args.remove("--edits")

        self.args = " ".join(args)

    def run_SE(self, threads: int = 3):
        televir_constants = Televir_Metadata_Constants()
        metaphlan_bidir = televir_constants.get_software_bin_directory("metaphlan")
        self.parse_args_db_edits()

        cmd = [
            os.path.join(metaphlan_bidir, "metaphlan"),
            self.query_path,
            "--input_type",
            "fastq",
            "--nproc",
            str(threads),
            "--index",
            os.path.basename(self.db_path).split(".")[0],
            "--bowtie2db",
            os.path.dirname(os.path.dirname(self.db_path)),
            self.args,
            "--bowtie2out",
            self.report_path.replace(".tsv", ".bowtie2.bam"),
            "--output_file",
            self.report_path,
        ]

        self.cmd.run_script(cmd, conda_env=os.path.dirname(metaphlan_bidir))

    def run_PE(self, threads: int = 3):
        televir_constants = Televir_Metadata_Constants()
        metaphlan_bidir = televir_constants.get_software_bin_directory("metaphlan")
        self.parse_args_db_edits()
        cmd_r1 = [
            os.path.join(metaphlan_bidir, "metaphlan"),
            self.query_path,
            "--input_type",
            "fastq",
            "--nproc",
            str(threads),
            self.args,
            "--index",
            os.path.basename(self.db_path).split(".")[0],
            "--bowtie2db",
            os.path.dirname(os.path.dirname(self.db_path)),
            "--bowtie2out",
            self.report_path.replace(".tsv", ".bowtie2.bz2"),
            ">",
            f"{self.report_path}.r1",
        ]

        cmd_r2 = [
            os.path.join(metaphlan_bidir, "metaphlan"),
            self.r2,
            "--input_type",
            "fastq",
            "--nproc",
            str(threads),
            self.args,
            "--index",
            os.path.basename(self.db_path).split(".")[0],
            "--bowtie2db",
            os.path.dirname(os.path.dirname(self.db_path)),
            "--bowtie2out",
            self.report_path.replace(".tsv", ".bowtie2.bz2"),
            ">",
            f"{self.report_path}.r2",
        ]

        self.cmd.run_script(cmd_r1, conda_env=os.path.dirname(metaphlan_bidir))
        self.cmd.run_script(cmd_r2, conda_env=os.path.dirname(metaphlan_bidir))

        merge_output_cmd = [
            os.path.join(metaphlan_bidir, "python3"),
            os.path.join(metaphlan_bidir, "merge_metaphlan_tables.py"),
            f"{self.report_path}.r1",
            f"{self.report_path}.r2",
            ">",
            self.report_path,
        ]

        self.cmd.run_python(merge_output_cmd)

    def get_report(self) -> pd.DataFrame:
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["description", "taxid", "abundance"])

        report = self.parse_metaphlan_output_as_tree(self.report_path)

        return report

    def get_report_simple(self) -> pd.DataFrame:
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        report = self.parse_metaphlan_output_as_tree(self.report_path)

        return report


class run_centrifuge(Classifier_init):
    method_name = "centrifuge"
    report_suffix = ".report.tsv"
    full_report_suffix = ".centrifuge"

    def run_SE(self, threads: int = 3, *args, **kwargs):
        """
        run single read file classification.
        """
        cmd = [
            "centrifuge",
            "-x",
            self.db_path,
            "--threads",
            f"{threads}",
            "-U",
            self.query_path,
            "-S",
            self.full_report_path,
            "--report-file",
            self.report_path,
            "--out-fmt",
            "sam",
            "--time",
            self.args,
        ]

        if kwargs:
            for key, value in kwargs.items():
                cmd.append(f"--{key}")
                cmd.append(f"{value}")

        self.cmd.run(cmd)

    def run_PE(self, threads: int = 3):
        """
        read classifier output. return pandas dataframe with standard query sequence id and accession column names.
        """
        cmd = [
            "centrifuge",
            "-x",
            self.db_path,
            "--threads",
            f"{threads}",
            "-1",
            self.query_path,
            "-2",
            self.r2,
            "-S",
            self.full_report_path,
            "--report-file",
            self.report_path,
            "--out-fmt",
            "sam",
            "--time",
            self.args,
        ]

        self.cmd.run(cmd)

    def get_report(self) -> pd.DataFrame:
        """
        read classifier output. return pandas dataframe with standard query sequence id and accession column names.
        """
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(self.report_path, sep="\t").rename(
            columns={"numReads": "counts"}
        )

    def get_report_simple(self) -> pd.DataFrame:
        """
        read classifier output, return only query and reference sequence id columns.
        """
        if check_report_empty(self.full_report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        report = pd.read_csv(
            self.full_report_path, sep="\t", header=None, usecols=[0, 2, 6]
        ).rename(columns={0: "qseqid", 2: "taxid", 6: "acc"})

        report = report[report.acc != "unclassified"][["qseqid", "taxid"]]
        report = report[report.taxid != 0][["qseqid", "taxid"]]

        return pd.DataFrame(report)


class run_deSamba(Classifier_init):
    method_name = "deSamba"
    report_suffix = ".sam"
    full_report_suffix = ".sam"

    def run_SE(self, threads: int = 3, **kwargs):
        """
        run single read file classification.
        """
        cmd = f"desamba classify  {self.query_path} -o {self.out_path}"
        cmd = [
            "deSAMBA",
            "classify",
            "-t",
            f"{threads}",
            self.db_path,
            self.query_path,
            "-o",
            self.report_path,
            self.args,
        ]

        if kwargs:
            for key, value in kwargs.items():
                cmd.append(f"--{key}")
                cmd.append(f"{value}")

        self.cmd.run(cmd)

    def run_PE(self, threads: int = 3, **kwargs):
        pass

    def get_report(self) -> pd.DataFrame:
        """
        read classifier output. return pandas dataframe with standard query sequence id and accession column names.
        """
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(self.report_path, sep="\t", header=None).rename(
            columns={
                0: "qseqid",
                2: "acc",
                3: "sstart",
                4: "length",
                8: "qual",
            }
        )

    def get_report_simple(self) -> pd.DataFrame:
        """
        read classifier output, return only query and reference sequence id columns.
        """
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        report = pd.read_csv(
            self.report_path, sep="\t", header=None, usecols=[0, 2, 4, 8], comment="@"
        ).rename(columns={0: "qseqid", 2: "acc", 4: "qual", 8: "length"})
        report = report[report.acc != "*"][["qseqid", "acc"]]

        return report


class run_kraken2(Classifier_init):
    method_name = "kraken2"
    report_suffix = ".tsv"
    full_report_suffix = ".tsv"

    def __init__(
        self,
        db_path: str,
        query_path: str,
        out_path: str,
        args="",
        r2: str = "",
        prefix: str = "",
        bin: str = "",
        log_dir="",
    ):
        super().__init__(db_path, query_path, out_path, args, r2, prefix, bin, log_dir)
        self.args = self.args.replace("--quick OFF", "")
        self.args = self.args.replace("--quick ON", "--quick")

    def run_SE(self, threads: int = 3, **kwargs):
        """
        run single read file classification.
        """
        cmd = [
            "kraken2",
            "--threads",
            f"{threads}",
            "--db",
            self.db_path,
            "--gzip-compressed",
            "--output",
            self.report_path,
            self.args,
        ]

        if kwargs:
            for k, g in kwargs.items():
                cmd.append(f"--{k}")
                cmd.append(g)

        cmd.append(self.query_path)

        self.cmd.run_script_software(cmd)

    def run_PE(self, threads: int = 3, **kwargs):
        """
        run paired read files classification.
        """
        cmd = f"kraken2 --threads 4 --db {self.db_path} --fastq-input --gzip-compressed --output {self.out_path} {self.query_path} {self.r2}"
        cmd = [
            "kraken2",
            "--threads",
            f"{threads}",
            "--db",
            self.db_path,
            "--fastq-input",
            "--gzip-compressed",
            "--output",
            self.report_path,
            self.args,
        ]

        if kwargs:
            for k, g in kwargs.items():
                cmd.append(f"--{k}")
                cmd.append(g)

        cmd.append(self.query_path)
        cmd.append(self.r2)

        self.cmd.run_script_software(cmd)

    def get_report(self) -> pd.DataFrame:
        """
        read classifier output. return pandas dataframe with standard query sequence id and accession column names.
        """
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(self.report_path, sep="\t", header=None).rename(
            columns={
                0: "CU",
                1: "seqid",
                2: "taxid",
                3: "length",
                4: "LCAmap",
            }
        )

    def get_report_simple(self) -> pd.DataFrame:
        """
        read classifier output, return only query and reference sequence id columns.
        """
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        report = pd.read_csv(
            self.report_path, sep="\t", header=None, usecols=[0, 1, 2, 3], comment="@"
        ).rename(columns={0: "status", 1: "qseqid", 2: "taxid", 3: "length"})
        report = report[report.status != "U"][
            ["qseqid", "taxid", "length"]
        ]  # remove unclassified

        return report


class run_diamond(Classifier_init):
    method_name = "diamond"
    report_suffix = ".daa"
    full_report_suffix = ".daa"

    def __init__(
        self,
        db_path: str,
        query_path: str,
        out_path: str,
        args="",
        r2: str = "",
        prefix: str = "",
        bin: str = "",
        log_dir="",
    ):
        super().__init__(db_path, query_path, out_path, args, r2, prefix, bin, log_dir)
        self.process_args()

    def process_args(self):
        """
        Process arguments to remove preset and mode option flags"""
        self.args = self.args.replace("sensitivity", "")

    def run_SE(self, threads: int = 3, *args, **kwargs):
        cmd = [
            "diamond",
            "blastx",
            "-p",
            f"{threads}",
            "-d",
            self.db_path,
            "-q",
            self.query_path,
            "-o",
            self.report_path,
            self.args,
        ]

        self.cmd.run(cmd)

    def run_PE(self, threads: int = 3, *args, **kwargs):
        cmd = [
            "diamond",
            "blastx",
            "-p",
            f"{threads}",
            "-d",
            self.db_path,
            "-q",
            self.query_path,
            "-o",
            self.report_path,
            self.args,
        ]

        self.cmd.run(cmd)

    @staticmethod
    def acc_name_simplify(acc: str):
        if "|" in acc:
            acc = acc.split("|")[2]

        return acc

    def get_report(self) -> pd.DataFrame:
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        report = pd.read_csv(self.report_path, sep="\t", header=None).rename(
            columns={
                0: "qseqid",
                1: "prot_acc",
                2: "pident",
                3: "length",
                4: "mismatch",
                5: "gapopen",
                6: "qstart",
                7: "qend",
                8: "sstart",
                9: "send",
                10: "evalue",
                11: "bitscore",
            }
        )
        report["prot_acc"] = report["prot_acc"].apply(self.acc_name_simplify)
        return report

    def get_report_simple(self) -> pd.DataFrame:
        """
        read classifier output, return only query and reference sequence id columns.
        """
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        report = pd.read_csv(
            self.report_path, sep="\t", header=None, usecols=[0, 1], comment="@"
        ).rename(columns={0: "qseqid", 1: "prot_acc"})
        report["prot_acc"] = report["prot_acc"].apply(self.acc_name_simplify)
        return report


class run_krakenuniq(Classifier_init):
    method_name = "krakenuniq"
    report_suffix = ".tsv"
    full_report_suffix = ".krakenuniq"

    def run_SE(self, threads: int = 3):
        cmd = [
            "krakenuniq",
            "--db",
            self.db_path,
            "--threads",
            f"{threads}",
            "--gzip-compressed",
            "--fastq-input",
            self.args,
            "--output",
            self.report_path,
            self.query_path,
        ]

        self.cmd.run(cmd)

    def run_PE(self, threads: int = 3):
        cmd = [
            "krakenuniq",
            "--db",
            self.db_path,
            "--threads",
            f"{threads}",
            "--gzip-compressed",
            "--fastq-input",
            self.args,
            "--output",
            self.report_path,
            self.query_path,
            self.r2,
        ]
        self.cmd.run(cmd)

    def get_report(self) -> pd.DataFrame:
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(self.report_path, sep="\t", header=None).rename(
            columns={
                0: "CU",
                1: "qseqid",
                2: "protid",
                3: "length",
                4: "LCAmap",
            }
        )

    def get_report_simple(self) -> pd.DataFrame:
        """
        read classifier output, return only query and reference sequence id columns.
        """

        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(
            self.report_path, sep="\t", header=None, usecols=[1, 2], comment="@"
        ).rename(columns={1: "qseqid", 2: "protid"})


class run_minimap2_illumina(Classifier_init):
    method_name = "minimap2_illumina"
    report_suffix = ".sam"
    full_report_suffix = ".minimap2_illumina"

    def run_SE(self, threads: int = 3):
        cmd = f"minimap2 -a -t {threads} {self.args} {self.db_path} {self.query_path} > {self.report_path}"
        self.cmd.run(cmd)

    def run_PE(self, threads: int = 3):
        cmd = f"minimap2 -a -t {threads} {self.args} {self.db_path} {self.query_path} {self.r2} > {self.report_path}"
        self.cmd.run(cmd)

    def get_report(self) -> pd.DataFrame:
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(self.report_path, sep="\t", header=None).rename(
            columns={
                0: "acc",
                1: "flag",
                2: "qseqid",
                3: "pos",
                4: "mapq",
                5: "cigar",
                6: "rnext",
                7: "pnext",
                8: "tlen",
                9: "seq",
                10: "qual",
            }
        )

    def get_report_simple(self) -> pd.DataFrame:
        """
        read classifier output, return only query and reference sequence id columns.
        """

        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        report = pd.read_csv(
            self.report_path, sep="\t", header=None, usecols=[0, 2], comment="@"
        ).rename(columns={0: "qseqid", 2: "acc"})

        report = report[report["acc"] != "*"]
        return report


class run_bowtie2(Classifier_init):
    method_name = "bowtie2"
    report_suffix = ".sam"
    full_report_suffix = ".bowtie2"

    def run_SE(self, threads: int = 3):
        index_name = self.index_reference()
        self.process_arguments()

        cmd = f"bowtie2 -t --sam-nohead --sam-nosq -x {index_name} -U {self.query_path} {self.args} -S {self.report_path}"
        self.cmd.run(cmd)

    def run_PE(self, threads: int = 3):
        index_name = self.index_reference()
        self.process_arguments()

        cmd = f"bowtie2 -t --sam-nohead --sam-nosq -x {index_name} -1 {self.query_path} -2 {self.r2} {self.args} -S {self.report_path}"
        self.cmd.run(cmd)

    def index_reference(self):
        """
        Index reference using bowtie2."""
        index_name = os.path.splitext(self.db_path)[0] + "_index"

        cmd = ["bowtie2-build", self.db_path, index_name]

        self.cmd.run(cmd)

        return index_name

    def process_arguments(self):
        """
        Process arguments to remove preset and mode option flags"""
        self.args = self.args.replace("[preset]", "").replace("[mode]", "")

    def get_report(self) -> pd.DataFrame:
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(self.report_path, sep="\t", header=None).rename(
            columns={
                0: "qseqid",
                1: "flag",
                2: "acc",
                3: "pos",
                4: "mapq",
                5: "cigar",
                6: "rnext",
                7: "pnext",
                8: "tlen",
                9: "seq",
                10: "qual",
            }
        )

    def get_report_simple(self) -> pd.DataFrame:
        """
        read classifier output, return only query and reference sequence id columns.
        """
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(
            self.report_path, sep="\t", header=None, usecols=[0, 2], comment="@"
        ).rename(columns={2: "acc", 0: "qseqid"})


class run_bwa_mem(Classifier_init):
    method_name = "bwa_illumina"
    report_suffix = ".sam"
    full_report_suffix = ".bwa_illumina"

    def run_SE(self, threads: int = 3):
        rundir = os.path.dirname(self.report_path)
        unzip_seq = f"gunzip -c {self.query_path} > {rundir}/seq.fq"

        self.cmd.run_bash(unzip_seq)

        cmd = f"bwa mem -t {threads} {self.args} {os.path.splitext(self.db_path)[0]} {rundir}/seq.fq > {self.report_path}"

        self.cmd.run(cmd)
        self.filter_secondary_alignments()

    def run_PE(self, threads: int = 3):
        rundir = os.path.dirname(self.report_path)
        unzip_seq = f"gunzip -c {self.query_path} > {rundir}/seq.fq"
        unzip_seq2 = f"gunzip -c {self.r2} > {rundir}/seq2.fq"

        self.cmd.run_bash(unzip_seq)
        self.cmd.run_bash(unzip_seq2)

        cmd = f"bwa mem -t {threads} {self.args} {os.path.splitext(self.db_path)[0]} {rundir}/seq.fq {rundir}/seq2.fq > {self.report_path}"

        self.cmd.run(cmd)
        self.filter_secondary_alignments()

    def get_report(self) -> pd.DataFrame:
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(self.report_path, sep="\t", header=None).rename(
            columns={
                0: "qseqid",
                1: "flag",
                2: "acc",
                3: "pos",
                4: "mapq",
                5: "cigar",
                6: "rnext",
                7: "pnext",
                8: "tlen",
                9: "seq",
                10: "qual",
            }
        )

    def get_report_simple(self) -> pd.DataFrame:
        """
        read classifier output, return only query and reference sequence id columns.
        """

        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        report = pd.read_csv(
            self.report_path,
            sep="\t",
            header=None,
            usecols=[0, 2],
            comment="@",
            quoting=csv.QUOTE_NONE,
        ).rename(columns={0: "qseqid", 2: "acc"})

        report = report[report["acc"] != "*"]
        return report


class run_bwa_mem_iterative(Classifier_init):
    method_name = "bwa_illumina_iterative"
    report_suffix = ".lst"
    full_report_suffix = ".bwa_illumina_iterative"

    def generate_tmp_file_name(self):
        """
        Generate a temporary file name
        """
        from random import randint

        seed = randint(1, 100000)
        tempdir = self.out_path
        tempname = f"temp_subsample_{seed}"

        temp1 = os.path.join(tempdir, tempname)

        return temp1

    def filter_mapped_fastq_using_bam(self, bam_file):
        """
        filter mapped fastq using bam file
        """

        tmp_read = self.generate_tmp_file_name()
        tmp_read1 = tmp_read + ".lst"

        cmd = [
            "samtools",
            "view",
            "-h",
            "-F",
            "4",
            bam_file,
            "|",
            "grep -v '^@'",
            "|",
            "cut -f1",
            "|",
            "sort",
            "|",
            "uniq",
            ">",
            tmp_read1,
        ]

        self.cmd.run_script_software(cmd)

        return tmp_read1

    def run_SE(self, threads: int = 3):
        rundir = os.path.dirname(self.report_path)
        unzip_seq = f"gunzip -c {self.query_path} > {rundir}/seq.fq"

        self.cmd.run_bash(unzip_seq)

        db_list = self.db_path.split(";")
        db_list = [os.path.splitext(db)[0] for db in db_list]
        reads_to_keep = set()
        for db in db_list:

            cmd = f"bwa mem -t {threads} {self.args} {db} {rundir}/seq.fq > {self.report_path}.sam"
            self.cmd.run(cmd)
            mapped_reads = self.filter_mapped_fastq_using_bam(f"{self.report_path}.sam")
            with open(mapped_reads, "r") as f:
                for line in f:
                    line = line.strip()
                    reads_to_keep.add(line)

        with open(self.report_path, "w") as f:
            for read in reads_to_keep:
                f.write(f"{read}\n")

    def run_PE(self, threads: int = 3):
        rundir = os.path.dirname(self.report_path)
        unzip_seq = f"gunzip -c {self.query_path} > {rundir}/seq.fq"
        unzip_seq2 = f"gunzip -c {self.r2} > {rundir}/seq2.fq"

        self.cmd.run_bash(unzip_seq)
        self.cmd.run_bash(unzip_seq2)

        db_list = self.db_path.split(";")
        db_list = [os.path.splitext(db)[0] for db in db_list]
        reads_to_keep = set()
        for db in db_list:
            cmd = f"bwa mem -t {threads} {self.args} {db} {rundir}/seq.fq {rundir}/seq2.fq > {self.report_path}.sam"
            self.cmd.run(cmd)
            mapped_reads = self.filter_mapped_fastq_using_bam(f"{self.report_path}.sam")
            with open(mapped_reads, "r") as f:
                for line in f:
                    line = line.strip()
                    reads_to_keep.add(line)
        with open(self.report_path, "w") as f:
            for read in reads_to_keep:
                f.write(f"{read}\n")

    def get_report(self) -> pd.DataFrame:

        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        report = pd.read_csv(self.report_path, sep="\t", header=None).rename(
            columns={0: "qseqid"}
        )

        report["qseqid"] = report["qseqid"].apply(lambda x: x.split("/")[0])
        report["qseqid"] = report["qseqid"].apply(lambda x: x.split(" ")[0])
        return report

    def get_report_simple(self) -> pd.DataFrame:

        return self.get_report()


class run_bowtie2_ONT(Classifier_init):
    method_name = "bowtie2_ONT"
    report_suffix = ".sam"
    full_report_suffix = ".bowtie2_ONT"

    def run_SE(self, threads: int = 3):
        cmd = f"bowtie2 -a --threads {threads} --sam-nohead --sam-nosq --no-unal {self.args} -x {self.db_path} -U {self.query_path} -S {self.report_path}"
        self.cmd.run(cmd)

    def run_PE(self, threads: int = 3):
        cmd = f"bowtie2 -a --threads {threads} --sam-nohead --sam-nosq --no-unal {self.args} -x {self.db_path} -1 {self.query_path} -2 {self.r2} -S {self.report_path}"
        self.cmd.run(cmd)

    def get_report(self) -> pd.DataFrame:
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(self.report_path, sep="\t", header=None).rename(
            columns={
                0: "acc",
                1: "flag",
                2: "qseqid",
                3: "pos",
                4: "mapq",
                5: "cigar",
                6: "rnext",
                7: "pnext",
                8: "tlen",
                9: "seq",
                10: "qual",
            }
        )

    def get_report_simple(self) -> pd.DataFrame:
        """
        read classifier output, return only query and reference sequence id columns.
        """

        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        report = pd.read_csv(
            self.report_path, sep="\t", header=None, usecols=[0, 2], comment="@"
        ).rename(columns={1: "acc", 0: "qseqid"})

        report = report[report["acc"] != "*"]
        return report


class run_minimap2_ONT(Classifier_init):
    method_name = "minimap2_ONT"
    report_suffix = ".paf"
    full_report_suffix = ".minimap2"

    def run_SE(self, threads: int = 3):
        cmd = f"minimap2 -t {threads} {self.args} {self.db_path} {self.query_path} > {self.report_path}"
        self.cmd.run(cmd)

    def run_PE(self, threads: int = 3):
        cmd = f"minimap2 -t {threads} {self.args} {self.db_path} {self.query_path} {self.r2} > {self.report_path} "
        self.cmd.run(cmd)

    def get_report(self) -> pd.DataFrame:
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(self.report_path, sep="\t", header=None).rename(
            columns={
                0: "qseqid",
                1: "flag",
                2: "acc",
                3: "tstart",
                4: "qual",
                5: "cigar",
                6: "rnext",
                7: "pnext",
                8: "tlen",
                9: "seq",
                10: "qual",
            }
        )

    def get_report_simple(self) -> pd.DataFrame:
        """
        read classifier output, return only query and reference sequence id columns.
        """

        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        # report = pd.read_csv(
        #    self.report_path, sep="\t", header=None, usecols=[0, 2], comment="@"
        # ).rename(columns={0: "qseqid", 2: "acc"})

        report = read_sam_file(
            self.report_path, sep="\t", columns_to_keep=[0, 2]
        ).rename(columns={0: "qseqid", 2: "acc"})

        report = report[report["acc"] != "*"]
        return report


class run_minimap2_asm(Classifier_init):
    method_name = "minimap2_asm"
    report_suffix = ".paf"
    full_report_suffix = ".minimap2"

    def run_SE(self, threads: int = 3):
        cmd = f"minimap2 -t {threads} -c {self.args} {self.db_path} {self.query_path} > {self.report_path}"
        self.cmd.run(cmd)

    def run_PE(self, threads: int = 3):
        pass

    def get_report(self) -> pd.DataFrame:
        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(self.report_path, sep="\t", header=None).rename(
            columns={
                0: "qseqid",
                1: "qlen",
                2: "qstart",
                3: "qend",
                4: "strand",
                5: "acc",
                6: "tlen",
                7: "tstart",
                8: "tend",
                9: "nmatches",
                10: "length",
                11: "qual",
            }
        )

    def get_report_simple(self) -> pd.DataFrame:
        """
        read classifier output, return only query and reference sequence id columns.
        """

        if check_report_empty(self.report_path):
            return pd.DataFrame(columns=["qseqid", "acc"])

        return pd.read_csv(
            self.report_path, sep="\t", header=None, usecols=[0, 5, 10], comment="@"
        ).rename(columns={0: "qseqid", 5: "acc", 10: "length"})


class Empty_classifier(Classifier_init):
    method_name = "None"
    report_suffix = ".blast_results.tsv"
    full_report_suffix = ".blast_full_results.tsv"

    def run_SE(self, threads: int = 3):
        pass

    def run_PE(self, threads: int = 3):
        pass

    def get_report(self) -> pd.DataFrame:
        return pd.DataFrame(columns=["qseqid", "acc"])

    def get_report_simple(self) -> pd.DataFrame:
        return pd.DataFrame(columns=["qseqid", "acc"])


class Classifier:
    """
    Classifier class.
    """

    classifier: Classifier_init

    available_software: dict = {
        "blastn": run_blast,
        "blastp": run_blast_p,
        "centrifuge": run_centrifuge,
        "voyager": run_voyager,
        "desamba": run_deSamba,
        "kraken2": run_kraken2,
        "minimap2_illumina": run_minimap2_illumina,
        "minimap2_illu": run_minimap2_illumina,
        "minimap2": run_minimap2_ONT,
        "minimap2_ont": run_minimap2_ONT,
        "minimap2_asm": run_minimap2_asm,
        "minimap2_remap": run_minimap2_illumina,
        "diamond": run_diamond,
        "metaphlan": run_metaphlan,
        "kaiju": run_kaiju,
        "krakenuniq": run_krakenuniq,
        "fastviromeexplorer": run_FastViromeExplorer,
        "clark": run_CLARK,
        "bowtie2_remap": run_bowtie2,
        "bowtie2": run_bowtie2,
        "bwa": run_bwa_mem,
        "bwa-filter": run_bwa_mem_iterative,
    }

    def __init__(
        self,
        classifier_method: SoftwareDetail,
        query_path: str = "",
        type: str = ConstantsSettings.SINGLE_END,
        r2: str = "",
        prefix: str = "",
        threads: int = 1,
        bin: str = "",
        logging_level: int = logging.INFO,
        log_dir: str = "",
    ):
        """
        Initialize classifier.

        :param classifier_method: classifier method name
        :param query_path: query file path
        :param type: SE or PE
        :param r2: r2 file path
        :param prefix: prefix for output file
        :param threads: threads number
        :param bin: bin path
        :param logging_level: logging level
        """
        self.logger = logging.getLogger(f"{__name__}_{classifier_method.name}_{prefix}")
        if self.logger.hasHandlers():
            self.logger.handlers.clear()
        self.logger.propagate = False
        self.logger.setLevel(logging_level)
        self.logger.addHandler(logging.StreamHandler())
        self.log_dir = log_dir
        self.cmd = RunCMD(
            bin,
            logdir=self.log_dir,
            prefix=prefix,
            task=f"classification_{classifier_method.name}_{prefix}",
        )
        self.prefix = prefix
        self.classifier_method = classifier_method
        self.r1 = query_path
        self.r2 = r2
        self.output_path = classifier_method.output_dir
        self.threads = str(threads)
        self.report_name = self.classifier_method.name + "_results.tsv"
        self.classifier = self.classifier_configure()
        self.type = type
        self.classification_report = pd.DataFrame(columns=["qseqid", "acc", "length"])
        self.classified_reads_list = []

        self.logger.info(f"Classifier: {self.classifier_method}")
        self.logger.info(f"Query: {self.r1}")
        self.logger.info(f"Output: {self.output_path}")
        self.logger.info(f"db: {self.classifier.db_path}")

        self.finished = False
        self.deployed = False

    def run(self):
        """
        deploy classifier method. read classifier output, return only query and reference sequence id columns.
        """

        if self.classifier.method_name == "None":
            print("No classifier method selected.")
            self.logger.info("No classifier method selected.")
            return

        if not self.check_r1():
            self.collect_report()
            self.finished = self.check_classifier_output()

        else:
            if not self.check_classifier_output():
                self.classify()

            self.collect_report()
            self.finished = self.check_classifier_output()

        self.deployed = True

    def classifier_configure(self):
        """
        configure classifier.
        """

        try:
            classifier = self.available_software[self.classifier_method.name.lower()]
            os.makedirs(self.output_path, exist_ok=True)
            return classifier(
                db_path=self.classifier_method.db,
                query_path=self.r1,
                out_path=self.output_path,
                args=self.classifier_method.args,
                r2=self.r2,
                prefix=self.prefix,
                bin=self.classifier_method.bin,
                log_dir=self.log_dir,
            )

        except KeyError:
            self.logger.info(f"{self.classifier_method.name} is not available")
            return Empty_classifier(
                db_path=self.classifier_method.db,
                query_path=self.r1,
                out_path=self.output_path,
                args=self.classifier_method.args,
                r2=self.r2,
                prefix=self.prefix,
                bin=self.classifier_method.bin,
                log_dir=self.log_dir,
            )

    def check_gz_file_not_empty(self, file):
        """
        check if gzipped file has lines starting with sequence or read indicators.
        """
        cmd = "gunzip -c {} | grep '^>\|^@' | wc -l".format(file)
        number_of_sequences = int(self.cmd.run_bash_return(cmd))
        if number_of_sequences > 0:
            return True
        else:
            return False

    def check_r1(self):
        """
        check r1 is not empty
        """
        print("CHECKING R1")
        print(self.r1)
        print(os.path.isfile(self.r1))
        print(self.check_gz_file_not_empty(self.r1))

        if os.path.isfile(self.r1):
            if self.check_gz_file_not_empty(self.r1):
                return True
            else:
                self.logger.info(f"{self.r1} is empty")
                return False
        else:
            self.logger.info(f"{self.r1} does not exist")
            return False

    def classify(self):
        """
        Classify a sequencem, deploy classifier method adjusted for the type of sequence.
        """
        if self.type == ConstantsSettings.SINGLE_END:
            self.classifier.run_SE(threads=self.threads)
        else:
            self.classifier.run_PE(threads=self.threads)

    def check_classifier_output(self) -> bool:
        """
        check if classifier report exists and is not empty.
        """

        if (
            os.path.isfile(self.classifier.report_path)
            and os.path.getsize(self.classifier.report_path) > 0
        ):
            return True
        else:
            return False

    def check_classifier_output_size(self) -> bool:
        if os.path.getsize(self.classifier.report_path) == 0:
            return False
        else:
            return True

    def get_report(self) -> pd.DataFrame:
        """
        check that classifier report exists and use pandas to read full report.
        """
        if self.check_classifier_output() and self.check_classifier_output_size():
            return self.classifier.get_report()
        else:
            return None

    def get_report_simple(self) -> pd.DataFrame:
        """
        return only query and reference sequence id columns from classifier output.
        """
        return self.classifier.get_report_simple()
        # try:
        #    return self.classifier.get_report_simple()
        # except Exception as e:
        #    return pd.DataFrame(columns=["qseqid", "acc"])

    def collect_report(self) -> pd.DataFrame:
        """
        set classification_report attribute query and reference sequence id columns from classifier output.
        set classified_reads_list attribute to list of classified reads from classifier report.
        """

        self.classification_report = self.get_report_simple()

        if "qseqid" in self.classification_report.columns:

            self.classified_reads_list = list(
                set(self.classification_report.qseqid.to_list())
            )
