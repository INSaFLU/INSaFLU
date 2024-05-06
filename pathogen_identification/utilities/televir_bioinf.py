import os
import re
import shutil
import subprocess
from typing import Dict, List, Optional

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from constants.constants import Televir_Metadata_Constants
from pathogen_identification.models import RawReference
from pathogen_identification.modules.object_classes import MappingStats


class DustMasker:

    def __init__(self, dustmasker_binary_dir, temp_dir=None, id: Optional[str] = 0):
        self.dustmasker_binary = "dustmasker"
        self.temp_dir = temp_dir
        if self.temp_dir is None:
            self.temp_dir = os.path.join(
                os.path.dirname(self.dustmasker_binary), "temp"
            )
        if not os.path.exists(self.temp_dir):
            os.mkdir(self.temp_dir)

        self.bedout = os.path.join(self.temp_dir, f"dust.{id}.bed")
        self.fasta_soft_mask = os.path.join(self.temp_dir, f"dust.{id}.fasta")
        self.fasta_hard_mask = os.path.join(self.temp_dir, f"dust.{id}.hard.fasta")

    def hardmask_fasta_output(self):
        """
        replace lower case bases with N
        """
        if not os.path.exists(self.fasta_soft_mask):
            return False

        with open(self.fasta_soft_mask, "r") as f:
            with open(self.fasta_hard_mask, "w") as g:
                for record in SeqIO.parse(f, "fasta"):
                    ## replace lower case with N
                    record.seq = Seq(re.sub("[a-z]", "N", str(record.seq)))

                    g.write(record.format("fasta"))

    def mask_sequence(self, fasta_file):
        """
        Run dustmasker on a fasta file
        """

        command = f"{self.dustmasker_binary} -in {fasta_file} -outfmt fasta -out {self.fasta_soft_mask}"
        print("RUNNING DUSTMASKER", command)
        subprocess.call(command, shell=True)

        return self.fasta_soft_mask

    def mask_sequence_hard(self, fasta_file):
        """
        Run dustmasker on a fasta file
        """

        self.mask_sequence(fasta_file)
        self.hardmask_fasta_output()

        return self.fasta_hard_mask

    def run_mask_hard(self, fasta_file):
        """
        Run dustmasker on a fasta file
        """
        self.mask_sequence_hard(fasta_file)
        os.remove(self.fasta_soft_mask)

        return self.fasta_hard_mask


class TelevirBioinf:

    def __init__(self):
        self.metadata_constants = Televir_Metadata_Constants()
        self.samtools_binary = self.metadata_constants.get_software_binary("samtools")
        self.bcf_tools_binary = self.metadata_constants.get_software_binary("bcftools")

    @staticmethod
    def virosaurus_formatting(accid):
        accid_simple = accid.split(".")[0]
        return f"{accid_simple}:{accid_simple};"

    @staticmethod
    def check_sourcefile_is_virosaurus(file_path):
        if "virosaurus" in file_path:
            return True

        return False

    def process_accid(self, accid, source_file):
        if self.check_sourcefile_is_virosaurus(source_file):
            return self.virosaurus_formatting(accid)

        return accid

    def check_file_exists_not_empty(self, file_path):
        return os.path.exists(file_path) and os.path.getsize(file_path) > 100

    def extract_reference(self, source_file, accid, output_file):
        accid_code = self.process_accid(accid, source_file)
        command = (
            f'{self.samtools_binary} faidx {source_file} "{accid_code}" > {output_file}'
        )
        subprocess.call(command, shell=True)

        return self.check_file_exists_not_empty(output_file)

    def merge_references(self, files: List[str], output_file):
        command = f"cat {' '.join(files)} > {output_file}"
        subprocess.call(command, shell=True)

        return self.check_file_exists_not_empty(output_file)

    def vcf_from_bam(self, bam_list: List[str], reference: str, output_vcf: str):
        command = f"{self.metadata_constants.get_software_binary('bcftools')} mpileup -f {reference} {' '.join(bam_list)} | {self.metadata_constants.get_software_binary('bcftools')} call -mv -Ov -o {output_vcf}"
        subprocess.call(command, shell=True)
        return self.check_file_exists_not_empty(output_vcf)

    def merge_vcf_files(self, files: List[str], output_file):

        if len(files) == 1:
            if files[0].endswith(".gz"):
                command = f"zcat {files[0]} > {output_file}"
                subprocess.call(command, shell=True)
            else:
                shutil.copy(files[0], output_file)
            return self.check_file_exists_not_empty(files[0])

        command = f"{self.bcf_tools_binary} merge {' '.join(files)} > {output_file}"

        subprocess.call(command, shell=True)

        return self.check_file_exists_not_empty(output_file)

    def create_igv_report(
        self, reference: str, vcf_file: str, tracks: Dict[int, dict], output_html: str
    ):

        create_report_binary = self.metadata_constants.get_software_binary(
            "create_report"
        )
        command = f"{create_report_binary} " + f"{vcf_file}" + f" --fasta {reference} "

        names = [x.get("name", None) for x in tracks.values()]
        command += f"--samples {' '.join(names)} "

        bam_files = [x.get("bam_file", None) for x in tracks.values()]
        command += f"--tracks {' '.join(bam_files)} "

        command += f"--output {output_html}"

        print(command)
        subprocess.call(command, shell=True)

    @staticmethod
    def read_paf(paf_file):

        paf_df = pd.read_csv(paf_file, sep="\t", header=None)
        if paf_df.shape[1] > 12:
            paf_df = paf_df.iloc[:, :12]
        paf_df.columns = [
            "query_name",
            "query_length",
            "query_start",
            "query_end",
            "strand",
            "target_name",
            "target_length",
            "target_start",
            "target_end",
            "num_matches",
            "alignment_length",
            "mapping_quality",
        ]
        return paf_df

    @staticmethod
    def read_sam(sam_file):
        sam_df = pd.read_csv(sam_file, sep="\t", header=None)
        sam_df.columns = [
            "query_name",
            "flag",
            "target_name",
            "position",
            "mapq",
            "cigar",
            "rnext",
            "pnext",
            "tlen",
            "sequence",
            "quality",
        ]
        return sam_df

    @staticmethod
    def read_alignment_file(align_file):
        if align_file.endswith(".paf"):
            return TelevirBioinf.read_paf(align_file)
        elif align_file.endswith(".sam"):
            return TelevirBioinf.read_sam(align_file)
        else:
            raise ValueError(f"File {align_file} is not a valid alignment file")

    def get_samtools_stats(self, alignment_file) -> str:
        """
        Get the stats from the alignment file, sam or bam"""
        samtools_binary = self.metadata_constants.get_software_binary("samtools")

        from utils.utils import Utils

        utils = Utils()

        tmp_file = utils.get_temp_file("mapping_error_rate", "txt")

        if tmp_file is None:
            raise ValueError("Could not create temp file")

        cmd = [
            samtools_binary,
            "stats",
            alignment_file,
            "|",
            "grep ^SN",
            "|",
            "cut -f 2-",
            ">",
            tmp_file,
        ]

        subprocess.call(" ".join(cmd), shell=True)

        return tmp_file

    def extract_mapping_stats(self, alignment_file) -> MappingStats:
        """
        read stats as pd data frame, pass to class"""

        samtools_stats = self.get_samtools_stats(alignment_file)

        stats_df = pd.read_csv(
            samtools_stats, sep="\t", header=None, index_col=0
        ).rename(columns={0: "stat", 1: "value", 2: "comment"})
        error_rate = stats_df.loc["error rate:", "value"]
        quality_avg = stats_df.loc["average quality:", "value"]

        return MappingStats(error_rate, quality_avg)

    def check_file_exists(self, file_path):

        return os.path.exists(file_path)

    def check_file_exists_not_empty(self, file_path):

        return os.path.exists(file_path) and os.path.getsize(file_path) > 100

    def alignment_agg_by_target(self, alignment_file):
        """
        Aggregate the alignment file by target name
        counts: number of reads aligned to the target
        tlen: total length of the reads aligned to the target

        """
        if not self.check_file_exists_not_empty(alignment_file):
            return pd.DataFrame(columns=["target_name", "query_name"])

        alignment_df = self.read_alignment_file(alignment_file)

        agg_df = alignment_df.groupby("target_name").agg({"query_name": "count"})

        return agg_df
