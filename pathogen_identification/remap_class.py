import logging
import os
import re
import shutil
from random import randint
from typing import List, Type

import numpy as np
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from result_display.plot_coverage import Bedgraph, plot_dotplot
from scipy.stats import kstest

from pathogen_detection.constants_settings import ConstantsSettings
from pathogen_detection.object_classes import (
    Read_class,
    Remap_Target,
    RunCMD,
    Software_detail,
)
from pathogen_detection.utilities import read_paf_coordinates

pd.options.mode.chained_assignment = None
np.warnings.filterwarnings("ignore")


class coverage_parse:
    def __init__(
        self,
        fastf: str,
        bedmap: str,
        output: str = "",
        Xm: int = 2,
        logging_level=logging.INFO,
    ):

        self.fastf = fastf
        self.bedmap = bedmap
        self.Xm = Xm
        self.output = output
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging_level)
        self.logger.addHandler(logging.StreamHandler())
        self.logger.propagate = False

    def fasta_segmentlength_extract(self):
        """
        calculate reference fasta total and contig lengths.
        :return: self.
        """
        ctg_lens = {}
        Tlen = 0

        with open(self.fastf) as fasta:
            for name, seq in SimpleFastaParser(fasta):
                seqLen = len(seq)
                ctg_lens[name] = seqLen
                Tlen += seqLen

        self.glen = Tlen

        self.ctgl = ctg_lens
        self.report = pd.DataFrame(
            [[x, 0, 0, 0, 0, 0, 0, 0, 0] for x in ctg_lens.keys()],
            columns=[
                "ID",
                "Hdepth",
                "HdepthR",
                "coverage",
                "nregions",
                "Rsize",
                "ngaps",
                "Gdist",
                "Gsize",
            ],
        )

        return self

    def read_bedfile(self):
        """
        read bedtools bedmap and store in self.bedm.
        :return: self.
        """
        try:
            bedm = pd.read_csv(self.bedmap, header=None, sep="\t")
            bedm.columns = ["contig", "i", "e", "x"]
            bedm["s"] = bedm.e - bedm.i
            self.bedm = bedm
        except ValueError:
            if os.path.isfile(self.bedmap):
                self.logger.error(f"bedfile {self.bedmap} is empty.")
            else:
                self.logger.error(f"bedfile {self.bedmap} is missing.")

            self.bedm = pd.DataFrame(columns=["contig", "i", "e", "x", "s"])

        return self

    def compress_bed(self, bed):
        """
        merge bedfile windows with contigous end and start points, 1 idx appart.
        :param bed: bedfile to compress.
        :return: compressed bedfile.
        """

        nrow = bed.loc[0]
        nbed = []

        for i in range(1, bed.shape[0]):
            row = bed.loc[i]

            if row.i == nrow.e or row.i == nrow.e + 1:
                nrow.e = row.e
            else:
                nbed.append(list(nrow))
                nrow = row

        nrow.e = row.e
        nbed.append(nrow)

        nbed = pd.DataFrame(nbed)
        nbed.columns = bed.columns

        nbed["s"] = bed.e - bed.i

        return nbed

    def bedstats(self, bedm, nwindows=50):
        """
        calculate average depth, % over self.Xm coverage, and gap number, size and distance.
        :param bedm: bedmap df to calculate stats on.
        :return:
        """

        depth = sum(bedm.s * bedm.x)
        overX = bedm[bedm.x >= self.Xm]
        depthR = sum(overX.s * overX.x)
        overX = sum(overX.s)
        results = [depth, depthR, overX]
        ### region operations
        bedp = bedm[bedm.x >= self.Xm].reset_index(drop=True)
        if bedp.shape[0] == 0:
            results.extend([0, 0])
        else:
            savg = sum(bedp.s) / bedp.shape[0]
            if len(bedp) > 1:
                bedp = self.compress_bed(bedp)
                pvals = []

                for ctg in bedp.contig.unique():

                    bp = bedp[bedp.contig == ctg].copy()
                    ctgsize = self.ctgl[ctg]
                    tdrange = list(np.linspace(0, ctgsize, nwindows, dtype=int))
                    td_windows = [
                        [tdrange[x], tdrange[x + 1] - 1]
                        for x in range(len(tdrange) - 1)
                    ]

                    present = []

                    for row in td_windows:
                        i, e = row
                        cmp = (bp.i <= e) & (bp.e >= i)
                        if sum(cmp) > 0:
                            present.append(i)

                    ks, pval = kstest(present, "uniform")
                    pvals.append(pval)

                pvals = sum(pvals) / len(pvals)
                savg = np.median(bedp.s)
                results.extend([bedp.shape[0], savg])

            else:
                results.extend([1, savg])
        ### gap operations.
        bedg = bedm[bedm.x < self.Xm].reset_index(drop=True)

        if bedg.shape[0] == 0:
            results.extend([0, 0, 0])
        else:
            savg = sum(bedg.s) / bedg.shape[0]
            if len(bedg) > 1:
                bedg = self.compress_bed(bedg)

                distances = []
                for ctg in bedg.contig.unique():
                    bg = bedg[bedg.contig == ctg].copy()
                    distances.append(np.sum(np.array(bg.i[1:]) - np.array(bg.e[:-1])))

                distances = np.sum(np.array(distances) / len(distances))
                savg = np.median(bedp.s)
                results.extend([bedg.shape[0], distances, savg])
            else:
                results.extend([len(bedg), 0, savg])

        return results

    def draft_report(self):
        """
        process bedstats output, join for total fasta and each sequencce individually.
        :return:
        """
        bed = self.bedm

        if bed.shape[0] == 0:
            return self
        tr = self.bedstats(bed)
        regions = bed[bed.x > self.Xm]
        regions["l"] = regions.e - regions.i
        tr[0] = tr[0] / self.glen
        sum_regions = sum(regions.l)
        if sum_regions == 0:
            tr[1] = 0
        else:
            tr[1] = tr[1] / sum(regions.l)

        tr[2] = tr[2] * 100 / self.glen
        report = []  # [["total"] + tr]

        for ctg in bed.contig.unique():
            nbed = bed[bed.contig == ctg].reset_index(drop=True)
            reggie = nbed[nbed.x > self.Xm]
            reggie["l"] = reggie.e - reggie.i
            tr = self.bedstats(nbed)
            tr[0] = tr[0] / self.ctgl[ctg]
            if sum(reggie.l) == 0:
                tr[1] = 0
            else:
                tr[1] = tr[1] / sum(reggie.l)
            tr[2] = tr[2] * 100 / self.ctgl[ctg]
            report.append([ctg] + tr)

        report = pd.DataFrame(report)

        new_columns = [
            "ID",
            "Hdepth",
            "HdepthR",
            "coverage",
            "nregions",
            "Rsize",
            "ngaps",
            "Gdist",
            "Gsize",
        ]
        report.columns = new_columns

        self.report = report

        return self

    def write(self):
        self.report.to_csv(self.output, sep="\t", index=False)


class Remapping:
    def __init__(
        self,
        r1: str,
        target: Type[Remap_Target],
        method: Type[Software_detail],
        assembly_path: str,
        type: str,
        prefix: str,
        rdir,
        threads: int = 3,
        r2: str = "",
        minimum_coverage: int = 1,
        bin: str = "",
        logging_level: int = logging.CRITICAL,
        cleanup: bool = False,
        log_dir="",
    ):
        """
        Args:
        :param r1:  path to read 1 fastq file
        :param target: target object to remap to.
        :param method: method object to use.
        :param assembly_path: path to assembly.
        :param type: type of assembly.
        :param prefix: prefix for output files.
        :param rdir: directory to write output files to.
        :param threads: number of threads to use.
        :param r2: path to read 2 fastq file.
        :param minimum_coverage: minimum coverage to use in remapping.
        :param bin: path to bin directory.
        :param logging_level: logging level to use.
        """
        self.method = method.name
        self.args = method.args
        self.rdir = rdir
        self.cleanup = cleanup

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging_level)
        self.logger.addHandler(logging.StreamHandler())
        self.logger.propagate = False
        self.logger.info("Starting remapping")

        self.target = target
        self.assembly_path = assembly_path
        self.type = type
        self.prefix = prefix
        self.threads = str(threads)
        self.r1 = r1
        self.r2 = r2
        self.minimum_coverage = minimum_coverage

        self.cmd = RunCMD(bin, logdir=log_dir, prefix=prefix, task="remapping")

        os.makedirs(self.rdir, exist_ok=True)

        self.reference_file = f"{self.rdir}/{self.prefix}_{target.acc_simple}_ref.fa"
        self.reference_fasta_index = (
            f"{self.rdir}/{self.prefix}_{target.acc_simple}_ref.fa.fai"
        )
        self.read_map_sam = f"{self.rdir}/{self.prefix}.sam"
        self.read_map_bam = f"{self.rdir}/{self.prefix}.bam"
        self.read_map_sorted_bam = (
            f"{self.rdir}/{self.prefix}.{target.acc_simple}.sorted.bam"
        )
        self.read_map_sorted_bam_index = (
            f"{self.rdir}/{self.prefix}.{target.acc_simple}.sorted.bam.bai"
        )

        self.genome_coverage = f"{self.rdir}/{self.prefix}.sorted.bedgraph"
        self.mapped_reads = f"{self.rdir}/{self.prefix}_reads_map.tsv"
        self.assembly_map_paf = f"{self.rdir}/{self.prefix}.paf"
        self.mapped_subset_r1 = f"{self.rdir}/{self.prefix}.R1.kept.fastq.gz"
        self.mapped_subset_r2 = f"{self.rdir}/{self.prefix}.R2.kept.fastq.gz"
        self.read_map_sam_rmdup = f"{self.rdir}/{self.prefix}.rmdup.sam"
        self.coverage_plot = (
            f"{self.rdir}/{self.prefix}.{target.acc_simple}.coverage.png"
        )
        self.dotplot = f"{self.rdir}/{self.prefix}.{target.acc_simple}.dotplot.png"

        self.read_map_sorted_bam_exists = os.path.isfile(self.read_map_sorted_bam)
        self.assembly_map_paf_exists = os.path.isfile(self.assembly_map_paf)
        self.mapped_subset_r1_exists = os.path.isfile(self.mapped_subset_r1)
        self.mapped_subset_r2_exists = os.path.isfile(self.mapped_subset_r2)
        self.genome_coverage_exists = os.path.isfile(self.genome_coverage)
        self.read_map_sam_rmdup_exists = False
        self.number_of_reads_mapped = 0
        self.coverage_plot_exists = os.path.isfile(self.coverage_plot)
        self.dotplot_exists = os.path.isfile(self.dotplot)

        self.logger.info("Remapping object created")
        self.logger.info("Starting remapping")
        self.logger.info("bin: %s", self.cmd.bin)
        self.logger.info("method: %s", self.method)
        self.logger.info("args: %s", self.args)
        self.logger.info("assembly: %s", self.assembly_path)

        self.output_analyser = coverage_parse(
            self.reference_file,
            self.genome_coverage,
            Xm=self.minimum_coverage,
            logging_level=self.logger.level,
        )

        self.reference_fasta_length = 0
        self.reference_fasta_string = "none"
        self.report = pd.DataFrame(
            [["None", 0, 0, 0, 0, 0, 0, 0, 0]],
            columns=[
                "ID",
                "Hdepth",
                "HdepthR",
                "coverage",
                "nregions",
                "Rsize",
                "ngaps",
                "Gdist",
                "Gsize",
            ],
        )

        self.mapped_contigs = []
        self.number_of_contigs_mapped = 0
        self.remapping_successful = False

    def cleanup_files(self):
        for file in [
            self.read_map_bam,
            self.read_map_sam,
            self.read_map_sorted_bam,
            self.read_map_sam_rmdup,
            self.assembly_map_paf,
            self.mapped_subset_r1,
            self.mapped_subset_r2,
        ]:
            if os.path.isfile(file):
                os.remove(file)

    def index_reference(self):
        cmd = "samtools faidx %s" % self.reference_file
        self.cmd.run(cmd)

    def check_remap_performed(self):
        """
        Checks if critical output files exist: bam, paf, subset reads and coverage bedgraph.
        """
        if (
            self.read_map_sorted_bam_exists
            and self.assembly_map_paf_exists
            and self.mapped_subset_r1_exists
            and self.genome_coverage_exists
        ):
            return True
        else:
            return False

    def get_input_read_number(self):
        """
        Get number of reads in current fastq files."""

        cmd = "zcat %s | wc -l" % self.r1
        rnumber = self.cmd.run_bash_return(cmd).decode("utf-8")
        rnumber = int(rnumber) // 4

        if self.type == "PE":
            cmd = "zcat %s | wc -l" % self.r2
            rnumber += int(self.cmd.run_bash_return(cmd).decode("utf-8")) // 4

        return rnumber

    def summarize(self):
        """
        Summarizes remapping results.
        generate report on read mapping to reference file.
        get number and name and length of mapped contigs.
        """
        self.report = self.calculate_mapping_statistics()
        self.get_mapped_contig_names()
        self.get_mapped_reads_number()
        self.get_reference_fasta_length()
        self.reference_fasta_string = self.get_reference_contig_name()
        self.plot_coverage()
        self.plot_dotplot_from_paf()
        self.remapping_successful = True

        if self.cleanup:
            self.cleanup_files()

    def get_reference_contig_name(self):
        cmd = [
            "grep",
            "'>'",
            self.reference_file,
            "|",
            "sed",
            "s/>//g",
        ]
        self.reference_fasta_string = self.cmd.run_bash_return(cmd)

    def run_remap(self):
        """
        Check if mapping hjas  already been perfored. if not,
        i) extract reference sequences from respective database.
        ii) map reads to reference,
        iii) map assembly to reference.

        """

        if self.check_remap_performed():
            self.logger.info("Remapping already performed")
            self.summarize()

            return self

        if len(self.target.accid_in_file) == 0:
            self.logger.info(f"No target contigs found for {self.target.accid}")
            return self

        os.makedirs(self.rdir, exist_ok=True)

        self.extract_reference_sequences()
        self.index_reference()

        if not self.check_mapping_output_exists():
            self.remap_deploy()

        if not self.check_mapping_output_exists():
            self.logger.error(
                f"Mapping output not found or unsuccesful after deploying on \
                    target(s): {self.target.accid_in_file}, file: {self.r1}, reference: {self.target.file}"
            )
            return

        self.remap_reads_post_process()
        self.assembly_to_reference_map()
        self.summarize()

    def remap_reads_post_process(self):
        """
        Process mapping output.
        1) filter read ids in bam.
        2) sort bam file.
        3) index bam file.
        4) get number of mapped reads."""

        self.filter_bamfile_read_names()
        self.sort_bam()
        self.index_sorted_bam()
        self.get_genomecoverage()
        self.get_mapped_reads_no_header()
        self.subset_mapped_reads()

    def check_mapping_output_exists(self):
        if (
            self.check_remap_status_bam()
            or self.check_remap_status_sam()
            or self.check_remap_status_paf()
        ):
            return True
        else:
            return False

    def check_remap_status_paf(self):
        if os.path.isfile(self.assembly_map_paf) and os.path.getsize(
            self.assembly_map_paf
        ):
            return True
        else:
            return False

    def check_assembly_exists(self):
        if not os.path.exists(self.assembly_path) or not os.path.getsize(
            self.assembly_path
        ):
            self.logger.error(
                f"Assembly file {self.assembly_path} does not exist or is empty"
            )
            return False
        else:
            return True

    def assembly_to_reference_map(self):

        if self.check_assembly_exists():
            self.minimap2_assembly_map()
        else:
            open(self.assembly_map_paf, "w").close()

        self.assembly_map_paf_exists = (
            os.path.isfile(self.assembly_map_paf)
            and os.path.getsize(self.assembly_map_paf) > 0
        )

    def get_mapped_contig_names(self):
        """
        Get names of contigs that are mapped to reference from paf file."""
        if not self.assembly_map_paf_exists:
            return

        self.mapped_contigs = pd.read_csv(
            self.assembly_map_paf,
            sep="\t",
            header=None,
            usecols=[5],
            names=["contig_name"],
        ).contig_name.unique()

        self.number_of_contigs_mapped = len(self.mapped_contigs)

    def extract_reference_sequences(self):
        """
        Extract reference sequences from respective database.
        """

        open(self.reference_file, "w").close()

        for accid in self.target.accid_in_file:
            cmd = f"samtools faidx {self.target.file} {accid} >> {self.reference_file}"
            self.cmd.run(cmd)

        self.reference_file_exists = (
            os.path.isfile(self.reference_file)
            and os.path.getsize(self.reference_file) > 0
        )

    def get_reference_fasta_length(self):
        """
        Get length of reference fasta file. assumes single sequence fasta.
        """

        with open(self.reference_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    continue
                self.reference_fasta_length += len(line.strip())

    def remap_bwa(self):
        """
        Remap reads to reference using bwa.
        """
        cmd = f"bwa mem -t {self.threads} {self.reference_file} {self.r1} {self.r2} > {self.read_map_sam}"
        self.cmd.run(cmd)
        self.logger.info("Finished remapping")

    def remap_snippy(self):
        """
        Remap reads to reference using snippy."""
        if self.type == "SE":
            self.remap_snippy_SE()
        elif self.type == "PE":
            self.remap_snippy_PE()
        else:
            self.logger.error(f"Remap type {self.type} not available")
            raise ValueError

    def remap_snippy_SE(self):
        """
        Remap reads to reference using snippy for single end reads."""
        cmd = [
            "snippy",
            self.args,
            "--cpus",
            self.threads,
            "--ram",
            "4",
            "--ref",
            self.reference_file,
            "--outdir",
            self.rdir,
            "--prefix",
            self.prefix,
            "--se",
            self.r1,
            "--force",
        ]
        self.cmd.run(cmd)

    def remap_snippy_PE(self):
        """
        Remap reads to reference using snippy for paired end reads."""
        cmd = [
            "snippy",
            self.args,
            "--cpus",
            self.threads,
            "--ram",
            "4",
            "--ref",
            self.reference_file,
            "--outdir",
            self.rdir,
            "--prefix",
            self.prefix,
            "--R1",
            self.r1,
            "--R2",
            self.r2,
            "--force",
        ]
        self.cmd.run(cmd)

    def remap_minimap2(self):
        """
        Remap reads to reference using minimap2. ONT data."""
        if self.type == "SE":
            self.remap_minimap2_SE()
        elif self.type == "PE":
            self.remap_minimap2_PE()

    def remap_minimap2_SE(self):

        cmd = f"minimap2 -t {self.threads} -ax map-ont {self.reference_file} {self.r1} > {self.read_map_sam}"
        self.cmd.run(cmd)

    def remap_minimap2_PE(self):

        cmd = f"minimap2 -t {self.threads} -ax map-ont {self.reference_file} {self.r1} {self.r2} > {self.read_map_sam}"
        self.cmd.run(cmd)

    def remap_minimap2_no_ref(self):
        cmd = f"minimap2 -t {self.threads} {self.r1} {self.r2} > {self.read_map_sam}"
        self.cmd.run(cmd)

    def remap_bowtie(self):
        if self.type == "SE":
            self.remap_bowtie_SE()
        elif self.type == "PE":
            self.remap_bowtie_PE()

    def remap_bowtie_SE(self):
        """
        Remap reads to reference using bowtie.
        """
        cmd = f"bowtie -p {self.threads} -x {self.reference_file} -U {self.r1} -S {self.read_map_sam}"
        self.cmd.run(cmd)

    def remap_bowtie_PE(self):
        """
        Remap reads to reference using bowtie.
        """
        cmd = f"bowtie2 -x {self.reference_file} -1 {self.r1} -2 {self.r2} -S {self.read_map_sam}"
        self.cmd.run(cmd)

    def minimap2_assembly_map(self):
        """
        Map assembly to reference using minimap2.
        """
        WHERETO = os.path.dirname(self.reference_file)
        tempfile = os.path.join(WHERETO, f"temp{randint(1,2000)}.fa")
        unzip_cmd = [
            "gunzip",
            "-c",
            self.assembly_path,
            ">",
            tempfile,
        ]
        self.cmd.run_bash(unzip_cmd)

        cmd = f"minimap2 -t {self.threads} -cx asm10 {self.assembly_path} {self.reference_file} > {self.assembly_map_paf}"
        self.cmd.run(cmd)

        os.remove(tempfile)

    def remap_deploy(self):
        """
        Configure which remapping method to use."""

        available_methods = {
            "bwa": self.remap_bwa,
            "snippy": self.remap_snippy,
            "minimap-rem": self.remap_minimap2,
            "minimap2_no_ref": self.remap_minimap2_no_ref,
            "bowtie": self.remap_bowtie,
        }

        try:
            available_methods[self.method]()
        except KeyError:
            self.logger.error(f"Remap Method {self.method} not available")
            raise KeyError

    def check_remap_status_sam(self):
        if os.path.exists(self.read_map_sam) and os.path.getsize(self.read_map_sam) > 0:
            return True
        else:
            return False

    def check_remap_status_bam(self):
        if os.path.exists(self.read_map_bam) and os.path.getsize(self.read_map_bam) > 0:
            return True
        else:
            return False

    def check_assembly_map_status(self):
        if (
            os.path.exists(self.assembly_map_paf)
            and os.path.getsize(self.assembly_map_paf) > 0
        ):
            return True
        else:
            return False

    def filter_samfile_read_names(self, same=True, output_sam=""):

        if not output_sam:
            output_sam = os.path.join(self.rdir, f"temp{randint(1,1999)}.sam")

        read_name_filter_regex = re.compile("^[A-Za-z0-9_-]*$")  # (r"@|=&$\t")

        with open(self.read_map_sam, "r") as f:
            with open(output_sam, "w") as f2:
                for line in f:
                    if line.startswith(tuple(["@HD", "@SQ", "@PG", "@RG", "@CO"])):
                        f2.write(line)

                    elif read_name_filter_regex.search(
                        line.split()[0].replace(":", "_")
                    ):
                        f2.write(line)

        if same:
            os.remove(self.read_map_sam)
            os.rename(output_sam, self.read_map_sam)

    def convert_bam_to_sam(self):
        cmd = [
            "samtools",
            "view",
            "-h",
            self.read_map_bam,
            "-o",
            self.read_map_sam,
        ]

        self.cmd.run(cmd)

    def filter_bamfile_read_names(self):
        """
        convert bam file to samfile,
        use python to filter read names in samfile.
        convert back to bam.
        """

        if not self.check_remap_status_sam():
            self.convert_bam_to_sam()

        self.filter_samfile_read_names()

        self.convert_sam_to_bam()

    def remove_duplicates_samfile(self, same=True):

        cmd = f"samtools rmdup -s {self.read_map_sam} {self.read_map_sam_rmdup}"

        self.cmd.run(cmd)

        if not os.path.isfile(self.read_map_sam_rmdup):
            self.logger.error(
                "Duplicate removal failed for file {}".format(self.read_map_sam)
            )
            return
        if same:
            os.remove(self.read_map_sam)
            os.rename(self.read_map_sam_rmdup, self.read_map_sam)

    def convert_sam_to_bam(self):

        if self.check_remap_status_sam():

            if os.path.isfile(self.read_map_bam):
                os.remove(self.read_map_bam)

            cmd = f"samtools view -bS {self.read_map_sam} > {self.read_map_bam}"
            self.cmd.run(cmd)
        else:

            self.logger.error("SAM file not found")
            raise FileNotFoundError

    def sort_bam(self):
        """
        sort bam file using samtools."""
        if self.check_remap_status_bam():
            cmd = f"samtools sort {self.read_map_bam} -o {self.read_map_sorted_bam}"
            self.cmd.run(cmd)
        else:
            self.logger.error("BAM file not found")
            raise FileNotFoundError

    def index_sorted_bam(self):
        """
        index sorted bam file using samtools."""
        cmd = f"samtools index {self.read_map_sorted_bam}"
        self.cmd.run(cmd)

    def get_genomecoverage(self):
        """
        Get genome coverage using bedtools."""
        cmd = f"bedtools genomecov -ibam {self.read_map_sorted_bam} -bga > {self.genome_coverage}"
        self.cmd.run(cmd)

    def get_mapped_reads_no_header(self):
        """
        Get number of mapped reads without header, use samtools."""
        cmd = f"samtools view -b -F 4 {self.read_map_sorted_bam} | samtools view -h | grep -v '^@' | cut -f1 | sort | uniq > {self.mapped_reads}"
        self.cmd.run(cmd)

    def get_mapped_reads_number(self):

        try:
            with open(self.mapped_reads, "r") as f:
                self.number_of_reads_mapped = len(f.readlines())
        except FileNotFoundError:
            self.number_of_reads_mapped = 0

    def subset_mapped_reads_r1(self, tempfile=""):
        """
        Subset mapped reads to R1, use seqtk."""

        cmd = f"seqtk subseq {self.r1} {tempfile} | gzip > {self.mapped_subset_r1}"
        self.cmd.run(cmd)

    def subset_mapped_reads_r2(self, tempfile=""):
        """
        Subset mapped reads to R2, use seqtk."""
        cmd = f"seqtk subseq {self.r2} {tempfile} | gzip > {self.mapped_subset_r2}"
        self.cmd.run(cmd)

    def subset_mapped_reads(self):
        """
        Subset mapped reads to R1 and R2, use seqtk."""

        tempfile = os.path.join(self.rdir, f"temp{randint(1,1999)}.rlst")
        self.cmd.run_bash(f"cat {self.mapped_reads} | cut -f1 > {tempfile}")

        if self.type == "SE":
            self.subset_mapped_reads_r1(tempfile)
        elif self.type == "PE":
            self.subset_mapped_reads_r1(tempfile)
            self.subset_mapped_reads_r2(tempfile)

        os.remove(tempfile)

    def calculate_mapping_statistics(self):
        """
        Calculate mapping statistics from bedtools bedgraph file.
        """

        self.output_analyser.fasta_segmentlength_extract()
        self.output_analyser.read_bedfile()
        self.output_analyser.draft_report()

        return self.output_analyser.report

    def plot_coverage(self):
        if os.path.getsize(self.genome_coverage):
            bedgraph = Bedgraph(self.genome_coverage)
            bedgraph.plot_coverage(self.coverage_plot, tlen=self.reference_fasta_length)

        self.coverage_plot_exists = os.path.exists(self.coverage_plot)

    def plot_dotplot_from_paf(self):
        if os.path.getsize(self.assembly_map_paf):

            df = read_paf_coordinates(self.assembly_map_paf)
            plot_dotplot(df, self.dotplot, "dotplot", xmax=self.reference_fasta_length)
            self.dotplot_exists = os.path.exists(self.dotplot)

    def move_coverage_plot(self, main_static, static_dir):
        """
        Move coverage plot to static directory."""
        new_coverage_plot = os.path.join(
            static_dir, os.path.basename(self.coverage_plot)
        )

        self.coverage_plot_exists = os.path.exists(
            os.path.join(
                ConstantsSettings.static_directory, main_static, self.coverage_plot
            )
        )

        if self.coverage_plot_exists:

            os.rename(
                self.coverage_plot,
                os.path.join(
                    ConstantsSettings.static_directory, main_static, new_coverage_plot
                ),
            )
            self.coverage_plot = new_coverage_plot

    def move_dotplot(self, main_static, static_dir):
        """
        Move dotplot to static directory."""

        new_coverage_plot = os.path.join(static_dir, os.path.basename(self.dotplot))

        self.dotplot_exists = os.path.exists(
            os.path.join(ConstantsSettings.static_directory, main_static, self.dotplot)
        )

        if self.dotplot_exists:

            os.rename(
                self.dotplot,
                os.path.join(
                    ConstantsSettings.static_directory, main_static, new_coverage_plot
                ),
            )
            self.dotplot = new_coverage_plot

    def move_igv_files(self, main_static, static_dir):
        """
        Move igv files to static directory."""
        new_bam = os.path.join(
            ConstantsSettings.static_directory,
            main_static,
            static_dir,
            os.path.basename(self.read_map_sorted_bam),
        )
        shutil.move(self.read_map_sorted_bam, new_bam)
        self.read_map_sorted_bam = new_bam

        new_bai = os.path.join(
            ConstantsSettings.static_directory,
            main_static,
            static_dir,
            os.path.basename(self.read_map_sorted_bam + ".bai"),
        )
        shutil.move(self.read_map_sorted_bam_index, new_bai)
        self.read_map_sorted_bam_index = new_bai

        new_reference_file = os.path.join(
            ConstantsSettings.static_directory,
            main_static,
            static_dir,
            os.path.basename(self.reference_file),
        )
        shutil.move(self.reference_file, new_reference_file)
        self.reference_file = new_reference_file

        new_reference_fasta_index = os.path.join(
            ConstantsSettings.static_directory,
            main_static,
            static_dir,
            os.path.basename(self.reference_fasta_index),
        )

        shutil.move(self.reference_fasta_index, new_reference_fasta_index)
        self.reference_fasta_index = new_reference_fasta_index


class Mapping_Instance:
    prefix: str
    reference: Type[Remapping] = None
    assembly: Type[Remapping] = None
    apres: bool = False
    rpres: bool = False
    success: str = "none"
    mapped: int = 0
    mapped_reads: int = 0
    original_reads: int = 0

    def __init__(
        self,
        reference: Type[Remapping],
        assembly: Type[Remapping],
        prefix,
        mapped_reads: int = 0,
        original_reads: int = 0,
    ):
        self.prefix = prefix
        self.reference = reference
        self.assembly = assembly
        self.mapped = self.reference.number_of_reads_mapped
        self.mapped_reads = mapped_reads
        self.original_reads = original_reads
        self.rpres = self.assert_reads_mapped()
        self.apres = self.assert_contigs_mapped()
        self.mapping_success = self.assert_mapping_success()
        self.classification_success = self.assert_classification_success()

        self.mapping_main_info = pd.DataFrame(
            [
                [
                    self.prefix,
                    self.reference.target.taxid,
                    self.reference.target.accid,
                    self.reference.target.description,
                    self.rpres,
                    self.apres,
                    self.mapping_success,
                    self.classification_success,
                ]
            ],
            columns=[
                "suffix",
                "taxid",
                "refseq",
                "description",
                "rclass",
                "aclass",
                "mapping_success",
                "classification_success",
            ],
        )

    def assert_reads_mapped(self):
        return self.reference.number_of_reads_mapped > 0

    def assert_contigs_mapped(self):
        return self.reference.number_of_contigs_mapped > 0

    def assert_mapping_success(self):
        ###
        if self.rpres and self.apres:
            success = "reads and contigs"
        elif self.rpres and not self.apres:
            success = "reads"
        elif self.apres and not self.rpres:
            success = "contigs"
        else:
            success = "none"

        return success

    def assert_classification_success(self):
        ###
        if self.reference.target.reads and self.reference.target.contigs:
            success = "reads and contigs"
        elif self.reference.target.reads and not self.reference.target.contigs:
            success = "reads"
        elif self.reference.target.contigs and not self.reference.target.reads:
            success = "contigs"
        else:
            success = "none"

        return success

    def generate_full_mapping_report_entry(self):

        ntax = pd.concat((self.mapping_main_info, self.reference.report), axis=1)

        def simplify_taxid(x):
            return (
                x.replace(";", "_")
                .replace(":", "_")
                .replace(".", "_")
                .replace("|", "_")
            )

        if len(self.reference.report) == 0:
            return pd.DataFrame()

        ntax["mapped"] = self.mapped
        if self.mapped_reads > 0:
            ntax["mapped_prop"] = 100 * (self.mapped / self.mapped_reads)
        else:
            self.mapped_reads = self.mapped
            ntax["mapped_prop"] = 1

        if self.original_reads > 0:
            ntax["ref_prop"] = 100 * (self.mapped / self.original_reads)

        else:
            ntax["ref_prop"] = 0

        ntax["refdb"] = self.reference.target.file
        ntax["ID"] = self.reference.target.accid
        ntax["simple_id"] = ntax["ID"].apply(simplify_taxid)
        ntax["unique_id"] = ntax["ID"].apply(simplify_taxid)

        ntax["contig_length"] = self.reference.reference_fasta_length
        ntax["contig_string"] = self.reference.reference_fasta_string

        ntax["refa_dotplot_exists"] = self.reference.dotplot_exists
        ntax["covplot_exists"] = self.reference.coverage_plot_exists
        ntax["refa_dotplot_path"] = self.reference.dotplot
        ntax["covplot_path"] = self.reference.coverage_plot
        ntax["bam_path"] = self.reference.read_map_sorted_bam
        ntax["bam_index_path"] = self.reference.read_map_sorted_bam_index
        ntax["reference_path"] = self.reference.reference_file
        ntax["reference_index_path"] = self.reference.reference_fasta_index
        ntax["reference_assembly_paf"] = self.reference.assembly_map_paf

        if self.reference.number_of_contigs_mapped > 0:
            ntax["mapped_scaffolds_path"] = self.assembly.reference_fasta_index
            ntax["mapped_scaffolds_index_path"] = self.assembly.assembly_map_paf
        else:
            ntax["mapped_scaffolds_path"] = ""
            ntax["mapped_scaffolds_index_path"] = ""

        return ntax


class Tandem_Remap:
    """
    Deep Remap class.
    deploys remap onto reference of reads and contigs.
    In a second step, if contigs mapped, maps mapped contigs to those.
    """

    reads_before_processing: int = 0
    reads_after_processing: int = 0
    logdir: str

    def __init__(
        self,
        r1,
        r2,
        remapping_method: Type[Software_detail],
        assembly_file: str,
        type: str,
        prefix,
        threads: int,
        minimum_coverage: int,
        bin: str,
        logging_level: int,
        cleanup: bool,
    ):

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging_level)
        self.logger.addHandler(logging.StreamHandler())
        self.logger.propagate = False
        self.logger.info("Reciprocal Remap started")

        self.remapping_method = remapping_method
        self.assembly_file = assembly_file
        self.type = type
        self.prefix = prefix
        self.threads = threads
        self.minimum_coverage = minimum_coverage
        self.bin = bin
        self.logging_level = logging_level
        self.cleanup = cleanup
        self.r1 = r1.current
        self.r2 = r2.current

    def reciprocal_map(self, remap_target):

        reference_remap_drone = self.reference_map(remap_target)
        assembly_map = self.assembly_map(reference_remap_drone)

        mapped_instance = Mapping_Instance(
            reference_remap_drone,
            assembly_map,
            self.prefix,
            mapped_reads=self.reads_after_processing,
            original_reads=self.reads_before_processing,
        )

        return mapped_instance

    def reference_map(self, remap_target: Type[Remap_Target]):
        rdir = os.path.join(
            self.remapping_method.dir,
            remap_target.name,
            "reference",
        )

        target_remap_drone = Remapping(
            self.r1,
            remap_target,
            self.remapping_method,
            self.assembly_file,
            self.type,
            self.prefix,
            rdir,
            self.threads,
            r2=self.r2,
            minimum_coverage=self.minimum_coverage,
            bin=self.bin,
            logging_level=self.logging_level,
            cleanup=self.cleanup,
            log_dir=self.logdir,
        )

        target_remap_drone.run_remap()

        return target_remap_drone

    def assembly_map(self, reference_remap: Type[Remapping]):

        if len(reference_remap.mapped_contigs) == 0:
            return None

        output_directory = os.path.join(
            self.remapping_method.dir,
            reference_remap.target.name,
            "assembly",
        )

        assembly_target = Remap_Target(
            "none",
            "assembly",
            "none",
            self.assembly_file,
            self.prefix,
            "description",
            reference_remap.mapped_contigs,
        )

        assembly_remap_drone = Remapping(
            reference_remap.mapped_subset_r1,
            assembly_target,
            self.remapping_method,
            self.assembly_file,
            self.type,
            self.prefix,
            output_directory,
            self.threads,
            r2=reference_remap.mapped_subset_r2,
            minimum_coverage=self.minimum_coverage,
            bin=self.bin,
            logging_level=self.logging_level,
            cleanup=self.cleanup,
            log_dir=self.logdir,
        )

        assembly_remap_drone.run_remap()

        return assembly_remap_drone


class Mapping_Manager(Tandem_Remap):
    """
    Deploy DeepRemap for a series of targets,
    generate mapping report.
    """

    mapped_instances: List[Mapping_Instance]

    def __init__(
        self,
        remap_targets: List[Remap_Target],
        r1: Type[Read_class],
        r2: Type[Read_class],
        remapping_method: Type[Software_detail],
        assembly_file: str,
        type: str,
        prefix,
        threads: int,
        minimum_coverage: int,
        bin: str,
        logging_level: int,
        cleanup: bool,
        logdir="",
    ):

        super().__init__(
            r1,
            r2,
            remapping_method,
            assembly_file,
            type,
            prefix,
            threads,
            minimum_coverage,
            bin,
            logging_level,
            cleanup,
        )

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging_level)
        self.logger.addHandler(logging.StreamHandler())
        self.logger.propagate = False
        self.logger.info("Mapping Manager started")
        self.logdir = logdir

        self.mapped_instances = []
        self.remap_targets = remap_targets
        self.reads_before_processing = r1.read_number_raw + r2.read_number_raw
        self.reads_after_processing = (
            r1.get_current_fastq_read_number() + r2.get_current_fastq_read_number()
        )
        self.max_gaps = 0
        self.max_prop = 0
        self.max_mapped = 0
        self.max_depth = 0
        self.max_depthR = 0
        self.mapped_instances = []
        self.report = pd.DataFrame(
            columns=[
                "suffix",
                "taxid",
                "refseq",
                "description",
                "rclass",
                "aclass",
                "ID",
                "Hdepth",
                "HdepthR",
                "coverage",
                "nregions",
                "Rsize",
                "ngaps",
                "Gdist",
                "Gsize",
            ],
        )

    def run_mappings(self):
        for target in self.remap_targets:
            mapped_instance = self.reciprocal_map(target)

            self.mapped_instances.append(mapped_instance)

    def move_plots_to_static(self, main_static, static_dir):
        for instance in self.mapped_instances:
            apres = instance.reference.number_of_contigs_mapped > 0
            rpres = instance.reference.number_of_reads_mapped > 0
            if rpres:
                instance.reference.move_coverage_plot(main_static, static_dir)
            if apres:
                instance.reference.move_dotplot(main_static, static_dir)

    def move_igv_to_static(self, main_static, static_dir):
        for instance in self.mapped_instances:

            if instance.reference.number_of_reads_mapped > 0:
                instance.reference.move_igv_files(main_static, static_dir)

    def merge_mapping_reports(self):

        full_report = []

        for instance in self.mapped_instances:

            ntax = instance.generate_full_mapping_report_entry()
            if len(ntax):
                full_report.append(ntax)

        if len(full_report) > 0:

            self.report = pd.concat(full_report, axis=0)
            self.clean_final_report()
            self.report = self.report.sort_values(
                ["coverage", "Hdepth"], ascending=False
            )

    def clean_final_report(self):

        self.report.ngaps = self.report.ngaps.fillna(0)

    def collect_final_report_summary_statistics(self):
        if self.report.shape[0] > 0:
            self.max_gaps = self.report.ngaps.max()
            if np.isnan(self.max_gaps):
                self.max_gaps = 0

            self.max_prop = self.report.ref_prop.max()
            if np.isnan(self.max_prop):
                self.max_prop = 0

            self.max_mapped = self.report.mapped.max()
            if np.isnan(self.max_mapped):
                self.max_mapped = 0

            self.max_depth = self.report.Hdepth.max()
            self.max_depthR = self.report.HdepthR.max()
