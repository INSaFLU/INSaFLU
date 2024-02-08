import logging
import os
import re
import shutil
from random import randint
from typing import List, Type

import numpy as np
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from scipy.stats import kstest

from constants.software_names import SoftwareNames
from pathogen_identification.constants_settings import ConstantsSettings as CS
from pathogen_identification.modules.object_classes import (Bedgraph,
                                                            MappingStats,
                                                            Read_class,
                                                            Remap_Target,
                                                            RunCMD,
                                                            SoftwareDetail,
                                                            SoftwareRemap)
from pathogen_identification.utilities.televir_parameters import RemapParams
from pathogen_identification.utilities.utilities_general import (
    plot_dotplot, read_paf_coordinates)

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
        if self.logger.hasHandlers():
            self.logger.handlers.clear()
        self.logger.propagate = False

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
            [[x, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0] for x in ctg_lens.keys()],
            columns=[
                "ID",
                "Hdepth",
                "HdepthR",
                "coverage",
                "nregions",
                "Rsize",
                "windows_covered",
                "pval_uniform",
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

    def compress_bed(self, bed):
        """
        merge bedfile windows with contigous end and start points, 1 idx appart.
        :param bed: bedfile to compress.
        :return: compressed bedfile.
        """

        nrow = bed.loc[0]
        nbed = []

        for i in range(bed.shape[0]):
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

    def calculate_windows(self, ctgsize, minw=3, maxw=10, wsize=2000):
        """calculate number of windows to use for coverage calculation."""
        ratio = ctgsize / wsize

        if ratio < minw:
            return minw
        elif ratio > maxw:
            return maxw
        else:
            return int(ratio)

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
        windows_covered = "NA"

        if bedp.shape[0] == 0:
            windows_covered = f"0/{nwindows}"
            results.extend([0, 0, windows_covered, 1])

        else:
            savg = sum(bedp.s) / bedp.shape[0]

            if len(bedp) > 0:
                bedp = self.compress_bed(bedp)
                pvals = []

                for ctg in bedp.contig.unique():
                    bp = bedp[bedp.contig == ctg].copy()
                    ctgsize = self.ctgl[ctg]
                    nwindows = self.calculate_windows(ctgsize)

                    tdrange = list(np.linspace(0, ctgsize, nwindows + 1, dtype=int))
                    td_windows = [
                        [
                            tdrange[x],
                            tdrange[x + 1] - 1 if tdrange[x + 1] < ctgsize else ctgsize,
                        ]
                        for x in range(len(tdrange) - 1)
                    ]

                    present = []

                    for row in td_windows:
                        i, e = row
                        cmp = (bp.i <= e) & (bp.e >= i)
                        if sum(cmp) > 0:
                            present.append(i)

                    windows_covered = f"{len(present)}/{nwindows}"

                    ks, pval = kstest(present, "uniform")
                    pvals.append(pval)

                pvals = sum(pvals) / len(pvals)
                savg = bedp.s.median()

                results.extend([bedp.shape[0], savg, windows_covered, pvals])

            else:
                results.extend([1, savg, "NA", 1])

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

                if len(distances) > 0:
                    distances = np.sum(np.array(distances) / len(distances))
                else:
                    distances = 0
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
            tr[1] = tr[1] / sum_regions

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
            tr[2] = tr[2] * 100 / self.ctgl[ctg]  # % of genome covered
            report.append([ctg] + tr)

        report = pd.DataFrame(report)

        new_columns = [
            "ID",
            "Hdepth",
            "HdepthR",
            "coverage",
            "nregions",
            "Rsize",
            "windows_covered",
            "pval_uniform",
            "ngaps",
            "Gdist",
            "Gsize",
        ]

        report.columns = new_columns

        self.report = report

    def write(self):
        self.report.to_csv(self.output, sep="\t", index=False)


class RemapMethod_init:
    def __init__(
        self,
        method: SoftwareDetail,
        r1,
        r2,
        args,
        type,
        prefix,
        reference,
        outdir,
        threads,
        force,
        logdir,
    ):
        self.r1 = r1
        self.r2 = r2
        self.args: str = args
        self.type: str = type
        self.prefix: str = prefix
        self.reference: str = reference
        self.outdir: str = outdir
        self.threads: str = threads
        self.force = force

        self.outbam = os.path.join(outdir, prefix + ".bam")
        self.outsam = os.path.join(outdir, prefix + ".sam")

        self.cmd = RunCMD(
            method.bin, logdir=logdir, prefix=prefix, task="remap_software"
        )


class Remap_Snippy(RemapMethod_init):
    def __init__(
        self,
        method: SoftwareDetail,
        r1,
        r2,
        args,
        type,
        prefix,
        reference,
        outdir,
        threads,
        force,
        logdir,
    ):
        super().__init__(
            method,
            r1,
            r2,
            args,
            type,
            prefix,
            reference,
            outdir,
            threads,
            force,
            logdir,
        )

        self.rundir = os.path.join(outdir, "snippy")
        os.makedirs(self.rundir, exist_ok=True)
        self.runbam = os.path.join(self.rundir, "snps.bam")

    def export(self):
        if os.path.exists(self.runbam):
            shutil.copy(self.runbam, self.outbam)

    def remap(self):
        """
        Remap reads to reference using snippy."""
        if self.type == CS.SINGLE_END:
            self.remap_SE()
        elif self.type == CS.PAIR_END:
            self.remap_PE()
        else:
            raise ValueError

    def remap_SE(self):
        """
        Remap reads to reference using snippy for single end reads."""
        cmd = [
            "snippy",
            self.args,
            "--cpus",
            self.threads,
            "--ref",
            self.reference,
            "--outdir",
            self.rundir,
            "--se",
            self.r1,
            "--force",
        ]
        self.cmd.run(cmd)
        self.export()

    def remap_PE(self):
        """
        Remap reads to reference using snippy for paired end reads."""
        cmd = [
            "snippy",
            self.args,
            "--cpus",
            self.threads,
            "--ref",
            self.reference,
            "--outdir",
            self.rundir,
            "--R1",
            self.r1,
            "--R2",
            self.r2,
            "--force",
        ]
        self.cmd.run(cmd)
        self.export()


class Remap_Bwa(RemapMethod_init):
    def remap(self):
        """
        Remap reads to reference using bwa."""
        if self.type == CS.SINGLE_END:
            self.remap_SE()
        elif self.type == CS.PAIR_END:
            self.remap_PE()
        else:
            raise ValueError

    def remap_SE(self):
        """

        Remap reads to reference using bwa for single end reads."""
        temp_sam = os.path.join(self.outdir, self.prefix + ".sam")
        cmd_01 = [
            "bwa",
            "mem",
            self.args,
            "-t",
            self.threads,
            self.reference,
            self.r1,
            ">",
            temp_sam,
        ]
        cmd_samtools = [
            "samtools",
            "view",
            "-b",
            "-o",
            self.outbam,
            temp_sam,
        ]
        self.cmd.run_script(cmd_01)
        self.cmd.run_script(cmd_samtools)

    def remap_PE(self):
        """
        Remap reads to reference using bwa for paired end reads."""
        temp_sam = os.path.join(self.outdir, self.prefix + ".sam")
        cmd = [
            "bwa",
            "mem",
            self.args,
            "-t",
            self.threads,
            self.reference,
            self.r1,
            self.r2,
            ">",
            temp_sam,
        ]
        cmd_samtools = [
            "samtools",
            "view",
            "-b",
            "-o",
            self.outbam,
            temp_sam,
        ]
        self.cmd.run_script(cmd)
        self.cmd.run_script(cmd_samtools)


class Remap_Minimap2(RemapMethod_init):
    def remap(self):
        """
        Remap reads to reference using minimap2."""
        if self.type == CS.SINGLE_END:
            self.remap_SE()
        elif self.type == CS.PAIR_END:
            self.remap_PE()
        else:
            raise ValueError

    def filter_secondary(self):
        """
        Filter secondary alignments from bam file."""
        cmd = [
            "samtools",
            "view",
            "-F0x900",
            self.outsam,
            ">",
            self.outsam + ".filtered",
        ]

        try:
            self.cmd.run(cmd)
            os.remove(self.outsam)
            shutil.copy(self.outsam + ".filtered", self.outsam)
        except:
            pass

    def remap_SE(self):
        """
        Remap reads to reference using minimap2 for single end reads."""
        cmd = [
            "minimap2",
            self.args,
            "-t",
            self.threads,
            "-ax",
            "map-ont",
            self.reference,
            self.r1,
            ">",
            self.outsam,
        ]
        self.cmd.run(cmd)
        # self.filter_secondary()

    def remap_PE(self):
        """
        Remap reads to reference using minimap2 for paired end reads."""
        cmd = [
            "minimap2",
            self.args,
            "-t",
            self.threads,
            "-ax",
            "map-ont",
            self.reference,
            self.r1,
            self.r2,
            ">",
            self.outsam,
        ]
        self.cmd.run(cmd)
        # self.filter_secondary()


class Remap_Bowtie2(RemapMethod_init):
    modes = ["--end-to-end", "--local"]
    preset_options = ["--very-fast", "--fast", "--sensitive", "--very-sensitive"]

    def remap(self):
        """
        Remap reads to reference using bowtie2."""
        if self.type == CS.SINGLE_END:
            self.remap_SE()
        elif self.type == CS.PAIR_END:
            self.remap_PE()
        else:
            raise ValueError

    def process_arguments(self):
        """
        Process arguments to remove preset and mode option flags"""
        self.args = self.args.replace("[preset]", "").replace("[mode]", "")

    def index_reference(self):
        """
        Index reference using bowtie2."""
        index_name = os.path.splitext(self.reference)[0] + "_index"

        cmd = ["bowtie2-build", self.reference, index_name]

        self.cmd.run(cmd)

        return index_name

    def remap_SE(self):
        """
        Remap reads to reference using bowtie2 for single end reads."""

        index_name = self.index_reference()
        self.process_arguments()

        cmd = [
            "bowtie2",
            self.args,
            "-p",
            self.threads,
            "-x",
            index_name,
            "-U",
            self.r1,
            "-S",
            self.outsam,
        ]
        self.cmd.run(cmd)

    def remap_PE(self):
        """
        Remap reads to reference using bowtie2 for paired end reads."""

        index_name = self.index_reference()
        self.process_arguments()

        cmd = [
            "bowtie2",
            self.args,
            "-p",
            self.threads,
            "-x",
            index_name,
            "-1",
            self.r1,
            "-2",
            self.r2,
            "-S",
            self.outsam,
        ]
        self.cmd.run(cmd)


class Remapping:
    remap_engine = RemapMethod_init

    def __init__(
        self,
        r1: str,
        target: Remap_Target,
        methods: SoftwareRemap,
        assembly_path: str,
        type: str,
        prefix: str,
        rdir: str,
        remap_params: RemapParams,
        threads: int = 3,
        r2: str = "",
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
        remap_method = methods.remap_software
        self.remap_filters = methods.remap_filters
        self.method = remap_method.name.split("_")[0]
        self.method_object = remap_method
        self.args = remap_method.args
        self.rdir = rdir
        self.cleanup = cleanup
        self.remap_params = remap_params

        self.logger = logging.getLogger(f"{__name__}_{prefix}")
        if self.logger.hasHandlers():
            self.logger.handlers.clear()
        self.logger.propagate = False

        self.logger.setLevel(logging.ERROR)
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
        # self.minimum_coverage = minimum_coverage
        self.logdir = log_dir

        self.cmd = RunCMD(bin, logdir=log_dir, prefix=prefix, task="remapping_instance")

        os.makedirs(self.rdir, exist_ok=True)

        self.reference_file = f"{self.rdir}/{self.prefix}_{target.acc_simple}_ref.fa"
        self.reference_fasta_index = (
            f"{self.rdir}/{self.prefix}_{target.acc_simple}_ref.fa.fai"
        )
        self.read_map_sam = f"{self.rdir}/{self.prefix}.sam"
        self.read_map_bam = f"{self.rdir}/{self.prefix}.bam"
        self.current = self.read_map_bam
        self.read_map_filtered_bam = f"{self.rdir}/{self.prefix}.filtered.bam"
        self.read_map_sorted_bam = (
            f"{self.rdir}/{self.prefix}.{target.acc_simple}.sorted.bam"
        )
        self.read_map_sorted_bam_index = (
            f"{self.rdir}/{self.prefix}.{target.acc_simple}.sorted.bam.bai"
        )
        self.read_map_sorted_bam_stats = (
            f"{self.rdir}/{self.prefix}.{target.acc_simple}.sorted.bam.stats"
        )

        self.genome_coverage = f"{self.rdir}/{self.prefix}.sorted.bedgraph"
        self.mapped_reads = f"{self.rdir}/{self.prefix}_reads_map.tsv"
        self.mapped_subset_r1_fasta = (
            f"{self.rdir}/{self.prefix}_{target.acc_simple}_reads_map_subset_r1.fasta"
        )
        self.mapped_subset_r2_fasta = (
            f"{self.rdir}/{self.prefix}_{target.acc_simple}_reads_map_subset_r2.fasta"
        )

        self.assembly_map_paf = f"{self.rdir}/{self.prefix}_{target.acc_simple}.paf"
        self.mapped_subset_r1 = f"{self.rdir}/{self.prefix}.R1.kept.fastq.gz"
        self.mapped_subset_r2 = f"{self.rdir}/{self.prefix}.R2.kept.fastq.gz"
        self.read_map_sam_rmdup = f"{self.rdir}/{self.prefix}.rmdup.sam"
        self.mapped_contigs_fasta = (
            f"{self.rdir}/{self.prefix}_{target.acc_simple}_mapped_contigs.fasta"
        )
        self.mapped_contigs_fasta_index = self.mapped_contigs_fasta + ".fai"

        self.vcf = f"{self.rdir}/{self.prefix}.vcf"

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
            Xm=remap_params.min_coverage,
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

    @staticmethod
    def relocate_file(filepath, destination):
        """move file in filepath to destination"""
        subdirectory = os.path.join(destination, "remapping")

        os.makedirs(subdirectory, exist_ok=True)
        final_file = os.path.join(subdirectory, os.path.basename(filepath))

        if os.path.exists(filepath) and final_file != filepath:
            if os.path.exists(final_file):
                os.remove(final_file)

            shutil.move(filepath, final_file)

        return final_file

    def relocate_reference_fasta(self, destination):
        """move reference fasta file to destination"""

        self.reference_file = self.relocate_file(self.reference_file, destination)
        self.reference_fasta_index = self.relocate_file(
            self.reference_fasta_index, destination
        )

    def relocate_mapping_files(self, destination):
        """move all files generated by remapping to destination
        possibly for future use:
        self.mapped_subset_r1 = self.relocate_file(self.mapped_subset_r1, destination)
        self.mapped_subset_r2 = self.relocate_file(self.mapped_subset_r2, destination)
        """

        self.mapped_subset_r1_fasta = self.relocate_file(
            self.mapped_subset_r1_fasta, destination
        )

        self.mapped_subset_r2_fasta = self.relocate_file(
            self.mapped_subset_r2_fasta, destination
        )

        self.reference_file = self.relocate_file(self.reference_file, destination)
        self.reference_fasta_index = self.relocate_file(
            self.reference_fasta_index, destination
        )

        self.read_map_sorted_bam = self.relocate_file(
            self.read_map_sorted_bam, destination
        )
        self.read_map_sorted_bam_index = self.relocate_file(
            self.read_map_sorted_bam_index, destination
        )
        self.assembly_map_paf = self.relocate_file(self.assembly_map_paf, destination)

        self.mapped_contigs_fasta = self.relocate_file(
            self.mapped_contigs_fasta, destination
        )
        self.mapped_contigs_fasta_index = self.relocate_file(
            self.mapped_contigs_fasta_index, destination
        )
        self.vcf = self.relocate_file(self.vcf, destination)

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
        rnumber = self.cmd.run_bash_return(cmd)
        rnumber = int(rnumber) // 4

        if self.type == CS.PAIR_END:
            cmd = "zcat %s | wc -l" % self.r2
            rnumber += int(self.cmd.run_bash_return(cmd)) // 4

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

        self.plot_coverage()
        self.plot_dotplot_from_paf()
        self.remapping_successful = True
        self.generate_mapped_contigs_fasta()
        self.index_mapped_contigs_fasta()

    def get_reference_contig_name(self):
        contig_names = []

        with open(self.reference_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    contig_names.append(line.strip().replace(">", ""))

        return "; ".join(contig_names)

    def sanitize_reference_fasta_contig_name(self):
        """
        Some fasta files have spaces in the contig names. This is not allowed for samtools.
        """
        with open(self.reference_file, "r") as f:
            lines = f.readlines()

        with open(self.reference_file, "w") as f:
            for line in lines:
                if line.startswith(">"):
                    line = line.replace(";", "_").replace(":", "_")
                f.write(line)

    def retrieve_reference(self):
        self.extract_reference_sequences()
        self.get_reference_fasta_length()
        self.reference_fasta_string = self.get_reference_contig_name()
        self.sanitize_reference_fasta_contig_name()

    def run_remap(self):
        """
        Check if mapping hjas  already been perfored. if not,
        i) extract reference sequences from respective database.
        ii) map reads to reference,
        iii) map assembly to reference.

        """

        if self.check_remap_performed():
            self.logger.info("Remapping already performed")
            self.reference_file_exists = True

            if self.reference_file_exists:
                self.retrieve_reference()
                self.index_reference()
                self.summarize()

            return self

        if len(self.target.accid_in_file) == 0:
            self.logger.info(f"No target contigs found for {self.target.accid}")
            return self

        os.makedirs(self.rdir, exist_ok=True)

        self.retrieve_reference()

        if not self.reference_file_exists:
            self.logger.info(f"No reference file found for {self.target.accid}")
            return self

        self.index_reference()

        if not self.check_mapping_output_exists():
            self.remap_deploy()

        if self.check_mapping_output_exists():
            self.remap_reads_post_process()
            self.assembly_to_reference_map()
            self.summarize()

        else:
            # self.logger.error(
            #    f"Mapping output not found or unsuccesful after deploying on \
            #        target(s): {self.target.accid_in_file}, file: {self.r1}, reference: {self.target.file}"
            # )
            return

    def remap_reads_post_process(self):
        """
        Process mapping output.
        1) filter read ids in bam.
        2) sort bam file.
        3) index bam file.
        4) get number of mapped reads."""

        self.process_bam()
        self.generate_vcf()
        self.get_genomecoverage()
        self.get_mapped_reads_no_header()
        self.filter_sam_file_mapped()
        self.subset_mapped_reads()
        self.mapped_reads_to_fasta()

    def process_bam(self):
        self.filter_bamfile_read_names()
        self.filter_bamfile()
        self.sort_bam()
        self.index_sorted_bam()

    def filter_bamfile(self):
        self.read_map_filtered_bam = self.read_map_bam

        for filter in self.remap_filters.software_list:
            if filter.name == SoftwareNames.SOFTWARE_BAMUTIL_name:
                self.filter_mapping_bamutil(filter)

            if filter.name == SoftwareNames.SOFTWARE_MSAMTOOLS_name:
                self.filter_mapping_msamtools(filter)

        self.filter_bam_unmapped()

    def filter_bam_unmapped(self):
        """
        filter reads marked as unmapped"""

        bam_path = self.read_map_filtered_bam

        filtered_bam_path = os.path.splitext(bam_path)[0] + ".filtered.bam"

        bash_cmd = f"samtools view -b -F 4 {bam_path} > {filtered_bam_path}"

        self.cmd.run_script_software(bash_cmd)

        if (
            os.path.isfile(filtered_bam_path)
            and os.path.getsize(filtered_bam_path) > 100
        ):
            os.remove(bam_path)
            shutil.move(filtered_bam_path, bam_path)

    def filter_mapping_bamutil(self, software: SoftwareDetail):
        """
        filter bam file by mapping quality.
        """

        temp_file = (
            os.path.splitext(self.read_map_bam)[0] + f"{np.random.randint(1000000)}.bam"
        )

        cmd = [
            "bam",
            "filter",
            "--in",
            self.read_map_filtered_bam,
            "--refFile",
            self.reference_file,
            software.args,
            "--out",
            temp_file,
            "--noPhoneHome",
        ]

        try:
            self.cmd.run(cmd)

        except Exception as e:
            self.logger.error("Bam filtering failed.")
            self.logger.error(e)
            if os.path.isfile(temp_file):
                os.remove(temp_file)
            return

        if os.path.isfile(temp_file) and os.path.getsize(temp_file) > 100:
            os.remove(self.read_map_filtered_bam)
            shutil.move(temp_file, self.read_map_filtered_bam)

        else:
            self.logger.error("Bam filtering failed, file missing.")

        return

    def filter_mapping_msamtools(self, software: SoftwareDetail):
        """
        filter bam file by mapping quality.
        """

        temp_file = (
            os.path.splitext(self.read_map_bam)[0] + f"{np.random.randint(1000000)}.bam"
        )

        cmd = [
            "msamtools",
            "filter",
            software.args,
            self.read_map_filtered_bam,
            ">",
            temp_file,
        ]

        try:
            self.cmd.run_script_software(cmd)

        except Exception as e:
            self.logger.error("Bam filtering failed.")
            self.logger.error(e)
            if os.path.isfile(temp_file):
                os.remove(temp_file)
            return

        if os.path.isfile(temp_file) and os.path.getsize(temp_file) > 100:
            os.remove(self.read_map_filtered_bam)
            shutil.move(temp_file, self.read_map_filtered_bam)

        return

    def generate_vcf(self):
        """
        Generate vcf file from bam file.
        """

        if not self.check_remap_status_bam():
            self.logger.error("Bam file not found.")
            return

        if self.check_vcf_exists():
            self.logger.info("Vcf file already exists.")
            return

        cmd = [
            "samtools",
            "mpileup",
            "-f",
            self.reference_file,
            self.read_map_sorted_bam,
            ">",
            self.vcf,
        ]

        try:
            self.cmd.run(cmd)

        except:
            self.logger.error("Vcf generation failed.")
            return

    def check_vcf_exists(self):
        """
        Check if vcf file exists.
        """

        if os.path.isfile(self.vcf):
            return True
        else:
            return False

    def filter_sam_file_mapped(self):
        """
        Filter sam file for mapped reads.
        """
        cmd = "samtools view -F 4 %s > %s" % (
            self.read_map_sorted_bam,
            self.read_map_sam_rmdup,
        )
        self.cmd.run_script_software(cmd)

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
            and os.path.getsize(self.assembly_map_paf) > 100
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
            usecols=[0],
            names=["contig_name"],
        ).contig_name.unique()

        self.number_of_contigs_mapped = len(self.mapped_contigs)

    def generate_mapped_contigs_fasta(self):
        """
        Generate fasta file of contigs that are mapped to reference."""
        if not self.assembly_map_paf_exists or self.number_of_contigs_mapped == 0:
            return

        for accid in self.mapped_contigs:
            cmd = f"samtools faidx {self.assembly_path} '{accid}' >> {self.mapped_contigs_fasta}"
            self.cmd.run(cmd)

    def index_mapped_contigs_fasta(self):
        """
        Index fasta file of contigs that are mapped to reference."""
        if not self.assembly_map_paf_exists or self.number_of_contigs_mapped == 0:
            return

        cmd = f"samtools faidx {self.mapped_contigs_fasta}"
        self.cmd.run(cmd)

    def extract_reference_sequences(self):
        """
        Extract reference sequences from respective database.
        """

        open(self.reference_file, "w").close()

        for accid in self.target.accid_in_file:
            cmd = (
                f"samtools faidx {self.target.file} '{accid}' >> {self.reference_file}"
            )
            self.cmd.run(cmd)

        self.reference_file_exists = (
            os.path.isfile(self.reference_file)
            and os.path.getsize(self.reference_file) > 100
        )

        if not self.reference_file_exists:
            self.logger.error(
                f"Reference file {self.reference_file} does not exist or is empty"
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

        cmd = f"minimap2 -t {self.threads} -cx asm10 {self.reference_file} {self.assembly_path} > {self.assembly_map_paf}"
        self.cmd.run(cmd)

        os.remove(tempfile)

    def remap_deploy(self):
        """
        Configure which remapping method to use."""

        available_methods = {
            "bwa": Remap_Bwa,
            "snippy": Remap_Snippy,
            "snippy_pi": Remap_Snippy,
            "minimap2": Remap_Minimap2,
            "bowtie2": Remap_Bowtie2,
            "bowtie2_remap": Remap_Bowtie2,
        }

        if self.method in available_methods:
            self.logger.info(f"Remapping with {self.method}")
            self.logger.info(f"reference {self.reference_file}")
            self.remap_engine = available_methods[self.method](
                self.method_object,
                self.r1,
                self.r2,
                self.args,
                self.type,
                self.prefix,
                self.reference_file,
                self.rdir,
                self.threads,
                False,
                self.logdir,
            )
            try:
                self.remap_engine.remap()
            except Exception as e:
                self.logger.error(f"Remapping failed with {e}")
                raise e

        else:
            self.logger.error(f"Remapping method {self.method} not available")
            raise ValueError

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

        read_name_filter_regex = re.compile("^[A-Za-z0-9_.-]*$")  # (r"@|=&$\t")

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
            cmd = f"samtools sort {self.read_map_filtered_bam} -o {self.read_map_sorted_bam}"
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

        cmd2 = f"samtools view -F 0x4 {self.read_map_sorted_bam} | cut -f 1 | sort | uniq > {self.mapped_reads}"
        self.cmd.run_script_software(cmd2)

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
        self.cmd.run_script_software(cmd)

    def subset_mapped_reads_r2(self, tempfile=""):
        """
        Subset mapped reads to R2, use seqtk."""
        cmd = f"seqtk subseq {self.r2} {tempfile} | gzip > {self.mapped_subset_r2}"
        self.cmd.run_script_software(cmd)

    def subset_mapped_reads(self):
        """
        Subset mapped reads to R1 and R2, use seqtk."""

        tempfile = os.path.join(self.rdir, f"temp{randint(1,1999)}.rlst")
        self.cmd.run_bash(f"cat {self.mapped_reads} | cut -f1 > {tempfile}")

        if self.type == CS.SINGLE_END:
            self.subset_mapped_reads_r1(tempfile)
        elif self.type == CS.PAIR_END:
            self.subset_mapped_reads_r1(tempfile)
            self.subset_mapped_reads_r2(tempfile)

        os.remove(tempfile)

    def mapped_reads_to_fasta(self):
        """convert fastq subsets to fasta"""
        if self.type == CS.SINGLE_END:
            cmd = (
                f"seqtk seq -a {self.mapped_subset_r1} > {self.mapped_subset_r1_fasta}"
            )
            self.cmd.run_script_software(cmd)
        elif self.type == CS.PAIR_END:
            cmd = (
                f"seqtk seq -a {self.mapped_subset_r1} > {self.mapped_subset_r1_fasta}"
            )
            self.cmd.run_script_software(cmd)
            cmd2 = (
                f"seqtk seq -a {self.mapped_subset_r2} > {self.mapped_subset_r2_fasta}"
            )
            self.cmd.run_script_software(cmd2)

    def calculate_mapping_statistics(self):
        """
        Calculate mapping statistics from bedtools bedgraph file.
        """

        self.output_analyser.fasta_segmentlength_extract()
        self.output_analyser.read_bedfile()
        self.output_analyser.draft_report()

        bedfile_report = self.output_analyser.report
        self.generate_samtools_stats()
        mapping_stats = self.extract_mapping_stats()
        bedfile_report["error_rate"] = mapping_stats.error_rate
        bedfile_report["quality_avg"] = mapping_stats.quality_avg

        return bedfile_report

    def generate_samtools_stats(self):
        """
        Generate samtools stats file.
        """
        cmd = [
            "samtools",
            "stats",
            self.read_map_sorted_bam,
            "|",
            "grep ^SN",
            "|",
            "cut -f 2-",
            ">",
            self.read_map_sorted_bam_stats,
        ]

        self.cmd.run_script_software(cmd)

    def extract_mapping_stats(self) -> MappingStats:
        """
        read stats as pd data frame, pass to class"""

        stats_df = pd.read_csv(
            self.read_map_sorted_bam_stats, sep="\t", header=None, index_col=0
        ).rename(columns={0: "stat", 1: "value", 2: "comment"})
        error_rate = stats_df.loc["error rate:", "value"]
        quality_avg = stats_df.loc["average quality:", "value"]

        return MappingStats(error_rate, quality_avg)

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

    def move_coverage_plot(self, static_dir_plots):
        """
        Move coverage plot to static directory."""
        new_coverage_plot = os.path.join(
            static_dir_plots, os.path.basename(self.coverage_plot)
        )

        self.coverage_plot_exists = os.path.exists(self.coverage_plot)

        self.full_path_coverage_plot = os.path.join(
            CS.static_directory, new_coverage_plot
        )

        if self.coverage_plot_exists:
            shutil.move(self.coverage_plot, self.full_path_coverage_plot)
            self.coverage_plot = new_coverage_plot

    def move_dotplot(self, static_dir_plots):
        """Move dotplot to static directory."""

        new_dotplot = os.path.join(static_dir_plots, os.path.basename(self.dotplot))

        self.full_path_dotplot = os.path.join(CS.static_directory, new_dotplot)

        self.dotplot_exists = os.path.exists(self.dotplot)

        if os.path.exists(self.dotplot):
            shutil.move(self.dotplot, self.full_path_dotplot)
            self.dotplot = new_dotplot

    def move_igv_files(self, static_dir):
        """
        Move igv files to static directory."""

        new_bam = os.path.join(
            static_dir,
            os.path.basename(self.read_map_sorted_bam),
        )
        shutil.move(self.read_map_sorted_bam, new_bam)
        self.read_map_sorted_bam = new_bam

        new_bai = os.path.join(
            static_dir,
            os.path.basename(self.read_map_sorted_bam + ".bai"),
        )
        shutil.move(self.read_map_sorted_bam_index, new_bai)
        self.read_map_sorted_bam_index = new_bai

        new_reference_file = os.path.join(
            static_dir,
            os.path.basename(self.reference_file),
        )
        shutil.move(self.reference_file, new_reference_file)
        self.reference_file = new_reference_file

        new_reference_fasta_index = os.path.join(
            static_dir,
            os.path.basename(self.reference_fasta_index),
        )

        shutil.move(self.reference_fasta_index, new_reference_fasta_index)
        self.reference_fasta_index = new_reference_fasta_index


class Mapping_Instance:
    prefix: str
    reference: Remapping = None
    assembly: Remapping = None
    apres: bool = False
    rpres: bool = False
    success: str = "none"
    mapped: int = 0
    mapped_reads: int = 0
    original_reads: int = 0

    def __init__(
        self,
        reference: Remapping,
        assembly: Remapping,
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
        self.produce_mapping_report()

    def produce_mapping_report(self):
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

    def export_mapping_files(self, destination):
        """move files to media directory"""

        if self.classification_success is not "none":
            # self.reference.move_igv_files(destination)
            self.reference.relocate_mapping_files(destination)

            if self.assembly:
                self.assembly.relocate_mapping_files(destination)

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
            ntax["mapped_scaffolds_path"] = self.reference.mapped_contigs_fasta
            ntax[
                "mapped_scaffolds_index_path"
            ] = self.reference.mapped_contigs_fasta_index
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
        remapping_methods: SoftwareRemap,
        remapping_params: RemapParams,
        assembly_file: str,
        type: str,
        prefix,
        threads: int,
        bin: str,
        logging_level: int,
        cleanup: bool,
    ):
        self.logger = logging.getLogger(__name__)
        if self.logger.hasHandlers():
            self.logger.handlers.clear()
        self.logger.propagate = False
        self.logger.setLevel(logging_level)
        self.logger.addHandler(logging.StreamHandler())
        self.logger.propagate = False
        self.logger.info("Reciprocal Remap started")

        self.remapping_methods = remapping_methods
        self.remapping_params = remapping_params
        self.assembly_file = assembly_file
        self.type = type
        self.prefix = prefix
        self.threads = threads
        # self.minimum_coverage = minimum_coverage
        self.bin = bin
        self.logging_level = logging_level
        self.cleanup = cleanup
        self.r1 = r1.current
        self.r2 = r2.current

    def reciprocal_map(self, remap_target: Remap_Target) -> Mapping_Instance:
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

    def reference_map(self, remap_target: Remap_Target):
        rdir = os.path.join(
            self.remapping_methods.remap_software.dir,
            remap_target.name,
            "reference",
        )

        target_remap_drone = Remapping(
            self.r1,
            remap_target,
            self.remapping_methods,
            self.assembly_file,
            self.type,
            self.prefix,
            rdir,
            self.remapping_params,
            self.threads,
            r2=self.r2,
            bin=self.bin,
            logging_level=self.logging_level,
            cleanup=self.cleanup,
            log_dir=self.logdir,
        )

        target_remap_drone.run_remap()

        return target_remap_drone

    def assembly_map(self, reference_remap: Remapping):
        if len(reference_remap.mapped_contigs) == 0:
            return None

        output_directory = os.path.join(
            self.remapping_methods.remap_software.dir,
            reference_remap.target.name,
            "assembly",
        )

        assembly_target = Remap_Target(
            "none",
            f"{reference_remap.target.name}_assembly",
            "none",
            self.assembly_file,
            self.prefix,
            "description",
            reference_remap.mapped_contigs,
        )

        assembly_remap_drone = Remapping(
            reference_remap.mapped_subset_r1,
            assembly_target,
            self.remapping_methods,
            self.assembly_file,
            self.type,
            self.prefix,
            output_directory,
            self.remapping_params,
            self.threads,
            r2=reference_remap.mapped_subset_r2,
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
    combined_fasta_filename: str = "reference_combined.fasta"

    def __init__(
        self,
        remap_targets: List[Remap_Target],
        r1: Read_class,
        r2: Read_class,
        remapping_methods: SoftwareRemap,
        assembly_file: str,
        type: str,
        prefix,
        threads: int,
        bin: str,
        logging_level: int,
        cleanup: bool,
        remap_params: RemapParams,
        logdir="",
    ):
        super().__init__(
            r1,
            r2,
            remapping_methods,
            remap_params,
            assembly_file,
            type,
            prefix,
            threads,
            bin,
            logging_level,
            cleanup,
        )

        self.logger = logging.getLogger(__name__)
        if self.logger.hasHandlers():
            self.logger.handlers.clear()
        self.logger.propagate = False
        self.logger.setLevel(logging_level)
        self.logger.addHandler(logging.StreamHandler())
        self.logger.propagate = False
        self.logger.info("Mapping Manager started")
        self.logdir = logdir

        self.mapped_instances = []
        self.remap_targets = remap_targets
        self.target_taxids = set([str(target.taxid) for target in remap_targets])

        self.reads_before_processing = r1.read_number_clean + r2.read_number_clean
        self.reads_after_processing = (
            r1.get_current_fastq_read_number() + r2.get_current_fastq_read_number()
        )
        self.max_gaps = 0
        self.max_prop = 0
        self.max_mapped = 0
        self.max_depth = 0
        self.max_depthR = 0
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
                "error_rate",
                "quality_avg",
            ],
        )
        self.combined_fasta_path = os.path.join(
            self.remapping_methods.output_dir, self.combined_fasta_filename
        )
        self.combined_fasta_gz_path = self.combined_fasta_path + ".gz"
        self.remap_params = remap_params
        self.cmd = RunCMD(bin, logdir=self.logdir, prefix=prefix, task="remapping_main")

    def check_targets_combined_fasta_exists(self):
        return os.path.exists(self.combined_fasta_gz_path)

    def generate_remap_targets_fasta(self):
        if self.check_targets_combined_fasta_exists():
            return
        open(self.combined_fasta_path, "w", encoding="utf-8").close()

        if os.path.exists(self.combined_fasta_gz_path):
            os.remove(self.combined_fasta_gz_path)

        for target in self.remap_targets:
            for accid in target.accid_in_file:
                accid_clean = (
                    accid.replace(";", "_")
                    .replace(":", "_")
                    .replace(".", "_")
                    .replace("|", "_")
                )
                tmp_fasta = os.path.join(
                    self.remapping_methods.output_dir, f"{accid_clean}.fasta"
                )
                cmd = f"samtools faidx {target.file} '{accid}' > {tmp_fasta}"

                self.cmd.run(cmd)

                if self.check_fasta_empty(tmp_fasta) is False:
                    self.append_fasta(tmp_fasta)

                os.remove(tmp_fasta)

        cmd_bgzip = f"bgzip {self.combined_fasta_path}"
        self.cmd.run(cmd_bgzip)

    def check_fasta_empty(self, fasta):
        """check if empty or only header"""

        with open(fasta, "r") as f:
            for line in f:
                if not line.startswith(">"):
                    return False
        return True

    def append_fasta(self, fasta):
        """append fasta to combined fasta"""

        with open(self.combined_fasta_path, "a") as f:
            with open(fasta, "r") as f2:
                for line in f2:
                    f.write(line)

    def run_mappings(self):
        """
        Run mappings for all targets."""
        for target in self.remap_targets:
            mapped_instance = self.reciprocal_map(target)

            self.mapped_instances.append(mapped_instance)

    def export_reference_fastas_if_failed(self, media_dir):
        for target in self.mapped_instances:
            target.reference.relocate_reference_fasta(media_dir)

    def run_mappings_move_clean(self, static_plots_dir, media_dir):
        for target in self.remap_targets:
            print(f"################## MAPPING ACC {target.accid} ##################")
            mapped_instance = self.reciprocal_map(target)

            self.mapped_instances.append(mapped_instance)

            apres = mapped_instance.reference.number_of_contigs_mapped > 0
            rpres = mapped_instance.reference.number_of_reads_mapped > 0

            print(f" Reads mapped: {rpres}")
            print(f" Contigs mapped: {apres}")

            if rpres:
                mapped_instance.reference.move_coverage_plot(static_plots_dir)
                mapped_instance.export_mapping_files(media_dir)
            else:
                print("No reads mapped, skipping coverage plot")
                mapped_instance.reference.cleanup_files()

            if apres:
                mapped_instance.reference.move_dotplot(static_plots_dir)
            else:
                print("No contigs mapped, skipping dotplot")
                if mapped_instance.assembly:
                    mapped_instance.assembly.cleanup_files()

    def move_igv_to_static(self, static_dir):
        print("Moving IGV files to static")
        for instance in self.mapped_instances:
            if instance.reference.number_of_reads_mapped > 0:
                instance.reference.move_igv_files(static_dir)

    def export_mapping_files(self, output_dir):
        for instance in self.mapped_instances:
            instance.export_mapping_files(output_dir)

    def merge_mapping_reports(self):
        full_report = []

        for instance in self.mapped_instances:
            ntax = instance.generate_full_mapping_report_entry()
            if len(ntax):
                full_report.append(ntax)

        if len(full_report) > 0:
            self.report = pd.concat(full_report, axis=0)

        if self.cleanup:
            self.clean_final_report()

    def clean_final_report(self):
        self.report.ngaps = self.report.ngaps.fillna(0)
        self.report = self.report[self.report.coverage > 0]
        self.report = self.report.sort_values(["coverage", "Hdepth"], ascending=False)

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

    def verify_mapped_instance(self, mapped_instance: Mapping_Instance):
        if (
            mapped_instance.reference.r1 == self.r1
            and mapped_instance.reference.r2 == self.r2
        ):
            return True
        else:
            return False

    def validate_mapped_instance_taxid(self, mapped_instance: Mapping_Instance):
        if str(mapped_instance.reference.target.taxid) in self.target_taxids:
            return True
        else:
            return False

    def update_mapped_instance_safe(self, mapped_instance: Mapping_Instance):
        if self.verify_mapped_instance(mapped_instance):
            if self.validate_mapped_instance_taxid(mapped_instance):
                self.mapped_instances.append(mapped_instance)

    def update_mapped_instance(self, mapped_instance: Mapping_Instance):
        if self.validate_mapped_instance_taxid(mapped_instance):
            self.mapped_instances.append(mapped_instance)

    def update_mapped_instances(self, mapped_instances: List[Mapping_Instance]):
        print(self.r1)
        print(self.target_taxids)

        self.mapped_instances = []
        for instance in mapped_instances:
            self.update_mapped_instance(instance)

    def get_mapped_instance(self, taxid):
        for instance in self.mapped_instances:
            if instance.reference.target.taxid == taxid:
                return instance
        return None
