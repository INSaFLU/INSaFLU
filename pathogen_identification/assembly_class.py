import gzip
import logging
import os
import shutil
from typing import Type

import pandas as pd
from Bio import SeqIO
from metagen_view.settings import STATICFILES_DIRS

from pathogen_detection.object_classes import RunCMD


class Assembly_init:
    def __init__(
        self,
        r1: str,
        assembly_dir: str,
        assembly_args: str,
        threads: int = 3,
        r2: str = "",
        bin: str = "",
        log_dir: str = "",
        prefix: str = "",
    ):

        self.r1 = r1
        self.r2 = r2
        self.assembly_dir = assembly_dir
        self.assembly_args = assembly_args
        self.threads = threads
        self.cmd = RunCMD(bin, logdir=log_dir, prefix=prefix, task="assembly")

        self.assembly_file_fasta = os.path.join(
            self.assembly_dir, self.assembly_file_name
        )

    def assembly_file_check(self):
        if not os.path.exists(self.assembly_file_fasta):
            return False
        else:
            return True


class Assembly_spades(Assembly_init):
    """
    Assembly with SPAdes
    """

    assembly_file_name = "scaffolds.fasta"

    def run_PE(self, threads: int = 3):
        """
        Assembly with SPAdes
        """
        if self.assembly_file_check():
            print("Assembly file already exists")
        else:
            print("Assembly file does not exist, running SPAdes")

            spades_cmd = "spades.py -1 {} -2 {} -o {} -t {} {}".format(
                self.r1, self.r2, self.assembly_dir, threads, self.assembly_args
            )

            self.cmd.run(spades_cmd)

    def run_SE(self, threads: int = 3):
        """
        Assembly with SPAdes
        """
        if self.assembly_file_check():
            print("Assembly file already exists")
        else:
            if "--meta" in self.assembly_args:
                print("Spades does not accept --meta for single end reads. Removing.")
                self.assembly_args = self.assembly_args.replace("--meta", "")

            print("Assembly file does not exist, running SPAdes")
            spades_cmd = "spades.py -s {} -o {} -t {} {}".format(
                self.r1, self.assembly_dir, threads, self.assembly_args
            )

            self.cmd.run(spades_cmd)


class Assembly_trinity(Assembly_init):
    """
    Assembly transcriptome from rna seq data with Trinity
    """

    assembly_file_name = "assembly"

    def run_PE(self, threads: int = 3):
        """
        Assembly with Trinity
        """
        if self.assembly_file_check():
            print("Assembly file already exists")
        else:
            print("Assembly file does not exist, running Trinity")
            trinity_cmd = "Trinity --seqType fq --max_memory {} --left {} --right {} --CPU {} --output {} {}".format(
                self.r1, self.r2, self.threads, self.assembly_dir, self.assembly_args
            )

            trinity_cmd = trinity_cmd.split(" ")
            trinity_cmd.extend(self.assembly_args)
            self.cmd.run(trinity_cmd)


class Assembly_flye(Assembly_init):
    """
    Assembly with Flye"""

    assembly_file_name = "assembly.fasta"

    def run_SE(self, threads: int = 3):
        """
        Assembly with Flye
        """
        if self.assembly_file_check():
            print("Assembly file already exists")
        else:
            print("Assembly file does not exist, running Flye")
            flye_cmd = "flye -t {} {} {} -o {}".format(
                threads, self.assembly_args, self.r1, self.assembly_dir
            )

            self.cmd.run_python(flye_cmd)


class Assembly_velvet(Assembly_init):
    """
    Assembly with Velvet
    """

    assembly_file_name = "contigs.fa"

    def run_PE(self, threads: int = 3):
        """
        Assembly with Velvet
        """
        if self.assembly_file_check():
            print("Assembly file already exists")
        else:
            print("Assembly file does not exist, running Velvet")
            velvet_cmd = "velveth {} {} -shortPaired -fastq -separate {} {} {}".format(
                self.assembly_dir, self.assembly_args, self.r1, self.r2
            )

            self.cmd.run(velvet_cmd)


class Assembly_raven(Assembly_init):
    """
    Assembly with Raven
    """

    assembly_file_name = "assembly.scaffolds.fasta"
    assembly_gfa = "assembly.gfa"

    def __init__(
        self,
        r1: str,
        assembly_dir: str,
        assembly_args: dict,
        threads: int = 3,
        r2: str = "",
        bin: str = "",
        log_dir: str = "",
        prefix: str = "",
    ):
        super().__init__(
            r1, assembly_dir, assembly_args, threads, r2, bin, log_dir, prefix
        )

        self.assembly_file_fasta = os.path.join(assembly_dir, self.assembly_file_name)
        self.assembly_file_gfa = os.path.join(assembly_dir, self.assembly_gfa)

    def gfa_to_assembly(self):
        """
        Convert gfa to fasta
        """

        if not os.path.isfile(self.assembly_file_gfa):
            open(self.assembly_file_gfa, "w").close()
            open(self.assembly_file_fasta, "w").close()

        with open(self.assembly_file_gfa, "r") as g:
            with open(self.assembly_file_fasta, "w") as f:
                for line in g:
                    if line.startswith("S"):
                        line = line.split()
                        f.write(f">{line[1]}\n{line[2]}\n")
                g.close()
                f.close()

    def run_SE(self, threads: int = 3):
        """
        Assembly with Raven
        """
        if self.assembly_file_check():
            print("Assembly file already exists")
        else:
            print("Assembly file does not exist, running Raven")
            raven_cmd = "raven -t {} --disable-checkpoints --graphical-fragment-assembly {} {} {}".format(
                threads, self.assembly_file_gfa, self.assembly_args, self.r1
            )
            raven_cmd = raven_cmd.split(" ")

            self.cmd.run(raven_cmd)
            self.gfa_to_assembly()


class Assembly_class:

    assemblers_available = {
        "spades": Assembly_spades,
        "raven": Assembly_raven,
        "trinity": Assembly_trinity,
        "velvet": Assembly_velvet,
        "flye": Assembly_flye,
    }

    def __init__(
        self,
        r1,
        assembly_method,
        assembly_type,
        min_scaffold_length: int = 50,
        r2: str = "",
        prefix: str = "assembly",
        threads: int = 3,
        bin: str = "",
        logging_level: int = logging.ERROR,
        log_dir: str = "",
    ):
        self.r1 = r1
        self.r2 = r2
        self.type = assembly_type
        self.cmd = RunCMD(bin, logdir=log_dir, prefix=prefix, task="assembly")

        self.assembly_args = assembly_method.args
        self.min_scaffold_length = min_scaffold_length
        self.threads = threads

        self.assembly_dir = assembly_method.output_dir
        self.output_dir = assembly_method.dir
        self.assembly_method = assembly_method
        self.assembler = Type[object]

        self.assembly_file_name = os.path.join(self.output_dir, prefix)
        self.assembly_file_fasta = self.assembly_file_name + ".fasta"
        self.assembly_file_fasta_gz = self.assembly_file_fasta + ".gz"
        self.assembly_file_fasta_gz_index = self.assembly_file_fasta_gz + ".fai"

        self.assembly_exists = False
        self.assembly_index_exists = False
        self.contig_names = []
        self.contig_summary = pd.DataFrame(
            [[0, 0]], columns=["contig_name", "contig_length"]
        )

        self.assembly_min = 0
        self.assembly_max = 0
        self.assembly_mean = 0
        self.assembly_number = 0

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging_level)
        self.logger.addHandler(logging.StreamHandler())
        self.logger.info("Assembly class initialized")
        self.logger.info("Assembly method: {}".format(self.assembly_method.name))
        self.logger.info("Assembly directory: {}".format(self.assembly_dir))
        self.logger.info("Assembly file fasta: {}".format(self.assembly_file_fasta))
        self.logger.info("Assembly args: {}".format(self.assembly_args))
        self.logger.info("Assembly bin: {}".format(self.assembly_method.bin))

    def run(self):
        self.assembly_configure()

        if not self.check_assembler_output():
            self.assembly_run()
        else:
            self.logger.info("Assembly file already exists")

        if not self.check_assembler_output():
            self.logger.error(f"Assembly failed. Exiting")
            return

        self.assembly_process()
        self.get_contig_summary()
        self.get_assembly_stats()

    def assembly_configure(self):
        """
        Configure assembly
        """

        try:
            assembly_class = self.assemblers_available[self.assembly_method.name]
        except KeyError:
            print("Assembly method not available")
            return False

        self.assembler = assembly_class(
            self.r1,
            self.assembly_dir,
            self.assembly_args,
            threads=self.threads,
            r2=self.r2,
            bin=self.assembly_method.bin,
            log_dir=self.cmd.logdir,
            prefix=self.cmd.prefix,
        )

        os.makedirs(self.assembly_dir, exist_ok=True)

    def assembly_run(self):
        """
        Run assembly
        """

        os.makedirs(self.output_dir, exist_ok=True)

        if self.type == "PE":
            self.assembler.run_PE(self.threads)
        elif self.type == "SE":
            self.assembler.run_SE(self.threads)

    def check_assembler_output(self):
        """
        Check if assembly output exists
        """
        if os.path.exists(self.assembler.assembly_file_fasta):
            return True
        else:
            return False

    def move_to_final_path(self):
        """
        move assembly output fasta to expected final path if it does not exist.
        """
        self.logger.info("Moving assembly file to final path")

        if self.assembly_file_fasta == self.assembler.assembly_file_fasta:
            print("Assembly file already exists")
        else:
            shutil.copy(self.assembler.assembly_file_fasta, self.assembly_file_fasta)

    def assembly_process(self):
        """
        process assembly output files: move to final fasta path, gzip, index.
        """

        if self.assembly_file_check_fasta_gz() and self.assembly_file_check_index():
            print("Assembly file already exists")
            self.assembly_exists = True
            self.assembly_index_exists = True
            return

        self.filter_fasta_move()
        self.assembly_bgzip()
        self.assembly_index()

        self.assembly_exists = self.assembly_file_check_fasta_gz()
        self.assembly_index_exists = self.assembly_file_check_index()

    def assembly_file_check_index(self):

        if os.path.isfile(self.assembly_file_fasta_gz_index):
            return True
        else:
            return False

    def assembly_file_check_fasta(self):
        if os.path.isfile(self.assembly_file_fasta):
            return True
        else:
            return False

    def check_gz_file_not_empty(self, file: str) -> bool:
        """
        check that gz file has contigs or reads. use grep to seach fir lines starting with > or @. return boolean.
        """
        cmd = "gunzip -c {} | grep '^>\|^@' | wc -l".format(file)
        number_of_sequences = int(self.cmd.run_bash_return(cmd))

        if number_of_sequences > 0:
            return True
        else:
            print("Number of sequences in {}: {}".format(file, number_of_sequences))
            return False

    def assembly_file_check_fasta_gz(self):
        if os.path.isfile(self.assembly_file_fasta_gz) and self.check_gz_file_not_empty(
            self.assembly_file_fasta_gz
        ):
            return True
        else:
            return False

    def assembly_bgzip(self):
        """
        Gzip assembly file, replace if it already exists.
        """
        if os.path.isfile(self.assembly_file_fasta_gz):
            os.remove(self.assembly_file_fasta_gz)

        self.cmd.run(
            f"bgzip -c {self.assembly_file_fasta} > {self.assembly_file_fasta_gz}"
        )

        if self.assembly_file_check_fasta_gz():
            os.remove(self.assembly_file_fasta)

    def assembly_index(self):
        """index assembly  file"""
        if self.assembly_file_check_index():
            print("Assembly file already indexed")
        else:
            self.cmd.run("samtools faidx %s" % self.assembly_file_fasta_gz)

    def filter_fasta_move(self):
        """
        Filter fasta file by length.
        """

        with open(self.assembler.assembly_file_fasta, "r") as handle:
            input_seq_iterator = SeqIO.parse(handle, "fasta")
            short_seq_iterator = (
                record
                for record in input_seq_iterator
                if len(record.seq) > int(self.min_scaffold_length)
            )

            with open(self.assembly_file_fasta, "w") as handle:
                for record in short_seq_iterator:

                    handle.write("{}\n".format(record.format("fasta")))

    def get_contig_summary(self):
        """
        read ocntig names and respective sequence lengths from final gzipped assembly fasta file.
        """

        contig_summary = []
        with gzip.open(self.assembly_file_fasta_gz, "rt") as handle:
            input_seq_iterator = SeqIO.parse(handle, "fasta")

            for record in input_seq_iterator:
                contig_summary.append([record.id, len(record.seq)])

        self.contig_summary = pd.DataFrame(
            contig_summary, columns=["contig_name", "contig_length"]
        )

        self.contig_names = self.contig_summary.contig_name.tolist()

    def get_assembly_stats(self):
        """Assembly contig summary_stats"""

        if self.contig_summary.shape[0] > 0:
            self.assembly_min = self.contig_summary["contig_length"].fillna(0).min()
            self.assembly_max = self.contig_summary["contig_length"].fillna(0).max()
            self.assembly_mean = self.contig_summary["contig_length"].fillna(0).mean()
            self.assembly_number = self.contig_summary.shape[0]
