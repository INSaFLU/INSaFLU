import argparse
import logging
import os
from dataclasses import dataclass

import numpy as np
import pandas as pd
from Bio import AlignIO


@dataclass
class Args:
    msa: str
    refname: str
    odir: str
    output: str


class VcfAlgn:
    def __init__(
        self,
        input,
        reference,
        nucl=True,
        output="test.vcf",
        odir="test/",
    ):
        """
        Class for aligning and processing alignments to reference sequence, map to existing look-up table, and write output to .tsv file.

        :param input: input alignment file
        :param gene: gene name
        :param reference: reference sequence
        :param nucl: boolean, True if nucleotide alignment, False if amino acid alignment
        :param output: output file name
        :param odir: output directory

        :type input: str
        :type gene: str
        :type reference: str
        :type nucl: bool
        :type output: str
        :type odir: str


        :return: summary geno-pheno look-up table

        """

        self.input = input
        self.reference = reference
        self.nucl = nucl
        if odir and odir[-1] != "/":
            odir += "/"
        self.odir = odir

        os.makedirs(self.odir, exist_ok=True)
        self.output = self.odir + output

    def read_input(self):
        """
        read in input (multifasta protein alignment) and geno-pheno database ("look-up table" resulting from gpdb2lut script)
        :return: self with alignment, reference sequence, ref seq index, and database
        """
        # read in alignment

        self.Alnmt = AlignIO.read(self.input, "fasta")
        ref_seq = self.reference
        # compares the full refseq identifier/description in the alnmt, to the ref_seq provided in argument.
        # This is a fix to deal with identifiers that (1) include spaces, (2) might consist of separate words of which the first is not unique.
        self.ref_seq_index = [
            i for i, s in enumerate(self.Alnmt) if s.description == ref_seq
        ]
        if len(self.ref_seq_index) == 0:
            raise ValueError(
                "Reference sequence not found in alignment. Please check alignment and reference sequence."
            )
        else:
            self.ref_seq_index = self.ref_seq_index[0]

        self.ref_seq = self.Alnmt[self.ref_seq_index]

    def read_ref(self):
        """
        read in reference sequence
        :return: self with reference sequence

        """
        samples = [s.description for _, s in enumerate(self.Alnmt)]

        refseq = [
            s.seq for _, s in enumerate(self.Alnmt) if s.description == self.reference
        ]
        if len(refseq) == 0:
            logging.info("Reference sequence not found in alignment")
            self.refseq = None

        else:
            self.refseq = refseq[0]

        self.samples = samples

    def msa2snp(self, merge_dels=False):
        """
        code from https://github.com/pinbo/msa2snp/blob/master/msa2snp.py
        Copyright 2019 Junli Zhang <zhjl86@gmail.com>

        iterates over mutations only.
        generate dataframe with (1) positions that vary, (2) aa at these positions in the refseq, (3) any alternatives in the dataset (sep by /), (4-end) for each sequence the aa at these positions
        :return: self with allele overview dataframe
        """
        fasta = {s.description: s.seq for _, s in enumerate(self.Alnmt)}
        seqnames = [s.description for _, s in enumerate(self.Alnmt)]
        refid = self.reference
        refseq = str(fasta[refid])
        self.refseq = refseq

        refseq_nogap = refseq.replace("-", "")

        n = -1
        seq2 = {}  # new dictionary of rotated sequences with no-gap postions as keys

        for i in range(len(refseq)):  # i is the position of seq poition
            if refseq[i] != "-":
                n += 1
            templist = []
            for j in seqnames:  # j is seq name
                seq = fasta[j]
                templist.append(seq[i])
            if n not in seq2:
                seq2[n] = templist
            else:
                seq2[n] = [s1 + s2 for s1, s2 in zip(seq2[n], templist)]

        # output
        outlist = []

        for k in range(len(refseq_nogap)):
            seq = [
                w[0] + w[1:].replace("-", "") if len(w) > 1 else w for w in seq2[k]
            ]  # remove "-" in alleles in "A--"
            # seq= seq2
            alleles = set(seq)  # set: unique values

            if len(alleles) != 1:
                ref_allele = refseq_nogap[k]  # string
                alt_allele_set = alleles - set(ref_allele)  # set
                alt_allele = "/".join(alt_allele_set)  # string
                outlist.append(
                    [str(k + 1), ref_allele, alt_allele] + seq
                )  # K+1: starting from 1

        # merge continuous deletions: ATTT to A---
        header = ["Pos", "ref", "alt"] + seqnames

        if merge_dels:

            def allpos(alist, avalue):
                return [pos for pos, element in enumerate(alist) if element == avalue]

            outlist2 = []  # merged continuous deletions
            mm = []  # temporary merging of lists
            n = 0  # the start position of continious deletion
            for i in range(len(outlist) - 1):
                p1 = outlist[i]
                p2 = outlist[i + 1]
                if (
                    int(p1[0]) + 1 == int(p2[0])
                    and "-" in p1
                    and "-" in p2
                    and allpos(p1, "-") == allpos(p2, "-")
                ):
                    if mm:  # if mm already have merged lists
                        mm = [s1 + s2 for s1, s2 in zip(mm, p1)]
                    else:
                        mm = p1
                        n = p1[0]
                else:
                    if mm:
                        mm = [s1 + s2 for s1, s2 in zip(mm, p1)]
                        mm[0] = n
                        outlist2.append(
                            ["-" if "-" in w else w for w in mm]
                        )  # replace "-----" to "-"
                        mm = []
                    else:
                        outlist2.append(p1)
                    if i + 1 == len(outlist) - 1:
                        outlist2.append(p2)

            outlist = outlist2

        outframe = pd.DataFrame(outlist)
        outframe.columns = header
        outframe.Pos = outframe.Pos.astype(int)

        self.vcf = outframe
        self.snps = outframe.Pos
        return self

    def find_stops(self):
        """
        determine "stop" positions for each sequence, to be able to identify sequences with early stops
        ## see https://github.com/pinbo/msa2snp/blob/master/msa2snp.py for reference.
        :return: self with an array of stop positions for each sequence
        """

        seq_matrix = [list(x.seq) for x in self.Alnmt]
        seq_matrix = np.array(seq_matrix, dtype=str)
        gapm = np.array(seq_matrix == "-", dtype=int)

        gapn = np.cumsum(gapm, axis=1)
        gapn = gapn[self.ref_seq_index]

        early_array = [seq_matrix.shape[1]] * seq_matrix.shape[0]
        # early_array = np.array(early_array)

        early_stops = np.where(seq_matrix == "*")
        early_finds = []

        for x, z in enumerate(early_stops[0]):
            est = early_stops[1][x]

            early_array[early_stops[0][x]] = est + 1 - gapn[est]
            early_finds += [est + 1 - gapn[est]]
        self.early_stops = early_array

        return self

    def write_vcf(self):
        """
        write vcf to file
        write vcf in format:
        ##fileformat=VCFv4.2
        ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
        ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
        #CHROM POS ID REF ALT QUAL FILTER INFO
        first:
            format self.vcf dataframe to standard vcf notation of allele relative to reference.
            exclude explicit reference column.
        """

        vcf = self.vcf.reset_index(drop=True)
        vcf = vcf[vcf.alt.str.contains("[ATCG]", regex=True)]
        # explode where alt contains /: split into separate rows
        vcf = vcf.assign(alt=vcf.alt.str.split("/")).explode("alt")
        vcf = vcf[vcf.alt.str.contains("[ATCG]", regex=True)]
        vcf = vcf.reset_index(drop=True)

        with open(self.output, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("##INFO=<ID=DP,Number=1,Type=Integer,Description='Total Depth'>\n")
            f.write(
                "##INFO=<ID=AF,Number=A,Type=Float,Description='Allele Frequency'>\n"
            )
            header = [
                "#CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                *self.samples,
            ]
            f.write("\t".join(header) + "\n")

            for i in range(len(vcf)):
                ref = vcf.loc[i, "ref"]
                alt = vcf.loc[i, "alt"]
                pos = vcf.loc[i, "Pos"]
                chrom = self.reference
                alleles = vcf.loc[i, self.samples]
                alleles.replace(ref, "0/0:1", inplace=True)
                alleles.replace(alt, "1/1:1", inplace=True)
                alleles.replace("N", "./.:0", inplace=True)
                alleles.replace("-", "./.:0", inplace=True)

                alleles = "\t".join(alleles)

                f.write(
                    f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tDP=1;AF=1\t{alleles}\n"
                )

        return self


def parse_args() -> Args:
    parser = argparse.ArgumentParser()

    parser.add_argument("--msa", help="Path to MSA")
    parser.add_argument("--refname", help="Name of reference sequence", required=True)
    parser.add_argument("--odir", default=".", help="Path to output directory")
    parser.add_argument("--output", default="test.vcf", help="Name of output file")
    args = parser.parse_args()

    return Args(msa=args.msa, refname=args.refname, odir=args.odir, output=args.output)


if __name__ == "__main__":

    args = parse_args()
    msa = args.msa
    odir = args.odir
    output = args.output
    refname = args.refname

    vcf_algn = VcfAlgn(msa, refname, odir=odir, output=output)
    vcf_algn.read_input()
    vcf_algn.read_ref()
    vcf_algn = vcf_algn.msa2snp()
    vcf_algn = vcf_algn.write_vcf()
