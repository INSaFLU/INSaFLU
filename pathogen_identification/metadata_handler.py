import logging
import os
from typing import List

import pandas as pd

from pathogen_detection.object_classes import Remap_Target
from pathogen_detection.utilities import merge_classes, scrape_description


class Metadata_handler:
    def __init__(self, config, sift_query: str = "phage", prefix: str = ""):
        """
        Initialize metadata handler.

        Args:
            metadata_paths: dictionary of paths to metadata files.
            sift_query: string to filter sift report.

        """
        self.prefix = prefix
        self.config = config
        self.metadata_paths = config["metadata"]
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.logger.addHandler(logging.StreamHandler())

        self.input_taxonomy_to_descriptor_path = self.metadata_paths[
            "input_taxonomy_to_descriptor_path"
        ]
        self.input_accession_to_taxid_path = self.metadata_paths[
            "input_accession_to_taxid_path"
        ]
        self.input_protein_accession_equivalent_path = self.metadata_paths[
            "input_protein_accession_to_protid_path"
        ]

        self.input_protein_accession_to_taxid_path = self.metadata_paths[
            "input_protein_accession_to_taxid_path"
        ]

        self.accession_to_taxid: pd.DataFrame
        self.taxonomy_to_description: pd.DataFrame
        self.protein_to_accession: pd.DataFrame
        self.protein_accession_to_csv: pd.DataFrame

        self.rclass: pd.DataFrame
        self.aclass: pd.DataFrame
        self.merged_targets: pd.DataFrame
        self.remap_targets: List[Remap_Target]
        self.remap_absent_taxid_list: List[str]
        self.remap_plan = pd.DataFrame
        self.sift_query = sift_query
        self.sift_report = pd.DataFrame(
            [[0, 0, 0]], columns=["input", "output", "removed"]
        )
        self.get_metadata()

    def match_and_select_targets(
        self,
        report_1: pd.DataFrame,
        report_2: pd.DataFrame,
        max_remap: int = 15,
        taxid_limit: int = 12,
    ):

        self.process_reports(
            report_1,
            report_2,
        )
        self.merge_reports_clean(
            max_remap=max_remap,
        )

        #######
        #######

        remap_targets, remap_absent = self.generate_mapping_targets(
            self.merged_targets,
            prefix=self.prefix,
            taxid_limit=taxid_limit,
            fasta_main_dir=self.config["source"]["REF_FASTA"],
        )

        self.remap_targets = remap_targets
        self.remap_absent_taxid_list = remap_absent

    def results_process(self, df: pd.DataFrame, sift: bool = True) -> pd.DataFrame:
        """
        Process results.
        merge df with metadata to create taxid columns.
        summarize merged dataframe to get counts per taxid.
        if sift is true, filter results to only include self.sift_query.
        """

        df = self.merge_report_to_metadata(df)

        df = self.map_hit_report(df)

        if sift:
            sifted_df = self.sift_report_filter(df, query=self.sift_query)
            self.sift_report = self.sift_summary(df, sifted_df)
            df = sifted_df

        return df

    def get_metadata(self):
        """
        Get metadata from files.
        """

        self.accession_to_taxid = pd.read_csv(
            self.input_accession_to_taxid_path, sep="\t", header=0
        )
        self.taxonomy_to_description = pd.read_csv(
            self.input_taxonomy_to_descriptor_path, sep="\t", header=0
        )
        self.protein_to_accession = pd.read_csv(
            self.input_protein_accession_equivalent_path, sep="\t", header=0
        )

        self.protein_accession_to_taxid = pd.read_csv(
            self.input_protein_accession_to_taxid_path, sep="\t", header=0
        )

        self.logger.info("Finished retrieving metadata")

    def merge_report_to_metadata(self, df: pd.DataFrame) -> pd.DataFrame:
        """

        Args:
            df: classifier output, possessing at least columns: acc, protid, prot_acc or taxid.

        Returns:
            df: classifier output, possessing original columns plus: description.

        """

        if df.shape[0] == 0:
            return pd.DataFrame(columns=["taxid", "description", "file"])

        if "taxid" not in df.columns:

            if "prot_acc" in df.columns:
                return self.merge_check_column_types(
                    df, self.protein_accession_to_taxid, "prot_acc"
                )

            if "protid" in df.columns:
                df = self.merge_check_column_types(
                    df, self.protein_to_accession, "protid"
                )

            if "acc" in df.columns:
                df = self.merge_check_column_types(
                    df, self.accession_to_taxid, column="acc", column_two="acc_in_file"
                )

            else:
                raise ValueError(
                    "No taxid, accid or protid in the dataframe, unable to retrieve description."
                )

        return self.merge_check_column_types(df, self.taxonomy_to_description, "taxid")

    @staticmethod
    def merge_check_column_types(
        df1: pd.DataFrame,
        df2: pd.DataFrame,
        column: str,
        column_two: str = "",
    ):
        """
        Merge df1 with df2 on given column, check that the column is of the same type.
        """

        if len(column_two) == 0:
            column_two = column

        if column not in df1.columns:
            df1[column] = ""
        if column_two not in df2.columns:
            df2[column_two] = ""
        if df1[column].dtype != str:
            df1[column] = df1[column].astype(str)
        if df2[column_two].dtype != str:
            df2[column_two] = df2[column_two].astype(str)

        return pd.merge(df1, df2, left_on=column, right_on=column_two)

    @staticmethod
    def sift_report_filter(df, query: str = "phage"):
        """
        Filter df to only include hits with query in description column.
        """
        if "description" not in df.columns:
            df["description"] = ""

        ntab = df[~df.description.str.contains(query)]
        ntab = ntab.drop_duplicates(subset=["qseqid"])

        return ntab

    @staticmethod
    def filter_files_to_map(nset) -> list:
        """return at most two files, give priority to refseq and virosaurus.

        args: pandas data frame of taxid, acc, files.
        """

        ref1 = "refseq_viral.genome.fna.gz"
        ref2 = "virosaurus90_vertebrate-20200330.fas.gz"

        files_to_map = []
        files_count = nset.file.value_counts()

        if ref1 in nset.file.unique():
            files_to_map.append(ref1)
        if ref2 in nset.file.unique():
            files_to_map.append(ref2)

        if len(files_to_map) == 0:
            files_to_map.append(files_count.index[0])
        if len(files_to_map) == 1 and files_count.shape[0] > 1:
            files_count = files_count[files_count.index != files_to_map[0]]
            files_to_map.append(files_count.index[0])

        return files_to_map

    def sift_summary(self, merged_report: pd.DataFrame, filtered_reads: pd.DataFrame):
        """
        generate report of the difference in sequence ids between merged_report and filtered_reads.
        """

        n_removed = merged_report.taxid.nunique() - filtered_reads.taxid.nunique()

        logging.info(f"{n_removed} reads sifted from the report")
        sift_report_df = pd.DataFrame(
            [
                [
                    merged_report.taxid.nunique(),
                    filtered_reads.taxid.nunique(),
                    n_removed,
                ]
            ],
            columns=["input", "output", "removed"],
        )

        return sift_report_df

    @staticmethod
    def map_hit_report(merged_table: pd.DataFrame):
        """
        create a report of the number of hits per taxid.
        """

        if merged_table.shape[0] == 0:
            return pd.DataFrame(columns=["taxid", "description", "file", "counts"])

        counts = merged_table.taxid.value_counts()
        counts = pd.DataFrame(counts).reset_index()
        counts.columns = ["taxid", "counts"]

        new_table = pd.merge(
            left=merged_table, right=counts, on="taxid"
        ).drop_duplicates(subset="taxid")

        new_table = new_table.sort_values("counts", ascending=False)

        return new_table

    def process_reports(
        self,
        report_1: pd.DataFrame,
        report_2: pd.DataFrame,
    ):
        self.rclass = self.results_process(report_1)
        self.aclass = self.results_process(report_2)

    def merge_reports_clean(
        self,
        max_remap: int = 15,
    ):
        """merge the reports and filter them."""

        targets = merge_classes(self.rclass, self.aclass, maxt=max_remap)
        targets.dropna(subset=["taxid"], inplace=True)
        targets["taxid"] = targets["taxid"].astype(int)

        self.merged_targets = targets

    def generate_mapping_targets(
        self,
        targets,
        prefix: str,
        taxid_limit: int = 9,
        fasta_main_dir: str = "",
    ):
        """
        check for presence of taxid in targets in self.accession_to_taxid.
        if present, find every accession ID associated, and create a ReferenceMap object
        for each accession ID.

        """

        remap_targets = []
        remap_absent = []
        taxf = self.accession_to_taxid
        remap_plan = []

        for taxid in targets.taxid.unique():

            if len(taxf[taxf.taxid == taxid]) == 0:
                remap_absent.append(taxid)

                nset = pd.DataFrame(columns=["taxid"])
                remap_plan.append([taxid, "none", "none"])
                continue

            nset = (
                taxf[taxf.taxid == taxid]
                .reset_index(drop=True)
                .drop_duplicates(subset=["acc"], keep="first")
            )
            ###
            files_to_map = self.filter_files_to_map(nset)

            ####

            for fileset in files_to_map:

                nsu = nset[nset.file == fileset]

                if nsu.shape[0] > taxid_limit:
                    nsu = nsu.drop_duplicates(
                        subset=["taxid"], keep="first"
                    ).reset_index()

                for pref in nsu.acc.unique():

                    nsnew = nsu[nsu.acc == pref].reset_index(drop=True)
                    pref_simple = (
                        pref.replace(".", "_")
                        .replace(";", "_")
                        .replace(":", "_")
                        .replace("|", "_")
                    )

                    self.taxonomy_to_description.taxid = (
                        self.taxonomy_to_description.taxid.astype(int)
                    )
                    description = self.taxonomy_to_description[
                        self.taxonomy_to_description.taxid == int(taxid)
                    ].description.unique()

                    if len(description) == 0:
                        description = [""]

                    if len(description) > 1:
                        description = sorted(description, key=len)

                    description = description[-1]
                    description = scrape_description(pref, description)

                    def determine_acc_in_file(acc, df: pd.DataFrame):
                        """
                        determine if an accession is in a dataframe.
                        """
                        if "acc" in df.columns:
                            return acc in df.acc.unique()
                        elif "acc_x" in df.columns and "acc_y" in df.columns:
                            return (acc in df.acc_x.unique()) or (
                                acc in df.acc_y.unique()
                            )
                        else:
                            return False

                    remap_targets.append(
                        Remap_Target(
                            pref,
                            pref_simple,
                            taxid,
                            os.path.join(fasta_main_dir, fileset),
                            prefix,
                            description,
                            [nsnew.acc_in_file[0]],
                            determine_acc_in_file(pref, self.rclass),
                            determine_acc_in_file(pref, self.aclass),
                        )
                    )
                    remap_plan.append([taxid, pref, fileset])

        self.remap_plan = pd.DataFrame(remap_plan, columns=["taxid", "acc", "file"])

        return remap_targets, remap_absent
