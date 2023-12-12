import http.client
import logging
import os
from typing import List, Optional

import pandas as pd

from pathogen_identification.constants_settings import ConstantsSettings as CS
from pathogen_identification.models import (
    PIProject_Sample,
    RawReference,
    ReferenceSourceFileMap,
    RunMain,
)
from pathogen_identification.modules.object_classes import Remap_Target
from pathogen_identification.utilities.entrez_wrapper import EntrezWrapper
from pathogen_identification.utilities.utilities_general import (
    description_fails_filter,
    merge_classes,
    scrape_description,
)
from pathogen_identification.utilities.utilities_pipeline import RawReferenceUtils


class RunMetadataHandler:
    remap_targets: List[Remap_Target] = []

    def __init__(
        self,
        username,
        config,
        sift_query: str = "phage",
        prefix: str = "",
        rundir: str = "",
    ):
        """
        Initialize metadata handler.

        Args:
            metadata_paths: dictionary of paths to metadata files.
            sift_query: string to filter sift report.

        """
        self.prefix = prefix
        self.rundir = rundir
        self.config = config
        self.metadata_paths = config["metadata"]
        self.logger = logging.getLogger(f"{__name__}_{self.prefix}")
        self.logger.setLevel(logging.INFO)
        self.logger.addHandler(logging.StreamHandler())
        if self.logger.hasHandlers():
            self.logger.handlers.clear()
        self.logger.propagate = False

        self.entrez_conn = EntrezWrapper(
            username,
            bindir=os.path.join(
                self.config["bin"]["ROOT"],
                self.config["bin"]["software"]["entrez_direct"],
                "bin",
            ),
            outdir=self.rundir,
            outfile="entrez_output.tsv",
        )

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
        self.merged_targets: pd.DataFrame = pd.DataFrame()
        self.remap_targets: List[Remap_Target] = []
        self.remap_absent_taxid_list: List[str] = []
        self.remap_plan = pd.DataFrame
        self.sift_query = sift_query
        self.sift_report = pd.DataFrame(
            [[0, 0, 0]], columns=["input", "output", "removed"]
        )
        self.get_metadata()

    def get_manual_references(self, sample: PIProject_Sample):
        references = RawReference.objects.filter(
            run__sample=sample,
            run__run_type=RunMain.RUN_TYPE_STORAGE,
        )
        fasta_main_dir = self.config["source"]["REF_FASTA"]

        for ref in references:
            refmap = ReferenceSourceFileMap.objects.filter(
                reference_source__taxid__taxid=ref.taxid,
                reference_source__accid=ref.accid,
            ).first()

            if refmap is None:
                continue

            accid_simple = (
                ref.accid.replace(".", "_")
                .replace(";", "_")
                .replace(":", "_")
                .replace("|", "_")
            )

            self.remap_targets.append(
                Remap_Target(
                    ref.accid,
                    accid_simple,
                    ref.taxid,
                    os.path.join(fasta_main_dir, refmap.reference_source_file.file),
                    self.prefix,
                    ref.description,
                    [ref.accid],
                    False,
                    False,
                )
            )

    def merge_sample_references(
        self, sample_registered: PIProject_Sample, max_taxids: int
    ):
        """
        Generate Remap Targets from all existing references for a given sample."""
        reference_utils = RawReferenceUtils(sample_registered)
        reference_utils.sample_reference_tables()
        reference_table = reference_utils.merge_ref_tables()

        proxy_rclass = reference_utils.reference_table_renamed(
            reference_table, {"read_counts": "counts"}
        )
        proxy_aclass = reference_utils.reference_table_renamed(
            reference_table, {"contig_counts": "counts"}
        )

        self.rclass = proxy_rclass
        self.aclass = proxy_aclass

        self.merge_reports_clean(
            max_taxids,
        )

        self.generate_targets_from_report(
            reference_table,
            max_taxids=max_taxids,
            skip_scrape=False,
        )

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
        if self.merged_targets.empty:
            self.merge_reports_clean(
                taxid_limit=taxid_limit,
            )

        #######
        #######

        self.generate_mapping_targets(
            self.merged_targets,
            prefix=self.prefix,
            max_remap=max_remap,
        )

    @staticmethod
    def prettify_reports(df: pd.DataFrame) -> pd.DataFrame:
        if "acc_x" in df.columns:
            if "accid" in df.columns:
                df = df.drop(columns=["acc_x"])
            else:
                df = df.rename(columns={"acc_x": "accid"})

        if "acc_y" in df.columns:
            if "accid" in df.columns:
                df = df.drop(columns=["acc_y"])
            else:
                df = df.rename(columns={"acc_y": "accid"})

        if "counts" in df.columns:
            if "counts_x" in df.columns:
                df = df.drop(columns=["counts_x"])
            if "counts_y" in df.columns:
                df = df.drop(columns=["counts_y"])

        return df

    @staticmethod
    def clean_report(df: pd.DataFrame):
        """
        Clean report.
        """

        if df.shape[0] > 0:
            for target_col in ["acc", "protid", "prot_acc", "taxid"]:
                if target_col in df.columns:
                    df = df.dropna(subset=[target_col])
                    # df = df.drop_duplicates(subset=[target_col])
                    df = df.reset_index(drop=True)

        return df

    def filter_references_table(self, references_table: pd.DataFrame) -> pd.DataFrame:
        references_table["taxid"] = references_table["taxid"].astype(str)

        references_table = references_table[references_table.taxid != "0"]
        references_table = references_table[references_table.taxid != "1"]
        references_table = references_table[
            ~references_table.description.isin(["root", "NA"])
        ]
        references_table = references_table[~references_table.accid.isin(["-"])]

        references_table["taxid"] = references_table["taxid"].astype(int)

        return references_table

    def generate_targets_from_report(
        self,
        df: pd.DataFrame,
        max_taxids: Optional[int] = None,
        skip_scrape: bool = True,
    ):
        references_table = self.filter_references_table(df)

        # references_table = references_table.drop_duplicates(subset=["taxid"])
        references_table.rename(columns={"accid": "acc"}, inplace=True)

        ## group by taxids
        if references_table.shape[0] > 0:
            references_table = (
                references_table.groupby(["taxid"])
                .agg(
                    {
                        "acc": "first",
                        "description": "first",
                        "read_counts": "first",
                        "standard_score": "first",
                        "contig_counts": "first",
                    }
                )
                .reset_index()
            )

        if "standard_score" in references_table.columns:
            references_table = references_table.sort_values(
                by="standard_score", ascending=False
            )
            print("##### standard score ######")
            print(references_table.head(10))

        if max_taxids is not None:
            references_table = references_table.iloc[:max_taxids, :]

        self.generate_mapping_targets(
            references_table,
            prefix=self.prefix,
            max_remap=1,
            skip_scrape=skip_scrape,
        )

    def results_collect_metadata(
        self, df: pd.DataFrame, sift: bool = True
    ) -> pd.DataFrame:
        """
        Process results.
        merge df with metadata to create taxid columns.
        summarize merged dataframe to get counts per taxid.
        if sift is true, filter results to only include self.sift_query.
        """

        df = self.clean_report(df)

        df = self.merge_report_to_metadata_taxid(df)

        df = self.map_hit_report(df)

        df = self.entrez_get_taxid_descriptions(df)
        # df = self.entrez_conn.entrez_get_taxid_descriptions(df)
        # df = self.merge_report_to_metadata_description(df)

        df = df.reset_index(drop=True)

        def get_acc(df: pd.DataFrame):
            if "acc_x" in df.columns:
                df["accid"] = df["acc_x"]
                df.drop(columns=["acc_x"])
            elif "acc_y" in df.columns:
                df["accid"] = df["acc_y"]
                df.drop(columns=["acc_y"])
            elif "acc" in df.columns:
                df["accid"] = df["acc"]
                df.drop(columns=["acc"])
            else:
                df["accid"] = df["taxid"].apply(self.get_taxid_representative_accid)

            return df

        if df.shape[0] > 0:
            df = get_acc(df)
            df["description"] = df["description"].fillna("NA")

        if sift:
            sifted_df = self.sift_report_filter(df, query=self.sift_query)
            self.sift_report = self.sift_summary(df, sifted_df)
            df = sifted_df

        df = self.prettify_reports(df)

        return df

    def get_metadata(self):
        """
        Get metadata from files.
        """
        try:
            self.accession_to_taxid = pd.read_csv(
                self.input_accession_to_taxid_path, sep="\t", header=0
            )
        except:
            self.accession_to_taxid = pd.DataFrame(columns=["acc", "taxid"])
            self.logger.info("No accession to taxid file found.")
            self.logger.info(
                "This file is required for mapping, check installation. Exiting."
            )
            exit()

        try:
            self.taxonomy_to_description = pd.read_csv(
                self.input_taxonomy_to_descriptor_path, sep="\t", header=0
            )
        except:
            self.taxonomy_to_description = pd.DataFrame(
                columns=["taxid", "description"]
            )
            self.logger.info("No taxonomy to description file found.")
            self.logger.info(
                "This file is required for mapping, check installation. Exiting."
            )
            exit()

        try:
            self.protein_to_accession = pd.read_csv(
                self.input_protein_accession_equivalent_path, sep="\t", header=0
            )
        except:
            self.protein_to_accession = pd.DataFrame(columns=["protid", "acc"])
            self.logger.info("No protein accession to protid file found.")

        try:
            self.protein_accession_to_taxid = pd.read_csv(
                self.input_protein_accession_to_taxid_path, sep="\t", header=0
            )
        except:
            self.protein_accession_to_taxid = pd.DataFrame(columns=["acc", "taxid"])
            self.logger.info("No protein accession to taxid file found.")

        self.logger.info("Finished retrieving metadata")

    def merge_report_to_metadata_taxid(self, df: pd.DataFrame) -> pd.DataFrame:
        """

        Args:
            df: classifier output, possessing at least columns: acc, protid, prot_acc or taxid.

        Returns:
            df: classifier output, possessing original columns plus: description.

        """

        if df.shape[0] == 0:
            return pd.DataFrame(columns=["taxid", "description", "file"])

        if "taxid" not in df.columns:
            if "prot_acc" in df.columns and "acc" not in df.columns:
                counts_df = df.groupby(["prot_acc"]).size().reset_index(name="counts")
                return self.merge_check_column_types(
                    counts_df, self.protein_accession_to_taxid, "prot_acc"
                )

            if "protid" in df.columns and "acc" not in df.columns:
                counts_df = df.groupby(["protid"]).size().reset_index(name="counts")
                df = self.merge_check_column_types(
                    counts_df, self.protein_to_accession, "protid"
                )

            if "acc" in df.columns:
                if "counts" in df.columns:
                    counts_df = df.groupby(["acc"]).agg({"counts": "sum"}).reset_index()

                else:
                    counts_df = df.groupby(["acc"]).size().reset_index(name="counts")

                df = self.merge_check_column_types(
                    counts_df,
                    self.accession_to_taxid,
                    column="acc",
                    column_two="acc_in_file",
                )

            else:
                raise ValueError(
                    "No taxid, accid or protid in the dataframe, unable to retrieve description."
                )

        df = df[(df.taxid != "0") | (df.taxid != 0)]

        return df

    def entrez_get_taxid_descriptions(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Get taxid descriptions from entrez.
        """

        return self.entrez_conn.entrez_get_taxid_descriptions(df)

    def merge_report_to_metadata_description(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Merge df with taxonomy to description file.
        """

        df["taxid"] = df["taxid"].astype(int)
        df = self.merge_check_column_types(df, self.taxonomy_to_description, "taxid")
        df = df.sort_values(by="taxid")
        df = df.drop_duplicates(subset=["taxid"], keep="first")

        df = df.dropna(subset=["taxid"])
        df.taxid = df.taxid.astype(float)
        df = df.dropna(subset=["taxid"])
        df.taxid = df.taxid.astype(int)
        df = df[df.taxid != 0]

        return df

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

        return pd.merge(df1, df2, left_on=column, right_on=column_two, how="left")

    @staticmethod
    def sift_report_filter(df, query: str = "phage"):
        """
        Filter df to only include hits with query in description column.
        """
        if "description" not in df.columns:
            df["description"] = ""

        ntab = df[~df.description.str.contains(query)]
        # ntab = ntab.drop_duplicates(subset=["qseqid"])

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
            return pd.DataFrame(columns=["taxid", "file", "counts"])

        if "counts" not in merged_table.columns:
            counts = merged_table.taxid.value_counts()
            counts = pd.DataFrame(counts).reset_index()
            counts.columns = ["taxid", "counts"]

            merged_table["taxid"] = merged_table["taxid"].astype(int)
            counts["taxid"] = counts["taxid"].astype(int)

            new_table = pd.merge(
                left=merged_table, right=counts, on="taxid"
            ).drop_duplicates(subset="taxid")

            new_table = new_table.sort_values("counts", ascending=False)
        else:
            new_table = (
                merged_table.groupby(["taxid"]).agg({"counts": "sum"}).reset_index()
            )
            new_table = new_table.sort_values("counts", ascending=False)

        return new_table

    def process_reports(
        self,
        report_1: pd.DataFrame,
        report_2: pd.DataFrame,
    ):
        self.rclass = self.results_collect_metadata(report_1)
        self.aclass = self.results_collect_metadata(report_2)

    def get_taxid_representative_accid(self, taxid: int) -> str:
        """
        Return representative accession for a given taxid.
        """
        if str(taxid) not in self.accession_to_taxid.taxid.astype(str).unique():
            return "-"

        accid_set = self.accession_to_taxid[
            self.accession_to_taxid.taxid.astype(str) == str(taxid)
        ].reset_index()

        accid_set = accid_set.dropna(subset=["acc"])

        if accid_set.shape[0] == 0:
            return "-"
        else:
            return accid_set.acc.iloc[0]

    def get_taxid_representative_accid_df(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Return representative accession for a given taxid.
        """

        new_df = df.copy()
        if "acc" in new_df.columns:
            new_df.drop(columns=["acc"], inplace=True)

        new_df = (
            new_df.merge(self.accession_to_taxid, on="taxid", how="left")
            .sort_values(["taxid", "acc"], ascending=[True, True])
            .drop_duplicates(subset=["taxid"], keep="first")
            .fillna("-")
        )

        return new_df

    def get_taxid_representative_description(self, taxid: int) -> str:
        """
        Return representative accession for a given taxid.
        """
        if str(taxid) not in self.taxonomy_to_description.taxid.astype(str).unique():
            return "-"

        desc_set = self.taxonomy_to_description[
            self.taxonomy_to_description.taxid.astype(str) == str(taxid)
        ].reset_index()

        desc_set = desc_set.dropna(subset=["description"])

        if desc_set.shape[0] == 0:
            return "-"
        else:
            return desc_set.description.iloc[0]

    def merge_reports_clean(
        self,
        taxid_limit: int = 15,
    ):
        """merge the reports and filter them."""

        targets, raw_targets = merge_classes(self.rclass, self.aclass, maxt=taxid_limit)

        raw_targets["accid"] = raw_targets["taxid"].apply(
            self.get_taxid_representative_accid
        )

        if "description" not in raw_targets.columns:
            taxid_descriptions = pd.concat(
                [
                    self.rclass[["taxid", "description"]],
                    self.aclass[["taxid", "description"]],
                ]
            )
            taxid_descriptions.dropna(subset=["description"], inplace=True)
            taxid_descriptions.drop_duplicates(subset=["taxid"], inplace=True)
            raw_targets["taxid"] = raw_targets["taxid"].astype(int)
            taxid_descriptions["taxid"] = taxid_descriptions["taxid"].astype(int)
            raw_targets = raw_targets.merge(taxid_descriptions, on="taxid", how="left")

        raw_targets[
            "status"
        ] = False  # raw_targets["taxid"].isin(targets["taxid"].to_list())

        self.raw_targets = raw_targets
        self.merged_targets = targets

    @staticmethod
    def metadata_from_taxid(taxid: str) -> pd.DataFrame:
        ref_db = pd.DataFrame(columns=["taxid", "acc", "description", "file"])
        reference_source = ReferenceSourceFileMap.objects.filter(
            reference_source__taxid__taxid=taxid
        )

        if len(reference_source) == 0:
            return ref_db

        ref_db = []
        for ref in reference_source:
            ref_db.append(
                [
                    ref.reference_source.taxid.taxid,
                    ref.reference_source.accid,
                    ref.reference_source.description,
                    ref.reference_source_file.file,
                ]
            )

        ref_db = pd.DataFrame(ref_db, columns=["taxid", "acc", "description", "file"])

        def acc_on_file(row: pd.Series):
            acc = row.acc
            file = row.file
            acc_on_file = acc
            if file == "virosaurus90_vertebrate-20200330.fas.gz":
                acc_on_file = acc.split(".")[0]
                acc_on_file = f"{acc_on_file}:{acc_on_file};"

            return acc_on_file

        ref_db["acc_in_file"] = ref_db.apply(acc_on_file, axis=1)

        return ref_db

    def generate_mapping_targets(
        self,
        targets,
        prefix: str,
        max_remap: int = 9,
        skip_scrape: bool = False,
    ):
        """
        check for presence of taxid in targets in self.accession_to_taxid.
        if present, find every accession ID associated, and create a ReferenceMap object
        for each accession ID.

        """
        fasta_main_dir = self.config["source"]["REF_FASTA"]

        remap_targets = []
        remap_absent = []
        taxf = self.accession_to_taxid
        remap_plan = []
        targets.taxid = targets.taxid.astype(int)

        for taxid in targets.taxid.unique():
            nset = self.metadata_from_taxid(taxid)

            if nset.empty:
                remap_absent.append(taxid)

                nset = pd.DataFrame(columns=["taxid"])
                remap_plan.append([taxid, "none", "none", "none"])
                continue

            # nset = (
            #    taxf[taxf.taxid == taxid]
            #    .reset_index(drop=True)
            #    .drop_duplicates(subset=["acc"], keep="first")
            # )
            ###

            files_to_map = self.filter_files_to_map(nset)

            ####

            for fileset in files_to_map:
                nsu = nset[nset.file == fileset]

                added_counts = 0

                if nsu.shape[0] > max_remap:
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
                        self.taxonomy_to_description.taxid.astype(int) == int(taxid)
                    ].description.unique()

                    if len(description) == 0:
                        description = [""]

                    if len(description) > 1:
                        description = sorted(description, key=len)

                    description = description[0]
                    if skip_scrape is False:
                        description = scrape_description(pref, description)

                    if description_fails_filter(description, CS.DESCRIPTION_FILTERS):
                        continue

                    def determine_taxid_in_file(taxid, df: pd.DataFrame):
                        """
                        determine if an accession is in a dataframe.
                        """
                        if "taxid" in df.columns:
                            return str(taxid) in df.taxid.astype(str).unique()

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
                            determine_taxid_in_file(taxid, self.rclass),
                            determine_taxid_in_file(taxid, self.aclass),
                        )
                    )

                    remap_plan.append([taxid, pref, fileset, description])

                    added_counts += 1
                    if added_counts > max_remap:
                        break

        self.remap_plan = pd.DataFrame(
            remap_plan, columns=["taxid", "acc", "file", "description"]
        )

        self.remap_targets.extend(remap_targets)
        self.remap_absent_taxid_list.extend(remap_absent)
