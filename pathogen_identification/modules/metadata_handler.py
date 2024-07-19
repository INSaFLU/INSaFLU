import http.client
import logging
import os
from typing import List, Optional

import pandas as pd

from pathogen_identification.constants_settings import ConstantsSettings as CS
from pathogen_identification.models import (PIProject_Sample, RawReference,
                                            RawReferenceCompoundModel,
                                            ReferenceSource,
                                            ReferenceSourceFileMap, RunMain)
from pathogen_identification.modules.object_classes import Remap_Target
from pathogen_identification.utilities.entrez_wrapper import EntrezWrapper
from pathogen_identification.utilities.utilities_general import (merge_classes,
                                                                 simplify_name)
from pathogen_identification.utilities.utilities_pipeline import \
    RawReferenceUtils


def determine_taxid_in_file(taxid, df: pd.DataFrame):
    """
    determine if an accession is in a dataframe.
    """
    if "taxid" in df.columns:
        return str(taxid) in df.taxid.astype(str).unique()

    return False


class RunMetadataHandler:
    # remap_targets: List[Remap_Target] = []

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

    def reset(self):
        self.remap_targets: List[Remap_Target] = []
        self.remap_absent_taxid_list: List[str] = []
        self.remap_plan = pd.DataFrame
        self.sift_report = pd.DataFrame(
            [[0, 0, 0]], columns=["input", "output", "removed"]
        )

    def get_manual_references(self, sample: PIProject_Sample, max_accids: int = 15):
        """
        Get manual references for a given sample. update map request with references.
        """

        references = RawReference.objects.filter(
            run__sample=sample,
            run__run_type=RunMain.RUN_TYPE_STORAGE,
        )

        self.update_map_request(references)

    def get_mapping_references(self, run_pk: int, max_accids: int = 15):
        """
        Get mapping references for a given run. update map request with references.
        """

        references = RawReference.objects.filter(run__pk=run_pk)

        self.update_map_request(references, max_accids=max_accids)

    def update_map_request(self, references: List[RawReference], max_accids: int = 15):
        """
        Update the remap_targets list with references from list"""

        for ref in references:
            refmaps = ReferenceSourceFileMap.objects.filter(
                reference_source__taxid__taxid=ref.taxid,
                reference_source__accid=ref.accid,
            )
            accids_replete = 0
            for refmap in refmaps:
                if accids_replete > max_accids:
                    break

                if ref.accid is None:
                    continue
                accid_simple = simplify_name(ref.accid)

                self.remap_targets.append(
                    Remap_Target(
                        ref.accid,
                        accid_simple,
                        ref.taxid,
                        refmap.reference_source_file.filepath,
                        self.prefix,
                        ref.description,
                        [refmap.accid_in_file],
                        False,
                        False,
                    )
                )

                accids_replete += 1

    def merge_sample_references_classic_compound(
        self, sample_registered: PIProject_Sample, max_taxids: int, max_remap: int = 15
    ):
        """
        Generate Remap Targets from all existing references for a given sample."""
        reference_utils = RawReferenceUtils(sample_registered)
        reference_utils.sample_reference_tables()
        reference_table = reference_utils.merged_table

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
            max_remap=max_remap,
            skip_scrape=False,
        )

    def merge_sample_references_ensemble(
        self,
        sample_registered: PIProject_Sample,
        max_taxids: Optional[int] = None,
        max_remap: int = 15,
    ):

        reference_utils = RawReferenceUtils(sample_registered)
        ### ############################################################# ###
        compound_refs: List[RawReferenceCompoundModel] = (
            reference_utils.query_sample_compound_references()
        )

        if max_taxids is not None:
            compound_refs = compound_refs[:max_taxids]

        remap_plan = []
        remap_targets = []
        remap_absent_taxid_list = []

        for ref in compound_refs:
            ref_in_file = ReferenceSourceFileMap.objects.filter(
                reference_source__taxid__taxid=ref.taxid,
                reference_source__accid=ref.accid,
            )

            if len(ref_in_file) == 0:
                remap_absent_taxid_list.append(ref.taxid)
                continue

            files_to_map = self.filter_query_set_files(ref_in_file)
            ref_in_file = ref_in_file.filter(
                reference_source_file__file__in=files_to_map
            )

            target = Remap_Target(
                ref.accid,
                simplify_name(ref.accid),
                ref.taxid,
                ref_in_file[0].filepath,
                self.prefix,
                ref.description,
                [ref_in_file[0].accid_in_file],
                False,
                False,
            )

            remap_targets.append(target)
            remap_plan.append(
                [
                    ref.taxid,
                    ref.accid,
                    ref_in_file[0].reference_source_file.file,
                    ref.description,
                ]
            )

        self.remap_plan = pd.DataFrame(
            remap_plan, columns=["taxid", "acc", "file", "description"]
        )

        self.remap_targets.extend(remap_targets)
        self.remap_absent_taxid_list.extend(remap_absent_taxid_list)

    def match_and_select_targets(
        self,
        report_1: pd.DataFrame,
        report_2: pd.DataFrame,
        max_remap: int = 15,
        taxid_limit: int = 12,
    ):
        
        print("### processing reports")
        self.process_reports(
            report_1,
            report_2,
        )
        
        print("### merging reports")
        if self.merged_targets.empty:
            self.merge_reports_clean(
                taxid_limit=taxid_limit,
            )

        print("### generating targets")
        #######
        #######

        self.generate_mapping_targets(
            self.merged_targets,
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
        if "description" not in references_table.columns:
            references_table["description"] = ""

        references_table = references_table[
            ~references_table.description.isin(["root", "NA"])
        ]
        if (
            "accid" not in references_table.columns
            and "acc" not in references_table.columns
        ):
            references_table["accid"] = references_table["taxid"].apply(
                self.get_taxid_representative_accid
            )

        references_table = references_table[~references_table.accid.isin(["-"])]

        references_table["taxid"] = references_table["taxid"].astype(int)

        return references_table

    def generate_targets_from_report(
        self,
        df: pd.DataFrame,
        max_taxids: Optional[int] = None,
        max_remap: int = 15,
        skip_scrape: bool = True,
    ):
        references_table = self.filter_references_table(df)

        # references_table = references_table.drop_duplicates(subset=["taxid"])
        references_table.rename(columns={"accid": "acc"}, inplace=True)

        if max_taxids is not None:
            references_table = references_table.iloc[:max_taxids, :]

        print("############# REFERENCES TABLE #############")
        print(references_table.head())
        print(references_table.shape)

        self.generate_mapping_targets(
            references_table,
            max_remap=max_remap,
        )

    @staticmethod
    def filter_taxids_not_in_db(df) -> pd.DataFrame:

        def get_refs_existing(taxid):
            refs_in_file = ReferenceSourceFileMap.objects.filter(
                reference_source__taxid__taxid=taxid,
            ).distinct("reference_source__accid")
            return len(refs_in_file) > 0

        df["has_refs"] = df["taxid"].apply(get_refs_existing)
        df = df[df["has_refs"] == True]
        df.drop(columns=["has_refs"], inplace=True)
        return df

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

        df = self.filter_taxids_not_in_db(df)

        df = self.map_hit_report(df)

        df = self.db_get_taxid_descriptions(df)
        # df = self.entrez_get_taxid_descriptions(df)
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

    def get_protacc_taxid(self, df: pd.DataFrame) -> pd.DataFrame:
        query_list = df.prot_acc.unique().tolist()
        self.entrez_conn.bin_query = self.entrez_conn.bin_query_factory.get_query(
            "fetch_protein_accession_taxon"
        )
        output = self.entrez_conn.run_entrez_query(query_list)
        self.entrez_conn.bin_query = self.entrez_conn.bin_query_factory.get_query(
            "fetch_taxid_description"
        )
        # merge with df
        df = df.merge(output, left_on="prot_acc", right_on="acc", how="left")
        df = df.drop(columns=["prot_acc"])
        return df

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
                return self.get_protacc_taxid(counts_df)

            elif "protid" in df.columns and "acc" not in df.columns:
                counts_df = df.groupby(["protid"]).size().reset_index(name="counts")
                df = self.merge_check_column_types(
                    counts_df, self.protein_to_accession, "protid"
                )

            if "acc" in df.columns and "taxid" not in df.columns:
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

            if "taxid" not in df.columns:
                raise ValueError(
                    "No taxid, accid or protid in the dataframe, unable to retrieve description."
                )

        df = df[(df.taxid != "0") | (df.taxid != 0)]

        return df

    def db_get_taxid_descriptions(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Get taxid descriptions from database.
        """

        def get_description(taxid: str):
            try:
                return (
                    ReferenceSource.objects.filter(taxid__taxid=taxid)
                    .first()
                    .description
                )
            except:
                return ""

        df["description"] = df["taxid"].apply(get_description)

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

    def filter_query_set_files(self, nset: List[ReferenceSourceFileMap]) -> list:

        ref1 = "refseq_viral.genome.fna.gz"
        ref2 = "virosaurus90_vertebrate-20200330.fas.gz"

        files_to_map = []
        files_count = [x.reference_source_file.file for x in nset]

        if ref1 in files_count:
            files_to_map.append(ref1)
        if ref2 in files_count:
            files_to_map.append(ref2)

        if len(files_to_map) == 0:
            files_to_map.append(files_count[0])
        if len(files_to_map) == 1 and len(files_count) > 1:
            files_count = [x for x in files_count if x != files_to_map[0]]
            files_to_map.append(files_count[0])

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

        sources = ReferenceSource.objects.filter(taxid__taxid=taxid)
        sources = [x for x in sources if x.accid != "-" and x.accid is not None]

        if len(sources) == 0:
            return "-"
        else:
            return sources[0].accid

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

        raw_targets["status"] = (
            False  # raw_targets["taxid"].isin(targets["taxid"].to_list())
        )

        self.raw_targets = raw_targets
        self.merged_targets = targets

    def generate_mapping_targets(
        self,
        targets: pd.DataFrame,
        max_remap: int = 9,
    ):

        print(
            "######################## GENERATING TARGETS ############################"
        )
        print(targets.head())
        print(targets.shape)
        remap_plan = []
        remap_targets = []
        remap_absent_taxid_list = []

        for taxid in targets.taxid.unique():

            refs_in_file = ReferenceSourceFileMap.objects.filter(
                reference_source__taxid__taxid=taxid,
            ).distinct("reference_source__accid")

            if len(refs_in_file) == 0:
                print("skipping taxid with no references", taxid)
                remap_absent_taxid_list.append(taxid)
                continue

            #
            refs_in_file = refs_in_file[:max_remap]
            print("#refs in file", len(refs_in_file))

            for ref_in_file_by_accid in refs_in_file:
                other_refs = ReferenceSourceFileMap.objects.filter(
                    reference_source__taxid__taxid=taxid,
                    reference_source__accid=ref_in_file_by_accid.reference_source.accid,
                )

                files_to_map = self.filter_query_set_files(other_refs)
                ref_in_file = other_refs.filter(
                    reference_source_file__file__in=files_to_map
                ).first()

                target = Remap_Target(
                    ref_in_file.reference_source.accid,
                    simplify_name(ref_in_file.reference_source.accid),
                    ref_in_file.taxid,
                    ref_in_file.reference_source_file.filepath,
                    self.prefix,
                    ref_in_file.description,
                    [ref_in_file.accid_in_file],
                    determine_taxid_in_file(taxid, self.rclass),
                    determine_taxid_in_file(taxid, self.aclass),
                )

                remap_targets.append(target)
                remap_plan.append(
                    [
                        ref_in_file.reference_source.taxid.taxid,
                        ref_in_file.reference_source.accid,
                        ref_in_file.reference_source_file.file,
                        ref_in_file.reference_source.description,
                    ]
                )

        self.remap_plan = pd.DataFrame(
            remap_plan, columns=["taxid", "acc", "file", "description"]
        )
        print("ABSENT TAXIDS")
        print(remap_absent_taxid_list)
        print(len(remap_targets))
        self.remap_targets.extend(remap_targets)
        self.remap_absent_taxid_list.extend(remap_absent_taxid_list)
