import os
from dataclasses import dataclass

import pandas as pd

from pathogen_identification.utilities.ncbi_connect import (
    entrez_fetch_taxid_from_org_description_curate,
)


@dataclass
class MergeArguments:
    file1: str
    file2: str
    file3: str
    output: str


def merge_files(args: MergeArguments, temp_dir: str):
    televir_reports = read_televir(temp_dir + args.file1)
    rpip_panel = read_panel(temp_dir + args.file2, panel="Microorganisms (RPIP)")
    upip_panel = read_panel(temp_dir + args.file3, panel="Microorganisms (UPIP)")

    illumina_found = get_illumina_found([rpip_panel, upip_panel], temp_dir)
    telebac_found = process_televir(televir_reports)

    merged_panel = merge_panels(illumina_found, telebac_found)
    merged_panel.to_csv(args.output, sep="\t", index=False)


def process_televir(televir_reports):
    columns_to_keep = [
        "Sample",
        "Taxid",
        "Description",
        "accID",
        "Cov (%)",
        "Depth",
        "DepthC",
        "Mapped reads",
        "Windows Covered",
        "Warning",
    ]
    if "Run" in televir_reports.columns:
        columns_to_keep.append("Run")

    telebac_found = televir_reports[columns_to_keep].drop_duplicates()

    telebac_found["Windows Covered"] = telebac_found["Windows Covered"].apply(
        lambda x: x.replace("/", "-")
    )
    telebac_found["Taxid"] = telebac_found["Taxid"].astype(str)

    return telebac_found


def read_televir(televir_filepath):
    """
    read tsv file from televir export
    """
    reports = pd.read_csv(televir_filepath, sep="\t")
    return reports


def read_panel(report, panel="Microorganisms"):
    """
    read excel extract spreadsheet name Microorganisms

    """
    panel = pd.read_excel(report, sheet_name=panel, engine="openpyxl")
    return panel


def get_illumina_found(panel_list: list, tmp_dir: str = "/tmp"):
    """
    merge rpip and upip panels, retrieve taxids
    """
    illumina_found_full = pd.concat(panel_list, axis=0, ignore_index=True)

    ## discard description= "None"
    illumina_found = illumina_found_full.drop_duplicates(
        subset=["Accession", "Microorganism Name"]
    )

    taxid_only = illumina_found[["Microorganism Name"]]

    def run_entrez_with_temp_dir(description: str):
        return entrez_fetch_taxid_from_org_description_curate(
            description, tmp_dir=tmp_dir
        )

    taxid_only["Taxid"] = taxid_only["Microorganism Name"].apply(
        run_entrez_with_temp_dir
    )

    illumina_found = pd.merge(
        illumina_found, taxid_only, on="Microorganism Name", how="left"
    )

    illumina_found = illumina_found[
        [
            "Accession",
            "Microorganism Name",
            "Class Type",
            "Coverage",
            "ANI",
            "Median Depth",
            "RPKM",
            "Taxid",
        ]
    ]

    # illumina_found["Taxid"] = illumina_found["Microorganism Name"].apply(
    #    entrez_fetch_taxid_from_org_description_curate
    # )

    illumina_found.rename(
        columns={"Accession": "Sample", "Microorganism Name": "Description"},
        inplace=True,
    )
    illumina_found = illumina_found.explode("Taxid").drop_duplicates()

    return illumina_found


def sample_in_televir(sample: str, televir_reports: pd.DataFrame) -> str:
    if sample in televir_reports["Sample"].values:
        return sample

    sample_safe = sample.replace("-", "_")
    if sample_safe in televir_reports["Sample"].values:
        return sample_safe
    for sample_televir in televir_reports["Sample"].values:
        if sample_safe in sample_televir:
            return sample_televir
    return sample


def merge_panels(illumina_found, telebac_found):
    """
    merge illumina_found and televir_reports
    """

    all_samples = illumina_found.Sample.unique()

    final_set = pd.DataFrame()
    print("Merging samples")
    for sample in all_samples:
        print(sample)

        illumina_sample = illumina_found[
            illumina_found["Sample"] == sample
        ].reset_index(drop=True)

        sample_televir = sample_in_televir(sample, telebac_found)

        telebac_sample = telebac_found[
            telebac_found["Sample"] == sample_televir
        ].reset_index(drop=True)

        new_set = pd.merge(
            illumina_sample, telebac_sample, on=["Sample", "Taxid"], how="outer"
        )

        new_set = new_set.drop_duplicates(
            subset=["Sample", "Description_x", "accID"]
        ).reset_index(drop=True)

        def support_for_taxid(taxid: str):
            in_illumina = taxid in illumina_sample["Taxid"].values
            in_telebac = taxid in telebac_sample["Taxid"].values

            if in_illumina and in_telebac:
                return "Both"
            elif in_illumina:
                return "Panels"
            elif in_telebac:
                return "TELEVir"
            else:
                return None

        def televir_run_support(taxid: str):
            in_telebac = taxid in telebac_sample["Taxid"].values
            if in_telebac:
                return telebac_sample[telebac_sample["Taxid"] == taxid]["Run"].nunique()
            else:
                return None

        if "Run" in new_set.columns:
            new_set["TELEVir Run Support"] = new_set["Taxid"].apply(televir_run_support)

        new_set["Support"] = new_set["Taxid"].apply(support_for_taxid)

        if "Run" in new_set.columns:
            new_set = new_set[
                [
                    "Sample",
                    "Support",
                    "Taxid",
                    "Description_y",
                    "accID",
                    "TELEVir Run Support",
                    "Cov (%)",
                    "Depth",
                    "DepthC",
                    "Mapped reads",
                    "Windows Covered",
                    "Warning",
                    "Description_x",
                    "Class Type",
                    "Coverage",
                    "ANI",
                    "Median Depth",
                    "RPKM",
                ]
            ]
        else:
            new_set = new_set[
                [
                    "Sample",
                    "Support",
                    "Taxid",
                    "Description_y",
                    "accID",
                    "Cov (%)",
                    "Depth",
                    "DepthC",
                    "Mapped reads",
                    "Windows Covered",
                    "Warning",
                    "Description_x",
                    "Class Type",
                    "Coverage",
                    "ANI",
                    "Median Depth",
                    "RPKM",
                ]
            ]

        new_set.rename(
            columns={
                "Description_x": "Description",
                "Description_y": "TELEVir Description",
            },
            inplace=True,
        )
        # sort by Mapped reads, then windows covered float
        new_set["Windows Covered"] = new_set["Windows Covered"].fillna("0-1")
        new_set["windows_covered_float"] = new_set["Windows Covered"].apply(
            lambda x: float(x.split("-")[0]) / float(x.split("-")[1])
        )

        new_set = new_set.sort_values(
            by=["Mapped reads", "windows_covered_float"], ascending=False
        ).reset_index(drop=True)
        # drop windows covered float
        new_set = new_set.drop(columns=["windows_covered_float"])

        ## drop duplicate taxids, keep the one with the highest Mapped reads
        new_set = new_set.drop_duplicates(subset=["Taxid"], keep="first")

        final_set = pd.concat([final_set, new_set], axis=0, ignore_index=True)

    # sort by sample, then support
    final_set = final_set.sort_values(
        by=["Sample", "Support"], ascending=False
    ).reset_index(drop=True)

    return final_set
