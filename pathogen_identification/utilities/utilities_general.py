import os
import zipfile

import matplotlib
import pandas as pd
import requests
from bs4 import BeautifulSoup

matplotlib.use("Agg")

import zipfile
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd

from fluwebvirus.settings import MEDIA_ROOT
from pathogen_identification.constants_settings import ConstantsSettings as CS
from pathogen_identification.models import Projects, RunMain


def generate_zip_file(file_list: list, zip_file_path: str) -> str:
    with zipfile.ZipFile(zip_file_path, "w") as zip_file:
        for file_path in file_list:
            zip_file.write(file_path, os.path.basename(file_path))

    return zip_file_path


def get_create_zip(file_list: list, outdir: str, zip_file_name: str) -> str:
    zip_file_path = os.path.join(outdir, zip_file_name)

    if os.path.exists(zip_file_path):
        os.unlink(zip_file_path)

    zip_file_path = generate_zip_file(file_list, zip_file_path)

    return zip_file_path


def description_fails_filter(description: str, filter_list: list):
    """
    Check if description contains any of the strings in filter_list
    """
    for filter_string in filter_list:
        if filter_string in description.lower():
            return True
    return False


def reverse_dict_of_lists(dict_of_lists: dict) -> dict:
    """
    Return dictionary of lists with keys as values and values as keys
    """
    return {value: key for key, values in dict_of_lists.items() for value in values}


def readname_from_fasta(fastafile) -> list:
    """
    Read in fasta file and return list of read names
    """
    read_names = []
    with open(fastafile) as f:
        for line in f:
            if line[0] == ">":
                read_names.append(line[1:].strip())
    return read_names


def simplify_name(name: str):
    """simplify sample name"""
    return name.replace(".", "_").replace(";", "_").replace(":", "_").replace("|", "_")


def simplify_name_lower(name: str):
    """simplify sample name"""
    return (
        name.replace("_", "_")
        .replace("-", "_")
        .replace(" ", "_")
        .replace(".", "_")
        .lower()
    )


def simplify_accid(accid):
    accid = accid.split("_")
    accid = accid.split("_")

    if len(accid) == 1:
        return accid[0]

    return "_".join(accid[:1])


def plot_dotplot(
    df: pd.DataFrame,
    out_file: str,
    title: str,
    inv_color: str = "red",
    xmax=0,
    borders=50,
):
    """Plot the dotplot from bamfile
    query and reference coordinates column names are ax, ay, bx, by.

    :param df: The dataframe with the dotplot.
    :param out_dir: The output directory.
    :param title: The title of the plot.
    :param inv_color: The color of the inverted regions.
    """
    fig, ax = plt.subplots(figsize=(11, 3))
    x_base = df.bx.min()

    for i, row in df.iterrows():
        x_coords = [row.ax, row.ay]
        y_coords = [row.bx + x_base, row.by + x_base]

        color = "black"
        if row.ay < row.ax:
            color = inv_color

        ax.plot(x_coords, y_coords, color=color)
        x_base += row.by

    ax.set_title(title)
    ax.set_xlabel("Reference")
    ax.set_ylabel("Contigs")
    if xmax:
        ax.set_xlim(0 - borders, xmax + borders)

    fig.savefig(out_file, bbox_inches="tight")
    ax.cla()
    fig.clf()
    plt.close("all")


def fastqc_parse(fastqc_path: str, stdin_fastqc_name: str = "stdin_fastqc"):
    """parse fastqc.
    returns a dict with the first 10 lines of the fastqc_data.txt file

    :param fastqc_path:
    :return: fastqc_data
    """

    if not os.path.isfile(fastqc_path):
        fqreads = pd.DataFrame(columns=["measure", "value"]).set_index("measure")
        return fqreads

    with zipfile.ZipFile(fastqc_path) as zf:
        fqreads = zf.read(f"{stdin_fastqc_name}/fastqc_data.txt").decode("utf-8")
        fqreads = fqreads.split("\n")[5:10]
        fqreads = [x.split("\t") for x in fqreads]
        fqreads = pd.DataFrame(fqreads, columns=["measure", "value"]).set_index(
            "measure"
        )
        fqreads.rename({"Total Sequences": "Total_Sequences"}, axis=0, inplace=True)

    return fqreads


def scrape_description(accid: str, existing_description: Optional[str] = None) -> str:
    """
    Scrape the description for the relevant information.
    """

    if accid.count("_") > 1:
        accid = accid.split("_")[:-1]
        accid = "_".join(accid)

    if accid in ["", "NA", "nan", "NAN", "N/A", "NAN", "-"]:
        return str(existing_description)

    url = f"https://www.ncbi.nlm.nih.gov/nuccore/{accid}"
    headers = {
        "Access-Control-Allow-Origin": "*",
        "Access-Control-Allow-Methods": "GET",
        "Access-Control-Allow-Headers": "Content-Type",
        "Access-Control-Max-Age": "3600",
        "User-Agent": "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:52.0) Gecko/20100101 Firefox/52.0",
    }

    try:
        req = requests.get(url, headers)
        soup = BeautifulSoup(req.content, "html.parser")
        title = soup.find_all("div", class_="rprtheader")[0].h1.text
    except:
        title = existing_description

    if title == "":
        return str(existing_description)
    else:
        return str(title)


def read_paf_coordinates(samfile: str) -> pd.DataFrame:
    """Read the paf file and return a dataframe with the coordinates."""

    try:
        df = pd.read_csv(samfile, sep="\t", header=None).rename(
            columns={0: "read_id", 7: "ax", 8: "ay", 2: "bx", 3: "by"}
        )

        df = df.sort_values("ax", ascending=True)

    except:
        df = pd.DataFrame(columns=["read_id", "ax", "ay", "bx", "by"])

    return df


def plotly_dotplot(df: pd.DataFrame):
    """Plot the dotplot using plotly."""

    if df.empty:
        return None

    import plotly.graph_objects as go
    from plotly.offline import plot

    x_base = df.bx.min()
    total_segdf = []

    for i, row in df.iterrows():
        x_coords = [row.ax, row.ay]
        y_coords = [row.bx + x_base, row.by + x_base]

        x_base += row.by
        segdf = pd.DataFrame(
            dict(x=x_coords, y=y_coords, label=[f"{i}_{row.read_id}"] * 2)
        )

        total_segdf.append(
            go.Scatter(x=x_coords, y=y_coords, mode="lines", name=f"{i}_{row.read_id}")
        )

    layout = {
        "xaxis_title": "Reference Sequence",
        "yaxis_title": "contigs",
        "legend": {"family": "Courier New, monospace", "size": 8},
        "font": {"family": "Courier New, monospace", "size": 8},
        "height": 250,
        "width": 760,
        "template": "none",
        # "colorway": px.colors.qualitative.Bold,
        "margin": dict(l=50, r=15, b=30, t=10, pad=4),
        "legend": dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    }

    plot_div = plot(
        {"data": total_segdf, "layout": layout},
        output_type="div",
    )

    return plot_div


def process_class(r2, maxt=6):
    """
    Process classification results.
    """
    r2 = r2.drop_duplicates(subset=["taxid"], keep="first")
    r2 = r2.reset_index(drop=True)
    # r2 = r2.sort_values("counts", ascending=False)

    taxids_tokeep = []
    nr2 = []

    if "length" in r2.columns:
        r2["length"] = r2["length"].astype(int)
        r2c = r2.copy().sort_values("length", ascending=False)

        for i in range(r2.shape[0]):
            if i < maxt:
                taxids_tokeep.append(r2c.taxid[i])
                nr2.append(r2c.loc[i])
                if r2.taxid.tolist()[i] not in taxids_tokeep:
                    taxids_tokeep.remove(r2.taxid[i])
                    nr2.append(r2.loc[i])
            else:
                break
    else:
        r2 = r2.head(maxt)

    if len(nr2):
        r2 = pd.concat(nr2, axis=1).T

    return r2


def merge_classes(r1: pd.DataFrame, r2: pd.DataFrame, maxt=6, exclude="phage"):
    """
    merge tables of taxids to columns.
    """

    ###
    if "description" in r1.columns:
        r1 = r1[~r1.description.str.contains(exclude)]
    if "description" in r2.columns:
        r2 = r2[~r2.description.str.contains(exclude)]

    ###

    r1 = r1[["taxid", "counts"]]
    r1_raw = r1.copy()
    r2_raw = r2.copy()

    ###
    if len(r2) > 0 and len(r1) > 0:
        full_descriptor = pd.merge(r1, r2, on="taxid", how="outer")
    elif len(r2) > 0:
        full_descriptor = r2.copy()
    else:
        full_descriptor = r1.copy()

    full_descriptor = full_descriptor.fillna(0)

    ###

    if len(r2) and len(r1):
        r1.taxid = r1.taxid.astype(str)
        r2.taxid = r2.taxid.astype(str)

        shared = pd.merge(r1, r2, on=["taxid"], how="inner").sort_values(
            "counts_x", ascending=False
        )

        maxt = maxt - shared.shape[0]

        if maxt < 0:
            r1 = shared

        else:
            r2 = (
                pd.merge(r2, shared, indicator=True, how="outer")
                .query('_merge=="left_only"')
                .drop("_merge", axis=1)
            )
            r2 = process_class(r2, maxt=maxt)
            r1p = r1[~r1.taxid.isin(shared.taxid)]

            r1 = (
                pd.concat([shared, r2, r1p.head(maxt - r2.shape[0])], axis=0)
                .drop_duplicates(subset=["taxid"], keep="first")
                .reset_index(drop=True)
            )

    elif len(r2) == 0:
        r1 = r1.head(maxt)

    elif len(r1) == 0:
        r1 = r2.head(maxt)

    ###

    def descriptor_description_remove(df: pd.DataFrame):
        if "description_x" in df.columns:
            df.drop("description_x", axis=1, inplace=True)
        if "description_y" in df.columns:
            df.drop("description_y", axis=1, inplace=True)
        if "description" in df.columns:
            df.drop("description", axis=1, inplace=True)

        return df

    def get_source(row):
        if row.counts_x > 0 and row.counts_y > 0:
            return 3
        elif row.counts_x > 0:
            return 1
        elif row.counts_y > 0:
            return 2

    def descriptor_sources(fd):
        if len(r1_raw) == 0:
            fd["source"] = 2
        elif len(r2_raw) == 0:
            fd["source"] = 1
        else:
            fd["source"] = fd.apply(get_source, axis=1)

        return fd

    def descriptor_counts(fd):
        if len(r1_raw) == 0 or len(r2_raw) == 0:
            if "counts" not in fd.columns:
                fd["counts"] = fd["counts_x"]
            if "counts_y" in fd.columns:
                fd["counts"] = fd["counts_y"] + fd["counts_x"]
            return fd
        else:
            fd["counts"] = fd.apply(get_counts, axis=1)

        return fd

    def get_counts(row):
        if row.counts_x > 0 and row.counts_y > 0:
            return f"{int(row.counts_x)} / {int(row.counts_y)}"
        elif row.counts_x > 0:
            return str(row.counts_x)
        elif row.counts_y > 0:
            return str(row.counts_y)

    full_descriptor["taxid"] = full_descriptor["taxid"].astype(int)
    full_descriptor = descriptor_description_remove(full_descriptor)
    full_descriptor = descriptor_sources(full_descriptor)
    full_descriptor = descriptor_counts(full_descriptor)

    r1["taxid"] = r1.taxid.astype(int)
    merged_final = full_descriptor[full_descriptor.taxid.isin(r1.taxid.to_list())]
    merged_final["source"] = merged_final.source.apply(
        lambda x: ["none", "reads", "contigs", "reads/contigs"][x]
    )
    merged_final = merged_final[["taxid", "counts", "source"]]

    return merged_final, full_descriptor


def get_project_dir(project: Projects):
    return os.path.join(
        MEDIA_ROOT, CS.televir_subdirectory, str(project.owner.pk), str(project.pk)
    )


def infer_run_media_dir(run_main: RunMain) -> Optional[str]:
    if run_main.params_file_path:
        params_exist = os.path.exists(run_main.params_file_path)
        media_classification_dir = os.path.dirname(run_main.params_file_path)
        media_dir = os.path.dirname(media_classification_dir)
        if os.path.isdir(media_dir):
            return media_dir

    elif run_main.processed_reads_r1:
        reads_r1_exist = os.path.exists(run_main.processed_reads_r1)

        media_dir = os.path.dirname(run_main.processed_reads_r1)
        if os.path.isdir(media_dir):
            return media_dir

    elif run_main.processed_reads_r2:
        reads_r2_exist = os.path.exists(run_main.processed_reads_r2)

        media_dir = os.path.dirname(run_main.processed_reads_r2)
        if os.path.isdir(media_dir):
            return media_dir

    return None
