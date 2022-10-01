import os
import zipfile

import matplotlib
import pandas as pd
import requests
from bs4 import BeautifulSoup

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd


def simplify_name(name):
    """simplify sample name"""
    return (
        name.replace("_", "_")
        .replace("-", "_")
        .replace(" ", "_")
        .replace(".", "_")
        .lower()
    )


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


def scrape_description(accid: str, existing_description: str = None) -> str:
    """
    Scrape the description for the relevant information.
    """

    if accid.count("_") > 1:
        accid = accid.split("_")[:-1]
        accid = "_".join(accid)

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
    """Read the sam file and return a dataframe with the coordinates."""

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
    r2 = r2.sort_values("counts", ascending=False)

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


def merge_classes(r1, r2, maxt=6, exclude="phage"):
    """
    merge tables of taxids to columns.
    """
    if "description" in r1.columns:
        r1 = (
            r1[~r1.description.str.contains(exclude)]
            .drop_duplicates(subset=["taxid"], keep="first")
            .sort_values("counts", ascending=False)
        )

    r1 = r1[["taxid", "counts"]]

    r2pres = 1

    if len(r2):
        r2pres = 2
        if "description" in r2.columns:
            r2 = r2[~r2.description.str.contains(exclude)]

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

            r1 = (
                pd.concat([shared, r2, r1.head(maxt)], axis=0)
                .drop_duplicates(subset=["taxid"], keep="first")
                .reset_index(drop=True)
            )

    return r1.head(maxt * r2pres)
