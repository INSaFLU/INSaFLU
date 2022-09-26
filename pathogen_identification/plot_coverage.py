import itertools as it
import os
import sys

import matplotlib

matplotlib.use("Agg")
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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


class Bedgraph:
    """Class to store and work with bedgraph files


    Methods:
    read_bedgraph: reads bedgraph wiith coverage column.
    get_coverage_array: returns numpy array of coverage column
    get_bins: generates unzipped bins from start & end positions.
    plot_coverage: barplot of coverage by window in bdgraph.
    """

    def __init__(self, bedgraph_file, max_bars=7000, nbins=300):
        self.max_bars = max_bars
        self.nbins = nbins
        self.bedgraph = self.read_bedgraph(bedgraph_file)
        self.reduce_number_bars()
        self.bar_to_histogram()

    def read_bedgraph(self, coverage_file) -> pd.DataFrame:
        coverage = pd.read_csv(coverage_file, sep="\t", header=None).rename(
            columns={0: "read_id", 1: "start", 2: "end", 3: "coverage"}
        )

        return coverage

    def get_bins(self) -> np.ndarray:
        """
        Get the bins.

        :param coverage_file: The coverage file.
        """
        self.bedgraph["width"] = self.bedgraph.end - self.bedgraph.start

    def get_coverage_array(self, coverage: pd.DataFrame) -> np.ndarray:
        """
        Get the coverage of the remapping.

        :param coverage_file: The coverage file.
        """
        coverage_values = np.array(coverage.coverage.to_list())

        return coverage_values

    def get_bar_coordinates(self):
        """
        Get the bar coordinates.
        """
        self.bedgraph["width"] = self.bedgraph.end - self.bedgraph.start

        self.bedgraph["x"] = self.bedgraph.start + self.bedgraph.end / 2
        self.bedgraph["y"] = self.bedgraph.coverage

        return self.bedgraph

    def reduce_number_bars(self):
        """
        Reduce the number of bars.
        """

        if self.bedgraph.shape > (self.max_bars,):
            self.bedgraph = self.bedgraph[self.bedgraph.coverage > 0]
            self.bedgraph = self.bedgraph.sample(self.max_bars)

    def merge_bedgraph_rows(self):
        """
        Merge the rows of the bedgraph.
        """

        for ix in range(1, self.bedgraph.shape[0]):
            if self.bedgraph.iloc[ix - 1].end < (self.bedgraph.iloc[ix].start - 1):
                self.bedgraph.iloc[ix].end = self.bedgraph.iloc[ix].start - 1

    def bar_to_histogram(self):
        """
        Bar to histogram.
        """

        self.bedgraph["coord"] = (self.bedgraph.start + self.bedgraph.end) / 2
        self.coverage = [
            [self.bedgraph.iloc[x]["coord"]] * self.bedgraph.iloc[x]["coverage"]
            for x in range(self.bedgraph.shape[0])
        ]
        self.coverage = list(it.chain.from_iterable(self.coverage))

    def plot_coverage(self, output_file, borders=50, tlen=0):
        """
        Plot the coverage of the remapping.

        :param coverage_file: The coverage file. bedgraph produced with samtools.
        :param output_file: The output file.
        """

        fig, ax = plt.subplots(figsize=(11, 3))

        if len(self.coverage) <= 1:
            return

        start_time = time.perf_counter()

        ax.hist(
            self.coverage,
            bins=self.nbins,
            color="skyblue",
            edgecolor="none",
        )

        ax.set_xlabel("Reference")
        ax.set_ylabel("Coverage")

        ##
        xmax = self.bedgraph.end.max()
        if tlen:
            xmax = tlen
        ax.set_xlim(0 - borders, xmax + borders)
        ##

        fig.savefig(output_file, bbox_inches="tight")
        ax.cla()
        fig.clf()
        plt.close("all")


def get_args():
    """
    Get arguments.
    """
    try:
        import argparse

        parser = argparse.ArgumentParser(
            description="Plot the coverage of the remapping."
        )
        parser.add_argument("--coverage_file", help="The coverage file.", required=True)
        parser.add_argument(
            "--output_file",
            help="The output file.",
            required=False,
            default="coverage.png",
        )
        args = parser.parse_args()

        return args

    except argparse.ArgumentError as e:
        print("ArgumentError: {}".format(e))
        print(e)
        sys.exit(1)


def main():
    """
    Main function.
    """

    args = get_args()

    bedgraph = Bedgraph(args.coverage_file)
    bedgraph.plot_coverage(args.output_file)


if __name__ == "__main__":
    main()
