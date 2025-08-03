from abc import ABC

import networkx as nx
import pandas as pd


class ClassifierOutputProcesseor(ABC):

    def __init__(self, class_output_path):
        self.output_path = class_output_path
        self.final_report: pd.DataFrame = pd.DataFrame(columns=["description", "taxID"])

    @classmethod
    def from_file(cls, class_output_path):
        pass

    @classmethod
    def process(self):
        pass

    def prep_final_report(self):
        self.final_report = self.final_report.drop_duplicates(subset=["description"])
        self.final_report["description"] = (
            self.final_report["description"].str.replace(" ", "_").str.lower()
        )
        return self

    def save(self, output_path):

        self.final_report.to_csv(output_path, sep="\t", index=False)


class CentrifugeOutputProcessor(ClassifierOutputProcesseor):

    def __init__(self, class_output_path, nuniq_threshold: int = 1):
        super().__init__(class_output_path)
        self.nuniq_threshold = nuniq_threshold

    def from_file(self):
        centrifuge_report = pd.read_csv(self.output_path, sep="\t")
        self.report = centrifuge_report
        return self

    def process(self):
        self.report = self.report[self.report["taxRank"] == "species"]
        self.report.sort_values("numUniqueReads", ascending=False, inplace=True)
        self.report = self.report[
            self.report["numUniqueReads"] >= self.nuniq_threshold
        ].reset_index(drop=True)
        self.final_report = self.report.rename(
            columns={
                "name": "description",
                "taxID": "taxid",
                "numUniqueReads": "counts",
            }
        )
        # self.final_report = self.report[["name", "taxID"]].rename(
        #    columns={"name": "description"}
        # )
        return self


def count_prefix_spaces(s):
    return len(s) - len(s.lstrip())


class KrakenOutputProcessor(ClassifierOutputProcesseor):

    def __init__(self, class_output_path, min_perc_reads: float = 10):
        super().__init__(class_output_path)
        self.min_perc_reads = min_perc_reads

    def from_file(self):

        kraken_report = pd.read_csv(self.output_path, sep="\t")
        kraken_report.columns = [
            "PercReads",
            "NumReadsRoot",
            "Nreads",
            "RankCode",
            "taxID",
            "name",
        ]
        kraken_report["prefix_spaces"] = kraken_report["name"].apply(
            count_prefix_spaces
        )

        nodes, edges = self.kraken_report_to_tree(kraken_report)
        self.nodes_dict = {node[0]: node for node in nodes}

        self.edges = edges
        self.report = kraken_report
        self.nodes = nodes

        return self

    def get_node_info(self, node):

        return self.nodes_dict[node]

    def process(self):
        leaves_simple = self.get_simplified_leaves(self.nodes, self.edges)
        leaves_simple = {
            self.get_node_info(parent): [self.get_node_info(leaf) for leaf in leaves]
            for parent, leaves in leaves_simple.items()
        }
        leaves_summary = self.summarize_leaves(leaves_simple)
        leaves_summary = leaves_summary[
            leaves_summary["perc_reads"] > self.min_perc_reads
        ]
        self.final_report = leaves_summary.rename(
            columns={
                "name": "description",
                "taxID": "taxid",
                "perc_reads": "PercReads",
                "Nreads": "counts",
            }
        )
        # self.final_report = self.final_report[["description", "taxID",

        return self

    @staticmethod
    def kraken_report_to_tree(report):
        taxid_dict = {}
        nodes_list = []
        edges_dict = {}
        edges = []

        for _, row in report.iterrows():
            name = row["name"].strip()
            tax_id = row["taxID"]
            prefix_spaces = count_prefix_spaces(row["name"])
            perc_reads = row["PercReads"]
            tax_rank = row["RankCode"]
            nreads = row["Nreads"]

            nodes_list.append(
                (tax_id, name, prefix_spaces, perc_reads, tax_rank, nreads)
            )
            taxid_dict[tax_id] = (name, prefix_spaces, perc_reads, tax_rank, nreads)

            # find parent node
            parent_node = None
            if prefix_spaces > 0:

                i = len(nodes_list) - 1
                while i > 0 and nodes_list[i][2] >= prefix_spaces:
                    i -= 1
                if i >= 0:
                    parent_node = nodes_list[i]

            if parent_node:
                parent_tax_id = parent_node[0]
                if parent_tax_id not in edges_dict:
                    edges_dict[parent_tax_id] = []
                edges_dict[parent_tax_id].append(tax_id)
                edges.append((parent_tax_id, tax_id))

        return nodes_list, edges

    @staticmethod
    def get_leaves(nodes, edges):
        G = nx.DiGraph()
        G.add_edges_from(edges)
        leaves = []
        for node in nodes:
            if G.out_degree(node[0]) == 0:
                leaves.append(node[0])
        return leaves

    @staticmethod
    def get_leaf_parents(nodes, edges):
        G = nx.DiGraph()
        G.add_edges_from(edges)
        leaf_parents = {}
        for node in nodes:
            if G.out_degree(node[0]) == 0:
                parent = list(G.predecessors(node[0]))
                if parent:
                    leaf_parents[node] = parent
        return leaf_parents

    def get_simplified_leaves(self, nodes, edges):
        """
        returns dict of {parent: [leaves]}
        """
        G = nx.DiGraph()
        G.add_edges_from(edges)
        ## starting from leaves:
        # - if leaf rank does not contain 'S', keep leaf.
        # - if leaf rank contains 'S', move up the until we find a parent without 'S' in its rank, keep the one before it.
        simplified_leaves = {}
        for node in nodes:
            if G.out_degree(node[0]) == 0:
                leaf = node[0]
                parent = list(G.predecessors(leaf))
                parent_detail = self.get_node_info(parent[0]) if parent else None
                while parent and "S" in parent_detail[4]:
                    leaf = parent[0]
                    parent = list(G.predecessors(leaf))
                    parent_detail = self.get_node_info(parent[0]) if parent else None
                if parent:
                    parent = parent[0]
                    if parent not in simplified_leaves:
                        simplified_leaves[parent] = []
                    simplified_leaves[parent].append(leaf)
        return simplified_leaves

    @staticmethod
    def summarize_leaves(leaves_simple) -> pd.DataFrame:
        """
        for each parent:
        if parent is species, return parent, else, return leaf with highest numUniqueReads
        """
        summary = []
        for parent, leaves in leaves_simple.items():
            if "S" in parent[4]:
                summary.append(parent)
            else:
                # find leaf with highest numUniqueReads
                best_leaf = max(leaves, key=lambda x: x[3])
                summary.append(best_leaf)
        return pd.DataFrame(
            summary,
            columns=[
                "taxID",
                "description",
                "prefix_spaces",
                "perc_reads",
                "rank_code",
                "Nreads",
            ],
        )
