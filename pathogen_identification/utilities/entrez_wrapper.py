import os
from abc import ABC, abstractmethod
from typing import List

import pandas as pd


class EntrezQuery(ABC):
    name: str

    def __init__(self, bindir: str):
        self.bindir = bindir

    @abstractmethod
    def query(self, query: List[str]) -> str:
        pass

    @abstractmethod
    def read_output(self, outptuf_path: str) -> List[str]:
        pass

    @property
    def epost(self) -> str:
        return os.path.join(self.bindir, "epost")

    @property
    def efetch(self) -> str:
        return os.path.join(self.bindir, "efetch")

    @property
    def esearch(self) -> str:
        return os.path.join(self.bindir, "esearch")

    @property
    def esummary(self) -> str:
        return os.path.join(self.bindir, "esummary")

    @property
    def xtract(self) -> str:
        return os.path.join(self.bindir, "xtract")


class EntrezFetchTaxidDescription(EntrezQuery):
    name: str = "fetch_taxid_description"
    db = "taxonomy"

    def query(self, query: List[str]) -> str:
        cmd = [
            self.efetch,
            "-db",
            self.db,
            "-id",
            ",".join(query),
            "-format",
            "docsum",
            "|",
            self.xtract,
            "-pattern",
            "DocumentSummary",
            "-element",
            "TaxId,ScientificName",
        ]

        return " ".join(cmd)

    def read_output(self, output_path: str) -> pd.DataFrame:
        """
        Read output from Entrez query using pandas
        """

        df = pd.read_csv(
            output_path, sep="\t", header=None, names=["taxid", "scientific_name"]
        )

        return df


class EntrezQueryFactory:
    def __init__(self, bindir: str):
        self.bindir = bindir

    def get_query(self, name: str) -> EntrezQuery:
        if name == "fetch_taxid_description":
            return EntrezFetchTaxidDescription(self.bindir)
        else:
            raise ValueError("Invalid query name")


class EntrezWrapper:
    bindir: str

    def __init__(
        self,
        bindir: str,
        outdir: str,
        outfile: str,
        query_type: str = "fetch_taxid_description",
        chunksize: int = 100,
    ):
        self.chunksize = chunksize
        self.bindir = bindir
        self.outfile = outfile
        self.outdir = outdir

        self.query_factory = EntrezQueryFactory(self.bindir)
        self.query = self.query_factory.get_query(query_type)

    @property
    def output_path(self) -> str:
        return os.path.join(self.outdir, self.outfile)

    def cmd_long(self, query: List[str]) -> str:
        cmd = self.query.query(query)

        cmd_long = [cmd, ">>", os.path.join(self.outdir, self.outfile)]

        print(" ".join(cmd_long))

        return " ".join(cmd_long)

    def cmd_chunks(self, query: List[str]) -> List[str]:
        chunks = [
            query[i : i + self.chunksize] for i in range(0, len(query), self.chunksize)
        ]

        return [self.cmd_long(chunk) for chunk in chunks]

    def run_queries(self, query: List[str]) -> None:
        cmds = self.cmd_chunks(query)

        print("CMD CHUNKS")

        for cmd in cmds:
            os.system(cmd)

        return None

    def read_output(self) -> pd.DataFrame:
        output_path = os.path.join(self.outdir, self.outfile)
        return self.query.read_output(output_path)

    def run(self, query: List[str]) -> pd.DataFrame:
        self.run_queries(query)
        df = self.read_output()

        return df
