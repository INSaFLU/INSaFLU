import os
from abc import ABC, abstractmethod
from typing import List

import pandas as pd
from Bio import Entrez
from django.contrib.auth.models import User


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
        username: str,
        bindir: str,
        outdir: str,
        outfile: str,
        query_type: str = "fetch_taxid_description",
        chunksize: int = 400,
    ):
        self.chunksize = chunksize
        self.bindir = bindir
        self.outfile = outfile
        self.outdir = outdir

        user = User.objects.get(username=username)

        Entrez.email = user.email
        Entrez.max_tries = 1
        Entrez.sleep_between_tries = 1

        self.query_factory = EntrezQueryFactory(self.bindir)
        self.query = self.query_factory.get_query(query_type)

    @property
    def output_path(self) -> str:
        return os.path.join(self.outdir, self.outfile)

    def cmd_long(self, query: List[str]) -> str:
        cmd = self.query.query(query)

        cmd_long = [cmd, ">>", os.path.join(self.outdir, self.outfile)]

        return " ".join(cmd_long)

    def split_query(self, query: List[str]) -> List[List[str]]:
        return [
            query[i : i + self.chunksize] for i in range(0, len(query), self.chunksize)
        ]

    def cmd_chunks(self, query: List[str]) -> List[str]:
        chunks = self.split_query(query)

        return [self.cmd_long(chunk) for chunk in chunks]

    def run_queries_binaries(self, query: List[str]) -> None:
        """
        run queries using entrez direct binaries"""
        cmds = self.cmd_chunks(query)

        if os.path.exists(self.output_path):
            os.remove(self.output_path)

        for cmd in cmds:
            os.system(cmd)

        output_path = os.path.join(self.outdir, self.outfile)
        df = pd.read_csv(output_path, sep="\t", header=None)
        df.columns = ["taxid", "scientific_name"]
        df.to_csv(self.output_path, sep="\t", index=False)

        return None

    def run_queries_biopy(self, query: List[str]) -> None:
        """
        run queries using Biopython Entrez
        """

        chunks = self.split_query(query)

        report = []

        for chunk in chunks:
            handle = Entrez.efetch(db="Taxonomy", id=",".join(chunk), retmode="xml")
            record = Entrez.read(handle)
            records = [[record["TaxId"], record["ScientificName"]] for record in record]
            report.extend(records)

        df = pd.DataFrame(report, columns=["taxid", "scientific_name"])

        df.to_csv(self.output_path, sep="\t", index=False)

    def read_output(self) -> pd.DataFrame:
        output_path = os.path.join(self.outdir, self.outfile)
        return pd.read_csv(output_path, sep="\t")

    def run(self, query: List[str]) -> pd.DataFrame:
        self.run_queries_biopy(query)
        df = self.read_output()

        return df
