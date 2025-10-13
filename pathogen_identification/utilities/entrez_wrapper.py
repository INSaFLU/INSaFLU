import http.client
import os
import urllib.error
from abc import ABC, abstractmethod
from typing import List

import pandas as pd
from Bio import Entrez
from django.contrib.auth.models import User


def split_query(query: List[str], chunksize: int) -> List[List[str]]:
    return [query[i : i + chunksize] for i in range(0, len(query), chunksize)]


class EntrezQuery(ABC):
    name: str
    output_columns: List[str]

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
    output_columns = ["taxid", "scientific_name"]

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


class EntrezFetchProteinAccession_Taxon(EntrezQuery):
    name: str = "fetch_protein_accession_taxon"
    db = "protein"
    output_columns = ["acc", "taxid"]

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
            "AccessionVersion,TaxId",
        ]

        return " ".join(cmd)

    def read_output(self, output_path: str) -> pd.DataFrame:
        """
        Read output from Entrez query using pandas
        """

        df = pd.read_csv(output_path, sep="\t", header=None, names=["acc", "taxid"])

        return df


class EntrezFetchAccessionDescription(EntrezQuery):
    name: str = "fetch_accession_description"
    db = "nuccore"
    output_columns = ["taxid", "accession", "description"]

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
            "TaxId,AccessionVersion,Title",
        ]

        return " ".join(cmd)

    def read_output(self, output_path: str) -> pd.DataFrame:
        """
        Read output from Entrez query using pandas
        """

        df = pd.read_csv(
            output_path,
            sep="\t",
            header=None,
            names=self.output_columns,
        )
        return df


class EntrezQueryFactory:
    def __init__(self, bindir: str):
        self.bindir = bindir

    def get_query(self, name: str) -> EntrezQuery:
        if name == "fetch_taxid_description":
            return EntrezFetchTaxidDescription(self.bindir)
        elif name == "fetch_accession_description":
            return EntrezFetchAccessionDescription(self.bindir)
        elif name == "fetch_protein_accession_taxon":
            return EntrezFetchProteinAccession_Taxon(self.bindir)
        else:
            raise ValueError("Invalid query name")


###########################################################


class BiopythonEntrezWrapper(ABC):
    output_path: str
    chuncksize: int
    output_column_names: List[str]

    def __init__(self, output_path: str, chuncksize: int = 400):
        self.output_path = output_path
        self.chuncksize = chuncksize

    @abstractmethod
    def run_query(self, query: List[str]) -> pd.DataFrame:
        pass

    def read_output(self) -> pd.DataFrame:
        try:
            return pd.read_csv(self.output_path, sep="\t")
        except FileNotFoundError:
            return pd.DataFrame()


class BiopyFetchTaxidDescription(BiopythonEntrezWrapper):
    def __init__(self, output_path: str, chuncksize: int = 400):
        super().__init__(output_path, chuncksize=chuncksize)
        self.output_column_names = ["taxid", "scientific_name"]

    def run_query(self, query: List[str]) -> None:
        chunks = split_query(query, self.chuncksize)

        report = []

        for chunk in chunks:
            handle = Entrez.efetch(db="Taxonomy", id=",".join(chunk), retmode="xml")
            record = Entrez.read(handle)
            records = [[record["TaxId"], record["ScientificName"]] for record in record]
            report.extend(records)

        df = pd.DataFrame(report, columns=self.output_column_names)

        df.to_csv(self.output_path, sep="\t", index=False)


class BiopyFetchAccessionDescription(BiopythonEntrezWrapper):
    def __init__(self, output_path: str, chuncksize: int = 400):
        super().__init__(output_path, chuncksize=chuncksize)
        self.output_column_names = ["acc", "description"]

    def run_query(self, query: List[str]) -> None:
        chunks = split_query(query, self.chuncksize)

        report = []

        for chunk in chunks:
            chunk = [str(i).split(".")[0] for i in chunk]
            handle = Entrez.efetch(db="Taxonomy", id=",".join(chunk), retmode="xml")
            record = Entrez.read(handle)

            records = [
                [record["AccessionVersion"], record["Title"]] for record in record
            ]
            report.extend(records)

        df = pd.DataFrame(report, columns=self.output_column_names)

        df.to_csv(self.output_path, sep="\t", index=False)


class BiopythonEntrezQueryFactory:
    def __init__(self, output_path: str, chuncksize: int = 400):
        self.output_path = output_path
        self.chuncksize = chuncksize

    def get_query(self, name: str) -> BiopythonEntrezWrapper:
        if name == "fetch_taxid_description":
            return BiopyFetchTaxidDescription(
                self.output_path, chuncksize=self.chuncksize
            )
        elif name == "fetch_accession_description":
            return BiopyFetchAccessionDescription(
                self.output_path, chuncksize=self.chuncksize
            )
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

        self.bin_query_factory = EntrezQueryFactory(self.bindir)
        self.bin_query = self.bin_query_factory.get_query(query_type)

        self.biopy_query_factory = BiopythonEntrezQueryFactory(
            self.output_path, chuncksize=self.chunksize
        )
        self.biopy_query = self.biopy_query_factory.get_query(query_type)

    def run_entrez_query(self, query_list: List[str]) -> pd.DataFrame:
        output = pd.DataFrame()
        try:
            # self.biopy_query.run_query(query_list)
            self.run_queries_binaries(query_list)
            output = self.read_output()

        except (
            urllib.error.URLError,
            http.client.RemoteDisconnected,
            http.client.IncompleteRead,
        ):
            self.run_queries_binaries(query_list)
            output = self.read_output()

        return output

    @property
    def output_path(self) -> str:
        return os.path.join(self.outdir, self.outfile)

    def cmd_long(self, query: List[str]) -> str:
        cmd = self.bin_query.query(query)

        cmd_long = [cmd, ">>", os.path.join(self.outdir, self.outfile)]

        return " ".join(cmd_long)

    def split_query(self, query: List[str]) -> List[List[str]]:
        return split_query(query, self.chunksize)

    def cmd_chunks(self, query: List[str]) -> List[str]:
        chunks = self.split_query(query)

        return [self.cmd_long(chunk) for chunk in chunks]

    def run_query_strategies(self, query: List[str]) -> None:
        try:
            self.run_taxid_description_queries_biopy(query)
        except (
            urllib.error.URLError,
            http.client.RemoteDisconnected,
            http.client.IncompleteRead,
        ):
            self.run_queries_binaries(query)
        except http.client.RemoteDisconnected:
            self.run_queries_binaries(query)
        except pd.errors.EmptyDataError:
            pass

    def entrez_get_taxid_descriptions(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Get taxid descriptions from entrez.
        """
        assert "taxid" in df.columns

        taxid_df = df.dropna(subset=["taxid"])
        taxid_list = taxid_df.taxid.unique().tolist()
        taxid_list = [str(int(i)) for i in taxid_list]

        if len(taxid_list) == 0:
            return pd.DataFrame(columns=["taxid", "counts", "description"])

        self.run_query_strategies(taxid_list)

        taxid_descriptions = self.read_output()

        if taxid_descriptions.shape[0] == 0:
            df["description"] = ""
            return df

        taxid_descriptions.rename(
            columns={"scientific_name": "description"}, inplace=True
        )

        df["taxid"] = df["taxid"].astype(int)
        taxid_descriptions["taxid"] = taxid_descriptions["taxid"].astype(int)

        df = df.merge(taxid_descriptions, on="taxid", how="left")

        df["taxid"] = df["taxid"].astype(float)
        df["taxid"] = df["taxid"].astype(int)

        return df

    def process_query_output(self, output_path: str) -> pd.DataFrame:
        """
        Process the output of the query. some rows have the taxid column repeated, ending wwith 4 columns instead of 3
        """

        tmp_duplicate_file = os.path.join(self.outdir, "tmp_duplicate_taxids.txt")
        tmp_file = os.path.join(self.outdir, "tmp_file.txt")
        os.system(f"awk -F'\t' 'NF==4' {output_path} > {tmp_duplicate_file}")
        os.system(f"awk -F'\t' 'NF==3' {output_path} > {tmp_file}")
        os.system("cut -f2,3,4 " + tmp_duplicate_file + " >> " + tmp_file)
        os.system("mv " + tmp_file + " " + output_path)
        os.system("rm " + tmp_duplicate_file)
        os.system("rm " + tmp_file)

    def filter_query_output(self, output_path: str) -> pd.DataFrame:
        tmp_file = os.path.join(self.outdir, "tmp_file.txt")
        os.system(f"awk -F'\t' 'NF==3' {output_path} > {tmp_file}")
        os.system("mv " + tmp_file + " " + output_path)

    def run_queries_binaries(self, query: List[str]) -> None:
        """
        run queries using entrez direct binaries"""
        cmds = self.cmd_chunks(query)

        if os.path.exists(self.output_path):
            os.remove(self.output_path)

        for cmd in cmds:
            os.system(cmd)

        output_path = os.path.join(self.outdir, self.outfile)
        self.process_query_output(output_path)
        try:
            df = pd.read_csv(output_path, sep="\t", header=None)
            df.columns = self.bin_query.output_columns
        except pd.errors.ParserError:
            self.filter_query_output(output_path)
            df = pd.read_csv(output_path, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            df = pd.DataFrame(columns=self.bin_query.output_columns)

        df.to_csv(self.output_path, sep="\t", index=False)

        return None

    def run_taxid_description_queries_biopy(self, query: List[str]) -> None:
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
        try:
            return pd.read_csv(output_path, sep="\t")
        except FileNotFoundError:
            return pd.DataFrame()

    def run(self, query: List[str]) -> pd.DataFrame:
        self.run_taxid_description_queries_biopy(query)
        df = self.read_output()

        return df
