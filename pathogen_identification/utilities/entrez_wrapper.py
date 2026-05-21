import http.client
import os
import urllib.error
from abc import ABC, abstractmethod
from typing import List, Dict, Tuple
from Bio import Entrez
from django.contrib.auth.models import User
from dataclasses import dataclass
from typing import Optional


def split_query(query: List[str], chunksize: int) -> List[List[str]]:
    return [query[i : i + chunksize] for i in range(0, len(query), chunksize)]

@dataclass
class LineageNode:
    taxid: str = ""
    name: str = ""
    rank: str = "no rank"


class EntrezQuery(ABC):
    name: str
    output_columns: List[str]

    def __init__(self, bindir: str):
        self.bindir = bindir

    @abstractmethod
    def query(self, query: List[str]) -> str:
        pass

    @abstractmethod
    def read_output(self, output_path: str) -> pd.DataFrame:
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

    def process_query_output(self, output_path: str) -> None:
        """
        Process the output of the query. some rows have the taxid column repeated, ending wwith 4 columns instead of 3
        """
        outdir = os.path.dirname(output_path)
        tmp_duplicate_file = os.path.join(outdir, "tmp_duplicate_taxids.txt")
        tmp_file = os.path.join(outdir, "tmp_file.txt")
        os.system(f"awk -F'\t' 'NF==4' {output_path} > {tmp_duplicate_file}")
        os.system(f"awk -F'\t' 'NF==3' {output_path} > {tmp_file}")
        os.system("cut -f2,3,4 " + tmp_duplicate_file + " >> " + tmp_file)
        os.system("mv " + tmp_file + " " + output_path)
        os.system("rm " + tmp_duplicate_file)
        os.system("rm " + tmp_file)



    def read_output(self, output_path: str) -> pd.DataFrame:
        """
        Read output from Entrez query using pandas
        """

        self.process_query_output(output_path)

        df = pd.read_csv(
            output_path,
            sep="\t",
            header=None,
            names=self.output_columns,
        )
        return df


class EntrezFetchTaxidLineage(EntrezQuery):
    """Fetch taxonomic lineage for TaxIDs via NCBI Taxonomy"""
    name: str = "fetch_taxid_lineage"
    db = "taxonomy"
    output_columns = ["taxid", "rank", "name"]

    def query(self, query: List[str]) -> str:
        # This is a Biopython-based query, not a binary command
        # Placeholder for consistency with interface
        return ""

    def read_output(self, output_path: str) -> pd.DataFrame:
        """Return empty DataFrame - this is handled by fetch_lineage method"""
        return pd.DataFrame()


class EntrezSearchOrganismByName(EntrezQuery):
    """Search NCBI Taxonomy database by organism name"""
    name: str = "search_organism_by_name"
    db = "taxonomy"
    output_columns = ["input_name", "taxid", "canonical_name", "rank", "confidence"]

    def query(self, query: List[str]) -> str:
        # This is a Biopython-based query, not a binary command
        # Placeholder for consistency with interface
        return ""

    def read_output(self, output_path: str) -> pd.DataFrame:
        """Return empty DataFrame - this is handled by search_organism_name method"""
        return pd.DataFrame()

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
        elif name == "fetch_taxid_lineage":
            return EntrezFetchTaxidLineage(self.bindir)
        elif name == "search_organism_by_name":
            return EntrezSearchOrganismByName(self.bindir)
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



    def run_queries_binaries(self, query: List[str]) -> None:
        """
        run queries using entrez direct binaries"""
        cmds = self.cmd_chunks(query)

        if os.path.exists(self.output_path):
            os.remove(self.output_path)

        for cmd in cmds:
            print(cmd)
            os.system(cmd)

        output_path = os.path.join(self.outdir, self.outfile)
        with open(output_path, "r") as f:
            print("########### read")
            print(f.read())
        import traceback
        try:
            df = self.bin_query.read_output(output_path)
            df.columns = self.bin_query.output_columns

        except pd.errors.EmptyDataError:
            traceback.print_exc()
            df = pd.DataFrame(columns=self.bin_query.output_columns)

        df.to_csv(self.output_path, sep="\t", index=False)

    # NEW METHODS FOR LINEAGE & NAME RESOLUTION

    def fetch_lineage(self, taxids: List[str], use_cache: bool = True) -> Dict[str, List[Dict]]:
        """
        Fetch taxonomic lineages for multiple taxids.

        Args:
            taxids: List of NCBI taxonomy IDs (as strings)
            use_cache: Whether to use cached lineages from database

        Returns:
            Dict mapping taxid → list of {rank, name, taxid} dicts
        """
        lineages = {}
        
        # Chunk taxids to respect NCBI rate limits
        chunks = split_query(taxids, self.chunksize)
        
        for chunk in chunks:
            try:
                # Use Biopython to fetch lineage
                handle = Entrez.efetch(db="Taxonomy", id=",".join(chunk), retmode="xml")
                records = Entrez.read(handle)
                
                # Extract lineage from each record
                for record in records:
                    taxid = str(record.get("TaxId", ""))
                    lineage_nodes = []
                    
                    # Parse LineageEx if present
                    if "LineageEx" in record:
                        for taxon in record["LineageEx"]:
                            lineage_nodes.append(LineageNode(
                                taxid=str(taxon.get("TaxId", "")),
                                name=taxon.get("ScientificName", ""),
                                rank=taxon.get("Rank", "no rank")
                            ))
                    
                    # Add the record itself as the leaf node
                    lineage_nodes.append(LineageNode(
                        taxid=taxid,
                        name=record.get("ScientificName", ""),
                        rank=record.get("Rank", "no rank")
                    ))
                    
                    lineages[taxid] = lineage_nodes
            except Exception as e:
                print(f"Error fetching lineage for chunk {chunk}: {e}")
                continue
        
        return lineages

    def persist_lineages(self, lineages: Dict[str, List[LineageNode]]) -> None:
        """
        Persist taxonomic lineages as Taxon objects in DB.

        Builds a deduplicated node registry and adjacency list, then BFS from
        roots to create/update Taxon objects. Updates ReferenceTaxid rank-level FKs.

        Args:
            lineages: taxid -> list of LineageNode from root to leaf
        """
        from collections import deque
        from pathogen_identification.models import Taxon, ReferenceTaxid
        from constants.constants_taxonomy import TaxonConstants

        RANK_TO_FIELD = {
            TaxonConstants.RANK_DOMAIN: "tax_domain",
            TaxonConstants.RANK_PHYLUM: "tax_phylum",
            TaxonConstants.RANK_CLASS: "tax_class",
            TaxonConstants.RANK_ORDER: "tax_order",
            TaxonConstants.RANK_FAMILY: "tax_family",
            TaxonConstants.RANK_GENUS: "tax_genus",
        }

        node_info = {}
        children = {}
        all_child_taxids = set()

        for leaf_taxid, lineage_list in lineages.items():
            prev = None
            for node in lineage_list:
                tid = node.taxid
                if not tid:
                    continue
                if tid not in node_info:
                    node_info[tid] = node
                if prev is not None:
                    children.setdefault(prev, set()).add(tid)
                    all_child_taxids.add(tid)
                prev = tid

        roots = sorted(set(node_info) - all_child_taxids)
        queue = deque()
        taxon_map = {}

        for root_tid in roots:
            queue.append((root_tid, None))

        while queue:
            tid, parent_taxon = queue.popleft()
            node = node_info[tid]

            taxon, _ = Taxon.objects.get_or_create(
                taxid=int(tid),
                defaults={
                    "name": node.name,
                    "rank": node.rank,
                    "parent": parent_taxon,
                },
            )

            changed = False
            if taxon.name != node.name:
                taxon.name = node.name
                changed = True
            if taxon.rank != node.rank:
                taxon.rank = node.rank
                changed = True
            if taxon.parent != parent_taxon:
                taxon.parent = parent_taxon
                changed = True
            if changed:
                taxon.save()

            taxon_map[tid] = taxon

            for child_tid in children.get(tid, set()):
                queue.append((child_tid, taxon))

        for leaf_taxid, lineage_list in lineages.items():
            if not leaf_taxid:
                continue
            try:
                ref_taxid_obj = ReferenceTaxid.objects.get(taxid=leaf_taxid)
            except ReferenceTaxid.DoesNotExist:
                continue

            changed = False
            for node in lineage_list:
                tid = node.taxid
                if not tid:
                    continue
                normalized_rank = TaxonConstants.normalize_rank(node.rank)
                field = RANK_TO_FIELD.get(normalized_rank)
                if field is not None:
                    taxon = taxon_map.get(tid)
                    if taxon is not None and getattr(ref_taxid_obj, field) != taxon:
                        setattr(ref_taxid_obj, field, taxon)
                        changed = True
            if changed:
                ref_taxid_obj.save()

    def search_organism_name(self, names: List[str], use_fuzzy: bool = True) -> Dict[str, Dict]:
        """
        Search for organisms by name. Tries exact match first, then NCBI search.

        Args:
            names: List of organism names (can be messy, abbreviated, etc.)
            use_fuzzy: Enable fuzzy matching for partial matches

        Returns:
            Dict mapping input_name → {taxid, canonical_name, confidence, source}
        """
        results = {}
        
        for name in names:
            # Try direct NCBI search
            try:
                handle = Entrez.esearch(db="Taxonomy", term=name, retmax=1)
                search_result = Entrez.read(handle)
                
                if search_result["IdList"]:
                    taxid = search_result["IdList"][0]
                    
                    # Fetch details for this taxid
                    handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
                    fetch_result = Entrez.read(handle)
                    
                    if fetch_result:
                        record = fetch_result[0]
                        results[name] = {
                            "taxid": taxid,
                            "canonical_name": record.get("ScientificName", ""),
                            "rank": record.get("Rank", ""),
                            "confidence": "high",
                            "source": "ncbi"
                        }
                    else:
                        results[name] = {
                            "taxid": None,
                            "canonical_name": None,
                            "rank": None,
                            "confidence": "low",
                            "source": "failed"
                        }
                else:
                    results[name] = {
                        "taxid": None,
                        "canonical_name": None,
                        "rank": None,
                        "confidence": "none",
                        "source": "not_found"
                    }
            except Exception as e:
                print(f"Error searching for organism name '{name}': {e}")
                results[name] = {
                    "taxid": None,
                    "canonical_name": None,
                    "rank": None,
                    "confidence": "error",
                    "source": "exception"
                }
        
        return results

    def enrich_references_dataframe(
        self,
        df: pd.DataFrame,
        lineages: Dict[str, List[LineageNode]] | None = None,
    ) -> pd.DataFrame:
        """
        Add lineage path to reference DataFrame.

        Assumes input df has columns: taxid, accession, description

        Args:
            df: Input DataFrame
            lineages: Optional pre-fetched lineages (avoids double fetch)

        Returns:
            DataFrame with additional column:
            - lineage_path (human-readable)
        """
        if df.empty:
            return df

        if lineages is None:
            taxids = df["taxid"].astype(str).unique().tolist()
            lineages = self.fetch_lineage(taxids)

        lineage_paths = {}

        for taxid_str, lineage_list in lineages.items():
            names = [node.name for node in lineage_list]
            lineage_paths[taxid_str] = " > ".join(filter(None, names))

        df["lineage_path"] = df["taxid"].astype(str).map(
            lambda x: lineage_paths.get(x, "")
        )

        return df

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
