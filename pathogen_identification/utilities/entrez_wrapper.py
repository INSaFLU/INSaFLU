import http.client
import os
import urllib.error
from abc import ABC, abstractmethod
from typing import List, Dict, Tuple, TYPE_CHECKING
from Bio import Entrez
from django.contrib.auth.models import User
from dataclasses import dataclass
from typing import Optional
import pandas as pd

if TYPE_CHECKING:
    from pathogen_identification.models import Taxon


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

        import traceback
        try:
            df = self.bin_query.read_output(output_path)
            df.columns = self.bin_query.output_columns

        except pd.errors.EmptyDataError:
            traceback.print_exc()
            df = pd.DataFrame(columns=self.bin_query.output_columns)

        df.to_csv(self.output_path, sep="\t", index=False)

    # NEW METHODS FOR LINEAGE & NAME RESOLUTION

    def fetch_lineage(self, taxid: str) -> List[LineageNode]:
        """
        Fetch taxonomic lineage for a single taxid.

        Args:
            taxid: Single NCBI taxonomy ID (as string)

        Returns:
            List of LineageNode from root to leaf
        """
        try:
            handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            
            if not records:
                return []
            
            record = records[0]
            lineage_nodes = []
            
            # Parse LineageEx if present (ancestors)
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
            
            return lineage_nodes
        except Exception as e:
            print(f"Error fetching lineage for taxid {taxid}: {e}")
            return []

    def fetch_lineages(self, taxids: List[str], strategy: str = "biopy") -> Dict[str, List[LineageNode]]:
        """
        Fetch taxonomic lineages for multiple taxids using specified strategy.

        Args:
            taxids: List of NCBI taxonomy IDs (as strings)
            strategy: "biopy" (Biopython batch - recommended) or "binary" (NCBI binaries)

        Returns:
            Dict mapping taxid -> list of LineageNode from root to leaf
        """
        if strategy == "binary":
            return self.fetch_lineages_binary(taxids)
        else:
            return self.fetch_lineages_biopy(taxids)

    def fetch_lineages_biopy(self, taxids: List[str]) -> Dict[str, List[LineageNode]]:
        """
        Fetch taxonomic lineages using Biopython (recommended).
        
        Batches requests to minimize API calls while respecting rate limits.
        Single Entrez.efetch call with comma-separated IDs is more efficient
        than per-taxid fetching.

        Args:
            taxids: List of NCBI taxonomy IDs (as strings)

        Returns:
            Dict mapping taxid -> list of LineageNode from root to leaf
        """
        lineages = {}
        
        # Chunk taxids to respect NCBI rate limits
        chunks = split_query(taxids, self.chunksize)
        
        for chunk in chunks:
            try:
                # Single batch call fetches multiple taxids at once
                handle = Entrez.efetch(
                    db="Taxonomy",
                    id=",".join(chunk),
                    retmode="xml"
                )
                records = Entrez.read(handle)
                
                # Extract lineage from each record
                for record in records:
                    taxid = str(record.get("TaxId", ""))
                    lineage_nodes = []
                    
                    # Parse LineageEx if present (ancestors)
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
                print(f"Error fetching lineage for chunk {chunk} via Biopython: {e}")
                continue
        
        return lineages

    def fetch_lineages_binary(self, taxids: List[str]) -> Dict[str, List[LineageNode]]:
        """
        Fetch taxonomic lineages using NCBI EDirect binaries.
        
        Uses binary utilities (esearch, efetch, xtract). More trustworthy
        but requires EDirect to be installed.

        Args:
            taxids: List of NCBI taxonomy IDs (as strings)

        Returns:
            Dict mapping taxid -> list of LineageNode from root to leaf
        """
        import subprocess
        
        lineages = {}
        
        try:
            # Build command using NCBI binaries
            taxid_list = ",".join(taxids)
            cmd = (
                f"echo '{taxid_list}' | "
                f"efetch -db taxonomy -format xml | "
                f"xtract -pattern 'Taxon' -element TaxId,ScientificName,Rank "
                f"    -pattern 'LineageEx' -element LineageEx/Taxon -block 'Taxon' "
                f"    -element TaxId,ScientificName,Rank"
            )
            
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                timeout=60
            )
            
            if result.returncode != 0:
                print(f"Error running binary command: {result.stderr}")
                print("Falling back to Biopython strategy")
                return self.fetch_lineages_biopy(taxids)
            
            # Parse output
            for line in result.stdout.strip().split("\n"):
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) >= 3:
                    taxid, name, rank = parts[0], parts[1], parts[2]
                    if taxid not in lineages:
                        lineages[taxid] = []
                    lineages[taxid].append(LineageNode(
                        taxid=taxid,
                        name=name,
                        rank=rank
                    ))
        except Exception as e:
            print(f"Error fetching lineage via binary: {e}")
            print("Falling back to Biopython strategy")
            return self.fetch_lineages_biopy(taxids)
        
        return lineages

    def persist_lineages(self, lineages: Dict[str, List[LineageNode]]) -> Dict[str, "Taxon"]:
        """
        Persist taxonomic lineages as Taxon objects in DB.

        Builds a deduplicated node registry and adjacency list, then BFS from
        roots to create/update Taxon objects.

        Args:
            lineages: taxid -> list of LineageNode from root to leaf

        Returns:
            Dict mapping taxid -> Taxon object for later use
        """
        from collections import deque
        from pathogen_identification.models import Taxon

        node_info = {}
        children = {}
        all_child_taxids = set()

        for leaf_taxid, lineage_list in lineages.items():
            print(f"  [DEBUG persist] lineage for taxid {leaf_taxid}: {len(lineage_list)} nodes")
            for n in lineage_list:
                print(f"    node: tid={n.taxid!r}, name={n.name!r}, rank={n.rank!r}")
            prev = None
            for node in lineage_list:
                tid = node.taxid
                if not tid:
                    print(f"    SKIP empty tid, prev stays {prev!r}")
                    continue
                if tid not in node_info:
                    node_info[tid] = node
                if prev is not None:
                    children.setdefault(prev, set()).add(tid)
                    all_child_taxids.add(tid)
                    print(f"    children[{prev}] += {{{tid}}}")
                prev = tid

        roots = sorted(set(node_info) - all_child_taxids)
        print(f"  [DEBUG persist] node_info keys: {sorted(node_info.keys())}")
        print(f"  [DEBUG persist] all_child_taxids: {sorted(all_child_taxids)}")
        print(f"  [DEBUG persist] roots: {roots}")
        queue = deque()
        taxon_map = {}

        for root_tid in roots:
            queue.append((root_tid, None))

        while queue:
            tid, parent_taxon = queue.popleft()
            node = node_info[tid]
            print(f"  [DEBUG persist] BFS create: tid={tid}, name={node.name!r}, rank={node.rank!r}, parent_taxid={parent_taxon.taxid if parent_taxon else None}")

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
                print(f"  [DEBUG persist]   enqueue child: {child_tid}")
                queue.append((child_tid, taxon))

        print(f"  [DEBUG persist] final taxon_map keys: {sorted(taxon_map.keys())}")
        return taxon_map

    def link_referencetaxid_to_lineage(
        self, 
        ref_taxid_str: str, 
        lineage: List[LineageNode],
        taxon_map: Dict[str, "Taxon"]
    ) -> None:
        """
        Update ReferenceTaxid rank-level fields from a lineage.

        Call this after ReferenceTaxid object has been created in the database.

        Args:
            ref_taxid_str: ReferenceTaxid taxid as string
            lineage: List of LineageNode from root to leaf
            taxon_map: Dict of taxid -> Taxon objects from persist_lineages()
        """
        from pathogen_identification.models import ReferenceTaxid
        from constants.constants_taxonomy import TaxonConstants

        RANK_TO_FIELD = {
            TaxonConstants.RANK_DOMAIN: "tax_domain",
            TaxonConstants.RANK_PHYLUM: "tax_phylum",
            TaxonConstants.RANK_CLASS: "tax_class",
            TaxonConstants.RANK_ORDER: "tax_order",
            TaxonConstants.RANK_FAMILY: "tax_family",
            TaxonConstants.RANK_GENUS: "tax_genus",
            TaxonConstants.RANK_SPECIES: "tax_species",
        }

        try:
            ref_taxid_obj = ReferenceTaxid.objects.get(taxid=ref_taxid_str)
        except ReferenceTaxid.DoesNotExist:
            return

        changed = False
        for node in lineage:
            tid = node.taxid
            if not tid:
                print(f"  [DEBUG] SKIP (empty tid): rank={node.rank!r}, name={node.name!r}")
                continue
            normalized_rank = TaxonConstants.normalize_rank(node.rank)
            field = RANK_TO_FIELD.get(normalized_rank)
            in_map = tid in taxon_map
            print(f"  [DEBUG] tid={tid!r}, rank={node.rank!r} -> norm={normalized_rank!r} -> field={field!r}, in_taxon_map={in_map}")
            if field is not None:
                taxon = taxon_map.get(tid)
                if taxon is not None and getattr(ref_taxid_obj, field) != taxon:
                    print(f"  [DEBUG]   -> SETTING {field} to Taxon(taxid={taxon.taxid})")
                    setattr(ref_taxid_obj, field, taxon)
                    changed = True
                else:
                    reason = "taxon is None" if taxon is None else f"field already {getattr(ref_taxid_obj, field)}"
                    print(f"  [DEBUG]   -> SKIP ({reason})")
            else:
                print(f"  [DEBUG]   -> SKIP (no matching field in RANK_TO_FIELD)")
        print(f"  [DEBUG] lineage processed, changed={changed}, final fields: dom={ref_taxid_obj.tax_domain}, phy={ref_taxid_obj.tax_phylum}, cla={ref_taxid_obj.tax_class}, ord={ref_taxid_obj.tax_order}, fam={ref_taxid_obj.tax_family}, gen={ref_taxid_obj.tax_genus}")
        if changed:
            ref_taxid_obj.save()

    def search_organism_name(self, names: List[str], use_fuzzy: bool = True) -> Dict[str, Dict]:
        """
        Search for organisms by name. Tries exact match first, then NCBI search.

        Args:
            names: List of organism names (can be messy, abbreviated, etc.)
            use_fuzzy: Enable fuzzy matching for partial matches

        Returns:
            Dict mapping input_name → {taxid, canonical_name, rank, confidence, source}
        """
        # Define result template with default values
        default_result = {
            "taxid": None,
            "canonical_name": None,
            "rank": None,
            "confidence": "error",
            "source": "failed"
        }
        
        results = {name: default_result.copy() for name in names}
        
        for name in names:
            try:
                # Try direct NCBI search
                handle = Entrez.esearch(db="Taxonomy", term=name, retmax=1)
                search_result = Entrez.read(handle)
                
                if not search_result["IdList"]:
                    results[name]["confidence"] = "none"
                    results[name]["source"] = "not_found"
                    continue
                
                # Found a taxid - fetch details
                taxid = search_result["IdList"][0]
                handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
                fetch_result = Entrez.read(handle)
                
                if fetch_result:
                    record = fetch_result[0]
                    # Update only the fields that were successfully fetched
                    results[name]["taxid"] = taxid
                    results[name]["canonical_name"] = record.get("ScientificName", "")
                    results[name]["rank"] = record.get("Rank", "")
                    results[name]["confidence"] = "high"
                    results[name]["source"] = "ncbi"
                else:
                    results[name]["confidence"] = "low"
                    results[name]["source"] = "failed_fetch"
                    
            except Exception as e:
                print(f"Error searching for organism name '{name}': {e}")
                results[name]["confidence"] = "error"
                results[name]["source"] = "exception"
        
        return results


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
