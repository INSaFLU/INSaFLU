# You need to install Biopython: pip install biopython
import json
import logging
import os
import subprocess
import time
from dataclasses import dataclass
from functools import wraps
from typing import Generator, List, Optional, Tuple

import dotenv
import numpy as np
import pandas as pd
from Bio import Entrez

dotenv.load_dotenv()
from decouple import config

Entrez.email = config("NCBI_EMAIL")
if Entrez.email is None:
    raise ValueError("NCBI_EMAIL environment variable not set. Please set it to your email address.")

NCBI_TAXONOMY_LEVELS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
NCBI_TAXONOMY_LEVELS_EXTENDED = ['acellular root', 'realm', 'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']


def retry_with_backoff(max_retries=3, initial_delay=1, backoff_factor=2):
    """
    Decorator for retrying functions with exponential backoff.
    Handles transient errors like rate limiting (429) and temporary network issues.
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            delay = initial_delay
            last_exception = None
            for attempt in range(max_retries):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    last_exception = e
                    error_msg = str(e).lower()
                    should_retry = any(x in error_msg for x in [
                        '429', 'rate limit', 'temporary failure',
                        'connection', 'timeout', 'service unavailable',
                        '500', '502', '503', '504'
                    ])
                    if not should_retry and attempt > 0:
                        should_retry = attempt < 2

                    if should_retry and attempt < max_retries - 1:
                        logging.warning(f"Attempt {attempt + 1}/{max_retries} failed for {func.__name__}: {e}. Retrying in {delay}s...")
                        time.sleep(delay)
                        delay *= backoff_factor
                    elif attempt == max_retries - 1:
                        logging.error(f"All {max_retries} attempts failed for {func.__name__}: {e}")
            raise last_exception
        return wrapper
    return decorator



@dataclass
class Passport:
    taxid: Optional[str]
    accession: Optional[str] = None
    lineage: Optional[str] = None
    description: Optional[str] = None

    def __init__(self, taxid: Optional[str], accession: Optional[str] = None, lineage: Optional[str] = None, description: Optional[str] = None):
        if taxid and "." in str(taxid):
            taxid = str(taxid).split(".")[0]
        self.taxid = taxid
        self.accession = accession
        self.lineage = lineage
        self.description = description

    def __str__(self):
        return f"TaxID: {self.taxid}, Accession: {self.accession}"

    @property
    def prefix(self):
        if self.accession:
            return f"{self.taxid}_{self.accession}"
        else:
            return f"{self.taxid}"


@dataclass
class LocalAssembly(Passport):
    file_path: Optional[str] = None

    def __str__(self):
        return f"TaxID: {self.taxid}, Accession: {self.accession}, File Path: {self.file_path}"


@dataclass
class ReferenceData(Passport):
    nucleotide_id: Optional[str] = None
    assembly_id: Optional[str] = None

    def __str__(self):
        return f"TaxID: {self.taxid}, Accession: {self.accession}, Description: {self.description}, Nucleotide ID: {self.nucleotide_id}, Assembly ID: {self.assembly_id}"





@retry_with_backoff(max_retries=3, initial_delay=1)
def retrieve_reference_sequence_id(accID: str, include_term=None, exclude_term=None) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    try:
        term = accID
        if exclude_term is not None:
            term += f" NOT {exclude_term}"
        if include_term is not None:
            term += f" AND {include_term}"

        handle = Entrez.esearch(db="nucleotide", term=term, retmax=20)
        record = Entrez.read(handle)
        handle.close()

        if not record['IdList']:
            raise ValueError(f"No sequence found for accession {accID}")

        sequence_id = record['IdList'][0]

        handle = Entrez.esummary(db="nucleotide", id=sequence_id)
        summary = Entrez.read(handle)
        handle.close()
        if summary is None:
            raise ValueError(f"No summary found for sequence ID {sequence_id}")
        accession = summary[0]['AccessionVersion']
        if accession != accID:
            print(f"Warning: Retrieved accession {accession} does not match requested {accID}")
        description = summary[0].get('Title', 'No description available')
        return accession, description, sequence_id

    except ValueError as e:
        print(e)
        return None, None, None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None, None, None


@retry_with_backoff(max_retries=3, initial_delay=1)
def get_reference_sequence_url(taxid, include_term=None, exclude_term=None) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    try:
        term = f"txid{taxid}[Organism:exp] AND refseq"
        if include_term is not None:
            term += f" AND {include_term}"
        if exclude_term is not None:
            term += f" NOT {exclude_term}"

        handle = Entrez.esearch(db="nucleotide", term=term, retmax=1)
        record = Entrez.read(handle)
        handle.close()
        if not record['IdList']:
            raise ValueError(f"No reference sequences found for taxid {taxid}")

        sequence_id = record['IdList'][0]

        handle = Entrez.esummary(db="nucleotide", id=sequence_id)
        summary = Entrez.read(handle)
        handle.close()

        if summary is None:
            raise ValueError(f"No summary found for sequence ID {sequence_id}")

        docsum = summary[0]
        accession = docsum['AccessionVersion']
        description = docsum.get('Title', 'No description available')

        return accession, description, sequence_id

    except ValueError as e:
        print(e)
        return None, None, None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None, None, None


@retry_with_backoff(max_retries=3, initial_delay=1)
def get_representative_assembly(taxid, include_term=None, exclude_term=None) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    try:
        term = f"txid{taxid}[Organism:exp]"
        if include_term is not None:
            term += f" AND {include_term}"
        if exclude_term is not None:
            term += f" NOT {exclude_term}"

        handle = Entrez.esearch(db="assembly", term=term, retmax=5)
        record = Entrez.read(handle)
        handle.close()
        if not record['IdList']:
            raise ValueError(f"No representative genomes found for taxid {taxid}")

        assembly_id = record['IdList'][0]

        handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
        summary = Entrez.read(handle)
        handle.close()
        docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
        accession = docsum['AssemblyAccession']
        description = docsum.get('SpeciesName', 'No description available')

        return accession, description, assembly_id

    except ValueError as e:
        print(e)
        return None, None, None
    except Exception as e:
        print(f"An error occurred while fetching assembly for taxid {taxid}: {e}")
        return None, None, None

@retry_with_backoff(max_retries=3, initial_delay=1)
def retrieve_reference_sequence(nucleotide_id, output_path, gzipped=True) -> bool:
    try:
        handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()
        if gzipped:
            import gzip
            with gzip.open(output_path, 'wt') as f:
                f.write(fasta_data)
        else:
            with open(output_path, 'w') as f:
                f.write(fasta_data)
        return True
    except Exception as e:
        print(f"An error occurred while downloading sequence: {e}")
        return False

@retry_with_backoff(max_retries=3, initial_delay=1)
def retrieve_assembly_sequence(assembly_id, output_path) -> bool:
    try:
        handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
        summary = Entrez.read(handle)
        handle.close()
        docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
        ftp_path = docsum['FtpPath_RefSeq'] or docsum['FtpPath_GenBank']
        if not ftp_path:
            print(f"No FTP path found for assembly ID {assembly_id}")
            return False
        asm_name = ftp_path.split('/')[-1]
        fasta_url = f"{ftp_path}/{asm_name}_genomic.fna.gz"
        result = subprocess.run(['wget', '-O', output_path, fasta_url], capture_output=True, check=False)
        if result.returncode != 0:
            raise RuntimeError(f"Failed to download file: {result.stderr.decode()}")
        return True
    except Exception as e:
        print(f"An error occurred while downloading assembly: {e}")
        return False


class NCBITools:
    def __init__(self):
        self.logger = logging.getLogger('NCBITools')
        self.logger.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)
        self.logger.propagate = False

    @retry_with_backoff(max_retries=3, initial_delay=1)
    def retrieve_passport_taxonomy(self, passport: Passport) -> Optional[str]:
        try:
            handle = Entrez.efetch(db="taxonomy", id=passport.taxid, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            if not records:
                raise ValueError(f"No taxonomy found for taxid {passport.taxid}")
            lineage = records[0]['Lineage']
            return lineage
        except Exception as e:
            self.logger.error(f"An error occurred while fetching taxonomy for taxid {passport.taxid}: {e}")
            return None
        
    @retry_with_backoff(max_retries=3, initial_delay=1)
    def retrieve_taxonomies_batch(self, taxids: list[int]) -> dict:
        taxid_list = ",".join(map(str, taxids))
        fetch_cmd = [
            "efetch",
            "-db",
            "taxonomy",
            "-id", taxid_list,
            "-format", "xml",
            "|", "xtract", "-pattern", "Taxon", "-element", "TaxId", "Lineage"
        ]
        try:
            result = subprocess.run(" ".join(fetch_cmd), shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                self.logger.error(f"Failed to retrieve taxonomies for taxids {taxids}: {result.stderr}")
                return {}
            taxonomies = {}
            for line in result.stdout.splitlines():
                line = line.split()
                taxids = []
                lineage = ""
                for part in line:
                    if part.isdigit():
                        taxids.append(int(part))
                    else:
                        lineage += part + " "

                if taxids:
                    taxonomies.update(dict.fromkeys(taxids, lineage.strip()))

            return taxonomies

        except Exception as e:
            import traceback
            
            self.logger.error(f"An error occurred while retrieving taxonomies for taxids {taxids}: {e}")
            self.logger.error(traceback.format_exc())
            return {}

    def query_sequence_databases(self, passport: Passport, include_term: Optional[str] = None, exclude_term: Optional[str] = None) -> ReferenceData:
        if passport.taxid is None or pd.isna(passport.taxid) or np.isnan(float(passport.taxid)):
            self.logger.error("No taxid provided in passport.")
            return ReferenceData(taxid=None, accession=None, lineage=None, description=None)

        lineage = self.retrieve_passport_taxonomy(passport)
        self.logger.info(f"Lineage for taxid {passport.taxid}: {lineage}")
        passport.lineage = lineage

        if passport.accession is not None:
            accession, description, nucleotide_id = retrieve_reference_sequence_id(passport.accession, include_term, exclude_term)
            if nucleotide_id is not None:
                return ReferenceData(
                    taxid=passport.taxid,
                    accession=passport.accession,
                    lineage=passport.lineage,
                    description=description,
                    nucleotide_id=nucleotide_id
                )
            else:
                self.logger.warning(f"No nucleotide ID found for accession {passport.accession}, falling back to taxid search.")

        accession, description, nucleotide_id = get_reference_sequence_url(passport.taxid, include_term, exclude_term)

        if accession is not None and nucleotide_id is not None:
            return ReferenceData(
                taxid=passport.taxid,
                accession=accession,
                lineage=passport.lineage,
                description=description,
                nucleotide_id=nucleotide_id
            )

        accession, description, assembly_id = get_representative_assembly(passport.taxid, include_term, exclude_term)

        return ReferenceData(
            taxid=passport.taxid,
            accession=accession,
            lineage=passport.lineage,
            description=description,
            assembly_id=assembly_id
        )

    def retrieve_sequence_databases(self, reference_data: ReferenceData, output_path: str, gzipped: bool = True) -> bool:
        if reference_data.nucleotide_id is not None:
            success = retrieve_reference_sequence(reference_data.nucleotide_id, output_path, gzipped)
            if not success:
                self.logger.warning(f"Failed to retrieve reference sequence for taxid {reference_data.taxid}")
            return success

        if reference_data.assembly_id is not None:
            success = retrieve_assembly_sequence(reference_data.assembly_id, output_path)
            if not success:
                self.logger.warning(f"Failed to retrieve assembly sequence for taxid {reference_data.taxid}")
            return success

        self.logger.warning(f"No sequence data found for taxid {reference_data.taxid}")
        return False
