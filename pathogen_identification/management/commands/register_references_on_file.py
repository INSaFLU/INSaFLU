import os
import time
from abc import ABC, abstractmethod
from datetime import date

import pandas as pd
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from django.db.models import Q

from constants.constants import Televir_Metadata_Constants
from managing_files.models import ProcessControler
from pathogen_identification.models import (
    ReferenceSource,
    ReferenceSourceFile,
    ReferenceSourceFileMap,
    ReferenceTaxid,
)
from pathogen_identification.utilities.entrez_wrapper import EntrezWrapper
from utils.process_SGE import ProcessSGE
from utils.utils import Utils


def return_zgrep(file, pattern, filter=None):

    if filter:
        cmd = f"zgrep {pattern} {file} | {filter}"
    else:
        cmd = f"zgrep {pattern} {file}"

    return cmd


def find_pattern_in_file(file, pattern, filter=None):
    """
    Return result of zgrep pattern file
    """

    cmd = return_zgrep(file, pattern, filter)

    result = os.popen(cmd).read()
    result = result.split("\n")

    return [line for line in result if line]


def find_pattern_multiple_files(files, pattern, filter=None):
    """
    Return result of zgrep pattern file
    """

    cmd = f"zgrep {pattern} {' '.join(files)}"

    result = os.popen(cmd).read()
    result = result.split("\n")
    return [line for line in result if line]


def extract_file_accids(file, output_file, pattern_include="", pattern_exclude=""):
    cmd = f"zgrep {pattern_include} {file} {pattern_exclude} | cut -f1 -d' ' | sort | uniq > {output_file}"
    os.system(cmd)
    # to dict
    with open(output_file, "r") as f:
        accids = f.readlines()
    return {accid.strip().split(":")[0].replace(">", ""): 0 for accid in accids}


class Command(BaseCommand):
    help = "deploy run"

    def add_arguments(self, parser):
        parser.add_argument(
            "--user_id",
            type=int,
            help="user id",
            default=1,
        )

        parser.add_argument(
            "-o",
            "--outdir",
            type=str,
            help="output directory",
        )

        parser.add_argument(
            "--curate",
            action="store_true",
            help="curate references",
        )

    def handle(self, *args, **options):
        ###
        # get user
        user = User.objects.get(pk=options["user_id"])
        utils: Utils = Utils()
        process_controler = ProcessControler()

        process_SGE = ProcessSGE()

        ### SETUP
        print(user.pk)
        process_SGE.set_process_controlers(
            user,
            process_controler.get_name_televir_reference_update(user_pk=user.pk),
            0,
        )

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_reference_update(
                user_pk=user.pk,
            ),
            ProcessControler.FLAG_RUNNING,
        )

        try:
            outdir = utils.get_temp_dir()
            os.makedirs(outdir, exist_ok=True)
            metadadata_constants = Televir_Metadata_Constants()

            # entrez direct interface
            entrez_connection = EntrezWrapper(
                user.username,
                bindir=metadadata_constants.get_software_bin_directory("entrez_direct"),
                outdir=outdir,
                outfile="entrez_output.tsv",
                query_type="fetch_accession_description",
                chunksize=500,
            )

            # accids on file
            accid_file_path = metadadata_constants.accession_to_taxid_path
            accid_file_df = pd.read_csv(accid_file_path, sep="\t")

            print(f"Number of accids on file: {len(accid_file_df)}")

            files = accid_file_df.file.unique().tolist()
            viros_file = [file for file in files if "virosaurus" in file]

            if len(viros_file) == 0:
                ignore_dict = {}
            else:
                viros_file = viros_file[0]

                ignore_dict = extract_file_accids(
                    os.path.join(
                        Televir_Metadata_Constants.SOURCE["REF_FASTA"],
                        viros_file,
                    ),
                    os.path.join(outdir, "ignore_accids.txt"),
                    "GENE",
                )

            if options["curate"] is False:
                entrez_descriptions = []
                for file_source, file_df in accid_file_df.groupby("file"):

                    file_descriptor = entrez_connection.run_entrez_query(
                        query_list=file_df.acc.unique().tolist(),
                    )
                    file_descriptor["file"] = file_source
                    entrez_descriptions.append(file_descriptor)

                entrez_descriptions = pd.concat(entrez_descriptions, ignore_index=True)

            else:
                entrez_descriptions = ReferenceSource.objects.all()
                entrez_descriptions = pd.DataFrame(
                    [
                        {
                            "accession": source.accid,
                            "description": source.description,
                            "taxid": source.taxid.taxid,
                        }
                        for source in entrez_descriptions
                    ]
                )

            print("Retrieved entrez descriptions")
            print(f"Number of entrez descriptions: {len(entrez_descriptions)}")
            print("Registering entrez descriptions")

            d = 0

            # raise Exception("Debugging point reached")

            for taxid_str, taxid_df in entrez_descriptions.groupby("taxid"):
                if pd.isna(taxid_str):
                    print("Skipping NaN taxid")
                    continue
                taxid_str = str(int(taxid_str))

                try:
                    ref_taxid = ReferenceTaxid.objects.get(taxid=taxid_str)
                except ReferenceTaxid.DoesNotExist:
                    ref_taxid = ReferenceTaxid.objects.create(taxid=taxid_str)

                for _, row in taxid_df.iterrows():
                    if pd.isna(row.accession):
                        print("Skipping NaN accession")
                        continue

                    if pd.isna(row.description):
                        print("Skipping NaN description")
                        continue

                    if pd.isna(row.file):
                        print("Skipping NaN file")
                        continue

                    ### register a log every 1000 accids
                    d += 1

                    if d % 1000 == 0:
                        print(f"Taxid: {taxid_str}")
                        print(f"Number of taxids processed: {d}")

                    accid_str = row.accession
                    simple_accid = accid_str.split(".")[0]

                    description = row.description

                    if len(description) > 300:
                        description = description[:300]

                    file_str = row.file

                    if (
                        ignore_dict.get(simple_accid, None) is not None
                        and options["curate"] is False
                        and "viros" in file_str
                    ):

                        ref_source = ReferenceSource.objects.filter(accid=accid_str)

                        viro_maps = ReferenceSourceFileMap.objects.filter(
                            reference_source__accid=accid_str,
                            reference_source_file__file=file_str,
                        )
                        viro_maps.delete()

                        any_left = ReferenceSourceFileMap.objects.filter(
                            reference_source__accid=accid_str,
                        )

                        if ref_source and not any_left.exists():
                            ref_source.delete()

                        continue

                    ref_source = ReferenceSource.objects.filter(accid=accid_str)

                    if ref_source.exists() is False:
                        ref_source = ReferenceSource.objects.create(
                            accid=accid_str, description=description, taxid=ref_taxid
                        )

                    elif ref_source.count() > 1:

                        ref_source.delete()
                        ref_source = ReferenceSource.objects.create(
                            accid=accid_str, description=description, taxid=ref_taxid
                        )

                    else:
                        ref_source = ref_source.first()

                    # get reference source file
                    try:
                        ref_source_file = ReferenceSourceFile.objects.get(file=file_str)
                    except ReferenceSourceFile.DoesNotExist:
                        ref_source_file = ReferenceSourceFile.objects.create(
                            file=file_str
                        )

                    # description = entrez_connection

                    try:
                        _ = ReferenceSourceFileMap.objects.get(
                            reference_source=ref_source,
                            reference_source_file=ref_source_file,
                        )

                    except ReferenceSourceFileMap.DoesNotExist:
                        _ = ReferenceSourceFileMap.objects.create(
                            reference_source=ref_source,
                            reference_source_file=ref_source_file,
                        )
        except Exception as e:
            print(e)
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_reference_update(
                    user_pk=user.pk,
                ),
                ProcessControler.FLAG_ERROR,
            )
            return

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_reference_update(
                user_pk=user.pk,
            ),
            ProcessControler.FLAG_FINISHED,
        )
