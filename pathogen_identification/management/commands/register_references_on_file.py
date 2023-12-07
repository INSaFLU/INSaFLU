import os
import time
from abc import ABC, abstractmethod
from datetime import date

import pandas as pd
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from django.db.models import Q

from constants.constants import Televir_Metadata_Constants
from pathogen_identification.models import (
    ReferenceSource,
    ReferenceSourceFile,
    ReferenceSourceFileMap,
    ReferenceTaxid,
)
from pathogen_identification.utilities.entrez_wrapper import EntrezWrapper


class Command(BaseCommand):
    help = "deploy run"

    def add_arguments(self, parser):
        parser.add_argument(
            "--user_id",
            type=int,
            help="user id",
            default=None,
        )

        parser.add_argument(
            "-o",
            "--outdir",
            type=str,
            help="output directory",
        )

    def handle(self, *args, **options):
        ###
        # get user
        if not options["user_id"]:
            username = "admin"

        user = User.objects.get(pk=options["user_id"])
        outdir = options["outdir"]
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

        entrez_descriptions = entrez_connection.run_entrez_query(
            query_list=accid_file_df.acc.unique()[:1000].tolist(),
        )

        for taxid_str, taxid_df in entrez_descriptions.groupby("taxid"):
            try:
                ref_taxid = ReferenceTaxid.objects.get(taxid=taxid_str)
            except ReferenceTaxid.DoesNotExist:
                ref_taxid = ReferenceTaxid.objects.create(taxid=taxid_str)

            for ix, row in taxid_df.iterrows():
                accid_str = row.accession
                # taxid_str = row.taxid
                description = row.description

                if len(description) > 300:
                    description = description[:300]

                files = accid_file_df[accid_file_df.acc == accid_str].file

                try:
                    ref_source = ReferenceSource.objects.get(accid=accid_str)
                except ReferenceSource.DoesNotExist:
                    ref_source = ReferenceSource.objects.create(
                        accid=accid_str, description=description, taxid=ref_taxid
                    )

                # get reference source file

                for file_str in files:
                    try:
                        ref_source_file = ReferenceSourceFile.objects.get(file=file_str)
                    except ReferenceSourceFile.DoesNotExist:
                        ref_source_file = ReferenceSourceFile.objects.create(
                            file=file_str
                        )

                    description = entrez_connection

                    try:
                        ref_source_file_map = ReferenceSourceFileMap.objects.get(
                            reference_source=ref_source,
                            reference_source_file=ref_source_file,
                        )

                    except ReferenceSourceFileMap.DoesNotExist:
                        ref_source_file_map = ReferenceSourceFileMap.objects.create(
                            reference_source=ref_source,
                            reference_source_file=ref_source_file,
                        )
