import os
from datetime import date
from typing import List

import pandas as pd
from Bio import SeqIO
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from managing_files.models import ProcessControler
from pathogen_identification.models import (
    ReferenceSource,
    ReferenceSourceFile,
    ReferenceSourceFileMap,
    ReferenceTaxid,
)
from pathogen_identification.utilities.reference_utils import raw_reference_to_insaflu
from pathogen_identification.utilities.televir_bioinf import TelevirBioinf
from utils.process_SGE import ProcessSGE


class Command(BaseCommand):
    help = "deploy televir sample reports sort"

    def add_arguments(self, parser):
        parser.add_argument(
            "--user_id",
            type=int,
            help="user deploying the run (pk)",
        )

        parser.add_argument(
            "--file_id",
            type=int,
            help="file_id",
        )

        parser.add_argument(
            "--fasta",
            type=str,
            help="fasta",
        )

        parser.add_argument(
            "--metadata",
            type=str,
            help="metadata",
        )

    def handle(self, *args, **options):
        ###
        # SETUP
        user_id = int(options["user_id"])
        user = User.objects.get(pk=user_id)

        file_id = int(options["file_id"])
        file = ReferenceSourceFile.objects.get(pk=file_id)
        filepath = file.filepath
        filepath_dir = os.path.dirname(filepath)
        os.makedirs(filepath_dir, exist_ok=True)

        bioinf_utils = TelevirBioinf()

        metadata_path = options["metadata"]
        fasta_path = options["fasta"]

        # PROCESS CONTROLER
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_file_upload(
                file_id=file_id,
            ),
            ProcessControler.FLAG_RUNNING,
        )

        process = ProcessControler.objects.filter(
            owner__id=user.pk,
            name=process_controler.get_name_televir_file_upload(
                file_id=file_id,
            ),
            is_finished=False,
        )

        if process.exists():
            process = process.first()

            print(f"Process {process.pk} has been submitted")

        else:
            print("Process does not exist")

        description = file.description

        # UTILITIES

        try:
            metadata_table = pd.read_csv(metadata_path, sep="\t")

            if "Description" not in metadata_table.columns:
                description_ref = f"Panel reference"
                if description != "":
                    description_ref = description_ref + f" - {description}"
                metadata_table["Description"] = description_ref

            with open(fasta_path) as handle_in:
                for record in SeqIO.parse(handle_in, "fasta"):
                    if (
                        len(metadata_table) > 0
                        and not record.id in metadata_table["Accession ID"].values
                    ):

                        continue

                    with open(filepath, "a") as handle_out:
                        SeqIO.write(record, handle_out, "fasta")

                    taxid = metadata_table.loc[
                        metadata_table["Accession ID"] == record.id
                    ]["TaxID"].values[0]

                    description = metadata_table.loc[
                        metadata_table["Accession ID"] == record.id
                    ]["Description"].values[0]

                    taxid = ReferenceTaxid.objects.get_or_create(
                        taxid=taxid,
                    )[0]

                    reference_source = ReferenceSource.objects.create(
                        accid=record.id,
                        taxid=taxid,
                        description=description,
                    )

                    ReferenceSourceFileMap.objects.create(
                        reference_source_file=file,
                        reference_source=reference_source,
                    )

                ## bgzip compress filepath and create index
                file_path_gz = bioinf_utils.bgzip(filepath)
                bioinf_utils.index_fasta(file_path_gz)

                ## update file path
                file.file = os.path.basename(file_path_gz)
                file.save()

                # delete files
                os.remove(metadata_path)
                os.remove(fasta_path)

        except Exception as e:
            print(e)
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_file_upload(
                    file_id=file_id,
                ),
                ProcessControler.FLAG_ERROR,
            )
            raise e

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_file_upload(
                file_id=file_id,
            ),
            ProcessControler.FLAG_FINISHED,
        )
