import os
import subprocess
from typing import List

from constants.constants import Televir_Metadata_Constants
from pathogen_identification.models import RawReference


class TelevirBioinf:
    def __init__(self):
        self.metadata_constants = Televir_Metadata_Constants()
        self.samtools_binary = self.metadata_constants.get_software_binary("samtools")

    @staticmethod
    def virosaurus_formatting(accid):
        accid_simple = accid.split(".")[0]
        return f"{accid_simple}:{accid_simple};"

    @staticmethod
    def check_sourcefile_is_virosaurus(file_path):
        if "virosaurus" in file_path:
            return True

        return False

    def process_accid(self, accid, source_file):
        if self.check_sourcefile_is_virosaurus(source_file):
            return self.virosaurus_formatting(accid)

        return accid

    def check_file_exists_not_empty(self, file_path):
        return os.path.exists(file_path) and os.path.getsize(file_path) > 100

    def extract_reference(self, source_file, accid, output_file):
        accid_code = self.process_accid(accid, source_file)
        command = (
            f'{self.samtools_binary} faidx {source_file} "{accid_code}" > {output_file}'
        )
        print(command)
        subprocess.call(command, shell=True)

        return self.check_file_exists_not_empty(output_file)

    def merge_references(self, files: List[str], output_file):
        command = f"cat {' '.join(files)} > {output_file}"
        subprocess.call(command, shell=True)

        return self.check_file_exists_not_empty(output_file)
