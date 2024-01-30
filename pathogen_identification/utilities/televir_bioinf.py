

from constants.constants import Televir_Metadata_Constants
import os
import subprocess

class TelevirBioinf:

    def __init__(self):
        self.metadata_constants = Televir_Metadata_Constants()
        self.samtools_binary = self.metadata_constants.get_software_binary("samtools")
    
    def check_file_exists_not_empty(self, file_path):
        return os.path.exists(file_path) and os.path.getsize(file_path) > 100
    
    
    def extract_reference(self, source_file, accid, output_file):
        command = f"{self.samtools_binary} faidx {source_file} {accid} > {output_file}"
        print(command)
        subprocess.call(command, shell=True)

        return self.check_file_exists_not_empty(output_file)