'''
Created on Nov 16, 2017

@author: mmp
'''
import unittest

from utils.constants import Constants
from utils.constants import FileType

class Test(unittest.TestCase):


	def testFileNames(self):
		constants = Constants()

		self.assertEqual("file_name.bam", constants.get_extensions_by_file_type("file_name", FileType.FILE_BAM))
		self.assertEqual("file_name.bam.bai", constants.get_extensions_by_file_type("file_name", FileType.FILE_BAM_BAI))
		self.assertEqual("file_name.consensus.fa", constants.get_extensions_by_file_type("file_name", FileType.FILE_CONSENSUS_FA))
		self.assertEqual("file_name.consensus.fasta", constants.get_extensions_by_file_type("file_name", FileType.FILE_CONSENSUS_FASTA))
		self.assertEqual("file_name.csv", constants.get_extensions_by_file_type("file_name", FileType.FILE_CSV))
		self.assertEqual("file_name.depth", constants.get_extensions_by_file_type("file_name", FileType.FILE_DEPTH))
		self.assertEqual("file_name.depth.gz", constants.get_extensions_by_file_type("file_name", FileType.FILE_DEPTH_GZ))
		self.assertEqual("file_name.depth.gz.tbi", constants.get_extensions_by_file_type("file_name", FileType.FILE_DEPTH_GZ_TBI))
		self.assertEqual("file_name.tab", constants.get_extensions_by_file_type("file_name", FileType.FILE_TAB))
		self.assertEqual("file_name.vcf", constants.get_extensions_by_file_type("file_name", FileType.FILE_VCF))
		self.assertEqual("file_name.vcf.gz", constants.get_extensions_by_file_type("file_name", FileType.FILE_VCF_GZ))
		self.assertEqual("file_name.vcf.gz.tbi", constants.get_extensions_by_file_type("file_name", FileType.FILE_VCF_GZ_TBI))

