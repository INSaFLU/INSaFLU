'''
Created on Oct 28, 2017

@author: mmp
'''
import unittest
from constants.constantsTestsCase import ConstantsTestsCase
from utils.utils import Utils
from constants.constants import Constants, FileType
from django.conf import settings 
import os

class Test(unittest.TestCase):


	### static
	constantsTestsCase = ConstantsTestsCase()
	constants = Constants()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), self.constantsTestsCase.MANAGING_TESTS)
		pass

	def tearDown(self):
		pass


	def testgetPathToReferenceFile(self):
		utils = Utils()
		self.assertTrue(getattr(settings, "MEDIA_ROOT", None) + "/" + Constants.DIR_PROCESSED_FILES_REFERENCE + "/userID_10/refID_20", 
			utils.get_path_to_reference_file(10, 20))


	def testFileNames(self):
		self.assertEqual("file_name.bam", self.constants.get_extensions_by_file_type("file_name", FileType.FILE_BAM))
		self.assertEqual("file_name.bam.bai", self.constants.get_extensions_by_file_type("file_name", FileType.FILE_BAM_BAI))
		self.assertEqual("file_name.consensus.fa", self.constants.get_extensions_by_file_type("file_name", FileType.FILE_CONSENSUS_FA))
		self.assertEqual("file_name.consensus.fasta", self.constants.get_extensions_by_file_type("file_name", FileType.FILE_CONSENSUS_FASTA))
		self.assertEqual("file_name.csv", self.constants.get_extensions_by_file_type("file_name", FileType.FILE_CSV))
		self.assertEqual("file_name.depth", self.constants.get_extensions_by_file_type("file_name", FileType.FILE_DEPTH))
		self.assertEqual("file_name.depth.gz", self.constants.get_extensions_by_file_type("file_name", FileType.FILE_DEPTH_GZ))
		self.assertEqual("file_name.depth.gz.tbi", self.constants.get_extensions_by_file_type("file_name", FileType.FILE_DEPTH_GZ_TBI))
		self.assertEqual("file_name.tab", self.constants.get_extensions_by_file_type("file_name", FileType.FILE_TAB))
		self.assertEqual("file_name.vcf", self.constants.get_extensions_by_file_type("file_name", FileType.FILE_VCF))
		self.assertEqual("file_name.vcf.gz", self.constants.get_extensions_by_file_type("file_name", FileType.FILE_VCF_GZ))
		self.assertEqual("file_name.vcf.gz.tbi", self.constants.get_extensions_by_file_type("file_name", FileType.FILE_VCF_GZ_TBI))