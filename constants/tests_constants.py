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
		
	def testConstants(self):
		self.assertEqual(self.constants.complement("AAATTTCCC"), "TTTAAAGGG", "must be equal")
		self.assertEqual(self.constants.complement("AAATTGGGTCCC"), "TTTAACCCAGGG", "must be equal")
		self.assertEqual(self.constants.reverse_complement("AAGGGATTTCCC"), "GGGAAATCCCTT", "must be equal")
		self.assertEqual(self.constants.reverse_complement("AAATTAGTCCC"), "GGGACTAATTT", "must be equal")
		self.assertEqual(self.constants.ambiguos_to_unambiguous("YARWATTTCCC"), "[TC]A[AG][AT]ATTTCCC", "must be equal")
		self.assertEqual(self.constants.ambiguos_to_unambiguous("RYKMSWBDHVN"), "[AG][TC][GT][AC][GC][AT][CGT][AGT][ACT][ACG][ACGT]", "must be equal")
		self.assertEqual(self.constants.complement("RYKMSWBDHVN"), "YRMKSWVHDBN", "must be equal")

	def testGet_diff_between_two_seq(self):
		self.assertEqual(self.constants.get_diff_between_two_seq("AAATTTCCC", "TTTAAAGGG"), 9)
		self.assertEqual(self.constants.get_diff_between_two_seq("AAATTTCCC", "TTTAAAGGGS"), 0)
		self.assertEqual(self.constants.get_diff_between_two_seq("AAATTTCCC", "AAATTTCCC"), 0)
		self.assertEqual(self.constants.get_diff_between_two_seq("GAATTTCCC", "AAATTTCCC"), 1)
		self.assertEqual(self.constants.get_diff_between_two_seq("AAATTTCCC", "AAATTTCCG"), 1)

	def test_is_poly_n(self):
		self.assertTrue(self.constants.is_poly_n("AAAAAAA"))
		self.assertFalse(self.constants.is_poly_n("AAAAAACA"))
		self.assertTrue(self.constants.is_poly_n("CCCCCCCCCCCCC"))
		self.assertTrue(self.constants.is_poly_n("TTTTTTTT"))
		self.assertFalse(self.constants.is_poly_n("TTCTTTTTT"))
		
	def test_short_name(self):
		"""
		test short names
		"""
		self.assertEquals('short_name_short_name_short_name', self.constants.short_name('short_name_short_name_short_name', 100))
		self.assertEquals('short..._name', self.constants.short_name('short_name_short_name_short_name', 10))
