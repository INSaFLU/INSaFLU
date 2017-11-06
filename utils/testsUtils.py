'''
Created on Oct 28, 2017

@author: mmp
'''
import unittest
from .constantsTestsCase import ConstantsTestsCase
from .utils import Utils
from django.conf import settings 
import os

class Test(unittest.TestCase):


	### static
	constantsTestsCase = ConstantsTestsCase()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), self.constantsTestsCase.MANAGING_TESTS)
		pass

	def tearDown(self):
		pass


	def testgetFileNameWithoutExtension(self):
		utils = Utils()
		self.assertEqual("fiels", utils.get_file_name_without_extension("/root/fiels.fasta"))
		self.assertEqual("fiels", utils.get_file_name_without_extension("/root/fiels"))
		self.assertEqual("fiels.fasta", utils.get_file_name_without_extension("/root/fiels.fasta.txt"))

	def test_is_fastq_gz(self):
		utils = Utils()
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		self.assertTrue(utils.is_fastq_gz(path_file))
		
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_ABRICATE, ConstantsTestsCase.MANAGING_TEST_ABRICATE)
		try:
			self.assertTrue(utils.is_fastq_gz(path_file))
		except Exception as e:
			self.assertEqual("File is not compressed in gzip format", e.args[0])
			
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		try:
			self.assertTrue(utils.is_fastq_gz(path_file))
		except Exception as e:
			self.assertEqual("File is not compressed in gzip format", e.args[0])
			
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA_FAKE_GZ)
		try:
			self.assertTrue(utils.is_fastq_gz(path_file))
		except Exception as e:
			self.assertEqual("File is not in fastq.gz format", e.args[0])
			