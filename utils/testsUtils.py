'''
Created on Oct 28, 2017

@author: mmp
'''
import unittest
from utils.constantsTestsCase import ConstantsTestsCase
from django.test.utils import override_settings
from utils.utils import Utils
from utils.software import Software
from django.conf import settings 
import os

class Test(unittest.TestCase):


	### static
	constantsTestsCase = ConstantsTestsCase()

	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), self.constantsTestsCase.MANAGING_TESTS)
		
	def tearDown(self):
		pass

	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_media_root(self):
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		
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
			self.assertEqual("File need to have suffix '.fastq.gz'", e.args[0])
			
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		try:
			self.assertTrue(utils.is_fastq_gz(path_file))
		except Exception as e:
			self.assertEqual("File need to have suffix '.fastq.gz'", e.args[0])
			
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA_FAKE_GZ)
		try:
			self.assertTrue(utils.is_fastq_gz(path_file))
		except Exception as e:
			self.assertEqual("File is not in fastq.gz format", e.args[0])
	
	def test_read_text_file(self):
		"""
		"""
		utils = Utils()
		txt_file = os.path.join(self.baseDirectory, 'abricate', ConstantsTestsCase.MANAGING_TEST_ABRICATE + ".xpto")
		self.assertFalse(os.path.exists(txt_file))
		try:
			utils.read_text_file(txt_file)
			self.fail("must raise error")
		except IOError as e:
			self.assertEqual("Error: file '/home/mmp/eclipse_oxygen/fluwebvirus/static_all/tests/abricate/abricate_out.txt.xpto' doens't exist.", e.args[0])
		
		txt_file = os.path.join(self.baseDirectory, 'abricate', ConstantsTestsCase.MANAGING_TEST_ABRICATE)
		self.assertTrue(os.path.exists(txt_file))
		
		vect_out = utils.read_text_file(txt_file)
		self.assertEquals(6, len(vect_out))
		self.assertEquals(0, vect_out[0].find("#FILE"))
		self.assertEquals(-1, vect_out[0].find("Victoria"))
		self.assertEquals(87, vect_out[5].find("Victoria"))
		
	def test_gzip_and_tabix(self):
		utils = Utils()
		vcf_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "temp.vcf")
		self.assertTrue(os.path.exists(vcf_file))
		
		temp_path = utils.get_temp_dir()
		temp_file = os.path.join(temp_path, os.path.basename(utils.get_temp_file("temp", ".vcf")))
		utils.copy_file(vcf_file, temp_file)
		self.assertTrue(os.path.exists(vcf_file))
		utils.compress_files(Software.SOFTWARE_BGZIP, temp_file)
		utils.compress_files(Software.SOFTWARE_BGZIP, temp_file)
		self.assertTrue(os.path.exists(temp_file + ".gz"))
		self.assertTrue(os.path.getsize(temp_file + ".gz") > 1000)
		utils.create_index_files(Software.SOFTWARE_TABIX, temp_file + ".gz")
		self.assertTrue(os.path.exists(temp_file + ".gz.tbi"))
		utils.remove_temp_file(temp_file + ".gz.tbi")
		try:
			utils.create_index_files(Software.SOFTWARE_TABIX, temp_file)
			self.fail("must raise exception")
		except Exception as e:
			self.assertEquals("Fail to create index", e.args[0])
		utils.remove_dir(temp_path)
		
		