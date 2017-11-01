'''
Created on Oct 28, 2017

@author: mmp
'''
from django.test import TestCase
from utils.ConstantsTestsCase import ConstantsTestsCase
from utils.Constants import Constants
from utils.Software import Software
from utils.Utils import Utils
from utils.ParseOutFiles import ParseOutFiles
from django.conf import settings 
import os, sys

class Test(TestCase):

	### static
	software = Software()
	utils = Utils()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass


	def tearDown(self):
		pass


	def testCreateFaiToFastaFile(self):
		"""
		Test samtools fai index
		"""
		## create an index file from 
		
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		if os.path.isfile(fasta_file + ".fai"): os.unlink(fasta_file + ".fai")
		try:
			self.software.createFaiToFastaFile(fasta_file)
			self.assertTrue(os.path.isfile(fasta_file + ".fai"))
			self.assertEquals(179, os.path.getsize(fasta_file + ".fai"))
			os.unlink(fasta_file + ".fai")
		except Exception as e:
			self.fail(e.args[0])
		
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FAIL_FASTA)
		if os.path.isfile(fasta_file + ".fai"): os.unlink(fasta_file + ".fai")
		try:
			self.software.createFaiToFastaFile(fasta_file)
			self.fail("Must throw exception")
		except Exception as e:
			self.assertEquals("Fail to run samtools", e.args[0])
			if os.path.isfile(fasta_file + ".fai"): os.unlink(fasta_file + ".fai")

	def test_is_exist_database_abricate(self):
		self.assertTrue(self.software.is_exist_database_abricate("ncbi"))
		self.assertFalse(self.software.is_exist_database_abricate("xpot"))
	
	def test_create_database_abricate(self):
		database_name = "xpto"
		if (self.software.is_exist_database_abricate(database_name)):
			cmd = "rm -r %s/%s*" % (Software.SOFTWARE_ABRICATE_DB, database_name)
			exist_status = os.system(cmd)
			self.assertTrue(exist_status == 0)
		self.assertFalse(self.software.is_exist_database_abricate(database_name))
		temp_file_file = os.path.join(self.baseDirectory, Constants.DIR_TYPE_IDENTIFICATION, ConstantsTestsCase.MANAGING_TEST_INFLUENZA_FILE)
		self.assertTrue(os.path.exists(temp_file_file))
		
		try:
			self.software.create_database_abricate(database_name, "fdssfd")
			self.fail("must throw error")
		except IOError:
			pass
		
		self.software.create_database_abricate(database_name, temp_file_file)
		self.assertTrue(self.software.is_exist_database_abricate(database_name))
		cmd = "rm -r %s/%s*" % (Software.SOFTWARE_ABRICATE_DB, database_name)
		exist_status = os.system(cmd)
		self.assertTrue(exist_status == 0)
		self.assertFalse(self.software.is_exist_database_abricate(database_name))
		
	def testRunSpadesAndAbricate(self):
		"""
		Test run spades and abricator
		"""
		b_run_this = getattr(settings, "RUN_TEST_IN_COMMAND_LINE", None)
		
		if (not b_run_this): return	## it only work from command line because of PYTHONPATH defined on the eclispe IDE 
		fastq1_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		fastq1_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)
		self.assertTrue(os.path.isfile(fastq1_1))
		self.assertTrue(os.path.isfile(fastq1_2))
		
		out_dir = self.utils.get_temp_dir()
		self.assertTrue(os.path.isdir(out_dir))
		cmd = self.software.run_spades(fastq1_1, fastq1_2, out_dir)
		file_out = os.path.join(out_dir, "contigs.fasta")

#		file_out = os.path.join("/tmp/insaFlu/insa_flu_85008740", "contigs.fasta")
		self.assertTrue(os.path.exists(file_out))
		self.assertTrue(os.path.getsize(file_out) > 100)
		
		database_name = "xpto"
		if (not self.software.is_exist_database_abricate(database_name)):
			temp_file_file = os.path.join(self.baseDirectory, Constants.DIR_TYPE_IDENTIFICATION, ConstantsTestsCase.MANAGING_TEST_INFLUENZA_FILE)
			self.assertTrue(os.path.exists(temp_file_file))
			self.software.create_database_abricate(database_name, temp_file_file)
		
		out_file = self.utils.get_temp_file("temp_abricate", ".txt")
		cmd = self.software.run_abricate(database_name, file_out, out_file)
		self.assertTrue(os.path.exists(out_file))

		parseOutFiles = ParseOutFiles()
		vect_data = parseOutFiles.parse_abricate_file(out_file)
		self.assertEqual('A', vect_data[0][ParseOutFiles.GENE])
		self.assertEqual(100.00, vect_data[0][ParseOutFiles.COVERAGE])
		self.assertEqual(99.69, vect_data[0][ParseOutFiles.IDENTITY])
		self.assertEqual('influenza_type', vect_data[0][ParseOutFiles.TYPE])
		self.assertEqual('XXXX', vect_data[0][ParseOutFiles.ACCESSION])
		
		self.assertEqual('N2', vect_data[3][ParseOutFiles.GENE])
		self.assertEqual(16.03, vect_data[3][ParseOutFiles.COVERAGE])
		self.assertEqual(96.46, vect_data[3][ParseOutFiles.IDENTITY])
		self.assertEqual('influenza_subtype', vect_data[3][ParseOutFiles.TYPE])
		self.assertEqual('XXXXXX', vect_data[3][ParseOutFiles.ACCESSION])

		## remove file
		os.unlink(out_file)
		
		## remove abricate db
		cmd = "rm -r %s/%s*" % (Software.SOFTWARE_ABRICATE_DB, database_name)
		exist_status = os.system(cmd)
		self.assertTrue(exist_status == 0)
		