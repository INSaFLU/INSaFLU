'''
Created on Mar 3, 2018

@author: mmp
'''
import unittest
from constants.constantsTestsCase import ConstantsTestsCase
from constants.software_names import SoftwareNames
from utils.utils import Utils
from django.conf import settings
from utils.software import Software, Contigs2Sequences
import os, filecmp


class Test(unittest.TestCase):

	utils = Utils()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass


	def test_get_most_recent_database_file(self):

		contigs_2_sequences = Contigs2Sequences(True)
		(version, file) = contigs_2_sequences.get_most_recent_database()
		self.assertEqual("22", version)
		self.assertEqual(os.path.join(self.baseDirectory, "db/contigs2sequences/sequences_v22.fasta"), file)
		
		
	def test_identify_contigs(self):
	
		software = Software()
		contigs_2_sequences = Contigs2Sequences(True)
	
		fastq1_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		fastq1_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)
		self.assertTrue(os.path.isfile(fastq1_1))
		self.assertTrue(os.path.isfile(fastq1_2))
	
		out_dir = self.utils.get_temp_dir()
		self.assertTrue(os.path.isdir(out_dir))
		cmd = software.run_spades(fastq1_1, fastq1_2, out_dir)
		file_out = os.path.join(out_dir, "contigs.fasta")
		self.assertTrue(os.path.exists(file_out))
		self.assertTrue(os.path.getsize(file_out) > 100)
	
		(out_file, clean_abricate_file) = contigs_2_sequences.identify_contigs(file_out, 'temp3.txt')
	
		expect_fasta = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_PV1_KX162693_spades_out_fasta)
		self.assertTrue(filecmp.cmp(out_file, expect_fasta))
	
		expected_abricate = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "expected_abricate_clean_2.txt")
		self.assertTrue(filecmp.cmp(expected_abricate, clean_abricate_file))
	
		## remove abricate db
		cmd = "rm -r %s/%s*" % (SoftwareNames.SOFTWARE_ABRICATE_DB, contigs_2_sequences.get_database_name())
		exist_status = os.system(cmd)
		self.assertTrue(exist_status == 0)
		self.utils.remove_dir(out_dir)
		os.unlink(out_file)
		os.unlink(clean_abricate_file)


