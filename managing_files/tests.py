from django.test import TestCase
from utils.ConstantsTestsCase import ConstantsTestsCase
from django.conf import settings 
from utils.Utils import Utils
import os

class testsReferenceFiles(TestCase):
	
	def test_fasta_and_gb_file(self):
		"""
		test fasta and genbank
		"""
		utils = Utils()
		constantsTestsCase = ConstantsTestsCase()
		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, constantsTestsCase.MANAGING_DIR, constantsTestsCase.MANAGING_FILES_FASTA)
		try:
			count_seqs = utils.is_fasta(fasta_file)
			self.assertEqual(8, count_seqs)
		except IOError as e:
			self.fail(e.args)
			
		try:
			self.assertFalse(utils.is_genbank(fasta_file))
			self.fail("Must throw exception")
		except IOError as e:
			pass
		
		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FAIL_FASTA)
		try:
			utils.is_fasta(fasta_file)
			self.fail("Must throw exception")
		except IOError as e:
			pass
			
		try:
			self.assertFalse(utils.is_genbank(fasta_file))
			self.fail("Must throw exception")
		except IOError as e:
			pass
		
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		try:
			utils.is_fasta(gb_file)
			self.fail("Must throw exception")
		except IOError as e:
			pass
			
		try:
			count_seqs = utils.is_genbank(gb_file)
			self.assertEqual(8, count_seqs)
		except IOError as e:
			self.fail(e.args)
		
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK_MISS_ONE)
		try:
			utils.is_fasta(gb_file)
			self.fail("Must throw exception")
		except IOError as e:
			pass
			
		try:
			count_seqs = utils.is_genbank(gb_file)
			self.assertEqual(7, count_seqs)
		except IOError as e:
			self.fail(e.args)


	def test_compare_locus_fasta_gb(self):
		"""
		Compare if fasta file has the same name of locus
		"""
		utils = Utils()
		constantsTestsCase = ConstantsTestsCase()
		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, constantsTestsCase.MANAGING_DIR, constantsTestsCase.MANAGING_FILES_FASTA)
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		
		## compare fasta and genbank
		try:
			count_seqs = utils.compare_locus_fasta_gb(fasta_file, gb_file)
			self.assertEqual(8, count_seqs)
		except ValueError as e:
			self.fail(e.message)
			
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK_MISS_ONE)
		try:
			utils.compare_locus_fasta_gb(fasta_file, gb_file)
			self.fail("Must throw exception")
		except ValueError as e:
			self.assertEqual("Number of locus are different from fasta to genbank.", e.args[0])
			
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK_DIFF_LENGTH)
		try:
			utils.compare_locus_fasta_gb(fasta_file, gb_file)
			self.fail("Must throw exception")
		except ValueError as e:
			self.assertEqual("Different length. Fasta seq: PB2 length: 2280; Fasta seq: PB2 length: 2277.", e.args[0])



