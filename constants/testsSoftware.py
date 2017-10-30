'''
Created on Oct 28, 2017

@author: mmp
'''
from django.test import TestCase
from constants.ConstantsTestsCase import ConstantsTestsCase
from constants.Software import Software
from django.conf import settings 
import os

class Test(TestCase):

	### static
	constantsTestsCase = ConstantsTestsCase()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), self.constantsTestsCase.MANAGING_TESTS)
		pass


	def tearDown(self):
		pass


	def testCreateFaiToFastaFile(self):
		"""
		Test samtools fai index
		"""
		## create an index file from 
		constantsTestsCase = ConstantsTestsCase()
		software = Software()
		fasta_file = os.path.join(self.baseDirectory, constantsTestsCase.MANAGING_DIR, constantsTestsCase.MANAGING_FILES_FASTA)
		if os.path.isfile(fasta_file + ".fai"): os.unlink(fasta_file + ".fai")
		try:
			software.createFaiToFastaFile(fasta_file)
			self.assertTrue(os.path.isfile(fasta_file + ".fai"))
			self.assertEquals(179, os.path.getsize(fasta_file + ".fai"))
			os.unlink(fasta_file + ".fai")
		except Exception as e:
			self.fail(e.args[0])
		
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FAIL_FASTA)
		if os.path.isfile(fasta_file + ".fai"): os.unlink(fasta_file + ".fai")
		try:
			software.createFaiToFastaFile(fasta_file)
			self.fail("Must throw exception")
		except Exception as e:
			self.assertEquals("Fail to run samtools", e.args[0])
			if os.path.isfile(fasta_file + ".fai"): os.unlink(fasta_file + ".fai")

