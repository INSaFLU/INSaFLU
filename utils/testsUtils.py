'''
Created on Oct 28, 2017

@author: mmp
'''
import unittest
from utils.ConstantsTestsCase import ConstantsTestsCase
from utils.Utils import Utils
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
