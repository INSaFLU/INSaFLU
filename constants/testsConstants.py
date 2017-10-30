'''
Created on Oct 28, 2017

@author: mmp
'''
import unittest
from constants.ConstantsTestsCase import ConstantsTestsCase
from constants.Constants import Constants
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


	def testgetPathToReferenceFile(self):
		constants = Constants()
		self.assertTrue(getattr(settings, "MEDIA_ROOT", None) + "/" + constants.DIR_PROCESSED_FILES_REFERENCE + "/userID_10/refID_20", 
			constants.get_path_to_reference_file(10, 20))

