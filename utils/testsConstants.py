'''
Created on Oct 28, 2017

@author: mmp
'''
import unittest
from .constantsTestsCase import ConstantsTestsCase
from .utils import Utils
from .constants import Constants
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
		utils = Utils()
		self.assertTrue(getattr(settings, "MEDIA_ROOT", None) + "/" + Constants.DIR_PROCESSED_FILES_REFERENCE + "/userID_10/refID_20", 
			utils.get_path_to_reference_file(10, 20))

