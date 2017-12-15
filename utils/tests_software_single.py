'''
Created on Oct 28, 2017

@author: mmp
'''

from django.test import TestCase
from django.conf import settings 
from constants.constantsTestsCase import ConstantsTestsCase
from constants.constants import Constants, TypePath, FileType, FileExtensions
from constants.meta_key_and_values import MetaKeyAndValue
from constants.tag_names_constants import TagNamesConstants
from utils.software import Software
from constants.software_names import SoftwareNames
from utils.utils import Utils
from utils.parse_out_files import ParseOutFiles
from utils.result import DecodeResultAverageAndNumberReads, DecodeObjects, Coverage, CountHits
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference, Statistics
from manage_virus.uploadFiles import UploadFiles
from manage_virus.models import UploadFile
from managing_files.manage_database import ManageDatabase
from django.test.utils import override_settings
from utils.tree import CreateTree
import os, filecmp


from managing_files.models import CountVariations

class Test(TestCase):

	### static
	software = Software()
	utils = Utils()
	constants = Constants()
	software_names = SoftwareNames()
	constants_tests_case = ConstantsTestsCase()

	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass


	def tearDown(self):
		pass

	
