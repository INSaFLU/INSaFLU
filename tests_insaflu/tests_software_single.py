'''
Created on Oct 28, 2017

@author: mmp
'''

from django.test import TestCase
from django.conf import settings 
from constants.constantsTestsCase import ConstantsTestsCase
from constants.constants_mixed_infection import ConstantsMixedInfection
from constants.constants import Constants, TypePath, FileType, FileExtensions
from constants.meta_key_and_values import MetaKeyAndValue
from constants.tag_names_constants import TagNamesConstants
from utils.software import Software, Contigs2Sequences
from constants.software_names import SoftwareNames
from utils.utils import Utils
from utils.parse_out_files import ParseOutFiles
from utils.result import DecodeObjects, Coverage, CountHits
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference, Statistics
from manage_virus.uploadFiles import UploadFiles
from manage_virus.models import UploadFile
from managing_files.manage_database import ManageDatabase
from django.test.utils import override_settings
from utils.tree import CreateTree
import os, filecmp, csv
from utils.parse_in_files import ParseInFiles
from utils.result import DecodeObjects, MixedInfectionMainVector
from managing_files.models import CountVariations, MixedInfections
from utils.mixed_infections_management import MixedInfectionsManagement
from manage_virus.constants_virus import ConstantsVirus
from manage_virus.models import Tags, SeqVirus, IdentifyVirus

class Test(TestCase):

	### static
	software = Software()
	utils = Utils()
	constants = Constants()
	constants_tests_case = ConstantsTestsCase()
	software_names = SoftwareNames()

	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass


	def tearDown(self):
		pass

	
	def test_parse_sample_files(self):
		"""
		Test input files
		"""
		## create an index file from 
		
		#MANAGING_TEMPLATE_INPUT = "template_input.tsv"
		#MANAGING_TEMPLATE_INPUT_FAIL_HEADER = "template_input_fail_header.tsv"
		#MANAGING_TEMPLATE_INPUT_DATA
	
		txt_file = os.path.join("/home/mmp/insa/import_cov/77053920_SARS_CoV_2_run_004_20200410_reads_headcrop30_.txt")
		self.assertTrue(os.path.exists(txt_file))
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.set_password(ConstantsTestsCase.TEST_USER_NAME)
			user.save()
			
		parse_in_files = ParseInFiles()
		b_test_char_encoding = False
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, ParseInFiles.STATE_READ_all)
		self.assertEquals(1, parse_in_files.get_errors().get_len_vect_results())
		self.assertEquals("'utf-8' codec can't decode byte 0xff in position 0: invalid start byte", parse_in_files.get_errors().get_error(0).message)
		self.assertFalse(parse_in_files.get_errors().get_error(0).is_success())
		self.assertEquals(0, len(parse_in_files.get_vect_samples()))



