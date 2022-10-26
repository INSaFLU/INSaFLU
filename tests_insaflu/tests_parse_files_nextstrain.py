'''
Created on Oct 28, 2017

@author: mmp
'''
from django.test import TestCase
from constants.constantsTestsCase import ConstantsTestsCase
from utils.parse_in_files_nextstrain import ParseNextStrainFiles
from django.contrib.auth.models import User
from django.conf import settings
from utils.utils import Utils
import os

class Test(TestCase):

	### static
	utils = Utils()
	
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
	
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_DATASET_FILES, "Nextstrain_metadata.tsv")
		self.assertTrue(os.path.exists(txt_file))
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.set_password(ConstantsTestsCase.TEST_USER_NAME)
			user.save()
			
		parse_nextstrin_files = ParseNextStrainFiles()
		b_test_char_encoding = True
		parse_nextstrin_files.parse_nextstrain_files(txt_file, user, b_test_char_encoding,
				ParseNextStrainFiles.STATE_READ_metadata_dont_detect_errors_and_chech_nexttrain)
		self.assertEquals(0, parse_nextstrin_files.get_errors().get_len_vect_results())
		self.assertEquals(9, len(parse_nextstrin_files.get_vect_samples()))
		self.assertEquals(21, len(parse_nextstrin_files.vect_header))
		self.assertEquals('strain', parse_nextstrin_files.vect_header[0])
		self.assertEquals('Scorpio Pangolin', parse_nextstrin_files.vect_header[-1])
		self.assertEquals('2022-06-20', parse_nextstrin_files.get_value('ERR4082026_covid_minion', 'date'))
		self.assertEquals('Europe', parse_nextstrin_files.get_value('ERR4082026_covid_minion', 'region'))
		self.assertEquals('B.1.93', parse_nextstrin_files.get_value('ERR4082025_covid_minion_1', 'Lineage Pangolin'))
		
	def test_parse_sample_files_fail(self):
		"""
		Test input files
		"""
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_DATASET_FILES, "Nextstrain_metadata_fail_date_and_age.tsv")
		self.assertTrue(os.path.exists(txt_file))
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.set_password(ConstantsTestsCase.TEST_USER_NAME)
			user.save()
			
		parse_nextstrin_files = ParseNextStrainFiles()
		b_test_char_encoding = True
		parse_nextstrin_files.parse_nextstrain_files(txt_file, user, b_test_char_encoding,
				ParseNextStrainFiles.STATE_READ_metadata_only_detect_errors_and_chech_nexttrain)
		self.assertEquals(12, parse_nextstrin_files.get_errors().get_len_vect_results())
		self.assertTrue(parse_nextstrin_files.get_errors().has_errors())
		self.assertEquals("The 'date' must have this format YYYY-MM-DD. Line: 7 Column: 2", str(parse_nextstrin_files.errors.vect_results[6].message))
		self.assertEquals(9, len(parse_nextstrin_files.get_vect_samples()))
		self.assertEquals(21, len(parse_nextstrin_files.vect_header))
		self.assertEquals('strain', parse_nextstrin_files.vect_header[0])
		self.assertEquals('Scorpio Pangolin', parse_nextstrin_files.vect_header[-1])
		self.assertEquals('2022-06-20', parse_nextstrin_files.get_value('ERR4082026_covid_minion', 'date'))
		self.assertEquals('Europe', parse_nextstrin_files.get_value('ERR4082026_covid_minion', 'region'))
		self.assertEquals('B.1.93', parse_nextstrin_files.get_value('ERR4082025_covid_minion_1', 'Lineage Pangolin'))


