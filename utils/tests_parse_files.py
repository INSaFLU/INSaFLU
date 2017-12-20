'''
Created on Oct 28, 2017

@author: mmp
'''
from django.test import TestCase
from constants.constantsTestsCase import ConstantsTestsCase
from constants.constants import TypeFile
from utils.parse_in_files import ParseInFiles
from django.contrib.auth.models import User
from managing_files.models import UploadFiles, MetaKey
from django.conf import settings
import os

class Test(TestCase):

	### static
	
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass


	def tearDown(self):
		pass


	def test_parse_abricate_file(self):
		"""
		Test input files
		"""
		## create an index file from 
		
		#MANAGING_TEMPLATE_INPUT = "template_input.tsv"
		#MANAGING_TEMPLATE_INPUT_FAIL_HEADER = "template_input_fail_header.tsv"
		#MANAGING_TEMPLATE_INPUT_DATA
	
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_INPUT_FILES, ConstantsTestsCase.MANAGING_TEMPLATE_INPUT)
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
		only_detect_errors = False
		b_test_char_encoding = False
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, only_detect_errors)
		self.assertEquals(1, parse_in_files.get_errors().get_len_vect_results())
		self.assertEquals("'utf-8' codec can't decode byte 0xff in position 0: invalid start byte", parse_in_files.get_errors().get_error(0).message)
		self.assertFalse(parse_in_files.get_errors().get_error(0).is_success())
		self.assertEquals(0, len(parse_in_files.get_vect_samples()))

	def test_parse_abricate_file_1(self):
		"""
		Test input files
		"""
		## create an index file from 
		
		#MANAGING_TEMPLATE_INPUT = "template_input.tsv"
		#MANAGING_TEMPLATE_INPUT_FAIL_HEADER = "template_input_fail_header.tsv"
		#MANAGING_TEMPLATE_INPUT_DATA
	
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_INPUT_FILES, ConstantsTestsCase.MANAGING_TEMPLATE_INPUT)
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
		only_detect_errors = False
		b_test_char_encoding = True
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, only_detect_errors)
		self.assertEquals(1, parse_in_files.get_errors().get_len_vect_results())
		self.assertEquals("There's no samples to process.", parse_in_files.get_errors().get_error(0).message)
		self.assertFalse(parse_in_files.get_errors().get_error(0).is_success())
		self.assertEquals(0, len(parse_in_files.get_vect_samples()))
		
		
	def test_parse_abricate_file_2(self):
		"""
		Test input files
		"""
		## create an index file from 
		#MANAGING_TEMPLATE_INPUT = "template_input.tsv"
		#MANAGING_TEMPLATE_INPUT_FAIL_HEADER = "template_input_fail_header.tsv"
		#MANAGING_TEMPLATE_INPUT_DATA
	
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_INPUT_FILES, ConstantsTestsCase.MANAGING_TEMPLATE_INPUT_FAIL_HEADER)
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
		only_detect_errors = False
		b_test_char_encoding = True
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, only_detect_errors)
		self.assertEquals(1, parse_in_files.get_errors().get_len_vect_results())
		self.assertEquals("Header not found in the file. Please, check the names in the header, must be equal and have the same order of the template file.", parse_in_files.get_errors().get_error(0).message)
		self.assertFalse(parse_in_files.get_errors().get_error(0).is_success())
		self.assertEquals(0, len(parse_in_files.get_vect_samples()))


	def test_parse_abricate_file_3(self):
		"""
		Test input files
		"""
		## create an index file from 
		#MANAGING_TEMPLATE_INPUT = "template_input.tsv"
		#MANAGING_TEMPLATE_INPUT_FAIL_HEADER = "template_input_fail_header.tsv"
		#MANAGING_TEMPLATE_INPUT_DATA
	
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_INPUT_FILES, ConstantsTestsCase.MANAGING_TEMPLATE_INPUT_DATA_TSV)
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
		only_detect_errors = False
		b_test_char_encoding = True
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, only_detect_errors)
		self.assertEquals(0, parse_in_files.get_errors().get_len_vect_results())
		self.assertEquals(2, len(parse_in_files.get_vect_samples()))
		self.assertEquals('xpto1', parse_in_files.get_vect_samples()[0][0].name)
		self.assertEquals('asds.fastq.gz', parse_in_files.get_vect_samples()[0][0].file_name_1)
		self.assertEquals('asdws.fastq.gz', parse_in_files.get_vect_samples()[0][0].file_name_2)
		self.assertEquals('Generic', parse_in_files.get_vect_samples()[0][0].data_set.name)
		self.assertEquals(12, parse_in_files.get_vect_samples()[0][0].week)
		self.assertEquals('SRID=4326;POINT (12.2 14.3)', str(parse_in_files.get_vect_samples()[0][0].geo_local))
		self.assertEquals('vaccine status', parse_in_files.get_vect_samples()[0][0].vaccine_status.name)
		self.assertEquals('2017-12-13 00:00:00', str(parse_in_files.get_vect_samples()[0][0].date_of_collection))
		self.assertEquals('2017-12-12 00:00:00', str(parse_in_files.get_vect_samples()[0][0].date_of_onset))
		self.assertEquals(2, len(parse_in_files.get_vect_samples()[0][1]))
		self.assertEquals('other fields from here', parse_in_files.get_vect_samples()[0][1][0].value)
		self.assertEquals('other fields from here', parse_in_files.get_vect_samples()[0][1][0].tag_name.name)
		self.assertEquals('other fields from here11', parse_in_files.get_vect_samples()[0][1][1].value)
		self.assertEquals('other fields from here1', parse_in_files.get_vect_samples()[0][1][1].tag_name.name)
		self.assertEquals('xpto2', parse_in_files.get_vect_samples()[1][0].name)
		self.assertEquals('2asds.fastq.gz', parse_in_files.get_vect_samples()[1][0].file_name_1)
		self.assertEquals('2asds3.fastq.gz', parse_in_files.get_vect_samples()[1][0].file_name_2)
		self.assertEquals(11, parse_in_files.get_vect_samples()[1][0].week)
		self.assertEquals('data set', parse_in_files.get_vect_samples()[1][0].data_set.name)
		self.assertEquals(None, parse_in_files.get_vect_samples()[1][0].vaccine_status)
		self.assertEquals('2017-12-23 00:00:00', str(parse_in_files.get_vect_samples()[1][0].date_of_collection))
		self.assertEquals('2017-12-22 00:00:00', str(parse_in_files.get_vect_samples()[1][0].date_of_onset))
		self.assertEquals(1, len(parse_in_files.get_vect_samples()[1][1]))
		self.assertEquals('other fields from here2', parse_in_files.get_vect_samples()[1][1][0].value)
		self.assertEquals('other fields from here', parse_in_files.get_vect_samples()[1][1][0].tag_name.name)
		self.assertEquals('SRID=4326;POINT (30.2 43.5)', str(parse_in_files.get_vect_samples()[1][0].geo_local))
		
	
	def test_parse_abricate_file_4(self):
		"""
		Test input files
		"""
		## create an index file from 
		#MANAGING_TEMPLATE_INPUT = "template_input.tsv"
		#MANAGING_TEMPLATE_INPUT_FAIL_HEADER = "template_input_fail_header.tsv"
		#MANAGING_TEMPLATE_INPUT_DATA
	
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_INPUT_FILES, ConstantsTestsCase.MANAGING_TEMPLATE_INPUT_DATA_CSV)
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
		only_detect_errors = False
		b_test_char_encoding = True
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, only_detect_errors)
		self.assertEquals(0, parse_in_files.get_errors().get_len_vect_results())
		self.assertEquals(2, len(parse_in_files.get_vect_samples()))
		self.assertEquals('xpto1', parse_in_files.get_vect_samples()[0][0].name)
		self.assertEquals('asds.fastq.gz', parse_in_files.get_vect_samples()[0][0].file_name_1)
		self.assertEquals('asdws.fastq.gz', parse_in_files.get_vect_samples()[0][0].file_name_2)
		self.assertEquals('Generic', parse_in_files.get_vect_samples()[0][0].data_set.name)
		self.assertEquals(12, parse_in_files.get_vect_samples()[0][0].week)
		self.assertEquals('SRID=4326;POINT (12.2 14.3)', str(parse_in_files.get_vect_samples()[0][0].geo_local))
		self.assertEquals('vaccine status', parse_in_files.get_vect_samples()[0][0].vaccine_status.name)
		self.assertEquals('2017-12-13 00:00:00', str(parse_in_files.get_vect_samples()[0][0].date_of_collection))
		self.assertEquals('2017-12-12 00:00:00', str(parse_in_files.get_vect_samples()[0][0].date_of_onset))
		self.assertEquals(2, len(parse_in_files.get_vect_samples()[0][1]))
		self.assertEquals('other fields from here', parse_in_files.get_vect_samples()[0][1][0].value)
		self.assertEquals('other fields from here', parse_in_files.get_vect_samples()[0][1][0].tag_name.name)
		self.assertEquals('other fields from here11', parse_in_files.get_vect_samples()[0][1][1].value)
		self.assertEquals('other fields from here1', parse_in_files.get_vect_samples()[0][1][1].tag_name.name)
		self.assertEquals('xpto2', parse_in_files.get_vect_samples()[1][0].name)
		self.assertEquals('2asds.fastq.gz', parse_in_files.get_vect_samples()[1][0].file_name_1)
		self.assertEquals('2asds3.fastq.gz', parse_in_files.get_vect_samples()[1][0].file_name_2)
		self.assertEquals(11, parse_in_files.get_vect_samples()[1][0].week)
		self.assertEquals('data set', parse_in_files.get_vect_samples()[1][0].data_set.name)
		self.assertEquals(None, parse_in_files.get_vect_samples()[1][0].vaccine_status)
		self.assertEquals('2017-12-23 00:00:00', str(parse_in_files.get_vect_samples()[1][0].date_of_collection))
		self.assertEquals('2017-12-22 00:00:00', str(parse_in_files.get_vect_samples()[1][0].date_of_onset))
		self.assertEquals(1, len(parse_in_files.get_vect_samples()[1][1]))
		self.assertEquals('other fields from here2', parse_in_files.get_vect_samples()[1][1][0].value)
		self.assertEquals('other fields from here', parse_in_files.get_vect_samples()[1][1][0].tag_name.name)
		self.assertEquals('SRID=4326;POINT (30.2 43.5)', str(parse_in_files.get_vect_samples()[1][0].geo_local))



	def test_parse_abricate_file_5_fail(self):
		"""
		Test input files
		"""
		## create an index file from 
		#MANAGING_TEMPLATE_INPUT = "template_input.tsv"
		#MANAGING_TEMPLATE_INPUT_FAIL_HEADER = "template_input_fail_header.tsv"
		#MANAGING_TEMPLATE_INPUT_DATA
	
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_INPUT_FILES, ConstantsTestsCase.MANAGING_TEMPLATE_INPUT_DATA_FAIL)
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
		only_detect_errors = False
		b_test_char_encoding = True
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, only_detect_errors)
		self.assertEquals(6, parse_in_files.get_errors().get_len_vect_results())
		self.assertEquals("Error - 'longitude' must have values between -180<long<180. Line: 7 Column: 3\n" +\
						"Error - Sample name 'xpto1' is repeated in the file. Line: 8 Column: 1\n" +\
						"Error - File 'asdws.fastq.gz' is repeated in the samples file. Line: 8 Column: 2\n" +\
						"Error - 'week' must be integer. Line: 8 Column: 3\n" +\
						"Error - The 'onset date' must have this format DD/MM/YYYY. Line: 8 Column: 7\n" +\
						"Error - 'latitude' must have values between -90<lat<90. Line: 8 Column: 3",\
						str(parse_in_files.get_errors()))
		self.assertEquals(0, len(parse_in_files.get_vect_samples()))
		
	
	def test_parse_abricate_file_5(self):
		"""
		Test input files
		"""
		## create an index file from 
		#MANAGING_TEMPLATE_INPUT = "template_input.tsv"
		#MANAGING_TEMPLATE_INPUT_FAIL_HEADER = "template_input_fail_header.tsv"
		#MANAGING_TEMPLATE_INPUT_DATA
	
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_INPUT_FILES, ConstantsTestsCase.MANAGING_TEMPLATE_INPUT_DATA_CSV)
		self.assertTrue(os.path.exists(txt_file))
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.set_password(ConstantsTestsCase.TEST_USER_NAME)
			user.save()
		
		try:
			type_file_fastaq = MetaKey.objects.get(name=TypeFile.TYPE_FILE_fastq_gz)
		except MetaKey.DoesNotExist:
			type_file_fastaq = MetaKey()
			type_file_fastaq.name = TypeFile.TYPE_FILE_fastq_gz
			type_file_fastaq.save()
			
		upload_files = UploadFiles()
		upload_files.file_name = 'asds.fastq.gz'
		upload_files.owner = user
		upload_files.is_processed = False
		upload_files.type_file = type_file_fastaq
		upload_files.save()
		
		upload_files = UploadFiles()
		upload_files.file_name = '2asds.fastq.gz'
		upload_files.owner = user
		upload_files.is_processed = True
		upload_files.type_file = type_file_fastaq
		upload_files.save()
		
		upload_files = UploadFiles()
		upload_files.file_name = '2asds3.fastq.gz'
		upload_files.owner = user
		upload_files.is_processed = False
		upload_files.type_file = type_file_fastaq
		upload_files.save()
		
		self.assertEquals(0, UploadFiles.objects.filter(file_name__iexact='2asds3.fastq.gz', owner=user,\
								is_processed=False, type_file__name=TypeFile.TYPE_FILE_sample_file).count())
		self.assertEquals(1, UploadFiles.objects.filter(file_name__iexact='2asds3.fastq.gz', owner=user,\
								is_processed=False, type_file__name=TypeFile.TYPE_FILE_fastq_gz).count())
		
		parse_in_files = ParseInFiles()
		only_detect_errors = False
		b_test_char_encoding = True
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, only_detect_errors)
		self.assertEquals(2, parse_in_files.get_errors().get_len_vect_results())
		self.assertEquals("Error - File 'asds.fastq.gz' is repeated in the database and it's not processed yet. Line: 7 Column: 2\n" +\
						"Error - File '2asds.fastq.gz' is repeated in the database and it's not processed yet. Line: 8 Column: 3",\
						str(parse_in_files.get_errors()))

