'''
Created on Oct 28, 2017

@author: mmp
'''
from django.test import TestCase
from constants.constantsTestsCase import ConstantsTestsCase
from constants.constants import TypeFile
from utils.parse_in_files import ParseInFiles, UploadFilesByDjangoQ
from django.contrib.auth.models import User
from managing_files.models import UploadFiles, MetaKey, Sample
from django.conf import settings
from utils.utils import Utils
from django.test.utils import override_settings
import os

class Test(TestCase):

	### static
	utils = Utils()
	
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
		b_test_char_encoding = False
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, ParseInFiles.STATE_READ_all)
		self.assertEquals(1, parse_in_files.get_errors().get_len_vect_results())
		self.assertEquals("'utf-8' codec can't decode byte 0xff in position 0: invalid start byte", parse_in_files.get_errors().get_error(0).message)
		self.assertFalse(parse_in_files.get_errors().get_error(0).is_success())
		self.assertEquals(0, len(parse_in_files.get_vect_samples()))

	def test_parse_abricate_file_error(self):
		"""
		Test input files
		"""
		## create an index file from 
		
		#MANAGING_TEMPLATE_INPUT = "template_input.tsv"
		#MANAGING_TEMPLATE_INPUT_FAIL_HEADER = "template_input_fail_header.tsv"
		#MANAGING_TEMPLATE_INPUT_DATA
	
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_INPUT_FILES, ConstantsTestsCase.MANAGING_TEMPLATE_INPUT_error)
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
		self.assertEquals("csv or tsv file. Could not determine delimiter", parse_in_files.get_errors().get_error(0).message)
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
		b_test_char_encoding = True
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, ParseInFiles.STATE_READ_all)
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
		b_test_char_encoding = True
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, ParseInFiles.STATE_READ_all)
		self.assertEquals(1, parse_in_files.get_errors().get_len_vect_results())
		self.assertEquals("Header not found in the file. Please, check the names in the header, must be equal and have the same order of the template file.", parse_in_files.get_errors().get_error(0).message)
		self.assertFalse(parse_in_files.get_errors().get_error(0).is_success())
		self.assertEquals(0, len(parse_in_files.get_vect_samples()))
		self.assertEquals(0, parse_in_files.get_number_samples())

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
		b_test_char_encoding = True
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, ParseInFiles.STATE_READ_all)
		self.assertEquals(0, parse_in_files.get_errors().get_len_vect_results())
		self.assertEquals(2, len(parse_in_files.get_vect_samples()))
		self.assertEquals(2, parse_in_files.get_number_samples())
		self.assertEquals('xpto1', parse_in_files.get_vect_samples()[0][0].name)
		self.assertEquals(None, parse_in_files.get_vect_samples()[0][0].file_name_1)
		self.assertEquals('asds.fastq.gz', parse_in_files.get_vect_samples()[0][0].candidate_file_name_1)
		self.assertEquals(None, parse_in_files.get_vect_samples()[0][0].file_name_2)
		self.assertEquals('asdws.fastq.gz', parse_in_files.get_vect_samples()[0][0].candidate_file_name_2)
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
		self.assertEquals(None, parse_in_files.get_vect_samples()[1][0].file_name_1)
		self.assertEquals('2asds.fastq.gz', parse_in_files.get_vect_samples()[1][0].candidate_file_name_1)
		self.assertEquals(None, parse_in_files.get_vect_samples()[1][0].file_name_2)
		self.assertEquals('2asds3.fastq.gz', parse_in_files.get_vect_samples()[1][0].candidate_file_name_2)
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
		b_test_char_encoding = True
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, ParseInFiles.STATE_READ_all)
		self.assertEquals(0, parse_in_files.get_errors().get_len_vect_results())
		self.assertEquals(2, len(parse_in_files.get_vect_samples()))
		self.assertEquals('xpto1', parse_in_files.get_vect_samples()[0][0].name)
		self.assertEquals('asds.fastq.gz', parse_in_files.get_vect_samples()[0][0].candidate_file_name_1)
		self.assertEquals('asdws.fastq.gz', parse_in_files.get_vect_samples()[0][0].candidate_file_name_2)
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
		self.assertEquals('2asds.fastq.gz', parse_in_files.get_vect_samples()[1][0].candidate_file_name_1)
		self.assertEquals('2asds3.fastq.gz', parse_in_files.get_vect_samples()[1][0].candidate_file_name_2)
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
		b_test_char_encoding = True
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, ParseInFiles.STATE_READ_all)
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
		
		UploadFiles.objects.all().delete()
		Sample.objects.all().delete()
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
		b_test_char_encoding = True
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, ParseInFiles.STATE_READ_all)
		self.assertEquals(0, parse_in_files.get_errors().get_len_vect_results())
		self.assertEquals(2, len(parse_in_files.get_vect_samples()))
# 		self.assertEquals("Error - File 'asds.fastq.gz' is repeated in the database and it's not processed yet. Line: 7 Column: 2\n" +\
# 						"Error - File '2asds3.fastq.gz' is repeated in the database and it's not processed yet. Line: 8 Column: 3",\
# 						str(parse_in_files.get_errors()))
		
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, ParseInFiles.STATE_READ_dont_detect_errors)
		self.assertEquals(0, parse_in_files.get_errors().get_len_vect_results())
		self.assertEquals(2, len(parse_in_files.get_vect_samples()))

		
	def test_parse_abricate_file_6(self):
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
		b_test_char_encoding = True
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, ParseInFiles.STATE_READ_all)
		self.assertEquals(0, parse_in_files.get_errors().get_len_vect_results())
		self.assertEquals(2, len(parse_in_files.get_vect_samples()))
		
		try:
			type_file_sample_files = MetaKey.objects.get(name=TypeFile.TYPE_FILE_sample_file)
		except MetaKey.DoesNotExist:
			type_file_sample_files = MetaKey()
			type_file_sample_files.name = TypeFile.TYPE_FILE_sample_file
			type_file_sample_files.save()
			
		### delete some upload files that can exist
		for upload_file in UploadFiles.objects.filter(owner=user):
			upload_file.delete()
			
		self.assertFalse(parse_in_files.has_samples_files_to_process(user))
		self.assertEquals(None, parse_in_files.get_upload_samples_file(user))
		upload_files = UploadFiles()
		upload_files.file_name = os.path.basename(txt_file)
		upload_files.owner = user
		upload_files.is_processed = False
		upload_files.is_valid = True
		upload_files.type_file = type_file_sample_files
		upload_files.save()

		self.assertTrue(parse_in_files.has_samples_files_to_process(user))
		
		### after read you can create samples
		parse_in_files.create_samples(upload_files, user)
		
		## start testing
		try:
			upload_files = UploadFiles.objects.get(file_name=os.path.basename(txt_file))
		except UploadFiles.DoesNotExist as e:
			self.fail("must exist")
			
		## count all samples
		self.assertEquals(2, upload_files.samples.all().count())
		self.assertFalse(False, upload_files.is_processed)
		self.assertEquals(2, upload_files.number_files_to_process)
		self.assertEquals(0, upload_files.number_files_processed)
		self.assertTrue(True, upload_files.is_valid)

		try:
			sample = Sample.objects.get(name='xpto1', owner=user)
			self.assertEquals(2, sample.tag_names.all().count())
			self.assertEquals(12, sample.week)
			self.assertEquals("asds.fastq.gz", sample.candidate_file_name_1)
			self.assertEquals("asdws.fastq.gz", sample.candidate_file_name_2)
			self.assertEquals("vaccine status", sample.vaccine_status.name)
			self.assertEquals("Generic", sample.data_set.name)
			self.assertEquals(False, sample.has_files)
			sample.delete()
		except Sample.DoesNotExist as e:
			self.fail("must exist")

		try:
			sample = Sample.objects.get(name='xpto2', owner=user)
			self.assertEquals(11, sample.week)
			self.assertEquals("2asds.fastq.gz", sample.candidate_file_name_1)
			self.assertEquals("2asds3.fastq.gz", sample.candidate_file_name_2)
			self.assertEquals(None, sample.vaccine_status)
			self.assertEquals("data set", sample.data_set.name)
			self.assertEquals(False, sample.has_files)
			sample.delete()
		except Sample.DoesNotExist as e:
			self.fail("must exist")
			
		try:
			sample = Sample.objects.get(name='assasaxpto2', owner=user)
			self.fail("must not exist")
		except Sample.DoesNotExist as e:
			pass
		
	def test_create_samples(self):
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
		b_test_char_encoding = True
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, ParseInFiles.STATE_READ_all)
		self.assertEquals(0, parse_in_files.get_errors().get_len_vect_results())
		self.assertEquals(2, len(parse_in_files.get_vect_samples()))
		
		try:
			type_file_sample_files = MetaKey.objects.get(name=TypeFile.TYPE_FILE_sample_file)
		except MetaKey.DoesNotExist:
			type_file_sample_files = MetaKey()
			type_file_sample_files.name = TypeFile.TYPE_FILE_sample_file
			type_file_sample_files.save()
			
		upload_files = UploadFiles()
		upload_files.file_name = os.path.basename(txt_file)
		upload_files.owner = user
		upload_files.is_processed = False
		upload_files.is_valid = True
		upload_files.type_file = type_file_sample_files
		upload_files.save()

		try:
			sample = Sample.objects.get(name='xpto1', owner=user)
			self.fail("must not exist")
		except Sample.DoesNotExist as e:
			pass

		try:
			sample = Sample.objects.get(name='xpto2', owner=user)
			self.fail("must not exist")
		except Sample.DoesNotExist as e:
			pass
			
		sample = Sample()
		sample.name = 'xpto1'
		sample.owner = user
		sample.save()
		
		## set relation
		upload_files.samples.add(sample)
		
		### after read you can create samples
		parse_in_files.create_samples(upload_files, user)
		
		## start testing
		try:
			upload_files = UploadFiles.objects.get(file_name=os.path.basename(txt_file))
		except UploadFiles.DoesNotExist as e:
			self.fail("must exist")
			
		## count all samples
		self.assertEquals(2, upload_files.samples.all().count())
		self.assertFalse(False, upload_files.is_processed)
		self.assertEquals(2, upload_files.number_files_to_process)
		self.assertEquals(0, upload_files.number_files_processed)


	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_link_files(self):
		"""
		link the files from samples
		"""
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT_TEST", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))
		
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_INPUT_FILES, ConstantsTestsCase.MANAGING_TEMPLATE_INPUT_DATA_2_CSV)
		file_1_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		file_1_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)
		file_2_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ2_1)
		file_2_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ2_2)
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
		b_test_char_encoding = True
		parse_in_files.parse_sample_files(txt_file, user, b_test_char_encoding, ParseInFiles.STATE_READ_all)
		
		try:
			type_file_sample_files = MetaKey.objects.get(name=TypeFile.TYPE_FILE_sample_file)
		except MetaKey.DoesNotExist:
			type_file_sample_files = MetaKey()
			type_file_sample_files.name = TypeFile.TYPE_FILE_sample_file
			type_file_sample_files.save()

		try:
			type_file_fastaq = MetaKey.objects.get(name=TypeFile.TYPE_FILE_fastq_gz)
		except MetaKey.DoesNotExist:
			type_file_fastaq = MetaKey()
			type_file_fastaq.name = TypeFile.TYPE_FILE_fastq_gz
			type_file_fastaq.save()
		
		UploadFiles.objects.all().delete()
		Sample.objects.all().delete()	
		upload_files_1_1 = self.temp_upload_files(user, type_file_fastaq, file_1_1)
		upload_files_1_2 = self.temp_upload_files(user, type_file_fastaq, file_1_2)
		upload_files_2_1 = self.temp_upload_files(user, type_file_fastaq, file_2_1)
		upload_files_2_2 = self.temp_upload_files(user, type_file_fastaq, file_2_2)
		self.assertFalse(upload_files_2_2.is_processed)
		self.assertEquals(0, upload_files_2_2.number_files_processed)
			
		## sample files
		upload_files = UploadFiles()
		upload_files.file_name = os.path.basename(txt_file)
		upload_files.owner = user
		upload_files.is_processed = False
		upload_files.is_valid = True
		upload_files.type_file = type_file_sample_files
		upload_files.save()

		self.assertTrue(parse_in_files.has_samples_files_to_process(user))
		self.assertEquals(upload_files.id, parse_in_files.get_upload_samples_file(user).id)
		### after read you can create samples
		parse_in_files.create_samples(upload_files, user)

		## link the files
		b_testing = True
		parse_in_files.link_files(user, b_testing)
		
		###
		try:
			upload_files = UploadFiles.objects.get(id=upload_files.id)
			self.assertTrue(upload_files.is_processed)
			self.assertEquals(2, upload_files.number_files_processed)
			self.assertEquals(2, upload_files.number_files_to_process)
		except UploadFiles.DoesNotExist as e:
			self.fail("must exist")

		###
		try:
			upload_files_temp = UploadFiles.objects.get(id=upload_files_1_1.id)
			self.assertTrue(upload_files_temp.is_processed)
			self.assertEquals(1, upload_files_temp.number_files_processed)
			self.assertEquals(upload_files_temp.upload_file.id, upload_files.id)
			self.assertEquals(1, upload_files_temp.samples.all().count())
		except UploadFiles.DoesNotExist as e:
			self.fail("must exist")
			
		###
		try:
			upload_files_temp = UploadFiles.objects.get(id=upload_files_1_2.id)
			self.assertTrue(upload_files_temp.is_processed)
			self.assertEquals(1, upload_files_temp.number_files_processed)
			self.assertEquals(upload_files_temp.upload_file.id, upload_files.id)
			self.assertEquals(1, upload_files_temp.samples.all().count())
		except UploadFiles.DoesNotExist as e:
			self.fail("must exist")
		
		###
		try:
			upload_files_temp = UploadFiles.objects.get(id=upload_files_2_1.id)
			self.assertTrue(upload_files_temp.is_processed)
			self.assertEquals(1, upload_files_temp.number_files_processed)
			self.assertEquals(upload_files_temp.upload_file.id, upload_files.id)
			self.assertEquals(1, upload_files_temp.samples.all().count())
		except UploadFiles.DoesNotExist as e:
			self.fail("must exist")
		
		###
		try:
			upload_files_temp = UploadFiles.objects.get(id=upload_files_2_2.id)
			self.assertTrue(upload_files_temp.is_processed)
			self.assertTrue(upload_files_temp.is_valid)
			self.assertEquals('fastq.gz', upload_files_temp.type_file.name)
			self.assertEquals(1, upload_files_temp.number_files_processed)
			self.assertEquals(upload_files_temp.upload_file.id, upload_files.id)
			self.assertEquals(1, upload_files_temp.samples.all().count())
		except UploadFiles.DoesNotExist as e:
			self.fail("must exist")
			
				
		try:
			sample = Sample.objects.get(name='xpto12', owner=user)
			self.assertEquals(2, sample.tag_names.all().count())
			self.assertEquals(12, sample.week)
			self.assertEquals("EVA001_S66_L001_R1_001.fastq.gz", sample.candidate_file_name_1)
			self.assertEquals("EVA001_S66_L001_R2_001.fastq.gz", sample.candidate_file_name_2)
			self.assertEquals("vaccine status", sample.vaccine_status.name)
			self.assertEquals("Generic", sample.data_set.name)
			self.assertEquals(True, sample.has_files)
			sample.delete()
		except Sample.DoesNotExist as e:
			self.fail("must exist")

		try:
			sample = Sample.objects.get(name='xpto23', owner=user)
			self.assertEquals(11, sample.week)
			self.assertEquals("EVA002_S52_L001_R1_001.fastq.gz", sample.candidate_file_name_1)
			self.assertEquals("EVA002_S52_L001_R2_001.fastq.gz", sample.candidate_file_name_2)
			self.assertEquals(None, sample.vaccine_status)
			self.assertEquals("data set", sample.data_set.name)
			self.assertEquals(True, sample.has_files)
			sample.delete()
		except Sample.DoesNotExist as e:
			self.fail("must exist")

	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_read_sample_and_link_files(self):
		"""
		link the files from samples
		"""
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT_TEST", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))
		
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_INPUT_FILES, ConstantsTestsCase.MANAGING_TEMPLATE_INPUT_DATA_2_CSV)
		file_1_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		file_1_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)
		file_2_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ2_1)
		file_2_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ2_2)
		self.assertTrue(os.path.exists(txt_file))
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.set_password(ConstantsTestsCase.TEST_USER_NAME)
			user.save()
		
		sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), self.utils.get_path_upload_file(user.id,\
													TypeFile.TYPE_FILE_sample_file), os.path.basename(txt_file))
		sz_file_to = self.utils.get_unique_file(sz_file_to)		## get unique file name, user can upload files with same name...
		self.utils.copy_file(txt_file, sz_file_to)
		
		try:
			type_file_sample_files = MetaKey.objects.get(name=TypeFile.TYPE_FILE_sample_file)
		except MetaKey.DoesNotExist:
			type_file_sample_files = MetaKey()
			type_file_sample_files.name = TypeFile.TYPE_FILE_sample_file
			type_file_sample_files.save()

		try:
			type_file_fastaq = MetaKey.objects.get(name=TypeFile.TYPE_FILE_fastq_gz)
		except MetaKey.DoesNotExist:
			type_file_fastaq = MetaKey()
			type_file_fastaq.name = TypeFile.TYPE_FILE_fastq_gz
			type_file_fastaq.save()
		
		UploadFiles.objects.all().delete()
		Sample.objects.all().delete()
		upload_files_1_1 = self.temp_upload_files(user, type_file_fastaq, file_1_1)
		upload_files_1_2 = self.temp_upload_files(user, type_file_fastaq, file_1_2)
		upload_files_2_1 = self.temp_upload_files(user, type_file_fastaq, file_2_1)
		upload_files_2_2 = self.temp_upload_files(user, type_file_fastaq, file_2_2)
		self.assertFalse(upload_files_2_2.is_processed)
		self.assertEquals(0, upload_files_2_2.number_files_processed)
			
		## sample files
		upload_files = UploadFiles()
		upload_files.file_name = os.path.basename(txt_file)
		upload_files.owner = user
		upload_files.is_processed = False
		upload_files.is_valid = True
		upload_files.type_file = type_file_sample_files
		upload_files.path_name.name = sz_file_to
		upload_files.save()

		upload_files_by_djangoq = UploadFilesByDjangoQ()
		b_testing = True
		self.assertTrue(upload_files_by_djangoq.read_sample_file(user, upload_files, b_testing))
		
		###
		try:
			upload_files = UploadFiles.objects.get(id=upload_files.id)
			self.assertTrue(upload_files.is_processed)
			self.assertEquals(2, upload_files.number_files_processed)
			self.assertEquals(2, upload_files.number_files_to_process)
		except UploadFiles.DoesNotExist as e:
			self.fail("must exist")

		###
		try:
			upload_files_temp = UploadFiles.objects.get(id=upload_files_1_1.id)
			self.assertTrue(upload_files_temp.is_processed)
			self.assertEquals(1, upload_files_temp.number_files_processed)
			self.assertEquals(upload_files_temp.upload_file.id, upload_files.id)
			self.assertEquals(1, upload_files_temp.samples.all().count())
		except UploadFiles.DoesNotExist as e:
			self.fail("must exist")
			
		###
		try:
			upload_files_temp = UploadFiles.objects.get(id=upload_files_1_2.id)
			self.assertTrue(upload_files_temp.is_processed)
			self.assertEquals(1, upload_files_temp.number_files_processed)
			self.assertEquals(upload_files_temp.upload_file.id, upload_files.id)
			self.assertEquals(1, upload_files_temp.samples.all().count())
		except UploadFiles.DoesNotExist as e:
			self.fail("must exist")
		
		###
		try:
			upload_files_temp = UploadFiles.objects.get(id=upload_files_2_1.id)
			self.assertTrue(upload_files_temp.is_processed)
			self.assertEquals(1, upload_files_temp.number_files_processed)
			self.assertEquals(upload_files_temp.upload_file.id, upload_files.id)
			self.assertEquals(1, upload_files_temp.samples.all().count())
		except UploadFiles.DoesNotExist as e:
			self.fail("must exist")
		
		###
		try:
			upload_files_temp = UploadFiles.objects.get(id=upload_files_2_2.id)
			self.assertTrue(upload_files_temp.is_processed)
			self.assertTrue(upload_files_temp.is_valid)
			self.assertEquals('fastq.gz', upload_files_temp.type_file.name)
			self.assertEquals(1, upload_files_temp.number_files_processed)
			self.assertEquals(upload_files_temp.upload_file.id, upload_files.id)
			self.assertEquals(1, upload_files_temp.samples.all().count())
		except UploadFiles.DoesNotExist as e:
			self.fail("must exist")
			
				
		try:
			sample = Sample.objects.get(name='xpto12', owner=user)
			self.assertEquals(2, sample.tag_names.all().count())
			self.assertEquals(12, sample.week)
			self.assertEquals("EVA001_S66_L001_R1_001.fastq.gz", sample.candidate_file_name_1)
			self.assertEquals("EVA001_S66_L001_R2_001.fastq.gz", sample.candidate_file_name_2)
			self.assertEquals("vaccine status", sample.vaccine_status.name)
			self.assertEquals("Generic", sample.data_set.name)
			self.assertEquals(True, sample.has_files)
			sample.delete()
		except Sample.DoesNotExist as e:
			self.fail("must exist")

		try:
			sample = Sample.objects.get(name='xpto23', owner=user)
			self.assertEquals(11, sample.week)
			self.assertEquals("EVA002_S52_L001_R1_001.fastq.gz", sample.candidate_file_name_1)
			self.assertEquals("EVA002_S52_L001_R2_001.fastq.gz", sample.candidate_file_name_2)
			self.assertEquals(None, sample.vaccine_status)
			self.assertEquals("data set", sample.data_set.name)
			self.assertEquals(True, sample.has_files)
			sample.delete()
		except Sample.DoesNotExist as e:
			self.fail("must exist")


	def temp_upload_files(self, user, type_file_fastq, file_name):
		try:
			upload_files = UploadFiles.objects.get(owner=user, type_file__name=type_file_fastq.name, file_name=os.path.basename(file_name), is_processed=False)
		except UploadFiles.DoesNotExist as e:
			upload_files = UploadFiles()
			upload_files.file_name = os.path.basename(file_name)
			upload_files.owner = user
			upload_files.is_processed = False
			upload_files.is_valid = True
			upload_files.type_file = type_file_fastq
			upload_files.path_name.name = file_name
			upload_files.save()
		return upload_files


