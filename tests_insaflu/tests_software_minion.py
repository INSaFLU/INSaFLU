'''
Created on 01/01/2021

@author: mmp
'''
import os
from django.test import TestCase
from django.conf import settings
from constants.constantsTestsCase import ConstantsTestsCase
from django.contrib.auth.models import User
from managing_files.models import Sample
from constants.constants import Constants, TypePath
from managing_files.manage_database import ManageDatabase
from utils.result import DecodeObjects
from constants.meta_key_and_values import MetaKeyAndValue
from utils.software_minion import SoftwareMinion
from utils.utils import Utils
from constants.software_names import SoftwareNames

class Test(TestCase):

	software_minion = SoftwareMinion()
	utils = Utils()
	constants = Constants()
	software_names = SoftwareNames()
	constants_tests_case = ConstantsTestsCase()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)

	def tearDown(self):
		pass


	def test_run_fastq_and_trimmomatic(self):
		"""
		Test run fastq and trimmomatic all together
 		"""
		file_name = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ_MINION_1)
		self.assertTrue(file_name)
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_name, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
			
		sample_name = "run_nanofilt"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = False
			sample.owner = user
			sample.save()
		
		### run software
		self.assertTrue(self.software_minion.run_nanofilt_and_stat(sample, user))
		
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_nanofilt_file(TypePath.MEDIA_ROOT)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_rabbitQC_output(TypePath.MEDIA_ROOT)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_rabbitQC_nanofilt(TypePath.MEDIA_ROOT)))))
		
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Number_And_Average_Reads, list_meta[0].meta_tag.name)
		
		### number of sequences
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Number_And_Average_Reads, list_meta[0].meta_tag.name)

		decodeResultAverageAndNumberReads = DecodeObjects()
		result_average = decodeResultAverageAndNumberReads.decode_result(list_meta[0].description)
		self.assertEqual('969', result_average.number_file_1)
		self.assertEqual('468.1', result_average.average_file_1)
		self.assertEqual(None, result_average.number_file_2)
		self.assertEqual(None, result_average.average_file_2)

		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, list_meta[0].meta_tag.name)
		self.assertEquals("Success, NanoStat(1.4.0), NanoFilt(2.6.0)", list_meta[0].description)
		
		## remove all files
		cmd = "rm -r %s*" % (temp_dir); os.system(cmd)
