'''
Created on 01/01/2021

@author: mmp
'''
import os, filecmp
from django.test import TestCase
from django.conf import settings
from constants.constantsTestsCase import ConstantsTestsCase
from django.contrib.auth.models import User
from managing_files.models import Sample, Reference, Project, ProjectSample, MixedInfectionsTag
from constants.constants import Constants, TypePath, FileType, FileExtensions
from managing_files.manage_database import ManageDatabase
from utils.result import DecodeObjects, Coverage, MaskingConsensus
from constants.meta_key_and_values import MetaKeyAndValue
from utils.software_minion import SoftwareMinion
from utils.software import Software
from utils.utils import Utils
from constants.software_names import SoftwareNames
from django.test.utils import override_settings
from settings.models import Parameter
from settings.models import Software as SoftwareSettings
from settings.default_software_project_sample import DefaultProjectSoftware
from manage_virus.uploadFiles import UploadFiles
from settings.constants_settings import ConstantsSettings
from settings.default_software import DefaultSoftware

class Test(TestCase):

	software_minion = SoftwareMinion()
	software = Software()
	utils = Utils()
	constants = Constants()
	software_names = SoftwareNames()
	constants_tests_case = ConstantsTestsCase()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)

	def tearDown(self):
		pass


	def test_run_clean_minion(self):
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
			sample.type_of_fastq = Sample.TYPE_OF_FASTQ_minion
			sample.owner = user
			sample.save()
		
		### run software
		self.assertTrue(self.software_minion.run_nanofilt_and_stat(sample, user))
		
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_nanofilt_file(TypePath.MEDIA_ROOT)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_rabbitQC_output(TypePath.MEDIA_ROOT)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_rabbitQC_nanofilt(TypePath.MEDIA_ROOT)))))
		
		### number of sequences
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Number_And_Average_Reads, list_meta[0].meta_tag.name)

		decodeResultAverageAndNumberReads = DecodeObjects()
		result_average = decodeResultAverageAndNumberReads.decode_result(list_meta[0].description)
		self.assertEqual('969', result_average.number_file_1)
		self.assertEqual('348.1', result_average.average_file_1)
		self.assertEqual(None, result_average.number_file_2)
		self.assertEqual(None, result_average.average_file_2)

		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
		self.assertTrue(len(list_meta) == 0)
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, list_meta[0].meta_tag.name)
		self.assertEquals("Success, NanoStat(1.4.0), NanoFilt(2.6.0)", list_meta[0].description)
		
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, list_meta[0].meta_tag.name)
		result_data = decodeResultAverageAndNumberReads.decode_result(list_meta[0].description)
		vect_soft = result_data.get_list_software_instance(self.software_names.get_NanoStat_name())
		self.assertEqual(2, len(vect_soft))
		
		self.assertEqual("Mean read length", vect_soft[0].get_vect_key_values()[0].key)
		self.assertEqual("489.3", vect_soft[0].get_vect_key_values()[0].value)
		self.assertEqual("Total bases", vect_soft[0].get_vect_key_values()[-1].key)
		self.assertEqual("517,186.0", vect_soft[0].get_vect_key_values()[-1].value)
		self.assertEqual("Mean read length", vect_soft[1].get_vect_key_values()[0].key)
		self.assertEqual("348.1", vect_soft[1].get_vect_key_values()[0].value)
		self.assertEqual("Total bases", vect_soft[1].get_vect_key_values()[-1].key)
		self.assertEqual("337,289.0", vect_soft[1].get_vect_key_values()[-1].value)
		
		### software
		meta_key_data = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, None)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, meta_key_data.value)
		self.assertEquals(MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, meta_key_data.meta_tag.name)
		
		decode_result = DecodeObjects()
		result = decode_result.decode_result(meta_key_data.description)
		self.assertTrue(not result is None)
		self.assertEquals("NanoStat-1.4.0",
				result.get_software(self.software_names.get_NanoStat_name()))
		self.assertEquals("NanoFilt-2.6.0; (-q 10 -l 50 --headcrop 70 --tailcrop 70)",\
 						result.get_software(self.software_names.get_NanoFilt_name()))
		self.assertEquals("RabbitQC-0.0.1; (-w 3 -D)",\
 						result.get_software(self.software_names.get_rabbitQC_name()))

		
		## remove all files
		cmd = "rm -r %s*" % (temp_dir); os.system(cmd)
		
	def test_run_clean_minion_2(self):
		"""
		Test run fastq and trimmomatic all together
 		"""
		uploadFiles = UploadFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		self.assertEqual("2", version)
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
		uploadFiles.upload_file(version, file)	## upload file
		
		file_name = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, "covid/test_minion_seq.fastq.gz")
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
			sample.type_of_fastq = Sample.TYPE_OF_FASTQ_minion
			sample.owner = user
			sample.save()
		
		### turn off abricate
		default_software = DefaultSoftware()
		default_software.test_all_defaults(user)
		self.assertTrue(default_software.is_software_to_run(SoftwareNames.SOFTWARE_ABRICATE_name,
							user, ConstantsSettings.TECHNOLOGY_minion))
		is_to_run = False
		default_software.set_software_to_run(SoftwareNames.SOFTWARE_ABRICATE_name, user,
							ConstantsSettings.TECHNOLOGY_minion, is_to_run)
		self.assertFalse(default_software.is_software_to_run(SoftwareNames.SOFTWARE_ABRICATE_name,
							user, ConstantsSettings.TECHNOLOGY_minion))

		b_make_identify_species = False
		self.assertTrue(self.software_minion.run_clean_minion(sample, user, b_make_identify_species))
		default_project_software = DefaultProjectSoftware()
		self.assertFalse(default_project_software.is_to_run_abricate(user, sample,
				ConstantsSettings.TECHNOLOGY_minion))
		
		try:
			software = SoftwareSettings.objects.get(name=SoftwareNames.SOFTWARE_NanoFilt_name, owner=user,\
							type_of_use = SoftwareSettings.TYPE_OF_USE_sample,
							technology__name = ConstantsSettings.TECHNOLOGY_minion,
							is_obsolete = False)
		except SoftwareSettings.DoesNotExist:
			self.fail("Must must exist")
		
		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(5, len(parameters))
		for parameter in parameters:
			if (parameter.name == "-q"):
				parameter.parameter = "5"
				parameter.save()
				break
		
		tag_name_mixed= "No"
		sample = Sample.objects.get(name=sample_name)
		if (b_make_identify_species):
			self.fail("Must have b_make_identify_species FALSE")
			self.assertEqual(tag_name_mixed, sample.mixed_infections_tag.name)
			self.assertEqual("A-H5N5", sample.type_subtype)
			
			sample.mixed_infections_tag = None
			sample.save()
			
			sample = Sample.objects.get(name=sample_name)
			self.assertEqual(None, sample.mixed_infections_tag)
			
			try:
				mixed_infections_tag = MixedInfectionsTag.objects.get(name=tag_name_mixed)
			except MixedInfectionsTag.DoesNotExist as e:
				self.fail("Must have tag mixed")
							
			if (not sample.mixed_infections_tag is None): sample.mixed_infections_tag.delete()
		else:
			self.assertEqual(Constants.EMPTY_VALUE_NA, sample.mixed_infections_tag.name)
			self.assertEqual(Constants.EMPTY_VALUE_NA, sample.type_subtype)
		
		self.assertEqual(0, sample.number_alerts)
		
		### run again
		default_project_software.set_abricate_to_run(user, sample,
				ConstantsSettings.TECHNOLOGY_minion, True)
		## global
		self.assertFalse(default_software.is_software_to_run(SoftwareNames.SOFTWARE_ABRICATE_name,
							user, ConstantsSettings.TECHNOLOGY_minion))
		### project sample
		self.assertTrue(default_project_software.is_to_run_abricate(user, 
					sample, ConstantsSettings.TECHNOLOGY_minion))
		self.assertTrue(self.software_minion.run_clean_minion(sample, user,
					default_project_software.is_to_run_abricate(user,
					sample, ConstantsSettings.TECHNOLOGY_minion) ))
		
		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(5, len(parameters))
		for parameter in parameters:
			if (parameter.name == "-q"):
				self.assertEqual("5", parameter.parameter)
				break;

		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_nanofilt_file(TypePath.MEDIA_ROOT)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_rabbitQC_output(TypePath.MEDIA_ROOT)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_rabbitQC_nanofilt(TypePath.MEDIA_ROOT)))))
		
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None)
		self.assertEqual(2, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Number_And_Average_Reads, list_meta[0].meta_tag.name)
		
		decodeResultAverageAndNumberReads = DecodeObjects()
		result_average = decodeResultAverageAndNumberReads.decode_result(list_meta[0].description)
		self.assertEqual('3000', result_average.number_file_1)
		self.assertEqual('1158.4', result_average.average_file_1)
		self.assertEqual(None, result_average.number_file_2)
		self.assertEqual(None, result_average.average_file_2)

		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
		self.assertTrue(len(list_meta) == 0)
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, None)
		self.assertEqual(2, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, list_meta[0].meta_tag.name)
		self.assertEquals("Success, NanoStat(1.4.0), NanoFilt(2.6.0)", list_meta[0].description)
		self.assertEquals("Success, NanoStat(1.4.0), NanoFilt(2.6.0)", list_meta[1].description)
		
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, None)
		self.assertEqual(2, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, list_meta[0].meta_tag.name)
		result_data = decodeResultAverageAndNumberReads.decode_result(list_meta[0].description)
		vect_soft = result_data.get_list_software_instance(self.software_names.get_NanoStat_name())
		self.assertEqual(2, len(vect_soft))

		self.assertEqual("Mean read length", vect_soft[0].get_vect_key_values()[0].key)
		self.assertEqual("1,298.4", vect_soft[0].get_vect_key_values()[0].value)
		self.assertEqual("Total bases", vect_soft[0].get_vect_key_values()[-1].key)
		self.assertEqual("3,895,144.0", vect_soft[0].get_vect_key_values()[-1].value)
		self.assertEqual("Mean read length", vect_soft[1].get_vect_key_values()[0].key)
		self.assertEqual("1,158.4", vect_soft[1].get_vect_key_values()[0].value)
		self.assertEqual("Total bases", vect_soft[1].get_vect_key_values()[-1].key)
		self.assertEqual("3,475,144.0", vect_soft[1].get_vect_key_values()[-1].value)
		
		stats_illumina, total_reads = self.software.get_stats_from_sample_reads(sample)
		self.assertEqual('Mean read length', stats_illumina[0][0])
		self.assertEqual('1,298.4', stats_illumina[0][1])
		self.assertEqual(3000, total_reads)
		
		sample = Sample.objects.get(name=sample_name)
		self.assertEqual("A-H5N5", sample.type_subtype)
	#	self.assertEqual("H5N5", sample.type_subtype)
		self.assertEqual("No", sample.mixed_infections_tag.name)
		self.assertEqual(0, sample.number_alerts)
		
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, None)
		self.assertEquals(2, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals("Success, NanoStat(1.4.0), NanoFilt(2.6.0)", list_meta[0].description)
		
		### software
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, None)
		self.assertEquals(2, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, list_meta[0].meta_tag.name)
		
		decode_result = DecodeObjects()
		result = decode_result.decode_result(list_meta[0].description)
		self.assertTrue(not result is None)
		self.assertEquals("NanoStat-1.4.0",
				result.get_software(self.software_names.get_NanoStat_name()))
		self.assertEquals("NanoFilt-2.6.0; (-q 5 -l 50 --headcrop 70 --tailcrop 70)",\
 						result.get_software(self.software_names.get_NanoFilt_name()))
		self.assertEquals("RabbitQC-0.0.1; (-w 3 -D)",\
 						result.get_software(self.software_names.get_rabbitQC_name()))
		
		## remove all files
		cmd = "rm -r %s*" % (temp_dir); os.system(cmd)

		
	def test_run_clean_minion_3(self):
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
			sample.type_of_fastq = Sample.TYPE_OF_FASTQ_minion
			sample.owner = user
			sample.save()
		
		### turn off nanostat
		default_software = DefaultSoftware()
		default_software.test_all_defaults(user)
		self.assertTrue(default_software.is_software_to_run(SoftwareNames.SOFTWARE_NanoFilt_name,
							user, ConstantsSettings.TECHNOLOGY_minion))
		is_to_run = False
		default_software.set_software_to_run(SoftwareNames.SOFTWARE_NanoFilt_name, user,
							ConstantsSettings.TECHNOLOGY_minion, is_to_run)
		self.assertFalse(default_software.is_software_to_run(SoftwareNames.SOFTWARE_NanoFilt_name,
							user, ConstantsSettings.TECHNOLOGY_minion))
		
		### run software
		self.assertTrue(self.software_minion.run_nanofilt_and_stat(sample, user))
		default_project_software = DefaultProjectSoftware()
		self.assertFalse(default_project_software.is_to_run_nanofilt(user, sample))
		
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
		self.assertFalse(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_nanofilt_file(TypePath.MEDIA_ROOT)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_rabbitQC_output(TypePath.MEDIA_ROOT)))))
		self.assertFalse(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_rabbitQC_nanofilt(TypePath.MEDIA_ROOT)))))
		
		### number of sequences
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Number_And_Average_Reads, list_meta[0].meta_tag.name)

		decodeResultAverageAndNumberReads = DecodeObjects()
		result_average = decodeResultAverageAndNumberReads.decode_result(list_meta[0].description)
		self.assertEqual('1057', result_average.number_file_1)
		self.assertEqual('489.3', result_average.average_file_1)
		self.assertEqual(None, result_average.number_file_2)
		self.assertEqual(None, result_average.average_file_2)

		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
		self.assertTrue(len(list_meta) == 0)
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, list_meta[0].meta_tag.name)
		self.assertEquals("Success, NanoStat(1.4.0)", list_meta[0].description)
		
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, list_meta[0].meta_tag.name)
		result_data = decodeResultAverageAndNumberReads.decode_result(list_meta[0].description)
		vect_soft = result_data.get_list_software_instance(self.software_names.get_NanoStat_name())
		self.assertEqual(1, len(vect_soft))
		
		self.assertEqual("Mean read length", vect_soft[0].get_vect_key_values()[0].key)
		self.assertEqual("489.3", vect_soft[0].get_vect_key_values()[0].value)
		self.assertEqual("Total bases", vect_soft[0].get_vect_key_values()[-1].key)
		self.assertEqual("517,186.0", vect_soft[0].get_vect_key_values()[-1].value)
		
		### software
		meta_key_data = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, None)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, meta_key_data.value)
		self.assertEquals(MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, meta_key_data.meta_tag.name)
		
		decode_result = DecodeObjects()
		result = decode_result.decode_result(meta_key_data.description)
		self.assertTrue(not result is None)
		self.assertEquals("NanoStat-1.4.0",
				result.get_software(self.software_names.get_NanoStat_name()))
		self.assertEquals("", result.get_software(self.software_names.get_NanoFilt_name()))
		self.assertEquals("RabbitQC-0.0.1; (-w 3 -D)",\
 						result.get_software(self.software_names.get_rabbitQC_name()))

		## remove all files
		cmd = "rm -r %s*" % (temp_dir); os.system(cmd)


	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_process_second_stage_medaka(self):
		"""
 		test global method
 		"""
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT_TEST", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))
		
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_COVID_GBK)
		reference_fasta = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_COVID_FASTA)
		file_fasta = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_nanopore)
		
		self.assertTrue(os.path.exists(gb_file))
		self.assertTrue(os.path.exists(reference_fasta))
		self.assertTrue(os.path.exists(file_fasta))
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		ref_name = "second_stagis_single_covid"
		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.reference_fasta.name = reference_fasta
			reference.reference_fasta_name = os.path.basename(reference_fasta)
			reference.reference_genbank.name = gb_file
			reference.reference_genbank_name = os.path.basename(gb_file)
			reference.owner = user
			reference.save()
		
		### set genetic generic element
		self.utils.get_elements_and_cds_from_db(reference, user)
		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_fasta, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
			
		sample_name = "run_ont_covid"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = False
			sample.file_name_2 = None
			sample.path_name_2.name = None
			sample.type_of_fastq = Sample.TYPE_OF_FASTQ_minion
			sample.owner = user
			sample.save()

		project_name = "file_naais_single_covid"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = project_name
			project.reference = reference
			project.owner = user
			project.save()
		
		## create project_sample
		project_sample = ProjectSample()
		project_sample.sample = sample
		project_sample.project = project
		project_sample.save()
		
		### copy to nanofilt file
		self.utils.copy_file(file_fasta, project_sample.sample.get_nanofilt_file(TypePath.MEDIA_ROOT))
		
		### set 60 to mask consensus
		default_software = DefaultProjectSoftware()
		default_software.test_all_defaults(user, project, None, None)
		try:
			software = SoftwareSettings.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, owner=user,\
							type_of_use = SoftwareSettings.TYPE_OF_USE_project,
							technology__name = ConstantsSettings.TECHNOLOGY_minion,
							is_obsolete = False)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertTrue(software.is_used_in_project())
			self.assertFalse(software.is_used_in_global())
		except SoftwareSettings.DoesNotExist:
			self.fail("Must exist this software name")
			
		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		
		### test set default
		self.assertEqual("70", parameters[0].parameter)
		parameter = parameters[0]
		parameter.parameter = "60"
		parameter.save()
		
		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		self.assertEqual("60", parameters[0].parameter)
		
		##################### limit coverage
		try:
			software = SoftwareSettings.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name, owner=user,\
							type_of_use = SoftwareSettings.TYPE_OF_USE_project,
							technology__name = ConstantsSettings.TECHNOLOGY_minion,
							is_obsolete = False)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertTrue(software.is_used_in_project())
			self.assertFalse(software.is_used_in_global())
		except SoftwareSettings.DoesNotExist:
			self.fail("Must exist this software name")
			
		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		
		### test set default
		self.assertEqual("30", parameters[0].parameter)
		parameter = parameters[0]
		parameter.parameter = "5"
		parameter.save()
		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		self.assertEqual("5", parameters[0].parameter)
		
		
		##################### set FREQ vcf limit
		try:
			software = SoftwareSettings.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name, owner=user,\
							type_of_use = SoftwareSettings.TYPE_OF_USE_project,
							technology__name = ConstantsSettings.TECHNOLOGY_minion,
							is_obsolete = False)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertTrue(software.is_used_in_project())
			self.assertFalse(software.is_used_in_global())
		except SoftwareSettings.DoesNotExist:
			self.fail("Must exist this software name")
			
		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		
		### test set default
		self.assertEqual("0.80", parameters[0].parameter)
		parameter = parameters[0]
		parameter.parameter = ".40"
		parameter.save()
		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		self.assertEqual(".40", parameters[0].parameter)
		
		## continue to process	
		manageDatabase = ManageDatabase()
		metaKeyAndValue = MetaKeyAndValue()
		meta_key_project_sample = metaKeyAndValue.get_meta_key_queue_by_project_sample_id(project_sample.id)
		manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue, "meta_sample.description")
		self.assertTrue(self.software_minion.process_second_stage_medaka(project_sample, user))
		try:
			project_sample = ProjectSample.objects.get(pk=project_sample.id)
		except ProjectSample.DoesNotExist:
			self.fail("Must exist")
		self.assertTrue(project_sample.is_finished)
		self.assertTrue(project_sample.is_mask_consensus_sequences)
		
		### test the files
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, SoftwareNames.SOFTWARE_Medaka_name)))
		self.assertFalse(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FA, SoftwareNames.SOFTWARE_Medaka_name)))
		### test sample name in the file name
		self.assertTrue(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, SoftwareNames.SOFTWARE_Medaka_name).index(sample_name) != -1)
	
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_BAM, SoftwareNames.SOFTWARE_Medaka_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_BAM_BAI, SoftwareNames.SOFTWARE_Medaka_name)))
		self.assertFalse(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH, SoftwareNames.SOFTWARE_Medaka_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ, SoftwareNames.SOFTWARE_Medaka_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ_TBI, SoftwareNames.SOFTWARE_Medaka_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_Medaka_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, SoftwareNames.SOFTWARE_Medaka_name)))
		self.assertFalse(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CSV, SoftwareNames.SOFTWARE_Medaka_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_Medaka_name)))
		self.assertFalse(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_REF_FASTA, SoftwareNames.SOFTWARE_Medaka_name)))
		
		### campare VCF
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "medaka_out.vcf")
		out_file = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_Medaka_name)
		temp_file = self.utils.get_temp_file("clean_vcf", ".txt")
		os.system('grep -v "##SnpEffCmd" {} > {}'.format(out_file, temp_file))
		self.assertTrue(os.path.exists(temp_file))
		self.assertTrue(filecmp.cmp(temp_file, expected_file_samples))
		
		### compare
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "medaka_output.tab")
		out_file = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_Medaka_name)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(os.path.exists(expected_file_samples))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		
		## test consensus file
		self.assertTrue(os.path.exists(project_sample.get_consensus_file(TypePath.MEDIA_ROOT)))
		expected_file_consensus = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "medaka_output_consensus.fasta")
		self.assertTrue(filecmp.cmp(project_sample.get_consensus_file(TypePath.MEDIA_ROOT), expected_file_consensus))
		self.assertTrue(filecmp.cmp(project_sample.get_consensus_file(TypePath.MEDIA_ROOT),
							project_sample.get_backup_consensus_file()))

		manageDatabase = ManageDatabase()
		meta_value = manageDatabase.get_project_sample_metakey_last(project_sample,
						MetaKeyAndValue.META_KEY_Masking_consensus_by_minfrac_VCF_medaka,
						MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(not meta_value is None)
		decode_masking_vcf = DecodeObjects()
		masking_consensus = decode_masking_vcf.decode_result(meta_value.description)
		self.assertTrue(masking_consensus.has_masking_data())
		mask_consensus = MaskingConsensus()
		mask_consensus.set_mask_sites("5788,11521,13642,13875,15272,21369,23816,24389,24643,29573")
		self.assertEqual(mask_consensus, masking_consensus.get_mask_consensus_element('MN908947'))
		self.assertEqual("Element:MN908947 -> (Mask sites:5788,11521,13642,13875,15272,21369,23816,24389,24643,29573)",
						masking_consensus.get_message_to_show_in_csv_file())
				
		### human file name, snippy tab
#		Check the ones with zero coverage		
#		print(project_sample.get_file_output_human(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_Medaka_name))
		self.assertTrue(os.path.exists(project_sample.get_file_output_human(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_Medaka_name)))
		self.assertTrue(os.path.getsize(project_sample.get_file_output_human(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_Medaka_name)) > 10)
		
		### set flag that is finished
		list_meta = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Count_Hits, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Count_Hits, list_meta[0].meta_tag.name)
		
		### get the hits value
		decode_coverage = DecodeObjects()
		count_hits = decode_coverage.decode_result(list_meta[0].description)
		self.assertEquals(0, count_hits.get_hits_50_90())
		self.assertEquals(3, count_hits.get_hits_less_50())
		self.assertEquals(3, count_hits.get_total_50_50_90())
		self.assertEquals(10, count_hits.get_total())
		
		self.assertEquals(3, project_sample.count_variations.var_less_50)
		self.assertEquals(0, project_sample.count_variations.var_bigger_50_90)
		self.assertEquals(7, project_sample.count_variations.var_bigger_90)
		self.assertEquals(0, project_sample.alert_first_level)
		self.assertEquals(1, project_sample.alert_second_level)
		
		### test mixed infections
		self.assertEquals(Constants.EMPTY_VALUE_NA, project_sample.mixed_infections.tag.name)
		self.assertFalse(project_sample.mixed_infections.has_master_vector)
		
		list_meta = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Medaka, None)
		self.assertEquals(1, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Medaka, list_meta[0].meta_tag.name)
		
		decode_result = DecodeObjects()
		result = decode_result.decode_result(list_meta[0].description)
		self.assertTrue(not result is None)
		self.assertEquals("MSA Masker-1.0; (--c 4; for coverages less than 5 in 40% of the regions.)", result.get_software(
			self.software_names.get_msa_masker_name()))
		self.assertEquals("Minimum depth of coverage per site to validate the sequence; (Threshold:5)", result.get_software(
			self.software_names.get_insaflu_parameter_limit_coverage_name()))
		self.assertEquals("Medaka-1.2.1; (consensus -m r941_min_high_g360)/Medaka-1.2.1; (variant --verbose)",\
 			result.get_software(self.software_names.get_medaka_name()))
		self.assertEquals("Samtools-1.3; (depth -aa)",\
 			result.get_software(self.software_names.get_samtools_name()))

		meta_value = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
		decode_coverage = DecodeObjects()
		coverage = decode_coverage.decode_result(meta_value.description)
		self.assertEqual(coverage.get_coverage('MN908947', Coverage.COVERAGE_ALL), "13.5")
		self.assertEqual(coverage.get_coverage('MN908947', Coverage.COVERAGE_MORE_9), "63.1")
		self.assertEqual(coverage.get_coverage('MN908947', Coverage.COVERAGE_MORE_0), "96.1")
		self.assertEqual(coverage.get_coverage('MN908947', Coverage.COVERAGE_PROJECT), "84.2") ### in this case is equal to 9
		
		lst_meta_sample = manageDatabase.get_project_sample_metakey(project_sample, meta_key_project_sample, None)
		self.assertEquals(2, len(lst_meta_sample))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Queue, lst_meta_sample[1].value)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, lst_meta_sample[0].value)
		
		### check if the images exist
		for gene in self.utils.get_elements_from_db(reference, user):
			output_image = project_sample.get_global_file_by_element(TypePath.MEDIA_ROOT, ProjectSample.PREFIX_FILE_COVERAGE, gene, FileExtensions.FILE_PNG)
			self.assertTrue(os.path.exists(output_image))

		meta_value = manageDatabase.get_project_sample_metakey_last(project_sample,
						MetaKeyAndValue.META_KEY_bam_stats,
						MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(not meta_value is None)
		decode_masking_vcf = DecodeObjects()
		result = decode_masking_vcf.decode_result(meta_value.description)
		self.assertEqual('1057', result.get_value_by_key(MetaKeyAndValue.SAMTOOLS_flagstat_mapped_reads))
		self.assertEqual('1057', result.get_value_by_key(MetaKeyAndValue.SAMTOOLS_flagstat_total_reads))
		self.assertEqual(None, result.get_value_by_key('xpto'))

		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))

