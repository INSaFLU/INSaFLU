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
from utils.software import Software, Contigs2Sequences
from constants.software_names import SoftwareNames
from utils.utils import Utils
from utils.parse_out_files import ParseOutFiles
from utils.result import DecodeObjects, Coverage
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference
from manage_virus.uploadFiles import UploadFiles
from manage_virus.models import UploadFile
from managing_files.manage_database import ManageDatabase
from django.test.utils import override_settings
from manage_virus.constants_virus import ConstantsVirus
import os, filecmp

class Test(TestCase):

	### static
	software = Software()
	utils = Utils()
	constants = Constants()
	software_names = SoftwareNames()
	constants_tests_case = ConstantsTestsCase()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)


	def tearDown(self):
		pass


	def testCreateFaiToFastaFile(self):
		"""
		Test samtools fai index
		"""
		## create an index file from 
		
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		temp_dir = os.path.join(self.utils.get_temp_dir())
		fasta_file_temp = self.utils.get_temp_file_from_dir(temp_dir, "FaiToFastaFile", FileExtensions.FILE_FASTA)
		self.utils.copy_file(fasta_file, fasta_file_temp)
		try:
			self.software.create_fai_fasta(fasta_file_temp)
			self.assertTrue(os.path.isfile(fasta_file_temp + ".fai"))
			self.assertEquals(179, os.path.getsize(fasta_file_temp + ".fai"))
		except Exception as e:
			self.fail(e.args[0])
		
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FAIL_FASTA)
		fasta_file_temp = self.utils.get_temp_file_from_dir(temp_dir, "FaiToFastaFile", FileExtensions.FILE_FASTA)
		self.utils.copy_file(fasta_file, fasta_file_temp)
		try:
			self.software.create_fai_fasta(fasta_file_temp)
			self.fail("Must throw exception")
		except Exception as e:
			self.assertEquals("Fail to run samtools", e.args[0])
			if os.path.isfile(fasta_file_temp + ".fai"): os.unlink(fasta_file_temp + ".fai")
		self.utils.remove_dir(temp_dir)
		
	def test_is_exist_database_abricate(self):
		self.assertTrue(self.software.is_exist_database_abricate("ncbi"))
		self.assertFalse(self.software.is_exist_database_abricate("xpot"))
	
	def test_create_database_abricate(self):
		database_name = "xpto"
		if (self.software.is_exist_database_abricate(database_name)):
			cmd = "rm -r %s/%s*" % (SoftwareNames.SOFTWARE_ABRICATE_DB, database_name)
			exist_status = os.system(cmd)
			self.assertTrue(exist_status == 0)
		self.assertFalse(self.software.is_exist_database_abricate(database_name))
		temp_file_file = os.path.join(self.baseDirectory, Constants.DIR_TYPE_IDENTIFICATION, ConstantsTestsCase.MANAGING_TEST_INFLUENZA_FILE)
		self.assertTrue(os.path.exists(temp_file_file))
		
		try:
			self.software.create_database_abricate(database_name, "fdssfd")
			self.fail("must throw error")
		except IOError:
			pass
		
		self.software.create_database_abricate(database_name, temp_file_file)
		self.assertTrue(self.software.is_exist_database_abricate(database_name))
		cmd = "rm -r %s/%s*" % (SoftwareNames.SOFTWARE_ABRICATE_DB, database_name)
		exist_status = os.system(cmd)
		self.assertTrue(exist_status == 0)
		self.assertFalse(self.software.is_exist_database_abricate(database_name))
	
	
	def test_run_spades_single(self):
		
		fastq1_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		self.assertTrue(os.path.isfile(fastq1_1))
		
		out_dir = self.utils.get_temp_dir()
		self.assertTrue(os.path.isdir(out_dir))
		cmd = self.software.run_spades(fastq1_1, "", out_dir)
		result_file = os.path.join(out_dir, "contigs.fasta")
		self.assertTrue(os.path.exists(result_file))
		self.assertTrue(os.path.getsize(result_file) > 100)
		self.utils.remove_dir(out_dir)
	
	##########################################################
	##
	##		Important, this only work from command line
	##
	## it only work from command line because of PYTHONPATH defined on the eclispe IDE 
	def testRunSpadesAndAbricate(self):
		"""
 		Test run spades and abricator
 		"""
# 		b_run_this = getattr(settings, "RUN_TEST_IN_COMMAND_LINE", None)
# 		
# 		if (not b_run_this): return	## it only work from command line because of PYTHONPATH defined on the eclispe IDE 
		fastq1_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		fastq1_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)
		self.assertTrue(os.path.isfile(fastq1_1))
		self.assertTrue(os.path.isfile(fastq1_2))
		
		out_dir = self.utils.get_temp_dir()
		self.assertTrue(os.path.isdir(out_dir))
		cmd = self.software.run_spades(fastq1_1, fastq1_2, out_dir)
		file_out = os.path.join(out_dir, "contigs.fasta")
		self.assertTrue(os.path.exists(file_out))
		self.assertTrue(os.path.getsize(file_out) > 100)
		
		database_name = "xpto"
		if (not self.software.is_exist_database_abricate(database_name)):
			temp_file_file = os.path.join(self.baseDirectory, Constants.DIR_TYPE_IDENTIFICATION, ConstantsTestsCase.MANAGING_TEST_INFLUENZA_FILE)
			self.assertTrue(os.path.exists(temp_file_file))
			self.software.create_database_abricate(database_name, temp_file_file)
		
		out_file = self.utils.get_temp_file("temp_abricate", ".txt")
		cmd = self.software.run_abricate(database_name, file_out, SoftwareNames.SOFTWARE_ABRICATE_PARAMETERS, out_file)
		self.assertTrue(os.path.exists(out_file))

		parseOutFiles = ParseOutFiles()
		(vect_data, clean_abricate_file) = parseOutFiles.parse_abricate_file(out_file, 'temp_2.txt', 1)
		self.assertEqual('A', vect_data[0][ParseOutFiles.GENE])
		self.assertEqual(100.00, vect_data[0][ParseOutFiles.COVERAGE])
		self.assertEqual(99.69, vect_data[0][ParseOutFiles.IDENTITY])
		self.assertEqual('influenza_type', vect_data[0][ParseOutFiles.TYPE])
		self.assertEqual('XXXX', vect_data[0][ParseOutFiles.ACCESSION])
		self.assertEqual('NODE_1_length_3445_cov_491.652019', vect_data[0][ParseOutFiles.SEQ_NAME])
		
		self.assertEqual('N2', vect_data[2][ParseOutFiles.GENE])
		self.assertEqual(100.0, vect_data[2][ParseOutFiles.COVERAGE])
		self.assertEqual(98.65, vect_data[2][ParseOutFiles.IDENTITY])
		self.assertEqual('influenza_subtype', vect_data[2][ParseOutFiles.TYPE])
		self.assertEqual('XXXXXX', vect_data[2][ParseOutFiles.ACCESSION])
		self.assertEqual('NODE_7_length_1478_cov_488.215560', vect_data[2][ParseOutFiles.SEQ_NAME])

		## test abricate
		self.assertFalse(filecmp.cmp(out_file, clean_abricate_file))
		
		## remove file
		os.unlink(out_file)
		os.unlink(clean_abricate_file)
		
		## remove abricate db
		cmd = "rm -r %s/%s*" % (SoftwareNames.SOFTWARE_ABRICATE_DB, database_name)
		exist_status = os.system(cmd)
		self.assertTrue(exist_status == 0)
		self.utils.remove_dir(out_dir)

	def test_get_lines_and_average_reads(self):
		"""
		test lines and average reads
		"""
		software = Software()
		fastq = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		(lines, average) = software.get_lines_and_average_reads(fastq)
		self.assertEquals("44425", lines)
		self.assertEquals("143.6", average)


	def test_run_fastq(self):
		"""
		run fastq
		"""
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_1, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
		self.utils.copy_file(file_2, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2))
			
		sample_name = "run_fastq_and_trimmomatic1"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = False
			sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
			sample.path_name_2.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2)
			sample.owner = user
			sample.save()

		out_put_path = self.software.run_fastq(sample.get_fastq(TypePath.MEDIA_ROOT, True), sample.get_fastq(TypePath.MEDIA_ROOT, False))
		out_file_1 = os.path.join(out_put_path, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, True)))
		out_file_2 = os.path.join(out_put_path, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, False)))
		self.assertTrue(os.path.exists(out_file_1))
		self.assertTrue(os.path.exists(out_file_2))
		self.assertTrue(os.path.getsize(out_file_1) > 1000)
		self.assertTrue(os.path.getsize(out_file_2) > 1000)
		self.utils.remove_dir(out_put_path)

		out_put_path = self.software.run_fastq(file_1, None)
		out_file_1 = os.path.join(out_put_path, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, True)))
		out_file_2 = os.path.join(out_put_path, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, False)))
		self.assertTrue(os.path.exists(out_file_1))
		self.assertFalse(os.path.exists(out_file_2))
		self.assertTrue(os.path.getsize(out_file_1) > 1000)
		self.utils.remove_dir(out_put_path)
		self.utils.remove_dir(temp_dir)


	def test_run_trimmomatic(self):
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_1, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
		self.utils.copy_file(file_2, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2))
			
		sample_name = "run_fastq_and_trimmomatic"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = False
			sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
			sample.path_name_2.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2)
			sample.owner = user
			sample.save()

		(out_put_path_trimmomatic, parameters) = self.software.run_trimmomatic(sample.get_fastq(TypePath.MEDIA_ROOT, True), sample.get_fastq(TypePath.MEDIA_ROOT, False), sample_name)
		out_file_1 = os.path.join(out_put_path_trimmomatic, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)))
		out_file_2 = os.path.join(out_put_path_trimmomatic, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False)))
		self.assertTrue(os.path.exists(out_file_1))
		self.assertTrue(os.path.exists(out_file_2))
		self.assertTrue(os.path.getsize(out_file_1) > 1000)
		self.assertTrue(os.path.getsize(out_file_2) > 1000)
		
		out_put_path = self.software.run_fastq(out_file_1, out_file_2)
		out_file_1 = os.path.join(out_put_path, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)))
		out_file_2 = os.path.join(out_put_path, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)))
		self.assertTrue(os.path.exists(out_file_1))
		self.assertTrue(os.path.exists(out_file_2))
		self.assertTrue(os.path.getsize(out_file_1) > 1000)
		self.assertTrue(os.path.getsize(out_file_2) > 1000)
		self.utils.remove_dir(out_put_path)
		self.utils.remove_dir(out_put_path_trimmomatic)
		self.utils.remove_dir(temp_dir);

		sample_name = "sample_name"
		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_1, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
			
		sample_name = "run_fastq_and_trimmomatic_222"
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
			sample.owner = user
			sample.save()
		
		(out_put_path, parameters) = self.software.run_trimmomatic(sample.get_fastq(TypePath.MEDIA_ROOT, True), sample.get_fastq(TypePath.MEDIA_ROOT, False), sample_name)
		self.assertTrue(sample.get_fastq(TypePath.MEDIA_ROOT, False) == None)
		out_file_1 = os.path.join(out_put_path, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)))
		self.assertTrue(os.path.exists(out_file_1))
		self.assertTrue(os.path.getsize(out_file_1) > 1000)
		self.utils.remove_dir(out_put_path)

		## remove all files
		self.utils.remove_dir(temp_dir);


	def test_run_fastq_and_trimmomatic(self):
		"""
		Test run fastq and trimmomatic all together
 		"""
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_1, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
		self.utils.copy_file(file_2, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2))
			
		sample_name = "run_fastq_and_trimmomatic"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = True
			sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
			sample.path_name_2.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2)
			sample.owner = user
			sample.save()
		
		### set the job
		taskID = "xpto_task" 
		manageDatabase = ManageDatabase()
		manageDatabase.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue, taskID)

		### run software
		self.assertTrue(self.software.run_fastq_and_trimmomatic(sample, user))
		
		meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_sample != None)
		self.assertEquals(taskID, meta_sample.description)

		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)))))
		
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
		self.assertEqual('41254', result_average.number_file_1)
		self.assertEqual('142.0', result_average.average_file_1)
		self.assertEqual('41254', result_average.number_file_2)
		self.assertEqual('139.1', result_average.average_file_2)

		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Fastq(0.11.5), Trimmomatic(0.27)", list_meta[0].description)
		
		## remove all files
		cmd = "rm -r %s*" % (temp_dir); os.system(cmd)


	def test_run_fastq_and_trimmomatic_single_file(self):
		"""
		Test run fastq and trimmomatic all together
 		"""
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_1, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
			
		sample_name = "run_fastq_and_trimmomatic"
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
		self.assertTrue(self.software.run_fastq_and_trimmomatic(sample, user))
		
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)))))
		
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
		self.assertEqual('43560', result_average.number_file_1)
		self.assertEqual('141.0', result_average.average_file_1)
		self.assertEqual(None, result_average.number_file_2)
		self.assertEqual(None, result_average.average_file_2)

		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Fastq(0.11.5), Trimmomatic(0.27)", list_meta[0].description)
		
		## remove all files
		cmd = "rm -r %s*" % (temp_dir); os.system(cmd)

	def test_identify_type_and_sub_type(self):
		"""
 		get type and sub_type
 		"""
		
		uploadFiles = UploadFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		self.assertEqual("2", version)
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
		uploadFiles.upload_file(version, file)	## upload file
		
		try:
			uploadFile = UploadFile.objects.order_by('-version')[0]
			self.assertEqual("test_db_influenza_typing_v2", uploadFile.abricate_name)
		except UploadFile.DoesNotExist:
			self.fail("must have values")
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()
			
		sample_name = "identify_type_and_sub_type"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = True
			sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
			sample.path_name_2.name = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)
			sample.owner = user
			sample.save()
			
		return_value = self.software.identify_type_and_sub_type(sample, sample.path_name_1.name, sample.path_name_2.name, user, True)
		self.assertTrue(return_value)
		
		vect_identify_virus = sample.identify_virus.all()
		self.assertEqual(3, len(vect_identify_virus))
		for identify_virus in vect_identify_virus:
			if (identify_virus.rank == 0):
				self.assertEquals("100.00", identify_virus.coverage)
				self.assertEquals("99.69", identify_virus.identity)
				self.assertEquals("A", identify_virus.seq_virus.name)
				self.assertEquals(ConstantsVirus.SEQ_VIRUS_TYPE, identify_virus.seq_virus.kind_type.name)
			elif (identify_virus.rank == 1):
				self.assertEquals("100.00", identify_virus.coverage)
				self.assertEquals("99.18", identify_virus.identity)
				self.assertEquals("H3", identify_virus.seq_virus.name)
				self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, identify_virus.seq_virus.kind_type.name)
			elif (identify_virus.rank == 2):
				self.assertEquals("100.00", identify_virus.coverage)
				self.assertEquals("98.65", identify_virus.identity)
				self.assertEquals("N2", identify_virus.seq_virus.name)
				self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, identify_virus.seq_virus.kind_type.name)
		file_abricate = sample.get_abricate_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(sample.get_abricate_output(TypePath.MEDIA_ROOT)))
		
		file_spades_contigs = sample.get_draft_contigs_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(file_spades_contigs))
		file_abricate_contigs = sample.get_draft_contigs_abricate_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(file_abricate_contigs))
		
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Identify_Sample, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Identify_Sample, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Spades(3.11.1), Abricate(0.8-dev)", list_meta[0].description)
		if (os.path.exists(file_abricate)): os.unlink(file_abricate)
		if (os.path.exists(file_spades_contigs)): os.unlink(file_spades_contigs)
		if (os.path.exists(file_abricate_contigs)): os.unlink(file_abricate_contigs)

		contigs_2_sequences = Contigs2Sequences(True)
		## remove abricate db
		cmd = "rm -r %s/%s*" % (SoftwareNames.SOFTWARE_ABRICATE_DB, contigs_2_sequences.get_database_name())
		exist_status = os.system(cmd)
		self.assertTrue(exist_status == 0)
	
	def test_identify_type_and_sub_type_2(self):
		"""
 		get type and sub_type
 		"""
		
		uploadFiles = UploadFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		self.assertEqual("2", version)
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
		uploadFiles.upload_file(version, file)	## upload file
		
		try:
			uploadFile = UploadFile.objects.order_by('-version')[0]
			self.assertEqual("test_db_influenza_typing_v2", uploadFile.abricate_name)
		except UploadFile.DoesNotExist:
			self.fail("must have values")
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()
			
		sample_name = "identify_type_and_sub_type_10"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ10_1
			sample.path_name_1.name = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ10_1)
			sample.is_valid_2 = True
			sample.file_name_2 = ConstantsTestsCase.FASTQ10_2
			sample.path_name_2.name = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ10_2)
			sample.owner = user
			sample.save()
			
		return_value = self.software.identify_type_and_sub_type(sample, sample.path_name_1.name, sample.path_name_2.name, user, True)
		self.assertTrue(return_value)
		
		vect_identify_virus = sample.identify_virus.all()
		self.assertEqual(3, len(vect_identify_virus))
		for identify_virus in vect_identify_virus:
			if (identify_virus.rank == 0):
				self.assertEquals("100.00", identify_virus.coverage)
				self.assertEquals("99.69", identify_virus.identity)
				self.assertEquals("A", identify_virus.seq_virus.name)
				self.assertEquals(ConstantsVirus.SEQ_VIRUS_TYPE, identify_virus.seq_virus.kind_type.name)
			elif (identify_virus.rank == 1):
				self.assertEquals("100.00", identify_virus.coverage)
				self.assertEquals("99.18", identify_virus.identity)
				self.assertEquals("H3", identify_virus.seq_virus.name)
				self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, identify_virus.seq_virus.kind_type.name)
			elif (identify_virus.rank == 2):
				self.assertEquals("100.00", identify_virus.coverage)
				self.assertEquals("98.37", identify_virus.identity)
				self.assertEquals("N2", identify_virus.seq_virus.name)
				self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, identify_virus.seq_virus.kind_type.name)
		file_abricate = sample.get_abricate_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(sample.get_abricate_output(TypePath.MEDIA_ROOT)))
		
		file_spades_contigs = sample.get_draft_contigs_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(file_spades_contigs))
		file_abricate_contigs = sample.get_draft_contigs_abricate_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(file_abricate_contigs))
		
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Identify_Sample, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Identify_Sample, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Spades(3.11.1), Abricate(0.8-dev)", list_meta[0].description)
		if (os.path.exists(file_abricate)): os.unlink(file_abricate)
		if (os.path.exists(file_spades_contigs)): os.unlink(file_spades_contigs)
		if (os.path.exists(file_abricate_contigs)): os.unlink(file_abricate_contigs)

		contigs_2_sequences = Contigs2Sequences(True)
		## remove abricate db
		cmd = "rm -r %s/%s*" % (SoftwareNames.SOFTWARE_ABRICATE_DB, contigs_2_sequences.get_database_name())
		exist_status = os.system(cmd)
		self.assertTrue(exist_status == 0)


	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_run_fastq_and_trimmomatic_and_identify_species(self):
		"""
		Test run fastq and trimmomatic all together
		"""

		uploadFiles = UploadFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		self.assertEqual("2", version)
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
		uploadFiles.upload_file(version, file)	## upload file
		
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_1, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
		self.utils.copy_file(file_2, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2))
			
		sample_name = "run_fastq_and_trimmomatic"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = True
			sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
			sample.path_name_2.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2)
			sample.owner = user
			sample.save()
		
		self.assertFalse(sample.is_ready_for_projects)
		
		### set the job
		taskID = "xpto_task" 
		manageDatabase = ManageDatabase()
		manageDatabase.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue, taskID)

		### run software
		self.assertTrue(self.software.run_fastq_and_trimmomatic_and_identify_species(sample, user))
		
		meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_sample != None)
		self.assertEquals(taskID, meta_sample.description)

		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)))))
		
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
		self.assertEqual('41254', result_average.number_file_1)
		self.assertEqual('142.0', result_average.average_file_1)
		self.assertEqual('41254', result_average.number_file_2)
		self.assertEqual('139.1', result_average.average_file_2)

		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Fastq(0.11.5), Trimmomatic(0.27)", list_meta[0].description)
		
		sample = Sample.objects.get(pk=sample.id)
		self.assertTrue(sample.is_ready_for_projects)
		self.assertTrue(sample.get_is_ready_for_projects())
		
		file_abricate = sample.get_abricate_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(sample.get_abricate_output(TypePath.MEDIA_ROOT)))
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Identify_Sample, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Identify_Sample, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Spades(3.11.1), Abricate(0.8-dev)", list_meta[0].description)
		self.assertEquals("A-H3N2", sample.type_subtype)
		self.assertEquals("A-H3N2", sample.get_type_sub_type())
		if (os.path.exists(file_abricate)): os.unlink(file_abricate)
		
		### mixed infections
		self.assertEquals(0, sample.number_alerts)
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, sample.mixed_infections_tag.name)
		
		manage_database = ManageDatabase()
		meta_data = manage_database.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE,\
								MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_data == None)
		
		meta_data = manage_database.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_TAG_MIXED_INFECTION_TYPE_SUBTYPE,\
								MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_data != None)
		self.assertTrue(meta_data.description == ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO)
		
		## remove all files
		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))


	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_run_fastq_and_trimmomatic_and_identify_species_2(self):
		"""
		Test run fastq and trimmomatic all together
		"""

		uploadFiles = UploadFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		self.assertEqual("2", version)
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
		uploadFiles.upload_file(version, file)	## upload file
		
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ5_1)
		file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ5_2)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_1, os.path.join(temp_dir, ConstantsTestsCase.FASTQ5_1))
		self.utils.copy_file(file_2, os.path.join(temp_dir, ConstantsTestsCase.FASTQ5_2))
			
		sample_name = "run_fastq_and_trimmomatic"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ5_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ5_1)
			sample.is_valid_2 = True
			sample.file_name_2 = ConstantsTestsCase.FASTQ5_2
			sample.path_name_2.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ5_2)
			sample.owner = user
			sample.save()
		
		self.assertFalse(sample.is_ready_for_projects)
		
		### set the job
		taskID = "xpto_task" 
		manageDatabase = ManageDatabase()
		manageDatabase.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue, taskID)

		### run software
		self.assertTrue(self.software.run_fastq_and_trimmomatic_and_identify_species(sample, user))
		
		meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_sample != None)
		self.assertEquals(taskID, meta_sample.description)

		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)))))
		
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
		self.assertEqual('128050', result_average.number_file_1)
		self.assertEqual('139.5', result_average.average_file_1)
		self.assertEqual('128050', result_average.number_file_2)
		self.assertEqual('136.8', result_average.average_file_2)

		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Fastq(0.11.5), Trimmomatic(0.27)", list_meta[0].description)
		
		sample = Sample.objects.get(pk=sample.id)
		self.assertTrue(sample.is_ready_for_projects)
		self.assertTrue(sample.get_is_ready_for_projects())
		
		file_abricate = sample.get_abricate_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(sample.get_abricate_output(TypePath.MEDIA_ROOT)))
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Identify_Sample, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Identify_Sample, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Spades(3.11.1), Abricate(0.8-dev)", list_meta[0].description)
		self.assertEquals("A-H1|H3|N1|N2", sample.type_subtype)
		self.assertEquals("A-H1|H3|N1|N2", sample.get_type_sub_type())
		if (os.path.exists(file_abricate)): os.unlink(file_abricate)
		
		### mixed infections
		self.assertEquals(1, sample.number_alerts)
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, sample.mixed_infections_tag.name)
		
		manage_database = ManageDatabase()
		meta_data = manage_database.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE,\
								MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_data != None)
		self.assertEquals("Warning: more than two subtypes were detected for this sample, suggesting that may represent a 'mixed infection'.", meta_data.description)
		
		meta_data = manage_database.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_TAG_MIXED_INFECTION_TYPE_SUBTYPE,\
								MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_data != None)
		self.assertTrue(meta_data.description == ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES)
		
		## remove all files
		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))
	
	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_run_fastq_and_trimmomatic_and_identify_species_3(self):
		"""
 		Test run fastq and trimmomatic all together
 		"""

		uploadFiles = UploadFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		self.assertEqual("2", version)
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
		uploadFiles.upload_file(version, file)	## upload file
		
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ6_1)
		file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ6_2)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_1, os.path.join(temp_dir, ConstantsTestsCase.FASTQ6_1))
		self.utils.copy_file(file_2, os.path.join(temp_dir, ConstantsTestsCase.FASTQ6_2))
			
		sample_name = "run_fastq_and_trimmomatic_new"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ6_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ6_1)
			sample.is_valid_2 = True
			sample.file_name_2 = ConstantsTestsCase.FASTQ6_2
			sample.path_name_2.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ6_2)
			sample.owner = user
			sample.save()
		
		self.assertFalse(sample.is_ready_for_projects)
		
		### set the job
		taskID = "xpto_task" 
		manageDatabase = ManageDatabase()
		manageDatabase.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue, taskID)

		### run software
		self.assertTrue(self.software.run_fastq_and_trimmomatic_and_identify_species(sample, user))
		
		meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_sample != None)
		self.assertEquals(taskID, meta_sample.description)

		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)))))
		
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
		self.assertEqual('152623', result_average.number_file_1)
		self.assertEqual('144.4', result_average.average_file_1)
		self.assertEqual('152623', result_average.number_file_2)
		self.assertEqual('142.3', result_average.average_file_2)

		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Fastq(0.11.5), Trimmomatic(0.27)", list_meta[0].description)
		
		sample = Sample.objects.get(pk=sample.id)
		self.assertTrue(sample.is_ready_for_projects)
		self.assertTrue(sample.get_is_ready_for_projects())
		
		file_abricate = sample.get_abricate_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(sample.get_abricate_output(TypePath.MEDIA_ROOT)))
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Identify_Sample, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Identify_Sample, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Spades(3.11.1), Abricate(0.8-dev)", list_meta[0].description)
		self.assertEquals("A-H3N2", sample.type_subtype)
		self.assertEquals("A-H3N2", sample.get_type_sub_type())
		if (os.path.exists(file_abricate)): os.unlink(file_abricate)
		
		### mixed infections
		self.assertEquals(0, sample.number_alerts)
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, sample.mixed_infections_tag.name)
		
		manage_database = ManageDatabase()
		meta_data = manage_database.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE,\
 								MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_data == None)
		
		meta_data = manage_database.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_TAG_MIXED_INFECTION_TYPE_SUBTYPE,\
								MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_data != None)
		self.assertTrue(meta_data.description == ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO)
		
		## remove all files
		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))


	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_run_fastq_and_trimmomatic_and_identify_species_4(self):
		"""
 		Test run fastq and trimmomatic all together
 		"""

		uploadFiles = UploadFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		self.assertEqual("2", version)
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
		uploadFiles.upload_file(version, file)	## upload file
		
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ8_1)
		file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ8_2)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_1, os.path.join(temp_dir, ConstantsTestsCase.FASTQ8_1))
		self.utils.copy_file(file_2, os.path.join(temp_dir, ConstantsTestsCase.FASTQ8_2))
			
		sample_name = "run_fastq_and_trimmomatic_new"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ8_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ8_1)
			sample.is_valid_2 = True
			sample.file_name_2 = ConstantsTestsCase.FASTQ8_2
			sample.path_name_2.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ8_2)
			sample.owner = user
			sample.save()
		
		self.assertFalse(sample.is_ready_for_projects)
		
		### set the job
		taskID = "xpto_task" 
		manageDatabase = ManageDatabase()
		manageDatabase.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue, taskID)

		### run software
		self.assertTrue(self.software.run_fastq_and_trimmomatic_and_identify_species(sample, user))
		
		meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_sample != None)
		self.assertEquals(taskID, meta_sample.description)

		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)))))
		
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
		self.assertEqual('83012', result_average.number_file_1)
		self.assertEqual('147.1', result_average.average_file_1)
		self.assertEqual('83012', result_average.number_file_2)
		self.assertEqual('145.5', result_average.average_file_2)

		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Fastq(0.11.5), Trimmomatic(0.27)", list_meta[0].description)
		
		sample = Sample.objects.get(pk=sample.id)
		self.assertTrue(sample.is_ready_for_projects)
		self.assertTrue(sample.get_is_ready_for_projects())
		
		file_abricate = sample.get_abricate_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(sample.get_abricate_output(TypePath.MEDIA_ROOT)))
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Identify_Sample, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Identify_Sample, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Spades(3.11.1), Abricate(0.8-dev)", list_meta[0].description)
		self.assertEquals("Yamagata", sample.type_subtype)
		self.assertEquals("Yamagata", sample.get_type_sub_type())
		
		vect_data = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, vect_data[0])
		self.assertEquals(1, vect_data[1])
		self.assertEquals("Warning: an incomplete type/subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).", vect_data[2])
		if (os.path.exists(file_abricate)): os.unlink(file_abricate)
		
		### mixed infections
		self.assertEquals(1, sample.number_alerts)
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, sample.mixed_infections_tag.name)
		
		manage_database = ManageDatabase()
		meta_data = manage_database.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE,\
 								MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_data != None)
		self.assertEquals("Warning: an incomplete type/subtype has been assigned (possible reasons: low number of influenza reads, " +\
						"same-subtype mixed infection, etc.).", meta_data.description)

		meta_data = manage_database.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_TAG_MIXED_INFECTION_TYPE_SUBTYPE,\
								MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_data != None)
		self.assertTrue(meta_data.description == ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO)
		
		## remove all files
		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))


	def test_run_snippy(self):
		"""
		run snippy
		"""

		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_1, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
		self.utils.copy_file(file_2, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2))
			
		sample_name = "run_snippy1_sdfs"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = False
			sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
			sample.path_name_2.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2)
			sample.owner = user
			sample.save()

		project_name = "file_name_3_dsf"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = sample_name
			project.owner = user
			project.save()
		
		project_sample = ProjectSample()
		project_sample.sample = sample
		project_sample.project = project
		project_sample.save()
		
		out_put_path = self.software.run_snippy(sample.get_fastq(TypePath.MEDIA_ROOT, True), sample.get_fastq(TypePath.MEDIA_ROOT, False), gb_file, sample.name)
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_BAM, "")))))
		file_size = os.path.getsize(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_BAM, ""))))
		self.assertTrue(file_size > 1000 )
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_CONSENSUS_FA, "")))))
		self.assertFalse(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_CONSENSUS_FASTA, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_CSV, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_DEPTH_GZ, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_DEPTH_GZ_TBI, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_TAB, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_VCF, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_VCF_GZ, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, "reference", os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_REF_FASTA, "")))))
		remove_path = os.path.dirname(out_put_path)
		if (len(remove_path.split('/')) > 2): self.utils.remove_dir(remove_path)
		else: self.utils.remove_dir(out_put_path)

		out_put_path = self.software.run_snippy(sample.get_fastq(TypePath.MEDIA_ROOT, True), None, gb_file, sample.name)
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_BAM, "")))))
		file_size_2 = os.path.getsize(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_BAM, ""))))
		self.assertTrue(file_size_2 > 1000)
		self.assertTrue(file_size_2 < file_size)
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_CONSENSUS_FA, "")))))
		self.assertFalse(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_CONSENSUS_FASTA, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_CSV, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_DEPTH_GZ, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_DEPTH_GZ_TBI, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_TAB, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_VCF, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_VCF_GZ, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, "reference", os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_REF_FASTA, "")))))
		remove_path = os.path.dirname(out_put_path)
		if (len(remove_path.split('/')) > 2): self.utils.remove_dir(remove_path)
		else: self.utils.remove_dir(out_put_path)
		self.utils.remove_dir(temp_dir)

	def test_run_freebayes(self):
		"""
 		test run freebayes
 		create a VCF
 		"""
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		genbank_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		bam_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_BAM, "run_snippy1.bam")	
		expect_tab_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "run_snippy1.tab")	
		vcf_expect_result = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "run_snippy1.vcf")	

		sample_name = "run_snippy1"
		out_put_path = self.software.run_freebayes(bam_file, fasta_file, genbank_file, sample_name)
		result_file = os.path.join(out_put_path, sample_name + ".vcf")
		result_file_tab = os.path.join(out_put_path, sample_name + ".tab")
		self.assertTrue(os.path.exists(result_file))
		self.assertTrue(os.path.getsize(result_file) > 1000)
		self.assertTrue(os.path.exists(result_file_tab))
		self.assertTrue(os.path.getsize(result_file_tab) > 1000)
		temp_file = self.utils.get_temp_file("file_name", ".txt")
		temp_file_1 = self.utils.get_temp_file("file_name", ".txt")
		cmd = "grep -E -v '##command|##reference|##fileDate|insaFlu' {} > {}".format(result_file, temp_file)
		os.system(cmd);
		cmd = "grep -E -v '##command|##reference|##fileDate|insaFlu' {} > {}".format(vcf_expect_result, temp_file_1)
		os.system(cmd);
		self.assertTrue(filecmp.cmp(temp_file_1, temp_file))
		self.assertTrue(filecmp.cmp(expect_tab_file, result_file_tab))
		self.utils.remove_temp_file(temp_file_1)
		self.utils.remove_temp_file(temp_file)
		self.utils.remove_dir(out_put_path)
		
	def test_run_freebayes_parallel(self):
		"""
 		test run freebayes
 		create a VCF
 		"""
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		genbank_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		bam_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_BAM, "run_snippy1.bam")	
		expect_tab_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "run_freebayes_parallel.tab")
		vcf_expect_result = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "run_freebayes_parallel.vcf")

		sample_name = "run_snippy1"
		out_put_path = self.software.run_freebayes_parallel(bam_file, fasta_file, genbank_file, sample_name)
		result_file = os.path.join(out_put_path, sample_name + ".vcf")
		result_file_tab = os.path.join(out_put_path, sample_name + ".tab")
		self.assertTrue(os.path.exists(result_file))
		self.assertTrue(os.path.getsize(result_file) > 1000)
		self.assertTrue(os.path.exists(result_file_tab))
		self.assertTrue(os.path.getsize(result_file_tab) > 1000)
		temp_file = self.utils.get_temp_file("file_name", ".txt")
		temp_file_1 = self.utils.get_temp_file("file_name", ".txt")
		cmd = "grep -E -v '##command|##reference|##fileDate|insaFlu' {} > {}".format(result_file, temp_file)
		os.system(cmd);
		cmd = "grep -E -v '##command|##reference|##fileDate|insaFlu' {} > {}".format(vcf_expect_result, temp_file_1)
		os.system(cmd);
		self.assertTrue(filecmp.cmp(temp_file_1, temp_file))
		self.assertTrue(filecmp.cmp(expect_tab_file, result_file_tab))
		self.utils.remove_temp_file(temp_file_1)
		self.utils.remove_temp_file(temp_file)
		self.utils.remove_dir(out_put_path)


	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_process_second_stage_analysis_single_file(self):
		"""
 		test global method
 		"""
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT_TEST", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))

		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		ref_name = "second_stagis_single_file2"
		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.reference_fasta.name = fasta_file
			reference.reference_fasta_name = os.path.basename(fasta_file)
			reference.reference_genbank.name = gb_file
			reference.reference_genbank_name = os.path.basename(gb_file)
			reference.owner = user
			reference.save()
			
		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_1, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
		self.utils.copy_file(file_2, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2))
			
		sample_name = "run_snippyis_single_file1"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = False
			sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
			sample.path_name_2.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2)
			sample.owner = user
			sample.save()

		project_name = "file_naais_single_filee_3"
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
		
		### copy to trimmomatic files
		self.utils.copy_file(file_1, project_sample.sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True))
		self.utils.copy_file(file_2, project_sample.sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False))
		
		manageDatabase = ManageDatabase()
		metaKeyAndValue = MetaKeyAndValue()
		meta_key_project_sample = metaKeyAndValue.get_meta_key_queue_by_project_sample_id(project_sample.id)
		manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue, "meta_sample.description")
		self.assertTrue(self.software.process_second_stage_snippy_coverage_freebayes(project_sample, user))
		try:
			project_sample = ProjectSample.objects.get(pk=project_sample.id)
		except ProjectSample.DoesNotExist:
			self.fail("Must exist")
		self.assertTrue(project_sample.is_finished)
		
		### test the files
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertFalse(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FA, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, SoftwareNames.SOFTWARE_SNIPPY_name).index(SoftwareNames.SOFTWARE_SNIPPY_name.lower()) != -1)
		
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_BAM, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_BAM_BAI, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertFalse(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ_TBI, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CSV, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_REF_FASTA, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_REF_FASTA, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_REF_FASTA_FAI, SoftwareNames.SOFTWARE_SNIPPY_name)))
		
		## freebayes
		self.assertTrue(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_FREEBAYES_name).index(SoftwareNames.SOFTWARE_FREEBAYES_name.lower()) != -1)
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, SoftwareNames.SOFTWARE_FREEBAYES_name)))

		## test freebayes clean
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		self.assertTrue(os.path.getsize(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)) > 10)

		## test consensus file
		self.assertTrue(os.path.exists(project_sample.get_consensus_file(TypePath.MEDIA_ROOT)))
		self.assertTrue(os.path.getsize(project_sample.get_consensus_file(TypePath.MEDIA_ROOT)) > 10)

		### human file name, snippy tab
		self.assertTrue(os.path.exists(project_sample.get_file_output_human(TypePath.MEDIA_ROOT, FileType.FILE_TAB,  self.software_names.get_snippy_name())))
		self.assertTrue(os.path.getsize(project_sample.get_file_output_human(TypePath.MEDIA_ROOT, FileType.FILE_TAB,  self.software_names.get_snippy_name())) > 10)
		
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
		self.assertEquals(125, count_hits.get_total())
		
		self.assertEquals(3, project_sample.count_variations.var_less_50)
		self.assertEquals(0, project_sample.count_variations.var_bigger_50_90)
		self.assertEquals(122, project_sample.count_variations.var_bigger_90)
		self.assertEquals(0, project_sample.alert_first_level)
		self.assertEquals(0, project_sample.alert_second_level)
		
		### test mixed infections
		self.assertEquals('0.815100969586423', '{}'.format(project_sample.mixed_infections.average_value))
		self.assertEquals('No', project_sample.mixed_infections.tag.name)
		self.assertFalse(project_sample.mixed_infections.has_master_vector)
		
		list_meta = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Snippy_Freebayes, None)
		self.assertEquals(1, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Snippy_Freebayes, list_meta[0].meta_tag.name)
		
		decode_result = DecodeObjects()
		result = decode_result.decode_result(list_meta[0].description)
		self.assertTrue(result is not None)
		self.assertEquals("Snippy-3.2-dev; (--mapqual 20 --mincov 10 --minfrac 0.51)", result.get_software(self.software_names.get_snippy_name()))
		self.assertEquals("Freebayes-v1.1.0-54-g49413aa; (-p 2 -q 20 -m 20 --min-coverage 100 --min-alternate-fraction 0.01 --min-alternate-count 10 -V)",\
 						result.get_software(self.software_names.get_freebayes_name()))

		meta_value = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
		decode_coverage = DecodeObjects()
		coverage = decode_coverage.decode_result(meta_value.description)
		self.assertEqual(coverage.get_coverage('PA', Coverage.COVERAGE_ALL), "527.4")
		self.assertEqual(coverage.get_coverage('MP', Coverage.COVERAGE_ALL), "2198.8")
		self.assertEqual(coverage.get_coverage('HA', Coverage.COVERAGE_ALL), "1449.3")
		self.assertEqual(coverage.get_coverage('PA', Coverage.COVERAGE_MORE_9), "100.0")
		self.assertEqual(coverage.get_coverage('MP', Coverage.COVERAGE_MORE_9), "100.0")
		self.assertEqual(coverage.get_coverage('HA', Coverage.COVERAGE_MORE_9), "100.0")
		self.assertEqual(coverage.get_coverage('PA', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('MP', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('HA', Coverage.COVERAGE_MORE_0), "100.0")
		
		lst_meta_sample = manageDatabase.get_project_sample_metakey(project_sample, meta_key_project_sample, None)
		self.assertEquals(2, len(lst_meta_sample))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Queue, lst_meta_sample[1].value)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, lst_meta_sample[0].value)
		
		### check if the images exist
		for gene in self.utils.get_elements_from_db(reference, user):
			output_image = project_sample.get_global_file_by_element(TypePath.MEDIA_ROOT, ProjectSample.PREFIX_FILE_COVERAGE, gene, FileExtensions.FILE_PNG)
			self.assertTrue(os.path.exists(output_image))

		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))


	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_process_second_stage_analysis_single_file_2(self):
		"""
 		test global method
 		"""
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT_TEST", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))

		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ2_1)
		file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ2_2)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		ref_name = "second_stagis_single_file2"
		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.reference_fasta.name = fasta_file
			reference.reference_fasta_name = os.path.basename(fasta_file)
			reference.reference_genbank.name = gb_file
			reference.reference_genbank_name = os.path.basename(gb_file)
			reference.owner = user
			reference.save()
			
		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_1, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
		self.utils.copy_file(file_2, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2))
			
		sample_name = "run_snippyis_single_file1"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = False
			sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
			sample.path_name_2.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2)
			sample.owner = user
			sample.save()

		project_name = "file_naais_single_filee_3"
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
		
		### copy to trimmomatic files
		self.utils.copy_file(file_1, project_sample.sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True))
		self.utils.copy_file(file_2, project_sample.sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False))

		manageDatabase = ManageDatabase()
		metaKeyAndValue = MetaKeyAndValue()
		meta_key_project_sample = metaKeyAndValue.get_meta_key_queue_by_project_sample_id(project_sample.id)
		manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue, "meta_sample.description")
		self.assertTrue(self.software.process_second_stage_snippy_coverage_freebayes(project_sample, user))
		try:
			project_sample = ProjectSample.objects.get(pk=project_sample.id)
		except ProjectSample.DoesNotExist:
			self.fail("Must exist")
		self.assertTrue(project_sample.is_finished)
		
		### test the files
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertFalse(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FA, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, SoftwareNames.SOFTWARE_SNIPPY_name).index(SoftwareNames.SOFTWARE_SNIPPY_name.lower()) != -1)
		
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_BAM, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_BAM_BAI, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertFalse(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ_TBI, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CSV, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_SNIPPY_name)))
		## freebayes
		self.assertTrue(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_FREEBAYES_name).index(SoftwareNames.SOFTWARE_FREEBAYES_name.lower()) != -1)
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, SoftwareNames.SOFTWARE_FREEBAYES_name)))


		### set flag that is finished
		list_meta = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Count_Hits, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Count_Hits, list_meta[0].meta_tag.name)
		
		### get the hits value
		decode_coverage = DecodeObjects()
		count_hits = decode_coverage.decode_result(list_meta[0].description)
		self.assertEquals(0, count_hits.get_hits_50_90())
		self.assertEquals(0, count_hits.get_hits_less_50())
		self.assertEquals(0, count_hits.get_total_50_50_90())
		self.assertEquals(4, count_hits.get_total())
		
		self.assertEquals(0, project_sample.count_variations.var_less_50)
		self.assertEquals(0, project_sample.count_variations.var_bigger_50_90)
		self.assertEquals(4, project_sample.count_variations.var_bigger_90)
		self.assertEquals(0, project_sample.alert_first_level)
		self.assertEquals(7, project_sample.alert_second_level)
		
		### test mixed infections
		self.assertEquals('0.0', '{}'.format(project_sample.mixed_infections.average_value))
		self.assertEquals('No', project_sample.mixed_infections.tag.name)
		self.assertFalse(project_sample.mixed_infections.has_master_vector)
		
		list_meta = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Snippy_Freebayes, None)
		self.assertEquals(1, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Snippy_Freebayes, list_meta[0].meta_tag.name)
		
		decode_result = DecodeObjects()
		result = decode_result.decode_result(list_meta[0].description)
		self.assertTrue(result is not None)
		self.assertEquals("Snippy-3.2-dev; (--mapqual 20 --mincov 10 --minfrac 0.51)", result.get_software(self.software_names.get_snippy_name()))
		self.assertEquals("Freebayes-v1.1.0-54-g49413aa; (-p 2 -q 20 -m 20 --min-coverage 100 --min-alternate-fraction 0.01 --min-alternate-count 10 -V)",\
 						result.get_software(self.software_names.get_freebayes_name()))

		meta_value = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
		decode_coverage = DecodeObjects()
		coverage = decode_coverage.decode_result(meta_value.description)

		self.assertEqual(coverage.get_coverage('PA', Coverage.COVERAGE_ALL), "11.0")
		self.assertEqual(coverage.get_coverage('MP', Coverage.COVERAGE_ALL), "139.5")
		self.assertEqual(coverage.get_coverage('HA', Coverage.COVERAGE_ALL), "20.2")
		self.assertEqual(coverage.get_coverage('PA', Coverage.COVERAGE_MORE_9), "44.1")
		self.assertEqual(coverage.get_coverage('MP', Coverage.COVERAGE_MORE_9), "100.0")
		self.assertEqual(coverage.get_coverage('HA', Coverage.COVERAGE_MORE_9), "85.4")
		self.assertEqual(coverage.get_coverage('PA', Coverage.COVERAGE_MORE_0), "99.3")
		self.assertEqual(coverage.get_coverage('MP', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('HA', Coverage.COVERAGE_MORE_0), "100.0")
		
		lst_meta_sample = manageDatabase.get_project_sample_metakey(project_sample, meta_key_project_sample, None)
		self.assertEquals(2, len(lst_meta_sample))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Queue, lst_meta_sample[1].value)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, lst_meta_sample[0].value)
		
		### check if the images exist
		for gene in self.utils.get_elements_from_db(reference, user):
			output_image = project_sample.get_global_file_by_element(TypePath.MEDIA_ROOT, ProjectSample.PREFIX_FILE_COVERAGE, gene, FileExtensions.FILE_PNG)
			self.assertTrue(os.path.exists(output_image))

		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))


	def test_run_snippy_vcf_to_tab(self):
		"""
  		test snippy_vcf_to_tab method
  		"""
		
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		vcf_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "temp_more_REF.vcf")
		result_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "resutl_vcf_to_tab.tab")

		out_file = self.utils.get_temp_file("snippy_vcf_to_tab", ".tab")
		out_file_2 = self.software.run_snippy_vcf_to_tab(fasta_file, gb_file, vcf_file, out_file)
		self.assertEquals(out_file, out_file_2)
		
		self.assertTrue(filecmp.cmp(out_file_2, result_file))
		os.unlink(out_file)


	def test_test_bgzip_and_tbi_in_vcf(self):
		"""
		test snippy_vcf_to_tab method
		"""
	
		vcf_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, 'temp.vcf')
		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(vcf_file, os.path.join(temp_dir, os.path.basename(vcf_file)))
		
		new_vcf_file = os.path.join(temp_dir, os.path.basename(vcf_file)) 
		self.software.test_bgzip_and_tbi_in_vcf(new_vcf_file)
		
		self.software.test_bgzip_and_tbi_in_vcf(new_vcf_file)
		
		self.assertTrue(os.path.exists(new_vcf_file))
		self.assertTrue(os.path.exists(new_vcf_file + '.gz'))
		self.assertTrue(os.path.getsize(new_vcf_file + '.gz') > 1000)
		self.assertTrue(os.path.getsize(new_vcf_file + '.gz.tbi') > 200)
		self.assertTrue(os.path.exists(new_vcf_file + '.gz.tbi'))
		self.utils.remove_dir(temp_dir)


	def test_run_genbank2gff3(self):
		"""
		test genbank2gff3 method
		"""
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		gff_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GFF)
		self.assertTrue(os.path.exists(gb_file))
		out_file = self.utils.get_temp_file("file_name", ".txt")
		out_file_2 = self.software.run_genbank2gff3(gb_file, out_file)
		self.assertEquals(out_file, out_file_2)
		
		self.assertTrue(filecmp.cmp(out_file_2, gff_file))
		os.unlink(out_file)

	def test_run_genbank2gff3_2(self):
		"""
		test genbank2gff3 method
		"""
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_RIM_sample_gbk)
		gff_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_RIM_sample_gff)
		self.assertTrue(os.path.exists(gb_file))
		out_file = self.utils.get_temp_file("file_name_rim", ".txt")
		out_file_2 = self.software.run_genbank2gff3(gb_file, out_file)
		self.assertEquals(out_file, out_file_2)
		
		self.assertTrue(filecmp.cmp(out_file_2, gff_file))
		os.unlink(out_file)

	def test_run_get_snpeff_config(self):
		"""
		test get_snpeff_config method
		"""
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		snpeff_config_expect = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_SNPEF_config)
		
		(base_name, out_file) = self.software.get_snpeff_config(fasta_file)
		file_name = os.path.basename(fasta_file)
		file_name = file_name[0:file_name.rfind('.')]
		self.assertEquals(base_name, file_name)
		self.assertTrue(filecmp.cmp(out_file, snpeff_config_expect))
		os.unlink(out_file)

	def test_run_snpEff(self):
		"""
		test snpEff method
		"""
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		genbank_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		freebayes_vcf = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, ConstantsTestsCase.MANAGING_FILES_FREEBAYES_VCF)
		freebayes_expect_vcf = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, ConstantsTestsCase.MANAGING_FILES_FREEBAYES_ANNOTATED_VCF)
		freebayes_expect_vcf_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, ConstantsTestsCase.MANAGING_FILES_FREEBAYES_ANNOTATED_VCF_2)
		
		out_file = self.utils.get_temp_file("file_name", ".vcf")
		out_file_2 = self.software.run_snpEff(fasta_file, genbank_file, freebayes_vcf, out_file)
		self.assertEquals(out_file, out_file_2)
		
		out_file_clean = self.utils.get_temp_file("file_name", ".vcf")
		cmd = "grep -v '{}' {} > {}".format(os.path.dirname(out_file), out_file, out_file_clean)
		os.system(cmd)
		self.assertTrue(filecmp.cmp(out_file_clean, freebayes_expect_vcf) or filecmp.cmp(out_file_clean, freebayes_expect_vcf_2))
		os.unlink(out_file_clean)
		os.unlink(out_file)

	def test_run_snpEff_2(self):
		"""
		test snpEff method
		"""
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		genbank_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		freebayes_vcf = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "run_snpeff.vcf")
		freebayes_expect_vcf = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "run_snpeff_expected.vcf")
		
		out_file = self.utils.get_temp_file("file_name", ".vcf")
		out_file_2 = self.software.run_snpEff(fasta_file, genbank_file, freebayes_vcf, out_file)
		self.assertEquals(out_file, out_file_2)
		
		out_file_clean = self.utils.get_temp_file("file_name", ".vcf")
		cmd = "grep -v '{}' {} > {}".format(os.path.dirname(out_file), out_file, out_file_clean)
		os.system(cmd)
		self.assertTrue(filecmp.cmp(out_file_clean, freebayes_expect_vcf))
		os.unlink(out_file_clean)
		os.unlink(out_file)
	
	def test_run_snpEff_3(self):
		"""
		test snpEff method
		"""
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_PV1_KX162693_fasta)
		genbank_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_PV1_KX162693_gbk)
		freebayes_vcf = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "run_snpeff_pv1_kx162693.vcf")
		freebayes_expect_vcf = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "run_snpeff_pv1_kx162693_expected.vcf")
		
		out_file = self.utils.get_temp_file("file_name", ".vcf")
		out_file_2 = self.software.run_snpEff(fasta_file, genbank_file, freebayes_vcf, out_file)
		self.assertEquals(out_file, out_file_2)
		
		out_file_clean = self.utils.get_temp_file("file_name", ".vcf")
		cmd = "grep -v '{}' {} > {}".format(os.path.dirname(out_file), out_file, out_file_clean)
		os.system(cmd)
		self.assertTrue(filecmp.cmp(out_file_clean, freebayes_expect_vcf))
		os.unlink(out_file_clean)
		os.unlink(out_file)
		
	def test_run_prokka(self):
		"""
		test prokka method
		"""
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, 'TwoGenesJoined.fasta')
		temp_out = self.software.run_prokka(fasta_file, 'TwoGenesJoined.fasta')

		name_strain = os.path.basename(fasta_file)
		name_strain = name_strain[:name_strain.rfind('.')]
		genbank_file = os.path.join(temp_out, name_strain + FileExtensions.FILE_GBK)
		self.assertTrue(os.path.exists(genbank_file))
		
		geneticElement = self.utils.get_elements_and_genes(genbank_file)
		self.assertEquals(2, len(geneticElement.get_sorted_elements()))
		self.assertTrue('PB2' in geneticElement.get_sorted_elements())
		self.assertTrue('MP' in geneticElement.get_sorted_elements())
		self.assertFalse('xptoNS' in geneticElement.get_sorted_elements())
		
		self.assertEquals(30, geneticElement.get_genes('PB2')[0].start)
		self.assertEquals(2280, geneticElement.get_genes('PB2')[0].end)
		self.assertEquals('PB2', geneticElement.get_genes('PB2')[0].name)
		self.assertTrue(geneticElement.get_genes('PB2')[0].is_forward())
		self.assertEquals(2325, geneticElement.get_genes('PB2')[1].start)
		self.assertEquals(4599, geneticElement.get_genes('PB2')[1].end)
		self.assertEquals('PB1', geneticElement.get_genes('PB2')[1].name)
		self.assertFalse(geneticElement.get_genes('PB2')[1].is_forward())
		
		self.assertEquals(0, geneticElement.get_genes('MP')[0].start)
		self.assertEquals(759, geneticElement.get_genes('MP')[0].end)
		self.assertEquals('M', geneticElement.get_genes('MP')[0].name)
		self.assertTrue(geneticElement.get_genes('MP')[0].is_forward())
		
		try:
			self.utils.compare_locus_fasta_gb(fasta_file, genbank_file)
		except ValueError as e:
			self.fail(e.args[0])
		self.utils.remove_dir(temp_out)

	def test_run_mauve(self):
		"""
 		test run mauve
 		create a VCF
 		"""
		## get all fastq files
		vect_files = self.constants_tests_case.get_all_consensus_files(self.baseDirectory)
		expect_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.FILE_OUT_MAUVE)
		self.assertTrue(os.path.exists(expect_file))
		
		temp_dir = self.utils.get_temp_dir()
		for file_ in vect_files:
			self.assertTrue(os.path.exists(file_))
			self.utils.copy_file(file_, os.path.join(temp_dir, os.path.basename(file_)))
			
		out_file = self.utils.get_temp_file("mauve", ".xfma")
		os.unlink(out_file)
		out_file = os.path.join(temp_dir, os.path.basename(out_file))
		output_file = self.software.run_mauve(temp_dir, out_file)
		self.assertEquals(output_file, out_file)
		temp_file = self.utils.get_temp_file("file_name", ".txt")
		temp_file_1 = self.utils.get_temp_file("file_name", ".txt")
		cmd = "grep -E -v '^>|#BackboneFile' {} > {}".format(output_file, temp_file)
		os.system(cmd);
		cmd = "grep -E -v '^>|#BackboneFile' {} > {}".format(expect_file, temp_file_1)
		os.system(cmd);
		self.assertTrue(filecmp.cmp(temp_file, temp_file_1))
		os.unlink(temp_file)
		os.unlink(temp_file_1)
		self.utils.remove_dir(temp_dir)
		
	def test_run_convert_mauve(self):
		"""
 		test run mauve
 		create a VCF
 		"""
		mauve_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.FILE_OUT_MAUVE)
		expect_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.FILE_OUT_CONVERT_MAUVE_RESULT)
		self.assertTrue(os.path.exists(mauve_file))
		self.assertTrue(os.path.exists(expect_file))
		
		out_file = self.utils.get_temp_file("mauve", ".fasta")
		output_file = self.software.run_convert_mauve(mauve_file, out_file)
		self.assertEquals(output_file, out_file)
		temp_file = self.utils.get_temp_file("file_name", ".txt")
		temp_file_1 = self.utils.get_temp_file("file_name", ".txt")
		cmd = "grep -E -v '^>' {} > {}".format(output_file, temp_file)
		os.system(cmd);
		cmd = "grep -E -v '^>' {} > {}".format(expect_file, temp_file_1)
		os.system(cmd);
		self.assertTrue(filecmp.cmp(temp_file, temp_file_1))
		os.unlink(temp_file)
		os.unlink(temp_file_1)
		os.unlink(out_file)

	def test_run_mafft(self):
		"""
 		test run mafft
 		"""
		in_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.FILE_OUT_CONVERT_MAUVE_RESULT)
		expect_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.FILE_OUT_MAFFT_RESULT)
		self.assertTrue(os.path.exists(in_file))
		self.assertTrue(os.path.exists(expect_file))
		
		out_file = self.utils.get_temp_file("mafft", ".fasta")
		output_file = self.software.run_mafft(in_file, out_file, SoftwareNames.SOFTWARE_MAFFT_PARAMETERS)
		self.assertEquals(output_file, out_file)
		self.assertTrue(filecmp.cmp(out_file, expect_file))
		os.unlink(out_file)

	def test_run_mafft_2(self):
		"""
 		test run mafft
 		"""
	
		in_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.MANAGING_MAFFT_IN_PROTEIN)
		expect_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.MANAGING_MAFFT_IN_PROTEIN_EXPECTED)
		self.assertTrue(os.path.exists(in_file))
		self.assertTrue(os.path.exists(expect_file))
		
		out_file = self.utils.get_temp_file("mafft", ".faa")
		output_file = self.software.run_mafft(in_file, out_file, SoftwareNames.SOFTWARE_MAFFT_PARAMETERS_PROTEIN)
		self.assertEquals(output_file, out_file)
		self.assertTrue(filecmp.cmp(out_file, expect_file))
		os.unlink(out_file)
		
	def test_run_mafft_3(self):
		"""
 		test run mafft
 		"""
		in_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.MANAGING_MAFFT_IN)
		expect_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.MANAGING_MAFFT_IN_EXPECTED)
		self.assertTrue(os.path.exists(in_file))
		self.assertTrue(os.path.exists(expect_file))
		
		out_file = self.utils.get_temp_file("mafft", ".fasta")
		output_file = self.software.run_mafft(in_file, out_file, SoftwareNames.SOFTWARE_MAFFT_PARAMETERS_TWO_SEQUENCES)
		self.assertEquals(output_file, out_file)
		self.assertTrue(filecmp.cmp(out_file, expect_file))
		os.unlink(out_file)
		
	def test_run_seqret(self):
		"""
 		test run seqret
 		"""
		in_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.FILE_OUT_MAFFT_RESULT)
		expect_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.FILE_OUT_NEX_RESULT)
		self.assertTrue(os.path.exists(in_file))
		self.assertTrue(os.path.exists(expect_file))
		
		out_file = self.utils.get_temp_file("seqret", ".nex")
		output_file = self.software.run_seqret_nex(in_file, out_file)
		self.assertEquals(output_file, out_file)
		temp_file = self.utils.get_temp_file("file_name", ".txt")
		temp_file_1 = self.utils.get_temp_file("file_name", ".txt")
		cmd = "grep -E -v 'TITLE: Written by EMBOSS' {} > {}".format(out_file, temp_file)
		os.system(cmd);
		cmd = "grep -E -v 'TITLE: Written by EMBOSS' {} > {}".format(expect_file, temp_file_1)
		os.system(cmd);
		self.assertTrue(filecmp.cmp(temp_file, temp_file_1))
		os.unlink(temp_file)
		os.unlink(temp_file_1)
		os.unlink(out_file)

	def test_run_fasttree(self):
		"""
 		test run mauve
 		create a VCF
 		"""
		in_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.FILE_OUT_MAFFT_RESULT)
		expect_file_nwk = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.FILE_FASTTREE_RESULT_NWK)
		self.assertTrue(os.path.exists(expect_file_nwk))
		self.assertTrue(os.path.exists(in_file))
		
		out_file = self.utils.get_temp_file("fasttree", ".nwk")
		output_file = self.software.run_fasttree(in_file, out_file, self.software_names.get_fasttree_parameters())
		self.assertTrue(self.software_names.get_fasttree().endswith(self.software_names.get_fasttree_name()))
		self.assertEquals(output_file, out_file)
		self.assertTrue(filecmp.cmp(out_file, expect_file_nwk))
		os.unlink(out_file)

	### test single file in sample
	def test_identify_type_and_sub_type_single_file(self):
		"""
		get type and sub_type
		"""
		uploadFiles = UploadFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		self.assertEqual("2", version)
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
		uploadFiles.upload_file(version, file)	## upload file
		
		try:
			uploadFile = UploadFile.objects.order_by('-version')[0]
			self.assertEqual("test_db_influenza_typing_v2", uploadFile.abricate_name)
		except UploadFile.DoesNotExist:
			self.fail("must have values")
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()
			
		sample_name = "identify_type_and_sub_type"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = False
			sample.owner = user
			sample.save()
			
		return_value = self.software.identify_type_and_sub_type(sample, sample.path_name_1.name, None, user, True)
		self.assertTrue(return_value)
		
		vect_identify_virus = sample.identify_virus.all()
		self.assertEqual(3, len(vect_identify_virus))
		for identify_virus in vect_identify_virus:
			if (identify_virus.rank == 0):
				self.assertEquals("100.00", identify_virus.coverage)
				self.assertEquals("99.69", identify_virus.identity)
				self.assertEquals("A", identify_virus.seq_virus.name)
				self.assertEquals(ConstantsVirus.SEQ_VIRUS_TYPE, identify_virus.seq_virus.kind_type.name)
			elif (identify_virus.rank == 1):
				self.assertEquals("100.00", identify_virus.coverage)
				self.assertEquals("99.18", identify_virus.identity)
				self.assertEquals("H3", identify_virus.seq_virus.name)
				self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, identify_virus.seq_virus.kind_type.name)
			elif (identify_virus.rank == 2):
				self.assertEquals("100.00", identify_virus.coverage)
				self.assertEquals("98.65", identify_virus.identity)
				self.assertEquals("N2", identify_virus.seq_virus.name)
				self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, identify_virus.seq_virus.kind_type.name)
		file_abricate = sample.get_abricate_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(file_abricate))
		
		file_spades_contigs = sample.get_draft_contigs_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(file_spades_contigs))
		
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Identify_Sample, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Identify_Sample, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Spades(3.11.1), Abricate(0.8-dev)", list_meta[0].description)
		if (os.path.exists(file_abricate)): os.unlink(file_abricate)
		if (os.path.exists(file_spades_contigs)): os.unlink(file_spades_contigs)

		contigs_2_sequences = Contigs2Sequences(True)
		## remove abricate db
		cmd = "rm -r %s/%s*" % (SoftwareNames.SOFTWARE_ABRICATE_DB, contigs_2_sequences.get_database_name())
		exist_status = os.system(cmd)
		self.assertTrue(exist_status == 0)

	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_process_second_stage_analysis_single_file_3(self):
		"""
 		test global method
 		"""
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT_TEST", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))

		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_PV1_KX162693_gbk)
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_PV1_KX162693_fasta)
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		ref_name = "second_stagis_single_file3_pk"
		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.reference_fasta.name = fasta_file
			reference.reference_fasta_name = os.path.basename(fasta_file)
			reference.reference_genbank.name = gb_file
			reference.reference_genbank_name = os.path.basename(gb_file)
			reference.owner = user
			reference.save()
			
		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_1, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
		self.utils.copy_file(file_2, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2))
			
		sample_name = "run_snippyis_single_file3_pk"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = False
			sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
			sample.path_name_2.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2)
			sample.owner = user
			sample.save()

		project_name = "file_naais_single_filee_4_pk"
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
		
		### copy to trimmomatic files
		self.utils.copy_file(file_1, project_sample.sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True))
		self.utils.copy_file(file_2, project_sample.sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False))
		
		manageDatabase = ManageDatabase()
		metaKeyAndValue = MetaKeyAndValue()
		meta_key_project_sample = metaKeyAndValue.get_meta_key_queue_by_project_sample_id(project_sample.id)
		manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue, "meta_sample.description")
		self.assertTrue(self.software.process_second_stage_snippy_coverage_freebayes(project_sample, user))
		try:
			project_sample = ProjectSample.objects.get(pk=project_sample.id)
		except ProjectSample.DoesNotExist:
			self.fail("Must exist")
		self.assertTrue(project_sample.is_finished)
		
		### test the files
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertFalse(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FA, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, SoftwareNames.SOFTWARE_SNIPPY_name).index(SoftwareNames.SOFTWARE_SNIPPY_name.lower()) != -1)
		
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_BAM, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_BAM_BAI, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertFalse(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ_TBI, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CSV, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_SNIPPY_name)))
		## freebayes
		self.assertTrue(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_FREEBAYES_name).index(SoftwareNames.SOFTWARE_FREEBAYES_name.lower()) != -1)
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, SoftwareNames.SOFTWARE_FREEBAYES_name)))

		## test freebayes clean
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		self.assertTrue(os.path.getsize(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)) > 10)

		## test consensus file
		self.assertFalse(os.path.exists(project_sample.get_consensus_file(TypePath.MEDIA_ROOT)))

		### human file name, snippy tab
		self.assertTrue(os.path.exists(project_sample.get_file_output_human(TypePath.MEDIA_ROOT, FileType.FILE_TAB,  self.software_names.get_snippy_name())))
		self.assertTrue(os.path.getsize(project_sample.get_file_output_human(TypePath.MEDIA_ROOT, FileType.FILE_TAB,  self.software_names.get_snippy_name())) > 10)
		
		### set flag that is finished
		list_meta = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Count_Hits, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Count_Hits, list_meta[0].meta_tag.name)
		
		### get the hits value
		decode_coverage = DecodeObjects()
		count_hits = decode_coverage.decode_result(list_meta[0].description)
		self.assertEquals(0, count_hits.get_hits_50_90())
		self.assertEquals(0, count_hits.get_hits_less_50())
		self.assertEquals(0, count_hits.get_total_50_50_90())
		self.assertEquals(0, count_hits.get_total())
		
		self.assertEquals(0, project_sample.count_variations.var_less_50)
		self.assertEquals(0, project_sample.count_variations.var_bigger_50_90)
		self.assertEquals(0, project_sample.count_variations.var_bigger_90)
		self.assertEquals(0, project_sample.alert_first_level)
		self.assertEquals(1, project_sample.alert_second_level)
		
		### test mixed infections
		self.assertEquals('0.0', '{}'.format(project_sample.mixed_infections.average_value))
		self.assertEquals('No', project_sample.mixed_infections.tag.name)
		self.assertFalse(project_sample.mixed_infections.has_master_vector)
		
		list_meta = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Snippy_Freebayes, None)
		self.assertEquals(1, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Snippy_Freebayes, list_meta[0].meta_tag.name)
		
		decode_result = DecodeObjects()
		result = decode_result.decode_result(list_meta[0].description)
		self.assertTrue(result is not None)
		self.assertEquals("Snippy-3.2-dev; (--mapqual 20 --mincov 10 --minfrac 0.51)", result.get_software(self.software_names.get_snippy_name()))
		self.assertEquals("Freebayes-v1.1.0-54-g49413aa; (-p 2 -q 20 -m 20 --min-coverage 100 --min-alternate-fraction 0.01 --min-alternate-count 10 -V)",\
 						result.get_software(self.software_names.get_freebayes_name()))

		meta_value = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
		decode_coverage = DecodeObjects()
		coverage = decode_coverage.decode_result(meta_value.description)
		self.assertEqual(coverage.get_coverage('KX162693', Coverage.COVERAGE_ALL), "0.0")
		self.assertEqual(coverage.get_coverage('KX162693', Coverage.COVERAGE_MORE_9), "0.0")
		self.assertEqual(coverage.get_coverage('KX162693', Coverage.COVERAGE_MORE_0), "0.0")
		
		lst_meta_sample = manageDatabase.get_project_sample_metakey(project_sample, meta_key_project_sample, None)
		self.assertEquals(2, len(lst_meta_sample))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Queue, lst_meta_sample[1].value)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, lst_meta_sample[0].value)
		
		### check if the images exist
		for gene in self.utils.get_elements_from_db(reference, user):
			output_image = project_sample.get_global_file_by_element(TypePath.MEDIA_ROOT, ProjectSample.PREFIX_FILE_COVERAGE, gene, FileExtensions.FILE_PNG)
			self.assertTrue(os.path.exists(output_image))

		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))

	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_identify_type_and_sub_type_corona(self):
		"""
 		get type and sub_type
 		"""
		
		uploadFiles = UploadFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		self.assertEqual("2", version)
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
		uploadFiles.upload_file(version, file)	## upload file
		
		try:
			uploadFile = UploadFile.objects.order_by('-version')[0]
			self.assertEqual("test_db_influenza_typing_v2", uploadFile.abricate_name)
		except UploadFile.DoesNotExist:
			self.fail("must have values")
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()
			
		sample_name = "identify_type_and_sub_type"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = "WHU02_SRR10903401_1.fastq.gz"
			sample.path_name_1.name = os.path.join("/home/mmp/insa/corona", "WHU02_SRR10903401_1.fastq.gz")
			sample.is_valid_2 = True
			sample.file_name_2 = "WHU02_SRR10903401_2.fastq.gz"
			sample.path_name_2.name = os.path.join("/home/mmp/insa/corona", "WHU02_SRR10903401_2.fastq.gz")
			sample.owner = user
			sample.save()
			
		return_value = self.software.identify_type_and_sub_type(sample, sample.path_name_1.name, sample.path_name_2.name, user, True)
		self.assertTrue(return_value)
		
		vect_identify_virus = sample.identify_virus.all()
		self.assertEqual(2, len(vect_identify_virus))
		for identify_virus in vect_identify_virus:
			if (identify_virus.rank == 0):
				self.assertEquals("100.00", identify_virus.coverage)
				self.assertEquals("100.00", identify_virus.identity)
				self.assertEquals("BetaCoV", identify_virus.seq_virus.name)
				self.assertEquals("M_gene_MN908947", identify_virus.seq_virus.accession)
				self.assertEquals(ConstantsVirus.SEQ_VIRUS_GENUS, identify_virus.seq_virus.kind_type.name)
			elif (identify_virus.rank == 1):
				self.assertEquals("100.00", identify_virus.coverage)
				self.assertEquals("100.00", identify_virus.identity)
				self.assertEquals("2019_nCoV", identify_virus.seq_virus.name)
				self.assertEquals("S_gene_MN908947", identify_virus.seq_virus.accession)
				self.assertEquals(ConstantsVirus.SEQ_VIRUS_HUMAN, identify_virus.seq_virus.kind_type.name)
		file_abricate = sample.get_abricate_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(sample.get_abricate_output(TypePath.MEDIA_ROOT)))
		
		file_spades_contigs = sample.get_draft_contigs_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(file_spades_contigs))
		file_abricate_contigs = sample.get_draft_contigs_abricate_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(file_abricate_contigs))

		##### get type and subtype
		self.assertEquals("BetaCoV-2019_nCoV", sample.get_type_sub_type())
		
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Identify_Sample, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Identify_Sample, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Spades(3.11.1), Abricate(0.8-dev)", list_meta[0].description)
		if (os.path.exists(file_abricate)): os.unlink(file_abricate)
		if (os.path.exists(file_spades_contigs)): os.unlink(file_spades_contigs)
		if (os.path.exists(file_abricate_contigs)): os.unlink(file_abricate_contigs)

		contigs_2_sequences = Contigs2Sequences(True)
		## remove abricate db
		cmd = "rm -r %s/%s*" % (SoftwareNames.SOFTWARE_ABRICATE_DB, contigs_2_sequences.get_database_name())
		exist_status = os.system(cmd)
		self.assertTrue(exist_status == 0)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))

	def test_make_downsize(self):
		
		fastq1_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		fastq1_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)
		self.assertTrue(os.path.isfile(fastq1_1))
		self.assertTrue(os.path.isfile(fastq1_2))
		
		(is_downsized, file_name_1, file_name_2) = self.software.make_downsize(fastq1_1, fastq1_2, 2000000)	## downsize to 2M
		
		self.assertEquals(is_downsized, True)
		self.assertTrue(os.path.exists(file_name_1))
		self.assertTrue(os.path.exists(file_name_2))
		
		self.assertEquals(os.path.getsize(file_name_1), 233201)
		self.assertEquals(os.path.getsize(file_name_2), 255802)
		self.utils.remove_dir(os.path.dirname(file_name_1))

		(is_downsized, file_name_1, file_name_2) = self.software.make_downsize(fastq1_1, None, 2000000)	## downsize to 2M
		
		self.assertEquals(is_downsized, True)
		self.assertTrue(os.path.exists(file_name_1))
		self.assertTrue(file_name_2 == None)
		
		self.assertEquals(os.path.getsize(file_name_1), 248769)
		self.utils.remove_dir(os.path.dirname(file_name_1))
		
		(is_downsized, file_name_1, file_name_2) = self.software.make_downsize(fastq1_1, fastq1_2, 20000000)	## downsize to 2M
		self.assertEquals(is_downsized, False)

