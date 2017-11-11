'''
Created on Oct 28, 2017

@author: mmp
'''

from django.test import TestCase
from django.conf import settings 
from utils.constantsTestsCase import ConstantsTestsCase
from utils.constants import Constants
from utils.meta_key_and_values import MetaKeyAndValue
from utils.software import Software
from utils.utils import Utils
from utils.parseOutFiles import ParseOutFiles
from utils.result import DecodeResultAverageAndNumberReads, DecodeResult
from django.contrib.auth.models import User
from managing_files.models import Sample
from manage_virus.uploadFiles import UploadFiles
from manage_virus.models import UploadFile
from managing_files.manage_database import ManageDatabase
import os

class Test(TestCase):

	### static
	software = Software()
	utils = Utils()
	constants = Constants()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass


	def tearDown(self):
		pass


	def testCreateFaiToFastaFile(self):
		"""
		Test samtools fai index
		"""
		## create an index file from 
		
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		if os.path.isfile(fasta_file + ".fai"): os.unlink(fasta_file + ".fai")
		try:
			self.software.createFaiToFastaFile(fasta_file)
			self.assertTrue(os.path.isfile(fasta_file + ".fai"))
			self.assertEquals(179, os.path.getsize(fasta_file + ".fai"))
			os.unlink(fasta_file + ".fai")
		except Exception as e:
			self.fail(e.args[0])
		
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FAIL_FASTA)
		if os.path.isfile(fasta_file + ".fai"): os.unlink(fasta_file + ".fai")
		try:
			self.software.createFaiToFastaFile(fasta_file)
			self.fail("Must throw exception")
		except Exception as e:
			self.assertEquals("Fail to run samtools", e.args[0])
			if os.path.isfile(fasta_file + ".fai"): os.unlink(fasta_file + ".fai")

	def test_is_exist_database_abricate(self):
		self.assertTrue(self.software.is_exist_database_abricate("ncbi"))
		self.assertFalse(self.software.is_exist_database_abricate("xpot"))
	
	def test_create_database_abricate(self):
		database_name = "xpto"
		if (self.software.is_exist_database_abricate(database_name)):
			cmd = "rm -r %s/%s*" % (Software.SOFTWARE_ABRICATE_DB, database_name)
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
		cmd = "rm -r %s/%s*" % (Software.SOFTWARE_ABRICATE_DB, database_name)
		exist_status = os.system(cmd)
		self.assertTrue(exist_status == 0)
		self.assertFalse(self.software.is_exist_database_abricate(database_name))
		
	def testRunSpadesAndAbricate(self):
		"""
		Test run spades and abricator
		"""
		b_run_this = getattr(settings, "RUN_TEST_IN_COMMAND_LINE", None)
		
		if (not b_run_this): return	## it only work from command line because of PYTHONPATH defined on the eclispe IDE 
		fastq1_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		fastq1_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)
		self.assertTrue(os.path.isfile(fastq1_1))
		self.assertTrue(os.path.isfile(fastq1_2))
		
		out_dir = self.utils.get_temp_dir()
		self.assertTrue(os.path.isdir(out_dir))
		cmd = self.software.run_spades(fastq1_1, fastq1_2, out_dir)
		file_out = os.path.join(out_dir, "contigs.fasta")

#		file_out = os.path.join("/tmp/insaFlu/insa_flu_85008740", "contigs.fasta")
		self.assertTrue(os.path.exists(file_out))
		self.assertTrue(os.path.getsize(file_out) > 100)
		
		database_name = "xpto"
		if (not self.software.is_exist_database_abricate(database_name)):
			temp_file_file = os.path.join(self.baseDirectory, Constants.DIR_TYPE_IDENTIFICATION, ConstantsTestsCase.MANAGING_TEST_INFLUENZA_FILE)
			self.assertTrue(os.path.exists(temp_file_file))
			self.software.create_database_abricate(database_name, temp_file_file)
		
		out_file = self.utils.get_temp_file("temp_abricate", ".txt")
		cmd = self.software.run_abricate(database_name, file_out, out_file)
		self.assertTrue(os.path.exists(out_file))

		parseOutFiles = ParseOutFiles()
		vect_data = parseOutFiles.parse_abricate_file(out_file)
		self.assertEqual('A', vect_data[0][ParseOutFiles.GENE])
		self.assertEqual(100.00, vect_data[0][ParseOutFiles.COVERAGE])
		self.assertEqual(99.69, vect_data[0][ParseOutFiles.IDENTITY])
		self.assertEqual('influenza_type', vect_data[0][ParseOutFiles.TYPE])
		self.assertEqual('XXXX', vect_data[0][ParseOutFiles.ACCESSION])
		
		self.assertEqual('N2', vect_data[3][ParseOutFiles.GENE])
		self.assertEqual(16.03, vect_data[3][ParseOutFiles.COVERAGE])
		self.assertEqual(96.46, vect_data[3][ParseOutFiles.IDENTITY])
		self.assertEqual('influenza_subtype', vect_data[3][ParseOutFiles.TYPE])
		self.assertEqual('XXXXXX', vect_data[3][ParseOutFiles.ACCESSION])

		## remove file
		os.unlink(out_file)
		
		## remove abricate db
		cmd = "rm -r %s/%s*" % (Software.SOFTWARE_ABRICATE_DB, database_name)
		exist_status = os.system(cmd)
		self.assertTrue(exist_status == 0)

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
		
		out_put_path = self.software.run_fastq(file_1, file_2)
		
		out_file_1 = os.path.join(out_put_path, os.path.basename(self.constants.get_fastq_output(file_1)))
		out_file_2 = os.path.join(out_put_path, os.path.basename(self.constants.get_fastq_output(file_2)))
		self.assertTrue(os.path.exists(out_file_1))
		self.assertTrue(os.path.exists(out_file_2))
		self.assertTrue(os.path.getsize(out_file_1) > 1000)
		self.assertTrue(os.path.getsize(out_file_2) > 1000)
		self.utils.remove_dir(out_put_path)

		out_put_path = self.software.run_fastq(file_1, None)
		out_file_1 = os.path.join(out_put_path, os.path.basename(self.constants.get_fastq_output(file_1)))
		out_file_2 = os.path.join(out_put_path, os.path.basename(self.constants.get_fastq_output(file_2)))
		self.assertTrue(os.path.exists(out_file_1))
		self.assertFalse(os.path.exists(out_file_2))
		self.assertTrue(os.path.getsize(out_file_1) > 1000)
		self.utils.remove_dir(out_put_path)
		
	def test_run_trimmomatic(self):
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)

		sample_name = "sample_name"
		out_put_path_trimmomatic = self.software.run_trimmomatic(file_1, file_2, sample_name)
		b_first_file = True
		out_file_1 = os.path.join(out_put_path_trimmomatic, os.path.basename(self.constants.get_trimmomatic_output(out_put_path_trimmomatic, sample_name, b_first_file)))
		b_first_file = False
		out_file_2 = os.path.join(out_put_path_trimmomatic, os.path.basename(self.constants.get_trimmomatic_output(out_put_path_trimmomatic, sample_name, b_first_file)))
		self.assertTrue(os.path.exists(out_file_1))
		self.assertTrue(os.path.exists(out_file_2))
		self.assertTrue(os.path.getsize(out_file_1) > 1000)
		self.assertTrue(os.path.getsize(out_file_2) > 1000)
		
		out_put_path = self.software.run_fastq(out_file_1, out_file_2)
		b_first_file = True
		out_file_1 = os.path.join(out_put_path, os.path.basename(self.constants.get_fastq_trimmomatic_output(out_put_path, sample_name, b_first_file)))
		b_first_file = False
		out_file_2 = os.path.join(out_put_path, os.path.basename(self.constants.get_fastq_trimmomatic_output(out_put_path, sample_name, b_first_file)))
		print(out_file_1)
		self.assertTrue(os.path.exists(out_file_1))
		self.assertTrue(os.path.exists(out_file_2))
		self.assertTrue(os.path.getsize(out_file_1) > 1000)
		self.assertTrue(os.path.getsize(out_file_2) > 1000)
		self.utils.remove_dir(out_put_path)
		self.utils.remove_dir(out_put_path_trimmomatic)

		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		sample_name = "sample_name"
		out_put_path = self.software.run_trimmomatic(file_1, None, sample_name)
		b_first_file = True
		out_file_1 = os.path.join(out_put_path, os.path.basename(self.constants.get_trimmomatic_output(out_put_path, sample_name, b_first_file)))
		b_first_file = False
		out_file_2 = os.path.join(out_put_path, os.path.basename(self.constants.get_trimmomatic_output(out_put_path, sample_name, b_first_file)))
		self.assertTrue(os.path.exists(out_file_1))
		self.assertFalse(os.path.exists(out_file_2))
		self.assertTrue(os.path.getsize(out_file_1) > 1000)
		self.utils.remove_dir(out_put_path)

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
			sample.is_rejected = False
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = True
			sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
			sample.path_name_2.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2)
			sample.owner = user
			sample.save()
		
		### run software
		self.assertTrue(self.software.run_fastq_and_trimmomatic(sample, user))
		
		self.assertTrue(os.path.exists(self.constants.get_fastq_output(sample.path_name_1)))
		self.assertTrue(os.path.exists(self.constants.get_fastq_output(sample.path_name_2)))
		self.assertTrue(os.path.exists(self.constants.get_trimmomatic_output(temp_dir, sample.name, True)))
		self.assertTrue(os.path.exists(self.constants.get_trimmomatic_output(temp_dir, sample.name, False)))
		self.assertTrue(os.path.exists(self.constants.get_fastq_trimmomatic_output(temp_dir, sample.name, True)))
		self.assertTrue(os.path.exists(self.constants.get_fastq_trimmomatic_output(temp_dir, sample.name, False)))
		
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_metakey(sample, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Number_And_Average_Reads, list_meta[0].meta_tag.name)
		
		### number of sequences
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_metakey(sample, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Number_And_Average_Reads, list_meta[0].meta_tag.name)

		decodeResultAverageAndNumberReads = DecodeResultAverageAndNumberReads()
		result_average = decodeResultAverageAndNumberReads.decode_result(list_meta[0].description)
		self.assertEqual('39845', result_average.number_file_1)
		self.assertEqual('144.0', result_average.average_file_1)
		self.assertEqual('39845', result_average.number_file_2)
		self.assertEqual('142.1', result_average.average_file_2)

		list_meta = manageDatabase.get_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
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
			sample.is_rejected = False
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = False
			sample.owner = user
			sample.save()
		
		### run software
		self.assertTrue(self.software.run_fastq_and_trimmomatic(sample, user))
		
		self.assertTrue(os.path.exists(self.constants.get_fastq_output(sample.path_name_1)))
		self.assertTrue(os.path.exists(self.constants.get_trimmomatic_output(temp_dir, sample.name, True)))
		self.assertTrue(os.path.exists(self.constants.get_fastq_trimmomatic_output(temp_dir, sample.name, True)))
		
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_metakey(sample, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Number_And_Average_Reads, list_meta[0].meta_tag.name)
		
		### number of sequences
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_metakey(sample, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Number_And_Average_Reads, list_meta[0].meta_tag.name)

		decodeResultAverageAndNumberReads = DecodeResultAverageAndNumberReads()
		result_average = decodeResultAverageAndNumberReads.decode_result(list_meta[0].description)
		self.assertEqual('42534', result_average.number_file_1)
		self.assertEqual('143.4', result_average.average_file_1)
		self.assertEqual(None, result_average.number_file_2)
		self.assertEqual(None, result_average.average_file_2)

		list_meta = manageDatabase.get_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
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
			sample.is_rejected = False
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = True
			sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
			sample.path_name_2.name = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)
			sample.owner = user
			sample.save()
			
		return_value = self.software.identify_type_and_sub_type(sample, sample.path_name_1.name, sample.path_name_2.name, user)
		self.assertTrue(return_value)
		
		vect_identify_virus = sample.identify_virus.all()
		self.assertEqual(3, len(vect_identify_virus))
		for identify_virus in vect_identify_virus:
			if (identify_virus.rank == 0):
				self.assertEquals("100.00", identify_virus.coverage)
				self.assertEquals("99.69", identify_virus.identity)
				self.assertEquals("A", identify_virus.seq_virus.name)
				self.assertEquals(Constants.SEQ_VIRUS_TYPE, identify_virus.seq_virus.kind_type.name)
			elif (identify_virus.rank == 1):
				self.assertEquals("100.00", identify_virus.coverage)
				self.assertEquals("99.18", identify_virus.identity)
				self.assertEquals("H3", identify_virus.seq_virus.name)
				self.assertEquals(Constants.SEQ_VIRUS_SUB_TYPE, identify_virus.seq_virus.kind_type.name)
			elif (identify_virus.rank == 2):
				self.assertEquals("100.00", identify_virus.coverage)
				self.assertEquals("98.65", identify_virus.identity)
				self.assertEquals("N2", identify_virus.seq_virus.name)
				self.assertEquals(Constants.SEQ_VIRUS_SUB_TYPE, identify_virus.seq_virus.kind_type.name)
		file_abricate = self.constants.get_abricate_output(sample.path_name_1)
		self.assertTrue(os.path.exists(file_abricate))
		
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_metakey(sample, MetaKeyAndValue.META_KEY_Identify_Sample, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Identify_Sample, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Spades(3.11.1), Abricate(0.8-dev)", list_meta[0].description)
		if (os.path.exists(file_abricate)): os.unlink(file_abricate)


###
###			THIS doesn't work because spades don't run with one single file
###
# 	def test_identify_type_and_sub_type_single_file(self):
# 		"""
# 		get type and sub_type
# 		"""
# 		uploadFiles = UploadFiles()
# 		to_test = True
# 		(version, file) = uploadFiles.get_file_to_upload(to_test)
# 		self.assertEqual("2", version)
# 		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
# 		uploadFiles.upload_file(version, file)	## upload file
# 		
# 		try:
# 			uploadFile = UploadFile.objects.order_by('-version')[0]
# 			self.assertEqual("test_db_influenza_typing_v2", uploadFile.abricate_name)
# 		except UploadFile.DoesNotExist:
# 			self.fail("must have values")
# 		
# 		try:
# 			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
# 		except User.DoesNotExist:
# 			user = User()
# 			user.username = ConstantsTestsCase.TEST_USER_NAME
# 			user.is_active = False
# 			user.password = ConstantsTestsCase.TEST_USER_NAME
# 			user.save()
# 			
# 		sample_name = "identify_type_and_sub_type"
# 		try:
# 			sample = Sample.objects.get(name=sample_name)
# 		except Sample.DoesNotExist:
# 			sample = Sample()
# 			sample.name = sample_name
# 			sample.is_rejected = False
# 			sample.is_valid_1 = True
# 			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
# 			sample.path_name_1.name = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
# 			sample.is_valid_2 = False
# 			sample.owner = user
# 			sample.save()
# 			
# 		return_value = self.software.identify_type_and_sub_type(sample, sample.path_name_1.name, None, user)
# 		self.assertTrue(return_value)
# 		
# 		vect_identify_virus = sample.identify_virus.all()
# 		self.assertEqual(3, len(vect_identify_virus))
# 		for identify_virus in vect_identify_virus:
# 			print(identify_virus.coverage)
# 			print(identify_virus.identity)
# 			print(identify_virus.seq_virus.name)
# 			if (identify_virus.rank == 0):
# 				self.assertEquals("100.00", identify_virus.coverage)
# 				self.assertEquals("99.69", identify_virus.identity)
# 				self.assertEquals("A", identify_virus.seq_virus.name)
# 				self.assertEquals(Constants.SEQ_VIRUS_TYPE, identify_virus.seq_virus.kind_type.name)
# 			elif (identify_virus.rank == 1):
# 				self.assertEquals("100.00", identify_virus.coverage)
# 				self.assertEquals("99.18", identify_virus.identity)
# 				self.assertEquals("H3", identify_virus.seq_virus.name)
# 				self.assertEquals(Constants.SEQ_VIRUS_SUB_TYPE, identify_virus.seq_virus.kind_type.name)
# 			elif (identify_virus.rank == 2):
# 				self.assertEquals("100.00", identify_virus.coverage)
# 				self.assertEquals("98.65", identify_virus.identity)
# 				self.assertEquals("N2", identify_virus.seq_virus.name)
# 				self.assertEquals(Constants.SEQ_VIRUS_SUB_TYPE, identify_virus.seq_virus.kind_type.name)
# 		file_abricate = self.constants.get_abricate_output(sample.path_name_1)
# 		self.assertTrue(os.path.exists(file_abricate))
# 		
# 		manageDatabase = ManageDatabase()
# 		list_meta = manageDatabase.get_metakey(sample, MetaKeyAndValue.META_KEY_Identify_Sample, None)
# 		self.assertTrue(len(list_meta) == 1)
# 		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
# 		self.assertEquals(MetaKeyAndValue.META_KEY_Identify_Sample, list_meta[0].meta_tag.name)
# 		self.assertEquals("Success, Spades(3.11.1), Abricate(0.8-dev)", list_meta[0].description)
# 		if (os.path.exists(file_abricate)): os.unlink(file_abricate)
