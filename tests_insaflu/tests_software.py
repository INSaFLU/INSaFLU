'''
Created on Oct 28, 2017

@author: mmp
'''

from pickle import TRUE
from django.test import TestCase
from django.conf import settings 
from constants.constantsTestsCase import ConstantsTestsCase
from constants.constants_mixed_infection import ConstantsMixedInfection
from constants.constants import Constants, TypePath, FileType, FileExtensions
from constants.meta_key_and_values import MetaKeyAndValue
from utils.software import Software, Contigs2Sequences
from utils.software_pangolin import SoftwarePangolin
from constants.software_names import SoftwareNames
from utils.utils import Utils
from utils.parse_out_files import ParseOutFiles
from utils.result import DecodeObjects, Coverage, MaskingConsensus, GeneticElement
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference
from manage_virus.uploadFiles import UploadFiles
from manage_virus.models import UploadFile
from managing_files.manage_database import ManageDatabase
from django.test.utils import override_settings
from manage_virus.constants_virus import ConstantsVirus
from settings.default_software import DefaultSoftware
from settings.default_software_project_sample import DefaultProjectSoftware
from settings.models import Software as SoftwareSettings, Parameter
from utils.parse_coverage_file import GetCoverage
from managing_files.models import Software as SoftwareModel
from settings.constants_settings import ConstantsSettings
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os, filecmp
import pandas

class Test(TestCase):

	### static
	software = Software()
	software_pangolin = SoftwarePangolin()
	utils = Utils()
	constants = Constants()
	software_names = SoftwareNames()
	constants_tests_case = ConstantsTestsCase()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)


	def tearDown(self):
		pass


	def test_get_species_tag(self):
		"""
		test SARS cov 
		"""
		
		ref_name = "File 1"
		reference_fasta = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_COVID_FASTA)
		ref_name_2 = "File 2"
		consensus_file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "SARSCoVDec200153.consensus.fasta")
		ref_name_3 = "File 3"
		consensus_file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "A_H3N2_A_Hong_Kong_4801_2014.fasta")
		ref_name_4 = "File 4"
		consensus_file_3 = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "monkeypox_MT903344_MPXV_UK_P2_wo_NN.fasta")
		
		self.assertTrue(os.path.exists(reference_fasta))
		self.assertTrue(os.path.exists(consensus_file_1))
		self.assertTrue(os.path.exists(consensus_file_2))
		self.assertTrue(os.path.exists(consensus_file_3))
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.reference_fasta.name = reference_fasta
			reference.reference_fasta_name = os.path.basename(reference_fasta)
			reference.owner = user
			reference.save()
			
		try:
			reference1 = Reference.objects.get(name=ref_name_2)
		except Reference.DoesNotExist:
			reference1 = Reference()
			reference1.name = ref_name_2
			reference1.reference_fasta.name = consensus_file_1
			reference1.reference_fasta_name = os.path.basename(consensus_file_1)
			reference1.owner = user
			reference1.save()
		
		try:
			reference2 = Reference.objects.get(name=ref_name_3)
		except Reference.DoesNotExist:
			reference2 = Reference()
			reference2.name = ref_name_3
			reference2.reference_fasta.name = consensus_file_2
			reference2.reference_fasta_name = os.path.basename(consensus_file_2)
			reference2.owner = user
			reference2.save()
			
		try:
			reference3 = Reference.objects.get(name=ref_name_4)
		except Reference.DoesNotExist:
			reference3 = Reference()
			reference3.name = ref_name_3
			reference3.reference_fasta.name = consensus_file_3
			reference3.reference_fasta_name = os.path.basename(consensus_file_3)
			reference3.owner = user
			reference3.save()
		
		uploadFiles = UploadFiles()
		version = 1
		file = os.path.join(self.baseDirectory, Constants.DIR_TYPE_IDENTIFICATION, "test_covid_typing.fasta")
		self.assertTrue(os.path.exists(file))
		uploadFiles.upload_file(version, file)	## upload file
		
		database_name = "xpto_sars_cov"
		if (not self.software.is_exist_database_abricate(database_name)):
			self.software.create_database_abricate(database_name, file)

		self.assertEquals(Reference.SPECIES_INFLUENZA, self.software.get_species_tag(reference2))		
		self.assertEquals(Reference.SPECIES_MPXV, self.software.get_species_tag(reference3))
		self.assertEquals(Reference.SPECIES_SARS_COV_2, self.software.get_species_tag(reference))
		self.assertEquals(Reference.SPECIES_SARS_COV_2, self.software.get_species_tag(reference1))
		
		##
		try:
			reference_test = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			self.assertFail("Record must exist")
		
		self.assertEqual(Reference.SPECIES_SARS_COV_2, reference_test.specie_tag)
		

	def test_fasta_2_upper(self):
		"""
		Test samtools fai index
		"""
		## create an index file from 
		
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "A_H3N2_reference_demo_lower_case.fasta")
		fasta_file_temp = self.utils.get_temp_file("fasta_lower", FileExtensions.FILE_FASTA)
		self.utils.copy_file(fasta_file, fasta_file_temp)
		self.software.fasta_2_upper(fasta_file_temp)
		fasta_file_upper = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "A_H3N2_reference_demo_upper_case.fasta")
		self.assertTrue(filecmp.cmp(fasta_file_temp, fasta_file_upper))
		self.utils.remove_file(fasta_file_temp)
		
	def test_get_first_sequence_fasta(self):
		"""
		Test samtools fai index
		"""
		## create an index file from 
		
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "A_H3N2_reference_demo_lower_case.fasta")
		fasta_file_temp = self.utils.get_temp_file("fasta_single_read", FileExtensions.FILE_FASTA)
		self.utils.copy_file(fasta_file, fasta_file_temp)
		self.software.set_first_sequence_fasta(fasta_file_temp)
		fasta_file_upper = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "A_H3N2_reference_demo_upper_case.fasta")
		self.assertTrue(filecmp.cmp(fasta_file_temp, fasta_file_upper))
		self.utils.remove_file(fasta_file_temp)

	def testCreateFaiToFastaFile(self):
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
		
		database_name = "xpto_insa"
		if (not self.software.is_exist_database_abricate(database_name)):
			temp_file_file = os.path.join(self.baseDirectory, Constants.DIR_TYPE_IDENTIFICATION, ConstantsTestsCase.MANAGING_TEST_INFLUENZA_FILE)
			self.assertTrue(os.path.exists(temp_file_file))
			self.software.create_database_abricate(database_name, temp_file_file)
		
		out_file = self.utils.get_temp_file("temp_abricate", ".txt")
		cmd = self.software.run_abricate(database_name, file_out, SoftwareNames.SOFTWARE_ABRICATE_PARAMETERS, out_file)
		self.assertTrue(os.path.exists(out_file))

		parseOutFiles = ParseOutFiles()
		(dict_out_abricate, clean_abricate_file) = parseOutFiles.parse_abricate_file(out_file, 'temp_2.txt', 1)
		self.assertEqual('A', dict_out_abricate['NODE_1_length_3445_cov_491.652019'][0][ParseOutFiles.GENE])
		self.assertEqual(100.00, dict_out_abricate['NODE_1_length_3445_cov_491.652019'][0][ParseOutFiles.COVERAGE])
		self.assertEqual(99.69, dict_out_abricate['NODE_1_length_3445_cov_491.652019'][0][ParseOutFiles.IDENTITY])
		self.assertEqual('influenza_type', dict_out_abricate['NODE_1_length_3445_cov_491.652019'][0][ParseOutFiles.TYPE])
		self.assertEqual('XXXX', dict_out_abricate['NODE_1_length_3445_cov_491.652019'][0][ParseOutFiles.ACCESSION])
		
		self.assertEqual('N2', dict_out_abricate['NODE_7_length_1478_cov_488.215560'][0][ParseOutFiles.GENE])
		self.assertEqual(100.0, dict_out_abricate['NODE_7_length_1478_cov_488.215560'][0][ParseOutFiles.COVERAGE])
		self.assertEqual(98.65, dict_out_abricate['NODE_7_length_1478_cov_488.215560'][0][ParseOutFiles.IDENTITY])
		self.assertEqual('influenza_subtype', dict_out_abricate['NODE_7_length_1478_cov_488.215560'][0][ParseOutFiles.TYPE])
		self.assertEqual('XXXXXX', dict_out_abricate['NODE_7_length_1478_cov_488.215560'][0][ParseOutFiles.ACCESSION])

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
		out_file_1 = os.path.join(out_put_path, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, True)))
		out_file_2 = os.path.join(out_put_path, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, False)))
		self.assertTrue(os.path.exists(out_file_1))
		self.assertTrue(os.path.exists(out_file_2))
		self.assertTrue(os.path.getsize(out_file_1) > 1000)
		self.assertTrue(os.path.getsize(out_file_2) > 1000)
		self.utils.remove_dir(out_put_path)

		out_put_path = self.software.run_fastq(file_1, None)
		out_file_1 = os.path.join(out_put_path, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, True)))
		out_file_2 = os.path.join(out_put_path, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, False)))
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

		(out_put_path_trimmomatic, filtering_result, parameters) = self.software.run_trimmomatic(sample.get_fastq(TypePath.MEDIA_ROOT, True),
						sample.get_fastq(TypePath.MEDIA_ROOT, False), sample)
		out_file_1 = os.path.join(out_put_path_trimmomatic, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)))
		out_file_2 = os.path.join(out_put_path_trimmomatic, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False)))
		self.assertTrue(os.path.exists(out_file_1))
		self.assertTrue(os.path.exists(out_file_2))
		self.assertTrue(os.path.getsize(out_file_1) > 1000)
		self.assertTrue(os.path.getsize(out_file_2) > 1000)
		self.assertEqual("SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33", parameters)
		self.assertEqual("Quality encoding detected:", filtering_result.get_key_value()[0].key)
		self.assertEqual("phred33", filtering_result.get_key_value()[0].value)
		self.assertEqual("Input Read Pairs:", filtering_result.get_key_value()[1].key)
		self.assertEqual("44425", filtering_result.get_key_value()[1].value)
		self.assertEqual("Dropped:", filtering_result.get_key_value()[-1].key)
		self.assertEqual("434 (0.98%)", filtering_result.get_key_value()[-1].value)
		self.assertNotEqual("434 (0.98%)_", filtering_result.get_key_value()[-1].value)
		
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
		
		(out_put_path, filtering_result, parameters) = self.software.run_trimmomatic(sample.get_fastq(TypePath.MEDIA_ROOT, True),
							sample.get_fastq(TypePath.MEDIA_ROOT, False), sample)
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
			sample.type_of_fastq = Sample.TYPE_OF_FASTQ_illumina
			sample.save()
		
		### set the job
		taskID = "xpto_task" 
		manageDatabase = ManageDatabase()
		manageDatabase.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue, taskID)

		### run software
		self.assertEqual((True, True), self.software.run_fastq_and_trimmomatic(sample, user))

		default_software = DefaultSoftware()
		self.assertTrue(default_software.is_software_to_run(SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
							user, ConstantsSettings.TECHNOLOGY_illumina))		
		
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)))))
		
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
		
		### get key values of data
		stats_illumina, total_reads = self.software.get_stats_from_sample_reads(sample)
		self.assertEqual('Number of reads R1', stats_illumina[0][0])
		self.assertEqual('44,425', stats_illumina[0][1])
		self.assertEqual(82508, total_reads)
		self.assertEqual('41,254', stats_illumina[0][2])
		self.assertEqual('-3,171.0', stats_illumina[0][3])
		self.assertEqual('92.9', stats_illumina[0][4])
		self.assertEqual('Average read length R1', stats_illumina[1][0])
		self.assertEqual('143.6', stats_illumina[1][1])
		self.assertEqual('142.0', stats_illumina[1][2])
		self.assertEqual('-1.6', stats_illumina[1][3])
		self.assertEqual('--', stats_illumina[1][4])
		self.assertEqual('Total reads', stats_illumina[-1][0])
		self.assertEqual('88,850', stats_illumina[-1][1])
		self.assertEqual('82,508', stats_illumina[-1][2])
		self.assertEqual('-6,342.0', stats_illumina[-1][3])
		self.assertEqual('92.9', stats_illumina[-1][4])
		
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Fastq(0.11.9), Trimmomatic(0.39)", list_meta[0].description)
		
		
		decodeResult = DecodeObjects()
		meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software, MetaKeyAndValue.META_VALUE_Success)
		result = decodeResult.decode_result(meta_sample.description)
		self.assertTrue(result.is_software_present(SoftwareNames.SOFTWARE_TRIMMOMATIC_name))
		self.assertTrue(result.is_software_present(SoftwareNames.SOFTWARE_FASTQ_name))
		## remove all files
		cmd = "rm -r %s*" % (temp_dir); os.system(cmd)

	def test_run_fastq_and_trimmomatic_2(self):
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
			
		sample_name = "run_fastq_and_trimmomatic_1"
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
			sample.type_of_fastq = Sample.TYPE_OF_FASTQ_illumina
			sample.owner = user
			sample.save()
		
		### set the job
		taskID = "xpto_task" 
		manageDatabase = ManageDatabase()
		manageDatabase.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue, taskID)

		### turn off trimmomatic
		default_software = DefaultSoftware()
		default_software.test_all_defaults(user)
		self.assertTrue(default_software.is_software_to_run(SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
							user, ConstantsSettings.TECHNOLOGY_illumina))
		is_to_run = False
		default_software.set_software_to_run(SoftwareNames.SOFTWARE_TRIMMOMATIC_name, user,
							ConstantsSettings.TECHNOLOGY_illumina, is_to_run)
		self.assertFalse(default_software.is_software_to_run(SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
							user, ConstantsSettings.TECHNOLOGY_illumina))
		
		### run software
		self.assertEqual((True, False), self.software.run_fastq_and_trimmomatic(sample, user))
		default_project_software = DefaultProjectSoftware()
		self.assertFalse(default_project_software.is_to_run_trimmomatic(user, sample))
		
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
		self.assertFalse(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False)))))
		self.assertFalse(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, True)))))
		self.assertFalse(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)))))
		self.assertFalse(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)))))
		
		manageDatabase = ManageDatabase()
		decodeResultAverageAndNumberReads = DecodeObjects()
		meta_value = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, meta_value.value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Number_And_Average_Reads, meta_value.meta_tag.name)
		result_average = decodeResultAverageAndNumberReads.decode_result(meta_value.description)
		self.assertEqual('44425', result_average.number_file_1)
		self.assertEqual('143.6', result_average.average_file_1)
		self.assertEqual('44425', result_average.number_file_2)
		self.assertEqual('143.6', result_average.average_file_2)

		### get key values of data
		stats_illumina, total_reads = self.software.get_stats_from_sample_reads(sample)
		self.assertEqual('Number of reads R1', stats_illumina[0][0])
		self.assertEqual('44,425', stats_illumina[0][1])
		self.assertEqual(88850, total_reads)
		self.assertEqual('--', stats_illumina[0][2])
		self.assertEqual('0.0', stats_illumina[0][3])
		self.assertEqual('100', stats_illumina[0][4])
		self.assertEqual('Average read length R1', stats_illumina[1][0])
		self.assertEqual('143.6', stats_illumina[1][1])
		self.assertEqual('--', stats_illumina[1][2])
		self.assertEqual('0.0', stats_illumina[1][3])
		self.assertEqual('--', stats_illumina[1][4])
		self.assertEqual('Total reads', stats_illumina[-1][0])
		self.assertEqual('88,850', stats_illumina[-1][1])
		self.assertEqual('--', stats_illumina[-1][2])
		self.assertEqual('0.0', stats_illumina[-1][3])
		self.assertEqual('100', stats_illumina[-1][4])
		
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
		self.assertEqual(1, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Fastq(0.11.9)", list_meta[0].description)
		
		decodeResult = DecodeObjects()
		meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software, MetaKeyAndValue.META_VALUE_Success)
		result = decodeResult.decode_result(meta_sample.description)
		self.assertFalse(result.is_software_present(SoftwareNames.SOFTWARE_TRIMMOMATIC_name))
		self.assertTrue(result.is_software_present(SoftwareNames.SOFTWARE_FASTQ_name))
		
		## remove all files
		cmd = "rm -r %s*" % (temp_dir); os.system(cmd)
	
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
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina)
			sample.owner = user
			sample.save()
		
		self.assertFalse(sample.is_ready_for_projects)
		
		### set the job
		taskID = "xpto_task" 
		manageDatabase = ManageDatabase()
		manageDatabase.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue, taskID)

		### turn off trimmomatic
		default_software = DefaultSoftware()
		default_software.test_all_defaults(user)
		self.assertTrue(default_software.is_software_to_run(SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
							user, ConstantsSettings.TECHNOLOGY_illumina))
		is_to_run = False
		default_software.set_software_to_run(SoftwareNames.SOFTWARE_TRIMMOMATIC_name, user,
							ConstantsSettings.TECHNOLOGY_illumina, is_to_run)
		self.assertFalse(default_software.is_software_to_run(SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
							user, ConstantsSettings.TECHNOLOGY_illumina))
		
		### run software
		self.assertTrue(self.software.run_fastq_and_trimmomatic_and_identify_species(sample, user))
		default_project_software = DefaultProjectSoftware()
		self.assertFalse(default_project_software.is_to_run_trimmomatic(user, sample))
		
		meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_sample != None)
		self.assertEquals(taskID, meta_sample.description)

		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
		self.assertFalse(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False)))))
		self.assertFalse(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, True)))))
		self.assertFalse(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)))))
		self.assertFalse(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)))))
		
		### number of sequences
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Number_And_Average_Reads, list_meta[0].meta_tag.name)

		decodeResultAverageAndNumberReads = DecodeObjects()
		result_average = decodeResultAverageAndNumberReads.decode_result(list_meta[0].description)
		self.assertEqual('44425', result_average.number_file_1)
		self.assertEqual('143.6', result_average.average_file_1)
		self.assertEqual('44425', result_average.number_file_2)
		self.assertEqual('143.6', result_average.average_file_2)

		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Fastq(0.11.9)", list_meta[0].description)
		
		meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software, MetaKeyAndValue.META_VALUE_Success)
		decodeResult = DecodeObjects()
		result = decodeResult.decode_result(meta_sample.description)
		soft_desc = result.get_software_instance(SoftwareNames.SOFTWARE_TRIMMOMATIC_name)
		self.assertTrue(soft_desc is None)
		
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
		self.assertEqual((True, True), self.software.run_fastq_and_trimmomatic(sample, user))
		
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)))))
		
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
		self.assertEquals("Success, Fastq(0.11.9), Trimmomatic(0.39)", list_meta[0].description)
		
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
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina)
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
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina)
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
	def test_run_fastq_and_trimmomatic_and_identify_species_other(self):
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
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina)
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
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)))))
		
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
		self.assertEquals("Success, Fastq(0.11.9), Trimmomatic(0.39)", list_meta[0].description)
		
		meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software, MetaKeyAndValue.META_VALUE_Success)
		decodeResult = DecodeObjects()
		result = decodeResult.decode_result(meta_sample.description)
		soft_desc = result.get_software_instance(SoftwareNames.SOFTWARE_TRIMMOMATIC_name)
		self.assertEqual(SoftwareNames.SOFTWARE_TRIMMOMATIC_name, soft_desc.name)
		self.assertEqual("Quality encoding detected:", soft_desc.get_vect_key_values()[0].key)
		self.assertEqual("phred33", soft_desc.get_vect_key_values()[0].value)
		self.assertEqual("Input Read Pairs:", soft_desc.get_vect_key_values()[1].key)
		self.assertEqual("44425", soft_desc.get_vect_key_values()[1].value)
		self.assertEqual("Dropped:", soft_desc.get_vect_key_values()[-1].key)
		self.assertEqual("434 (0,98%)", soft_desc.get_vect_key_values()[-1].value)
		self.assertNotEqual("434 (0,98%)_", soft_desc.get_vect_key_values()[-1].value)
		
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
	def test_run_fastq_and_trimmomatic_and_identify_species_no_abricate(self):
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
			
		sample_name = "run_fastq_and_trimmomatic_no_abricate"
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
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina)
			sample.owner = user
			sample.save()
		
		self.assertFalse(sample.is_ready_for_projects)
		
		### set the job
		taskID = "xpto_task" 
		manageDatabase = ManageDatabase()
		manageDatabase.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue, taskID)

		### turn off abricate
		default_software = DefaultSoftware()
		default_software.test_all_defaults(user)
		
		default_project_software = DefaultProjectSoftware()
		self.assertTrue(default_software.is_software_to_run(SoftwareNames.SOFTWARE_ABRICATE_name,
							user, ConstantsSettings.TECHNOLOGY_illumina))
		is_to_run = False
		default_software.set_software_to_run(SoftwareNames.SOFTWARE_ABRICATE_name, user,
							ConstantsSettings.TECHNOLOGY_illumina, is_to_run)
		
		### run software
		self.assertTrue(self.software.run_fastq_and_trimmomatic_and_identify_species(sample, user))
		self.assertFalse(default_project_software.is_to_run_abricate(user, sample,
				ConstantsSettings.TECHNOLOGY_illumina))
		
		meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_sample != None)
		self.assertEquals(taskID, meta_sample.description)

		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)))))
		
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
		self.assertEquals("Success, Fastq(0.11.9), Trimmomatic(0.39)", list_meta[0].description)
		
		meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software, MetaKeyAndValue.META_VALUE_Success)
		decodeResult = DecodeObjects()
		result = decodeResult.decode_result(meta_sample.description)
		soft_desc = result.get_software_instance(SoftwareNames.SOFTWARE_TRIMMOMATIC_name)
		self.assertEqual(SoftwareNames.SOFTWARE_TRIMMOMATIC_name, soft_desc.name)
		self.assertEqual("Quality encoding detected:", soft_desc.get_vect_key_values()[0].key)
		self.assertEqual("phred33", soft_desc.get_vect_key_values()[0].value)
		self.assertEqual("Input Read Pairs:", soft_desc.get_vect_key_values()[1].key)
		self.assertEqual("44425", soft_desc.get_vect_key_values()[1].value)
		self.assertEqual("Dropped:", soft_desc.get_vect_key_values()[-1].key)
		self.assertEqual("434 (0,98%)", soft_desc.get_vect_key_values()[-1].value)
		self.assertNotEqual("434 (0,98%)_", soft_desc.get_vect_key_values()[-1].value)
		
		sample = Sample.objects.get(pk=sample.id)
		self.assertTrue(sample.is_ready_for_projects)
		self.assertTrue(sample.get_is_ready_for_projects())
		
		self.assertFalse(os.path.exists(sample.get_abricate_output(TypePath.MEDIA_ROOT)))
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Identify_Sample, None)
		self.assertTrue(len(list_meta) == 0)
		
		### mixed infections
		self.assertEquals(0, sample.number_alerts)
		self.assertEqual(Constants.EMPTY_VALUE_NA, sample.mixed_infections_tag.name)
		self.assertEqual(Constants.EMPTY_VALUE_NA, sample.type_subtype)
		
		manage_database = ManageDatabase()
		meta_data = manage_database.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE,\
								MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_data != None)
		self.assertEqual("Info: Abricate turned OFF by the user.", meta_data.description)
		
		meta_data = manage_database.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_TAG_MIXED_INFECTION_TYPE_SUBTYPE,\
								MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_data is None)
		
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
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina)
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
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)))))
		
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
		self.assertEquals("Success, Fastq(0.11.9), Trimmomatic(0.39)", list_meta[0].description)
		
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
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina)
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
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)))))
		
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
		self.assertEquals("Success, Fastq(0.11.9), Trimmomatic(0.39)", list_meta[0].description)
		
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
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina)
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
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastqc_output(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)))))
		
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
		self.assertEquals("Success, Fastq(0.11.9), Trimmomatic(0.39)", list_meta[0].description)
		
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

		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
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
		
		software_names = SoftwareNames()
		out_put_path = self.software.run_snippy(sample.get_fastq(TypePath.MEDIA_ROOT, True), sample.get_fastq(TypePath.MEDIA_ROOT, False),
				fasta_file, gb_file, sample.name, software_names.get_snippy_parameters())
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

		out_put_path = self.software.run_snippy(sample.get_fastq(TypePath.MEDIA_ROOT, True), None, fasta_file, gb_file, sample.name, software_names.get_snippy_parameters())
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
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina) 
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
		self.assertTrue(project_sample.is_mask_consensus_sequences)
		
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
		
		expected_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "run_snippyis_single_file1.tab")
		self.assertTrue(filecmp.cmp(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_SNIPPY_name),
				expected_file))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_SNIPPY_name)))
		
		## freebayes
		self.assertTrue(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_FREEBAYES_name).index(SoftwareNames.SOFTWARE_FREEBAYES_name.lower()) != -1)
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		expected_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "run_freebayes_single_file1.tab")
		self.assertTrue(filecmp.cmp(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name),
				expected_file))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, SoftwareNames.SOFTWARE_FREEBAYES_name)))

		## test freebayes clean
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		self.assertTrue(os.path.getsize(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)) > 10)

		## test consensus file
		self.assertTrue(os.path.exists(project_sample.get_consensus_file(TypePath.MEDIA_ROOT)))
		self.assertTrue(os.path.getsize(project_sample.get_consensus_file(TypePath.MEDIA_ROOT)) > 10)
		self.assertTrue(filecmp.cmp(project_sample.get_consensus_file(TypePath.MEDIA_ROOT),
							project_sample.get_backup_consensus_file()))

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
		self.assertEquals("Freebayes-v1.1.0-54-g49413aa; (--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10 --min-alternate-fraction 0.01 --ploidy 2 -V)",\
 						result.get_software(self.software_names.get_freebayes_name()))
		self.assertTrue(result.is_software_present(self.software_names.get_freebayes_name()))
		self.assertFalse(result.is_software_present(self.software_names.get_trimmomatic_name()))
		
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

		meta_value = manageDatabase.get_project_sample_metakey_last(project_sample,
						MetaKeyAndValue.META_KEY_bam_stats,
						MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(not meta_value is None)
		decode_masking_vcf = DecodeObjects()
		result = decode_masking_vcf.decode_result(meta_value.description)
		self.assertEqual('87818', result.get_value_by_key(MetaKeyAndValue.SAMTOOLS_flagstat_mapped_reads))
		self.assertEqual('87818', result.get_value_by_key(MetaKeyAndValue.SAMTOOLS_flagstat_total_reads))
		self.assertEqual(None, result.get_value_by_key('xpto'))
		
		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))


	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_process_second_stage_analysis_single_file_no_freebayes(self):
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
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina) 
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
		
		### copy to trimmomatic files, don't copy trimmommatic files to project
		##self.utils.copy_file(file_1, project_sample.sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True))
		##self.utils.copy_file(file_2, project_sample.sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False))
		
		### turn off freebayes
		default_software = DefaultSoftware()
		default_software.test_all_defaults(user)
		self.assertTrue(default_software.is_software_to_run(SoftwareNames.SOFTWARE_FREEBAYES_name,
							user, ConstantsSettings.TECHNOLOGY_illumina))
		### turn off freebayes
		is_to_run = False
		default_software.set_software_to_run(SoftwareNames.SOFTWARE_FREEBAYES_name, user,
							ConstantsSettings.TECHNOLOGY_illumina, is_to_run)
		self.assertFalse(default_software.is_software_to_run(SoftwareNames.SOFTWARE_FREEBAYES_name,
							user, ConstantsSettings.TECHNOLOGY_illumina))
		### set mincov 30 for snippy
		try:
			software = SoftwareSettings.objects.get(name=SoftwareNames.SOFTWARE_SNIPPY_name, owner=user,
							type_of_use = SoftwareSettings.TYPE_OF_USE_global,
							technology__name = ConstantsSettings.TECHNOLOGY_illumina)
			self.assertTrue(software.is_used_in_global())
		except SoftwareSettings.DoesNotExist:
			self.fail("Must exist this software name")
		
		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(3, len(parameters))
		
		### test set default
		self.assertEqual("10", parameters[1].parameter)
		parameter = parameters[1]
		parameter.parameter = 5
		parameter.save()
		self.assertEqual("--mapqual 20 --mincov 5 --minfrac 0.51", default_software.get_snippy_parameters(user))
		
		default_software_project = DefaultProjectSoftware()
		default_software_project.test_all_defaults(user, project, None, None) ## the user can have defaults yet		
		default_software_project.test_default_db(SoftwareNames.SOFTWARE_FREEBAYES_name, user,
	 			SoftwareSettings.TYPE_OF_USE_project_sample,
	 			None, project_sample, None, ConstantsSettings.TECHNOLOGY_illumina)
	
		## test if it is necessary to run freebayes
		self.assertFalse(default_software_project.is_to_run_freebayes(user, project_sample))
				
		manageDatabase = ManageDatabase()
		metaKeyAndValue = MetaKeyAndValue()
		meta_key_project_sample = metaKeyAndValue.get_meta_key_queue_by_project_sample_id(project_sample.id)
		manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue, "meta_sample.description")
		self.assertTrue(self.software.process_second_stage_snippy_coverage_freebayes(project_sample, user))
		## test if freebays run
		default_project_software = DefaultProjectSoftware()
		self.assertFalse(default_project_software.is_to_run_freebayes(user, project_sample))
		
		try:
			project_sample = ProjectSample.objects.get(pk=project_sample.id)
		except ProjectSample.DoesNotExist:
			self.fail("Must exist")
		self.assertTrue(project_sample.is_finished)
		self.assertTrue(project_sample.is_mask_consensus_sequences)
		
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
		
		expected_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "run_snippyis_single_file1.tab")
		self.assertTrue(filecmp.cmp(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_SNIPPY_name),
				expected_file))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_SNIPPY_name)))
		
		## freebayes
		self.assertTrue(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_FREEBAYES_name).index(SoftwareNames.SOFTWARE_FREEBAYES_name.lower()) != -1)
		self.assertFalse(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		self.assertFalse(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		expected_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "run_freebayes_single_file1.tab")
		self.assertFalse(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, SoftwareNames.SOFTWARE_FREEBAYES_name)))

		## test consensus file
		self.assertTrue(os.path.exists(project_sample.get_consensus_file(TypePath.MEDIA_ROOT)))
		self.assertTrue(os.path.getsize(project_sample.get_consensus_file(TypePath.MEDIA_ROOT)) > 10)
		self.assertTrue(filecmp.cmp(project_sample.get_consensus_file(TypePath.MEDIA_ROOT),
							project_sample.get_backup_consensus_file()))

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
		self.assertEquals(0, project_sample.alert_second_level)
		
		### test mixed infections
		self.assertEquals(Constants.EMPTY_VALUE_NA, project_sample.mixed_infections.tag.name)
		
		list_meta = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Snippy_Freebayes, None)
		self.assertEquals(1, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Snippy_Freebayes, list_meta[0].meta_tag.name)
		
		decode_result = DecodeObjects()
		result = decode_result.decode_result(list_meta[0].description)
		self.assertTrue(result is not None)
		self.assertEquals("Snippy-3.2-dev; (--mapqual 20 --mincov 5 --minfrac 0.51)", result.get_software(self.software_names.get_snippy_name()))
		self.assertEquals("MSA Masker-1.0; (--c 4; for coverages less than 5 in 30% of the regions.)", result.get_software(self.software_names.get_msa_masker_name()))
		self.assertEquals("", result.get_software(self.software_names.get_freebayes_name()))
		self.assertTrue(result.is_software_present(self.software_names.get_snippy_name()))
		self.assertFalse(result.is_software_present(self.software_names.get_freebayes_name()))
		self.assertFalse(result.is_software_present(self.software_names.get_trimmomatic_name()))
		
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

		meta_value = manageDatabase.get_project_sample_metakey_last(project_sample,
						MetaKeyAndValue.META_KEY_bam_stats,
						MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(not meta_value is None)
		decode_masking_vcf = DecodeObjects()
		result = decode_masking_vcf.decode_result(meta_value.description)
		self.assertEqual('87818', result.get_value_by_key(MetaKeyAndValue.SAMTOOLS_flagstat_mapped_reads))
		self.assertEqual('87818', result.get_value_by_key(MetaKeyAndValue.SAMTOOLS_flagstat_total_reads))
		self.assertEqual(None, result.get_value_by_key('xpto'))
		
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
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina)
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
		self.assertEquals("Freebayes-v1.1.0-54-g49413aa; (--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10 --min-alternate-fraction 0.01 --ploidy 2 -V)",\
 						result.get_software(self.software_names.get_freebayes_name()))

		meta_value = manageDatabase.get_project_sample_metakey_last(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
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

		###############################################
		###################################
		### run again
		default_software = DefaultSoftware()
		default_software.test_all_defaults(user)
		try:
			software = SoftwareSettings.objects.get(name=SoftwareNames.SOFTWARE_SNIPPY_name, owner=user,\
							type_of_use = SoftwareSettings.TYPE_OF_USE_global)
			self.assertTrue(software.is_used_in_global())
		except SoftwareSettings.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(3, len(parameters))
		
		### test set default
		self.assertEqual("20", parameters[0].parameter)
		parameter = parameters[0]
		parameter.parameter = "15"
		parameter.save()
		self.assertEqual("10", parameters[1].parameter)
		parameter = parameters[1]
		parameter.parameter = "5"
		parameter.save()
		self.assertEqual("0.51", parameters[2].parameter)
		parameter = parameters[2]
		parameter.parameter = "0.40"
		parameter.save()
		self.assertEqual("--mapqual 15 --mincov 5 --minfrac 0.40", default_software.get_snippy_parameters(user))

		try:
			project_sample = ProjectSample.objects.get(pk=project_sample.id)
		except ProjectSample.DoesNotExist:
			self.fail("Must exist")
		project_sample.is_finished = False
		project_sample.save()
		
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
		self.assertTrue(len(list_meta) == 2)
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
		self.assertEquals(2, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Snippy_Freebayes, list_meta[0].meta_tag.name)
		
		decode_result = DecodeObjects()
		result = decode_result.decode_result(list_meta[0].description)
		self.assertTrue(result is not None)
		self.assertEquals("Snippy-3.2-dev; (--mapqual 20 --mincov 10 --minfrac 0.51)", result.get_software(self.software_names.get_snippy_name()))
		self.assertEquals("Freebayes-v1.1.0-54-g49413aa; (--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10 --min-alternate-fraction 0.01 --ploidy 2 -V)",\
 						result.get_software(self.software_names.get_freebayes_name()))

		meta_value = manageDatabase.get_project_sample_metakey_last(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
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
		self.assertEquals(4, len(lst_meta_sample))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Queue, lst_meta_sample[3].value)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, lst_meta_sample[2].value)
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

		#### add the transform p.Val423Glu to p.V423G
		### this run a on SNP_EFF
		parse_out_files = ParseOutFiles()
		vcf_file_with_pp = parse_out_files.add_amino_single_letter_code(vcf_file)
		
		out_file = self.utils.get_temp_file("snippy_vcf_to_tab", ".tab")
		out_file_2 = self.software.run_snippy_vcf_to_tab(fasta_file, gb_file, vcf_file_with_pp, out_file)
		self.assertEquals(out_file, out_file_2)
		self.assertTrue(filecmp.cmp(out_file_2, result_file))
		os.unlink(out_file)
		os.unlink(vcf_file_with_pp)


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


	def test_run_genbank2gff3_positions_1(self):
		"""
		test genbank2gff3 method
		"""
		## covid
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_COVID_GBK)
		gff_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "covid.gff3")
		self.assertTrue(os.path.exists(gb_file))
		out_file = self.utils.get_temp_file("file_name", ".txt")
		out_file_2 = self.software.run_genbank2gff3(gb_file, out_file)
		self.assertFalse(out_file_2 is None)
		self.assertTrue(filecmp.cmp(out_file_2, gff_file))
		os.unlink(out_file)
		
	def test_run_genbank2gff3(self):
		"""
		test genbank2gff3 method
		"""
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		gff_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GFF)
		self.assertTrue(os.path.exists(gb_file))
		out_file = self.utils.get_temp_file("file_name", ".txt")
		out_file_2 = self.software.run_genbank2gff3(gb_file, out_file)
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
		self.assertTrue(filecmp.cmp(out_file_2, gff_file))
		os.unlink(out_file)

	def test_run_genbank2gff3_for_nextclade(self):
		"""
		test genbank2gff3 method
		"""
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_COVID_GBK)
		gff_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "covid_for_nextclade.gff3")
		self.assertTrue(os.path.exists(gb_file))
		out_file = self.utils.get_temp_file("file_name_rim", ".txt")
		for_nextclade = True
		out_file_2 = self.software.run_genbank2gff3(gb_file, out_file, for_nextclade)
		self.assertTrue(filecmp.cmp(out_file_2, gff_file))
		os.unlink(out_file)

	def test_run_genbank2gff3_for_nextclade_2(self):
		"""
		test genbank2gff3 method
		"""
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		gff_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "covid_for_nextclade_2.gff3")
		self.assertTrue(os.path.exists(gb_file))
		out_file = self.utils.get_temp_file("file_name", ".txt")
		run_for_nextstrain = True
		out_file_2 = self.software.run_genbank2gff3(gb_file, out_file, run_for_nextstrain)
		self.assertTrue(filecmp.cmp(out_file_2, gff_file))
		os.unlink(out_file)
		
	def test_run_genbank2gff3_positions_comulative(self):
		"""
		test genbank2gff3 method
		"""
		## covid
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_COVID_GBK)
		gff_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "covid_comulative_gene_position.gff")
		self.assertTrue(os.path.exists(gb_file))
		out_file = self.utils.get_temp_file("file_name", ".txt")
		out_file_2 = self.software.run_genbank2gff3_positions_comulative(gb_file, out_file)
		self.assertFalse(out_file_2 is None)
		self.assertEquals(out_file, out_file_2)
		self.assertTrue(filecmp.cmp(out_file_2, gff_file))
		
		### flu
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		gff_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "flu_comulative_gene_position.gff")
		self.assertTrue(os.path.exists(gb_file))
		out_file = self.utils.get_temp_file("file_name", ".txt")
		out_file_2 = self.software.run_genbank2gff3_positions_comulative(gb_file, out_file)
		self.assertFalse(out_file_2 is None)
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
		
		out_file = self.utils.get_temp_file("file_name", ".vcf")
		out_file_2 = self.software.run_snpEff(fasta_file, genbank_file, freebayes_vcf, out_file)
		self.assertEquals(out_file, out_file_2)
		
		out_file_clean = self.utils.get_temp_file("file_name", ".vcf")
		cmd = "grep -v '{}' {} > {}".format(os.path.dirname(out_file), out_file, out_file_clean)
		os.system(cmd)
		self.assertTrue(filecmp.cmp(out_file_clean, freebayes_expect_vcf))
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
		
		### reroot the tree with fake leaf name		
		reroot_leaf = "EVA011_S54_"
		output_file = self.software.run_fasttree(in_file, out_file, self.software_names.get_fasttree_parameters(), reroot_leaf)
		self.assertTrue(filecmp.cmp(out_file, expect_file_nwk))
		
		### reroot the tree		
		reroot_leaf = "EVA011_S54"
		output_file = self.software.run_fasttree(in_file, out_file, self.software_names.get_fasttree_parameters(), reroot_leaf)
		expect_file_nwk = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "2_" + ConstantsTestsCase.FILE_FASTTREE_RESULT_NWK)
		self.assertTrue(filecmp.cmp(out_file, expect_file_nwk))
		os.unlink(out_file)


	def test_run_nextalign(self):
		
		reference_fasta = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		sequence_fasta = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "run_snippy1_sdfs.consensus.fa")
		
		temp_dir = self.utils.get_temp_dir()
		gff_file = self.utils.get_temp_file_from_dir(temp_dir, "file_name", ".gff3")
		self.software.run_genbank2gff3(gb_file, gff_file, True)
		genes_to_process = "PA"
		
		ref_file = self.utils.get_temp_file_from_dir(temp_dir, "reference", ".fasta")
		records = []
		with open(reference_fasta) as handle_ref: 
			record_dict_ref = SeqIO.to_dict(SeqIO.parse(handle_ref, "fasta"))

			if genes_to_process in record_dict_ref:
				records.append(record_dict_ref[genes_to_process])
				records[-1].id = "reference"
				### save file
				with open(ref_file, 'w') as handle_write:
					SeqIO.write(records, handle_write, "fasta")
			self.assertTrue(len(records) > 0)
					
		seq_file = self.utils.get_temp_file_from_dir(temp_dir, "sequence", ".fasta")
		records = []
		with open(sequence_fasta) as handle_ref: 
			record_dict_ref = SeqIO.to_dict(SeqIO.parse(handle_ref, "fasta"))

			if genes_to_process in record_dict_ref:
				records.append(record_dict_ref[genes_to_process])
				### save file
				with open(seq_file, 'w') as handle_write:
					SeqIO.write(records, handle_write, "fasta")
			self.assertTrue(len(records) > 0)
		
		
		alignment_file, vect_protein_file, insert_positions_file = \
			self.software.run_nextalign(ref_file, seq_file, gff_file, genes_to_process, temp_dir)
		self.assertEqual(1, len(vect_protein_file))
		expect_result = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "nextalign.aligned.fasta")
		self.assertTrue(filecmp.cmp(alignment_file, expect_result))
		expect_result = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "nextalign.PA.fasta")
		self.assertTrue(filecmp.cmp(vect_protein_file[0], expect_result))
		expect_result = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "nextalign.insertions.csv")
		self.assertTrue(filecmp.cmp(insert_positions_file, expect_result))
		self.utils.remove_dir(temp_dir)


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
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina)
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
		
		### files illumina covid
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, "covid/WHU02_SRR10903401_1.fastq.gz")
		file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, "covid/WHU02_SRR10903401_2.fastq.gz")
		self.assertTrue(os.path.exists(file_1))
		self.assertTrue(os.path.exists(file_2))
		
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
			sample.file_name_1 = os.path.basename(file_1)
			sample.path_name_1.name = file_1
			sample.is_valid_2 = True
			sample.file_name_2 = os.path.basename(file_2)
			sample.path_name_2.name = file_2
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina)
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
		
		self.assertEquals(os.path.getsize(fastq1_1), 5314574)
		self.assertEquals(os.path.getsize(fastq1_2), 5715245)
		(is_downsized, file_name_1, file_name_2) = self.software.make_downsize(fastq1_1, fastq1_2, 2000000)	## downsize to 2M
		
		self.assertEquals(is_downsized, True)
		self.assertTrue(os.path.exists(file_name_1))
		self.assertTrue(os.path.exists(file_name_2))
		
		self.assertEquals(os.path.getsize(file_name_1), 1724181)
		self.assertEquals(os.path.getsize(file_name_2), 1924367)
		self.utils.remove_dir(os.path.dirname(file_name_1))

		(is_downsized, file_name_1, file_name_2) = self.software.make_downsize(fastq1_1, None, 2000000)	## downsize to 2M
		
		self.assertEquals(is_downsized, True)
		self.assertTrue(os.path.exists(file_name_1))
		self.assertTrue(file_name_2 == None)
		
		self.assertEquals(os.path.getsize(file_name_1), 1839658)
		self.utils.remove_dir(os.path.dirname(file_name_1))
		
		(is_downsized, file_name_1, file_name_2) = self.software.make_downsize(fastq1_1, fastq1_2, 20000000)	## downsize to 2M
		self.assertEquals(is_downsized, False)


	def test_make_mask_consensus_by_deep(self):
		"""
		run make mask consensus
		"""

		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
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

		ref_name = "second_stage_2"
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
			
		project_name = "file_name_3_dsf"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = sample_name
			project.reference = reference
			project.owner = user
			project.save()
		
		project_sample = ProjectSample()
		project_sample.sample = sample
		project_sample.project = project
		project_sample.save()
		
		software_names = SoftwareNames()
		out_put_path = self.software.run_snippy(sample.get_fastq(TypePath.MEDIA_ROOT, True), sample.get_fastq(TypePath.MEDIA_ROOT, False),
				fasta_file, gb_file, sample.name, software_names.get_snippy_parameters())

		### set coverage statistics
		get_coverage = GetCoverage()
		coverage = get_coverage.get_coverage(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_DEPTH_GZ, ""))),\
						project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT),
						700)

		self.assertFalse(coverage.is_100_more_9("NP"))
		self.assertEqual(0.0, coverage.ratio_value_defined_by_user)
		self.assertFalse(coverage.is_100_more_9("PA"))
		self.assertEqual(30.0, coverage.ratio_value_defined_by_user)
		self.assertFalse(coverage.is_100_more_9("NA"))
		self.assertEqual(86.0, coverage.ratio_value_defined_by_user)
		self.assertFalse(coverage.is_100_more_9("PB2"))
		self.assertEqual(49.0, coverage.ratio_value_defined_by_user)	

		limit_make_mask = 70
		msa_parameters = self.software.make_mask_consensus_by_deep(os.path.join(out_put_path,
			os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_CONSENSUS_FA, ""))),
 			project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT),
 			os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_DEPTH_GZ, ""))), 
 			coverage, sample_name, limit_make_mask)
		self.assertEqual("--c 699", msa_parameters)

		self.assertTrue(filecmp.cmp(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_CONSENSUS_FA, ""))),
					os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, 'run_snippy1_sdfs.consensus.fa' )) )
		
		remove_path = os.path.dirname(out_put_path)
		if (len(remove_path.split('/')) > 2): self.utils.remove_dir(remove_path)
		else: self.utils.remove_dir(out_put_path)
		self.utils.remove_dir(temp_dir)

	def test_run_snippy_vcf_to_tab_freq_and_evidence(self):
		
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		vcf_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "run_snippy_vcf_to_tab_freq_and_evidence.vcf")
		self.assertTrue(os.path.exists(fasta_file))
		self.assertTrue(os.path.exists(gb_file))
		self.assertTrue(os.path.exists(vcf_file))
		
		out_file_vcf = self.utils.get_temp_file("vcf_file", ".vcf")
		out_file_2 = self.utils.get_temp_file("parse_vcf_file", ".tab")
		tab_expect_result = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "run_snippy_vcf_to_tab_freq_and_evidence.tab")
		
		#### add the transform p.Val423Glu to p.V423G
		### this run a on SNP_EFF
		parse_out_files = ParseOutFiles()
		vcf_file_with_pp = parse_out_files.add_amino_single_letter_code(vcf_file)
		
		### first need to add FREQ
		out_file_vcf = self.utils.add_freq_to_vcf(vcf_file_with_pp, out_file_vcf)
		self.assertTrue(os.path.exists(out_file_vcf))
		
		self.software.run_snippy_vcf_to_tab_freq_and_evidence(fasta_file, gb_file, out_file_vcf, out_file_2)
		self.assertTrue(os.path.exists(out_file_2))
		self.assertTrue(os.path.exists(tab_expect_result))
		
		self.assertTrue(filecmp.cmp(tab_expect_result, out_file_2))
		if (os.path.exists(out_file_vcf)): os.unlink(out_file_vcf)
		if (os.path.exists(out_file_2)): os.unlink(out_file_2)
		if (os.path.exists(vcf_file_with_pp)): os.unlink(vcf_file_with_pp)


	def test_run_pangolin(self):
		"""
		test genbank2gff3 method
		"""
		## covid
		
		consensus_file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "SARSCoVDec200153.consensus.fasta")
		consensus_file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "SARSCoVDec200234.consensus.fasta")
		out_file_consensus = self.utils.get_temp_file("all_file_name", ".fasta")
		cmd = "cat {} {} > {}".format(consensus_file_1, consensus_file_2, out_file_consensus)
		os.system(cmd)
		
		#print("Testing " + SoftwareNames.SOFTWARE_Pangolin_name)
		#try:
		#	software = SoftwareModel.objects.get(name=SoftwareNames.SOFTWARE_Pangolin_name)
		#	self.fail("Must not exist software name")
		#except SoftwareModel.DoesNotExist:	## need to create with last version
		#	pass

		self.software_pangolin.run_pangolin_update()

		#try:
		#	software = SoftwareModel.objects.get(name=SoftwareNames.SOFTWARE_Pangolin_name)
		#	software.is_updated_today()
		#	dt_software = software.get_version_long()
		#	self.assertTrue(len(dt_software) > 0)
		#	self.assertTrue(len(software.version) > 0)
		#except SoftwareModel.DoesNotExist:	## need to create with last version
		#	self.fail("Must not exist software name")
		
		out_file = self.utils.get_temp_file("file_name", ".txt")
		self.software_pangolin.run_pangolin(out_file_consensus, out_file)
		
		vect_data = self.utils.read_text_file(out_file)
		self.assertEqual(3, len(vect_data))

		#try:
		#	software = SoftwareModel.objects.get(name=SoftwareNames.SOFTWARE_Pangolin_name)
		#	dt_versions = software.get_version_long()
		#	self.assertTrue(len(dt_versions) > 0)
		#	self.assertTrue(len(software.version) > 0)							
		#except SoftwareModel.DoesNotExist:	## need to create with last version
		#	self.fail("Must exist software name")
		
		os.unlink(out_file)
		os.unlink(out_file_consensus)

	def test_is_ref_sars_cov(self):
		"""
		test SARS cov 
		"""
		uploadFiles = UploadFiles()
		version = 1
		file = os.path.join(self.baseDirectory, Constants.DIR_TYPE_IDENTIFICATION, "test_covid_typing.fasta")
		uploadFiles.upload_file(version, file)	## upload file
		
		database_name = "xpto_sars_cov"
		if (not self.software.is_exist_database_abricate(database_name)):
			self.software.create_database_abricate(database_name, file)
			
		consensus_file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "SARSCoVDec200153.consensus.fasta")
		self.assertTrue(self.software_pangolin.is_ref_sars_cov_2(consensus_file_1))
		
		consensus_file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "A_H3N2_A_Hong_Kong_4801_2014.fasta")
		self.assertFalse(self.software_pangolin.is_ref_sars_cov_2(consensus_file_1))
	
	def test_align_two_sequences(self):
		
		seq_ref, seq_consensus = self.software.align_two_sequences(
			"ATGGAAGATTTTGTGCGACAATGCTTCAACCCGATGATTGTCGAACTTGCAGAAAAAGCAATGAAAGAGTATGGGGAGGATCGAAAATTGAAACCAACAAATTT",
			"ATGGAAGATTTTGTGCGACAATGCTTCAACCGATGATTGTCGAACTTGCAGACGAAAAGCAATGAAAGAGTATGGGGAGGATCTGAAAATTGAAACCAACAAATTT")
		self.assertEqual("ATGGAAGATTTTGTGCGACAATGCTTCAACCCGATGATTGTCGAACTTGCAGA--AAAAGCAATGAAAGAGTATGGGGAGGATC-GAAAATTGAAACCAACAAATTT",
						seq_ref)
		self.assertEqual("ATGGAAGATTTTGTGCGACAATGCTTCAA-CCGATGATTGTCGAACTTGCAGACGAAAAGCAATGAAAGAGTATGGGGAGGATCTGAAAATTGAAACCAACAAATTT",
						seq_consensus)


	def test_mask_sequences_by_position(self):
		
		seq_ref =       SeqRecord(Seq("AACA-AC--AAA"), id="xpto")
		seq_consensus = SeqRecord(Seq("AA-AAACAC--C"), id="xpto")
		mask_sites = "5,6"
		mask_from_beginning = "20"
		mask_from_end = "20"
		mask_range = "1-20"
		sequence_consensus = self.software.mask_sequence(seq_ref, seq_consensus, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual(12, len(str(sequence_consensus.seq)))
		self.assertEqual("AAANNNNNNNNN", str(sequence_consensus.seq))
		
		seq_ref =       SeqRecord(Seq("AACA-AC--AAA"), id="xpto")
		seq_consensus = SeqRecord(Seq("AA-AAACAC--C"), id="xpto")
		mask_sites = "1,2"
		mask_from_beginning = "-1"
		mask_from_end = "-1"
		mask_range = ""
		sequence_consensus = self.software.mask_sequence(seq_ref, seq_consensus, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual(9, len(str(sequence_consensus.seq)))
		self.assertEqual("AAANNCACC", str(sequence_consensus.seq))

		seq_ref =       SeqRecord(Seq("AACAACAAA"), id="xpto")
		seq_consensus = SeqRecord(Seq("AAAAACACC"), id="xpto")
		mask_sites = ""
		mask_from_beginning = "2"
		mask_from_end = "2"
		mask_range = ""
		sequence_consensus = self.software.mask_sequence(seq_ref, seq_consensus, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual(11, len(str(sequence_consensus.seq)))
		self.assertEqual("AAANNCACCNN", str(sequence_consensus.seq))
		
		seq_ref =       SeqRecord(Seq("AACAACAAA"), id="xpto")
		seq_consensus = SeqRecord(Seq("AACAACAAA"), id="xpto")
		mask_sites = ""
		mask_from_beginning = "2"
		mask_from_end = "2"
		mask_range = ""
		sequence_consensus = self.software.mask_sequence(seq_ref, seq_consensus, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual(9, len(str(sequence_consensus.seq)))
		self.assertEqual("NNCAACANN", str(sequence_consensus.seq))
			
		seq_ref =       SeqRecord(Seq("AACAACAAA"), id="xpto")
		seq_consensus = SeqRecord(Seq("AACAACAAA"), id="xpto")
		mask_sites = ""
		mask_from_beginning = "2"
		mask_from_end = "2"
		mask_range = "0-40"
		sequence_consensus = self.software.mask_sequence(seq_ref, seq_consensus, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual(9, len(str(sequence_consensus.seq)))
		self.assertEqual("NNNNNNNNN", str(sequence_consensus.seq))
		
		seq_ref =       SeqRecord(Seq("AACAACCAAA"), id="xpto")
		seq_consensus = SeqRecord(Seq("AACAACAAA"), id="xpto")
		mask_sites = "5,6"
		mask_from_beginning = "0"
		mask_from_end = "0"
		mask_range = "0"
		sequence_consensus = self.software.mask_sequence(seq_ref, seq_consensus, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual(9, len(str(sequence_consensus.seq)))
		self.assertEqual("AACANNAAA", str(sequence_consensus.seq))
		
		seq_ref =       SeqRecord(Seq("AACAAAAA"), id="xpto")
		seq_consensus = SeqRecord(Seq("AACAACAAA"), id="xpto")
		mask_sites = "5,6"
		mask_from_beginning = "0"
		mask_from_end = "0"
		mask_range = "0"
		sequence_consensus = self.software.mask_sequence(seq_ref, seq_consensus, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual(9, len(str(sequence_consensus.seq)))
		self.assertEqual("AACANNAAA", str(sequence_consensus.seq))
		
		seq_ref =       SeqRecord(Seq("AACAAAAAAAAAAAAAAA"), id="xpto")
		seq_consensus = SeqRecord(Seq("AACAACAAAAAAAAAAAAA"), id="xpto")
		mask_sites = "5,6"
		mask_from_beginning = "0"
		mask_from_end = "0"
		mask_range = "0"
		sequence_consensus = self.software.mask_sequence(seq_ref, seq_consensus, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual(19, len(str(sequence_consensus.seq)))
		self.assertEqual("AACANNAAAAAAAAAAAAA", str(sequence_consensus.seq))

		temp_ref_file = self.utils.get_temp_file("ref_test", ".fasta")
		temp_consensus_file = self.utils.get_temp_file("consensus_test", ".fasta")
		vect_record_out = [SeqRecord(Seq("AACAAAAAAAAAAAAAAA"), id="xpto")]
		with open(temp_ref_file, "w") as handle_fasta_out:
			SeqIO.write(vect_record_out, handle_fasta_out, "fasta")
		vect_record_out = [SeqRecord(Seq("AACAACAAAAAAAAAAAAA"), id="xpto")]
		with open(temp_consensus_file, "w") as handle_fasta_out:
			SeqIO.write(vect_record_out, handle_fasta_out, "fasta")
		
		masking_consensus = MaskingConsensus()
		masking_consensus.set_mask_sites("7")
		masking_consensus.set_mask_from_beginning("2")
		masking_consensus.set_mask_from_ends("2")
		masking_consensus.set_mask_regions("10-11")
		genetic_element = GeneticElement()
		genetic_element.set_mask_consensus_element("xpto", masking_consensus)
		self.software.mask_sequence_by_sites(temp_ref_file, temp_consensus_file, genetic_element)
		dt_records = SeqIO.to_dict(SeqIO.parse(temp_consensus_file, "fasta"))
		self.assertEqual("NNCAACNAANNAAAAANNA", str(dt_records["xpto"].seq))
		
		vect_record_out = [SeqRecord(Seq("AACAACAAAAAAAAAAAA"), id="xpto")]
		with open(temp_consensus_file, "w") as handle_fasta_out:
			SeqIO.write(vect_record_out, handle_fasta_out, "fasta")
		masking_consensus = MaskingConsensus()
		masking_consensus.set_mask_sites("7")
		masking_consensus.set_mask_from_beginning("2")
		masking_consensus.set_mask_from_ends("3")
		masking_consensus.set_mask_regions("10-14")
		genetic_element = GeneticElement()
		genetic_element.set_mask_consensus_element("xpto", masking_consensus)
		self.software.mask_sequence_by_sites(temp_ref_file, temp_consensus_file, genetic_element)
		dt_records = SeqIO.to_dict(SeqIO.parse(temp_consensus_file, "fasta"))
		self.assertEqual("NNCAACNAANNNNNANNN", str(dt_records["xpto"].seq))
		
		os.unlink(temp_consensus_file)
		os.unlink(temp_ref_file)


	def test_get_statistics_bam(self):
		""" test stasts in bam file """
		bam_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_BAM, "run_snippy1.bam")
		self.assertTrue(os.path.exists(bam_file))
		result = self.software.get_statistics_bam(bam_file)
		self.assertEqual('43909', result.get_value_by_key(MetaKeyAndValue.SAMTOOLS_flagstat_mapped_reads))
		self.assertEqual('43909', result.get_value_by_key(MetaKeyAndValue.SAMTOOLS_flagstat_total_reads))
		self.assertEqual(None, result.get_value_by_key('xpto'))


	def test_run_nextstrain_ncov(self):
		""" test running nexstrain ncov build"""

		alignments_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "sequences_ncov.fasta")
		metadata_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "metadata_ncov.tsv")

		[tree_file, alignments_file, zip_file] = self.software.run_nextstrain_ncov(alignments=alignments_file, metadata=metadata_file)

		self.assertEqual(os.path.exists(zip_file),True)
		self.assertEqual(os.path.exists(alignments_file),True)
		self.assertEqual(os.path.exists(tree_file),True)

		self.utils.remove_file(zip_file)
		self.utils.remove_file(alignments_file)
		self.utils.remove_file(tree_file)

		# TODO look inside the results to see if we get what we expect to get
		#self.assertEqual(os.path.exists(temp_dir + "/auspice"),True)
		#self.assertEqual(os.path.exists(temp_dir + "/auspice/ncov_current.json"),True)
		#self.assertEqual(os.path.exists(temp_dir + "/auspice/ncov_current_root-sequence.json"),True)
		#self.assertEqual(os.path.exists(temp_dir + "/auspice/ncov_current_tip-frequencies.json"),True)




	def test_run_nextstrain_generic(self):
		""" test running nexstrain generic (default) build"""
		
		alignments_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "sequences_generic.fasta")
		metadata_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "metadata_generic.tsv")

		reference_fasta = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "reference_generic.fasta")
		reference_gb = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "reference_generic.gb")

		[tree_file, alignments_file, zip_file] = self.software.run_nextstrain_generic(alignments=alignments_file, metadata=metadata_file, 
														 ref_fasta=reference_fasta, ref_genbank=reference_gb)

		self.assertEqual(os.path.exists(zip_file),True)
		self.assertEqual(os.path.exists(alignments_file),True)
		self.assertEqual(os.path.exists(tree_file),True)

		self.utils.remove_file(zip_file)
		self.utils.remove_file(alignments_file)
		self.utils.remove_file(tree_file)

		# TODO look inside the results to see if we get what we expect to get
		#self.assertEqual(os.path.exists(temp_dir + "/auspice"),True)
		#self.assertEqual(os.path.exists(temp_dir + "/auspice/generic.json"),True)


	def test_run_nextstrain_mpx(self):
		""" test running nexstrain mpx build"""
		
		alignments_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "sequences_mpx.fasta")
		metadata_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "metadata_mpx.tsv")

		[tree_file, alignments_file, zip_file] = self.software.run_nextstrain_mpx(alignments=alignments_file, metadata=metadata_file)

		self.assertEqual(os.path.exists(zip_file),True)
		self.assertEqual(os.path.exists(alignments_file),True)
		self.assertEqual(os.path.exists(tree_file),True)

		self.utils.remove_file(zip_file)
		self.utils.remove_file(alignments_file)
		self.utils.remove_file(tree_file)

		# TODO look inside the results to see if we get what we expect to get, eg. test if they are empty!!
		#self.assertEqual(os.path.exists(temp_dir + "/auspice"),True)
		#self.assertEqual(os.path.exists(temp_dir + "/auspice/monkeypox.json"),True)
		#self.assertEqual(os.path.exists(temp_dir + "/auspice/monkeypox_root-sequence.json"),True)

		#self.utils.remove_dir(temp_dir)	


	def test_run_nextstrain_flu(self):
		""" test running nexstrain flu build"""

		alignments_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "sequences_h3n2_ha.fasta")
		metadata_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "metadata_h3n2_ha.tsv")

		# TODO make tests for other strains and periods...
		[tree_file, alignments_file, zip_file] = self.software.run_nextstrain_flu(alignments=alignments_file, metadata=metadata_file, strain="h3n2", period="12y")

		self.assertEqual(os.path.exists(zip_file),True)
		self.assertEqual(os.path.exists(alignments_file),True)
		self.assertEqual(os.path.exists(tree_file),True)

		self.utils.remove_file(zip_file)
		self.utils.remove_file(alignments_file)
		self.utils.remove_file(tree_file)

		# TODO look inside the results to see if we get what we expect to get
		#self.assertEqual(os.path.exists(temp_dir + "/auspice"),True)
		#self.assertEqual(os.path.exists(temp_dir + "/auspice/flu_h3n2_ha_12y.json"),True)
		#self.assertEqual(os.path.exists(temp_dir + "/auspice/flu_h3n2_ha_12y_root-sequence.json"),True)
		#self.assertEqual(os.path.exists(temp_dir + "/auspice/flu_h3n2_ha_12y_tip-frequencies.json"),True)

		#self.utils.remove_dir(temp_dir)	



	def test_run_aln2pheno(self):
		""" test running aln2pheno """
		alignments_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "test_Alignment_aa_SARS_CoV_2_S.fasta")

		temp_dir = os.path.join(self.utils.get_temp_dir())
		report = os.path.join(temp_dir, '/final_report.tsv')
		flagged = os.path.join(temp_dir, '/flagged_mutation_report.tsv')

		exit_status = self.software.run_aln2pheno(sequences=alignments_file,  reference="SARS_CoV_2_Wuhan_Hu_1_MN908947_SARS_CoV_2_S", gene="S", report=report, flagged=flagged)

		self.assertEqual(os.path.exists(report),True)
		self.assertEqual(os.path.exists(flagged),True)

		report_data = pandas.read_csv(report, delimiter=Constants.SEPARATOR_TAB)

		# create a function in utils, or within DataColumns or something (not sure if this is sufficiently generic...)
		pangolin_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "test_Alignment_aa_SARS_CoV_2_S_pangolin_lineage.csv")
		pangolin_data = pandas.read_csv(pangolin_file, delimiter=Constants.SEPARATOR_COMMA)
		pangolin_data = pangolin_data[['taxon','lineage']]
		pangolin_data['taxon'] = pangolin_data['taxon'].str.replace('__SARS_CoV_2','')
		pangolin_data.rename(columns = {'taxon':'Sequence'}, inplace = True)

		report_data = report_data.merge(pangolin_data, on=["Sequence"])
		self.assertEqual(report_data['lineage'][0], "BA.2")

		final_report = os.path.join(temp_dir, '/flagged_mutation_report_lineage.tsv')
		report_data.to_csv(final_report, sep=Constants.SEPARATOR_TAB, index=False)
		self.assertEqual(os.path.exists(final_report),True)

		# TODO collect results and check content
		self.utils.remove_dir(temp_dir)

		

