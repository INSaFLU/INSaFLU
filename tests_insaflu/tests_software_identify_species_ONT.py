'''
Created on Oct 28, 2017

@author: mmp
'''

import os
from django.test import TestCase
from django.conf import settings 
from constants.constantsTestsCase import ConstantsTestsCase
from constants.constants import Constants, TypePath
from constants.meta_key_and_values import MetaKeyAndValue
from utils.software import Software, Contigs2Sequences
from constants.software_names import SoftwareNames
from utils.utils import Utils
from django.contrib.auth.models import User
from managing_files.models import Sample
from manage_virus.uploadFiles import UploadFiles
from manage_virus.models import UploadFile
from managing_files.manage_database import ManageDatabase
from manage_virus.constants_virus import ConstantsVirus

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

	
	def test_identify_type_and_sub_type(self):
		"""
 		get type and sub_type
 		
 		need to check line 1802 software.py
 		## NEED to check coverage for CANU, it is checked for Spades
 		## Need to define genome size in CANU, need to check this with Vitor 
 		"""
		
		uploadFiles = UploadFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		self.assertEqual("2", version)
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
		uploadFiles.upload_file(version, file)	## upload file
		
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, "covid/test_minion_seq.fastq.gz")
		self.assertTrue(os.path.exists(file_1))
		
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
			sample.is_valid_2 = False
			sample.file_name_2 = ""
			sample.path_name_2.name = ""
			sample.owner = user
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_ont)
			sample.save()
			
		return_value = self.software.identify_type_and_sub_type(sample, sample.path_name_1.name, None, user, True)
		self.assertTrue(return_value)
		
		vect_identify_virus = sample.identify_virus.all()
		self.assertEqual(3, len(vect_identify_virus))
		for identify_virus in vect_identify_virus:
			if (identify_virus.rank == 0):
				self.assertEquals("99.19", identify_virus.coverage)
				self.assertEquals("85.71", identify_virus.identity)
				self.assertEquals("A", identify_virus.seq_virus.name)
				self.assertEquals(ConstantsVirus.SEQ_VIRUS_TYPE, identify_virus.seq_virus.kind_type.name)
			elif (identify_virus.rank == 1):
				self.assertEquals("98.80", identify_virus.coverage)
				self.assertEquals("85.83", identify_virus.identity)
				self.assertEquals("N5", identify_virus.seq_virus.name)
				self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, identify_virus.seq_virus.kind_type.name)
			elif (identify_virus.rank == 2):
				self.assertEquals("98.59", identify_virus.coverage)
				self.assertEquals("87.10", identify_virus.identity)
				self.assertEquals("H5", identify_virus.seq_virus.name)
				self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, identify_virus.seq_virus.kind_type.name)
		
		### full descriptiom		
		self.assertEquals("A-H5N5", sample.get_type_sub_type())
		
		file_abricate = sample.get_abricate_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(sample.get_abricate_output(TypePath.MEDIA_ROOT)))
		
		file_spades_contigs = sample.get_draft_contigs_output(TypePath.MEDIA_ROOT)
		self.assertFalse(os.path.exists(file_spades_contigs))
		file_abricate_contigs = sample.get_draft_contigs_abricate_output(TypePath.MEDIA_ROOT)
		self.assertFalse(os.path.exists(file_abricate_contigs))
		file_abricate_contigs = sample.get_draft_reads_abricate_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(file_abricate_contigs))
		
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Identify_Sample, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Identify_Sample, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Abricate(0.8-dev)", list_meta[0].description)
		if (os.path.exists(file_abricate)): os.unlink(file_abricate)
		if (os.path.exists(file_spades_contigs)): os.unlink(file_spades_contigs)
		if (os.path.exists(file_abricate_contigs)): os.unlink(file_abricate_contigs)

		contigs_2_sequences = Contigs2Sequences(True)
		## remove abricate db
		cmd = "rm -r %s/%s*" % (SoftwareNames.SOFTWARE_ABRICATE_DB, contigs_2_sequences.get_database_name())
		exist_status = os.system(cmd)
		self.assertTrue(exist_status == 0)
		
				
		