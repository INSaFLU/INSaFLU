'''
Created on Nov 5, 2017

@author: mmp
'''
import unittest, os 
from django.conf import settings 
from django_q.cluster import Cluster
from manage_virus.uploadFiles import UploadFiles
from manage_virus.models import UploadFile
from django.contrib.auth.models import User
from constants.constantsTestsCase import ConstantsTestsCase
from managing_files.models import Sample
from constants.constants import Constants, TypePath
from manage_virus.constants_virus import ConstantsVirus
from constants.meta_key_and_values import MetaKeyAndValue
from utils.software import Software
from managing_files.manage_database import ManageDatabase

class Test(unittest.TestCase):

	software = Software()
	constants = Constants()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		self.c = Cluster()
		self.c.start()

	def tearDown(self):
		self.c.stop()
	
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
		
		task_id = async_task(self.software.identify_type_and_sub_type, sample, sample.path_name_1.name, sample.path_name_2.name, user, sync=True)
		# the task will then be available immediately
		return_value = fetch(task_id)
		self.assertTrue(return_value)
		
		vect_identify_virus = sample.identify_virus.all()
		self.assertEqual(3, len(vect_identify_virus))
		for identify_virus in vect_identify_virus:
			if (identify_virus.rank == 0):
				self.assertEquals("100.00", identify_virus.coverage)
				self.assertEquals("99.69", identify_virus.identity)
				self.assertEquals("A", identify_virus.seq_virus.name)
				self.assertEquals("XXXX", identify_virus.seq_virus.accession)
				self.assertEquals(ConstantsVirus.SEQ_VIRUS_TYPE, identify_virus.seq_virus.kind_type.name)
			elif (identify_virus.rank == 1):
				self.assertEquals("100.00", identify_virus.coverage)
				self.assertEquals("99.18", identify_virus.identity)
				self.assertEquals("XXXXX", identify_virus.seq_virus.accession)
				self.assertEquals("H3", identify_virus.seq_virus.name)
				self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, identify_virus.seq_virus.kind_type.name)
			elif (identify_virus.rank == 2):
				self.assertEquals("100.00", identify_virus.coverage)
				self.assertEquals("98.65", identify_virus.identity)
				self.assertEquals("N2", identify_virus.seq_virus.name)
				self.assertEquals("XXXXXX", identify_virus.seq_virus.accession)
				self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, identify_virus.seq_virus.kind_type.name)
		file_abricate = sample.get_abricate_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(file_abricate))
		
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Identify_Sample, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Identify_Sample, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Spades(3.11.1), Abricate(0.8-dev)", list_meta[0].description)
		if (os.path.exists(file_abricate)): os.unlink(file_abricate)


