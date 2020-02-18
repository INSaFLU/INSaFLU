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
from constants.tag_names_constants import TagNamesConstants
from utils.software import Software, Contigs2Sequences
from constants.software_names import SoftwareNames
from utils.utils import Utils
from utils.parse_out_files import ParseOutFiles
from utils.result import DecodeObjects, Coverage, CountHits
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference, Statistics
from manage_virus.uploadFiles import UploadFiles
from manage_virus.models import UploadFile
from managing_files.manage_database import ManageDatabase
from django.test.utils import override_settings
from utils.tree import CreateTree
import os, filecmp, csv
from utils.parse_in_files import ParseInFiles
from utils.result import DecodeObjects, MixedInfectionMainVector
from managing_files.models import CountVariations, MixedInfections
from utils.mixed_infections_management import MixedInfectionsManagement
from manage_virus.constants_virus import ConstantsVirus
from manage_virus.models import Tags, SeqVirus, IdentifyVirus

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
			sample.file_name_1 = "Ct-LGV-Gc-6477-G2A-y2019_S12_L001_R2_001.fastq.gz"
			sample.path_name_1.name = os.path.join("/home/mmp/insa", "Ct-LGV-Gc-6477-G2A-y2019_S12_L001_R2_001.fastq.gz")
			sample.is_valid_2 = False
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