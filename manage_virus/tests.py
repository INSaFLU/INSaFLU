from django.test import TestCase

# Create your tests here.
from utils.constantsTestsCase import ConstantsTestsCase
from django.conf import settings 
from .uploadFiles import UploadFiles
from utils.constants import Constants
from utils.parseOutFiles import ParseOutFiles
from .models import UploadFile, Tags, SeqVirus
import os

class Test(TestCase):

	### static
	constantsTestsCase = ConstantsTestsCase()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), self.constantsTestsCase.MANAGING_TESTS)
		pass


	def tearDown(self):
		pass


	def test_get_file_to_upload(self):
		uploadFiles = UploadFiles()
		
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		self.assertEqual("2", version)
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)

	def test_parse_title_name(self):
		uploadFiles = UploadFiles()

		dt_return = uploadFiles.parse_title_name("influenza_type~~~B~~~AF100378")
		self.assertTrue(Constants.SEQ_VIRUS_TYPE in dt_return)
		self.assertEqual("B", dt_return[Constants.SEQ_VIRUS_TYPE])
		self.assertFalse(Constants.SEQ_VIRUS_SUB_TYPE in dt_return)
		self.assertFalse(UploadFiles.DESCRIPTION in dt_return)
		self.assertEqual("AF100378", dt_return[UploadFiles.ACCESSION])
		
		dt_return = uploadFiles.parse_title_name("influenza_subtype~~~sB~~~AF100378 adadadd")
		self.assertTrue(Constants.SEQ_VIRUS_SUB_TYPE in dt_return)
		self.assertEqual("sB", dt_return[Constants.SEQ_VIRUS_SUB_TYPE])
		self.assertEqual("adadadd", dt_return[UploadFiles.DESCRIPTION])
		self.assertFalse(Constants.SEQ_VIRUS_TYPE in dt_return)
		self.assertEqual("AF100378", dt_return[UploadFiles.ACCESSION])
		
		try:
			dt_return = uploadFiles.parse_title_name("influenza_ype~~~sB~~~AF100378 adadadd")
			self.fail("must throw exception")
		except Exception as e:
			pass


	def test_upload_file(self):
		
		uploadFiles = UploadFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		
		uploadFiles.upload_file(version, file)
		base_name = os.path.basename(file)
		try:
			uploadFile = UploadFile.objects.get(name=base_name)
		except UploadFile.DoesNotExist:	## not exist
			self.fail("Must have value")
			
		try:
			tag = Tags.objects.get(name=Constants.SEQ_VIRUS_TYPE)
		except Tags.DoesNotExist:	## not exist
			self.fail("Must have value")
		
		try:
			tag = Tags.objects.get(name=Constants.SEQ_VIRUS_SUB_TYPE)
		except Tags.DoesNotExist:	## not exist
			self.fail("Must have value")
		
		try:
			seqVirus = SeqVirus.objects.get(name='B')
			self.assertEqual("B", seqVirus.name)
			self.assertEqual(Constants.SEQ_VIRUS_TYPE, seqVirus.kind_type.name)
			self.assertEqual("AF100378", seqVirus.accession)
		except SeqVirus.DoesNotExist:	## not exist
			self.fail("Must have value")
			
		try:
			seqVirus = SeqVirus.objects.get(name='H9')
			self.assertEqual("H9", seqVirus.name)
			self.assertEqual(Constants.SEQ_VIRUS_SUB_TYPE, seqVirus.kind_type.name)
			self.assertEqual("KX879589", seqVirus.accession)
			self.assertEqual("test_db_influenza_typing_v2.fasta", seqVirus.file.name)
			self.assertEqual(2, seqVirus.file.version)
		except SeqVirus.DoesNotExist:	## not exist
			self.fail("Must have value")
		
		
	def test_upload_file_and_results(self):
		
		uploadFiles = UploadFiles()
		parseOutFiles = ParseOutFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		uploadFiles.upload_file(version, file)	## upload file
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
		self.assertEquals('2', version)
		
		uploadFile = UploadFile.objects.order_by('-version')[0]
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_ABRICATE, ConstantsTestsCase.MANAGING_TEST_ABRICATE)
		self.assertTrue(os.path.exists(txt_file))
		vect_data = parseOutFiles.parse_abricate_file(txt_file)
		
		### create an IdentifyVirus instance
		vect_identify_virus = uploadFiles.uploadIdentifyVirus(vect_data, uploadFile.abricate_name)
				
		self.assertEquals(0, vect_identify_virus[0].rank)
		self.assertEquals("100.00", vect_identify_virus[0].coverage)
		self.assertEquals("99.69", vect_identify_virus[0].identity)
		self.assertEquals(Constants.SEQ_VIRUS_TYPE, vect_identify_virus[0].seq_virus.kind_type.name)
		self.assertEquals(Constants.SEQ_VIRUS_SUB_TYPE, vect_identify_virus[1].seq_virus.kind_type.name)
		self.assertEquals("100.00", vect_identify_virus[1].coverage)
		self.assertEquals("99.18", vect_identify_virus[1].identity)
		self.assertEquals("H3", vect_identify_virus[1].seq_virus.name)
		self.assertEquals(Constants.SEQ_VIRUS_SUB_TYPE, vect_identify_virus[1].seq_virus.kind_type.name)
		self.assertEquals("100.00", vect_identify_virus[2].coverage)
		self.assertEquals("98.65", vect_identify_virus[2].identity)
		self.assertEquals("N2", vect_identify_virus[2].seq_virus.name)
		self.assertEquals(Constants.SEQ_VIRUS_SUB_TYPE, vect_identify_virus[2].seq_virus.kind_type.name)
		self.assertEquals("80.00", vect_identify_virus[3].coverage)
		self.assertEquals("78.65", vect_identify_virus[3].identity)
		self.assertEquals("Victoria", vect_identify_virus[3].seq_virus.name)
		self.assertEquals(Constants.SEQ_VIRUS_LINEAGE, vect_identify_virus[3].seq_virus.kind_type.name)
		
		
		
		
		
		
	
		
		