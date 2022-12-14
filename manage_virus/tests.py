from django.test import TestCase

# Create your tests here.
from constants.constantsTestsCase import ConstantsTestsCase
from django.conf import settings 
from utils.utils import Utils
from .uploadFiles import UploadFiles
from constants.constants import Constants, TypePath, FileExtensions
from utils.parse_out_files import ParseOutFiles
from django.contrib.auth.models import User
from .models import UploadFile, Tags, SeqVirus
from managing_files.models import Reference
from django.test.utils import override_settings
from manage_virus.constants_virus import ConstantsVirus
import os, filecmp

class Test(TestCase):

	### static
	constantsTestsCase = ConstantsTestsCase()
	utils = Utils()
	
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
		self.assertTrue(ConstantsVirus.SEQ_VIRUS_TYPE in dt_return)
		self.assertEqual("B", dt_return[ConstantsVirus.SEQ_VIRUS_TYPE])
		self.assertFalse(ConstantsVirus.SEQ_VIRUS_SUB_TYPE in dt_return)
		self.assertFalse(UploadFiles.DESCRIPTION in dt_return)
		self.assertEqual("AF100378", dt_return[UploadFiles.ACCESSION])
		
		dt_return = uploadFiles.parse_title_name("influenza_subtype~~~sB~~~AF100378 adadadd")
		self.assertTrue(ConstantsVirus.SEQ_VIRUS_SUB_TYPE in dt_return)
		self.assertEqual("sB", dt_return[ConstantsVirus.SEQ_VIRUS_SUB_TYPE])
		self.assertEqual("adadadd", dt_return[UploadFiles.DESCRIPTION])
		self.assertFalse(ConstantsVirus.SEQ_VIRUS_TYPE in dt_return)
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
			tag = Tags.objects.get(name=ConstantsVirus.SEQ_VIRUS_TYPE)
		except Tags.DoesNotExist:	## not exist
			self.fail("Must have value")
		
		try:
			tag = Tags.objects.get(name=ConstantsVirus.SEQ_VIRUS_SUB_TYPE)
		except Tags.DoesNotExist:	## not exist
			self.fail("Must have value")
		
		try:
			seqVirus = SeqVirus.objects.get(name='B')
			self.assertEqual("B", seqVirus.name)
			self.assertEqual(ConstantsVirus.SEQ_VIRUS_TYPE, seqVirus.kind_type.name)
			self.assertEqual("AF100378", seqVirus.accession)
		except SeqVirus.DoesNotExist:	## not exist
			self.fail("Must have value")
			
		try:
			seqVirus = SeqVirus.objects.get(name='H9')
			self.assertEqual("H9", seqVirus.name)
			self.assertEqual(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, seqVirus.kind_type.name)
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
		(dict_out_abricate, clean_abricate_file) = parseOutFiles.parse_abricate_file(txt_file, 'test.txt', 1)
		
		self.assertTrue(filecmp.cmp(txt_file, clean_abricate_file))
		os.unlink(clean_abricate_file)
		
		### create an IdentifyVirus instance
		vect_identify_virus = uploadFiles.uploadIdentifyVirus(dict_out_abricate, uploadFile.abricate_name)
				
		self.assertEquals(0, vect_identify_virus[0].rank)
		self.assertEquals("100.00", vect_identify_virus[0].coverage)
		self.assertEquals("99.69", vect_identify_virus[0].identity)
		self.assertEquals(ConstantsVirus.SEQ_VIRUS_TYPE, vect_identify_virus[0].seq_virus.kind_type.name)
		self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, vect_identify_virus[1].seq_virus.kind_type.name)
		self.assertEquals("100.00", vect_identify_virus[1].coverage)
		self.assertEquals("99.18", vect_identify_virus[1].identity)
		self.assertEquals("H3", vect_identify_virus[1].seq_virus.name)
		self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, vect_identify_virus[1].seq_virus.kind_type.name)
		self.assertEquals("100.00", vect_identify_virus[2].coverage)
		self.assertEquals("98.65", vect_identify_virus[2].identity)
		self.assertEquals("N2", vect_identify_virus[2].seq_virus.name)
		self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, vect_identify_virus[2].seq_virus.kind_type.name)
		self.assertEquals("80.00", vect_identify_virus[3].coverage)
		self.assertEquals("78.65", vect_identify_virus[3].identity)
		self.assertEquals("Victoria", vect_identify_virus[3].seq_virus.name)
		self.assertEquals(ConstantsVirus.SEQ_VIRUS_LINEAGE, vect_identify_virus[3].seq_virus.kind_type.name)
		
	def test_upload_file_and_results_2(self):
		
		uploadFiles = UploadFiles()
		parseOutFiles = ParseOutFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		uploadFiles.upload_file(version, file)	## upload file
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
		self.assertEquals('2', version)
		
		uploadFile = UploadFile.objects.order_by('-version')[0]
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_ABRICATE, ConstantsTestsCase.MANAGING_TEST_ABRICATE_2)
		self.assertTrue(os.path.exists(txt_file))
		(dict_out_abricate, clean_abricate_file) = parseOutFiles.parse_abricate_file(txt_file, 'test.txt', 1)
		
		self.assertTrue(filecmp.cmp(txt_file, clean_abricate_file))
		os.unlink(clean_abricate_file)
		
		### create an IdentifyVirus instance
		vect_identify_virus = uploadFiles.uploadIdentifyVirus(dict_out_abricate, uploadFile.abricate_name)
				
		self.assertEquals(0, vect_identify_virus[0].rank)
		self.assertEquals("100.00", vect_identify_virus[0].coverage)
		self.assertEquals("99.59", vect_identify_virus[0].identity)
		self.assertEquals("A", vect_identify_virus[0].seq_virus.name)
		self.assertEquals(ConstantsVirus.SEQ_VIRUS_TYPE, vect_identify_virus[0].seq_virus.kind_type.name)
		self.assertEquals("47.26", vect_identify_virus[1].coverage)
		self.assertEquals("95.47", vect_identify_virus[1].identity)
		self.assertEquals("B", vect_identify_virus[1].seq_virus.name)
		self.assertEquals(ConstantsVirus.SEQ_VIRUS_TYPE, vect_identify_virus[1].seq_virus.kind_type.name)
		
		self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, vect_identify_virus[2].seq_virus.kind_type.name)
		self.assertEquals("100.00", vect_identify_virus[2].coverage)
		self.assertEquals("98.82", vect_identify_virus[2].identity)
		self.assertEquals("H3", vect_identify_virus[2].seq_virus.name)
		self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, vect_identify_virus[2].seq_virus.kind_type.name)
		self.assertEquals("90.92", vect_identify_virus[3].coverage)
		self.assertEquals("97.66", vect_identify_virus[3].identity)
		self.assertEquals("N2", vect_identify_virus[3].seq_virus.name)
		self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, vect_identify_virus[3].seq_virus.kind_type.name)
		self.assertEquals("12.75", vect_identify_virus[4].coverage)
		self.assertEquals("92.34", vect_identify_virus[4].identity)
		self.assertEquals("Yamagata", vect_identify_virus[4].seq_virus.name)
		self.assertEquals(ConstantsVirus.SEQ_VIRUS_LINEAGE, vect_identify_virus[4].seq_virus.kind_type.name)
		
	def test_upload_file_and_results_3(self):
		
		uploadFiles = UploadFiles()
		parseOutFiles = ParseOutFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		uploadFiles.upload_file(version, file)	## upload file
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
		self.assertEquals('2', version)
		
		uploadFile = UploadFile.objects.order_by('-version')[0]
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_ABRICATE, ConstantsTestsCase.MANAGING_TEST_ABRICATE_2)
		self.assertTrue(os.path.exists(txt_file))
		(dict_out_abricate, clean_abricate_file) = parseOutFiles.parse_abricate_file(txt_file, 'test.txt', 10)
		
		## clean abricate
		expected_abricate = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "expected_abricate_clean.txt")
		
		self.assertTrue(filecmp.cmp(expected_abricate, clean_abricate_file))
		os.unlink(clean_abricate_file)
		
		### create an IdentifyVirus instance
		vect_identify_virus = uploadFiles.uploadIdentifyVirus(dict_out_abricate, uploadFile.abricate_name)
			
		self.assertEqual(3, len(vect_identify_virus))
		self.assertEquals(0, vect_identify_virus[0].rank)
		self.assertEquals("100.00", vect_identify_virus[0].coverage)
		self.assertEquals("99.59", vect_identify_virus[0].identity)
		self.assertEquals("A", vect_identify_virus[0].seq_virus.name)
		self.assertEquals(ConstantsVirus.SEQ_VIRUS_TYPE, vect_identify_virus[0].seq_virus.kind_type.name)
		
		self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, vect_identify_virus[1].seq_virus.kind_type.name)
		self.assertEquals("100.00", vect_identify_virus[1].coverage)
		self.assertEquals("98.82", vect_identify_virus[1].identity)
		self.assertEquals("H3", vect_identify_virus[1].seq_virus.name)
		self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, vect_identify_virus[1].seq_virus.kind_type.name)
		self.assertEquals("90.92", vect_identify_virus[2].coverage)
		self.assertEquals("97.66", vect_identify_virus[2].identity)
		self.assertEquals("N2", vect_identify_virus[2].seq_virus.name)
		self.assertEquals(ConstantsVirus.SEQ_VIRUS_SUB_TYPE, vect_identify_virus[2].seq_virus.kind_type.name)
		
	def test_upload_file_and_results_4(self):
		
		uploadFiles = UploadFiles()
		parseOutFiles = ParseOutFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		uploadFiles.upload_file(version, file)	## upload file
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
		self.assertEquals('2', version)
		
		uploadFile = UploadFile.objects.order_by('-version')[0]
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_ABRICATE, "abricate_out_3.txt")
		self.assertTrue(os.path.exists(txt_file))
		(dict_out_abricate, clean_abricate_file) = parseOutFiles.parse_abricate_file(txt_file, 'test.txt', 10)
		
		### create an IdentifyVirus instance
		vect_identify_virus = uploadFiles.uploadIdentifyVirus(dict_out_abricate, uploadFile.abricate_name)
			
		self.assertEqual(1, len(vect_identify_virus))
		self.assertEquals(0, vect_identify_virus[0].rank)
		self.assertEquals("100.00", vect_identify_virus[0].coverage)
		self.assertEquals("100.00", vect_identify_virus[0].identity)
		self.assertEquals("RSV_A_GA1", vect_identify_virus[0].seq_virus.name)
		self.assertEquals(ConstantsVirus.SEQ_VIRUS_GENOTYPE, vect_identify_virus[0].seq_virus.kind_type.name)
		
	
	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_upload_references(self):
		"""
		"""
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT_TEST", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))
		uploadFiles = UploadFiles()
		
		try:
			User.objects.get(username=Constants.DEFAULT_USER)
			### great, the default user exist 
		except User.DoesNotExist:
			
			### need to create it
			user = User()
			user.username = Constants.DEFAULT_USER
			user.password = Constants.DEFAULT_USER_PASS
			user.first_name = Constants.DEFAULT_USER
			user.is_active = False
			user.is_staff = False
			user.is_superuser = False
			user.save()
			
		b_test = True
		uploadFiles.upload_default_references(user, b_test)
		
		path_to_find = os.path.join(getattr(settings, "STATIC_ROOT", None), Constants.DIR_TEST_TYPE_REFERENCES)
		n_files = 0
		for file in self.utils.get_all_files(path_to_find):
			n_files += 1
			
			name = self.utils.clean_extension(os.path.basename(file))
			try:
				reference = Reference.objects.get(owner=user, is_obsolete=False, is_deleted=False, name__iexact=name)
			except Reference.DoesNotExist as e:
				self.fail("must exist")
			
			self.assertTrue(os.path.exists(reference.get_reference_fasta_index(TypePath.MEDIA_ROOT)))
			self.assertTrue(os.path.exists(reference.get_reference_fasta(TypePath.MEDIA_ROOT)))
			self.assertTrue(os.path.exists(reference.get_reference_gbk(TypePath.MEDIA_ROOT)))
			self.assertTrue(reference.get_gff3(TypePath.MEDIA_ROOT).endswith(FileExtensions.FILE_GFF3))
			self.assertTrue(reference.get_gff3_comulative_positions(TypePath.MEDIA_ROOT).
						endswith(".comulative_positions" + FileExtensions.FILE_GFF3))
			self.assertTrue(os.path.exists(reference.get_gff3(TypePath.MEDIA_ROOT)))
			self.assertTrue(os.path.exists(reference.get_gff3_comulative_positions(TypePath.MEDIA_ROOT)))
		self.assertEqual(4, n_files)	## it has a gbk file
	#	self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))

	
		
		