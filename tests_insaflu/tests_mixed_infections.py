'''
Created on Dec 14, 2017

@author: mmp
'''
import unittest, os
from utils.mixed_infections_management import MixedInfectionsManagement
from managing_files.models import MixedInfections
from utils.result import MixedInfectionMainVector, DecodeObjects, CountHits
from utils.utils import Utils
from django.conf import settings 
from django.test.utils import override_settings
from constants.constantsTestsCase import ConstantsTestsCase
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference
from constants.meta_key_and_values import MetaKeyAndValue
from constants.constants import Constants
from managing_files.manage_database import ManageDatabase
from constants.constants_mixed_infection import ConstantsMixedInfection
from manage_virus.constants_virus import ConstantsVirus
from manage_virus.models import Tags, SeqVirus, IdentifyVirus

class Test(unittest.TestCase):

	utils = Utils()

	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass

	def test_get_mixed_infection_empty_value(self):
	
		mixed_infections_management = MixedInfectionsManagement()
		mixed_infections = mixed_infections_management.get_mixed_infections_empty_value()
		self.assertEqual(Constants.EMPTY_VALUE_NA, mixed_infections.tag.name)
		
	def test_get_mixed_infection_main_vector(self):
		
		try:
			mixed_infections_ = MixedInfections.objects.get(has_master_vector=True)
			mixed_infections_.delete()
		except MixedInfections.DoesNotExist as e:
			pass
		
		mixed_infection_main_vector = MixedInfectionMainVector()
		mixed_infections_management = MixedInfectionsManagement()
		self.assertEquals(mixed_infection_main_vector, mixed_infections_management.get_mixed_infection_main_vector())
		
		try:
			mixed_infections_ = MixedInfections.objects.get(has_master_vector=True)
			decodeResult = DecodeObjects()
			self.assertEquals(decodeResult.decode_result(mixed_infections_.description), mixed_infection_main_vector)
		except MixedInfections.DoesNotExist as e:
			self.fail('must have value')
		
		mixed_infection_main_vector = mixed_infections_management.get_mixed_infection_main_vector()
		mixed_infection_main_vector.add_vector([10, 20])
		mixed_infections_.description = mixed_infection_main_vector.to_json()
		mixed_infections_.save()
		
		try:
			mixed_infections_ = MixedInfections.objects.get(has_master_vector=True)
			decodeResult = DecodeObjects()
			self.assertEquals(decodeResult.decode_result(mixed_infections_.description), mixed_infection_main_vector)
		except MixedInfections.DoesNotExist as e:
			self.fail('must have value')
			
		mixed_infection_main_vector.add_vector([10, 22220])
		self.assertNotEqual(decodeResult.decode_result(mixed_infections_.description), mixed_infection_main_vector)
	
	
	def test_get_value_mixed_infection(self):
		
		try:
			mixed_infections_ = MixedInfections.objects.get(has_master_vector=True)
			mixed_infections_.delete()
		except MixedInfections.DoesNotExist as e:
			pass
		
		count_hits = CountHits()
		count_hits.set_hits_less_50(30)
		count_hits.set_hits_50_90(10)
		self.assertEquals([30, 10], count_hits.get_vect_mixed_infections())
		
		mixed_infections_management = MixedInfectionsManagement()
		self.assertEquals('0.9542543216935453', '{}'.format(mixed_infections_management.get_value_mixed_infection(count_hits)))
		count_hits.set_hits_less_50(75)
		count_hits.set_hits_50_90(66)
		self.assertEquals('0.9899940504994065', '{}'.format(mixed_infections_management.get_value_mixed_infection(count_hits)))
		count_hits.set_hits_less_50(66)
		count_hits.set_hits_50_90(75)
		self.assertEquals('0.9681224546782758', '{}'.format(mixed_infections_management.get_value_mixed_infection(count_hits)))
		count_hits.set_hits_less_50(4)
		count_hits.set_hits_50_90(2)
		self.assertEquals('0.984995168387904', '{}'.format(mixed_infections_management.get_value_mixed_infection(count_hits)))
		count_hits.set_hits_less_50(0)
		count_hits.set_hits_50_90(0)
		self.assertEquals('0.0', '{}'.format(mixed_infections_management.get_value_mixed_infection(count_hits)))
		count_hits.set_hits_less_50(4)
		count_hits.set_hits_50_90(2)
		self.assertEquals('0.984995168387904', '{}'.format(mixed_infections_management.get_value_mixed_infection(count_hits)))

	def test_constants_mixed_infection(self):
		constants_mixed_infection = ConstantsMixedInfection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, constants_mixed_infection.get_tag_by_value(0.99))
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, constants_mixed_infection.get_tag_by_value(0.98))
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_MAY_BE, constants_mixed_infection.get_tag_by_value(0.971))
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_MAY_BE, constants_mixed_infection.get_tag_by_value(0.970))
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, constants_mixed_infection.get_tag_by_value(0.96))
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, constants_mixed_infection.get_tag_by_value(0.95))
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, constants_mixed_infection.get_tag_by_value(0.969))
		
	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_get_mixed_infections(self):
		"""
 		test global method
 		"""
		
		try:
			mixed_infections_ = MixedInfections.objects.get(has_master_vector=True)
			mixed_infections_.delete()
		except MixedInfections.DoesNotExist as e:
			pass
		
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))

		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
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

		ref_name = "second_statest_get_mixed_infections"
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
			
		sample_name = "run_snippytest_get_mixed_infections"
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
			sample.type_of_fastq = Sample.TYPE_OF_FASTQ_illumina
			sample.save()

		project_name = "file_nametest_get_mixed_infections"
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

		## set count hits		
		count_hits = CountHits()
		count_hits.set_hits_less_50(70)
		count_hits.set_hits_50_90(50)
		
		decodeResult = DecodeObjects()
		mixed_infection_main_vector = MixedInfectionMainVector()
		mixed_infections_management = MixedInfectionsManagement()
		mixed_infections = mixed_infections_management.get_mixed_infections(project_sample, user, count_hits)
		self.assertEquals('0.995925821845502', '{}'.format(mixed_infections.average_value))
		self.assertEquals(mixed_infection_main_vector, decodeResult.decode_result(mixed_infections.description))
		self.assertFalse(mixed_infections.has_master_vector)
		self.assertEquals('Yes', mixed_infections.tag.name)
		
		project_sample_ = ProjectSample.objects.get(pk=project_sample.id)
		self.assertEquals(1, project_sample_.alert_first_level)
		
		self.assertEquals(0.7142857142857143, count_hits.get_mixed_infection_ratio())
		self.assertEquals("0.7", count_hits.get_mixed_infection_ratio_str())
		self.assertTrue(count_hits.is_mixed_infection_ratio_test())
		
		### 
# 		manage_database = ManageDatabase()
# 		meta_project = manage_database.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_COSINE_DISTANCE, MetaKeyAndValue.META_VALUE_Success)
# 		self.assertFalse(meta_project == None)
# 		self.assertEquals("Warning, this sample has an average cosine distance of '0.995925821845502'.\nSuggest mixed infection.", meta_project.description)
# 		meta_project = manage_database.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_RATIO_TEST, MetaKeyAndValue.META_VALUE_Success)
# 		self.assertFalse(meta_project == None)
# 		self.assertEquals("Warning, this sample has a ratio of '0.7' and a total of '120' variations.\nSuggest mixed infection.", meta_project.description)
		
		## remove all files
		self.utils.remove_dir(temp_dir)
		
	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_get_mixed_infections_2(self):
		"""
 		test global method
 		"""
		
		try:
			mixed_infections_ = MixedInfections.objects.get(has_master_vector=True)
			mixed_infections_.delete()
		except MixedInfections.DoesNotExist as e:
			pass
		
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))

		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
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

		ref_name = "second_statest_get_mixed_infections"
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
			
		sample_name = "run_snippytest_get_mixed_infections_2323"
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
			sample.type_of_fastq = Sample.TYPE_OF_FASTQ_illumina
			sample.owner = user
			sample.save()

		project_name = "file_nametest_get_mixed_infections"
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

		manage_database = ManageDatabase()
		manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE,\
								MetaKeyAndValue.META_VALUE_Success, 'this is a mixed infection of type and subtype')
		
		## set count hits		
		count_hits = CountHits()
		count_hits.set_hits_less_50(50)
		count_hits.set_hits_50_90(70)
		
		decodeResult = DecodeObjects()
		mixed_infection_main_vector = MixedInfectionMainVector()
		mixed_infections_management = MixedInfectionsManagement()
		mixed_infections = mixed_infections_management.get_mixed_infections(project_sample, user, count_hits)
		self.assertEquals('0.9394790479330498', '{}'.format(mixed_infections.average_value))
		self.assertEquals(mixed_infection_main_vector, decodeResult.decode_result(mixed_infections.description))
		self.assertFalse(mixed_infections.has_master_vector)
		self.assertEquals('Yes', mixed_infections.tag.name)
		
		project_sample_ = ProjectSample.objects.get(pk=project_sample.id)
		self.assertEquals(1, project_sample_.alert_first_level)
		
		self.assertEquals(1.4, count_hits.get_mixed_infection_ratio())
		self.assertEquals("1.4", count_hits.get_mixed_infection_ratio_str())
		self.assertTrue(count_hits.is_mixed_infection_ratio_test())
		
		meta_data = manage_database.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_SUM_TEST,\
								MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_data == None)
		meta_data = manage_database.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_RATIO_TEST,\
								MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_data != None)
		self.assertEquals("Warning: this sample has a ratio of the number of iSNVs at frequency 1-50% (minor iSNVs) and 50-90% of '1.4' " +\
						"(within the range 0.5-2.0) and a total number of iSNVs from the two categories of '120' (i.e., above 20) suggesting " +\
						"that may represent a 'mixed infection'.", meta_data.description)

		meta_data = manage_database.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE,\
								MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_data != None)
		self.assertEquals('this is a mixed infection of type and subtype', meta_data.description)

		### 
# 		manage_database = ManageDatabase()
# 		meta_project = manage_database.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_COSINE_DISTANCE, MetaKeyAndValue.META_VALUE_Success)
# 		self.assertFalse(meta_project == None)
# 		self.assertEquals("Warning, this sample has an average cosine distance of '0.9394790479330498'.\nSuggest mixed infection.", meta_project.description)
# 		meta_project = manage_database.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_RATIO_TEST, MetaKeyAndValue.META_VALUE_Success)
# 		self.assertFalse(meta_project == None)
# 		self.assertEquals("Warning, this sample has a ratio of '1.4' and a total of '120' variations.\nSuggest mixed infection.", meta_project.description)

		## remove all files
		self.utils.remove_dir(temp_dir)

	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_get_mixed_infections_3(self):
		"""
 		test global method
 		"""
		
		try:
			mixed_infections_ = MixedInfections.objects.get(has_master_vector=True)
			mixed_infections_.delete()
		except MixedInfections.DoesNotExist as e:
			pass
		
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))

		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
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

		ref_name = "second_statest_get_mixed_infections"
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
			
		sample_name = "run_snippytest_get_mixed_infections"
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
			sample.type_of_fastq = Sample.TYPE_OF_FASTQ_illumina
			sample.owner = user
			sample.save()

		project_name = "file_nametest_get_mixed_infections"
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

		## set count hits		
		count_hits = CountHits()
		count_hits.set_hits_less_50(1)
		count_hits.set_hits_50_90(3)
		
		decodeResult = DecodeObjects()
		mixed_infection_main_vector = MixedInfectionMainVector()
		mixed_infections_management = MixedInfectionsManagement()
		mixed_infections = mixed_infections_management.get_mixed_infections(project_sample, user, count_hits)
		self.assertEquals('0.8007024955568505', '{}'.format(mixed_infections.average_value))
		self.assertEquals(mixed_infection_main_vector, decodeResult.decode_result(mixed_infections.description))
		self.assertFalse(mixed_infections.has_master_vector)
		self.assertEquals('No', mixed_infections.tag.name)
		
		project_sample_ = ProjectSample.objects.get(pk=project_sample.id)
		self.assertEquals(0, project_sample_.alert_first_level)
		
		self.assertEquals(3.0, count_hits.get_mixed_infection_ratio())
		self.assertEquals("3.0", count_hits.get_mixed_infection_ratio_str())
		self.assertFalse(count_hits.is_mixed_infection_ratio_test())
		
		### 
		manage_database = ManageDatabase()
		meta_project = manage_database.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_COSINE_DISTANCE, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_project == None)
		meta_project = manage_database.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_RATIO_TEST, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_project == None)
		meta_project = manage_database.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_SUM_TEST, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_project == None)
		
		## remove all files
		self.utils.remove_dir(temp_dir)
		
		
	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_get_mixed_infections_4(self):
		"""
 		test global method
 		"""
		
		try:
			mixed_infections_ = MixedInfections.objects.get(has_master_vector=True)
			mixed_infections_.delete()
		except MixedInfections.DoesNotExist as e:
			pass
		
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))

		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
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

		ref_name = "second_statest_get_mixed_infections"
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
			
		sample_name = "run_snippytest_get_mixed_infections"
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
			sample.type_of_fastq = Sample.TYPE_OF_FASTQ_illumina
			sample.owner = user
			sample.save()

		project_name = "file_nametest_get_mixed_infections"
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

		## set count hits		
		count_hits = CountHits()
		count_hits.set_hits_less_50(200)
		count_hits.set_hits_50_90(50)
		
		decodeResult = DecodeObjects()
		mixed_infection_main_vector = MixedInfectionMainVector()
		mixed_infections_management = MixedInfectionsManagement()
		mixed_infections = mixed_infections_management.get_mixed_infections(project_sample, user, count_hits)
		self.assertEquals('0.9295706783216093', '{}'.format(mixed_infections.average_value))
		self.assertEquals(mixed_infection_main_vector, decodeResult.decode_result(mixed_infections.description))
		self.assertFalse(mixed_infections.has_master_vector)
		self.assertEquals('Yes', mixed_infections.tag.name)
		
		project_sample_ = ProjectSample.objects.get(pk=project_sample.id)
		self.assertEquals(1, project_sample_.alert_first_level)
		
		self.assertEquals(0.25, count_hits.get_mixed_infection_ratio())
		self.assertEquals("0.2", count_hits.get_mixed_infection_ratio_str())
		self.assertFalse(count_hits.is_mixed_infection_ratio_test())
		
		### 
		manage_database = ManageDatabase()
		meta_project = manage_database.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_COSINE_DISTANCE, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_project == None)
		meta_project = manage_database.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_RATIO_TEST, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_project == None)
		meta_project = manage_database.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_SUM_TEST, MetaKeyAndValue.META_VALUE_Success)
		self.assertFalse(meta_project == None)
		
		## remove all files
		self.utils.remove_dir(temp_dir)

	
	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_get_mixed_infections_sentences(self):
		"""
 		test global method
 		"""
		
		try:
			mixed_infections_ = MixedInfections.objects.get(has_master_vector=True)
			mixed_infections_.delete()
		except MixedInfections.DoesNotExist as e:
			pass
		
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))

		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
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

		ref_name = "second_statest_get_mixed_infections"
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
		
		### test with none
		sample = self.get_sample("run_snippytest_get_mixed_infections_11", user, temp_dir)
		self.assertEquals('Not assigned', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: no typing data was obtained (possible reason: low number of influenza reads).", data_result[2])

		### test normal A-H1N1
		sample = self.get_sample("run_snippytest_get_mixed_infections_12", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A, ConstantsVirus.TYPE_A], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1', 'N1']})
		self.assertEquals('A-H1N1', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
		self.assertEquals(0, data_result[1])
		self.assertEquals(None, data_result[2])
		
		### test  A-H1N1H3
		sample = self.get_sample("run_snippytest_get_mixed_infections_13", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1', 'N1', 'N3']})
		self.assertEquals('A-H1|N1|N3', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: more than two subtypes were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])

		### test  A-H1N1 YAMAGATA
		sample = self.get_sample("run_snippytest_get_mixed_infections_14", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1', 'N1'],\
									ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata']})
		self.assertEquals('A-H1N1; Yamagata', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])

		### test  A YAMAGATA
		sample = self.get_sample("run_snippytest_get_mixed_infections_15", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A], ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata']})
		self.assertEquals('A; Yamagata', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])

		### test  A-H1
		sample = self.get_sample("run_snippytest_get_mixed_infections_16", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1']})
		self.assertEquals('A-H1', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: an incomplete subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).", data_result[2])
		
		### test  A-N1
		sample = self.get_sample("run_snippytest_get_mixed_infections_16b", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['N1']})
		self.assertEquals('A-N1', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: an incomplete subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).", data_result[2])
		
		### test  A
		sample = self.get_sample("run_snippytest_get_mixed_infections_16a", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A]})
		self.assertEquals('A', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: an incomplete subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).", data_result[2])

		### test  B-H1 Yamagata
		sample = self.get_sample("run_snippytest_get_mixed_infections_17", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_B], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1'],\
									ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata']})
		self.assertEquals('H1; B-Yamagata', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: more than one type/subtypes were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])
		
		### test normal B-Yamagata
		sample = self.get_sample("run_snippytest_get_mixed_infections_18", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_B], ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata']})
		self.assertEquals('B-Yamagata', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
		self.assertEquals(0, data_result[1])
		self.assertEquals(None, data_result[2])
		
		### test  B-Yamagata-Xpto
		sample = self.get_sample("run_snippytest_get_mixed_infections_19", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_B], ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata', 'Xpto']})
		self.assertEquals('B-XptoYamagata', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: more than one lineage were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])
		
		### test  B
		sample = self.get_sample("run_snippytest_get_mixed_infections_20", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_B]})
		self.assertEquals('B', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: an incomplete lineage has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).", data_result[2])

		### test  A-B
		sample = self.get_sample("run_snippytest_get_mixed_infections_21", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A, ConstantsVirus.TYPE_B]})
		self.assertEquals('A; B', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])

		### test  A-B H1
		sample = self.get_sample("run_snippytest_get_mixed_infections_22", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A, ConstantsVirus.TYPE_B], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1']})
		self.assertEquals('A-H1; B', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])

		### test  A-B Yamagata
		sample = self.get_sample("run_snippytest_get_mixed_infections_23", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A, ConstantsVirus.TYPE_B], ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata']})
		self.assertEquals('A; B-Yamagata', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])

		### test  A-B H1 Yamagata
		sample = self.get_sample("run_snippytest_get_mixed_infections_24", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A, ConstantsVirus.TYPE_B], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1'],\
							ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata']})
		self.assertEquals('A-H1; B-Yamagata', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])


		### test H1 Yamagata
		sample = self.get_sample("run_snippytest_get_mixed_infections_25", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1'], ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata']})
		self.assertEquals('H1; Yamagata', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])

		### test  H1 
		sample = self.get_sample("run_snippytest_get_mixed_infections_26", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1'] })
		self.assertEquals('H1', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: an incomplete type/subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).", data_result[2])

		### test  Yamagata
		sample = self.get_sample("run_snippytest_get_mixed_infections_27", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata']})
		self.assertEquals('Yamagata', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: an incomplete type/subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).", data_result[2])
	
		### test  H1N2 
		sample = self.get_sample("run_snippytest_get_mixed_infections_28", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1', 'N2'] })
		self.assertEquals('H1N2', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: an incomplete type/subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).", data_result[2])

		### test  N1N2 
		sample = self.get_sample("run_snippytest_get_mixed_infections_29", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['N1', 'N2'] })
		self.assertEquals('N1N2', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])

		### test  N1N2  H1H2 
		sample = self.get_sample("run_snippytest_get_mixed_infections_30", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['N1', 'N2', 'H1', 'H2'] })
		self.assertEquals('H1|H2|N1|N2', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])

		### test  N1Yamagata
		sample = self.get_sample("run_snippytest_get_mixed_infections_31", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata'], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['N1'] })
		self.assertEquals('N1; Yamagata', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])


		#### corona virus
		sample = self.get_sample("run_snippytest_get_mixed_infections_32", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_GENUS : ['BetaCoV'], ConstantsVirus.SEQ_VIRUS_HUMAN : ['SARS_CoV_2'] })
		self.assertEquals('BetaCoV-SARS_CoV_2', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
		self.assertEquals(0, data_result[1])
		self.assertEquals(None, data_result[2])
		
		sample = self.get_sample("run_snippytest_get_mixed_infections_33", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_GENUS : ['BetaCoV'], ConstantsVirus.SEQ_VIRUS_HUMAN : ['SARS_CoV_2', 'SARS_CoV_24'] })
		self.assertEquals('BetaCoV-SARS_CoV_2|SARS_CoV_24', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: more than one human BetaCoV virus are likely present in this sample, suggesting that may represent a 'mixed infection'", data_result[2])

		sample = self.get_sample("run_snippytest_get_mixed_infections_34", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_GENUS : ['BetaCoV'] })
		self.assertEquals('BetaCoV', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: an incomplete human BetaCoV identification has been obtained (possible reasons: low number of  human BetaCoV reads, mixed infection, etc)", data_result[2])

		sample = self.get_sample("run_snippytest_get_mixed_infections_35", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_HUMAN : ['SARS_CoV_2'] })
		self.assertEquals('SARS_CoV_2', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: an incomplete human BetaCoV identification has been obtained (possible reasons: low number of  human BetaCoV reads, mixed infection, etc)", data_result[2])

		sample = self.get_sample("run_snippytest_get_mixed_infections_36", user, temp_dir)
		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_HUMAN : ['SARS_CoV_2', 'SARS_CoV_24'] })
		self.assertEquals('SARS_CoV_2|SARS_CoV_24', sample.get_type_sub_type())
		data_result = sample.get_mixed_infection()
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
		self.assertEquals(1, data_result[1])
		self.assertEquals("Warning: more than one human BetaCoV virus are likely present in this sample, suggesting that may represent a 'mixed infection'", data_result[2])


		## remove all files
		self.utils.remove_dir(temp_dir)



	def get_sample(self, name, user, temp_dir):
		"""
		return a test sample
		"""
		sample_name = name
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
			sample.type_of_fastq = Sample.TYPE_OF_FASTQ_illumina
			sample.owner = user
			sample.save()
		return sample


	def add_type_sub_type(self, sample, dict_data):
		"""
		add type and sub type
		"""
		n_rank = 0
		for key_type in dict_data:
			
			try:
				tag = Tags.objects.get(name=key_type)
			except Tags.DoesNotExist as e:
				tag = Tags()
				tag.name = key_type;
				tag.save()
			
			for name in dict_data[key_type]:
				try:
					seq_virus = SeqVirus.objects.get(name=name, kind_type__name=key_type)
				except SeqVirus.DoesNotExist as e:
					seq_virus = SeqVirus()
					seq_virus.name = name
					seq_virus.kind_type = tag
					seq_virus.save()
					
				identify_virus = IdentifyVirus()
				identify_virus.seq_virus = seq_virus
				identify_virus.rank = n_rank
				identify_virus.save()
				n_rank += 1
				
				sample.identify_virus.add(identify_virus)
		sample.save()
		return sample


