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
from managing_files.manage_database import ManageDatabase
from constants.constants_mixed_infection import ConstantsMixedInfection

class Test(unittest.TestCase):

	utils = Utils()

	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass

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
			sample.is_rejected = False
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = False
			sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
			sample.path_name_2.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2)
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
		
		### 
		manage_database = ManageDatabase()
		meta_project = manage_database.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION, MetaKeyAndValue.META_VALUE_Success)
		self.assertFalse(meta_project == None)
		self.assertEquals("Warning, this sample has an average cosine distance of '0.995925821845502'\nSuggest mixed infection", meta_project.description)
		
		
		


