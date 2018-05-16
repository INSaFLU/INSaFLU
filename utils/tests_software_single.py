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

# 	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
# 	def test_get_mixed_infections_sentences(self):
# 		"""
#  		test global method
#  		"""
# 		
# 		try:
# 			mixed_infections_ = MixedInfections.objects.get(has_master_vector=True)
# 			mixed_infections_.delete()
# 		except MixedInfections.DoesNotExist as e:
# 			pass
# 		
# 		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
# 		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))
# 
# 		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
# 		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
# 		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
# 		file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)
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
# 		ref_name = "second_statest_get_mixed_infections"
# 		try:
# 			reference = Reference.objects.get(name=ref_name)
# 		except Reference.DoesNotExist:
# 			reference = Reference()
# 			reference.name = ref_name
# 			reference.reference_fasta.name = fasta_file
# 			reference.reference_fasta_name = os.path.basename(fasta_file)
# 			reference.reference_genbank.name = gb_file
# 			reference.reference_genbank_name = os.path.basename(gb_file)
# 			reference.owner = user
# 			reference.save()
# 			
# 		temp_dir = self.utils.get_temp_dir()
# 		self.utils.copy_file(file_1, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
# 		self.utils.copy_file(file_2, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2))
# 		
# 		### test with none
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_11", user, temp_dir)
# 		self.assertEquals('Not assigned', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: no type/subtype has been assigned (possible reason: low number of influenza reads).", data_result[2])
# 
# 		### test normal A-H1N1
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_12", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1', 'N1']})
# 		self.assertEquals('A-H1N1', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
# 		self.assertEquals(0, data_result[1])
# 		self.assertEquals(None, data_result[2])
# 		
# 		### test  A-H1N1H3
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_13", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1', 'N1', 'N3']})
# 		self.assertEquals('A-H1|N1|N3', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: more than two subtypes were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])
# 
# 		### test  A-H1N1 YAMAGATA
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_14", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1', 'N1'],\
# 									ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata']})
# 		self.assertEquals('A-H1N1; Yamagata', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])
# 
# 		### test  A YAMAGATA
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_15", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A], ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata']})
# 		self.assertEquals('A; Yamagata', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])
# 
# 		### test  A-H1
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_16", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1']})
# 		self.assertEquals('A-H1', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: an incomplete subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).", data_result[2])
# 		
# 		### test  A-N1
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_16b", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['N1']})
# 		self.assertEquals('A-N1', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: an incomplete subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).", data_result[2])
# 		
# 		### test  A
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_16a", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A]})
# 		self.assertEquals('A', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: an incomplete subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).", data_result[2])
# 
# 		### test  B-H1 Yamagata
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_17", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_B], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1'],\
# 									ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata']})
# 		self.assertEquals('H1; B-Yamagata', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: more than one type/subtypes were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])
# 		
# 		### test normal B-Yamagata
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_18", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_B], ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata']})
# 		self.assertEquals('B-Yamagata', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
# 		self.assertEquals(0, data_result[1])
# 		self.assertEquals(None, data_result[2])
# 		
# 		### test  B-Yamagata-Xpto
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_19", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_B], ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata', 'Xpto']})
# 		self.assertEquals('B-XptoYamagata', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: more than one lineage were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])
# 		
# 		### test  B
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_20", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_B]})
# 		self.assertEquals('B', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: an incomplete lineage has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).", data_result[2])
# 
# 		### test  A-B
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_21", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A, ConstantsVirus.TYPE_B]})
# 		self.assertEquals('A; B', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])
# 
# 		### test  A-B H1
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_22", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A, ConstantsVirus.TYPE_B], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1']})
# 		self.assertEquals('A-H1; B', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])
# 
# 		### test  A-B Yamagata
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_23", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A, ConstantsVirus.TYPE_B], ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata']})
# 		self.assertEquals('A; B-Yamagata', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])
# 
# 		### test  A-B H1 Yamagata
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_24", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_TYPE : [ConstantsVirus.TYPE_A, ConstantsVirus.TYPE_B], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1'],\
# 							ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata']})
# 		self.assertEquals('A-H1; B-Yamagata', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])
# 
# 
# 		### test H1 Yamagata
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_25", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1'], ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata']})
# 		self.assertEquals('H1; Yamagata', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])
# 
# 		### test  H1 
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_26", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1'] })
# 		self.assertEquals('H1', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: an incomplete type/subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).", data_result[2])
# 
# 		### test  Yamagata
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_27", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata']})
# 		self.assertEquals('Yamagata', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: an incomplete type/subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).", data_result[2])
# 	
# 		### test  H1N2 
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_28", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['H1', 'N2'] })
# 		self.assertEquals('H1N2', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: an incomplete type/subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).", data_result[2])
# 
# 		### test  N1N2 
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_29", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['N1', 'N2'] })
# 		self.assertEquals('N1N2', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])
# 
# 		### test  N1Yamagata
# 		sample = self.get_sample("run_snippytest_get_mixed_infections_30", user, temp_dir)
# 		self.add_type_sub_type(sample, { ConstantsVirus.SEQ_VIRUS_LINEAGE : ['Yamagata'], ConstantsVirus.SEQ_VIRUS_SUB_TYPE : ['N1'] })
# 		self.assertEquals('N1; Yamagata', sample.get_type_sub_type())
# 		data_result = sample.get_mixed_infection()
# 		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, data_result[0])
# 		self.assertEquals(1, data_result[1])
# 		self.assertEquals("Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.", data_result[2])
# 
# 		## remove all files
# 		self.utils.remove_dir(temp_dir)
