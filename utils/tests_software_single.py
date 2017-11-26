'''
Created on Oct 28, 2017

@author: mmp
'''

from django.test import TestCase
from django.conf import settings 
from utils.constantsTestsCase import ConstantsTestsCase
from utils.constants import Constants, TypePath, FileType, FileExtensions
from utils.meta_key_and_values import MetaKeyAndValue
from utils.software import Software
from utils.software_names import SoftwareNames
from utils.utils import Utils
from utils.parseOutFiles import ParseOutFiles
from utils.result import DecodeResultAverageAndNumberReads, DecodeResult, DecodeCoverage, Coverage
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference
from manage_virus.uploadFiles import UploadFiles
from manage_virus.models import UploadFile
from managing_files.manage_database import ManageDatabase
from django.test.utils import override_settings
from utils.software_names import SoftwareNames
from utils.tree import CreateTree
import os, filecmp

class Test(TestCase):

	### static
	software = Software()
	utils = Utils()
	constants = Constants()
	software_names = SoftwareNames()
	constants_tests_case = ConstantsTestsCase()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass


	def tearDown(self):
		pass

		



# 	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
# 	def test_process_second_stage_analysis_to_produce_data(self):
# 		"""
#  		test global method
#  		"""
# 		manageDatabase = ManageDatabase()
# 		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
# 		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))
# 	
# 		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
# 		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
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
# 		ref_name = "second_stage_2"
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
# 		project_name = "several_names"
# 		try:
# 			project = Project.objects.get(name=project_name)
# 		except Project.DoesNotExist:
# 			project = Project()
# 			project.name = project_name
# 			project.reference = reference
# 			project.owner = user
# 			project.save()
# 		
# 		vect_files = self.constants_tests_case.get_all_fastq_files(self.baseDirectory)
# 		
# 		temp_dir = self.utils.get_temp_dir()
# 		for vect_file in vect_files:
# 			self.utils.copy_file(vect_file[0], os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
# 			self.utils.copy_file(vect_file[1], os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2))
# 				
# 			sample_name = "_".join(os.path.basename(vect_file[0]).split('_')[0:2])
# 			try:
# 				sample = Sample.objects.get(name=sample_name)
# 			except Sample.DoesNotExist:
# 				sample = Sample()
# 				sample.name = sample_name
# 				sample.is_rejected = False
# 				sample.is_valid_1 = True
# 				sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
# 				sample.path_name_1.name = os.path.join(temp_dir, vect_file[0])
# 				sample.is_valid_2 = False
# 				sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
# 				sample.path_name_2.name = os.path.join(temp_dir, vect_file[1])
# 				sample.owner = user
# 				
# 				sample.is_ready_for_projects = True
# 				sample.is_obsolete = False
# 				sample.is_rejected = False
# 				sample.save()
# 
# 			## create project_sample
# 			project_sample = ProjectSample()
# 			project_sample.sample = sample
# 			project_sample.project = project
# 			project_sample.is_finished = True
# 			project_sample.is_deleted = False
# 			project_sample.is_error = False
# 			project_sample.save()
# 		
# ## IF exist the directory
# # 			self.assertTrue(self.software.process_second_stage_snippy_coverage_freebayes(project_sample, user))
# # 			try:
# # 				project_sample = ProjectSample.objects.get(pk=project_sample.id)
# # 			except ProjectSample.DoesNotExist:
# # 				self.fail("Must exist")
# # 			self.assertTrue(project_sample.is_finished)
# 
# 			meta_value = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
# 			if (meta_value == None):
# 				coverage = Coverage()
# 				coverage.add_coverage('MP', Coverage.COVERAGE_ALL, '2198.8')
# 				coverage.add_coverage('MP', Coverage.COVERAGE_MORE_0, '100.0')
# 				coverage.add_coverage('MP', Coverage.COVERAGE_MORE_9, '100.0')
# 				coverage.add_coverage('PA', Coverage.COVERAGE_ALL, '527.8')
# 				coverage.add_coverage('PA', Coverage.COVERAGE_MORE_0, '100.0')
# 				coverage.add_coverage('PA', Coverage.COVERAGE_MORE_9, '100.0')
# 				coverage.add_coverage('HA', Coverage.COVERAGE_ALL, '1449.8')
# 				coverage.add_coverage('HA', Coverage.COVERAGE_MORE_0, '100.0')
# 				coverage.add_coverage('HA', Coverage.COVERAGE_MORE_9, '100.0')
# 				coverage.add_coverage('PB1', Coverage.COVERAGE_ALL, '618.8')
# 				coverage.add_coverage('PB1', Coverage.COVERAGE_MORE_0, '100.0')
# 				coverage.add_coverage('PB1', Coverage.COVERAGE_MORE_9, '100.0')
# 				coverage.add_coverage('NP', Coverage.COVERAGE_ALL, '439.8')
# 				coverage.add_coverage('NP', Coverage.COVERAGE_MORE_0, '100.0')
# 				coverage.add_coverage('NP', Coverage.COVERAGE_MORE_9, '100.0')
# 				coverage.add_coverage('NS', Coverage.COVERAGE_ALL, '1214.8')
# 				coverage.add_coverage('NS', Coverage.COVERAGE_MORE_0, '100.0')
# 				coverage.add_coverage('NS', Coverage.COVERAGE_MORE_9, '100.0')
# 				coverage.add_coverage('PB2', Coverage.COVERAGE_ALL, '690.8')
# 				coverage.add_coverage('PB2', Coverage.COVERAGE_MORE_0, '100.0')
# 				coverage.add_coverage('PB2', Coverage.COVERAGE_MORE_9, '100.0')
# 				coverage.add_coverage('NA', Coverage.COVERAGE_ALL, '1092.8')
# 				coverage.add_coverage('NA', Coverage.COVERAGE_MORE_0, '100.0')
# 				coverage.add_coverage('NA', Coverage.COVERAGE_MORE_9, '100.0')
# 				meta_value = manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success, coverage.to_json())
# 			
# 		##### make the trees
# 		create_tree = CreateTree()
# 		create_tree.create_tree_and_alignments_all(project, user)
# 		
# 		self.utils.remove_dir(temp_dir)
# #		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))

