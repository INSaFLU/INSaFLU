'''
Created on Nov 25, 2017

@author: mmp
'''
import unittest, os
from django.conf import settings 
from utils.software import Software
from utils.utils import Utils
from constants.constantsTestsCase import ConstantsTestsCase
from constants.constants import Constants, TypePath
from django.test.utils import override_settings
from managing_files.manage_database import ManageDatabase
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference
from utils.result import DecodeObjects
from constants.meta_key_and_values import MetaKeyAndValue
from utils.result import Coverage
from utils.tree import CreateTree
from constants.software_names import SoftwareNames

class Test(unittest.TestCase):


	### static
	software = Software()
	utils = Utils()
	constants = Constants()
	constants_tests_case = ConstantsTestsCase()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass


	def tearDown(self):
		pass


	def test_Name(self):
		pass


	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_create_tree_and_alignments(self):
		"""
 		test global method
 		"""
		manageDatabase = ManageDatabase()
		metaKeyAndValue = MetaKeyAndValue()
		
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))
	
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME + '5000')
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME + '5000'
			user.id = 5000
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		ref_name = "second_stage_2_ test_create_tree"
		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.display_name = ref_name
			reference.reference_fasta.name = fasta_file
			reference.reference_fasta_name = os.path.basename(fasta_file)
			reference.reference_genbank.name = gb_file
			reference.reference_genbank_name = os.path.basename(gb_file)
			reference.owner = user
			reference.save()
		
		project_name = "several_names_test_create_tree"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.id= 5000
			project.name = project_name
			project.reference = reference
			project.owner = user
			project.save()
		
		## get all fastq files
		vect_files = self.constants_tests_case.get_all_fastq_files(self.baseDirectory)
		
		ProjectSample.objects.all().delete()
		Sample.objects.all().delete()
		temp_dir = self.utils.get_temp_dir()
		count = 0
		b_at_least_one_less_than_100 = False
		for vect_file in vect_files:
			self.utils.copy_file(vect_file[0], os.path.join(temp_dir, os.path.basename(vect_file[0])))
			self.utils.copy_file(vect_file[1], os.path.join(temp_dir, os.path.basename(vect_file[1])))
				
			sample_name = "_".join(os.path.basename(vect_file[0]).split('_')[0:2])
			sample = Sample()
			sample.id = 5000 + count + 1
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = vect_file[0]
			sample.path_name_1.name = os.path.join(temp_dir, os.path.basename(vect_file[0]))
			sample.is_valid_2 = False
			sample.file_name_2 = vect_file[1]
			sample.path_name_2.name = os.path.join(temp_dir, os.path.basename(vect_file[1]))
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina)
			sample.owner = user
			sample.is_ready_for_projects = True
			sample.is_obsolete = False
			sample.save()

			## create project_sample
			project_sample = ProjectSample()
			project_sample.id = 5000 + count + 1
			project_sample.sample = sample
			project_sample.project = project
			project_sample.is_finished = True
			project_sample.is_deleted = False
			project_sample.is_error = False
			project_sample.save()
			count += 1

			meta_value = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
			self.assertTrue(meta_value == None)
			
			coverage = Coverage()
			coverage.add_coverage('MP', Coverage.COVERAGE_ALL, '2198.8')
			coverage.add_coverage('MP', Coverage.COVERAGE_MORE_0, '100.0')
			coverage.add_coverage('MP', Coverage.COVERAGE_MORE_9, '100.0')
			coverage.add_coverage('PA', Coverage.COVERAGE_ALL, '527.8')
			coverage.add_coverage('PA', Coverage.COVERAGE_MORE_0, '100.0')
			coverage.add_coverage('PA', Coverage.COVERAGE_MORE_9, '100.0')
			coverage.add_coverage('HA', Coverage.COVERAGE_ALL, '1449.8')
			coverage.add_coverage('HA', Coverage.COVERAGE_MORE_0, '100.0')
			coverage.add_coverage('HA', Coverage.COVERAGE_MORE_9, '100.0')
			coverage.add_coverage('PB1', Coverage.COVERAGE_ALL, '618.8')
			coverage.add_coverage('PB1', Coverage.COVERAGE_MORE_0, '100.0')
			if (sample_name == 'EVA003_S91'):
				coverage.add_coverage('PB1', Coverage.COVERAGE_MORE_9, '90.0')
				b_at_least_one_less_than_100 = True
			else: coverage.add_coverage('PB1', Coverage.COVERAGE_MORE_9, '100.0')
			coverage.add_coverage('NP', Coverage.COVERAGE_ALL, '439.8')
			coverage.add_coverage('NP', Coverage.COVERAGE_MORE_0, '100.0')
			coverage.add_coverage('NP', Coverage.COVERAGE_MORE_9, '100.0')
			coverage.add_coverage('NS', Coverage.COVERAGE_ALL, '1214.8')
			coverage.add_coverage('NS', Coverage.COVERAGE_MORE_0, '100.0')
			coverage.add_coverage('NS', Coverage.COVERAGE_MORE_9, '100.0')
			coverage.add_coverage('PB2', Coverage.COVERAGE_ALL, '690.8')
			coverage.add_coverage('PB2', Coverage.COVERAGE_MORE_0, '100.0')
			coverage.add_coverage('PB2', Coverage.COVERAGE_MORE_9, '100.0')
			coverage.add_coverage('NA', Coverage.COVERAGE_ALL, '1092.8')
			coverage.add_coverage('NA', Coverage.COVERAGE_MORE_0, '100.0')
			coverage.add_coverage('NA', Coverage.COVERAGE_MORE_9, '100.0')
			meta_value = manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success, coverage.to_json())

		### test if it has the less 100 in EVA003_S91 sample
		self.assertTrue(b_at_least_one_less_than_100)
		self.assertEquals(len(vect_files), count)
		
## print all coverage values
# 		try:
# 			project = Project.objects.get(name=project_name)
# 			for project_sample in project.project_samples.all():
# 				if (not project_sample.get_is_ready_to_proccess()): continue
# 				meta_value = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
# 				decode_coverage = DecodeObjects()
# 				coverage = decode_coverage.decode_result(meta_value.description)
# 				print(project_sample.sample.name)
# 				print(coverage)
# 		except User.DoesNotExist:
# 			self.fail("there's no project name: " + project_name)
		#########
		
		## copy all the data for getattr(settings, "MEDIA_ROOT", None)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))
		path_destination = os.path.join(getattr(settings, "MEDIA_ROOT", None), ConstantsTestsCase.DIR_PROJECTS)
		self.utils.make_path(os.path.join(getattr(settings, "MEDIA_ROOT", None), ConstantsTestsCase.DIR_PROJECTS))
		cmd = 'cp -r {}/{}/* {}'.format(self.baseDirectory, ConstantsTestsCase.DIR_PROJECTS, path_destination)
		os.system(cmd)
		
		create_tree = CreateTree()
		create_tree.create_tree_and_alignments(project, user)
		
		#### test all type of files for global sequences
		for type_files in project.vect_clean_file:
			self.assertTrue(os.path.exists(project.get_global_file_by_project(TypePath.MEDIA_ROOT, type_files)))
		
		decode_result = DecodeObjects()
		software_names = SoftwareNames()
		meta_sample = manageDatabase.get_project_metakey(project, MetaKeyAndValue.META_KEY_Run_Tree_All_Sequences, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_sample != None)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, meta_sample.value)
		result = decode_result.decode_result(meta_sample.description)
#		self.assertEquals('Mauve-2.4.0, Feb 13 2015', result.get_software(software_names.get_mauve_name()))
		self.assertEquals('Mafft-7.453; (--preservecase --quiet)', result.get_software(software_names.get_mafft_name()))
		self.assertEquals(3, result.get_number_softwares())
		
		meta_sample = manageDatabase.get_project_metakey(project, MetaKeyAndValue.META_KEY_Tree_Count_All_Sequences, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_sample != None)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, meta_sample.value)
		self.assertEquals('4', meta_sample.description)
		self.assertTrue(95, os.path.getsize(project.get_global_file_by_project(TypePath.MEDIA_ROOT, project.PROJECT_FILE_NAME_FASTTREE)))
		
		### get all elements and gene names
		for sequence_name in self.utils.get_elements_from_db(reference, user):
			for type_files in project.vect_clean_file:
				self.assertTrue(os.path.exists(project.get_global_file_by_element(TypePath.MEDIA_ROOT, sequence_name, type_files)))
			self.assertTrue(os.path.getsize(project.get_global_file_by_element(TypePath.MEDIA_ROOT, sequence_name, project.PROJECT_FILE_NAME_FASTTREE)) > 98)
			self.assertTrue(os.path.getsize(project.get_global_file_by_element(TypePath.MEDIA_ROOT, sequence_name, project.PROJECT_FILE_NAME_FASTA)) > 100)

			meta_key = metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_Run_Tree_By_Element, sequence_name)
			meta_sample = manageDatabase.get_project_metakey(project, meta_key, MetaKeyAndValue.META_VALUE_Success)
			self.assertTrue(meta_sample != None)
			self.assertEquals(MetaKeyAndValue.META_VALUE_Success, meta_sample.value)
			result = decode_result.decode_result(meta_sample.description)
#			self.assertEquals('Mauve-2.4.0, Feb 13 2015', result.get_software(software_names.get_mauve_name()))
			self.assertEquals('Mafft-7.453; (--preservecase --quiet)', result.get_software(software_names.get_mafft_name()))
			self.assertEquals('seqret (EMBOSS)-6.6.0.0; (-sformat fasta -osformat2 nexusnon)', result.get_software(software_names.get_seqret_name()))
			self.assertEquals(3, result.get_number_softwares())
		
			meta_key = metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_Tree_Count_By_Element, sequence_name)
			meta_sample = manageDatabase.get_project_metakey(project, meta_key, MetaKeyAndValue.META_VALUE_Success)
			self.assertTrue(meta_sample != None)
			self.assertEquals(MetaKeyAndValue.META_VALUE_Success, meta_sample.value)
			if (sequence_name == 'PB1'): self.assertEquals('4', meta_sample.description)
			else: self.assertEquals('5', meta_sample.description)
		
		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))


