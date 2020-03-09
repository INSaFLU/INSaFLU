'''
Created on Nov 27, 2017

@author: mmp
'''
import unittest, os, filecmp
from django.conf import settings 
from utils.utils import Utils
from constants.constantsTestsCase import ConstantsTestsCase
from utils.result import CountHits 
from managing_files.manage_database import ManageDatabase
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference, TagName, TagNames, CountVariations
from constants.meta_key_and_values import MetaKeyAndValue
from utils.collect_extra_data import CollectExtraData
from django.test.utils import override_settings
from utils.parse_coverage_file import GetCoverage
from constants.constants import TypePath, FileType, Constants
from constants.software_names import SoftwareNames

class Test(unittest.TestCase):

	constants_tests_case = ConstantsTestsCase()
	utils = Utils()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass

	def tearDown(self):
		pass


	def test_create_graph_minor_variants(self):
		manageDatabase = ManageDatabase()
	
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		expect_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_IMAGES, ConstantsTestsCase.MANAGING_FILE_GRAPH_VAR_HTML)
		expect_file_png = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_IMAGES, ConstantsTestsCase.MANAGING_FILE_GRAPH_VAR_PNG)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		ref_name = "second_stage_2"
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
		
		project_name = "several_names"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = project_name
			project.reference = reference
			project.owner = user
			project.save()
		
		## get all fastq files
		vect_files = self.constants_tests_case.get_all_fastq_files(self.baseDirectory)
		
		count = 1
		for vect_file in vect_files:
			sample_name = "_".join(os.path.basename(vect_file[0]).split('_')[0:2])
			try:
				sample = Sample.objects.get(name=sample_name)
			except Sample.DoesNotExist:
				sample = Sample()
				sample.name = sample_name
				sample.is_valid_1 = True
				sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
				sample.is_valid_2 = False
				sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
				sample.owner = user
				
				sample.is_ready_for_projects = True
				sample.is_obsolete = False
				sample.save()

			## create project_sample
			project_sample = ProjectSample()
			project_sample.sample = sample
			project_sample.project = project
			project_sample.is_finished = True
			project_sample.is_deleted = False
			project_sample.is_error = False
			project_sample.save()
			
			meta_value = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Count_Hits, MetaKeyAndValue.META_VALUE_Success)
			if (meta_value == None):
				count_hits = CountHits()
				count_hits.set_hits_50_90(count*10)
				count_hits.set_hits_less_50(count*8)
				manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Count_Hits, MetaKeyAndValue.META_VALUE_Success, count_hits.to_json())
			count += 1

		## test group tag_names
		self.assertTrue(project_sample.sample.get_tag_names() == None)
		
		### test the graph
		collectExtraData = CollectExtraData()
		(file_out_html, file_out_png) = collectExtraData.create_graph_minor_variants(project, user)
		self.assertTrue(os.path.exists(file_out_html))
		self.assertTrue(file_out_png == None)
		self.assertTrue(os.path.getsize(file_out_html), os.path.getsize(expect_file))
		# self.assertTrue(filecmp.cmp(file_out_png, expect_file_png))

		meta_project = manageDatabase.get_project_metakey(project, MetaKeyAndValue.META_KEY_Count_Samples_Var_Graph, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_project != None)
		self.assertEquals('4', meta_project.description)
		os.unlink(file_out_html)
		if (file_out_png != None): os.unlink(file_out_png)


	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_create_tree_and_alignments(self):
		"""
 		test global method
 		"""
		software_names = SoftwareNames()
		manageDatabase = ManageDatabase()
		
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT_TEST", None))
		path_destination = os.path.join(getattr(settings, "MEDIA_ROOT_TEST", None), ConstantsTestsCase.DIR_PROJECTS)
		self.utils.make_path(os.path.join(getattr(settings, "MEDIA_ROOT_TEST", None), ConstantsTestsCase.DIR_PROJECTS))
		cmd = 'cp -r {}/{}/* {}'.format(self.baseDirectory, ConstantsTestsCase.DIR_PROJECTS, path_destination)
		os.system(cmd)
		
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		expected_file_coverage = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_coverage_output.tsv")
		
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
		
		tag_name_name = 'xpto'
		try:
			tag_name = TagName.objects.get(name='xpto')
		except TagName.DoesNotExist as e:
			tag_name = TagName()
			tag_name.owner = user
			tag_name.name = tag_name_name
			tag_name.is_meta_data = False
			tag_name.save()
		
		get_coverage = GetCoverage()
		ProjectSample.objects.all().delete()
		Sample.objects.all().delete()
		temp_dir = self.utils.get_temp_dir()
		count = 0
		for vect_file in vect_files:
			self.utils.copy_file(vect_file[0], os.path.join(temp_dir, os.path.basename(vect_file[0])))
			self.utils.copy_file(vect_file[1], os.path.join(temp_dir, os.path.basename(vect_file[1])))
			
			n_id = 5000 + count + 1
			sample_name = "_".join(os.path.basename(vect_file[0]).split('_')[0:2])
			try:
				
				sample = Sample.objects.get(pk = n_id)
			except Sample.DoesNotExist as e:
				sample = Sample()
				sample.id = n_id
				sample.name = sample_name
				sample.is_valid_1 = True
				sample.file_name_1 = os.path.basename(vect_file[0])
				sample.path_name_1.name = os.path.join(temp_dir, os.path.basename(vect_file[0]))
				sample.is_valid_2 = False
				sample.file_name_2 = os.path.basename(vect_file[1])
				sample.path_name_2.name = os.path.join(temp_dir, os.path.basename(vect_file[1]))
				sample.owner = user
				sample.is_ready_for_projects = True
				sample.is_obsolete = False
				sample.type_subtype = 'xpto, zpto'
				sample.save()

			## add tag names to sample
			tag_names = TagNames()
			tag_names.value = tag_name_name + " " + tag_name_name
			tag_names.tag_name = tag_name
			tag_names.sample = sample
			tag_names.save()
			
			n_id = 5000 + count + 1
			try:
				project_sample = ProjectSample.objects.get(pk = n_id)
			except ProjectSample.DoesNotExist as e:
				## create project_sample
				project_sample = ProjectSample()
				project_sample.id = n_id
				project_sample.sample = sample
				project_sample.project = project
				project_sample.is_finished = True
				project_sample.is_deleted = False
				project_sample.is_error = False
			
				count_variations = CountVariations()
				count_variations.var_bigger_50_90 = n_id + 1
				count_variations.var_bigger_90 = n_id + 12
				count_variations.var_less_50 = 12 + count
				count_variations.total = n_id + 1 + n_id + 12
				count_variations.save()
				project_sample.count_variations = count_variations
				project_sample.save()
				
			coverage = get_coverage.get_coverage(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ,\
						software_names.get_snippy_name()), project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT))
			meta_value = manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success, coverage.to_json())
			count += 1
		
		### test group set tag names
		self.assertTrue(project_sample.sample.get_tag_names().count() == 1) 
		
		collect_extra_data = CollectExtraData();
		out_file = collect_extra_data.create_coverage_file(project, user)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_coverage))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### samples test
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output.csv")
		out_file = collect_extra_data.collect_sample_table(project, Constants.SEPARATOR_COMMA)
		
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output.tsv")
		out_file = collect_extra_data.collect_sample_table(project, Constants.SEPARATOR_TAB)
		
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### collect variations from snippy
		out_file = collect_extra_data.collect_variations_snippy(project)
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_variations_snippy.tsv")
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### collect variations from freebayes
		out_file = collect_extra_data.collect_variations_freebayes(project)
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_variations_freebayes.tsv")
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### count variations in project
		out_file = collect_extra_data.calculate_count_variations(project)
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_variations_total.tsv")
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### json file
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output.json")
		out_file = collect_extra_data.create_json_file_from_sample_csv(project)
		
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(os.path.exists(expected_file_samples))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### get sample result file
		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))


