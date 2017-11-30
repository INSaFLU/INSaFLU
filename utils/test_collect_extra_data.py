'''
Created on Nov 27, 2017

@author: mmp
'''
import unittest, os, filecmp
from django.conf import settings 
from utils.software import Software
from utils.utils import Utils
from utils.constantsTestsCase import ConstantsTestsCase
from utils.result import CountHits 
from managing_files.manage_database import ManageDatabase
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference
from utils.meta_key_and_values import MetaKeyAndValue
from utils.collect_extra_data import CollectExtraData
from django.test.utils import override_settings

class Test(unittest.TestCase):

	constants_tests_case = ConstantsTestsCase()

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
				sample.is_rejected = False
				sample.is_valid_1 = True
				sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
				sample.is_valid_2 = False
				sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
				sample.owner = user
				
				sample.is_ready_for_projects = True
				sample.is_obsolete = False
				sample.is_rejected = False
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


# 	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
# 	def test_collect_extra_data(self):
# 
# 		manageDatabase = ManageDatabase()
# 		metaKeyAndValue = MetaKeyAndValue()	
# 		
# 		### set the taskID
# 		manageDatabase.set_project_metakey(project, user, MetaKeyAndValue.META_KEY_Queue_TaskID,\
# 								MetaKeyAndValue.META_VALUE_Queue, 'meta_sample.description')
# 
# 		
# 		lst_meta_sample = manageDatabase.get_project_metakey(project, MetaKeyAndValue.META_KEY_Queue_TaskID, None)
# 		self.assertEquals(2, len(lst_meta_sample))
# 		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, lst_meta_sample[0].value)
# 		self.assertEquals('meta_sample.description', lst_meta_sample[0].description)


