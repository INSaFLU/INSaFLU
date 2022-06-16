'''
Created on Nov 27, 2017

@author: mmp
'''
import unittest, os, filecmp
from datetime import datetime
from django.conf import settings 
from utils.utils import Utils
from constants.constantsTestsCase import ConstantsTestsCase
from utils.result import CountHits
from managing_files.manage_database import ManageDatabase
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference, TagName, TagNames, CountVariations
from constants.meta_key_and_values import MetaKeyAndValue
from utils.collect_extra_data import CollectExtraData, ParsePangolinResult
from django.test.utils import override_settings
from utils.parse_coverage_file import GetCoverage
from constants.constants import TypePath, FileType, Constants, FileExtensions
from constants.software_names import SoftwareNames
from utils.result import Result, SoftwareDesc, KeyValues, KeyValue, GeneticElement, Gene, MaskingConsensus

class Test(unittest.TestCase):

	constants_tests_case = ConstantsTestsCase()
	utils = Utils()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass

	def tearDown(self):
		pass

	def test_pangolin_file(self):
		"""
		'taxon,lineage,conflict,pangoLEARN_version,status,note',
		'MN908947 SARSCoVDec200153,B.1.177,1.0,2021-04-01,passed_qc,'
		'MN908947 SARSCoVDec200234,B.1.1.7,1.0,2021-04-01,passed_qc,17/17 B.1.1.7 SNPs'
		"""
		
		pangolin_results = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "pangolin_results.csv")
		parse_pangolin_result = ParsePangolinResult(pangolin_results)
		self.assertTrue(parse_pangolin_result.has_data())
		self.assertEqual("B.1.177", parse_pangolin_result.get_value("MN908947 SARSCoVDec200153", ParsePangolinResult.KEY_LINEAGE))
		self.assertEqual("", parse_pangolin_result.get_value("MN908947 SARSCoVDec200153", ParsePangolinResult.KEY_SCORPIO))
		self.assertEqual("", parse_pangolin_result.get_value("MN908947 SARSCoVDec200153", "xpto"))
		self.assertEqual("B.1.1.7", parse_pangolin_result.get_value("MN908947 SARSCoVDec200234", ParsePangolinResult.KEY_LINEAGE))
		self.assertEqual("Alpha (B.1.1.7-like)", parse_pangolin_result.get_value("MN908947 SARSCoVDec200234", ParsePangolinResult.KEY_SCORPIO))
		self.assertEqual("B.1.177;B.1.1.7", parse_pangolin_result.get_value("MN908947 SARSCoVDec200", ParsePangolinResult.KEY_LINEAGE))
		self.assertEqual("Alpha (B.1.1.7-like)", parse_pangolin_result.get_value("MN908947 SARSCoVDec20", ParsePangolinResult.KEY_SCORPIO))
		self.assertEqual("", parse_pangolin_result.get_value("MN908947_SARSCoVDec200___", ParsePangolinResult.KEY_LINEAGE))
		self.assertEqual("", parse_pangolin_result.get_value("MN908947_SARSCo200___", ParsePangolinResult.KEY_SCORPIO))
		self.assertEqual("", parse_pangolin_result.get_value(None, ParsePangolinResult.KEY_LINEAGE))
		self.assertEqual("", parse_pangolin_result.get_value(None, ParsePangolinResult.KEY_SCORPIO))
		self.assertEqual("", parse_pangolin_result.get_value(None, "xpto"))
		
		parse_pangolin_result = ParsePangolinResult("xpto.xpt")
		self.assertFalse(parse_pangolin_result.has_data())
		

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
				sample.type_of_fastq = Sample.TYPE_OF_FASTQ_illumina
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
		manage_database = ManageDatabase()
		
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
				sample.type_of_fastq = Sample.TYPE_OF_FASTQ_illumina
				sample.type_subtype = 'xpto, zpto'
				sample.save()

				## KEY META_KEY_Identify_Sample_Software and META_KEY_Fastq_Trimmomatic_Software
				parameters = "let's go again";
				result_all_2 = Result()
				result_all_2.add_software(SoftwareDesc(software_names.get_fastq_name(), 
						software_names.get_fastq_version(), parameters))
				manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Identify_Sample_Software,
					MetaKeyAndValue.META_VALUE_Success, result_all_2.to_json())
				
				result_all_2.add_software(SoftwareDesc(software_names.get_trimmomatic_name(), 
						software_names.get_trimmomatic_version(), parameters))
				result_all_2.add_software(SoftwareDesc(software_names.get_trimmomatic_name(), 
						software_names.get_trimmomatic_version(), parameters + "second version trimmomatic"))
				manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software,
					MetaKeyAndValue.META_VALUE_Success, result_all_2.to_json())
			## add tag names to sample
			tag_names = TagNames()
			tag_names.value = tag_name_name + " " + tag_name_name
			tag_names.tag_name = tag_name
			tag_names.sample = sample
			tag_names.save()
			
			result_all = Result()
			n_id = 5000 + count + 1
			try:
				project_sample = ProjectSample.objects.get(pk = n_id)
			except ProjectSample.DoesNotExist as e:
				## create project_sample
				project_sample = ProjectSample()
				project_sample.id = n_id
				project_sample.sample = sample
				project_sample.seq_name_all_consensus = project_sample.sample.name
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
				
				### KEY META_KEY_Snippy_Freebayes
				parameters = "let's go";
				result_all.add_software(SoftwareDesc(software_names.get_spades_name(), 
						software_names.get_spades_version(), parameters))
				parameters = "let's go";
				result_all.add_software(SoftwareDesc(software_names.get_snippy_name(), 
						software_names.get_snippy_version(), parameters + "_2222"))
				### don't save this one to stay empty in the file
				if (not vect_file[0].endswith("EVA011_S54_L001_R1_001.fastq.gz")):
					manage_database.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Snippy_Freebayes, 
						MetaKeyAndValue.META_VALUE_Success, result_all.to_json())
				
			coverage = get_coverage.get_coverage(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ,\
						software_names.get_snippy_name()), project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT))
			meta_value = manage_database.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage,
					MetaKeyAndValue.META_VALUE_Success, coverage.to_json())
			count += 1
		
		### masking data
		manageDatabase = ManageDatabase()
		geneticElement = GeneticElement()
		geneticElement.add_gene('element_name', 100, Gene('name', 12, 45, 1))
		masking_consensus = MaskingConsensus()
		masking_consensus.set_mask_sites("2,5,-5,5,cf")
		masking_consensus.set_mask_from_beginning("2")
		masking_consensus.set_mask_from_ends("2")
		masking_consensus.set_mask_regions("[2-4],[4-30],[1-2], [500-320], [50-32]")
		geneticElement.dt_elements_mask['element_name'] = masking_consensus
		manageDatabase.set_project_metakey(project, project.owner, MetaKeyAndValue.META_KEY_Masking_consensus,
			MetaKeyAndValue.META_VALUE_Success, geneticElement.to_json())
		
		### test group set tag names
		self.assertTrue(project_sample.sample.get_tag_names().count() == 1) 
		
		collect_extra_data = CollectExtraData();
		out_file = collect_extra_data.create_coverage_file(project, user)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_coverage))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### samples test
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output.csv")
		out_file = collect_extra_data.collect_sample_table(project, Constants.SEPARATOR_COMMA, CollectExtraData.SAMPLE_LIST_list)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		self.utils.copy_file(out_file,
			project.get_global_file_by_project(TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV))
		if (os.path.exists(out_file)): os.unlink(out_file)

		### samples test
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output_settings.csv")
		out_file = collect_extra_data.collect_sample_table(project, Constants.SEPARATOR_COMMA, CollectExtraData.SAMPLE_LIST_list_settings)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		self.utils.copy_file(out_file,
			project.get_global_file_by_project(TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_SAMPLE_RESULT_SETTINGS_CSV))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### samples test CSV simple
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output_simple.csv")
		out_file = collect_extra_data.collect_sample_table(project, Constants.SEPARATOR_COMMA, CollectExtraData.SAMPLE_LIST_simple)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		self.utils.copy_file(out_file, 
			project.get_global_file_by_project(TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV_simple))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### samples test TSV
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output.tsv")
		out_file = collect_extra_data.collect_sample_table(project, Constants.SEPARATOR_TAB, CollectExtraData.SAMPLE_LIST_list)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### collect variations from snippy
		out_file = collect_extra_data.collect_variations_snippy(project)
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_variations_snippy.tsv")
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### collect variations from freebayes
		vect_type_remove = ['ins', 'del']
		out_file = collect_extra_data.collect_variations_freebayes(project, vect_type_remove)
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_variations_freebayes.tsv")
		self.utils.copy_file(out_file,
			project.get_global_file_by_project(TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_TAB_VARIATIONS_FREEBAYES))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### count variations in project
		out_file = collect_extra_data.calculate_count_variations(project)
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_variations_total.tsv")
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### json file
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output_1.json")
		out_file = collect_extra_data.create_json_file_from_sample_csv(project)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(os.path.exists(expected_file_samples))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### zip file, zip files change every time that are produced
		collect_extra_data.zip_several_files(project)
		out_file = project.get_global_file_by_project(TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_all_files_zipped)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(os.path.getsize(out_file) > 1616)
		
		### masking data
		manageDatabase = ManageDatabase()
		geneticElement = GeneticElement()
		geneticElement.add_gene('PB2', 100, Gene('name', 12, 45, 1))
		masking_consensus = MaskingConsensus()
		masking_consensus.set_mask_sites("2,5,-5,5,cf")
		masking_consensus.set_mask_from_beginning("2")
		masking_consensus.set_mask_from_ends("3")
		geneticElement.dt_elements_mask['PB2'] = masking_consensus
		manageDatabase.set_project_metakey(project, project.owner, MetaKeyAndValue.META_KEY_Masking_consensus,
			MetaKeyAndValue.META_VALUE_Success, geneticElement.to_json())
		
		collect_extra_data.mask_all_consensus_files(project)
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "Consensus_EVA001_S66.fasta")
		self.assertTrue(os.path.exists(expected_file_samples))
		for project_sample in project.project_samples.all():
			self.assertTrue(filecmp.cmp(project_sample.get_consensus_file(TypePath.MEDIA_ROOT), expected_file_samples))
			self.assertTrue(os.path.exists(project_sample.get_backup_consensus_file()))
			self.assertFalse(filecmp.cmp(project_sample.get_consensus_file(TypePath.MEDIA_ROOT), project_sample.get_backup_consensus_file()))
		
		masking_consensus = MaskingConsensus()
		masking_consensus.set_mask_sites("2,5,-5,5,cf")
		masking_consensus.set_mask_from_beginning("10")
		geneticElement.dt_elements_mask['PB2'] = masking_consensus
		geneticElement.add_gene('MP', 100, Gene('MP', 12, 45, 1))
		masking_consensus = MaskingConsensus()
		masking_consensus.set_mask_regions("20-30,35-40")
		masking_consensus.set_mask_from_beginning("10")
		masking_consensus.set_mask_from_ends("8")
		geneticElement.dt_elements_mask['MP'] = masking_consensus
		manageDatabase.set_project_metakey(project, project.owner, MetaKeyAndValue.META_KEY_Masking_consensus,
			MetaKeyAndValue.META_VALUE_Success, geneticElement.to_json())
		collect_extra_data.mask_all_consensus_files(project)
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "Consensus_EVA001_S66_second_mask.fasta")
		self.assertTrue(os.path.exists(expected_file_samples))
		for project_sample in project.project_samples.all():
			self.assertTrue(filecmp.cmp(project_sample.get_consensus_file(TypePath.MEDIA_ROOT), expected_file_samples))
			self.assertTrue(os.path.exists(project_sample.get_backup_consensus_file()))
			self.assertFalse(filecmp.cmp(project_sample.get_consensus_file(TypePath.MEDIA_ROOT), project_sample.get_backup_consensus_file()))
		
		### get sample result file
		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))

	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_create_tree_and_alignments_2(self):
		"""
 		test global method
 		"""
		software_names = SoftwareNames()
		manage_database = ManageDatabase()
		
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
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME + '6000')
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME + '6000'
			user.id = 6000
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
		
		project_name = "several_names_test_create_tree_6000"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.id= 6000
			project.name = project_name
			project.reference = reference
			project.owner = user
			project.creation_date = datetime(2002, 12, 4, 20, 30, 40)
			project.last_change_date = datetime(2002, 12, 14, 20, 30, 40)
			project.save()
		
		## get all fastq files
		vect_files = self.constants_tests_case.get_all_fastq_files(self.baseDirectory)
		
		tag_name_name = 'xpto_1'
		try:
			tag_name = TagName.objects.get(name='xpto_1')
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
			
			n_id = 6000 + count + 1
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
				sample.type_of_fastq = Sample.TYPE_OF_FASTQ_illumina
				sample.is_obsolete = False
				sample.type_subtype = 'xpto, zpto'
				sample.save()

				## KEY META_KEY_Identify_Sample_Software and META_KEY_Fastq_Trimmomatic_Software
				parameters = "let's go again";
				result_all_2 = Result()
				result_all_2.add_software(SoftwareDesc(software_names.get_fastq_name(), 
						software_names.get_fastq_version(), parameters))
				manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Identify_Sample_Software,
					MetaKeyAndValue.META_VALUE_Success, result_all_2.to_json())
				
				### set statistics of number of reads
				key_values = KeyValues()
				for _, key in enumerate(SoftwareNames.SOFTWARE_ILLUMINA_stat_collect):
					key_values.add_key_value(KeyValue(key, str(_ + 10000)))
				result_all_2.add_software(SoftwareDesc(SoftwareNames.SOFTWARE_ILLUMINA_stat, "", "",
											key_values))
					
				key_values = None
				if (n_id < 6004):
					key_values = KeyValues()
					key_values.add_key_value(KeyValue("Input Read Pairs:", "10"))
					key_values.add_key_value(KeyValue("Both Surviving:", "5"))
					key_values.add_key_value(KeyValue("Forward Only Surviving:", "2"))
					if (n_id == 6001): key_values.add_key_value(KeyValue("Reverse Only Surviving:", "33"))
					key_values.add_key_value(KeyValue("Dropped:", "44"))
				result_all_2.add_software(SoftwareDesc(software_names.get_trimmomatic_name(), 
						software_names.get_trimmomatic_version(), parameters, key_values))
				result_all_2.add_software(SoftwareDesc(software_names.get_trimmomatic_name(), 
						software_names.get_trimmomatic_version(), parameters + "second version trimmomatic"))
				manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software,
					MetaKeyAndValue.META_VALUE_Success, result_all_2.to_json())
			## add tag names to sample
			tag_names = TagNames()
			tag_names.value = tag_name_name + " " + tag_name_name
			tag_names.tag_name = tag_name
			tag_names.sample = sample
			tag_names.save()
			
			result_all = Result()
			n_id = 6000 + count + 1
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
				if (n_id == 6001): project_sample.is_mask_consensus_sequences = False
				else: project_sample.is_mask_consensus_sequences = True
			
				count_variations = CountVariations()
				count_variations.var_bigger_50_90 = n_id + 1
				count_variations.var_bigger_90 = n_id + 12
				count_variations.var_less_50 = 12 + count
				count_variations.total = n_id + 1 + n_id + 12
				count_variations.save()
				project_sample.count_variations = count_variations
				project_sample.save()
				
				### KEY META_KEY_Snippy_Freebayes
				parameters = "let's go";
				result_all.add_software(SoftwareDesc(software_names.get_spades_name(), 
						software_names.get_spades_version(), parameters))
				parameters = "let's go";
				result_all.add_software(SoftwareDesc(software_names.get_snippy_name(), 
						software_names.get_snippy_version(), parameters + "_2222"))
				### don't save this one to stay empty in the file
				if (not vect_file[0].endswith("EVA011_S54_L001_R1_001.fastq.gz")):
					manage_database.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Snippy_Freebayes, 
						MetaKeyAndValue.META_VALUE_Success, result_all.to_json())
				
			coverage = get_coverage.get_coverage(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ,\
						software_names.get_snippy_name()), project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT))
			manage_database.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage,
					MetaKeyAndValue.META_VALUE_Success, coverage.to_json())
			count += 1
		
		### test group set tag names
		self.assertTrue(project_sample.sample.get_tag_names().count() == 1) 
		
		collect_extra_data = CollectExtraData();
		out_file = collect_extra_data.create_coverage_file(project, user)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_coverage))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### samples test
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output_2.csv")
		out_file = collect_extra_data.collect_sample_table(project, Constants.SEPARATOR_COMMA, CollectExtraData.SAMPLE_LIST_list)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)

		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output_2_settings.csv")
		out_file = collect_extra_data.collect_sample_table(project, Constants.SEPARATOR_COMMA, CollectExtraData.SAMPLE_LIST_list_settings)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output_2_simple.csv")
		out_file = collect_extra_data.collect_sample_table(project, Constants.SEPARATOR_COMMA, CollectExtraData.SAMPLE_LIST_simple)
		self.assertTrue(os.path.exists(out_file))
		self.utils.copy_file(out_file, 
			project.get_global_file_by_project(TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV_simple))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output_2.tsv")
		out_file = collect_extra_data.collect_sample_table(project, Constants.SEPARATOR_TAB, CollectExtraData.SAMPLE_LIST_list)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### collect variations from snippy
		out_file = collect_extra_data.collect_variations_snippy(project)
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_variations_snippy.tsv")
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### collect variations from freebayes
		vect_type_remove = ['ins', 'del']
		out_file = collect_extra_data.collect_variations_freebayes(project, vect_type_remove)
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_variations_freebayes.tsv")
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### count variations in project
		out_file = collect_extra_data.calculate_count_variations(project)
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_variations_total_2.tsv")
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### json file
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output.json")
		out_file = collect_extra_data.create_json_file_from_sample_csv(project)
		
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(os.path.exists(expected_file_samples))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		#################################3
		## test collect samples
		collect_extra_data.collect_sample_list(user, True)
		csv_file = self.utils.get_sample_list_by_user(user.id, "MEDIA_ROOT", FileExtensions.FILE_CSV)
		self.assertTrue(os.path.exists(csv_file))
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "AllSamples_1.csv")
		self.assertTrue(os.path.exists(expected_file_samples))
		self.assertTrue(filecmp.cmp(csv_file, expected_file_samples))
		
		tsv_file = self.utils.get_sample_list_by_user(user.id, "MEDIA_ROOT", FileExtensions.FILE_TSV)
		self.assertTrue(os.path.exists(tsv_file))
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "AllSamples_1.tsv")
		self.assertTrue(os.path.exists(expected_file_samples))
		self.assertTrue(filecmp.cmp(tsv_file, expected_file_samples))
		
		#########################################
		## test collect projects
		collect_extra_data.collect_project_list(user, True)
		csv_file = self.utils.get_project_list_by_user(user.id, "MEDIA_ROOT", FileExtensions.FILE_CSV)
		self.assertTrue(os.path.exists(csv_file))
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "AllProjects_1.csv")
		self.assertTrue(os.path.exists(expected_file_samples))
		self.assertTrue(filecmp.cmp(csv_file, expected_file_samples))
		
		tsv_file = self.utils.get_project_list_by_user(user.id, "MEDIA_ROOT", FileExtensions.FILE_TSV)
		self.assertTrue(os.path.exists(tsv_file))
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "AllProjects_1.tsv")
		self.assertTrue(os.path.exists(expected_file_samples))
		self.assertTrue(filecmp.cmp(tsv_file, expected_file_samples))
		
		### get sample result file
		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))
		
	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_create_tree_and_alignments_3(self):
		"""
 		test global method
 		mmp@cs-nb0008:~/git/INSaFLU$ more /tmp/tests_insa_flu/projects/result/user_6000/project_6000/main_result/proportions_iSNVs_graph.tsv
 		Sample	Less 50	Between 90 and 50
		EVA011_S54	15	6005
		EVA003_S91	14	6004
		mmp@cs-nb0008:~/git/INSaFLU$ more /home/mmp/git/INSaFLU/static_all/tests/global_project/insa_flu_variations_total_2_ont.tsv
		Sample	Less 50	Between 90 and 50
		EVA011_S54	15	6005
		EVA003_S91	14	6004
		EVA002_S52	13	6003
		EVA001_S66	12	6002
 		"""
		software_names = SoftwareNames()
		manage_database = ManageDatabase()
		
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT_TEST", None))
		path_destination = os.path.join(getattr(settings, "MEDIA_ROOT_TEST", None), ConstantsTestsCase.DIR_PROJECTS)
		self.utils.make_path(os.path.join(getattr(settings, "MEDIA_ROOT_TEST", None), ConstantsTestsCase.DIR_PROJECTS))
		cmd = 'cp -r {}/{}/* {}'.format(self.baseDirectory, ConstantsTestsCase.DIR_PROJECTS, path_destination)
		os.system(cmd)
		
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME + '6000')
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME + '6000'
			user.id = 6000
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
		
		project_name = "several_names_test_create_tree_6000"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.id= 6000
			project.name = project_name
			project.reference = reference
			project.owner = user
			project.creation_date = datetime(2002, 12, 4, 20, 30, 40)
			project.last_change_date = datetime(2002, 12, 14, 20, 30, 40)
			project.save()
		
		## get all fastq files
		vect_files = self.constants_tests_case.get_all_fastq_files(self.baseDirectory)
		
		tag_name_name = 'xpto_1'
		try:
			tag_name = TagName.objects.get(name='xpto_1')
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
		
		### illumina files
		for vect_file in vect_files:
			self.utils.copy_file(vect_file[0], os.path.join(temp_dir, os.path.basename(vect_file[0])))
			self.utils.copy_file(vect_file[1], os.path.join(temp_dir, os.path.basename(vect_file[1])))
			
			n_id = 6000 + count + 1
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
				if (n_id > 6002): sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina)
				sample.is_obsolete = False
				sample.type_subtype = 'xpto, zpto'
				sample.save()

				## KEY META_KEY_Identify_Sample_Software and META_KEY_Fastq_Trimmomatic_Software
				parameters = "let's go again";
				result_all_2 = Result()
				result_all_2.add_software(SoftwareDesc(software_names.get_fastq_name(), 
						software_names.get_fastq_version(), parameters))
				manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Identify_Sample_Software,
					MetaKeyAndValue.META_VALUE_Success, result_all_2.to_json())
				
				key_values = None
				if (n_id < 6004):
					key_values = KeyValues()
					key_values.add_key_value(KeyValue("Input Read Pairs:", "xpto"))
					key_values.add_key_value(KeyValue("Both Surviving:", "xpto1"))
					key_values.add_key_value(KeyValue("Forward Only Surviving:", "xpto2"))
					if (n_id == 6001): key_values.add_key_value(KeyValue("Reverse Only Surviving:", "xpto3"))
					key_values.add_key_value(KeyValue("Dropped:", "xpto4"))
				result_all_2.add_software(SoftwareDesc(software_names.get_trimmomatic_name(), 
						software_names.get_trimmomatic_version(), parameters, key_values))
				result_all_2.add_software(SoftwareDesc(software_names.get_trimmomatic_name(), 
						software_names.get_trimmomatic_version(), parameters + "second version trimmomatic"))
				manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software,
					MetaKeyAndValue.META_VALUE_Success, result_all_2.to_json())
			## add tag names to sample
			tag_names = TagNames()
			tag_names.value = tag_name_name + " " + tag_name_name
			tag_names.tag_name = tag_name
			tag_names.sample = sample
			tag_names.save()
			
			result_all = Result()
			n_id = 6000 + count + 1
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
				if (n_id == 6001): project_sample.is_mask_consensus_sequences = False
				else: project_sample.is_mask_consensus_sequences = True
			
				count_variations = CountVariations()
				count_variations.var_bigger_50_90 = n_id + 1
				count_variations.var_bigger_90 = n_id + 12
				count_variations.var_less_50 = 12 + count
				count_variations.total = n_id + 1 + n_id + 12
				count_variations.save()
				project_sample.count_variations = count_variations
				project_sample.save()
				
				### KEY META_KEY_Snippy_Freebayes
				parameters = "let's go";
				result_all.add_software(SoftwareDesc(software_names.get_spades_name(), 
						software_names.get_spades_version(), parameters))
				parameters = "let's go";
				result_all.add_software(SoftwareDesc(software_names.get_snippy_name(), 
						software_names.get_snippy_version(), parameters + "_2222"))
				### don't save this one to stay empty in the file
				if (not vect_file[0].endswith("EVA011_S54_L001_R1_001.fastq.gz")):
					manage_database.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Snippy_Freebayes, 
						MetaKeyAndValue.META_VALUE_Success, result_all.to_json())
				
			coverage = get_coverage.get_coverage(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ,\
						software_names.get_snippy_name()), project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT))
			meta_value = manage_database.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage,
					MetaKeyAndValue.META_VALUE_Success, coverage.to_json())
			count += 1
		
		### ONT files
		count_before = n_id
		for vect_file in vect_files:
			self.utils.copy_file(vect_file[0], os.path.join(temp_dir, os.path.basename(vect_file[0])))
			
			n_id = 6000 + count + 1
			sample_name = "_".join(os.path.basename(vect_file[0]).split('_')[0:2]) + "_ONT"
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
				sample.owner = user
				sample.is_ready_for_projects = True
				sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_ont)
				sample.is_obsolete = False
				sample.type_subtype = 'xpto, zpto'
				sample.save()

				## KEY META_KEY_Identify_Sample_Software and META_KEY_Fastq_Trimmomatic_Software
				parameters = "let's go again";
				result_all_2 = Result()

				key_values = None
				if (count_before < 6007):
					key_values = KeyValues()
					for _, key in enumerate(SoftwareNames.SOFTWARE_NANOSTAT_vect_info_to_collect):
						key_values.add_key_value(KeyValue(key, str(_)))
						if (n_id == 6005): key_values.add_key_value(KeyValue(key, str(_ * 2)))
					result_all_2.add_software(SoftwareDesc(software_names.get_NanoStat_name(), software_names.get_NanoStat_version(),
						software_names.get_NanoStat_parameters(), key_values))
					result_all_2.add_software(SoftwareDesc(software_names.get_rabbitQC_name(), software_names.get_rabbitQC_version(),
						software_names.get_rabbitQC_parameters()))
					result_all_2.add_software(SoftwareDesc(software_names.get_NanoFilt_name(), software_names.get_NanoFilt_version(), parameters))
					result_all_2.add_software(SoftwareDesc(software_names.get_NanoStat_name(), software_names.get_NanoStat_version(),
						software_names.get_NanoStat_parameters(), key_values))
				else:
					result_all_2.add_software(SoftwareDesc(software_names.get_NanoStat_name(), software_names.get_NanoStat_version(),
						software_names.get_NanoStat_parameters()))
					result_all_2.add_software(SoftwareDesc(software_names.get_rabbitQC_name(), software_names.get_rabbitQC_version(),
						software_names.get_rabbitQC_parameters()))
					result_all_2.add_software(SoftwareDesc(software_names.get_NanoFilt_name(), software_names.get_NanoFilt_version(), parameters))
					result_all_2.add_software(SoftwareDesc(software_names.get_NanoStat_name(), software_names.get_NanoStat_version(),
						software_names.get_NanoStat_parameters() + "second"))
				manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software,
					MetaKeyAndValue.META_VALUE_Success, result_all_2.to_json())
			## add tag names to sample
			tag_names = TagNames()
			tag_names.value = tag_name_name + " " + tag_name_name
			tag_names.tag_name = tag_name
			tag_names.sample = sample
			tag_names.save()
			
			result_all = Result()
			n_id = 6000 + count + 1
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
				if (n_id == 6001): project_sample.is_mask_consensus_sequences = False
				else: project_sample.is_mask_consensus_sequences = True
			
				count_variations = CountVariations()
				count_variations.var_bigger_50_90 = n_id + 1
				count_variations.var_bigger_90 = n_id + 12
				count_variations.var_less_50 = 12 + count
				count_variations.total = n_id + 1 + n_id + 12
				count_variations.save()
				project_sample.count_variations = count_variations
				project_sample.save()
				
				### KEY META_KEY_Snippy_Freebayes
				parameters = "let's go";
				parameters_depth = "-q 10"
				parameters_medaka_consensus = ""
				coverage_limit = n_id
				result_all.add_software(SoftwareDesc(software_names.get_medaka_name(),\
								software_names.get_medaka_version(), "consensus " + parameters_medaka_consensus))
				result_all.add_software(SoftwareDesc(software_names.get_samtools_name(),\
								software_names.get_samtools_version(), "depth -aa {}".format(parameters_depth)))
				result_all.add_software(SoftwareDesc(software_names.get_medaka_name(),\
								software_names.get_medaka_version(), "variant --verbose"))
				result_all.add_software(SoftwareDesc(software_names.get_insaflu_parameter_limit_coverage_name(),\
								"", "Threshold:{}".format(coverage_limit)))
				limit_to_mask_consensus = 60
				msa_parameters = "MSA parameters"
				result_all.add_software(SoftwareDesc(software_names.get_msa_masker_name(), software_names.get_msa_masker_version(),\
					"{}; for coverages less than {} in {}% of the regions.".format(msa_parameters,\
					coverage_limit,									
					100 - limit_to_mask_consensus) ))
				### don't save this one to stay empty in the file
				if (not vect_file[0].endswith("EVA011_S54_L001_R1_001.fastq.gz")):
					manage_database.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Medaka, 
						MetaKeyAndValue.META_VALUE_Success, result_all.to_json())
			meta_value = manage_database.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage,
					MetaKeyAndValue.META_VALUE_Success, coverage.to_json())
			count += 1

		### test group set tag names
		self.assertTrue(project_sample.sample.get_tag_names().count() == 1) 
		
		collect_extra_data = CollectExtraData();
		expected_file_coverage = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_coverage_output_ont.tsv")
		out_file = collect_extra_data.create_coverage_file(project, user)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_coverage))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### samples test
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output_2_ont.csv")
		out_file = collect_extra_data.collect_sample_table(project, Constants.SEPARATOR_COMMA, CollectExtraData.SAMPLE_LIST_list)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)

		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output_2_ont_settings.csv")
		out_file = collect_extra_data.collect_sample_table(project, Constants.SEPARATOR_COMMA, CollectExtraData.SAMPLE_LIST_list_settings)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### samples test
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output_2_ont_simple.csv")
		out_file = collect_extra_data.collect_sample_table(project, Constants.SEPARATOR_COMMA, CollectExtraData.SAMPLE_LIST_simple)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		self.utils.copy_file(out_file,
			project.get_global_file_by_project(TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV_simple))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output_2_ont.tsv")
		out_file = collect_extra_data.collect_sample_table(project, Constants.SEPARATOR_TAB, CollectExtraData.SAMPLE_LIST_list)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### collect variations from snippy
		out_file = collect_extra_data.collect_variations_snippy(project)
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_variations_snippy_ont.tsv")
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### collect variations from freebayes
		vect_type_remove = ['ins', 'del']
		out_file = collect_extra_data.collect_variations_freebayes(project, vect_type_remove)
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_variations_freebayes_ont.tsv")
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### count variations in project
		out_file = collect_extra_data.calculate_count_variations(project)
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_variations_total_2_ont.tsv")
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### json file
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output_ont.json")
		out_file = collect_extra_data.create_json_file_from_sample_csv(project)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(os.path.exists(expected_file_samples))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### get sample result file
		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))

