'''
Created on Nov 21, 2017

@author: mmp
'''
import unittest
from utils.coverage import DrawCoverage, DrawAllCoverage
import os, filecmp
from django.conf import settings 
from constants.constantsTestsCase import ConstantsTestsCase
from utils.utils import Utils
from django.test.utils import override_settings
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference
from constants.constants import FileType, TypePath, FileExtensions
from constants.meta_key_and_values import MetaKeyAndValue
from managing_files.manage_database import ManageDatabase
from utils.parse_coverage_file import GetCoverage
from constants.software_names import SoftwareNames
from utils.result import Gene

class Test(unittest.TestCase):

	utils = Utils()
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass
	
	def test_coverage(self):
		coverage = DrawCoverage(None)
		vect_coverage = [1,2,4,5,6,7]
		vect_coverage.extend([2] * 100)
		vect_coverage.extend([20] * 100)
		vect_coverage.extend([200] * 100)
		vect_coverage.extend([500] * 100)
		vect_coverage.extend([300] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([3] * 100)
		vect_coverage.extend([33] * 100)
		vect_genes = [Gene('xpto', 10, 300, 1), Gene('xpto2', 500, 800, -1), Gene('xpto4', 700, 1000, 1)]
		var_more_90 = [20, 50, 700]
		var_more_50 = [2, 100, 500, 800]
		var_less_50 = [6, 200, 350]
		
		output_image = self.utils.get_temp_file("image_test_coverage", ".png")
		average_coverage = "23432"
		ratio_more_zero = "100"
		rati_more_nine = "98"
		sample_name = "xpto"
		sequence_name = "1"
		coverage.create_coverage(vect_coverage, vect_genes, var_more_90, var_more_50, var_less_50, output_image,
					average_coverage, ratio_more_zero, rati_more_nine, sample_name, sequence_name)

		self.assertTrue(os.path.exists(output_image))
		self.assertTrue(filecmp.cmp(output_image, os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_IMAGES, 'first.png' )))
		os.unlink(output_image)
		
		
	def test_coverage_2(self):
		coverage = DrawCoverage(None)
		vect_coverage = [1,2,4,5,6,7]
		vect_coverage.extend([2] * 100)
		vect_coverage.extend([20] * 100)
		vect_coverage.extend([200] * 100)
		vect_coverage.extend([500] * 100)
		vect_coverage.extend([300] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([20] * 100)
		vect_coverage.extend([200] * 100)
		vect_coverage.extend([500] * 100)
		vect_coverage.extend([300] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([20] * 100)
		vect_coverage.extend([200] * 100)
		vect_coverage.extend([500] * 100)
		vect_coverage.extend([300] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([3] * 100)
		vect_coverage.extend([33] * 100)
		vect_genes = [Gene('xpto', 10, 300, 1), Gene('xpto2', 500, 800, -1)]
		var_more_90 = [20, 50, 700]
		var_more_50 = [2, 100, 500, 800]
		var_less_50 = [6, 200, 350]
		
		output_image = self.utils.get_temp_file("image_test_coverage", ".png")
		average_coverage = "23432"
		ratio_more_zero = "100"
		rati_more_nine = "98"
		sample_name = "xpto"
		sequence_name = "1"
		coverage.create_coverage(vect_coverage, vect_genes, var_more_90, var_more_50, var_less_50, output_image,
					average_coverage, ratio_more_zero, rati_more_nine, sample_name, sequence_name)

		self.assertTrue(os.path.exists(output_image))
		self.assertTrue(filecmp.cmp(output_image, os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_IMAGES, 'second.png' )))
		os.unlink(output_image)

	def test_coverage_3(self):
		coverage = DrawCoverage(20)
		vect_coverage = [1,2,4,5,6,7]
		vect_coverage.extend([2] * 100)
		vect_coverage.extend([20] * 100)
		vect_coverage.extend([200] * 100)
		vect_coverage.extend([500] * 100)
		vect_coverage.extend([300] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([20] * 100)
		vect_coverage.extend([200] * 100)
		vect_coverage.extend([500] * 100)
		vect_coverage.extend([300] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([20] * 100)
		vect_coverage.extend([200] * 100)
		vect_coverage.extend([500] * 100)
		vect_coverage.extend([300] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([3] * 100)
		vect_coverage.extend([33] * 100)
		vect_genes = [Gene('xpto', 10, 300, 1), Gene('xpto2', 500, 800, -1)]
		var_more_90 = [20, 50, 700]
		var_more_50 = [2, 100, 500, 800]
		var_less_50 = [6, 200, 350]
		
		output_image = self.utils.get_temp_file("image_test_coverage", ".png")
		average_coverage = "23432"
		ratio_more_zero = "100"
		rati_more_nine = "98"
		sample_name = "xpto"
		sequence_name = "1"
		coverage.create_coverage(vect_coverage, vect_genes, var_more_90, var_more_50, var_less_50, output_image,
					average_coverage, ratio_more_zero, rati_more_nine, sample_name, sequence_name)

		self.assertTrue(os.path.exists(output_image))
		self.assertTrue(filecmp.cmp(output_image, os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_IMAGES, 'third.png' )))
		os.unlink(output_image)


	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_coverage_all(self):

		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))

		tab_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "resutl_vcf_to_tab2.tab")
		coverage_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, "A_H3N2_A_Hong_Kong_4801_2014.depth.gz")
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		
		self.assertTrue(os.path.exists(coverage_file))
		self.assertTrue(os.path.exists(tab_file))
		self.assertTrue(os.path.exists(gb_file))
		self.assertTrue(os.path.exists(fasta_file))

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
			
		sample_name = "run_snippy2_1"
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
			sample.save()

		project_name = "file_name_3"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = sample_name
			project.reference = reference
			project.owner = user
			project.save()
		
		## create project_sample
		project_sample = ProjectSample()
		project_sample.sample = sample
		project_sample.project = project
		project_sample.save()

		tab_file_from_freebayes = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)
		self.utils.copy_file(tab_file, tab_file_from_freebayes)
		coverage_file_test = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ, SoftwareNames.SOFTWARE_SNIPPY_name)
		self.utils.copy_file(coverage_file, coverage_file_test)

		### set coverage statistics
		manageDatabase = ManageDatabase()
		get_coverage = GetCoverage()
		coverage = get_coverage.get_coverage(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ,\
						SoftwareNames.SOFTWARE_SNIPPY_name), project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT))
		meta_sample = manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage,\
								MetaKeyAndValue.META_VALUE_Success, coverage.to_json())

		draw_all_coverage = DrawAllCoverage()
		draw_all_coverage.draw_all_coverages(project_sample)
		self.assertEqual(['HA', 'MP', 'NA', 'NP', 'NS', 'PA', 'PB1', 'PB2'], self.utils.get_elements_from_db(reference, user))
		self.assertEqual(['HA', 'MP', 'NA', 'NP', 'NS', 'PA', 'PB1', 'PB2'], self.utils.get_elements_with_CDS_from_db(reference, user))
		for gene in self.utils.get_elements_from_db(reference, user):
			output_image = project_sample.get_global_file_by_element(TypePath.MEDIA_ROOT, ProjectSample.PREFIX_FILE_COVERAGE, gene, FileExtensions.FILE_PNG)
			self.assertTrue(os.path.exists(output_image))
			self.assertTrue(filecmp.cmp(output_image, os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_IMAGES, 'coverage_{}.png'.format(gene) )))
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))


	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_coverage_all_1(self):

		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))

		tab_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_COVERAGE, "MEA_420468.tab")
		coverage_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_COVERAGE, "MEA_420468.depth.gz")
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_COVERAGE, "Measles_final.gbk")
		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_COVERAGE, "Measles_final.fasta")
		
		self.assertTrue(os.path.exists(coverage_file))
		self.assertTrue(os.path.exists(tab_file))
		self.assertTrue(os.path.exists(gb_file))
		self.assertTrue(os.path.exists(fasta_file))

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		ref_name = "second_stage_2_new"
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
			
		sample_name = "run_snippy2_1_12"
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
			sample.save()

		project_name = "file_name_3_ds"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = sample_name
			project.reference = reference
			project.owner = user
			project.save()
		
		## create project_sample
		project_sample = ProjectSample()
		project_sample.sample = sample
		project_sample.project = project
		project_sample.save()

		tab_file_from_freebayes = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)
		self.utils.copy_file(tab_file, tab_file_from_freebayes)
		coverage_file_test = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ, SoftwareNames.SOFTWARE_SNIPPY_name)
		self.utils.copy_file(coverage_file, coverage_file_test)

		### set coverage statistics
		manageDatabase = ManageDatabase()
		get_coverage = GetCoverage()
		coverage = get_coverage.get_coverage(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ,\
						SoftwareNames.SOFTWARE_SNIPPY_name), project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT))
		meta_sample = manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage,\
								MetaKeyAndValue.META_VALUE_Success, coverage.to_json())

		draw_all_coverage = DrawAllCoverage()
		draw_all_coverage.draw_all_coverages(project_sample)
		for gene in self.utils.get_elements_from_db(reference, user):
			output_image = project_sample.get_global_file_by_element(TypePath.MEDIA_ROOT, ProjectSample.PREFIX_FILE_COVERAGE, gene, FileExtensions.FILE_PNG)
			self.assertTrue(os.path.exists(output_image))
			self.assertTrue(filecmp.cmp(output_image, os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_IMAGES, 'coverage_{}.png'.format(gene) )))
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))


