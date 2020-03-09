'''
Created on Jan 7, 2018

@author: mmp
'''
import unittest
from constants.constantsTestsCase import ConstantsTestsCase
from utils.utils import Utils
from utils.software import Software
from utils.result import Coverage
from django.conf import settings
from utils.parse_out_files import ParseOutFiles
from constants.constants import Constants, FileExtensions
import os, filecmp, csv

class Test(unittest.TestCase):

	constantsTestsCase = ConstantsTestsCase()
	utils = Utils()
	software = Software()

	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), self.constantsTestsCase.MANAGING_TESTS)


	def test_parse_snippy_tab_files(self):
		"""
		test ParseOutFiles for tab files
		"""
		parse_out_files = ParseOutFiles()
		tab_file_to_process = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "have_multiple_var.tab")
		expected_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "out_aggregate_tab_files.tab")

		temp_out_file = self.utils.get_temp_file('test_tab_out', '.tab')
		vect_type_out = ['snp']
		vect_type_remove = []
		b_add_header = True
		with open(temp_out_file, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=Constants.SEPARATOR_TAB, quotechar='"', quoting=csv.QUOTE_ALL)
			parse_out_files.parse_tab_files('sample_name', tab_file_to_process, csv_writer, vect_type_out, vect_type_remove, 101, b_add_header)
			
			vect_type_out = ['del']
			b_add_header = False
			parse_out_files.parse_tab_files('sample_name', tab_file_to_process, csv_writer, vect_type_out, vect_type_remove, 101, b_add_header)

		self.assertTrue(filecmp.cmp(expected_file, temp_out_file))
		if (os.path.exists(temp_out_file)): os.unlink(temp_out_file)


	def test_parse_snippy_tab_files_2(self):
		"""
		test ParseOutFiles for tab files
		"""
		parse_out_files = ParseOutFiles()
		tab_file_to_process = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "have_multiple_var.tab")
		expected_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "out_aggregate_tab_files_2.tab")

		temp_out_file = self.utils.get_temp_file('test_tab_out', '.tab')
		vect_type_out = ['snp', 'del']
		vect_type_remove = []
		b_add_header = True
		with open(temp_out_file, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=Constants.SEPARATOR_TAB, quotechar='"', quoting=csv.QUOTE_ALL)
			parse_out_files.parse_tab_files('sample_name', tab_file_to_process, csv_writer, vect_type_out, vect_type_remove, 101, b_add_header)
			
		self.assertTrue(os.path.exists(temp_out_file))
		self.assertTrue(filecmp.cmp(expected_file, temp_out_file))
		if (os.path.exists(temp_out_file)): os.unlink(temp_out_file)

		
	def test_parse_freebays_tab_files_2(self):
		"""
		test ParseOutFiles for tab files
		"""
		parse_out_files = ParseOutFiles()
		tab_file_to_process = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "run_freebayes_parallel.tab")
		expected_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "out_aggregate_tab_files_freebays.tab")

		temp_out_file = self.utils.get_temp_file('test_tab_out', '.tab')
		vect_type_out = ['snp', 'del']
		vect_type_remove = []
		b_add_header = True
		with open(temp_out_file, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=Constants.SEPARATOR_TAB, quotechar='"', quoting=csv.QUOTE_ALL)
			parse_out_files.parse_tab_files('sample_name', tab_file_to_process, csv_writer, vect_type_out, vect_type_remove, 50, b_add_header)
			
		self.assertTrue(os.path.exists(temp_out_file))
		self.assertTrue(filecmp.cmp(expected_file, temp_out_file))
		if (os.path.exists(temp_out_file)): os.unlink(temp_out_file)

	def test_parse_freebays_tab_files_3(self):
		"""
		test ParseOutFiles for tab files
		"""
		parse_out_files = ParseOutFiles()
		tab_file_to_process = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "run_freebayes_parallel.tab")
		expected_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "out_aggregate_tab_files_freebays_2.tab")

		temp_out_file = self.utils.get_temp_file('test_tab_out', '.tab')
		vect_type_out = []
		vect_type_remove = ['snp']
		b_add_header = True
		n_out_lines = 0
		with open(temp_out_file, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=Constants.SEPARATOR_TAB, quotechar='"', quoting=csv.QUOTE_ALL)
			n_out_lines += parse_out_files.parse_tab_files('sample_name', tab_file_to_process, csv_writer, vect_type_out, vect_type_remove, 101, b_add_header)
		
		self.assertEquals(4, n_out_lines)
		self.assertTrue(os.path.exists(temp_out_file))
		self.assertTrue(filecmp.cmp(expected_file, temp_out_file))
		if (os.path.exists(temp_out_file)): os.unlink(temp_out_file)

	def test_add_variants_in_incomplete_locus(self):
		"""
		test ParseOutFiles for tab files
		"""
		parse_out_files = ParseOutFiles()
		tab_file_to_process = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "have_multiple_var_snippy.tab")
		expected_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "out_snippy_translated.tsv")

		temp_file = self.utils.get_temp_file('test_add_variants_in_incomplete_locus', FileExtensions.FILE_TSV)
		self.utils.copy_file(tab_file_to_process, temp_file)
		
		coverage = Coverage()
		coverage.add_coverage('2', Coverage.COVERAGE_ALL, '1400.0')
		coverage.add_coverage('2', Coverage.COVERAGE_MORE_0, '100.0')
		coverage.add_coverage('2', Coverage.COVERAGE_MORE_9, '98.0')
		coverage.add_coverage('PB2', Coverage.COVERAGE_ALL, '1400.0')
		coverage.add_coverage('PB2', Coverage.COVERAGE_MORE_0, '100.0')
		coverage.add_coverage('PB2', Coverage.COVERAGE_MORE_9, '100.0')
		
		###
		parse_out_files.add_variants_in_incomplete_locus(temp_file, coverage)
			
		self.assertTrue(os.path.exists(temp_file))
		self.assertTrue(filecmp.cmp(expected_file, temp_file))
		
		self.utils.copy_file(tab_file_to_process, temp_file)
		coverage = Coverage()
		coverage.add_coverage('2', Coverage.COVERAGE_ALL, '1400.0')
		coverage.add_coverage('2', Coverage.COVERAGE_MORE_0, '100.0')
		coverage.add_coverage('2', Coverage.COVERAGE_MORE_9, '98.0')
		coverage.add_coverage('PB2', Coverage.COVERAGE_ALL, '1400.0')
		coverage.add_coverage('PB2', Coverage.COVERAGE_MORE_0, '98.0')
		coverage.add_coverage('PB2', Coverage.COVERAGE_MORE_9, '100.0')
		
		###
		parse_out_files.add_variants_in_incomplete_locus(temp_file, coverage)
		
		expected_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "out_snippy_translated_2.tsv")
		self.assertTrue(os.path.exists(temp_file))
		self.assertTrue(filecmp.cmp(expected_file, temp_file))
		if (os.path.exists(temp_file)): os.unlink(temp_file)


