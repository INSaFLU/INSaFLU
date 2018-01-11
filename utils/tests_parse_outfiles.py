'''
Created on Jan 7, 2018

@author: mmp
'''
import unittest
from constants.constantsTestsCase import ConstantsTestsCase
from utils.utils import Utils
from utils.software import Software
from django.conf import settings
from utils.parse_out_files import ParseOutFiles
from utils.collect_extra_data import CollectExtraData
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
		b_add_header = True
		with open(temp_out_file, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=CollectExtraData.SEPARATOR_TAB, quotechar='"', quoting=csv.QUOTE_ALL)
			parse_out_files.parse_tab_files('sample_name', tab_file_to_process, csv_writer, vect_type_out, 101, b_add_header)
			
			vect_type_out = ['del']
			b_add_header = False
			parse_out_files.parse_tab_files('sample_name', tab_file_to_process, csv_writer, vect_type_out, 101, b_add_header)

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
		b_add_header = True
		with open(temp_out_file, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=CollectExtraData.SEPARATOR_TAB, quotechar='"', quoting=csv.QUOTE_ALL)
			parse_out_files.parse_tab_files('sample_name', tab_file_to_process, csv_writer, vect_type_out, 101, b_add_header)
			
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
		b_add_header = True
		with open(temp_out_file, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=CollectExtraData.SEPARATOR_TAB, quotechar='"', quoting=csv.QUOTE_ALL)
			parse_out_files.parse_tab_files('sample_name', tab_file_to_process, csv_writer, vect_type_out, 50, b_add_header)
			
		self.assertTrue(os.path.exists(temp_out_file))
		self.assertTrue(filecmp.cmp(expected_file, temp_out_file))
		if (os.path.exists(temp_out_file)): os.unlink(temp_out_file)


