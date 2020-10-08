'''
Created on Nov 17, 2017

@author: mmp
'''
import unittest
import os
from utils.parse_coverage_file import GetCoverage, ParseFile
from django.conf import settings 
from utils.result import Coverage
from constants.constantsTestsCase import ConstantsTestsCase

class Test(unittest.TestCase):

	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_COVERAGE)

	def testReference(self):
		reference_file = os.path.join(self.baseDirectory, "files_2/reference_zeros.fasta")
		get_coverage = GetCoverage()
		get_coverage.read_reference_fasta(reference_file)
		
		self.assertTrue("1" in get_coverage.get_dict_reference())
		self.assertTrue("2" in get_coverage.get_dict_reference())
		self.assertFalse("3" in get_coverage.get_dict_reference())
		self.assertEquals(20, get_coverage.get_dict_reference()["1"])
		self.assertEquals(20, get_coverage.get_dict_reference()["2"])
		
	def testReferenceGz(self):
		reference_file = os.path.join(self.baseDirectory, "files_2/reference_zeros.fasta.gz")
		get_coverage = GetCoverage()
		get_coverage.read_reference_fasta(reference_file)
		
		self.assertTrue("1" in get_coverage.get_dict_reference())
		self.assertTrue("2" in get_coverage.get_dict_reference())
		self.assertFalse("3" in get_coverage.get_dict_reference())
		self.assertEquals(10, get_coverage.get_dict_reference()["1"])
		self.assertEquals(10, get_coverage.get_dict_reference()["2"])
		
	def testFile(self):
		
		reference_file = os.path.join(self.baseDirectory, "ref/ref_H3.fasta")
		input_file = os.path.join(self.baseDirectory, "EVA001_S66.depth")
		
		get_coverage = GetCoverage()
		parse_file = ParseFile()
		data_file = parse_file.parse_file(input_file)
		get_coverage.read_reference_fasta(reference_file)

		self.assertEqual(len(data_file.get_vect_chromosomes()), 8)
		self.assertEqual(data_file.get_vect_chromosomes()[0], "1")
		self.assertEqual(data_file.get_vect_chromosomes()[-1], "8")
		self.assertEqual(get_coverage.get_dict_reference()["1"], 2280)
		self.assertEqual(get_coverage.get_dict_reference()["8"], 838)
		self.assertEqual("%.2f" % data_file.get_coverage(data_file.get_vect_chromosomes()[0], get_coverage.get_dict_reference()["1"]), "7.78")
		self.assertEqual("%.2f" % data_file.get_coverage(data_file.get_vect_chromosomes()[-1], get_coverage.get_dict_reference()["8"]), "33.57")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than(data_file.get_vect_chromosomes()[0], get_coverage.get_dict_reference()["1"], 9), "0.32")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than(data_file.get_vect_chromosomes()[-1], get_coverage.get_dict_reference()["8"], 9), "0.96")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than(data_file.get_vect_chromosomes()[0], get_coverage.get_dict_reference()["1"], 1), "0.98")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than(data_file.get_vect_chromosomes()[-1], get_coverage.get_dict_reference()["8"], 1), "1.00")
		self.assertEqual(data_file.get_file_name(), "EVA001_S66")

	def testFileRef(self):
		reference_file = os.path.join(self.baseDirectory, "ref/ompA_amplicon_14_sero.fasta")
		
		get_coverage = GetCoverage()
		get_coverage.read_reference_fasta(reference_file)
		self.assertEqual(len(get_coverage.get_dict_reference()), 14)
		self.assertEqual(len(get_coverage.get_vect_reference()), 14)
		self.assertTrue('A_Har13_NC007429' in get_coverage.get_dict_reference())
		self.assertTrue('F_IC_Cal3_X52080' in get_coverage.get_dict_reference())
		self.assertTrue('L2b_UCH_1_NC_010280' in get_coverage.get_dict_reference())
		self.assertTrue(get_coverage.get_vect_reference()[0] in get_coverage.get_dict_reference())
		self.assertTrue(get_coverage.get_vect_reference()[-1] in get_coverage.get_dict_reference())
		self.assertTrue(get_coverage.get_vect_reference()[3] in get_coverage.get_dict_reference())


	def testFile_1(self):
		reference_file = os.path.join(self.baseDirectory, "ref/ref_H3.fasta")
		input_file = os.path.join(self.baseDirectory, "EVA003_S91.depth.gz")
		
		get_coverage = GetCoverage()
		parse_file = ParseFile()
		data_file = parse_file.parse_file(os.path.join(os.getcwd(), input_file))
		get_coverage.read_reference_fasta(reference_file)
				
		self.assertEqual(len(data_file.get_vect_chromosomes()), 8)
		self.assertEqual(data_file.get_vect_chromosomes()[0], "1")
		self.assertEqual(data_file.get_vect_chromosomes()[-1], "8")
		self.assertEqual(get_coverage.get_dict_reference()["1"], 2280)
		self.assertEqual(get_coverage.get_dict_reference()["8"], 838)
		self.assertEqual("%.2f" % data_file.get_coverage(data_file.get_vect_chromosomes()[0], get_coverage.get_dict_reference()["1"]), "905.41")
		self.assertEqual("%.2f" % data_file.get_coverage(data_file.get_vect_chromosomes()[-1], get_coverage.get_dict_reference()["8"]), "1752.53")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than(data_file.get_vect_chromosomes()[0], get_coverage.get_dict_reference()["1"], 9), "1.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than(data_file.get_vect_chromosomes()[-1], get_coverage.get_dict_reference()["8"], 9), "1.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than(data_file.get_vect_chromosomes()[0], get_coverage.get_dict_reference()["1"], 1), "1.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than(data_file.get_vect_chromosomes()[-1], get_coverage.get_dict_reference()["8"], 1), "1.00")
		self.assertEqual(data_file.get_file_name(), "EVA003_S91")

	def testFile_2(self):
		reference_file = os.path.join(self.baseDirectory, "ref/ref_H3.fasta")
		input_file = os.path.join(self.baseDirectory, "EVA003_S91.depth")

		get_coverage = GetCoverage()
		parse_file = ParseFile()
		data_file = parse_file.parse_file(os.path.join(os.getcwd(), input_file))
		get_coverage.read_reference_fasta(reference_file)
		
		self.assertEqual(len(data_file.get_vect_chromosomes()), 8)
		self.assertEqual(data_file.get_vect_chromosomes()[0], "1")
		self.assertEqual(data_file.get_vect_chromosomes()[-1], "8")
		self.assertEqual(get_coverage.get_dict_reference()["1"], 2280)
		self.assertEqual(get_coverage.get_dict_reference()["8"], 838)
		self.assertEqual("%.2f" % data_file.get_coverage(data_file.get_vect_chromosomes()[0], get_coverage.get_dict_reference()["1"]), "905.41")
		self.assertEqual("%.2f" % data_file.get_coverage(data_file.get_vect_chromosomes()[-1], get_coverage.get_dict_reference()["8"]), "1752.53")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than(data_file.get_vect_chromosomes()[0], get_coverage.get_dict_reference()["1"], 9), "1.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than(data_file.get_vect_chromosomes()[-1], get_coverage.get_dict_reference()["8"], 9), "1.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than(data_file.get_vect_chromosomes()[0], get_coverage.get_dict_reference()["1"], 1), "1.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than(data_file.get_vect_chromosomes()[-1], get_coverage.get_dict_reference()["8"], 1), "1.00")
		self.assertEqual(data_file.get_file_name(), "EVA003_S91")

	def testFile_with_zeros(self):
		reference_file = os.path.join(self.baseDirectory, "files_3/reference_zeros.fasta")
		input_file = os.path.join(self.baseDirectory, "files_3/EVA001_S67_zeros.depth")
		
		get_coverage = GetCoverage()
		parse_file = ParseFile()
		data_file = parse_file.parse_file(os.path.join(os.getcwd(), input_file))
		get_coverage.read_reference_fasta(reference_file)

		self.assertEqual(len(data_file.get_vect_chromosomes()), 2)
		self.assertEqual(data_file.get_vect_chromosomes()[0], "1")
		self.assertEqual(data_file.get_vect_chromosomes()[-1], "2")
		self.assertEqual(get_coverage.get_dict_reference()["1"], 20)
		self.assertEqual(get_coverage.get_dict_reference()["2"], 20)
		self.assertEqual("%.2f" % data_file.get_coverage("1", get_coverage.get_dict_reference()["1"]), "5.00")
		self.assertEqual("%.2f" % data_file.get_coverage("2", get_coverage.get_dict_reference()["2"]), "5.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than("1", get_coverage.get_dict_reference()["1"], 9), "0.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than("2", get_coverage.get_dict_reference()["2"], 9), "0.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than("1", get_coverage.get_dict_reference()["1"], 1), "1.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than("2", get_coverage.get_dict_reference()["2"], 1), "1.00")

	def testFile_with_zeros_1(self):
		reference_file = os.path.join(self.baseDirectory, "files_2/reference_zeros.fasta")
		input_file = os.path.join(self.baseDirectory, "files_2/EVA001_S67_zeros_2.depth")
		
		get_coverage = GetCoverage()
		parse_file = ParseFile()
		data_file = parse_file.parse_file(os.path.join(os.getcwd(), input_file))
		get_coverage.read_reference_fasta(reference_file)

		self.assertEqual(len(data_file.get_vect_chromosomes()), 2)
		self.assertEqual(data_file.get_vect_chromosomes()[0], "1")
		self.assertEqual(data_file.get_vect_chromosomes()[-1], "2")
		self.assertEqual(len(get_coverage.get_vect_reference()), 2)
		self.assertEqual(get_coverage.get_dict_reference()["1"], 20)
		self.assertEqual(get_coverage.get_dict_reference()["2"], 20)
		self.assertEqual("%.2f" % data_file.get_coverage("1", get_coverage.get_dict_reference()["1"]), "4.50")
		self.assertEqual("%.2f" % data_file.get_coverage("2", get_coverage.get_dict_reference()["2"]), "5.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than("1", get_coverage.get_dict_reference()["1"], 9), "0.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than("2", get_coverage.get_dict_reference()["2"], 9), "0.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than("1", get_coverage.get_dict_reference()["1"], 1), "0.90")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than("2", get_coverage.get_dict_reference()["2"], 1), "1.00")
		
	def testFile_with_zeros_fault(self):
		reference_file = os.path.join(self.baseDirectory, "files_2/reference_zeros_fault.fasta")
		input_file = os.path.join(self.baseDirectory, "files_2/EVA001_S67_zeros.depth")
		
		get_coverage = GetCoverage()
		parse_file = ParseFile()
		data_file = parse_file.parse_file(os.path.join(os.getcwd(), input_file))
		get_coverage.read_reference_fasta(reference_file)

		self.assertEqual(len(data_file.get_vect_chromosomes()), 2)
		self.assertEqual(data_file.get_vect_chromosomes()[0], "1")
		self.assertEqual(data_file.get_vect_chromosomes()[-1], "2")
		self.assertEqual(get_coverage.get_dict_reference()["1"], 10)
		self.assertEqual(get_coverage.get_dict_reference()["2"], 10)
		try:
			self.assertEqual("%.2f" % data_file.get_coverage("1", get_coverage.get_dict_reference()["1"]), "10.00")
			self.fail("Must raise exception")
		except Exception as e:
			self.assertEqual("Chromosome '1' has different sizes. Coverage: 20; Reference: 10", e.args[0])
			
		try:
			self.assertEqual("%.2f" % data_file.get_ratio_more_than("1", get_coverage.get_dict_reference()["1"], 4), "2.00")
			self.fail("Must raise exception")
		except Exception as e:
			self.assertEqual("Chromosome '1' has different sizes. Coverage: 20; Reference: 10", e.args[0])
		
	def testFile_with_zeros_2(self):
		reference_file = os.path.join(self.baseDirectory, "files_3/reference_zeros.fasta")
		input_file = os.path.join(self.baseDirectory, "files_3/EVA001_S67_zeros.depth")
		
		get_coverage = GetCoverage()
		parse_file = ParseFile()
		data_file = parse_file.parse_file(os.path.join(os.getcwd(), input_file))
		get_coverage.read_reference_fasta(reference_file)

		self.assertEqual(len(data_file.get_vect_chromosomes()), 2)
		self.assertEqual(data_file.get_vect_chromosomes()[0], "1")
		self.assertEqual(data_file.get_vect_chromosomes()[-1], "2")
		self.assertEqual(get_coverage.get_dict_reference()["1"], 20)
		self.assertEqual(get_coverage.get_dict_reference()["2"], 20)
		self.assertEqual(get_coverage.get_dict_reference()["3"], 20)
		self.assertEqual("%.2f" % data_file.get_coverage("1", get_coverage.get_dict_reference()["1"]), "5.00")
		self.assertEqual("%.2f" % data_file.get_coverage("2", get_coverage.get_dict_reference()["2"]), "5.00")
		self.assertEqual("%.2f" % data_file.get_coverage("3", get_coverage.get_dict_reference()["3"]), "0.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than("1", get_coverage.get_dict_reference()["1"], 9), "0.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than("2", get_coverage.get_dict_reference()["2"], 9), "0.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than("3", get_coverage.get_dict_reference()["3"], 9), "0.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than("1", get_coverage.get_dict_reference()["1"], 1), "1.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than("2", get_coverage.get_dict_reference()["2"], 1), "1.00")
		self.assertEqual("%.2f" % data_file.get_ratio_more_than("3", get_coverage.get_dict_reference()["3"], 1), "0.00")


	def test_coverage(self):
		reference_file = os.path.join(self.baseDirectory, "files_3/reference_zeros.fasta")
		input_file = os.path.join(self.baseDirectory, "files_3/EVA001_S67_zeros.depth")
		
		get_coverage = GetCoverage()
		coverage = get_coverage.get_coverage(input_file, reference_file)

		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_ALL), "5.0")
		self.assertEqual(coverage.get_coverage('2', Coverage.COVERAGE_ALL), "5.0")
		self.assertEqual(coverage.get_coverage('3', Coverage.COVERAGE_ALL), "0.0")
		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_MORE_9), "0.0")
		self.assertEqual(coverage.get_coverage('2', Coverage.COVERAGE_MORE_9), "0.0")
		self.assertEqual(coverage.get_coverage('3', Coverage.COVERAGE_MORE_9), "0.0")
		self.assertFalse(coverage.is_100_more_9('3'))
		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('2', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('3', Coverage.COVERAGE_MORE_0), "0.0")

	def test_coverage_2(self):
		reference_file = os.path.join(self.baseDirectory, "ref/ref_H3.fasta")
		input_file = os.path.join(self.baseDirectory, "EVA003_S91.depth")

		get_coverage = GetCoverage()
		coverage = get_coverage.get_coverage(input_file, reference_file)
		
		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_ALL), "905.4")
		self.assertEqual(coverage.get_coverage('8', Coverage.COVERAGE_ALL), "1752.5")
		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_MORE_9), "100.0")
		self.assertEqual(coverage.get_coverage('8', Coverage.COVERAGE_MORE_9), "100.0")
		self.assertTrue(coverage.is_100_more_9('8'))
		self.assertTrue(coverage.is_100_more_0('8'))
		self.assertTrue(coverage.is_100_more_9('8'))
		self.assertEqual("Fail, locus '8': the % of locus size covered by at least 10-fold is '100.0%' (below 100%)", coverage.get_fault_message_9('8'))
		self.assertEqual("Fail, locus '8': the % of locus size covered by at least 1-fold is '100.0%' (below 100%)", coverage.get_fault_message_0('8'))
		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('8', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertTrue(coverage.is_100_more_9('8'))
		self.assertTrue(coverage.is_100_more_9('8'))
		self.assertTrue(coverage.is_100_more_0('8'))
		

	def test_coverage_define_by_user_1(self):
		reference_file = os.path.join(self.baseDirectory, "files_3/reference_zeros.fasta")
		input_file = os.path.join(self.baseDirectory, "files_3/EVA001_S67_zeros.depth")
		
		get_coverage = GetCoverage()
		defined_by_user_limit = 3
		coverage = get_coverage.get_coverage(input_file, reference_file, defined_by_user_limit)

		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_ALL), "5.0")
		self.assertTrue(coverage.is_100_more_9('1'))
		self.assertEqual(coverage.get_coverage('2', Coverage.COVERAGE_ALL), "5.0")
		self.assertEqual(coverage.get_coverage('3', Coverage.COVERAGE_ALL), "0.0")
		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_MORE_9), "0.0")
		self.assertEqual(coverage.get_coverage('2', Coverage.COVERAGE_MORE_9), "0.0")
		self.assertEqual(coverage.get_coverage('3', Coverage.COVERAGE_MORE_9), "0.0")
		self.assertFalse(coverage.is_100_more_9('3'))
		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('2', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('3', Coverage.COVERAGE_MORE_0), "0.0")
		
	def test_coverage_define_by_user_2(self):
		reference_file = os.path.join(self.baseDirectory, "files_3/reference_zeros.fasta")
		input_file = os.path.join(self.baseDirectory, "files_3/EVA001_S67_zeros.depth")
		
		get_coverage = GetCoverage()
		defined_by_user_limit = 4
		coverage = get_coverage.get_coverage(input_file, reference_file, defined_by_user_limit)

		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_ALL), "5.0")
		self.assertTrue(coverage.is_100_more_9('1'))
		self.assertEqual(coverage.get_coverage('2', Coverage.COVERAGE_ALL), "5.0")
		self.assertEqual(coverage.get_coverage('3', Coverage.COVERAGE_ALL), "0.0")
		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_MORE_9), "0.0")
		self.assertEqual(coverage.get_coverage('2', Coverage.COVERAGE_MORE_9), "0.0")
		self.assertEqual(coverage.get_coverage('3', Coverage.COVERAGE_MORE_9), "0.0")
		self.assertFalse(coverage.is_100_more_9('3'))
		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('2', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('3', Coverage.COVERAGE_MORE_0), "0.0")
	
	def test_coverage_define_by_user_3(self):
		reference_file = os.path.join(self.baseDirectory, "files_3/reference_zeros.fasta")
		input_file = os.path.join(self.baseDirectory, "files_3/EVA001_S67_zeros.depth")
		
		get_coverage = GetCoverage()
		defined_by_user_limit = 5
		coverage = get_coverage.get_coverage(input_file, reference_file, defined_by_user_limit)

		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_ALL), "5.0")
		self.assertTrue(coverage.is_100_more_9('1'))
		self.assertEqual(coverage.get_coverage('2', Coverage.COVERAGE_ALL), "5.0")
		self.assertEqual(coverage.get_coverage('3', Coverage.COVERAGE_ALL), "0.0")
		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_MORE_9), "0.0")
		self.assertEqual(coverage.get_coverage('2', Coverage.COVERAGE_MORE_9), "0.0")
		self.assertEqual(coverage.get_coverage('3', Coverage.COVERAGE_MORE_9), "0.0")
		self.assertFalse(coverage.is_100_more_9('3'))
		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('2', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('3', Coverage.COVERAGE_MORE_0), "0.0")
	
	def test_coverage_define_by_user_4(self):
		reference_file = os.path.join(self.baseDirectory, "files_3/reference_zeros.fasta")
		input_file = os.path.join(self.baseDirectory, "files_3/EVA001_S67_zeros.depth")
		
		get_coverage = GetCoverage()
		defined_by_user_limit = 6
		coverage = get_coverage.get_coverage(input_file, reference_file, defined_by_user_limit)

		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_ALL), "5.0")
		self.assertFalse(coverage.is_100_more_9('1'))
		self.assertEqual(coverage.get_coverage('2', Coverage.COVERAGE_ALL), "5.0")
		self.assertEqual(coverage.get_coverage('3', Coverage.COVERAGE_ALL), "0.0")
		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_MORE_9), "0.0")
		self.assertEqual(coverage.get_coverage('2', Coverage.COVERAGE_MORE_9), "0.0")
		self.assertEqual(coverage.get_coverage('3', Coverage.COVERAGE_MORE_9), "0.0")
		self.assertFalse(coverage.is_100_more_9('3'))
		self.assertEqual(coverage.get_coverage('1', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('2', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('3', Coverage.COVERAGE_MORE_0), "0.0")

