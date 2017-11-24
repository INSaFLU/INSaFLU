'''
Created on Oct 28, 2017

@author: mmp
'''
from django.test import TestCase
from utils.result import Output, SoftwareDesc, DecodeResult, Result, ResultAverageAndNumberReads, CountHits
from utils.result import DecodeResultAverageAndNumberReads, Coverage, DecodeCoverage
from utils.software import Software

class Test(TestCase):

	def test_result(self):
		result = Result()
		result.set_error("xpto")
		result.add_output(Output("name", 'path'))
		result.add_software(SoftwareDesc("name_s", '11212', 'ads'))
		
		sz_return = result.to_json()
		self.assertTrue(sz_return.find('"version": "11212"') != 0)
		self.assertTrue(sz_return.find('"name": "name_s"') != 0)
		
		decodeResult = DecodeResult()
		result_2 = decodeResult.decode_result(sz_return)
		self.assertEqual(result_2.result, result.result)
		self.assertEqual(len(result_2.outputs.list_output), len(result.outputs.list_output))
		self.assertEqual(result_2.outputs.list_output[0].file_name, result.outputs.list_output[0].file_name)
		self.assertEqual(len(result_2.softwares.list_software), len(result.softwares.list_software))
		self.assertEqual(result_2.softwares.list_software[0].name, result.softwares.list_software[0].name)
		
	def test_result_software(self):
		result = Result()
		result.set_error("xpto")
		result.add_output(Output("name", 'path'))
		result.outputs.list_output[0].set_software(SoftwareDesc("name44_s", '4444', 'sdfd'))
		result.add_software(SoftwareDesc("name_s", '11212', "sdd"))
		
		sz_return = result.to_json()
		self.assertTrue(sz_return.find('"version": "11212"') != 0)
		self.assertTrue(sz_return.find('"name": "name_s"') != 0)
		
		decodeResult = DecodeResult()
		result_2 = decodeResult.decode_result(sz_return)
		self.assertEqual(result_2.result, result.result)
		self.assertEqual(len(result_2.outputs.list_output), len(result.outputs.list_output))
		self.assertEqual(result_2.outputs.list_output[0].file_name, result.outputs.list_output[0].file_name)
		self.assertEqual(result_2.outputs.list_output[0].software.name, result.outputs.list_output[0].software.name)
		self.assertEqual(result_2.outputs.list_output[0].software.parameters, result.outputs.list_output[0].software.parameters)
		self.assertEqual(len(result_2.softwares.list_software), len(result.softwares.list_software))
		self.assertEqual(result_2.softwares.list_software[0].name, result.softwares.list_software[0].name)

	def test_ResultAverageAndNumberReads(self):
		resultAverageAndNumberReads = ResultAverageAndNumberReads(21, 43, 53, 12)
		self.assertEqual(21, resultAverageAndNumberReads.number_file_1)
		self.assertEqual(43, resultAverageAndNumberReads.average_file_1)
		self.assertEqual(53, resultAverageAndNumberReads.number_file_2)
		self.assertEqual(12, resultAverageAndNumberReads.average_file_2)
		
		sz_return = resultAverageAndNumberReads.to_json()
		self.assertTrue(sz_return.find('"number_file_1": "21"') != 0)
		self.assertTrue(sz_return.find('"average_file_1": "43"') != 0)
		
		decodeResultAverageAndNumberReads = DecodeResultAverageAndNumberReads()
		result_2 = decodeResultAverageAndNumberReads.decode_result(sz_return)
		self.assertEqual(result_2, resultAverageAndNumberReads)
		
	def test_ResultSoftware(self):
		result = Result()
		result.set_error("xpto")
		result.add_software(SoftwareDesc(Software.SOFTWARE_SPAdes_name, Software.SOFTWARE_SPAdes_VERSION, Software.SOFTWARE_SPAdes_PARAMETERS))
		result.add_software(SoftwareDesc(Software.SOFTWARE_TRIMMOMATIC_name, Software.SOFTWARE_TRIMMOMATIC_VERSION, Software.SOFTWARE_TRIMMOMATIC_PARAMETERS))
		
		self.assertEqual("Trimmomatic-0.27; (SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:55 TOPHRED33)", result.get_software(Software.SOFTWARE_TRIMMOMATIC_name))
		self.assertEqual("SPAdes-3.11.1", result.get_software(Software.SOFTWARE_SPAdes_name))
		
		
	def test_Coverage(self):
		
		coverage = Coverage()	
		coverage.add_coverage('3', Coverage.COVERAGE_ALL, '10')
		coverage.add_coverage('3', Coverage.COVERAGE_MORE_0, '0')
		coverage.add_coverage('3', Coverage.COVERAGE_MORE_9, '9')
		coverage.add_coverage('2', Coverage.COVERAGE_ALL, '101')
		coverage.add_coverage('2', Coverage.COVERAGE_MORE_0, '10')
		coverage.add_coverage('2', Coverage.COVERAGE_MORE_9, '19')
		
		json = coverage.to_json()
		decodeCoverage = DecodeCoverage()
		coverage_2 = decodeCoverage.decode_result(json)
		self.assertTrue("10", coverage_2.get_coverage('3', Coverage.COVERAGE_ALL))
		self.assertTrue("0", coverage_2.get_coverage('3', Coverage.COVERAGE_MORE_0))
		self.assertTrue("9", coverage_2.get_coverage('3', Coverage.COVERAGE_MORE_9))
		self.assertTrue("101", coverage_2.get_coverage('2', Coverage.COVERAGE_ALL))
		self.assertTrue("10", coverage_2.get_coverage('2', Coverage.COVERAGE_MORE_0))
		self.assertTrue("19", coverage_2.get_coverage('2', Coverage.COVERAGE_MORE_9))
		try:
			self.assertTrue("19", coverage_2.get_coverage('4', Coverage.COVERAGE_MORE_9))
			self.fail("must raise exception")
		except Exception as e:
			self.assertEquals("Error: there's no key like this: 4", e.args[0])
		
		try:
			self.assertTrue("19", coverage_2.get_coverage('2', "xpto"))
			self.fail("must raise exception")
		except Exception as e:
			self.assertEquals("Error: there's no key like this: xpto", e.args[0])


	def test_count_hits(self):
		count_hits = CountHits()	
		count_hits.set_hits_50_90(40)
		count_hits.set_hits_less_50(140)
		
		json = count_hits.to_json()
		decodeCoverage = DecodeCoverage()
		count_hits_2 = decodeCoverage.decode_result(json)
		self.assertEquals(count_hits_2, count_hits)

		count_hits_2.add_one_hits_less_50()
		count_hits.add_one_hits_50_90()
		self.assertNotEqual(count_hits_2, count_hits)
		self.assertEquals(count_hits_2.get_hits_less_50(), 141)
		self.assertEquals(count_hits.get_hits_50_90(), 41)

		count_hits.add_one_hits_less_50()
		count_hits_2.add_one_hits_50_90()
		self.assertEquals(count_hits_2, count_hits)



