'''
Created on Oct 28, 2017

@author: mmp
'''
from django.test import TestCase
from utils.result import Output, SoftwareDesc, Result, ResultAverageAndNumberReads, CountHits, MixedInfectionMainVector
from utils.result import Coverage, DecodeObjects, TasksToProcess, GeneticElement, Gene
from utils.result import ProcessResults, SingleResult, Coverage, DecodeObjects, TasksToProcess, GeneticElement, Gene
from constants.software_names import SoftwareNames

class Test(TestCase):

	def test_result(self):
		result = Result()
		result.set_error("xpto")
		result.add_output(Output("name", 'path'))
		result.add_software(SoftwareDesc("name_s", '11212', 'ads'))
		
		sz_return = result.to_json()
		self.assertTrue(sz_return.find('"version": "11212"') != 0)
		self.assertTrue(sz_return.find('"name": "name_s"') != 0)
		
		decodeResult = DecodeObjects()
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
		
		self.assertFalse(result.is_success())
		sz_return = result.to_json()
		self.assertTrue(sz_return.find('"version": "11212"') != 0)
		self.assertTrue(sz_return.find('"name": "name_s"') != 0)
		
		decodeResult = DecodeObjects()
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
		
		decodeResultAverageAndNumberReads = DecodeObjects()
		result_2 = decodeResultAverageAndNumberReads.decode_result(sz_return)
		self.assertEqual(result_2, resultAverageAndNumberReads)
		
	def test_ResultSoftware(self):
		result = Result()
		result.set_error("xpto")
		result.add_software(SoftwareDesc(SoftwareNames.SOFTWARE_SPAdes_name, SoftwareNames.SOFTWARE_SPAdes_VERSION, SoftwareNames.SOFTWARE_SPAdes_PARAMETERS))
		result.add_software(SoftwareDesc(SoftwareNames.SOFTWARE_TRIMMOMATIC_name, SoftwareNames.SOFTWARE_TRIMMOMATIC_VERSION, SoftwareNames.SOFTWARE_TRIMMOMATIC_PARAMETERS))
		
		self.assertEqual("Trimmomatic-0.27; (SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:55 TOPHRED33)", result.get_software(SoftwareNames.SOFTWARE_TRIMMOMATIC_name))
		self.assertEqual("SPAdes-3.11.1", result.get_software(SoftwareNames.SOFTWARE_SPAdes_name))
		
		
	def test_Coverage(self):
		
		coverage = Coverage()	
		coverage.add_coverage('3', Coverage.COVERAGE_ALL, '10')
		coverage.add_coverage('3', Coverage.COVERAGE_MORE_0, '0')
		coverage.add_coverage('3', Coverage.COVERAGE_MORE_9, '9')
		coverage.add_coverage('2', Coverage.COVERAGE_ALL, '101')
		coverage.add_coverage('2', Coverage.COVERAGE_MORE_0, '10')
		coverage.add_coverage('2', Coverage.COVERAGE_MORE_9, '19')
		
		json = coverage.to_json()
		decodeCoverage = DecodeObjects()
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
		count_hits.set_hits_more_90(10)
		
		json = count_hits.to_json()
		decodeCoverage = DecodeObjects()
		count_hits_2 = decodeCoverage.decode_result(json)
		self.assertEquals(count_hits_2, count_hits)

		count_hits_2.add_one_hits_less_50()
		count_hits.add_one_hits_50_90()
		self.assertNotEqual(count_hits_2, count_hits)
		self.assertEquals(count_hits_2.get_hits_less_50(), 141)
		self.assertEquals(count_hits.get_hits_50_90(), 41)
		self.assertEquals(count_hits.get_hits_more_90(), 10)

		count_hits.add_one_hits_less_50()
		count_hits_2.add_one_hits_50_90()
		count_hits_2.add_one_hits_more_90()
		self.assertNotEqual(count_hits_2, count_hits)
		
		count_hits.add_one_hits_more_90()
		self.assertEquals(count_hits_2, count_hits)

	def test_tasks_to_process(self):
		tasksToProcess = TasksToProcess()	
		tasksToProcess.add_taskd_id('sdsdds')
		tasksToProcess.add_taskd_id('sdsdds1')
		tasksToProcess.add_taskd_id('sdsdds3')
		
		json = tasksToProcess.to_json()
		decodeCoverage = DecodeObjects()
		tasksToProcess_2 = decodeCoverage.decode_result(json)
		self.assertEquals(tasksToProcess, tasksToProcess_2)


	def test_elements_genes(self):
		
		geneticElement = GeneticElement()
		self.assertTrue(geneticElement.add_gene('element_name', Gene('name', 12, 45, 1)))
		self.assertFalse(geneticElement.add_gene('element_name', Gene('name', 12, 45, 1)))
		self.assertTrue(geneticElement.add_gene('element_name', Gene('name2', 35, 55, -1)))
		self.assertTrue(geneticElement.add_gene('element_name2', Gene('name', 12, 45, 1)))
		geneticElement.add_gene('element_name2', Gene('name2', 212, 245, 1))
		geneticElement.add_gene('element_name2', Gene('name3', 412, 445, -1))
		geneticElement.add_gene('element', Gene('name', 412, 445, -1))

		json = geneticElement.to_json()
		decodeCoverage = DecodeObjects()
		geneticElement_2 = decodeCoverage.decode_result(json)
		self.assertEquals(geneticElement, geneticElement_2)
		self.assertEquals("element,element_name,element_name2", ",".join(geneticElement_2.get_sorted_elements()))
		self.assertEquals(2, len(geneticElement_2.get_genes('element_name')))


	def test_mixed_infection_main_vector(self):
		
		mixedInfectionMainVector = MixedInfectionMainVector()

		self.assertEquals(8, len(mixedInfectionMainVector.get_vector()))
		mixedInfectionMainVector.add_vector([10, 34])
		mixedInfectionMainVector.add_vector([78, 56])
		self.assertEquals(9, len(mixedInfectionMainVector.get_vector()))
		json = mixedInfectionMainVector.to_json()
		
		decodeCoverage = DecodeObjects()
		mixedInfectionMainVector_2 = decodeCoverage.decode_result(json)
		self.assertEquals(9, len(mixedInfectionMainVector_2.get_vector()))


		mixedInfectionMainVector_3 = MixedInfectionMainVector()
		self.assertEquals(8, len(mixedInfectionMainVector_3.get_vector()))

	def test_process_results(self):
		
		process_results = ProcessResults()
		process_results.add_single_result(SingleResult(SingleResult.ERROR, 'xprto errroer'))
		process_results.add_single_result(SingleResult(SingleResult.SUCCESS, 'xprto'))
		self.assertEquals(2, len(process_results.get_vect_results()))
		
		json = process_results.to_json()
		
		decodeCoverage = DecodeObjects()
		process_results_2 = decodeCoverage.decode_result(json)
		self.assertEquals(process_results, process_results_2)
		self.assertEquals("Error - xprto errroer\nSuccess - xprto", str(process_results_2))

