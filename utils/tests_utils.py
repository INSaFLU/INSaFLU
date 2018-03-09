'''
Created on Oct 28, 2017

@author: mmp
'''
import unittest
from constants.constantsTestsCase import ConstantsTestsCase
from django.test.utils import override_settings
from utils.utils import Utils
from utils.software import Software
from django.conf import settings 
import os, filecmp
from constants.software_names import SoftwareNames
from utils.result import Coverage
from constants.constants import FileExtensions, Constants
from utils.result import Gene

class Test(unittest.TestCase):


	### static
	constantsTestsCase = ConstantsTestsCase()
	utils = Utils()
	software = Software()

	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), self.constantsTestsCase.MANAGING_TESTS)
		
	def tearDown(self):
		pass

	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_media_root(self):
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		
	def testgetFileNameWithoutExtension(self):
		
		self.assertEqual("fiels", self.utils.get_file_name_without_extension("/root/fiels.fasta"))
		self.assertEqual("fiels", self.utils.get_file_name_without_extension("/root/fiels"))
		self.assertEqual("fiels.fasta", self.utils.get_file_name_without_extension("/root/fiels.fasta.txt"))

	def test_is_fastq_gz(self):
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		self.assertTrue(self.utils.is_fastq_gz(path_file))
		
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_ABRICATE, ConstantsTestsCase.MANAGING_TEST_ABRICATE)
		try:
			self.assertTrue(self.utils.is_fastq_gz(path_file))
		except Exception as e:
			self.assertEqual("File need to have suffix '.fastq.gz'", e.args[0])
			
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		try:
			self.assertTrue(self.utils.is_fastq_gz(path_file))
		except Exception as e:
			self.assertEqual("File need to have suffix '.fastq.gz'", e.args[0])
			
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA_FAKE_GZ)
		try:
			self.assertTrue(self.utils.is_fastq_gz(path_file))
		except Exception as e:
			self.assertEqual("File is not in fastq.gz format.", e.args[0])
	
	def test_read_text_file(self):
		"""
		"""
		txt_file = os.path.join(self.baseDirectory, 'abricate', ConstantsTestsCase.MANAGING_TEST_ABRICATE + ".xpto")
		self.assertFalse(os.path.exists(txt_file))
		try:
			self.utils.read_text_file(txt_file)
			self.fail("must raise error")
		except IOError as e:
			self.assertTrue(e.args[0].endswith("static_all/tests/abricate/abricate_out.txt.xpto' doens't exist."))
		
		txt_file = os.path.join(self.baseDirectory, 'abricate', ConstantsTestsCase.MANAGING_TEST_ABRICATE)
		self.assertTrue(os.path.exists(txt_file))
		
		vect_out = self.utils.read_text_file(txt_file)
		self.assertEquals(6, len(vect_out))
		self.assertEquals(0, vect_out[0].find("#FILE"))
		self.assertEquals(-1, vect_out[0].find("Victoria"))
		self.assertEquals(87, vect_out[5].find("Victoria"))
		
	def test_gzip_and_tabix(self):
		vcf_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "temp.vcf")
		self.assertTrue(os.path.exists(vcf_file))
		
		temp_path = self.utils.get_temp_dir()
		test_file_to_remove = self.utils.get_temp_file("test_gzip_and_tabix", ".vcf")
		temp_file = os.path.join(temp_path, os.path.basename(test_file_to_remove))
		self.utils.copy_file(vcf_file, temp_file)
		self.assertTrue(os.path.exists(vcf_file))
		self.utils.compress_files(SoftwareNames.SOFTWARE_BGZIP, temp_file)
		self.utils.compress_files(SoftwareNames.SOFTWARE_BGZIP, temp_file)
		self.assertTrue(os.path.exists(temp_file + ".gz"))
		self.assertTrue(os.path.getsize(temp_file + ".gz") > 1000)
		self.software.create_index_files(temp_file + ".gz")
		self.assertTrue(os.path.exists(temp_file + ".gz.tbi"))
		self.utils.remove_temp_file(temp_file + ".gz.tbi")
		try:
			self.software.create_index_files(temp_file)
			self.fail("must raise exception")
		except Exception as e:
			self.assertEquals("File doesn't exist", e.args[0])
		os.unlink(test_file_to_remove)
		print("These two faults are expected: 'Not a BGZF file' and 'tbx_index_build failed'")
		self.utils.remove_dir(temp_path)
	

	def test_add_vcf_FREQ(self):
		vcf_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "temp.vcf")
		vcf_file_expected = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "temp_more_REF.vcf")
		temp_dir = self.utils.get_temp_dir()
		temp_vcf_file = os.path.join(temp_dir, os.path.basename(vcf_file))
		
		self.utils.copy_file(vcf_file, temp_vcf_file)
		
		self.software.test_bgzip_and_tbi_in_vcf(temp_vcf_file)
		vcf_file_out = self.utils.get_temp_file("temp_vcf", ".vcf")
		vcf_file_out = self.utils.add_freq_to_vcf(vcf_file, vcf_file_out)
		self.assertTrue(filecmp.cmp(vcf_file_expected, vcf_file_out))
	
		## remove file
		os.unlink(vcf_file_out)
		self.utils.remove_dir(temp_dir)

	def test_count_hits_from_tab(self):
		"""
		test cout hits
		"""
		tab_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "resutl_vcf_to_tab2.tab")
		vect_count_type = ['snp']
		count_hits = self.utils.count_hits_from_tab(tab_file, vect_count_type)
		self.assertEqual(4, count_hits.get_hits_less_50())
		self.assertEqual(2, count_hits.get_hits_50_90())
		self.assertEqual(118, count_hits.get_hits_more_90())
		self.assertEqual(6, count_hits.get_total_50_50_90())
		self.assertEqual(124, count_hits.get_total())

	def test_count_hits_from_tab_2(self):
		"""
		test cout hits
		"""
		tab_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "run_snippy2_other.tab")
		vect_count_type = ['snp']
		count_hits = self.utils.count_hits_from_tab(tab_file, vect_count_type)
		self.assertEqual(3, count_hits.get_hits_less_50())
		self.assertEqual(0, count_hits.get_hits_50_90())
		self.assertEqual(122, count_hits.get_hits_more_90())
		self.assertEqual(3, count_hits.get_total_50_50_90())
		self.assertEqual(125, count_hits.get_total())
		
	def test_count_hits_from_tab_3(self):
		"""
		test cout hits
		"""
		tab_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "run_snippy2_other.tab")
		vect_count_type = ['snp', 'del']
		count_hits = self.utils.count_hits_from_tab(tab_file, vect_count_type)
		self.assertEqual(4, count_hits.get_hits_less_50())
		self.assertEqual(0, count_hits.get_hits_50_90())
		self.assertEqual(122, count_hits.get_hits_more_90())
		self.assertEqual(4, count_hits.get_total_50_50_90())
		self.assertEqual(126, count_hits.get_total())
		
	def test_count_hits_from_tab_4(self):
		"""
		test cout hits
		"""
		tab_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "have_multiple_var.tab")
		vect_count_type = ['snp']
		count_hits = self.utils.count_hits_from_tab(tab_file, vect_count_type)
		self.assertEqual(1, count_hits.get_hits_less_50())
		self.assertEqual(1, count_hits.get_hits_50_90())
		self.assertEqual(2, count_hits.get_hits_more_90())
		self.assertEqual(2, count_hits.get_total_50_50_90())
		self.assertEqual(4, count_hits.get_total())
		
	def test_count_hits_from_tab_5(self):
		"""
		test cout hits
		"""
		tab_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "have_multiple_var.tab")
		vect_count_type = ['del']
		count_hits = self.utils.count_hits_from_tab(tab_file, vect_count_type)
		self.assertEqual(1, count_hits.get_hits_less_50())
		self.assertEqual(0, count_hits.get_hits_50_90())
		self.assertEqual(0, count_hits.get_hits_more_90())
		self.assertEqual(1, count_hits.get_total_50_50_90())
		self.assertEqual(1, count_hits.get_total())
		
	def test_get_variations_by_freq_from_tab(self):
		"""
		test get variations by freq
		"""
		tab_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "resutl_vcf_to_tab2.tab")
		vect_count_type = ['snp']
		(dict_less_50, dict_more_50, dict_more_90) = self.utils.get_variations_by_freq_from_tab(tab_file, vect_count_type)
		self.assertEqual(2, len(dict_less_50))
		self.assertEqual(817, dict_less_50['NS'][0])
		self.assertEqual(237, dict_less_50['PB2'][0])
		self.assertEqual(774, dict_less_50['PB2'][1])
		self.assertEqual(896, dict_less_50['PB2'][2])
		self.assertEqual(1, len(dict_more_50))
		self.assertEqual(1554, dict_more_50['PB2'][0])
		self.assertEqual(8, len(dict_more_90))
		self.assertEqual(324, dict_more_90['MP'][0])
		self.assertEqual(99, dict_more_90['NS'][0])
		self.assertEqual(117, dict_more_90['PA'][0])

	def test_get_elements_and_genes(self):
		"""
		test cout hits
		"""
		genbank_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS,\
					ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		geneticElement = self.utils.get_elements_and_genes(genbank_file)
		self.assertEquals(8, len(geneticElement.get_sorted_elements()))
		self.assertTrue('NS' in geneticElement.get_sorted_elements())
		self.assertFalse('xptoNS' in geneticElement.get_sorted_elements())
		self.assertTrue('PB2' in geneticElement.get_sorted_elements())
		self.assertEquals('HA,MP,NA,NP,NS,PA,PB1,PB2', ','.join(geneticElement.get_sorted_elements()))
		self.assertEquals(30, geneticElement.get_genes('PB2')[0].start)
		self.assertEquals(2280, geneticElement.get_size_element('PB2'))
		self.assertEquals(2280, geneticElement.get_genes('PB2')[0].end)
		self.assertEquals('PB2', geneticElement.get_genes('PB2')[0].name)
		self.assertEquals(1, geneticElement.get_genes('PB2')[0].strand)
		self.assertTrue(geneticElement.get_genes('PB2')[0].is_forward())
		self.assertEquals(1410, geneticElement.get_size_element('NA'))

		
	def test_filter_fasta_all_sequences(self):
		"""
		test filter fasta
		"""
		
		"""
		EVA011_S54
		MP - All 2198.8; 0 100.0; 9 100.0
		PA - All 527.4; 0 100.0; 9 100.0
		HA - All 1449.3; 0 100.0; 9 100.0
		PB1 - All 618.4; 0 100.0; 9 100.0
		NP - All 439.3; 0 100.0; 9 100.0
		NS - All 1214.7; 0 100.0; 9 100.0
		PB2 - All 690.6; 0 100.0; 9 100.0
		NA - All 1092.1; 0 100.0; 9 100.0
		
		EVA003_S91
		MP - All 2198.8; 0 100.0; 9 100.0
		PA - All 527.4; 0 100.0; 9 100.0
		HA - All 1449.3; 0 100.0; 9 100.0
		PB1 - All 618.4; 0 100.0; 9 100.0
		NP - All 439.3; 0 100.0; 9 100.0
		NS - All 1214.7; 0 100.0; 9 100.0
		PB2 - All 690.6; 0 100.0; 9 100.0
		NA - All 1092.1; 0 100.0; 9 100.0
		
		EVA002_S52
		MP - All 2198.8; 0 100.0; 9 100.0
		PA - All 527.4; 0 100.0; 9 100.0
		HA - All 1449.3; 0 100.0; 9 100.0
		PB1 - All 618.4; 0 100.0; 9 100.0
		NP - All 439.3; 0 100.0; 9 100.0
		NS - All 1214.7; 0 100.0; 9 100.0
		PB2 - All 690.6; 0 100.0; 9 100.0
		NA - All 1092.1; 0 100.0; 9 100.0
		
		EVA001_S66
		MP - All 2198.8; 0 100.0; 9 100.0
		PA - All 527.4; 0 100.0; 9 100.0
		HA - All 1449.3; 0 100.0; 9 100.0
		PB1 - All 618.4; 0 100.0; 9 100.0
		NP - All 439.3; 0 100.0; 9 100.0
		NS - All 1214.7; 0 100.0; 9 100.0
		PB2 - All 690.6; 0 100.0; 9 100.0
		NA - All 1092.1; 0 100.0; 9 100.0
		"""

		consensus_EVA001_S66 = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.FILES_EVA001_S66_consensus)
		self.assertTrue(os.path.exists(consensus_EVA001_S66))
		
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
		coverage.add_coverage('PB1', Coverage.COVERAGE_MORE_9, '100.0')
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

		sample_name = 'xpto' 
		utils = Utils()
		temp_dir = utils.get_temp_dir()
		result_fasta = utils.filter_fasta_all_sequences(consensus_EVA001_S66, sample_name, coverage, temp_dir)
		self.assertTrue(result_fasta != None)
		self.assertTrue(os.path.exists(result_fasta))
		locus_fasta = self.utils.is_fasta(result_fasta)
		self.assertEquals(8, locus_fasta)
		self.assertEquals(8, self.utils.get_max_length_fasta(result_fasta))
		
		coverage.add_coverage('NA', Coverage.COVERAGE_MORE_9, '99.9')
		result_fasta = utils.filter_fasta_all_sequences(consensus_EVA001_S66, sample_name, coverage, temp_dir)
		self.assertTrue(result_fasta == None)

		## remove dir
		self.utils.remove_dir(temp_dir)


	def test_filter_by_sequence_sequences(self):
		"""
		test filter fasta
		"""
		
		"""
		EVA011_S54
		MP - All 2198.8; 0 100.0; 9 100.0
		PA - All 527.4; 0 100.0; 9 100.0
		HA - All 1449.3; 0 100.0; 9 100.0
		PB1 - All 618.4; 0 100.0; 9 100.0
		NP - All 439.3; 0 100.0; 9 100.0
		NS - All 1214.7; 0 100.0; 9 100.0
		PB2 - All 690.6; 0 100.0; 9 100.0
		NA - All 1092.1; 0 100.0; 9 100.0
		
		EVA003_S91
		MP - All 2198.8; 0 100.0; 9 100.0
		PA - All 527.4; 0 100.0; 9 100.0
		HA - All 1449.3; 0 100.0; 9 100.0
		PB1 - All 618.4; 0 100.0; 9 100.0
		NP - All 439.3; 0 100.0; 9 100.0
		NS - All 1214.7; 0 100.0; 9 100.0
		PB2 - All 690.6; 0 100.0; 9 100.0
		NA - All 1092.1; 0 100.0; 9 100.0
		
		EVA002_S52
		MP - All 2198.8; 0 100.0; 9 100.0
		PA - All 527.4; 0 100.0; 9 100.0
		HA - All 1449.3; 0 100.0; 9 100.0
		PB1 - All 618.4; 0 100.0; 9 100.0
		NP - All 439.3; 0 100.0; 9 100.0
		NS - All 1214.7; 0 100.0; 9 100.0
		PB2 - All 690.6; 0 100.0; 9 100.0
		NA - All 1092.1; 0 100.0; 9 100.0
		
		EVA001_S66
		MP - All 2198.8; 0 100.0; 9 100.0
		PA - All 527.4; 0 100.0; 9 100.0
		HA - All 1449.3; 0 100.0; 9 100.0
		PB1 - All 618.4; 0 100.0; 9 100.0
		NP - All 439.3; 0 100.0; 9 100.0
		NS - All 1214.7; 0 100.0; 9 100.0
		PB2 - All 690.6; 0 100.0; 9 100.0
		NA - All 1092.1; 0 100.0; 9 100.0
		"""

		consensus_EVA001_S66 = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.FILES_EVA001_S66_consensus)
		self.assertTrue(os.path.exists(consensus_EVA001_S66))
		
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
		coverage.add_coverage('PB1', Coverage.COVERAGE_MORE_9, '100.0')
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

		sample_name = 'xpto' 
		utils = Utils()
		temp_dir = utils.get_temp_dir()
		result_fasta = utils.filter_fasta_by_sequence_names(consensus_EVA001_S66, sample_name, 'NP', coverage, temp_dir)
		self.assertTrue(result_fasta != None)
		self.assertTrue(os.path.exists(result_fasta))
		locus_fasta = self.utils.is_fasta(result_fasta)
		self.assertEquals(1, locus_fasta)
		
		
		coverage.add_coverage('NA', Coverage.COVERAGE_MORE_9, '99.9')
		result_fasta = utils.filter_fasta_by_sequence_names(consensus_EVA001_S66, sample_name, 'NA', coverage, temp_dir)
		self.assertTrue(result_fasta == None)
		
		coverage.add_coverage('PB2', Coverage.COVERAGE_MORE_9, '19.9')
		result_fasta = utils.filter_fasta_by_sequence_names(consensus_EVA001_S66, sample_name, 'PB2', coverage, temp_dir)
		self.assertTrue(result_fasta == None)
		
		result_fasta = utils.filter_fasta_by_sequence_names(consensus_EVA001_S66, sample_name, 'PA', coverage, temp_dir)
		self.assertTrue(result_fasta != None)
		self.assertTrue(os.path.exists(result_fasta))
		locus_fasta = self.utils.is_fasta(result_fasta)
		self.assertEquals(1, locus_fasta)

		result_fasta = utils.filter_fasta_by_sequence_names(consensus_EVA001_S66, sample_name, 'PA', None, temp_dir)
		self.assertTrue(result_fasta != None)
		self.assertTrue(os.path.exists(result_fasta))
		locus_fasta = self.utils.is_fasta(result_fasta)
		self.assertEquals(1, locus_fasta)

		## remove dir
		self.utils.remove_dir(temp_dir)


	def test_clean_fasta_names(self):
		"""
		test clean name
		"""
		vect_names_to_clean = [ FileExtensions.FILE_FASTA, FileExtensions.FILE_FA, FileExtensions.FILE_CONSENSUS_FASTA]
		temp_in_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_GLOBAL_PROJECT,\
					ConstantsTestsCase.FILE_TO_CLEAN_NAME)
		expect_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_GLOBAL_PROJECT,\
					ConstantsTestsCase.FILE_RESULT_TO_CLEAN_NAME)

		temp_out_file = self.utils.get_temp_file('clean_name_test', FileExtensions.FILE_CONSENSUS_FASTA)
		self.utils.clean_fasta_names(vect_names_to_clean, temp_in_file, temp_out_file)
		self.assertTrue(filecmp.cmp(expect_file, temp_out_file))
		os.unlink(temp_out_file)

	def test_clean_extension(self):
		"""
		clean extension in file names
		"""
		self.assertEquals('sdfsdf', self.utils.clean_extension('sdfsdf.as'))
		self.assertEquals('sdfsdf', self.utils.clean_extension('sdfsdf'))
		self.assertEquals('sdfsdf.as', self.utils.clean_extension('sdfsdf.as.fasta'))
		self.assertEquals('', self.utils.clean_extension(''))


	def test_str2bool(self):
		self.assertTrue(self.utils.str2bool('true'))
		self.assertTrue(self.utils.str2bool('True'))
		self.assertTrue(self.utils.str2bool('yes'))
		self.assertTrue(self.utils.str2bool('t'))
		self.assertTrue(self.utils.str2bool('y'))
		self.assertTrue(self.utils.str2bool('1'))
		self.assertFalse(self.utils.str2bool('false'))
		self.assertFalse(self.utils.str2bool('0'))
		
	def test_read_file_to_string(self):
		
		self.assertEquals(None, self.utils.read_file_to_string("/root/fiels.fasta"))
		self.assertEquals("(EVA003_S91:0.0,EVA001_S66:0.0,EVA002_S52:0.0,EVA011_S54:0.0);\n",\
				self.utils.read_file_to_string(os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT,\
				ConstantsTestsCase.FILE_FASTTREE_RESULT_NWK)))


	def test_get_sequence_from_genbank(self):
		
		genbank_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		self.assertTrue(os.path.join(genbank_file))
		
		gene = Gene('PB2', 10, 20, 1)
		seq = self.utils.get_sequence_from_genbank('PB2', gene, genbank_file)
		self.assertEquals('ATGTCGCAGTCTCGCACTCGCGAGATACTGACAAAAACCACAGTGGACCATATGGCCATAATTAAGAAGTACACATCGGGGAGACAGGAAAAGAACCCGTCACTTAGGATGAAATGGATGATGGCAATGAAATATCCAATCACTG' +\
						'CTGACAAAAGGGTAACAGAAATGGTTCCGGAGAGAAATGAACAAGGACAAACTCTATGGAGCAAAATGAGTGATGCTGGATCAGATAGAGTGATGGTATCACCTTTGGCTGTAACATGGTGGAATAGGAATGGACCCGTGACAAGTA' +\
						'CGGTCCATTACCCAAAAGTGTACAAAACTTATTTCGACAAAGTCGAAAGGTTAAAACATGGAACCTTTGGCCCTGTCCATTTTAGAAATCAAGTCAAGATACGCAGAAGAGTAGACATAAACCCTGGTCATGCAGACCTCAGTGCCA' +\
						'AAGAGGCACAAGATGTAATTATGGAAGTTGTTTTTCCCAATGAAGTGGGAGCCAGAATACTAACATCAGAATCACAACTAACAATAACTAAAGAGAAAAAAGAAGAACTCCGAGATTGCAAAATTTCTCCCTTGATGGTCGCATACA' +\
						'TGTTAGAGAGAGAACTTGTGCGGAAAACAAGATTTCTCCCAGTTGCTGGCGGAACAAGCAGTATATACATTGAAGTTTTACATTTGACTCAAGGAACGTGTTGGGAACAAATGTACACTCCAGGTGGAGGAGTGAGGAATGACGATG' +\
						'TTGACCAAAGCCTAATTATTGCGGCCAGGAACATAGTAAGAAGAGCCGCAGTATCAGCAGATCCATTAGCATCTTTATTGGAGATGTGCCACAGCACGCAAATTGGCGGAACAAGGATGGTGGACATTCTTAGACAGAACCCGACTG' +\
						'AAGAACAAGCTGTGGATATATGCAAGGCTGCAATGGGATTGAGAATCAGCTCATCCTTCAGCTTTGGTGGCTTTACATTTAAAAGAACAAGCGGGTCGTCAGTCAAAAAAGAAGAAGAGGTGCTTACAGGCAATCTCCAAACATTGA' +\
						'GAATAAGAGTACATGAGGGGTATGAGGAGTTCACAATGGTGGGGAAAAGAGCAACAGCTATACTAAGAAAAGCAACCAGAAGATTGGTTCAACTCATAGTGAGTGGAAGAGACGAACAGTCAATAGCCGAAGCAATAATCGTGGCCA' +\
						'TGGTGTTTTCACAAGAAGATTGCATGATAAAAGCAGTTAGAGGTGACCTGAATTTTGTCAACAGAGCAAATCAGCGGTTGAACCCCATGCATCAGCTTTTAAGGCATTTTCAGAAAGATGCGAAAGTACTCTTTCAAAATTGGGGAG' +\
						'TTGAACACATCGACAGTGTGATGGGAATGGTTGGAGTATTACCAGATATGACTCCAAGCACAGAGATGTCAATGAGAGGAATAAGAGTCAGCAAAATGGGTGTGGATGAATACTCCAGTACAGAGAGGGTGGTGGTTAGCATTGATC' +\
						'GGTTTTTGAGAGTTCGAGACCAACGTGGGAATGTATTATTATCTCCTGAGGAGGTTAGTGAAACACAGGGAACTGAGAGACTGACAATAACTTATTCATCGTCGATGATGTGGGAGATTAACGGTCCTGAGTCAGTCTTGGTCAATA' +\
						'CCTATCAATGGATCATCAGGAATTGGGAAGCTGTTAAAATTCAATGGTCTCAGAATCCTGCAATGTTGTACAACAAAATGGAATTTGAACCATTTCAATCTTTAGTCCCCAAGGCCATTAGAAGCCAATACAGTGGGTTTGTCAGAA' +\
						'CTCTATTCCAACAAATGAGAGACGTACTTGGGACATTTGACACTGCCCAGATAATAAAGCTTCTCCCTTTTGCAGCTGCTCCACCGAAGCAAAGCAGAATGCAGTTCTCTTCACTGACTGTGAATGTGAGGGGATCAGGGATGAGAA' +\
						'TACTTGTAAGGGGCAATTCTCCTGTATTCAACTACAACAAGACCACTAAAAGGCTAACAATTCTCGGAAAAGATGCCGGCACTTTAATTGAAGACCCAGATGAAAGCACATCCGGAGTGGAGTCCGCCGTCTTGAGAGGGTTCCTCA' +\
						'TTATAGGTAAAGAAGACAGAAGATACGGACCTGCATTAAGCATCAATGAACTGAGTAACCTTGCAAAAGGAGAAAAGGCTAATGTGCTAATTGGGCAAGGAGACGTGGTGTTGGTAATGAAACGAAAACGGGACTCTAGTATACTTA' +\
						'CTGACAGCCAGACAGCGACCAAAAGAATTCGGATGGCCATCAATTAA', str(seq))
		self.assertEquals('MSQSRTREILTKTTVDHMAIIKKYTSGRQEKNPSLRMKWMMAMK' +\
                    'YPITADKRVTEMVPERNEQGQTLWSKMSDAGSDRVMVSPLAVTWWNRNGPVTSTVHYP' +\
                    'KVYKTYFDKVERLKHGTFGPVHFRNQVKIRRRVDINPGHADLSAKEAQDVIMEVVFPN' +\
                    'EVGARILTSESQLTITKEKKEELRDCKISPLMVAYMLERELVRKTRFLPVAGGTSSIY' +\
                    'IEVLHLTQGTCWEQMYTPGGGVRNDDVDQSLIIAARNIVRRAAVSADPLASLLEMCHS' +\
                    'TQIGGTRMVDILRQNPTEEQAVDICKAAMGLRISSSFSFGGFTFKRTSGSSVKKEEEV' +\
                    'LTGNLQTLRIRVHEGYEEFTMVGKRATAILRKATRRLVQLIVSGRDEQSIAEAIIVAM' +\
                    'VFSQEDCMIKAVRGDLNFVNRANQRLNPMHQLLRHFQKDAKVLFQNWGVEHIDSVMGM' +\
                    'VGVLPDMTPSTEMSMRGIRVSKMGVDEYSSTERVVVSIDRFLRVRDQRGNVLLSPEEV' +\
                    'SETQGTERLTITYSSSMMWEINGPESVLVNTYQWIIRNWEAVKIQWSQNPAMLYNKME' +\
                    'FEPFQSLVPKAIRSQYSGFVRTLFQQMRDVLGTFDTAQIIKLLPFAAAPPKQSRMQFS' +\
                    'SLTVNVRGSGMRILVRGNSPVFNYNKTTKRLTILGKDAGTLIEDPDESTSGVESAVLR' +\
                    'GFLIIGKEDRRYGPALSINELSNLAKGEKANVLIGQGDVVLVMKRKRDSSILTDSQTA' +\
                    'TKRIRMAIN*', str(seq.translate(table=Constants.TRANSLATE_TABLE_NUMBER)))
		
		
	def test_get_sequence_from_genbank_2(self):
		
		genbank_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_TWO_GENES_JOINED_GBK)
		self.assertTrue(os.path.join(genbank_file))
		
		gene = Gene('PB2', 10, 20, 1)
		seq = self.utils.get_sequence_from_genbank('PB2', gene, genbank_file)
		self.assertEquals('MSQSRTREILTKTTVDHMAIIKKYTSGRQEKNPSLRMKWMMAMK' +\
                    'YPITADKRVTEMVPERNEQGQTLWSKMSDAGSDRVMVSPLAVTWWNRNGPVTSTVHYP' +\
                    'KVYKTYFDKVERLKHGTFGPVHFRNQVKIRRRVDINPGHADLSAKEAQDVIMEVVFPN' +\
                    'EVGARILTSESQLTITKEKKEELRDCKISPLMVAYMLERELVRKTRFLPVAGGTSSIY' +\
                    'IEVLHLTQGTCWEQMYTPGGGVRNDDVDQSLIIAARNIVRRAAVSADPLASLLEMCHS' +\
                    'TQIGGTRMVDILRQNPTEEQAVDICKAAMGLRISSSFSFGGFTFKRTSGSSVKKEEEV' +\
                    'LTGNLQTLRIRVHEGYEEFTMVGKRATAILRKATRRLVQLIVSGRDEQSIAEAIIVAM' +\
                    'VFSQEDCMIKAVRGDLNFVNRANQRLNPMHQLLRHFQKDAKVLFQNWGVEHIDSVMGM' +\
                    'VGVLPDMTPSTEMSMRGIRVSKMGVDEYSSTERVVVSIDRFLRVRDQRGNVLLSPEEV' +\
                    'SETQGTERLTITYSSSMMWEINGPESVLVNTYQWIIRNWEAVKIQWSQNPAMLYNKME' +\
                    'FEPFQSLVPKAIRSQYSGFVRTLFQQMRDVLGTFDTAQIIKLLPFAAAPPKQSRMQFS' +\
                    'SLTVNVRGSGMRILVRGNSPVFNYNKTTKRLTILGKDAGTLIEDPDESTSGVESAVLR' +\
                    'GFLIIGKEDRRYGPALSINELSNLAKGEKANVLIGQGDVVLVMKRKRDSSILTDSQTA' +\
                    'TKRIRMAIN*', str(seq.translate(table=Constants.TRANSLATE_TABLE_NUMBER)))
		
		gene = Gene('PB1', 10, 20, 1)
		seq = self.utils.get_sequence_from_genbank('PB2', gene, genbank_file)
		self.assertEquals('MDVNPTLLFLKVPAQNAISTTFPYTGDPPYSHGTGTGYTMDTVN' +\
                    'RTHQYSERGKWTTNTETGAPQLNPIDGPLPEDNEPSGYAQTDCVLEAMAFLEESHPGI' +\
                    'FENSCLETMEAVQQTRVDKLTQGRQTYDWTLNRNQPAATALANTIEVFRTNGLTANES' +\
                    'GRLIDYLKDVMESMDKEEMEITTHFQRKRRVRDNMTKKMVTQRTIGKKKQRVNKRGYL' +\
                    'IRALTLNTMTKDAERGKLKRRAIATPGMQIRGFVYFVETLARSICEKLEQSGLPVGGN' +\
                    'EKKAKLANVVRKMMTNSQDTELSFTITGDNTKWNENQNPRMFLAMITYITKNQPEWFR' +\
                    'NILSIAPIMFSNKMARLGKGYMFESKRMKLRTQIPAEMLASIDLKYFNESTRKKIEKI' +\
                    'RPLLIDGTASLSPGMMMGMFNMLSTVLGVSILNLGQKKYTKTTYWWDGLQSSDDFALI' +\
                    'VNAPNHEGIQAGVDRFYRTCKLVGINMSKKKSYINKTGTFEFTSFFYRYGFVANFSME' +\
                    'LPSFGVSGINESADMSIGVTVIKNNMINNDLGPATAQMALQLFIKDYRYTYRCHRGDT' +\
                    'QIQTRRSFEIKKLWDQTQSRTGLLVSDGGPNLYNIRNLHIPEVCLKWELMDENYRGRL' +\
                    'CNPLNPFVSHKEIESVNNAVVMPAHGPAKSMEYDAVATTHSWIPKRNRSILNTSQRGI' +\
                    'LEDEQMYQKCCNLFEKFFPSSSYRRPIGISSMVEAMVSRARIDARIDFESGRIKKEEF' +\
                    'SEIMKICSTIEELRRQK*', str(seq.translate(table=Constants.TRANSLATE_TABLE_NUMBER)))


	def test_get_sequence_from_genbank_reversed(self):
		"""
		If the sequence is complement the Bio.Seq() return in the right order
		"""
		genbank_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_pb2_reversed)
		self.assertTrue(os.path.join(genbank_file))
		
		gene = Gene('PB2', 10, 20, 1)
		seq = self.utils.get_sequence_from_genbank('PB2_reversed', gene, genbank_file)
		self.assertEquals('MSQSRTREILTKTTVDHMAIIKKYTSGRQEKNPSLRMKWMMAMK' +\
                    'YPITADKRVTEMVPERNEQGQTLWSKMSDAGSDRVMVSPLAVTWWNRNGPVTSTVHYP' +\
                    'KVYKTYFDKVERLKHGTFGPVHFRNQVKIRRRVDINPGHADLSAKEAQDVIMEVVFPN' +\
                    'EVGARILTSESQLTITKEKKEELRDCKISPLMVAYMLERELVRKTRFLPVAGGTSSIY' +\
                    'IEVLHLTQGTCWEQMYTPGGGVRNDDVDQSLIIAARNIVRRAAVSADPLASLLEMCHS' +\
                    'TQIGGTRMVDILRQNPTEEQAVDICKAAMGLRISSSFSFGGFTFKRTSGSSVKKEEEV' +\
                    'LTGNLQTLRIRVHEGYEEFTMVGKRATAILRKATRRLVQLIVSGRDEQSIAEAIIVAM' +\
                    'VFSQEDCMIKAVRGDLNFVNRANQRLNPMHQLLRHFQKDAKVLFQNWGVEHIDSVMGM' +\
                    'VGVLPDMTPSTEMSMRGIRVSKMGVDEYSSTERVVVSIDRFLRVRDQRGNVLLSPEEV' +\
                    'SETQGTERLTITYSSSMMWEINGPESVLVNTYQWIIRNWEAVKIQWSQNPAMLYNKME' +\
                    'FEPFQSLVPKAIRSQYSGFVRTLFQQMRDVLGTFDTAQIIKLLPFAAAPPKQSRMQFS' +\
                    'SLTVNVRGSGMRILVRGNSPVFNYNKTTKRLTILGKDAGTLIEDPDESTSGVESAVLR' +\
                    'GFLIIGKEDRRYGPALSINELSNLAKGEKANVLIGQGDVVLVMKRKRDSSILTDSQTA' +\
                    'TKRIRMAIN*', str(seq.translate(table=Constants.TRANSLATE_TABLE_NUMBER)))

	def test_validate_date(self):
		
		try:
			self.utils.validate_date('12/12/2017')
		except ValueError as e:
			self.fail('Real date')
		
		try:
			self.utils.validate_date('12-12-2017')
		except ValueError as e:
			self.fail('Real date')

		try:
			self.utils.validate_date('12.12.2017')
			self.fail('must throw exception')
		except ValueError as e:
			pass
		
		try:
			self.utils.validate_date('2017-12-12')
			self.fail('must throw exception')
		except ValueError as e:
			pass
		
		try:
			self.utils.validate_date('12-21-2017')
			self.fail('must throw exception')
		except ValueError as e:
			pass


	def test_clean_fasta_file(self):
		"""
		test clean file
		"""
		to_clean_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "to_clean_sequences.fasta")
		cleanned_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "cleanned_sequences.fasta")
		self.assertTrue(os.path.exists(to_clean_file))
		self.assertTrue(os.path.exists(cleanned_file))
		out_file = self.utils.get_temp_file('clean_sequences', '.fasta')
		self.utils.clean_fasta_file(to_clean_file, out_file)
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, cleanned_file))
		os.unlink(out_file)

	def test_from_genbank_to_bed(self):
		"""
		test create bed from genbank
		"""
		
		file_expected = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_BED)
		file_in = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		file_out = self.utils.get_temp_file("test_genbank_to_bed", FileExtensions.FILE_BED)
		self.utils.from_genbank_to_bed(file_in, file_out)
		self.assertTrue(os.path.exists(file_out))
		self.assertTrue(os.path.getsize(file_out) > 10)
		self.assertTrue(filecmp.cmp(file_expected, file_out))
		
		### create .tbi
		self.software.create_index_files_from_igv_tools(file_out)
		self.assertTrue(os.path.exists(file_out + FileExtensions.FILE_IDX))
		self.assertTrue(os.path.getsize(file_out + FileExtensions.FILE_IDX) > 100)
		os.unlink(file_out + FileExtensions.FILE_IDX)
		
		file_expected = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_TWO_GENES_JOINED_BED)
		file_in = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_TWO_GENES_JOINED_GBK)
		self.utils.from_genbank_to_bed(file_in, file_out)
		self.assertTrue(os.path.exists(file_out))
		self.assertTrue(os.path.getsize(file_out) > 10)
		self.assertTrue(filecmp.cmp(file_expected, file_out))
		os.unlink(file_out)
	
		file_in =  os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_PV1_KX162693_gbk)
		file_expected = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_PV1_KX162693_bed)
		self.utils.from_genbank_to_bed(file_in, file_out)
		self.assertTrue(os.path.exists(file_out))
		self.assertTrue(os.path.getsize(file_out) > 10)
		self.assertTrue(filecmp.cmp(file_expected, file_out))
		os.unlink(file_out)


	def test_fasta_files(self):
		"""
		test fasta files to degenerated bases
		"""
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTA, "fasta_1.fasta")
		self.assertTrue(os.path.exists(txt_file))
		try:
			self.assertFalse(self.utils.has_degenerated_bases(txt_file))
		except Exception as e:
			self.fail("must pass")
		
		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTA, "fasta_fail_1.fasta")
		self.assertTrue(os.path.exists(txt_file))
		try:
			self.utils.has_degenerated_bases(txt_file)
			self.fail("must raise error")
		except Exception as e:
			self.assertEquals(e.args[0], "Error: file is not in FASTA format.")

		txt_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTA, "fasta_fail_2.fasta")
		self.assertTrue(os.path.exists(txt_file))
		try:
			self.utils.has_degenerated_bases(txt_file)
			self.fail("must raise error")
		except Exception as e:
			self.assertEquals(e.args[0], "Error: sequence 'EVA011_S54' must have only 'A', 'C', 'T' and 'G' bases.")

	def test_clean_name(self):
		"""
		clean name
		"""
		self.assertEquals('cpto_', self.utils.clean_name('cpto ()', { ' ' : '_', '(' : '' , ')' : '' }))
	
	
	def test_clean_genbank_version_name(self):
		"""
		clean genbank
		"""
		genbank_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_PV1_KX162693_gbk)
		expect_genbank_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_PV1_KX162693_clean_version_gbk)
		out_file = self.utils.get_temp_file("file_name", ".gb")
		self.utils.clean_genbank_version_name(genbank_file, out_file)
		self.assertTrue(filecmp.cmp(out_file, expect_genbank_file))
		os.unlink(out_file)
		
	def test_clean_genbank_version_name_2(self):
		"""
		clean genbank
		"""
		genbank_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_PV_Sabin1_V01150_gbk)
		expect_genbank_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_PV_Sabin1_V01150_clean_version_gbk)
		out_file = self.utils.get_temp_file("file_name", ".gb")
		self.utils.clean_genbank_version_name(genbank_file, out_file)
		self.assertTrue(filecmp.cmp(out_file, expect_genbank_file))
		os.unlink(out_file)
		
		
		
	def test_get_number_seqs_names_bigger_than(self):
		
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_PV1_KX162693_spades_out_fasta)
		
		n_count = self.utils.get_number_seqs_names_bigger_than(fasta_file, Constants.MAX_LENGTH_SEQ_NAME)
		self.assertEquals(0, n_count)
		
		n_count = self.utils.get_number_seqs_names_bigger_than(fasta_file, 10)
		self.assertEquals(7, n_count)
		
		n_count = self.utils.get_number_seqs_names_bigger_than(fasta_file, 15)
		self.assertEquals(7, n_count)
		
		n_count = self.utils.get_number_seqs_names_bigger_than(fasta_file, 17)
		self.assertEquals(6, n_count)



