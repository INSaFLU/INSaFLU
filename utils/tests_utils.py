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
from constants.constants import FileExtensions
from Bio import SeqIO

class Test(unittest.TestCase):


	### static
	constantsTestsCase = ConstantsTestsCase()
	utils = Utils()

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
			self.assertEqual("File is not in fastq.gz format", e.args[0])
	
	def test_read_text_file(self):
		"""
		"""
		txt_file = os.path.join(self.baseDirectory, 'abricate', ConstantsTestsCase.MANAGING_TEST_ABRICATE + ".xpto")
		self.assertFalse(os.path.exists(txt_file))
		try:
			self.utils.read_text_file(txt_file)
			self.fail("must raise error")
		except IOError as e:
			self.assertEqual("Error: file '/home/mmp/eclipse_oxygen/fluwebvirus/static_all/tests/abricate/abricate_out.txt.xpto' doens't exist.", e.args[0])
		
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
		self.utils.create_index_files(SoftwareNames.SOFTWARE_TABIX, temp_file + ".gz")
		self.assertTrue(os.path.exists(temp_file + ".gz.tbi"))
		self.utils.remove_temp_file(temp_file + ".gz.tbi")
		try:
			self.utils.create_index_files(SoftwareNames.SOFTWARE_TABIX, temp_file)
			self.fail("must raise exception")
		except Exception as e:
			self.assertEquals("Fail to create index", e.args[0])
		os.unlink(test_file_to_remove)
		print("These two faults are expected: 'Not a BGZF file' and 'tbx_index_build failed'")
		self.utils.remove_dir(temp_path)
	

	def test_add_vcf_FREQ(self):
		vcf_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "temp.vcf")
		vcf_file_expected = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "temp_more_REF.vcf")
		temp_dir = self.utils.get_temp_dir()
		temp_vcf_file = os.path.join(temp_dir, os.path.basename(vcf_file))
		
		self.utils.copy_file(vcf_file, temp_vcf_file)
		
		software = Software()
		software.test_bgzip_and_tbi_in_vcf(temp_vcf_file)
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
		self.assertEqual(3, count_hits.get_hits_50_90())
		self.assertEqual(7, count_hits.get_total())

	def test_count_hits_from_tab_2(self):
		"""
		test cout hits
		"""
		tab_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "run_snippy2_other.tab")
		vect_count_type = ['snp']
		count_hits = self.utils.count_hits_from_tab(tab_file, vect_count_type)
		self.assertEqual(3, count_hits.get_hits_less_50())
		self.assertEqual(0, count_hits.get_hits_50_90())
		self.assertEqual(3, count_hits.get_total())
		
	def test_count_hits_from_tab_3(self):
		"""
		test cout hits
		"""
		tab_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "run_snippy2_other.tab")
		vect_count_type = ['snp', 'del']
		count_hits = self.utils.count_hits_from_tab(tab_file, vect_count_type)
		self.assertEqual(4, count_hits.get_hits_less_50())
		self.assertEqual(0, count_hits.get_hits_50_90())
		self.assertEqual(4, count_hits.get_total())
		
	def test_get_variations_by_freq_from_tab(self):
		"""
		test get variations by freq
		"""
		tab_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.DIR_VCF, "resutl_vcf_to_tab2.tab")
		vect_count_type = ['snp']
		(dict_less_50, dict_more_50) = self.utils.get_variations_by_freq_from_tab(tab_file, vect_count_type)
		self.assertEqual(2, len(dict_less_50))
		self.assertEqual(817, dict_less_50['NS'][0])
		self.assertEqual(237, dict_less_50['PB2'][0])
		self.assertEqual(774, dict_less_50['PB2'][1])
		self.assertEqual(896, dict_less_50['PB2'][2])
		self.assertEqual(8, len(dict_more_50))
		self.assertEqual(324, dict_more_50['MP'][0])
		self.assertEqual(99, dict_more_50['NS'][0])
		self.assertEqual(117, dict_more_50['PA'][0])

	def test_get_elements_and_genes(self):
		"""
		test cout hits
		"""
		genbank_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS,\
					ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		dt_data = self.utils.get_elements_and_genes(genbank_file)
		self.assertEquals(8, len(dt_data))
		self.assertTrue('NS' in dt_data)
		self.assertFalse('xptoNS' in dt_data)
		self.assertTrue('PB2' in dt_data)
		self.assertEquals(30, dt_data['PB2'][0][0])
		self.assertEquals(2280, dt_data['PB2'][0][1])
		self.assertEquals('PB2', dt_data['PB2'][0][2])
		self.assertEquals(1, dt_data['PB2'][0][3])

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
		
		
		coverage.add_coverage('NA', Coverage.COVERAGE_MORE_9, '99.9')
		result_fasta = utils.filter_fasta_all_sequences(consensus_EVA001_S66, sample_name, coverage, temp_dir)
		self.assertTrue(result_fasta != None)
		self.assertTrue(os.path.exists(result_fasta))
		locus_fasta = self.utils.is_fasta(result_fasta)
		self.assertEquals(7, locus_fasta)
		record_dict = SeqIO.to_dict(SeqIO.parse(result_fasta, "fasta"))
		self.assertFalse('NA' in record_dict)
		self.assertTrue('PB2' in record_dict)
		
		coverage.add_coverage('PB2', Coverage.COVERAGE_MORE_9, '19.9')
		result_fasta = utils.filter_fasta_all_sequences(consensus_EVA001_S66, sample_name, coverage, temp_dir)
		self.assertTrue(result_fasta != None)
		self.assertTrue(os.path.exists(result_fasta))
		locus_fasta = self.utils.is_fasta(result_fasta)
		self.assertEquals(6, locus_fasta)
		record_dict = SeqIO.to_dict(SeqIO.parse(result_fasta, "fasta"))
		self.assertFalse('NA' in record_dict)
		self.assertFalse('PB2' in record_dict)
		self.assertTrue('NP' in record_dict)

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




