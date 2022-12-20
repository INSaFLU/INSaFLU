'''
Created on Oct 28, 2017

@author: mmp
'''
import unittest, time, os, filecmp
from constants.constantsTestsCase import ConstantsTestsCase
from django.test.utils import override_settings
from utils.utils import Utils
from utils.software import Software
from django.conf import settings
from constants.software_names import SoftwareNames
from utils.result import Coverage, FeatureLocationSimple
from constants.constants import Constants, FileExtensions
from managing_files.models import ProjectSample
from utils.result import Gene
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
		self.assertEqual((True, Constants.FORMAT_FASTQ_illumina), self.utils.is_fastq_gz(path_file))
		
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_ABRICATE, ConstantsTestsCase.MANAGING_TEST_ABRICATE)
		try:
			self.assertEqual((True, Constants.FORMAT_FASTQ_illumina), self.utils.is_fastq_gz(path_file))
		except Exception as e:
			self.assertEqual("File need to have suffix '.fastq.gz'/'.fq.gz'", e.args[0])
			
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		try:
			self.assertEqual((True, Constants.FORMAT_FASTQ_illumina), self.utils.is_fastq_gz(path_file))
		except Exception as e:
			self.assertEqual("File need to have suffix '.fastq.gz'/'.fq.gz'", e.args[0])
			
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA_FAKE_GZ)
		try:
			self.assertEqual((True, Constants.FORMAT_FASTQ_ont), self.utils.is_fastq_gz(path_file))
		except Exception as e:
			self.assertEqual("File is not in fastq.gz format.", e.args[0])
			
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ10_1_test)
		try:
			self.assertEqual((True, Constants.FORMAT_FASTQ_illumina), self.utils.is_fastq_gz(path_file))
		except Exception as e:
			self.assertEqual("File is not in fastq.gz format.", e.args[0])
		
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_nanopore)
		try:
			self.assertEqual((True, Constants.FORMAT_FASTQ_ont), self.utils.is_fastq_gz(path_file))
			self.assertEqual(True, self.utils.is_fastq_gz(path_file)[0])
			self.assertEqual(Constants.FORMAT_FASTQ_ont, self.utils.is_fastq_gz(path_file)[1])
		except Exception as e:
			self.assertEqual("Can not detect file format. Ensure Illumina fastq file.", e.args[0])
	
		### IONtorrent
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_ion_torrent_1)
		self.assertTrue(os.path.exists(path_file))
		try:
			self.assertEqual((True, Constants.FORMAT_FASTQ_illumina), self.utils.is_fastq_gz(path_file))
		except Exception as e:
			self.assertEqual("File is not in fastq.gz format.", e.args[0])
			
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_ion_torrent_2)
		self.assertTrue(os.path.exists(path_file))
		try:
			self.assertEqual((True, Constants.FORMAT_FASTQ_illumina), self.utils.is_fastq_gz(path_file))
		except Exception as e:
			self.assertEqual("File is not in fastq.gz format.", e.args[0])
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_ion_torrent_3)
		self.assertTrue(os.path.exists(path_file))
		try:
			self.assertEqual((True, Constants.FORMAT_FASTQ_illumina), self.utils.is_fastq_gz(path_file))
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
		self.assertEquals(51, vect_out[5].find("Victoria"))
		
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
		
	def test_get_elements_and_genes_2(self):
		"""
		test cout hits
		"""
		genbank_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS,\
					ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_COVID_GBK)
		geneticElement = self.utils.get_elements_and_genes(genbank_file)
		self.assertEquals(1, len(geneticElement.get_sorted_elements()))
		self.assertEquals(265, geneticElement.get_genes('MN908947')[0].start)
		self.assertEquals(21555, geneticElement.get_genes('MN908947')[0].end)
		self.assertEquals(21562, geneticElement.get_genes('MN908947')[1].start)
		self.assertEquals(25384, geneticElement.get_genes('MN908947')[1].end)
		for gene in geneticElement.get_genes('MN908947'):
			if (gene.name == 'orf1ab'):
				self.assertEqual(2, len(gene.get_feature_locations()))
				self.assertEqual(FeatureLocationSimple(265, 13468, 1), gene.get_feature_locations()[0])
				self.assertEqual(FeatureLocationSimple(13467, 21555, 1), gene.get_feature_locations()[1])
			else:
				self.assertEqual(0, len(gene.get_feature_locations()))
		
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
		result_fasta = utils.filter_fasta_all_sequences(consensus_EVA001_S66, sample_name, coverage, -1, temp_dir)
		self.assertTrue(result_fasta != None)
		self.assertTrue(os.path.exists(result_fasta))
		locus_fasta = self.utils.is_fasta(result_fasta)
		self.assertEquals(8, locus_fasta)
		self.assertEquals(8, self.utils.get_max_length_fasta(result_fasta))
		
		coverage.add_coverage('NA', Coverage.COVERAGE_MORE_9, '99.9')
		result_fasta = utils.filter_fasta_all_sequences(consensus_EVA001_S66, sample_name, coverage, -1, temp_dir)
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
		result_fasta = utils.filter_fasta_by_sequence_names(consensus_EVA001_S66, sample_name, 'NP', coverage, None, 70, temp_dir)
		self.assertTrue(result_fasta != None)
		self.assertTrue(os.path.exists(result_fasta))
		locus_fasta = self.utils.is_fasta(result_fasta)
		self.assertEquals(1, locus_fasta)
		
		
		coverage.add_coverage('NA', Coverage.COVERAGE_MORE_9, '99.9')
		result_fasta = utils.filter_fasta_by_sequence_names(consensus_EVA001_S66, sample_name, 'NA', coverage, None, 100, temp_dir)
		self.assertTrue(result_fasta == None)
		
		coverage.add_coverage('PB2', Coverage.COVERAGE_MORE_9, '19.9')
		result_fasta = utils.filter_fasta_by_sequence_names(consensus_EVA001_S66, sample_name, 'PB2', coverage, None, 70, temp_dir)
		self.assertTrue(result_fasta == None)
		
		result_fasta = utils.filter_fasta_by_sequence_names(consensus_EVA001_S66, sample_name, 'PA', coverage, None, 70, temp_dir)
		self.assertTrue(result_fasta != None)
		self.assertTrue(os.path.exists(result_fasta))
		locus_fasta = self.utils.is_fasta(result_fasta)
		self.assertEquals(1, locus_fasta)

		result_fasta = utils.filter_fasta_by_sequence_names(consensus_EVA001_S66, sample_name, 'PA', None, None, 70, temp_dir)
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
		
		n_count = self.utils.get_number_seqs_names_bigger_than(fasta_file, Constants.MAX_LENGTH_SEQ_NAME + 12)
		self.assertEquals(0, n_count)
		
		n_count = self.utils.get_number_seqs_names_bigger_than(fasta_file, Constants.MAX_LENGTH_SEQ_NAME)
		self.assertEquals(2, n_count)
		
		n_count = self.utils.get_number_seqs_names_bigger_than(fasta_file, 10)
		self.assertEquals(9, n_count)
		
		n_count = self.utils.get_number_seqs_names_bigger_than(fasta_file, 15)
		self.assertEquals(9, n_count)
		
		n_count = self.utils.get_number_seqs_names_bigger_than(fasta_file, 17)
		self.assertEquals(8, n_count)
		
		n_count = self.utils.get_number_seqs_names_bigger_than(fasta_file, 17, 0)
		self.assertEquals(8, n_count)
		
		n_count = self.utils.get_number_seqs_names_bigger_than(fasta_file, 17, 10)
		self.assertEquals(9, n_count)
		
		n_count = self.utils.get_number_seqs_names_bigger_than(fasta_file, 7, 0)
		self.assertEquals(9, n_count)
		

	def test_sequences_same_length(self):
		
		utils = Utils()

		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_PV1_KX162693_spades_out_fasta)
		self.assertTrue(utils.test_sequences_same_length(fasta_file))
		
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_TEST_LENGTH_fasta_1)
		self.assertTrue(utils.test_sequences_same_length(fasta_file))

		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_TEST_LENGTH_fasta_2)
		self.assertFalse(utils.test_sequences_same_length(fasta_file))
		
	def test_count_gff_seq(self):

		utils = Utils()
		gff_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR_GFF, "insa_flu_temp_empty.gff")
		self.assertTrue(os.path.exists(gff_file))
		self.assertEqual(0, utils.get_number_sequeces_in_gff_file(gff_file))
		
		gff_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR_GFF, "insa_flu_temp_2.gff")
		self.assertTrue(os.path.exists(gff_file))
		self.assertEqual(2, utils.get_number_sequeces_in_gff_file(gff_file))


	def test_remove(self):
		
		utils = Utils()
		main_path = os.path.join(Constants.TEMP_DIRECTORY, Constants.COUNT_DNA_TEMP_DIRECTORY)
		if (not os.path.exists(main_path)): os.makedirs(main_path)
		self.assertTrue(os.path.exists(main_path))
		utils.remove_dir(main_path)
		self.assertTrue(os.path.exists(main_path))
		utils.remove_dir(main_path + "/")
		self.assertTrue(os.path.exists(main_path))
		
		main_path_2 = os.path.join(Constants.TEMP_DIRECTORY, Constants.COUNT_DNA_TEMP_DIRECTORY, "xpto")
		if (not os.path.exists(main_path_2)): os.makedirs(main_path_2)
		self.assertTrue(os.path.exists(main_path_2))
		utils.remove_dir(main_path_2)
		self.assertFalse(os.path.exists(main_path_2))
		self.assertTrue(os.path.exists(main_path))
		
		### test 
		last_value = int(os.path.getmtime(main_path))
		time.sleep(2)
		utils.touch_file(main_path)
		self.assertEqual(2, int(os.path.getmtime(main_path)) - last_value)
		utils.remove_dir(main_path)
		
	def test_link_file(self):
		
		utils = Utils()
		main_path = os.path.join(Constants.TEMP_DIRECTORY, Constants.COUNT_DNA_TEMP_DIRECTORY, "temp_link" )
		utils.remove_dir(main_path)
		if (not os.path.exists(main_path)): os.makedirs(main_path)
		self.assertTrue(os.path.exists(main_path))
		file_1 = os.path.join(main_path, "temp2.txt")
		file_2 = os.path.join(main_path, "temp3.txt")
		
		cmd = "touch {}".format(file_1)
		os.system(cmd)
		self.assertTrue(os.path.exists(file_1))
		self.assertFalse(os.path.exists(file_2))
		utils.link_file(file_1, file_2)
		
		last_value = int(os.path.getmtime(file_2))
		time.sleep(2)
		utils.touch_file(file_2)
		self.assertTrue((int(os.path.getmtime(file_2)) - last_value) in [2,3])
		self.assertTrue(os.path.exists(file_2))
		utils.remove_dir(main_path)


	def test_merge_fasta_first_sequence(self):
		"""
		testing merging fasta
		"""
		
		utils = Utils()
		path_to_fasta = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTA)
		out_file = self.utils.get_temp_file("file_name", ".fasta")
		self.assertTrue(os.path.exists(out_file))
		utils.merge_fasta_first_sequence(path_to_fasta , out_file)
		self.assertTrue(os.path.exists(out_file))
		vect_read_file = utils.read_text_file(out_file)
		self.assertEqual(50, len(vect_read_file))
		
		self.assertTrue(">fasta_1 EVA003_S91" in vect_read_file)
		os.unlink(out_file)
	
	def test_merge_fasta_files(self):
		"""
		testing merging fasta
		"""
		
		utils = Utils()
		path_to_fasta = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTA)
		out_file = self.utils.get_temp_file("file_name", ".fasta")
		self.assertTrue(os.path.exists(out_file))
		vest_data = [[os.path.join(path_to_fasta, "fasta_1.fasta"), "test_sample_1", -1],
					[os.path.join(path_to_fasta, "fasta_2.fasta"), "test_sample_2", -1]]
		utils.merge_fasta_files(vest_data , out_file)
		self.assertTrue(os.path.exists(out_file))
		vect_read_file = utils.read_text_file(out_file)
		self.assertEqual(78, len(vect_read_file))
		
		self.assertTrue(">test_sample_1__EVA003_S91" in vect_read_file)
		self.assertTrue(">test_sample_2__EVA003_S91" in vect_read_file)
		self.assertTrue(">test_sample_1__EVA002_S52_1" in vect_read_file)
		os.unlink(out_file)


	def test_merge_fasta_files_1(self):
		"""
		testing merging fasta
		"""
		
		project_sample = ProjectSample()
		project_sample.save()
		
		project_sample_1 = ProjectSample()
		project_sample_1.save()
		
		utils = Utils()
		path_to_fasta = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTA)
		out_file = self.utils.get_temp_file("file_name", ".fasta")
		self.assertTrue(os.path.exists(out_file))
		vest_data = [[os.path.join(path_to_fasta, "fasta_1.fasta"), "test_sample_1", project_sample.id],
					[os.path.join(path_to_fasta, "fasta_2.fasta"), "test_sample_2", project_sample_1.id]]
		utils.merge_fasta_files(vest_data , out_file)
		self.assertTrue(os.path.exists(out_file))
		vect_read_file = utils.read_text_file(out_file)
		self.assertEqual(78, len(vect_read_file))
		
		self.assertTrue(">test_sample_1__EVA003_S91" in vect_read_file)
		self.assertTrue(">test_sample_2__EVA003_S91" in vect_read_file)
		
		try:
			project_sample = ProjectSample.objects.get(id=project_sample.id)
			self.assertEqual("test_sample_1{}".format(Constants.SEPARATOR_sample_record_id),
				project_sample.seq_name_all_consensus)
		except ProjectSample.DoesNotExist:	## need to create with last version
			self.fail("Must exist...")
			
		try:
			project_sample = ProjectSample.objects.get(id=project_sample_1.id)
			self.assertEqual("test_sample_2{}".format(Constants.SEPARATOR_sample_record_id),
				project_sample.seq_name_all_consensus)
		except ProjectSample.DoesNotExist:	## need to create with last version
			self.fail("Must exist...")
		os.unlink(out_file)


	def test_merge_fasta_files_and_join_multifasta(self):
		"""
		testing merging fasta
		"""
		utils = Utils()
		path_to_fasta = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTA)
		out_file = self.utils.get_temp_file("file_name", ".fasta")
		self.assertTrue(os.path.exists(out_file))
		vest_data = [[os.path.join(path_to_fasta, "fasta_1.fasta"), "test_sample_1", -1],
					 [os.path.join(path_to_fasta, "fasta_2.fasta"), "test_sample_1", -1],
					 [os.path.join(path_to_fasta, "fasta_3.fasta"), "test_sample_2", -1]]
		utils.merge_fasta_files_and_join_multifasta(vest_data , out_file)
		self.assertTrue(os.path.exists(out_file))
		
		with open(out_file, "rU") as handle_fasta:
			count = 0
			for record in SeqIO.parse(handle_fasta, "fasta"):
				if (count == 0):
					self.assertEqual("test_sample_1", record.id)
					self.assertEqual(1982, len(str(record.seq)))
				elif (count == 1):
					self.assertEqual("test_sample_1_1", record.id)
					self.assertEqual(1982, len(str(record.seq)))
				elif (count == 2):
					self.assertEqual("test_sample_2", record.id)
					self.assertEqual(2042, len(str(record.seq)))
				count += 1
				self.assertTrue(count < 4)
		os.unlink(out_file)


	def test_merge_fasta_and_join_sequences(self):
		"""
		testing merging fasta and join sequences
		"""
		
		utils = Utils()
		path_to_fasta = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_MERGE_FASTA)
		vect_read_file = utils.read_text_file(os.path.join(path_to_fasta, "fasta_1.fasta"))
		self.assertEqual(39, len(vect_read_file))
		
		out_file = self.utils.get_temp_file("file_name", ".fasta")
		self.assertTrue(os.path.exists(out_file))
		vect_elements = ['EVA003_S91', 'EVA001_S66', 'EVA002_S52', 'EVA002_S53', 'EVA011_S54']
		utils.merge_fasta_and_join_sequences(path_to_fasta, vect_elements, out_file)
		self.assertTrue(os.path.exists(out_file))
		vect_read_file = utils.read_text_file(out_file)
		self.assertEqual(70, len(vect_read_file))
		
		self.assertTrue(">fasta_1" in vect_read_file)
		self.assertTrue(">fasta_2" in vect_read_file)
		os.unlink(out_file)

	def test_get_last_name_from_fasta(self):
		
		utils = Utils()
		path_to_fasta = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_MERGE_FASTA)
		file_fasta = os.path.join(path_to_fasta, "fasta_1.fasta")
		sz_name = utils.get_last_name_from_fasta(file_fasta)
		self.assertEqual("EVA011_S54", sz_name)
		
	def test_get_type_files(self):
		
		utils = Utils()
		path_to_fasta = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_MERGE_FASTA)
		self.assertTrue(os.path.exists(os.path.join(path_to_fasta, "fasta_1.fasta")))
		self.assertEqual(Constants.FORMAT_FASTA, utils.get_type_file(os.path.join(path_to_fasta, "fasta_1.fasta")))
	
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		self.assertTrue(os.path.exists(path_file))
		self.assertEqual(Constants.FORMAT_FASTQ_illumina, utils.get_type_file(path_file))

		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_nanopore)
		self.assertTrue(os.path.exists(path_file))
		self.assertEqual(Constants.FORMAT_FASTQ_ont, utils.get_type_file(path_file))

		gff_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR_GFF, "insa_flu_temp_empty.gff")
		self.assertTrue(os.path.exists(gff_file))
		try:
			self.assertTrue(utils.get_type_file(gff_file))
		except Exception as e:
			self.assertEqual("File is not in fastq.gz format.", e.args[0])
			
		gff_file_fake_fastq = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR_GFF, "insa_flu_temp_2.gff.fastq.gz")
		self.assertTrue(os.path.exists(gff_file_fake_fastq))
		try:
			self.assertTrue(utils.get_type_file(gff_file_fake_fastq))
		except Exception as e:
			self.assertEqual("File is not in fastq.gz format.", e.args[0])
			
	def test_parse_amino_HGVS_code(self):
		""" test HGVS parse code """
		
		utils = Utils()
		self.assertEqual("p.N292N", utils.parse_amino_HGVS_code("p.Asn292Asn"))
		self.assertEqual("p.S299S", utils.parse_amino_HGVS_code("p.Ser299Ser"))
		self.assertEqual("p.S299fs", utils.parse_amino_HGVS_code("p.Ser299fs"))
		self.assertEqual("p.S299*", utils.parse_amino_HGVS_code("p.Ser299*"))
		self.assertEqual("p.X299*", utils.parse_amino_HGVS_code("p.Ter299*"))
		self.assertEqual("p.A299D", utils.parse_amino_HGVS_code("p.Ala299Asp"))
		self.assertEqual("p.G299Dir", utils.parse_amino_HGVS_code("p.Gly299Dir"))
		self.assertEqual("p.SV534SV", utils.parse_amino_HGVS_code("p.SerVal534SerVal"))
		self.assertEqual("p.S534SV", utils.parse_amino_HGVS_code("p.Ser534SerVal"))
		self.assertEqual("", utils.parse_amino_HGVS_code("xpto"))

	def test_add_freq_ao_ad_and_type_to_vcf(self):

		utils = Utils()
		coverage_limit = 20
		freq_vcf_limit = 0.51
		
		vcf_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "run_snippyis_single_covid.vcf")
		depth_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "run_snippyis_single_covid.depth.gz")
		expecteded_vcf_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "add_tags_from_medaka_expected_2.vcf")
		expecteded_vcf_file_removed = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "add_tags_from_medaka_expected_2_removed.vcf")
		self.assertTrue(os.path.exists(vcf_file))
		self.assertTrue(os.path.exists(depth_file))
		self.assertTrue(os.path.exists(expecteded_vcf_file))
		
		vcf_file_out = utils.get_temp_file("vcf_medaka", ".vcf")
		vcf_file_removed_out = utils.get_temp_file("vcf_medaka", "_2.vcf")
		utils.add_freq_ao_ad_and_type_to_vcf(vcf_file, depth_file, vcf_file_out, vcf_file_removed_out, coverage_limit,
								freq_vcf_limit)
		
		self.assertTrue(os.path.exists(vcf_file_out))
		self.assertTrue(os.path.exists(expecteded_vcf_file))
		self.assertTrue(filecmp.cmp(vcf_file_out, expecteded_vcf_file))
		self.assertTrue(filecmp.cmp(vcf_file_removed_out, expecteded_vcf_file_removed))
		if (os.path.exists(vcf_file_out)): os.unlink(vcf_file_out)
		if (os.path.exists(vcf_file_removed_out)): os.unlink(vcf_file_removed_out)

	def test_add_freq_ao_ad_and_type_to_vcf_2(self):

		utils = Utils()
		coverage_limit = 20
		freq_vcf_limit = 0.38
		
		vcf_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "run_snippyis_single_covid.vcf")
		depth_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "run_snippyis_single_covid.depth.gz")
		expecteded_vcf_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "add_tags_from_medaka_expected.vcf")
		expecteded_vcf_file_removed = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "add_tags_from_medaka_expected_3_removed.vcf")
		self.assertTrue(os.path.exists(vcf_file))
		self.assertTrue(os.path.exists(depth_file))
		self.assertTrue(os.path.exists(expecteded_vcf_file))
		
		vcf_file_out = utils.get_temp_file("vcf_medaka", ".vcf")
		vcf_file_removed_out = utils.get_temp_file("vcf_medaka", "_2.vcf")
		utils.add_freq_ao_ad_and_type_to_vcf(vcf_file, depth_file, vcf_file_out, vcf_file_removed_out, coverage_limit,
								freq_vcf_limit)
		
		self.assertTrue(os.path.exists(vcf_file_out))
		self.assertTrue(os.path.exists(expecteded_vcf_file))
		self.assertTrue(filecmp.cmp(vcf_file_out, expecteded_vcf_file))
		self.assertTrue(filecmp.cmp(vcf_file_removed_out, expecteded_vcf_file_removed))
		if (os.path.exists(vcf_file_out)): os.unlink(vcf_file_out)
		if (os.path.exists(vcf_file_removed_out)): os.unlink(vcf_file_removed_out)

	def test_add_freq_ao_ad_and_type_to_vcf_3(self):

		utils = Utils()
		coverage_limit = 20
		freq_vcf_limit = 0.39
		
		vcf_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "run_snippyis_single_covid.vcf")
		depth_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "run_snippyis_single_covid.depth.gz")
		expecteded_vcf_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "add_tags_from_medaka_expected_2.vcf")
		expecteded_vcf_file_removed = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "add_tags_from_medaka_expected_4_removed.vcf")
		self.assertTrue(os.path.exists(vcf_file))
		self.assertTrue(os.path.exists(depth_file))
		self.assertTrue(os.path.exists(expecteded_vcf_file))
		
		vcf_file_out = utils.get_temp_file("vcf_medaka", ".vcf")
		vcf_file_removed_out = utils.get_temp_file("vcf_medaka", "_2.vcf")
		utils.add_freq_ao_ad_and_type_to_vcf(vcf_file, depth_file, vcf_file_out, vcf_file_removed_out, coverage_limit,
								freq_vcf_limit)
		
		self.assertTrue(os.path.exists(vcf_file_out))
		self.assertTrue(os.path.exists(expecteded_vcf_file))
		self.assertTrue(filecmp.cmp(vcf_file_out, expecteded_vcf_file))
		self.assertTrue(filecmp.cmp(vcf_file_removed_out, expecteded_vcf_file_removed))
		if (os.path.exists(vcf_file_out)): os.unlink(vcf_file_out)
		if (os.path.exists(vcf_file_removed_out)): os.unlink(vcf_file_removed_out)

	def test_add_freq_ao_ad_and_type_to_vcf_4(self):

		utils = Utils()
		coverage_limit = 20
		freq_vcf_limit = 0.39
		
		vcf_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "UK_ERR4082154.vcf")
		depth_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "UK_ERR4082154.depth.gz")
		expecteded_vcf_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "add_tags_from_medaka_expected_3.vcf")
		self.assertTrue(os.path.exists(vcf_file))
		self.assertTrue(os.path.exists(depth_file))
		self.assertTrue(os.path.exists(expecteded_vcf_file))
		
		vcf_file_out = utils.get_temp_file("vcf_medaka", ".vcf")
		vcf_file_removed_out = utils.get_temp_file("vcf_medaka", "_2.vcf")
		utils.add_freq_ao_ad_and_type_to_vcf(vcf_file, depth_file, vcf_file_out, vcf_file_removed_out, coverage_limit,
								freq_vcf_limit)
		
		self.assertTrue(os.path.exists(vcf_file_out))
		self.assertTrue(os.path.exists(expecteded_vcf_file))
		self.assertTrue(filecmp.cmp(vcf_file_out, expecteded_vcf_file))
		if (os.path.exists(vcf_file_out)): os.unlink(vcf_file_out)
		if (os.path.exists(vcf_file_removed_out)): os.unlink(vcf_file_removed_out)
		
	def test_add_freq_ao_ad_and_type_to_vcf_5(self):

		utils = Utils()
		coverage_limit = 20
		freq_vcf_limit = 0.39
		
		vcf_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "UK_ERR4248992.vcf")
		depth_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "UK_ERR4248992.depth.gz")
		expecteded_vcf_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "add_tags_from_medaka_expected_4.vcf")
		self.assertTrue(os.path.exists(vcf_file))
		self.assertTrue(os.path.exists(depth_file))
		self.assertTrue(os.path.exists(expecteded_vcf_file))
		
		vcf_file_out = utils.get_temp_file("vcf_medaka", ".vcf")
		vcf_file_removed_out = utils.get_temp_file("vcf_medaka", "_2.vcf")
		utils.add_freq_ao_ad_and_type_to_vcf(vcf_file, depth_file, vcf_file_out, vcf_file_removed_out, coverage_limit,
								freq_vcf_limit)
		
		self.assertTrue(os.path.exists(vcf_file_out))
		self.assertTrue(os.path.exists(expecteded_vcf_file))
		self.assertTrue(filecmp.cmp(vcf_file_out, expecteded_vcf_file))
		if (os.path.exists(vcf_file_out)): os.unlink(vcf_file_out)
		if (os.path.exists(vcf_file_removed_out)): os.unlink(vcf_file_removed_out)
		
	def test_medaka_models(self):
		
		utils = Utils()
		vect_data = utils.get_all_medaka_models()
		self.assertTrue(len(vect_data) > 0)
		self.assertTrue(SoftwareNames.SOFTWARE_Medaka_default_model in vect_data)

	
	def test_difference_files(self):

		utils = Utils()
		
		fasta_file_out = utils.get_temp_file("fasta_data", ".fasta")
		with open(fasta_file_out, 'w') as handle_write:
			handle_write.write(">1\nAAAAAAAAAA\n")
			handle_write.write(">2\nAAAAAAAAAA\n")
		self.assertFalse(utils.is_differente_fasta_size(fasta_file_out, 1))
		self.assertFalse(utils.is_differente_fasta_size(fasta_file_out, 10))
		
		with open(fasta_file_out, 'w') as handle_write:
			handle_write.write(">1\nAAAAAAAAAA\n")
			handle_write.write(">2\nAAAAAAAAAAAAA\n")
		self.assertTrue(utils.is_differente_fasta_size(fasta_file_out, 1))
		self.assertTrue(utils.is_differente_fasta_size(fasta_file_out, 5))
		self.assertTrue(utils.is_differente_fasta_size(fasta_file_out, 10))
		self.assertFalse(utils.is_differente_fasta_size(fasta_file_out, 30))
		
		with open(fasta_file_out, 'w') as handle_write:
			handle_write.write(">2\nAAAAAAAAAAAAAAAAAAAAAA\n")
			handle_write.write(">1\nAAAAAAAAAA\n")
		self.assertTrue(utils.is_differente_fasta_size(fasta_file_out, 1))
		self.assertTrue(utils.is_differente_fasta_size(fasta_file_out, 10))
		self.assertTrue(utils.is_differente_fasta_size(fasta_file_out, 50))
		self.assertFalse(utils.is_differente_fasta_size(fasta_file_out, 60))

		### remove file
		utils.remove_file(fasta_file_out)

	def test_get_number_sequences_fastq(self):
		
		utils = Utils()
		path_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		self.assertTrue(os.path.exists(path_file))
		self.assertEqual((44425, 143.6, 20.9), utils.get_number_sequences_fastq(path_file))
		
		fasta_file_out = utils.get_temp_file("fasta_data", ".fasta")
		with open(fasta_file_out, 'w') as handle_write:
			handle_write.write(">2\nAAAAAAAAAAAAAAAAAAAAAA\n")
			handle_write.write(">1\nAAAAAAAAAA\n")
		try:
			self.assertEqual((44425, 10.5, 20.9), utils.get_number_sequences_fastq(path_file))
			self.fail("must throw exception")
		except Exception as e:
			pass

	def test_mask_sequences(self):
		""" mask sequences """
		
		utils = Utils()
		sequence = SeqRecord(Seq("AACCTTTAATTTAATTTATTTAATTT"), id="xpto")
		mask_sites = "5,6"
		mask_from_beginning = "3"
		mask_from_end = "2"
		mask_range = "10-12,14-16"
		seq_changed = utils.mask_sequence(sequence, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual("NNNCNNTAANNNANNNTATTTAATNN", str(seq_changed.seq))
		
		sequence = SeqRecord(Seq("AACCTTTAATTTAATTTATTTAATTT"), id="xpto")
		mask_sites = ""
		mask_from_beginning = ""
		mask_from_end = ""
		mask_range = ""
		seq_changed = utils.mask_sequence(sequence, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual("AACCTTTAATTTAATTTATTTAATTT", str(seq_changed.seq))
		
		sequence = SeqRecord(Seq("AACCTTTAATTTAATTTATTTAATTT"), id="xpto")
		mask_sites = None
		mask_from_beginning = None
		mask_from_end = None
		mask_range = None
		seq_changed = utils.mask_sequence(sequence, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual("AACCTTTAATTTAATTTATTTAATTT", str(seq_changed.seq))




