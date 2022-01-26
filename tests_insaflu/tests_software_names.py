'''
Created on May 17, 2018

@author: mmp
'''
import os
import unittest
from constants.software_names import SoftwareNames

"""
Used to check all 
"""
class Test(unittest.TestCase):


	def setUp(self):
		pass


	def tearDown(self):
		pass


	def test_software_names(self):
		
		software_names = SoftwareNames()
		
		self.assertTrue(os.path.exists(software_names.get_samtools()), "Fail to test - " + software_names.get_samtools_name())
		self.assertTrue(os.path.exists(software_names.get_spades()), "Fail to test - " + software_names.get_spades_name())
		self.assertTrue(os.path.exists(software_names.get_abricate()), "Fail to test - " + software_names.get_abricate_name())
		self.assertTrue(os.path.exists(software_names.get_fastq()), "Fail to test - " + software_names.get_fastq_name())
		self.assertTrue(os.path.exists(software_names.get_trimmomatic()), "Fail to test - " + software_names.get_trimmomatic_name())
		self.assertTrue(os.path.exists(software_names.get_snippy()), "Fail to test - " + software_names.get_snippy_name())
		self.assertTrue(os.path.exists(software_names.get_snippy_vcf_to_tab()), "Fail to test - " + software_names.get_snippy_vcf_to_tab_name())
		self.assertTrue(os.path.exists(software_names.get_snp_eff()), "Fail to test - " + software_names.get_snp_eff_name())
		self.assertTrue(os.path.exists(software_names.get_freebayes()), "Fail to test - " + software_names.get_freebayes_name())
		self.assertTrue(os.path.exists(software_names.get_freebayes_parallel()), "Fail to test - " + software_names.get_freebayes_name() + " parallel")
		self.assertTrue(os.path.exists(software_names.get_fasta_generate_regions()), "Fail to test - " + software_names.get_fasta_generate_regions_name())
		self.assertTrue(os.path.exists(software_names.get_coverage_to_regions()), "Fail to test - " + software_names.get_coverage_to_regions_name(), )
		self.assertTrue(os.path.exists(software_names.get_bamtools()), "Fail to test - " + software_names.get_bamtools_name())
		self.assertTrue(os.path.exists(software_names.get_bgzip()), "Fail to test - " + software_names.get_bgzip_name())
		self.assertTrue(os.path.exists(software_names.get_tabix()), "Fail to test - " + software_names.get_tabix_name())
		self.assertTrue(os.path.exists(software_names.get_prokka()), "Fail to test - " + software_names.get_prokka_name())
		self.assertTrue(os.path.exists(software_names.get_mauve()), "Fail to test - " + software_names.get_mauve_name())
		self.assertTrue(os.path.exists(software_names.get_convert_mauve()), "Fail to test - " + software_names.get_convert_mauve_name())
		self.assertTrue(os.path.exists(software_names.get_mafft()), "Fail to test - " + software_names.get_mafft_name())
		self.assertTrue(os.path.exists(software_names.get_seqret()), "Fail to test - " + software_names.get_seqret_name())
		self.assertTrue(os.path.exists(software_names.get_fasttree()), "Fail to test - " + software_names.get_fasttree_name(), )
		self.assertTrue(os.path.exists(software_names.get_fastqtools_sample()), "Fail to test - " + software_names.get_fastqtools_sample_name())

		### pangolin names
		dt_software_names_pangolin = software_names.get_pangolin_all_names_version()
		self.assertTrue(SoftwareNames.SOFTWARE_Pangolin_name in dt_software_names_pangolin)
		self.assertTrue(SoftwareNames.SOFTWARE_Pangolin_learn_name in dt_software_names_pangolin)
		self.assertTrue(SoftwareNames.SOFTWARE_Pangolin_designation_name in dt_software_names_pangolin)
		self.assertTrue(SoftwareNames.SOFTWARE_Pangolin_constellations_name in dt_software_names_pangolin)
		self.assertTrue(SoftwareNames.SOFTWARE_Pangolin_scorpio_name in dt_software_names_pangolin)


if __name__ == "__main__":
	#import sys;sys.argv = ['', 'Test.testName']
	unittest.main()