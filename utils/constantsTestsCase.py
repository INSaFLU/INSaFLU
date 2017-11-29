'''
Created on Oct 13, 2017

@author: mmp
'''

import os

class ConstantsTestsCase(object):
	'''
	classdocs
	'''
	MANAGING_TESTS = "tests"
	MANAGING_DIR = "managing_files"
	MANAGING_FILES_FASTA = "A_H3N2_A_Hong_Kong_4801_2014.fasta"
	MANAGING_FILES_SNPEF_config = "snpeff.config"
	MANAGING_FILES_FAIL_FASTA = "A_H3N2_A_Hong_Kong_4801_2014_fail.fasta"
	MANAGING_FILES_GBK = "A_H3N2_A_Hong_Kong_4801_2014.gbk"
	MANAGING_FILES_GFF = "A_H3N2_A_Hong_Kong_4801_2014.gff"
	MANAGING_FILES_FREEBAYES_VCF = "freebayes.vcf"
	MANAGING_FILES_FREEBAYES_ANNOTATED_VCF = "freebayes_annotated.vcf"
	MANAGING_FILES_GBK_MISS_ONE = "A_H3N2_A_Hong_Kong_4801_2014_miss_one.gbk"
	MANAGING_FILES_GBK_DIFF_LENGTH = "A_H3N2_A_Hong_Kong_4801_2014_diff_length.gbk"
	MANAGING_FILES_FASTA_FAKE_GZ = "test.fasta.gz"
	MANAGING_TEST_INFLUENZA_FILE = "test_influenza_typing.fasta"
	MANAGING_FILES_FASTA_2 = "temp_2.fasta"
	MANAGING_FILES_GBK_2 = "temp_2.gbk"
	MANAGING_TEST_ABRICATE = "abricate_out.txt"

	### consensus files
	FILES_EVA001_S66_consensus = "EVA001_S66.consensus.fasta"
	FILES_EVA002_S52_consensus = "EVA002_S52.consensus.fasta"
	FILES_EVA003_S91_consensus = "EVA003_S91.consensus.fasta"
	FILES_EVA011_S54_consensus = "EVA011_S54.consensus.fasta"
	FILE_OUT_MAUVE = 'mauve.xfma'
	FILE_OUT_CONVERT_MAUVE_RESULT = 'convert_mauve_result.fasta'
	FILE_OUT_MAFFT_RESULT = 'mafft_result.fasta'
	FILE_FASTTREE_RESULT_NWK = 'mafft_result.nwk'
	FILE_FASTTREE_RESULT_TREE = 'mafft_result.tree'
	FILE_TO_CLEAN_NAME = "to_clean_name.fasta"
	FILE_RESULT_TO_CLEAN_NAME = "result_to_clean_name.fasta"
	MANAGING_FILE_GRAPH_VAR_HTML = 'file_graph_var.html'
	
	DIR_ABRICATE = "abricate"
	DIR_COVERAGE = "coverage"
	DIR_BAM = "bam"
	DIR_VCF = "vcf"
	DIR_IMAGES = "images"
	DIR_PROJECTS = "projects"
	DIR_GLOBAL_PROJECT = "global_project"
	
	## fastq files
	DIR_FASTQ = "fastq"
	FASTQ1_1 = "EVA001_S66_L001_R1_001.fastq.gz"
	FASTQ1_2 = "EVA001_S66_L001_R2_001.fastq.gz"
	FASTQ2_1 = "EVA002_S52_L001_R1_001.fastq.gz"
	FASTQ2_2 = "EVA002_S52_L001_R2_001.fastq.gz"
	FASTQ3_1 = "EVA003_S91_L001_R1_001.fastq.gz"
	FASTQ3_2 = "EVA003_S91_L001_R2_001.fastq.gz"
	FASTQ4_1 = "EVA011_S54_L001_R1_001.fastq.gz"
	FASTQ4_2 = "EVA011_S54_L001_R2_001.fastq.gz"

	META_KEY_TEST = "meta_key_test"
	VALUE_TEST = "value_test"
	VALUE_TEST_2 = "value_test_2"
	TEST_USER_NAME = "test_user"
	
	
	def get_all_fastq_files(self, base_dir):
		"""
		get all fastq files
		"""
		return [[os.path.join(base_dir, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1),\
						os.path.join(base_dir, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)],\
					[os.path.join(base_dir, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ2_1),\
						os.path.join(base_dir, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ2_2)],\
					[os.path.join(base_dir, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ3_1),\
						os.path.join(base_dir, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ3_2)],\
					[os.path.join(base_dir, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ4_1),\
						os.path.join(base_dir, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ4_2)]]
		
	def get_all_consensus_files(self, base_dir):
		"""
		get all consensus files
		"""
		return [os.path.join(base_dir, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.FILES_EVA001_S66_consensus),\
				os.path.join(base_dir, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.FILES_EVA002_S52_consensus),\
				os.path.join(base_dir, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.FILES_EVA003_S91_consensus),\
				os.path.join(base_dir, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.FILES_EVA011_S54_consensus)]




