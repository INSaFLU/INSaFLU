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
	MANAGING_DIR_GFF = "gff_files"
	MANAGING_FILES_FASTA = "A_H3N2_A_Hong_Kong_4801_2014.fasta"
	MANAGING_FILES_SNPEF_config = "snpeff.config"
	MANAGING_FILES_FAIL_FASTA = "A_H3N2_A_Hong_Kong_4801_2014_fail.fasta"
	MANAGING_FILES_GBK = "A_H3N2_A_Hong_Kong_4801_2014.gbk"
	MANAGING_FILES_BED = "A_H3N2_A_Hong_Kong_4801_2014.bed"
	MANAGING_FILES_GFF = "A_H3N2_A_Hong_Kong_4801_2014.gff"
	MANAGING_FILES_FREEBAYES_VCF = "freebayes.vcf"
	MANAGING_FILES_FREEBAYES_ANNOTATED_VCF = "freebayes_annotated.vcf"
	MANAGING_FILES_FREEBAYES_ANNOTATED_VCF_2 = "freebayes_annotated_2.vcf"
	MANAGING_FILES_GBK_MISS_ONE = "A_H3N2_A_Hong_Kong_4801_2014_miss_one.gbk"
	MANAGING_FILES_GBK_DIFF_LENGTH = "A_H3N2_A_Hong_Kong_4801_2014_diff_length.gbk"
	MANAGING_FILES_FASTA_FAKE_GZ = "test.fasta.gz"
	MANAGING_FILES_pb2_reversed = "pb2_reversed.gbk"
	MANAGING_TEST_INFLUENZA_FILE = "test_influenza_typing.fasta"
	MANAGING_FILES_FASTA_2 = "temp_2.fasta"
	MANAGING_FILES_GBK_2 = "temp_2.gbk"
	MANAGING_TEST_ABRICATE = "abricate_out.txt"
	MANAGING_TEST_ABRICATE_2 = "abricate_2.txt"
	MANAGING_TWO_GENES_JOINED_GBK = "TwoGenesJoined.gbk"
	MANAGING_TWO_GENES_JOINED_BED = "TwoGenesJoined.bed"
	MANAGING_CONSENSUS_ALIGNMENT_PROTEIN = "concensus_alignment_protein.fasta"
	MANAGING_OUT_PROTEIN = "out_protein.faa"
	MANAGING_MAFFT_IN = 'mafft_2.fasta'
	MANAGING_MAFFT_IN_EXPECTED = 'mafft_2_expected.fasta'
	MANAGING_MAFFT_IN_PROTEIN = 'mafft_protein.fasta'
	MANAGING_MAFFT_IN_PROTEIN_EXPECTED = 'mafft_protein_expected.fasta'
	MANAGING_TREE_out_protein = 'Alignment_MP_M.nex'
	MANAGING_TEMPLATE_INPUT = "template_input.tsv"
	MANAGING_TEMPLATE_INPUT_error = "vncServer.txt"
	MANAGING_TEMPLATE_INPUT_FAIL_HEADER = "template_input_fail_header.tsv"
	MANAGING_TEMPLATE_INPUT_DATA_TSV = "template_input_data.tsv"
	MANAGING_TEMPLATE_INPUT_DATA_UPDATE_TSV = "template_input_data_update.tsv"
	MANAGING_TEMPLATE_INPUT_DATA_CSV = "template_input_data.csv"
	MANAGING_TEMPLATE_INPUT_DATA_3_CSV = "Cyril_input_test1_minus_heavy_files.txt"
	MANAGING_TEMPLATE_INPUT_DATA_2_CSV = "template_input_data_2.csv"
	MANAGING_TEMPLATE_INPUT_DATA_ASSIGN2CONTIGS_TSV = "test_update_assign2contigs.tsv"
	MANAGING_TEMPLATE_INPUT_DATA_FAIL = "template_input_data_fail.csv"
	MANAGING_TEMPLATE_INPUT_big_data_csv = "sample_data_season_2016_17_.csv"
	MANAGING_FILES_RIM_sample = "Rim1_haplogroup_A1a1a1_Italy.fasta"
	MANAGING_FILES_RIM_sample_gbk = "Rim1_haplogroup_A1a1a1_Italy.gbk"
	MANAGING_FILES_RIM_sample_gff = "Rim1_haplogroup_A1a1a1_Italy.gff"
	MANAGING_FILES_PV1_KX162693_fasta = "PV1_KX162693.fasta"
	MANAGING_FILES_PV1_KX162693_gbk = "PV1_KX162693.gb"
	MANAGING_FILES_PV1_KX162693_clean_version_gbk = "PV1_KX162694_clean_version.gb"
	MANAGING_FILES_PV1_KX162693_spades_out_fasta = "PV1_KX162693_spades_out.fasta"
	MANAGING_FILES_TEST_LENGTH_fasta_1 = "test_length_fasta_1.fasta"
	MANAGING_FILES_TEST_LENGTH_fasta_2 = "test_length_fasta_2.fasta"
	MANAGING_FILES_PV1_KX162693_bed = "PV1_KX162693.bed"
	
	MANAGING_FILES_PV_Sabin1_V01150_fasta = "PV_Sabin1_V01150.fasta"
	MANAGING_FILES_PV_Sabin1_V01150_gbk = "PV_Sabin1_V01150.gb"
	MANAGING_FILES_PV_Sabin1_V01150_clean_version_gbk = "PV_Sabin1_V01150_clean_version.gb"
	MANAGING_FILES_PV_Sabin1_V01150_spades_out_fasta = "PV1_KX162693_spades_out.fasta"
	MANAGING_FILES_PV_Sabin1_V01150_bed = "PV1_KX162693.bed"
	
	MANAGING_FILES_TEMPLATE_MULTIPLE_FILES_data_csv = "test_multiple_files.csv"

	
	### consensus files
	FILES_EVA001_S66_consensus = "EVA001_S66.consensus.fasta"
	FILES_EVA002_S52_consensus = "EVA002_S52.consensus.fasta"
	FILES_EVA003_S91_consensus = "EVA003_S91.consensus.fasta"
	FILES_EVA011_S54_consensus = "EVA011_S54.consensus.fasta"
	FILE_OUT_MAUVE = 'mauve.xfma'
	FILE_OUT_CONVERT_MAUVE_RESULT = 'convert_mauve_result.fasta'
	FILE_OUT_MAFFT_RESULT = 'mafft_result.fasta'
	FILE_OUT_NEX_RESULT = 'mafft_result.nex'
	FILE_FASTTREE_RESULT_NWK = 'mafft_result.nwk'
	FILE_FASTTREE_RESULT_TREE = 'mafft_result.tree'
	FILE_TO_CLEAN_NAME = "to_clean_name.fasta"
	FILE_RESULT_TO_CLEAN_NAME = "result_to_clean_name.fasta"
	MANAGING_FILE_GRAPH_VAR_HTML = 'file_graph_var.html'
	MANAGING_FILE_GRAPH_VAR_PNG = 'file_graph_var.png'
	
	
	DIR_ABRICATE = "abricate"
	DIR_COVERAGE = "coverage"
	DIR_BAM = "bam"
	DIR_TREE = "tree"
	DIR_VCF = "vcf"
	DIR_IMAGES = "images"
	DIR_PROJECTS = "projects"
	DIR_GLOBAL_PROJECT = "global_project"
	DIR_INPUT_FILES = "input_files"
	DIR_FASTA = "fasta"
	
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
	FASTQ5_1 = "x_mix_H3N2_H1N1_1P.fastq.gz"
	FASTQ5_2 = "x_mix_H3N2_H1N1_2P.fastq.gz"
	FASTQ6_1 = "RL017_1P.fastq.gz"
	FASTQ6_2 = "RL017_2P.fastq.gz"
	FASTQ7_1 = "Ae-MIT-602_S25_L001_R1_001.fastq.gz"
	FASTQ7_2 = "Ae-MIT-602_S25_L001_R2_001.fastq.gz"
	FASTQ8_1 = "SU104_S10_L001_R1_001.fastq.gz"
	FASTQ8_2 = "SU104_S10_L001_R2_001.fastq.gz"
	FASTQ9_1 = "CLEVB-1_S2_L001_R1_001.fastq.gz"
	FASTQ9_2 = "CLEVB-1_S2_L001_R2_001.fastq.gz"
	FASTQ10_1 = "EVA082_S12_L001_R1_001.fastq.gz"
	FASTQ10_2 = "EVA082_S12_L001_R2_001.fastq.gz"
	FASTQ10_1_test = "EVA082_S12_L001_R2_001_test.fastq.gz"
	FASTQ1_nanopore = "nanopore.fastq.gz"
	
	FASTQ11_1 = "EVA174_S64_L001_R1_001.fastq.gz"
	FASTQ11_2 = "EVA174_S64_L001_R2_001.fastq.gz"
	
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




