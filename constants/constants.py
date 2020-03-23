'''
Created on Oct 13, 2017

@author: mmp
'''
from enum import Enum

class Constants(object):
	'''
	classdocs
	'''
	
	### default user that has the default references to be used in mapping
	DEFAULT_USER = "system"
	DEFAULT_USER_PASS = "default_user_123_$%_2"
	## DEFAULT_USER_EMAIL = "insaflu@insa.min-saude.pt"		### it's defined in .env
	
	### user anonymous
	USER_ANONYMOUS = "demo"
	USER_ANONYMOUS_PASS = "demo_user"
	## USER_ANONYMOUS_EMAIL = "insaflu@insa.min-saude.pt"		### it's defined in .env
	
	META_KEY_VALUE_NOT_NEED = "value not needed"
	
	## MAX LOCUS FROM FASTA
	MAX_SEQUENCES_FROM_FASTA = 20		### update this values
	MAX_LENGTH_SEQ_NAME = 20			###  it must be less than 20 because of prokka constrainments 
	
	## MAX LENGTH_SEQUENCE_FROM_FASTA
	MAX_LENGTH_SEQUENCE_FROM_FASTA = 10000 
	MAX_LENGTH_SEQUENCE_TOTAL_FROM_FASTA = 50000
# 	MAX_FASTQ_FILE = 50971520			## 50M
# 	MAX_FASTQ_FILE_WITH_DOWNSIZE = 250971520			## 250M
	MAX_REF_FASTA_FILE = 100000			## 100k
	MAX_REF_GENBANK_FILE = 150000		## 150k
	
	### has the minimun number of files to calculate global files
	MINIMUN_NUMER_SAMPLES_CACULATE_GLOBAL_FILES = 2
	
	### Session variables
	NUMBER_LOCUS_FASTA_FILE = "number_locus_fasta_file"
	
	## https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
	## translate table number
	TRANSLATE_TABLE_NUMBER = 11
	
	## start expand tag name rows
	START_EXPAND_SAMPLE_TAG_NAMES_ROWS = 4
	
	DIR_PROCESSED_FILES_UPLOADS = "uploads"
	DIR_PROCESSED_PROCESSED = 'processed'
	
	### queue names
	QUEUE_SGE_NAMES = ['queue_1.q', 'queue_2.q']
	QUEUE_SGE_NAME_GLOBAL = 'all.q'
	QUEUE_SGE_NAME_FAST = 'fast.q'						## jobs that are fast to run
	QUEUE_SGE_NAME_INSA = 'insa.q'						## insa queue
##	QUEUE_SGE_NAME_EMAIL = 'email.q' direct for now...
	
	### separators
	SEPARATOR_COMMA = ','
	SEPARATOR_TAB = '\t'

	## DIR_PROCESSED_FILES_FROM_WEB/userId_<id>/refId_<id>
	DIR_PROCESSED_FILES_REFERENCE = DIR_PROCESSED_FILES_UPLOADS + "/references"
	DIR_PROCESSED_FILES_FASTQ = DIR_PROCESSED_FILES_UPLOADS + "/fastq"
	DIR_PROCESSED_FILES_PROJECT = DIR_PROCESSED_FILES_UPLOADS + "/project"
	DIR_PROCESSED_FILES_MULTIPLE_SAMPLES = DIR_PROCESSED_FILES_UPLOADS + "/multiple_samples"

	DIR_ICONS = "icons"
	DIR_TEMPLATE_INPUT = "template_input"
	TEMP_DIRECTORY = "/tmp"
	COUNT_DNA_TEMP_DIRECTORY = "insaFlu"

	FILE_TEMPLATE_INPUT_csv = "template_input.csv"
	FILE_TEMPLATE_INPUT_tsv = "template_input.tsv"
	FILE_TEMPLATE_INPUT_data_csv = "template_input_data.csv"
	FILE_TEMPLATE_INPUT_data_tsv = "template_input_data.tsv"

	FILE_TEMPLATE_INPUT_METADATA_csv = "template_metadata_input.csv"
	FILE_TEMPLATE_INPUT_METADATA_tsv = "template_metadata_input.tsv"
	FILE_TEMPLATE_INPUT_METADATA_data_csv = "template_input_metadata_data.csv"
	FILE_TEMPLATE_INPUT_METADATA_data_tsv = "template_input_metadata_data.tsv"
	
	FORMAT_FASTA = "fasta"
	FORMAT_FASTQ = "fastq"
	EXTENSION_ZIP = ".gz"
	MAX_LENGHT_ILLUMINA_FASQC_SEQ = 301
	
	## vect with keys to get the ID
	VECT_GENBANK_TAG_NAME = ['gene', 'CDS', 'locus_tag', 'protein_id']
	
	## Has all the versions of type influenza typing
	DIR_TYPE_CONTIGS_2_SEQUENCES = "db/contigs2sequences/"
	DIR_TYPE_IDENTIFICATION = "db/type_identification/"
	DIR_TYPE_REFERENCES = "db/references/"
	DIR_TEST_TYPE_REFERENCES = "tests/db/references/"
	
	INSAFLU_NAME = 'insaflu'
	
	### key for a session with project name
	PROJECT_NAME_SESSION = 'project_name_session'
	
	#####
	DIR_STATIC = "static"
	DIR_ICONS = "icons"
	ICON_GREEN_16_16 = "bullet_ball_glass_green.png"
	ICON_YELLOW_16_16 = "bullet_ball_glass_yellow.png"
	ICON_RED_16_16 = "bullet_ball_glass_red.png"
	
	AJAX_LOADING_GIF = "ajax-loading-gif.gif"
	AJAX_LOADING_GIF_13 = "ajax-loading-gif-13.gif"
	
	### data_set 
	DATA_SET_GENERIC = "Generic"	## default name for a dataset
	
	## NUMBER OF SETs to paginate
	PAGINATE_NUMBER = 12
	PAGINATE_NUMBER_SMALL = 2

	## tag for check box all in the tables
	CHECK_BOX_ALL = 'check_box_all'
	CHECK_BOX = 'check_box'
	GET_CHECK_BOX_SINGLE = 'get_check_box_single'
	GET_CHANGE_CHECK_BOX_SINGLE = 'get_change_check_box_single'
	COUNT_CHECK_BOX = 'count_check_boxes'
	CHECK_BOX_VALUE = 'value'
	CHECK_BOX_not_show_processed_files = 'check_box_not_show_processed_files'
	
	### empty value used in tables
	EMPTY_VALUE_TABLE = "-"
	
	### 
	EMPTY_VALUE_TYPE_SUBTYPE = "Not assigned"
	
	## session values
	SESSION_KEY_USER_ID = 'session_key_user_id'
	
	## errors
	PROJECT_NAME = 'project_name'
	ERROR_REFERENCE = 'error_reference'
	ERROR_PROJECT_NAME = 'error_project_name'

	vect_ambigous = ['R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N', '*']
	dt_ambigous = { 'R':'[AG]', 'Y':'[TC]', 'K':'[GT]', 'M':'[AC]', 'S':'[GC]', 
			'W':'[AT]', 'B':'[CGT]', 'D':'[AGT]', 'H':'[ACT]', 'V':'[ACG]',
			'N':'[ACGT]', '*':'.' }
	dict_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
					'R': 'Y', 
					'Y': 'R', 
					'K': 'M', 
					'M': 'K', 
					'S': 'S', 
					'W': 'W', 
					'B': 'V', 
					'D': 'H',
					'H': 'D',
					'V': 'B',
					'N': 'N'}
	
	def get_extensions_by_file_type(self, file_name, file_type):
		"""
		get extensions by file type
		"""
		if (file_type == FileType.FILE_BAM): return "{}.bam".format(file_name)
		if (file_type == FileType.FILE_BAM_BAI): return "{}.bam.bai".format(file_name)
		if (file_type == FileType.FILE_CONSENSUS_FA): return "{}.consensus.fa".format(file_name)
		if (file_type == FileType.FILE_CONSENSUS_FASTA): return "{}{}".format(file_name, FileExtensions.FILE_CONSENSUS_FASTA)
		if (file_type == FileType.FILE_CSV): return "{}.csv".format(file_name)
		if (file_type == FileType.FILE_DEPTH): return "{}.depth".format(file_name)
		if (file_type == FileType.FILE_DEPTH_GZ): return "{}.depth.gz".format(file_name)
		if (file_type == FileType.FILE_DEPTH_GZ_TBI): return "{}.depth.gz.tbi".format(file_name)
		if (file_type == FileType.FILE_TAB): return "{}.tab".format(file_name)
		if (file_type == FileType.FILE_VCF): return "{}.vcf".format(file_name)
		if (file_type == FileType.FILE_VCF_GZ): return "{}.vcf.gz".format(file_name)
		if (file_type == FileType.FILE_VCF_GZ_TBI): return "{}.vcf.gz.tbi".format(file_name)
		if (file_type == FileType.FILE_REF_FASTA): return "ref.fa"
		if (file_type == FileType.FILE_REF_FASTA_FAI): return "ref.fa.fai"
		return ""

	### complement
	def complement(self, seq):  
		complseq = [self.dict_complement[base] if base in self.dict_complement else base for base in seq]  
		return ''.join(complseq)
	
	#reverse
	def reverse_complement(self, seq):  
		seq = list(seq)  
		seq.reverse()   
		return self.complement(''.join(seq))
	
	
	###
	def ambiguos_to_unambiguous(self, sequence):
		for ambig in self.vect_ambigous:
			sequence = sequence.replace(ambig, self.dt_ambigous[ambig])
		#print sequence
		return sequence
	
	def get_diff_between_two_seq(self, seq1, seq2):
		if (len(seq1) != len(seq2)): return 0
		n_diff = 0
		for i in range(0, len(seq1)):
			if (seq1[i] != seq2[i]): n_diff += 1
		return n_diff

	def is_poly_n(self, sequence):
		"""
		if it has the same letter is poly N
		"""
		letter_previous = ""
		for letter in sequence:
			if (letter_previous == ""):
				letter_previous = letter
				continue
			if (letter_previous != letter): return False
		return True
	

	def short_name(self, name, max_size):
		"""
		short the name for the size of max_size
		"""
		if (len(name) > max_size): return "{}...{}".format(name[:int(max_size/2)], name[int(len(name) - (max_size/2)):])
		return name
	
class TypePath(Enum):
	"""
	Has the type of paths you can get from file paths
	"""
	MEDIA_ROOT = 0
	MEDIA_URL = 1


class FileType(Enum):
	"""
	Has the type of files
	[06:29:16] * /tmp/insafli/xpto/xpto.bam
	[06:29:16] * /tmp/insafli/xpto/xpto.bam.bai
	[06:29:16] * /tmp/insafli/xpto/xpto.bed
	[06:29:16] * /tmp/insafli/xpto/xpto.consensus.fa
	[06:29:16] * /tmp/insafli/xpto/xpto.consensus.subs.fa
	[06:29:16] * /tmp/insafli/xpto/xpto.csv
	[06:29:16] * /tmp/insafli/xpto/xpto.depth.gz
	[06:29:16] * /tmp/insafli/xpto/xpto.depth.gz.tbi
	[06:29:16] * /tmp/insafli/xpto/xpto.filt.subs.vcf
	[06:29:16] * /tmp/insafli/xpto/xpto.filt.subs.vcf.gz
	[06:29:16] * /tmp/insafli/xpto/xpto.filt.subs.vcf.gz.tbi
	[06:29:16] * /tmp/insafli/xpto/xpto.filt.vcf
	[06:29:16] * /tmp/insafli/xpto/xpto.gff
	[06:29:16] * /tmp/insafli/xpto/xpto.html
	[06:29:16] * /tmp/insafli/xpto/xpto.log
	[06:29:16] * /tmp/insafli/xpto/xpto.raw.vcf
	[06:29:16] * /tmp/insafli/xpto/xpto.tab
	[06:29:16] * /tmp/insafli/xpto/xpto.txt
	[06:29:16] * /tmp/insafli/xpto/xpto.vcf
	[06:29:16] * /tmp/insafli/xpto/xpto.vcf.gz
	[06:29:16] * /tmp/insafli/xpto/xpto.vcf.gz.tbi
				/tmp/insafli/xpto/ref/ref.fa
	"""
	FILE_BAM = 0
	FILE_BAM_BAI = 1
	FILE_CONSENSUS_FA = 2
	FILE_CONSENSUS_FASTA = 3
	FILE_DEPTH = 4
	FILE_DEPTH_GZ = 5
	FILE_DEPTH_GZ_TBI = 6
	FILE_TAB = 7
	FILE_VCF = 8
	FILE_VCF_GZ = 9
	FILE_VCF_GZ_TBI = 10
	FILE_CSV = 11
	FILE_REF_FASTA = 12		## ref/ref.fa
	FILE_REF_FASTA_FAI = 13		## ref/ref.fa.fai

class TypeFile(object):
	
	TYPE_FILE_fastq_gz = "fastq.gz"								## fastq.fz files
	TYPE_FILE_sample_file = "sample-file imported" 				## file that the user import with sample descriptions
	TYPE_FILE_sample_file_metadata = "sample-file metadata" 	## file that the user import with new metadata

	
class FileExtensions(object):
	"""
	file extensions
	"""
	FILE_TSV = '.tsv'
	FILE_TXT = '.txt'
	FILE_TAB = '.tab'
	FILE_VCF = '.vcf'
	FILE_VCF_BGZ = 'vcf.bgz'
	FILE_VCF_BGZ_TBI = 'vcf.bgz.tbi'
	FILE_VCF_GZ = 'vcf.gz'
	FILE_CSV = '.csv'
	FILE_PNG = '.png'
	FILE_GBK = '.gbk'
	FILE_GB = '.gb'
	FILE_FASTA = '.fasta'
	FILE_FASTQ = '.fastq'
	FILE_FNA = '.fna'
	FILE_FAA = '.faa'
	FILE_FA = '.fa'
	FILE_FAI = '.fai'
	FILE_CONSENSUS_FASTA = '.consensus.fasta'
	FILE_TREE = '.tree'
	FILE_NWK = '.nwk'
	FILE_GZ = '.gz'
	FILE_BED = '.bed'
	FILE_TBI = '.tbi'	### create with tabix
	FILE_IDX = '.idx'	### created from igvtools
	FILE_JSON = '.json'

	### all GBK
	VECT_ALL_GBK_EXTENSIONS = [FILE_GBK, FILE_GB] 
