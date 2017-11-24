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
	DEFAULT_USER = "default_user"
	DEFAULT_USER_PASS = "default_user_123_$%_2"
	
	META_KEY_VALUE_NOT_NEED = "value not needed"
	
	### Session variables
	NUMBER_LOCUS_FASTA_FILE = "number_locus_fasta_file"
	
	## main path for all paths
	MAIN_PATH = "/usr/local/insaflu"
	

	DIR_PROCESSED_FILES_UPLOADS = "uploads"
	
	## DIR_PROCESSED_FILES_FROM_WEB/userId_<id>/refId_<id>
	DIR_PROCESSED_FILES_REFERENCE = DIR_PROCESSED_FILES_UPLOADS + "/references"
	DIR_PROCESSED_FILES_FASTQ = DIR_PROCESSED_FILES_UPLOADS + "/fastq"
	DIR_PROCESSED_FILES_PROJECT = DIR_PROCESSED_FILES_UPLOADS + "/project"

	DIR_ICONS = "icons"
	TEMP_DIRECTORY = "/tmp"
	COUNT_DNA_TEMP_DIRECTORY = "insaFlu"

	FORMAT_FASTA = "fasta"
	FORMAT_FASTQ = "fastq"
	EXTENSION_ZIP = ".gz"
	
	## Has all the versions of type influenza typing
	DIR_TYPE_IDENTIFICATION = "db/type_identification/"
	
	INSAFLU_NAME = 'insaflu'
	
	####
	SEQ_VIRUS_TYPE = "Type"
	SEQ_VIRUS_SUB_TYPE = "Subtype"
	SEQ_VIRUS_LINEAGE = "Lineage"
	
	### data_set 
	DATA_SET_GENERIC = "Generic"	## default name for a dataset
	
	## NUMBER OF SETs to paginate
	PAGINATE_NUMBER = 15
	PAGINATE_NUMBER_SMALL = 2

	## tag for check box all in the tables
	CHECK_BOX_ALL = 'check_box_all'
	CHECK_BOX = 'check_box'
	GET_CHECK_BOX_SINGLE = 'get_check_box_single'
	COUNT_CHECK_BOX = 'count_check_boxes'
	CHECK_BOX_VALUE = 'value'
		
	def get_extensions_by_file_type(self, file_name, file_type):
		"""
		get extensions by file type
		"""
		if (file_type == FileType.FILE_BAM): return "{}.bam".format(file_name)
		if (file_type == FileType.FILE_BAM_BAI): return "{}.bam.bai".format(file_name)
		if (file_type == FileType.FILE_CONSENSUS_FA): return "{}.consensus.fa".format(file_name)
		if (file_type == FileType.FILE_CONSENSUS_FASTA): return "{}.consensus.fasta".format(file_name)
		if (file_type == FileType.FILE_CSV): return "{}.csv".format(file_name)
		if (file_type == FileType.FILE_DEPTH): return "{}.depth".format(file_name)
		if (file_type == FileType.FILE_DEPTH_GZ): return "{}.depth.gz".format(file_name)
		if (file_type == FileType.FILE_DEPTH_GZ_TBI): return "{}.depth.gz.tbi".format(file_name)
		if (file_type == FileType.FILE_TAB): return "{}.tab".format(file_name)
		if (file_type == FileType.FILE_VCF): return "{}.vcf".format(file_name)
		if (file_type == FileType.FILE_VCF_GZ): return "{}.vcf.gz".format(file_name)
		if (file_type == FileType.FILE_VCF_GZ_TBI): return "{}.vcf.gz.tbi".format(file_name)
		return ""
		

			
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
	




