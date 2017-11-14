'''
Created on Oct 13, 2017

@author: mmp
'''
import os
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
	
	####
	SEQ_VIRUS_TYPE = "Type"
	SEQ_VIRUS_SUB_TYPE = "Subtype"
	SEQ_VIRUS_LINEAGE = "Lineage"
	
	### data_set 
	DATA_SET_GENERIC = "Generic"	## default name for a dataset
	
	## NUMBER OF SETs to paginate
	PAGINATE_NUMBER = 15
	PAGINATE_NUMBER_SMALL = 6
	


class TypePath(Enum):
	"""
	Has the type of paths you can get from file paths
	"""
	MEDIA_ROOT = 0
	MEDIA_URL = 1
