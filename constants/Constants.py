'''
Created on Oct 13, 2017

@author: mmp
'''
from Bio import SeqIO
from django.utils.translation import ugettext_lazy as _
import os, random, gzip
import logging

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
	
	## DIR_PROCESSED_FILES_PROCESSED/user/fasta/day/month/year/
	## DIR_PROCESSED_FILES_PROCESSED/user/project_<id_db>/
	DIR_PROCESSED_FILES_PROCESSED =  MAIN_PATH + "/processed"
	## DIR_PROCESSED_FILES_FROM_WEB/user/fasta/day/month/year/
	DIR_PROCESSED_FILES_UPLOADS = "uploads"
	
	## DIR_PROCESSED_FILES_FROM_WEB/userId_<id>/refId_<id>
	DIR_PROCESSED_FILES_REFERENCE = DIR_PROCESSED_FILES_UPLOADS + "/references"

	DIR_ICONS = "icons"
	TEMP_DIRECTORY = "/tmp"
	COUNT_DNA_TEMP_DIRECTORY = "insaFlu"

	FORMAT_FASTA = "fasta"
	FORMAT_FASTQ = "fastq"
	EXTENSION_ZIP = ".gz"
	
	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	def __init__(self):
		'''
		Constructor
		'''
		pass

	def get_path_to_reference_file(self, user_id, ref_id):
		"""
		get the path and reference file
		"""
		return os.path.join(self.DIR_PROCESSED_FILES_REFERENCE, "userId_{0}".format(user_id), "refId_{0}".format(ref_id))

	def get_temp_file(self, file_name, sz_type):
		"""
		return a temp file name
		"""
		main_path = os.path.join(self.TEMP_DIRECTORY, self.COUNT_DNA_TEMP_DIRECTORY)
		if (not os.path.exists(main_path)): os.makedirs(main_path)
		while 1:
			return_file = os.path.join(main_path, "count_dna_" + file_name + "_" + str(random.randrange(10000000, 99999999, 10)) + "_file" + sz_type)
			if (not os.path.exists(return_file)): return return_file


	def remove_temp_file(self, sz_file_name):
		"""
		prevent to remove files outside of temp directory
		"""
		if os.path.exists(sz_file_name) and len(sz_file_name) > 0 and sz_file_name.find(self.TEMP_DIRECTORY) == 0:
			cmd = "rm " + sz_file_name
			exist_status = os.system(cmd)
			if (exist_status != 0):
				self.logger_production.error('Fail to run: ' + cmd)
				self.logger_debug.error('Fail to run: ' + cmd)
				raise Exception("Fail to remove a file") 

	def move_file(self, sz_file_from, sz_file_to):
		if os.path.exists(sz_file_from):
			self.make_path(os.path.dirname(sz_file_to))
			cmd = "mv " + sz_file_from + " " + sz_file_to
			exist_status = os.system(cmd)
			if (exist_status != 0):
				self.logger_production.error('Fail to run: ' + cmd)
				self.logger_debug.error('Fail to run: ' + cmd)
				raise Exception("Fail to make a move a file") 
			
	def make_path(self, path_name):
		if (not os.path.isdir(path_name) and not os.path.isfile(path_name)):
			cmd = "mkdir -p " + path_name
			os.system(cmd)
			exist_status = os.system(cmd)
			if (exist_status != 0):
				self.logger_production.error('Fail to run: ' + cmd)
				self.logger_debug.error('Fail to run: ' + cmd)
				raise Exception("Fail to make a path") 
		
	def is_integer(self, n_value):
		try:
			int(n_value)
			return True
		except ValueError: 
			return False


	def is_float(self, n_value):
		try:
			float(n_value)
			return True
		except ValueError: 
			return False

	def is_gzip(self, file_name):
		"""
		test if the file name ends in gzip
		""" 
		return True if (file_name.rfind(".gz") == len(file_name) - 3) else False
	
	
	def get_type_file(self, file_name):
		"""
		return 'fasta' or 'fastq' 
		raise exception if can't detected
		"""
		if (self.is_gzip(file_name)): handle = gzip.open(file_name)
		else: handle = open(file_name)
		for record in SeqIO.parse(handle, self.FORMAT_FASTQ):
			handle.close() 
			return self.FORMAT_FASTQ
		handle.close()
		
		if (self.is_gzip(file_name)): handle = gzip.open(file_name)
		else: handle = open(file_name)
		for record in SeqIO.parse(handle, self.FORMAT_FASTA):
			handle.close() 
			return self.FORMAT_FASTA
		handle.close()
		
		raise Exception("Can't detect file format for the file '" + file_name + "'")
	

	def is_fasta(self, sz_file_name):
		"""
		Test Fata file
		"""
		handle = open(sz_file_name)
		b_pass = False
		for line in handle:
			sz_temp = line.strip()
			if (len(sz_temp) == 0): continue
			if (sz_temp[0] == ">"): 
				b_pass = True
				break
			else: 
				handle.close()
				raise IOError(_("Error: the file is not in FASTA format."))
		handle.close()
		if (not b_pass): raise IOError(_("Error: the file is not in FASTA format."))

		record_dict = SeqIO.index(sz_file_name, "fasta")
		if (len(record_dict) > 0): return len(record_dict)
		raise IOError(_("Error: the file is not in FASTA format."))

	
	def is_genbank(self, sz_file_name):
		"""
		Test GenBank file
		"""
		handle = open(sz_file_name)
		b_pass = False
		for line in handle:
			sz_temp = line.strip()
			if (len(sz_temp) == 0): continue
			if (sz_temp.startswith("LOCUS", 0)): 
				b_pass = True
				break
			else: 
				handle.close()
				raise IOError(_("Error: the file is not in GenBank format."))
		handle.close()
		if (not b_pass): raise IOError(_("Error: the file is not in GenBank format."))

		n_number_locus = 0
		for record in SeqIO.parse(sz_file_name, "genbank"):
			n_number_locus += 1
		if (n_number_locus > 0): return n_number_locus
		raise IOError(_("Error: the file is not in GenBank format."))


	
	def compare_locus_fasta_gb(self, fasta_file, gb_file):
		"""
		Test if all fasta locus are in gb file
		"""
		locus_fasta = self.is_fasta(fasta_file)
		locus_gb = self.is_genbank(gb_file)
		if (locus_gb != locus_fasta): raise ValueError(_("Number of locus are different from fasta to genbank."))
		
		record_dict = SeqIO.index(fasta_file, "fasta")
		for record in SeqIO.parse(gb_file, "genbank"):
			b_found = False
			for seq in record_dict:
				if (seq == record.name):
					if (len(record_dict[seq].seq) != len(record.seq)):
						raise ValueError(_("Different length. Fasta seq: %s length: %d; Fasta seq: %s length: %d." 
									% (seq, len(record_dict[seq].seq), record.name, len(record.seq))) )
					b_found = True
					break
			if (not b_found): raise ValueError(_("This locus '" + record.locus + "' is not in fasta file."))
		return locus_fasta


