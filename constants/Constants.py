'''
Created on Oct 13, 2017

@author: mmp
'''
from Bio import SeqIO

class Constants(object):
	'''
	classdocs
	'''
	META_KEY_VALUE_NOT_NEED = "value not needed"
	
	### Session variables
	NUMBER_LOCUS_FASTA_FILE = "number_locus_fasta_file"

	def __init__(self):
		'''
		Constructor
		'''
		pass
	
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
				raise IOError("Error: the file is not in FASTA format.")
		handle.close()
		if (not b_pass): raise IOError("Error: the file is not in FASTA format.")

		record_dict = SeqIO.index(sz_file_name, "fasta")
		if (len(record_dict) > 0): return len(record_dict)
		raise IOError("Error: the file is not in FASTA format.")
	
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
				raise IOError("Error: the file is not in GenBank format.")
		handle.close()
		if (not b_pass): raise IOError("Error: the file is not in GenBank format.")

		record_dict = SeqIO.index(sz_file_name, "genbank")
		if (len(record_dict) > 0): return len(record_dict)
		raise IOError("Error: the file is not in GenBank format.")

