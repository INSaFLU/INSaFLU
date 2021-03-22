'''
Created on Nov 16, 2017

@author: mmp
'''
import os, gzip
from utils.utils import Utils
from utils.result import Coverage

from Bio import SeqIO

class ParseFile(object):
	'''
	classdocs
	'''
	utils = Utils()

	def __init__(self):
		'''
		Constructor
		'''
		self.data_file = None
		self.reference_dict = {}
		self.vect_reference = None
			
	def is_gzip(self, file_name): return True if (file_name.rfind(".gz") == len(file_name) - 3) else False
	
	def parse_file(self, file_name):
		"""
		"""
		self.data_file = DataFile(file_name)
		if (self.is_gzip(file_name)): handle = gzip.open(file_name, mode='rt')
		else: handle = open(file_name)
		for line in handle:
			sz_temp = line.strip().lower()
			if (len(sz_temp) == 0 or sz_temp[0] == '#'): continue
			self.data_file.add_data(line)
		handle.close()
		return self.data_file
	
		
	def read_reference_fasta(self, reference_file):
		"""
		test if the reference_file and ge the handle
		"""
		if (not os.path.exists(reference_file)): raise Exception("Can't locate the reference file: '" + reference_file + "'")

		### set temp file name
		temp_file_name = reference_file
		
		## create temp file
		b_temp_file = False
		if self.utils.is_gzip(reference_file):
			b_temp_file = True
			temp_file_name = self.utils.get_temp_file("reference_file_", ".fasta")
			cmd = "gzip -cd " + reference_file + " > " + temp_file_name
			os.system(cmd)
			
		for rec in SeqIO.parse(temp_file_name, 'fasta'):
			self.reference_dict[rec.id] = len(str(rec.seq))
			self.vect_reference.append(rec.id)
		
		### remove temp file if necessary
		if (b_temp_file): os.remove(temp_file_name)
		

class DataFile(object):
	'''
	classdocs
	'''
	util = Utils();

	def __init__(self, file_name):
		'''
		Constructor
		'''
		self.file_name = file_name
		self.vect_chromosomes = []
		self.dict_data = {}
		self.dict_data_coverage = {}
		self.previous_position = -1
		
	def get_vect_chromosomes(self): return self.vect_chromosomes
	def get_dict_data(self): return self.dict_data
	
	def add_data(self, line):
		if (len(line) == 0 or line[0] == '#'): return
		vect_data = line.split()
		if (len(vect_data) != 3): raise Exception("File: " + self.file_name + "\nThis line must have three values '" + line + "'")
		if (not self.util.is_integer(vect_data[1])): raise Exception("File: " + self.file_name + "\nLine: '" + line + "'\nThe locus need to be integer")
		if (not self.util.is_integer(vect_data[2])): raise Exception("File: " + self.file_name + "\nLine: '" + line + "'\nThe coverage need to be integer")
		if (vect_data[0] in self.dict_data): 
			if (int(vect_data[1]) <= (self.previous_position)): raise Exception("File: " + self.file_name + "\nLine: '" + line + "'\nThe locus need to be greater than the predecessor in the file")
			self.dict_data[vect_data[0]].append([vect_data[1], vect_data[2]])
			self.previous_position = int(vect_data[1])
		else:
			self.vect_chromosomes.append(vect_data[0])
			self.dict_data[vect_data[0]] = [[vect_data[1], vect_data[2]]]
			self.previous_position = int(vect_data[1])
			
		
	def get_coverage(self, sz_chromosome, length_chromosome):
		if (sz_chromosome not in self.dict_data): return 0
		if (sz_chromosome in self.dict_data_coverage): return self.dict_data_coverage[sz_chromosome]
		if (length_chromosome == 0): return 0
#		medaka sometimes creates bigger references than the original, difference 2 or 3 number of bases
		if (len(self.dict_data[sz_chromosome]) > (length_chromosome * 1.10) or 
			len(self.dict_data[sz_chromosome]) < (length_chromosome - (length_chromosome * 0.10))): 
			raise Exception("Chromosome '%s' has different sizes. Coverage: %d; Reference: %d" % (sz_chromosome, len(self.dict_data[sz_chromosome]), length_chromosome))
		sum_total = 0
		for data_ in self.dict_data[sz_chromosome]: sum_total += int(data_[1])
		self.dict_data_coverage[sz_chromosome] = sum_total / float(length_chromosome)
		return self.dict_data_coverage[sz_chromosome]
	
	def get_ratio_more_than(self, sz_chromosome, length_chromosome, value):
		if (sz_chromosome not in self.dict_data): return 0
		if (length_chromosome == 0): return 0
#		medaka sometimes creates bigger references than the original, difference 2 or 3 number of bases
		if (len(self.dict_data[sz_chromosome]) > (length_chromosome * 1.10) or 
			len(self.dict_data[sz_chromosome]) < (length_chromosome - (length_chromosome * 0.10))):
			raise Exception("Chromosome '%s' has different sizes. Coverage: %d; Reference: %d" % (sz_chromosome, len(self.dict_data[sz_chromosome]), length_chromosome))
		sum_total = 0
		for data_ in self.dict_data[sz_chromosome]: sum_total += (1 if (int(data_[1]) > value) else 0)
		return sum_total / float(length_chromosome)
	
	def get_file_name(self):
		sz_return = os.path.basename(self.file_name)
		if (sz_return.rfind(".gz") == len(sz_return) - 3): sz_return = sz_return[:-3]
		if (sz_return.rfind(".") != -1): sz_return = sz_return[:-1 * (len(sz_return) - sz_return.rfind("."))]
		return sz_return

class GetCoverage(object):
	"""
	get coverage from deep.gz file
	need deep.gz file and reference
	"""
	utils = Utils()

	def __init__(self):
		self.reference_dict = {}
		self.vect_reference = []

	def get_dict_reference(self): return self.reference_dict
	def get_vect_reference(self): return self.vect_reference

	def get_dict_with_coverage(self, deep_file):
		"""
		get a dictonary of elements with coverage 
		"""
		self.reference_dict = {}
		self.vect_reference = []

		parse_file = ParseFile()
		data_file = parse_file.parse_file(deep_file)
		dt_out = {}
		for key in data_file.get_dict_data().keys():
			dt_out[key] = [int(value_[1]) for value_ in data_file.get_dict_data()[key] ]
		return dt_out

	def get_coverage(self, deep_file, reference, limit_defined_by_user = None, limit_defined_to_project = None):
		"""
		get an instance of coverage 
		"""
		self.reference_dict = {}
		self.vect_reference = []

		parse_file = ParseFile()
		data_file = parse_file.parse_file(deep_file)
		self.read_reference_fasta(reference)

		coverage = Coverage(limit_defined_by_user, limit_defined_to_project)
		for chromosome in self.vect_reference:
			if (chromosome not in self.reference_dict): raise Exception("Can't locate the chromosome '" + chromosome + "' in reference file")
			coverage.add_coverage(chromosome, Coverage.COVERAGE_ALL, "%.1f" % (data_file.get_coverage(chromosome, self.reference_dict[chromosome])))
			coverage.add_coverage(chromosome, Coverage.COVERAGE_MORE_0, "%.1f" % (data_file.get_ratio_more_than(chromosome, self.reference_dict[chromosome], 0) * 100))
			coverage.add_coverage(chromosome, Coverage.COVERAGE_MORE_9, "%.1f" % (data_file.get_ratio_more_than(chromosome, self.reference_dict[chromosome], 9) * 100))
			if (not limit_defined_by_user is None):
				## need to decrease one value because is more than something
				coverage.add_coverage(chromosome, Coverage.COVERAGE_MORE_DEFINED_BY_USER,\
					"%.1f" % (data_file.get_ratio_more_than(chromosome, self.reference_dict[chromosome], limit_defined_by_user - 1) * 100))
			if (not limit_defined_to_project is None):
				## need to decrease one value because is more than something
				coverage.add_coverage(chromosome, Coverage.COVERAGE_PROJECT,\
					"%.1f" % (data_file.get_ratio_more_than(chromosome, self.reference_dict[chromosome], limit_defined_to_project - 1) * 100))
		return coverage


	def read_reference_fasta(self, reference_file):
		"""
		test if the reference_file and ge the handle
		"""
		if (not os.path.exists(reference_file)): raise Exception("Can't locate the reference file: '" + reference_file + "'")

		### set temp file name
		temp_file_name = reference_file
		
		## create temp file
		b_temp_file = False
		if self.utils.is_gzip(reference_file):
			b_temp_file = True
			temp_file_name = self.utils.get_temp_file("reference_file_", ".fasta")
			cmd = "gzip -cd " + reference_file + " > " + temp_file_name
			os.system(cmd)
			
		for rec in SeqIO.parse(temp_file_name, 'fasta'):
			self.reference_dict[rec.id] = len(str(rec.seq))
			self.vect_reference.append(rec.id)
		
		###
		if (b_temp_file): os.remove(temp_file_name)
