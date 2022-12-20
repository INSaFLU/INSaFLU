'''
Created on 13/06/2022

@author: mmp
'''

import csv
from datetime import date
from constants.software_names import SoftwareNames
from utils.collect_extra_data import CollectExtraData
from constants.constants import Constants
from django.conf import settings

"""
#*Mandatory metadata fields:
#strain	Sample ID (Characters “()[]{}|#><” are disallowed)
#date	YEAR-MONTH-DAY ex: 2021-02-19
#virus	ncov
#region	Africa, Asia, Europe, North America, Oceania or South America
#gisaid_epi_isl	GISAID ID; if not available needs to be “?” 
#genbank_accession	Genbank accession #; if not available needs to be “?” 
#length	Genome length; can be filled with “?” 
#segment		Filled with “genome”
#sex	host sex; if not available needs to be “?” 
#age	host age; if not available needs to be “?” 
#host	host; if not available needs to be “?”  - from ncov apparently it is not mandatory??

# id in the Sample_list.tsv corresponds to strain in nextstrain metadata
# fastq1 and fastq2 in Sample_list is to ignore
# "data set" in Sample_list is to ignore

# Lineage Pangolin (for ncov) corresponds to lineage in nextstrain metadata
"""


NEXTSTRAIN_strain = "strain"	##	Sample ID (Characters “()[]{}|#><” are disallowed)
NEXTSTRAIN_date = "date"		##	YEAR-MONTH-DAY ex: 2021-02-19
NEXTSTRAIN_virus = "virus"	  ## ncov
NEXTSTRAIN_region = "region"	##   Africa, Asia, Europe, North America, Oceania or South America
NEXTSTRAIN_gisaid_epi_isl = "gisaid_epi_isl" #	GISAID ID; if not available needs to be “?” 
NEXTSTRAIN_genbank_accession = "genbank_accession"  #  Genbank accession #; if not available needs to be “?” 
NEXTSTRAIN_accession = "accession"  #  Genbank accession #; if not available needs to be “?” 
NEXTSTRAIN_length = "length"	##   Genome length; can be filled with “?” 
NEXTSTRAIN_segment = "segment"  ##   Filled with “genome”
NEXTSTRAIN_sex = "sex"		  ## host sex; if not available needs to be “?” 
NEXTSTRAIN_age = "age"		  ## host age; if not available needs to be “?” 
NEXTSTRAIN_host = "host"		## host; if not available needs to be “?”  - from ncov apparently it is not mandatory??
NEXTSTRAIN_clade = "clade"		## host; if not available needs to be “?”  - from ncov apparently it is not mandatory??

DATASET_LIST_INSAFLU_project_name = "project_name"	## project name or empty, if not from project

VECT_NEXTSTRAIN_mandatory_ncov = [
		NEXTSTRAIN_strain,
		NEXTSTRAIN_date,
		NEXTSTRAIN_virus,
		NEXTSTRAIN_region,
		NEXTSTRAIN_gisaid_epi_isl, 
		NEXTSTRAIN_genbank_accession, 
		NEXTSTRAIN_length, 
		NEXTSTRAIN_segment,
		NEXTSTRAIN_sex, 
		NEXTSTRAIN_age, 
		NEXTSTRAIN_host,
	]

VECT_NEXTSTRAIN_mandatory_mpx = [
		NEXTSTRAIN_accession,
		NEXTSTRAIN_genbank_accession,
		NEXTSTRAIN_strain,
		NEXTSTRAIN_date,
		NEXTSTRAIN_region,
		NEXTSTRAIN_host, 
		NEXTSTRAIN_clade,
	]	

VECT_NEXTSTRAIN_mandatory_generic = [
		NEXTSTRAIN_strain,
		NEXTSTRAIN_date,
		NEXTSTRAIN_genbank_accession,		 
		NEXTSTRAIN_region,
		NEXTSTRAIN_host,
	]

DICT_MANDATORY_FIELDS = {
	SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_ncov : VECT_NEXTSTRAIN_mandatory_ncov,
	SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_mpx : VECT_NEXTSTRAIN_mandatory_mpx,
	SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_generic : VECT_NEXTSTRAIN_mandatory_generic,
	SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_flu_h3n2_12y : VECT_NEXTSTRAIN_mandatory_generic,
	SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_flu_h1n1pdm_12y : VECT_NEXTSTRAIN_mandatory_generic,
	SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_flu_vic_12y : VECT_NEXTSTRAIN_mandatory_generic,
	SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_flu_yam_12y : VECT_NEXTSTRAIN_mandatory_generic
}


DICT_NEXTSTRAIN_default_ncov = {
		NEXTSTRAIN_virus : "ncov",
		NEXTSTRAIN_region : "Europe",
		NEXTSTRAIN_gisaid_epi_isl : "?",
		NEXTSTRAIN_genbank_accession : "?",
		NEXTSTRAIN_segment : "genome",
		NEXTSTRAIN_sex : "?",
		NEXTSTRAIN_age : "?",
		NEXTSTRAIN_host : "?",
	}

DICT_NEXTSTRAIN_default_mpx = {
		NEXTSTRAIN_region : "Europe",
		NEXTSTRAIN_genbank_accession : "?",
		NEXTSTRAIN_host : "Homo sapiens",
		NEXTSTRAIN_clade : "hMPXV-1",
	}

DICT_NEXTSTRAIN_default_generic = {
		NEXTSTRAIN_genbank_accession : "?",		 
		NEXTSTRAIN_region : "Europe",
		NEXTSTRAIN_host : "?",
	}

DICT_MANDATORY_FIELDS_DEFAULTS = {
	SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_ncov : DICT_NEXTSTRAIN_default_ncov,
	SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_mpx : DICT_NEXTSTRAIN_default_mpx,
	SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_flu_h3n2_12y : DICT_NEXTSTRAIN_default_generic,
	SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_flu_h1n1pdm_12y : DICT_NEXTSTRAIN_default_generic,
	SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_flu_vic_12y : DICT_NEXTSTRAIN_default_generic,
	SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_flu_yam_12y : DICT_NEXTSTRAIN_default_generic,
	SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_generic : DICT_NEXTSTRAIN_default_generic,	
}

### if None pass
DICT_INSAFLU_to_NEXTSTRAIN = {
		'id' : NEXTSTRAIN_strain,
		'collection date': NEXTSTRAIN_date,
		'onset date': NEXTSTRAIN_date,
		'lab reception date': NEXTSTRAIN_date,
		'data set' : None,
		'fastq2' : None,
		'fastq1' : None,
	}


DICT_NEXTSTRAIN_to_INSAFLU = {
		NEXTSTRAIN_strain : ['id'],
		NEXTSTRAIN_accession : ['id'],
		NEXTSTRAIN_date : ['collection date', 'onset date', 'lab reception date'],
}


##### Important, array with exceptions for nextstrain TSV file 
NEXTSTRAIN_exception_nextstrain_file = [
	'Wuhan/Hu-1/2019',
	'MK783032',
	'MK783030',
]

class MetaRow(object):
	
	def __init__(self, project_name, seq_name_consensus, row, consensus_length):
		"""
		:param project_name, has the project name for a consensus that came from a project
							otherwise has the name of the reference or the consensus
		"""
		self.seq_name_consensus = seq_name_consensus
		self.row = row
		self.consensus_length = consensus_length
		self.project_name = project_name if not project_name is None and len(project_name) > 0 else "?" 

class Metadata(object):
	
	def __init__(self, header):
		'''
		Constructor
		'''
		self.header = header
		self.dt_header = dict(zip(header, list(range(0, len(header)))))
		self.vect_rows_id = []
		self.dt_rows_id = {}
		
	def add_metadata(self, project_sample_pk, project_name, seq_name_consensus, row, consensus_length):
		"""
		:param project_name, has the project name for a consensus that came from a project
							otherwise has the name of the reference or the consensus
		"""
		
		if not project_sample_pk in self.dt_rows_id:
			self.dt_rows_id[project_sample_pk] = MetaRow(project_name, seq_name_consensus, row, consensus_length)
	
	def get_vect_out(self, vect_header_out, csv_writer):
	
		dt_out_id_project_sample = {}
		count = 0
		for project_sample_pk in self.dt_rows_id:
			if project_sample_pk in dt_out_id_project_sample: continue
			dt_out_id_project_sample[project_sample_pk] = 1
			vect_out = [self.dt_rows_id[project_sample_pk].seq_name_consensus]
			vect_out.append(self.dt_rows_id[project_sample_pk].project_name)
				
			for column in vect_header_out:
				if column == CollectExtraData.HEADER_SAMPLE_OUT_ID: continue 
				if column == DATASET_LIST_INSAFLU_project_name: continue
				if column in self.dt_header: vect_out.append(self.dt_rows_id[project_sample_pk].row[self.dt_header[column]])
				else: vect_out.append("")
			count += 1 
			csv_writer.writerow(vect_out)
		return count
	
	def get_vect_out_nextstrain(self, vect_header_out, dt_header_normal_out, csv_writer, build, parse_in_files):
		"""
		:param vect_header_out all headers out
		:param dt_header_normal_out keys that are present only in INSAFLu list files
		:param csv_writer 
		"""
	
		dt_out_id_project_sample = {}
		count = 0
		for project_sample_pk in self.dt_rows_id:

			# Avoid repeating samples
			if project_sample_pk in dt_out_id_project_sample: continue

			dt_out_id_project_sample[project_sample_pk] = 1

			## NEXTSTRAIN_strain
			vect_out = [self.dt_rows_id[project_sample_pk].seq_name_consensus]

			dt_out_header = {}
			for column in vect_header_out:

				if column == NEXTSTRAIN_strain: continue	## already out
				if column in dt_out_header: continue		## already out, can be synonymous

				### set nextstrain metadata, uploaded by file
				value = parse_in_files.get_value(vect_out[0], column)
				if not value is None:
					vect_out.append(value)
					continue
				
				## exception
				if column == NEXTSTRAIN_length:
					if self.dt_header.get(column, -1) > 0 and \
						len(self.dt_rows_id[project_sample_pk].row[self.dt_header[column]]) > 0:
						vect_out.append(self.dt_rows_id[project_sample_pk].row[self.dt_header[column]])
					elif not self.dt_rows_id[project_sample_pk].consensus_length is None:
						vect_out.append(str(self.dt_rows_id[project_sample_pk].consensus_length))
					else: vect_out.append('?')
					continue

				### date column, has data
				if (column == NEXTSTRAIN_date and self.dt_header.get(column, -1) > 0 and \
						len(self.dt_rows_id[project_sample_pk].row[self.dt_header[column]]) > 0):
					vect_out.append(self.dt_rows_id[project_sample_pk].row[self.dt_header[column]])
					continue
				
				### test synonymous, And try date synonymous
				b_found = False
				for column_insaflu in DICT_NEXTSTRAIN_to_INSAFLU.get(column, []):	 
					if self.dt_header.get(column, -1) == -1 and column_insaflu in self.dt_header and \
						len(self.dt_rows_id[project_sample_pk].row[self.dt_header[column_insaflu]]) > 0:
						vect_out.append(self.dt_rows_id[project_sample_pk].row[self.dt_header[column_insaflu]])
						for column_insaflu in DICT_NEXTSTRAIN_to_INSAFLU[column]: dt_out_header[column_insaflu] = 1
						b_found = True
						break
				
				## found synonymous before
				if b_found: continue
				if (column == NEXTSTRAIN_date):	 ## need to add default date
					vect_out.append(date.today().strftime(settings.DATE_FORMAT_FOR_SHOW))
					continue
				
				## test default NEXTstrain columns names
				if column in DICT_MANDATORY_FIELDS[build]:
					if self.dt_header.get(column, -1) > 0 and \
						len(self.dt_rows_id[project_sample_pk].row[self.dt_header[column]]) > 0:
						vect_out.append(self.dt_rows_id[project_sample_pk].row[self.dt_header[column]])
					else:
						vect_out.append(DICT_MANDATORY_FIELDS_DEFAULTS[build][column])
					dt_out_header[column] = 1
					continue
				
				### regular INSAFLU
				if column in self.dt_header: vect_out.append(self.dt_rows_id[project_sample_pk].row[self.dt_header[column]])
				else: vect_out.append("?")
			count += 1
			
			csv_writer.writerow(['?' if len(_) == 0 else _ for _ in vect_out])

		return count
	
	

class Reference(object):
	'''
	This is to represent build-specific references which have their own metadata tsv
	'''

	def __init__(self, file_name):
		self.file_name = file_name
		self._read_file()

	def _read_file(self):

		self.dt_out_rows = {}
		with open(self.file_name) as handle_in: 
			csv_reader = csv.reader(handle_in, delimiter=Constants.SEPARATOR_TAB)
			for line, row in enumerate(csv_reader):
				if line == 0: self.dt_header = dict(zip(row, list(range(0, len(row)))))
				else: self.dt_out_rows[row[0]] = row
		## done

	def save_out_nextstrain(self, csv_writer, vect_header_out, build, parse_in_files):
		"""
		"""
		
		count = 0
		for ref_id in self.dt_out_rows:
			vect_out = []
			for index, column in enumerate(DICT_MANDATORY_FIELDS[build]):
				if (self.dt_header.get(column, -1) < 0): vect_out.append(DICT_MANDATORY_FIELDS_DEFAULTS[build].get(column, '?'))
				else: vect_out.append(self.dt_out_rows[ref_id][self.dt_header[column]])
			vect_out += ['?'] * (len(vect_header_out) - index - 1)   
			count += 1
			csv_writer.writerow(vect_out)
		return count


class DataColumns(object):
	'''
	Create table for nextstrain metadata. 
	Depends on the nextstrain build
	'''

	def __init__(self, build = None):
		'''
		Constructor
		'''
		self.build = build
		self.dt_project = {}
	
	def add_header(self, project_pk, header):
		"""
		add header for this project
		"""
	
		if not project_pk in self.dt_project: 
			self.dt_project[project_pk] = Metadata(header)
		
	def add_metadata(self, project_pk, project_sample_pk, project_name, seq_name_consensus, row, consensus_length = None):
		"""
		add metadata for a specific project,project_sample
		"""
		if project_pk in self.dt_project:
			self.dt_project[project_pk].add_metadata(project_sample_pk, project_name, seq_name_consensus, row, consensus_length)
		
	def get_type_header_nextstrain_and_check_repeated(self, vect_header):
		"""
		get type header and repeated
		:out TypeHeader, Dict with repeated
		:
		"""
		## check repeated
		dict_out_value = {}
		for value in vect_header:
			if value in dict_out_value: dict_out_value[value] += 1
			else: dict_out_value[value] = 1
		dict_repeated = { key:dict_out_value[key] for key in dict_out_value if dict_out_value[key] > 1 }
		
		if self.__test_type_header(vect_header, VECT_NEXTSTRAIN_mandatory_ncov):
			return SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_ncov, dict_repeated
		if self.__test_type_header(vect_header, VECT_NEXTSTRAIN_mandatory_mpx):
			return SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_mpx, dict_repeated
		if self.__test_type_header(vect_header, VECT_NEXTSTRAIN_mandatory_generic):
			return SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_generic, dict_repeated
		return None, dict_repeated

	def __test_type_header(self, vect_header, vect_reference):
		"""
		:param vect_header		-> pass the header to test
		:param vect_reference	-> vector to find match
		Test if the header passed is equal to a specific vect, like:
		1) VECT_NEXTSTRAIN_mandatory_ncov
		2) VECT_NEXTSTRAIN_mandatory_mpx
		3) VECT_NEXTSTRAIN_mandatory_generic
		"""
		
		dict_to_test = dict(zip(vect_reference, [0] * len(vect_reference)))
		for value in vect_header:
			if value in dict_to_test: dict_to_test[value] += 1
		return len([value for value in dict_to_test if dict_to_test[value] > 0]) == len(vect_reference)
		
	def _get_header(self):
		"""
		return header
		## CollectExtraData.HEADER_SAMPLE_OUT_CSV_simple
		"""
		
		self.vect_header_out = []
		dt_header_out = {}
		
		if len(self.dt_project) > 0:
			for key_metadata in self.dt_project:
				if len(self.vect_header_out) == 0: 
					self.vect_header_out = self.dt_project[key_metadata].header
					## set CollectExtraData.HEADER_SAMPLE_OUT_ID be the first
					try:
						self.vect_header_out.pop(self.vect_header_out.index(CollectExtraData.HEADER_SAMPLE_OUT_ID))
					except ValueError as e:
						pass
					self.vect_header_out = [CollectExtraData.HEADER_SAMPLE_OUT_ID, DATASET_LIST_INSAFLU_project_name] +\
						self.vect_header_out
					dt_header_out = dict(zip(self.vect_header_out, [1] * len(self.vect_header_out)))
				else:
					for column in self.dt_project[key_metadata].header:
						if not column in dt_header_out:
							self.vect_header_out.append(column)
							dt_header_out[column] = 1
		return self.vect_header_out


	def _get_header_nextstrain(self, parse_in_files):
		"""
		return header
		## CollectExtraData.HEADER_SAMPLE_OUT_CSV_simple
		"""

		regular_header = self._get_header()
		self.vect_header_out = []
		self.dt_header_normal_out = {}

		## mandatory fields
		self.vect_header_out = DICT_MANDATORY_FIELDS[self.build].copy()

		## try to find others in the other file
		for regular_name in regular_header:
			if regular_name in DATASET_LIST_INSAFLU_project_name: continue	## exclude project_name of regular header
			if regular_name in DICT_INSAFLU_to_NEXTSTRAIN: continue
			if regular_name in self.vect_header_out: continue
			self.vect_header_out.append(regular_name)
			self.dt_header_normal_out[regular_name] = 1
		
		### check remove or add columns
		if not parse_in_files is None and len(parse_in_files.vect_header) > 0:
			for value in parse_in_files.vect_header:
				if value in DICT_MANDATORY_FIELDS[self.build]: continue
				
				### add if not there
				if not value in self.vect_header_out: 
					self.vect_header_out.append(value)
					self.dt_header_normal_out[value] = 1
			
			### NOT to remove because cna bring the case that a column never appears		
			### remove if not in 
			vect_index_to_remove = []
			for index, value in enumerate(self.vect_header_out):
				if value in DICT_MANDATORY_FIELDS[self.build]: continue
				if not value in parse_in_files.vect_header:
					vect_index_to_remove.append(index)
			
			## reverse
			#vect_index_to_remove = vect_index_to_remove[::-1]
			#for index in vect_index_to_remove:
			#	del self.dt_header_normal_out[self.vect_header_out[index]]
			#	self.vect_header_out.pop(index)
		return self.vect_header_out

	
	def save_rows(self, csv_writer):
		"""
		save all data
		"""
		count = 0 ## number of rows saved
		
		## header save
		csv_writer.writerow(self._get_header())
		
		## save data
		for key_metadata in self.dt_project:
			count += self.dt_project[key_metadata].get_vect_out(
				self.vect_header_out, csv_writer) 
		return count


	def save_rows_nextstrain(self, csv_writer, reference_tsv, parse_in_files):
		"""
		:param reference_tsv has the data that has the template to save this
		:param parse_in_files has thedata that came from metadata files imported by users
		create file for nextstrain
		"""
		count = 0 ## number of rows saved
		
		## save header
		csv_writer.writerow(self._get_header_nextstrain(parse_in_files))

		# read NextStrain reference (obsolete)
		if(not reference_tsv is None):
			reference = Reference(reference_tsv)
			### save reference
			count += reference.save_out_nextstrain(csv_writer, self.vect_header_out, self.build, parse_in_files)

		## save data
		for key_metadata in self.dt_project:
			
			count += self.dt_project[key_metadata].get_vect_out_nextstrain(
					self.vect_header_out, self.dt_header_normal_out, csv_writer, self.build, parse_in_files)			

		return count
		
			