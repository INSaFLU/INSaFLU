'''
Created on Nov 1, 2017

@author: mmp
'''

from utils.utils import Utils
from utils.result import ProcessResults, SingleResult
from constants.constants import Constants, TypeFile
from managing_files.models import Sample, Reference
from datasets.models import Consensus
from datasets.models import UploadFiles
from django.utils.translation import ugettext_lazy as _
import chardet, csv, os
from utils.data_columns import DataColumns
from utils.data_columns import NEXTSTRAIN_age
from utils.data_columns import NEXTSTRAIN_length
from utils.data_columns import NEXTSTRAIN_date
from utils.data_columns import NEXTSTRAIN_exception_nextstrain_file

class ParseNextStrainFiles(object):
	'''
	classdocs
	'''
	STATE_READ_only_detect_errors = 'only_detect_errors'
	STATE_READ_dont_detect_errors = 'dont_detect_errors'
	STATE_READ_all = 'all'

	## for nextstrain
	STATE_READ_metadata_only_detect_errors_and_chech_nexttrain = 'metadata_only_detect_errors_check_nexttrain'
	STATE_READ_metadata_dont_detect_errors_and_chech_nexttrain = 'metadata_dont_detect_errors_check_nexttrain'

	## only for metadata
	vect_only_read_metadata = [STATE_READ_metadata_only_detect_errors_and_chech_nexttrain,
				STATE_READ_metadata_dont_detect_errors_and_chech_nexttrain]
	
	utils = Utils()

	def __init__(self):
		'''
		Constructor
		'''
		self.clean_data()

	def get_errors(self):
		return self.errors
	
	def get_vect_samples(self):
		return self.vect_samples

	def get_number_samples(self):
		return self.number_samples

	def clean_data(self):
		"""
		clean data
		"""
		self.errors = ProcessResults()
		self.vect_header = []					## header, ordered
		self.dict_header = {}					## {'sample_name' : 0, 'sample_name_1' : 1, }
		self.vect_samples = []					## [[sample, vect_tag_names], [sample1, vect_tag_names1], ... ]
		
		self.dict_file_names = {}				## to test if it is a file names repeated
		self.dict_other_fields_repeated = {}	## to test if it is other fields repeated
		self.dict_samples_out = {}				## to test if it is other samples in the file
		self.dict_samples = {}
		
		self.number_samples = 0					## has the number of samples
		
	def parse_nextstrain_files(self, file_name, user, b_test_char_encoding, read_state):
		"""
		Test possible type of build:
			1) SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_ncov
			2) SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_mpx
			3) SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_generic
		if something is wrong, add the errors 
		in: dont_detect_sample_names sometimes is necessary to read the file whitout
		"""
		
		### clean data first
		self.clean_data()
		if (not os.path.exists(file_name)): return
		
		try:
			if (b_test_char_encoding):
				with open(file_name, 'rb') as f:
					result = chardet.detect(f.read())

			if (b_test_char_encoding and result['confidence'] == 1.0):
				with open(file_name, 'rt', encoding=result['encoding']) as f:
					self.process_file(f, user, read_state)
			else:
				with open(file_name) as f:
					self.process_file(f, user, read_state)
		except UnicodeDecodeError as e:
			self.errors.add_single_result(SingleResult(SingleResult.ERROR, str(e)))
		except csv.Error as e:
			self.errors.add_single_result(SingleResult(SingleResult.ERROR, 'csv or tsv file. ' + str(e)))

	def process_file(self, f, user, read_state):
		"""
		processing file of samples
		"""
		sniffer = csv.Sniffer()
		n_lines = 0
		dict_delimiter = {}
		## get biggest delimiter
		for line in f:
			if (line.startswith("#") or line.startswith("\"#")): continue
			dialect = sniffer.sniff(line)
			delimiter = dialect.delimiter
			if (delimiter in dict_delimiter): dict_delimiter[delimiter] += 1
			else: dict_delimiter[delimiter] = 1
			if (n_lines > 20): break
			n_lines += 1
			
		### get delimiter
		data_columns = DataColumns()
		delimiter = ','
		for key in dict_delimiter:
			delimiter = key
			f.seek(0)
			n_lines = 0
			b_found = False
			reader = csv.reader(f, delimiter=delimiter)
			for row in reader:
				
				(header, dict_repeated) = data_columns.get_type_header_nextstrain_and_check_repeated(row)
				if not header is None: 
					b_found = True
					break
			## header was founded				
			if (b_found): break

		## start reading
		f.seek(0)
		reader = csv.reader(f, delimiter=delimiter)
		header = None
		count_row = 1
		for row in reader:
			if (header is None):
				(header, dict_repeated) = data_columns.get_type_header_nextstrain_and_check_repeated(row)

				## test repetition in header				
				if (not header is None):
					self.vect_header = row
					self.dict_header = dict(zip((row), range(len(row))))
					for repeated_column in dict_repeated:
						self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("Column name '{}' is repeated in the header. Line: {}".\
											format(repeated_column, count_row))))

			elif (not header is None): ## line with data
				self.process_row_metadata(row, count_row, self.vect_header, user, read_state)
			count_row += 1
		if (header is None):
			self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("Header not found in the file. Please, check the names in the header, must be equal and have the same order of the template file.")))
		elif (self.number_samples == 0 and not self.errors.has_errors()):
			if (read_state in ParseNextStrainFiles.vect_only_read_metadata):
				self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("There's no samples to update.")))
			else:
				self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("There's no samples to process.")))


	def process_row_metadata(self, row, count_row, header, user, read_state):
		"""
		process the lines for update metadata
		in: read_state -> metadata_only_detect_errors_check_samples; metadata_dont_detect_errors_check_samples;
		in: only_detect_errors -> only to detect errors, not to add samples
		
				# 	STATE_READ_metadata_only_detect_errors_and_chech_samples = 'metadata_only_detect_errors_check_samples'
				# 	STATE_READ_metadata_dont_detect_errors_and_chech_samples = 'metadata_dont_detect_errors_check_samples'
				
		process header
		"""

		### to check the errors
		if (read_state == ParseNextStrainFiles.STATE_READ_metadata_only_detect_errors_and_chech_nexttrain):
			if len(row) > 0 and len(row[0].strip()) > 0:
				sample_name = row[0].strip()
				
				if not sample_name in NEXTSTRAIN_exception_nextstrain_file:
					## test clean name
					try:
						sample = Sample.objects.get(name__iexact=sample_name, owner=user, is_deleted=False)
					except Sample.DoesNotExist as error:
						
						### try repeated names of samples
						b_found = False
						if (self.utils.is_integer(sample_name.split('_')[-1])):
							repeated_name = "_".join(sample_name.split('_')[:-1])
							if len(repeated_name) > 0:
								try:
									sample = Sample.objects.get(name__iexact=repeated_name, owner=user, is_deleted=False)
									b_found = True
								except Sample.DoesNotExist as error:
									pass
						
						if (not b_found):
							try:
								consensus = Consensus.objects.get(name__iexact=sample_name, owner=user, is_deleted=False)
							except Consensus.DoesNotExist as error:
								try:
									reference = Reference.objects.get(name__iexact=sample_name, owner=user, is_deleted=False)
								except Reference.DoesNotExist as error:
									try:
										reference = Reference.objects.get(name__iexact=sample_name, owner__username=Constants.DEFAULT_USER, is_deleted=False)
									except Reference.DoesNotExist as error:
										self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("Name '{}' doesn't exists in database. Line: {} Column: {}".format(sample_name, count_row, 1))))
										pass
				
				## test repeated samples
				if (sample_name in self.dict_samples_out):
					self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("Name '{}' is repeated in the file. Line: {} Column: {}".\
											format(sample_name, count_row, 1))))
				else:
					self.dict_samples_out[sample_name] = 1
			else:
				self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("There's no sample name Line: {} Column: {}".format(count_row, 1))))

			### check if fastq1 file as gz
			if len(row) > 1:
				for i in range(1, len(header)):
					if(header[i].strip().lower() == NEXTSTRAIN_length):
						if (row[i] != "?" and len(row[i]) > 0 and not self.utils.is_integer(row[i])):
							self.errors.add_single_result(SingleResult(SingleResult.ERROR, 
							_("'{}' must be integer. Line: {} Column: {}".format(NEXTSTRAIN_length, count_row, i+1))))
					elif(header[i].strip().lower() == NEXTSTRAIN_age):
						if (row[i] != "?" and len(row[i]) > 0 and not self.utils.is_integer(row[i])):
							self.errors.add_single_result(SingleResult(SingleResult.ERROR,
							_("'{}' must be integer. Line: {} Column: {}".format(NEXTSTRAIN_age, count_row, i+1))))
					elif(header[i].strip().lower() == NEXTSTRAIN_date):
						self.validate_date(row, i, count_row, header)	## validate onset date

			else:	## there's no data
				self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("There's no data to process Line: {} Column: {}".format(count_row, 1))))
			
		self.number_samples += 1
		
		self.vect_samples.append(row[0].strip())	## add samples
		self.dict_samples[row[0].strip()] = row

	def get_value(self, sample_name, header_name):
		"""
		:out return a value from a structure, based on sample name and value
		"""
		if not header_name in self.dict_header: return None
		if not sample_name in self.dict_samples: return None
		index = self.dict_header[header_name]
		if(index < len(self.dict_samples[sample_name])): return self.dict_samples[sample_name][index]
		return None
		
		
	def validate_date(self, row, column, count_row, vect_header):
		"""
		validate date
		"""
		if len(row[column].strip()) > 0 and row[column].strip() != "?":
			try:
				return self.utils.validate_date_format(row[column].strip(), '%Y-%m-%d', 'YYYY-mm-dd')
			except ValueError as e:
				self.errors.add_single_result(SingleResult(SingleResult.ERROR, 
					_("The '{}' must have this format YYYY-MM-DD. Line: {} Column: {}".\
					format(vect_header[column], count_row, column + 1))))
