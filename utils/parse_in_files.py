'''
Created on Nov 1, 2017

@author: mmp
'''

from utils.utils import Utils
from utils.result import ProcessResults, SingleResult
from constants.constants import Constants, TypeFile
from managing_files.models import Sample, DataSet, VaccineStatus, TagName, TagNames, UploadFiles
from django.utils.translation import ugettext_lazy as _
from django.contrib.gis.geos import Point
from constants.constants import FileExtensions
import chardet
import csv

class ParseInFiles(object):
	'''
	classdocs
	'''

	## this header must be present
	vect_header = ['sample name', 'fastq1', 'fastq2', 'data set', 'vaccine status', 'week', 'onset date', 'collection date', 'lab reception date', 'latitude', 'longitude']
	vect_madatory_header = ['sample name', 'fastq1']
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

	def clean_data(self):
		"""
		clean data
		"""
		self.errors = ProcessResults()
		self.vect_samples = []					## [[sample, vect_tag_names], [sample1, vect_tag_names1], ... ]
		
		self.dict_file_names = {}				## to test if it is a file names repeated
		self.dict_other_fields_repeated = {}	## to test if it is other fields repeated
		self.dict_samples_out = {}				## to test if it is other samples in the file

	def parse_sample_files(self, file_name, user, b_test_char_encoding, only_detect_errors):
		"""
		#The first eleven header fields are mandatory, then you can add what you want as a field											
		#Mandatory fields: 'sample name', 'fastq1'										
		#Date format: dd/mm/yyyy											
		#Fastq file names must be equal to the files that you are going to upload later											
		#Sample names must be unique even with the ones that are in database											
		sample name 	fastq1	fastq2	data set	vaccine status	week	onset date	collection date	lab reception date	latitude	longitude	other fields from here
		
		return a instance of a class with all lines
		if something is wrong, add the errors 
		"""
		
		### clean data first
		self.clean_data()
		
		try:
			if (b_test_char_encoding):
				with open(file_name, 'rb') as f:
					result = chardet.detect(f.read())

			if (b_test_char_encoding and result['confidence'] == 1.0):
				with open(file_name, 'rt', encoding=result['encoding']) as f:
					self.process_file(f, user, only_detect_errors)
			else:
				with open(file_name) as f:
					self.process_file(f, user, only_detect_errors)
		except UnicodeDecodeError as e:
			self.errors.add_single_result(SingleResult(SingleResult.ERROR, str(e)))

	def process_file(self, f, user, only_detect_errors):
		"""
		"""
		sniffer = csv.Sniffer()
		n_lines = 0
		delimiter = ','
		for line in f:
			dialect = sniffer.sniff(line)
			delimiter = dialect.delimiter
			if (n_lines > 15): break

		f.seek(0)
		reader = csv.reader(f, delimiter=delimiter)
		header = None
		(count_row, count_column) = (1, 1)
		for row in reader:
			if (header == None):
				count_column = 0
				n_columns_ok = 0
				for col in row:
					if (count_column >= len(self.vect_header) or col.replace(' ', '').lower() != self.vect_header[count_column].replace(' ', '').lower()): break
					n_columns_ok += 1
					count_column += 1
				if (n_columns_ok == len(self.vect_header)):
					header = row
					for i in range(0, len(row)):
						if (row[i].strip() in self.dict_other_fields_repeated):
							self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("Column name '{}' is repeated in the header. Line: {} Column: {}".\
										format(row[i].strip(), count_row, i+1))))
						else:
							self.dict_other_fields_repeated[row[i].strip()] = 1

			elif (header != None): ## line with data
				self.process_row(row, count_row, header, user, only_detect_errors)
			count_row += 1
		if (header == None):
			self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("Header not found in the file. Please, check the names in the header, must be equal and have the same order of the template file.")))
		elif (len(self.vect_samples) == 0 and not self.errors.has_errors()):
			self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("There's no samples to process.")))

	def process_row(self, row, count_row, header, user, only_detect_errors):
		"""
		in: only_detect_errors -> only to detect errors, not to add samples
		process header
		"""
		### check sample name
		n_errors = len(self.errors.get_vect_results())
		if len(row) > 0:
			sample_name = row[0].strip()
			try:
				sample = Sample.objects.get(name__iexact=sample_name, owner=user)
				self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("Sample name '{}' exists in database. Line: {} Column: {}".format(sample_name, count_row, 1))))
			except Sample.DoesNotExist as e:
				pass
			
			## test repeated samples
			if (sample_name in self.dict_samples_out):
				self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("Sample name '{}' is repeated in the file. Line: {} Column: {}".\
										format(sample_name, count_row, 1))))
			else:
				self.dict_samples_out[sample_name] = 1
		else:
			self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("There's no sample name Line: {} Column: {}".format(count_row, 1))))
			
		### check if fastq1 file as gz
		if len(row) > 1:
			fastq1 = row[1]
			if (not fastq1.strip().endswith(FileExtensions.FILE_GZ)):
				self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("File '{}' not ends with '{}'. Fastq gzip compressed file is needed. Line: {} Column: {}".\
										format(fastq1.strip(), FileExtensions.FILE_GZ, count_row, 2))))
			elif (fastq1.strip() in self.dict_file_names):
				self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("File '{}' is repeated in the samples file. Line: {} Column: {}".\
										format(fastq1.strip(), count_row, 2))))
			else:
				number_files = UploadFiles.objects.filter(file_name__iexact=fastq1.strip(), owner=user,\
									is_processed=False, type_file__name=TypeFile.TYPE_FILE_fastq_gz).count()
				if (number_files > 0):
					self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("File '{}' is repeated in the database and it's not processed yet. Line: {} Column: {}".\
																format(fastq1.strip(), count_row, 2))))
				else:
					self.dict_file_names[fastq1.strip()] = 1
		else:
			self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("There's no fastq1 file name. Line: {} Column: {}".format(count_row, 2))))

		### check if fastq2 file as gz
		if len(row) > 2 and len(row[2].strip()) > 0:
			fastq2 = row[2]
			if (not fastq2.endswith(FileExtensions.FILE_GZ)):
				self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("File '' not ends with '{}'. Fastq gzip compressed file is needed. Line: {} Column: {}".\
										format(fastq1, FileExtensions.FILE_GZ, count_row, 3))))
			elif (fastq2.strip() in self.dict_file_names):
				self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("File '{}' is repeated in the samples file. Line: {} Column: {}".\
										format(fastq2.strip(), count_row, 3))))
			else:
				number_files = UploadFiles.objects.filter(file_name__iexact=fastq2.strip(), owner=user,\
									is_processed=False, type_file__name=TypeFile.TYPE_FILE_fastq_gz).count()
				if (number_files > 0):
					self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("File '{}' is repeated in the database and it's not processed yet. Line: {} Column: {}".\
																format(fastq1.strip(), count_row, 3))))
				else:
					self.dict_file_names[fastq2.strip()] = 1

		### anything to do with data set and vaccine status, are free fields
		
		### week
		if len(row) > 5:
			week = row[5]
			if (len(week) > 0 and not self.utils.is_integer(week)):
				self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("'week' must be integer. Line: {} Column: {}".format(count_row, 3))))

		self.validate_date(row, 6, count_row)	## validate onset date
		self.validate_date(row, 7, count_row)	## validate collection date
		self.validate_date(row, 8, count_row)	## validate lab reception date
		
		### latitude
		if len(row) > 9 and len(row[9].strip()) > 0:
			latitude = row[9].strip()
			if (not self.utils.is_float(latitude) or float(latitude) > 90 or float(latitude) < -90):
				self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("'latitude' must have values between -90<lat<90. Line: {} Column: {}".format(count_row, 3))))

		### longitude
		if len(row) > 10 and len(row[10].strip()) > 0:
			longitude = row[10].strip()
			if (not self.utils.is_float(longitude) or float(longitude) > 100 or float(longitude) < -100):
				self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("'longitude' must have values between -180<long<180. Line: {} Column: {}".format(count_row, 3))))

		### there's no errors, process this sample
		if (n_errors == len(self.errors.get_vect_results()) and not only_detect_errors):
			
			### add samples
			sample = Sample()
			sample.id = count_row
			sample.owner = user
			sample.name = row[0]
			sample.file_name_1 = row[1].strip()
			if (len(row[2].strip()) > 0): sample.file_name_2 = row[2].strip()
			sample.data_set = DataSet()
			sample.data_set.owner = user
			if (len(row[3].strip()) > 0): sample.data_set.name = row[3].strip()
			else: sample.data_set.name = Constants.DATA_SET_GENERIC	
			if (len(row[4].strip()) > 0):
				sample.vaccine_status = VaccineStatus()
				sample.vaccine_status.name = row[4].strip()
				sample.vaccine_status.owner = user
			if (len(row[5].strip()) > 0): sample.week = int(row[5].strip())
			
			## dates
			if (len(row[6].strip()) > 0): sample.date_of_onset = self.utils.validate_date(row[6].strip())
			if (len(row[7].strip()) > 0): sample.date_of_collection = self.utils.validate_date(row[7].strip())
			if (len(row[8].strip()) > 0): sample.date_of_receipt_lab = self.utils.validate_date(row[8].strip())
			
			if len(row[9].strip()) > 0:
				sample.geo_local = Point(float(row[9].strip()), float(row[10].strip()))
			
			vect_tag_names = []	## this trick need to be done because 'Extra fields on many-to-many'
			for i in range(11, len(header)):
				if len(row[i].strip()) > 0:
					tag_name = TagName()
					tag_name.id = count_row
					tag_name.owner = user
					tag_name.name = header[i].strip()
					tag_name.is_meta_data = False
					
					tag_names = TagNames()
					tag_names.id = count_row
					tag_names.value = row[i].strip()
					tag_names.tag_name = tag_name
					tag_names.sample = sample
					vect_tag_names.append(tag_names) ## this trick need to be done because 'Extra fields on many-to-many'
			self.vect_samples.append([sample, vect_tag_names])	## add samples


	def validate_date(self, row, column, count_row):
		"""
		validate date
		"""
		if len(row) > 6 and len(row[column].strip()) > 0:
			try:
				return self.utils.validate_date(row[column].strip())
			except ValueError as e:
				self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("The '{}' must have this format DD/MM/YYYY. Line: {} Column: {}".\
						format(self.vect_header[column], count_row, column + 1))))

