'''
Created on Nov 1, 2017

@author: mmp
'''

from django.db import transaction
from utils.utils import Utils
from utils.result import ProcessResults, SingleResult
from constants.constants import Constants, TypeFile, TypePath
from managing_files.models import Sample, DataSet, VaccineStatus, TagName, TagNames, UploadFiles, ProcessControler
from django.utils.translation import ugettext_lazy as _
from django.contrib.gis.geos import Point
from constants.constants import FileExtensions
import chardet, csv, os, re, time
from django.conf import settings
from django_q.tasks import async
from utils.software import Software
from constants.meta_key_and_values import MetaKeyAndValue
from managing_files.manage_database import ManageDatabase
from datetime import datetime
from utils.process_SGE import ProcessSGE

class ParseInFiles(object):
	'''
	classdocs
	'''
	STATE_READ_only_detect_errors = 'only_detect_errors'
	STATE_READ_dont_detect_errors = 'dont_detect_errors'
	STATE_READ_all = 'all'
	
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

	def get_number_samples(self):
		return self.number_samples

	def clean_data(self):
		"""
		clean data
		"""
		self.errors = ProcessResults()
		self.vect_samples = []					## [[sample, vect_tag_names], [sample1, vect_tag_names1], ... ]
		
		self.dict_file_names = {}				## to test if it is a file names repeated
		self.dict_other_fields_repeated = {}	## to test if it is other fields repeated
		self.dict_samples_out = {}				## to test if it is other samples in the file

		self.number_samples = 0					## has the number of samples
		
	def parse_sample_files(self, file_name, user, b_test_char_encoding, read_state):
		"""
		#The first eleven header fields are mandatory, then you can add what you want as a field											
		#Mandatory fields: 'sample name', 'fastq1'										
		#Date format: dd/mm/yyyy											
		#Fastq file names must be equal to the files that you are going to upload later											
		#Sample names must be unique even with the ones that are in database											
		sample name 	fastq1	fastq2	data set	vaccine status	week	onset date	collection date	lab reception date	latitude	longitude	other fields from here
		
		return a instance of a class with all lines
		if something is wrong, add the errors 
		in: dont_detect_sample_names sometimes is necessary to read the file witout
		"""
		
		### clean data first
		self.clean_data()
		
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
			dialect = sniffer.sniff(line)
			delimiter = dialect.delimiter
			if (delimiter in dict_delimiter): dict_delimiter[delimiter] += 1
			else: dict_delimiter[delimiter] = 1
			if (n_lines > 20): break
			n_lines += 1
			
		### get delimiter
		delimiter = ','
		for key in dict_delimiter:
			delimiter = key
			f.seek(0)
			n_lines = 0
			b_found = False
			reader = csv.reader(f, delimiter=delimiter)
			for row in reader:
				(count_column, n_columns_ok) = (0, 0)
				for col in row:
					if (count_column >= len(self.vect_header) or col.replace(' ', '').lower() != self.vect_header[count_column].replace(' ', '').lower()): break
					n_columns_ok += 1
					count_column += 1
				if (n_columns_ok == len(self.vect_header)):
					b_found = True
					break

				## after XX lines stop find				
				if (n_lines > 40): break
				n_lines += 1

			## header was founded				
			if (b_found): break

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
					for i in range(0, len(row)):	## test all header names
						if (len(row[i].strip()) == 0): continue
						if (row[i].strip() in self.dict_other_fields_repeated):
							self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("Column name '{}' is repeated in the header. Line: {} Column: {}".\
										format(row[i].strip(), count_row, i+1))))
						else:
							self.dict_other_fields_repeated[row[i].strip()] = 1

			elif (header != None): ## line with data
				self.process_row(row, count_row, header, user, read_state)
			count_row += 1
		if (header == None):
			self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("Header not found in the file. Please, check the names in the header, must be equal and have the same order of the template file.")))
		elif (self.number_samples == 0 and not self.errors.has_errors()):
			self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("There's no samples to process.")))

	def process_row(self, row, count_row, header, user, read_state):
		"""
		in: read_state -> STATE_READ_only_detect_errors; STATE_READ_detect_errors_and_not_test_samples_name; STATE_READ_all
		in: only_detect_errors -> only to detect errors, not to add samples
		process header
		"""
		### check sample name
		n_errors = len(self.errors.get_vect_results())
		if (read_state != ParseInFiles.STATE_READ_dont_detect_errors):
			if len(row) > 0 and len(row[0].strip()) > 0:
				sample_name = row[0].strip()
				
				## test clean name
				result_filer_sample_name = re.sub('[^A-Za-z0-9_]+', '', sample_name)
				if (len(result_filer_sample_name) != len(sample_name)):
					self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("Sample name '{}' only letters, numbers and underscores are allowed. Line: {} Column: {}".format(sample_name, count_row, 1))))
				else:
					try:
						sample = Sample.objects.get(name__iexact=sample_name, owner=user, is_deleted=False)
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
# 					number_files = UploadFiles.objects.filter(file_name__iexact=fastq1.strip(), owner=user,\
# 										is_processed=False, is_deleted=False, type_file__name=TypeFile.TYPE_FILE_fastq_gz).count()
# 					if (number_files > 0):
# 						self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("File '{}' is repeated in the database and it's not processed yet. Line: {} Column: {}".\
# 																	format(fastq1.strip(), count_row, 2))))
# 					else:
					self.dict_file_names[fastq1.strip()] = 1
			else:
				self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("There's no fastq1 file name. Line: {} Column: {}".format(count_row, 2))))
	
			### check if fastq2 file as gz
			if len(row) > 2 and len(row[2].strip()) > 0:
				fastq2 = row[2]
				if (not fastq2.endswith(FileExtensions.FILE_GZ)):
					self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("File '' not ends with '{}'. Fastq gzip compressed file is needed. Line: {} Column: {}".\
											format(fastq2, FileExtensions.FILE_GZ, count_row, 3))))
				elif (fastq2.strip() in self.dict_file_names):
					self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("File '{}' is repeated in the samples file. Line: {} Column: {}".\
											format(fastq2.strip(), count_row, 3))))
				else:
# 					number_files = UploadFiles.objects.filter(file_name__iexact=fastq2.strip(), owner=user,\
# 										is_processed=False, is_deleted=False, type_file__name=TypeFile.TYPE_FILE_fastq_gz).count()
# 					if (number_files > 0):
# 						self.errors.add_single_result(SingleResult(SingleResult.ERROR, _("File '{}' is repeated in the database and it's not processed yet. Line: {} Column: {}".\
# 																	format(fastq2.strip(), count_row, 3))))
# 					else:
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

		self.number_samples += 1
		### there's no errors, process this sample
		if (n_errors == len(self.errors.get_vect_results()) and\
				(read_state == ParseInFiles.STATE_READ_all or\
				read_state == ParseInFiles.STATE_READ_dont_detect_errors)):
			
			### add samples
			sample = Sample()
			sample.id = count_row
			sample.owner = user
			sample.name = row[0]
			sample.candidate_file_name_1 = row[1].strip()
			sample.candidate_file_name_2 = ""
			if (len(row[2].strip()) > 0): sample.candidate_file_name_2 = row[2].strip()
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


	@transaction.atomic
	def create_samples(self, upload_files, user):
		"""
		create samples in database
		run this after the file sample_files is uploaded
		Only run if there's no errors
		"""
		## you can do anything with errors
		if (self.errors.get_len_vect_results() > 0): return
		
		for vect_sample in self.vect_samples:
			
			try:
				sample = Sample.objects.get(name__iexact=vect_sample[0].name, owner=user, is_deleted=False)
				## if exist don't add
				upload_files.samples.add(sample)
				continue 	## already exist
			except Sample.DoesNotExist as e:
				pass
			
			
			## data set
			data_set = None
			if (vect_sample[0].data_set != None):
				try:
					data_set = DataSet.objects.get(name__iexact=vect_sample[0].data_set.name, owner=user)
				except DataSet.DoesNotExist as e:
					data_set = DataSet()
					data_set.name = vect_sample[0].data_set.name
					data_set.owner = user
					data_set.save()
			
			vaccine_status = None
			if (vect_sample[0].vaccine_status != None):
				try:
					vaccine_status = VaccineStatus.objects.get(name__iexact=vect_sample[0].vaccine_status.name, owner=user)
				except VaccineStatus.DoesNotExist as e:
					vaccine_status = VaccineStatus()
					vaccine_status.name = vect_sample[0].vaccine_status.name
					vaccine_status.owner = user
					vaccine_status.save()
					
			## first set sample
			sample = Sample()
			sample.name = vect_sample[0].name
			sample.owner = vect_sample[0].owner
			sample.candidate_file_name_1 = vect_sample[0].candidate_file_name_1
			sample.candidate_file_name_2 = vect_sample[0].candidate_file_name_2
			sample.data_set = data_set
			sample.vaccine_status = vaccine_status
			sample.week = vect_sample[0].week
			sample.date_of_onset = vect_sample[0].date_of_onset
			sample.date_of_collection = vect_sample[0].date_of_collection
			sample.date_of_receipt_lab = vect_sample[0].date_of_receipt_lab
			sample.geo_local = vect_sample[0].geo_local
			sample.save()

			### now save the tag names
			for tag_samples_temp in vect_sample[1]:
				try:
					tag_name = TagName.objects.get(name__iexact=tag_samples_temp.tag_name.name, owner=user)
				except TagName.DoesNotExist as e:
					tag_name = TagName()
					tag_name.owner = user
					tag_name.name = tag_samples_temp.tag_name.name
					tag_name.is_meta_data = False
					tag_name.save()
					
				tag_names = TagNames(tag_name=tag_name, sample=sample, value=tag_samples_temp.value)
				tag_names.save()

			upload_files.samples.add(sample)
			
		## save uplaod files
		upload_files.is_valid = True
		upload_files.is_processed = False			## True when all samples are set
		upload_files.number_files_to_process = len(self.vect_samples)
		upload_files.number_files_processed = 0		## has the number of files linked 
		upload_files.save()

	def has_samples_files_to_process(self, user):
		"""
		Test if there are any files to process
		"""
		upload_files = self.get_upload_samples_file(user)
		return not (upload_files == None)

	def get_upload_samples_file(self, user):
		"""
		test if there are any samples files to process
		"""
		upload_files = UploadFiles.objects.filter(owner=user, is_valid=True, is_deleted=False, is_processed=False,\
						type_file__name=TypeFile.TYPE_FILE_sample_file)
		if (upload_files.count() == 0): return None
		return upload_files[0]

	def link_files(self, user, b_testing = False):
		
		### make it running 
		process_controler = ProcessControler()
		process_SGE = ProcessSGE()
		process_SGE.set_process_controler(user, process_controler.get_name_link_files_user(user), ProcessControler.FLAG_RUNNING)
		
		## need to add a delay for the test in command line
		if (settings.RUN_TEST_IN_COMMAND_LINE): time.sleep(4)
		
		try:
			self._link_files(user, b_testing)
		except:
			## finished with error
			process_SGE.set_process_controler(user, process_controler.get_name_link_files_user(user), ProcessControler.FLAG_ERROR)
			return
		
		### finished
		process_SGE.set_process_controler(user, process_controler.get_name_link_files_user(user), ProcessControler.FLAG_FINISHED)
		
	@transaction.atomic
	def _link_files(self, user, b_testing = False):
		"""
		after a fastq.gz is uploaded you can run this to link the fastq files to samples.
		"""
		
		software = Software()
		utils = Utils()
		upload_files = self.get_upload_samples_file(user)
		if (upload_files == None): return	## there's no files to match

		n_files_processed = 0
		for sample in upload_files.samples.all():
			if (sample.has_files):
				n_files_processed += 1
				continue
			
			### if not files
			upload_files_1 = None
			upload_files_2 = None
			if (len(sample.candidate_file_name_1) > 0):
				try:
					upload_files_1 = UploadFiles.objects.get(file_name=sample.candidate_file_name_1, owner=user, is_processed=False,\
										is_deleted=False, is_valid=True, type_file__name=TypeFile.TYPE_FILE_fastq_gz)
					path_to_file = os.path.join(getattr(settings, "MEDIA_ROOT", None), upload_files_1.path_name.name)
					if (not os.path.exists(path_to_file) or not utils.is_fastq_gz(path_to_file)): upload_files_1 = None
				except UploadFiles.DoesNotExist as e:
					pass
			if (upload_files_1 != None and len(sample.candidate_file_name_2) > 0):
				try:
					upload_files_2 = UploadFiles.objects.get(file_name=sample.candidate_file_name_2, owner=user, is_processed=False,\
										is_deleted=False, is_valid=True, type_file__name=TypeFile.TYPE_FILE_fastq_gz)
					path_to_file = os.path.join(getattr(settings, "MEDIA_ROOT", None), upload_files_2.path_name.name)
					if (not os.path.exists(path_to_file) or not utils.is_fastq_gz(path_to_file)): upload_files_2 = None
				except UploadFiles.DoesNotExist as e:
					pass
			
			## can set the data
			if ((len(sample.candidate_file_name_2) > 0 and upload_files_1 != None and upload_files_2 != None) or\
				(len(sample.candidate_file_name_2) == 0 and upload_files_1 == None)):
				
				## link the files to the right place
				sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), utils.get_path_to_fastq_file(user.id, sample.id), sample.candidate_file_name_1)
				utils.link_file(os.path.join(getattr(settings, "MEDIA_ROOT", None), upload_files_1.path_name.name), sz_file_to)
				sample.path_name_1.name = os.path.join(utils.get_path_to_fastq_file(user.id, sample.id), sample.candidate_file_name_1)
				sample.file_name_1 = sample.candidate_file_name_1
				sample.is_valid_1 = True
				
				upload_files_1.is_processed = True
				upload_files_1.upload_file = upload_files
				upload_files_1.number_files_processed = 1
				upload_files_1.samples.add(sample)
				upload_files_1.attached_date = datetime.now()
				upload_files_1.save()
				
				if (upload_files_2 != None):
					sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), utils.get_path_to_fastq_file(user.id, sample.id), sample.candidate_file_name_2)
					utils.link_file(os.path.join(getattr(settings, "MEDIA_ROOT", None), upload_files_2.path_name.name), sz_file_to)
					sample.path_name_2.name = os.path.join(utils.get_path_to_fastq_file(user.id, sample.id), sample.candidate_file_name_2)
					sample.file_name_2 = sample.candidate_file_name_2
					sample.is_valid_2 = True
					
					upload_files_2.is_processed = True
					upload_files_2.upload_file = upload_files
					upload_files_2.number_files_processed = 1
					upload_files_2.samples.add(sample)
					upload_files_2.attached_date = datetime.now()
					upload_files_2.save()
					
				### sample file
				sample.has_files = True
				sample.save()
				n_files_processed += 1
				
				if (not b_testing):
					### send signal to process the type and sub type
					### create a task to perform the analysis of fastq and trimmomatic
					try:
						## here can be direct because came from a djangoq
						if (settings.RUN_SGE):
							process_SGE = ProcessSGE()
							if (settings.RUN_SGE_INTO_DJANGOQ):
								taskID = async(process_SGE.set_run_trimmomatic_species, sample, user)
							else:
								taskID = process_SGE.set_run_trimmomatic_species(sample, user)
						else:
							taskID = async(software.run_fastq_and_trimmomatic_and_identify_species, sample, user)
						
						### 
						manageDatabase = ManageDatabase()
						manageDatabase.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue, taskID)
					except:
						pass
					
		### set the files already processed
		upload_files.number_files_processed = n_files_processed
		if (upload_files.number_files_processed == upload_files.number_files_to_process): upload_files.is_processed = True
		upload_files.save()
		
			
	
class UploadFilesByDjangoQ(object):
	
	def __init__(self):
		pass
	
	def read_sample_file(self, user, upload_files, b_testing = False):
		"""
		read samples csv file, and link files if they exists		
		"""
		### make it running 
		process_controler = ProcessControler()
		process_SGE = ProcessSGE()
		process_SGE.set_process_controler(user, process_controler.get_name_upload_files(upload_files), ProcessControler.FLAG_RUNNING)
		
		## need to add a delay for the test in command line
		if (settings.RUN_TEST_IN_COMMAND_LINE): time.sleep(4)
		
		try:
			parse_in_files = ParseInFiles()
			b_test_char_encoding = True
			parse_in_files.parse_sample_files(upload_files.get_path_to_file(TypePath.MEDIA_ROOT), user, b_test_char_encoding, ParseInFiles.STATE_READ_all)
			if (parse_in_files.get_errors().has_errors()): return False
			
			parse_in_files.create_samples(upload_files, user)
			parse_in_files._link_files(user, b_testing)	### without passing the messages 
		except:
			## finished with error
			process_SGE.set_process_controler(user, process_controler.get_name_upload_files(upload_files), ProcessControler.FLAG_ERROR)
			return
		
		### finished
		process_SGE.set_process_controler(user, process_controler.get_name_upload_files(upload_files), ProcessControler.FLAG_FINISHED)
		return True

