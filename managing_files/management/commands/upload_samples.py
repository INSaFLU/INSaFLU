'''
Created on Aug 1, 2022

@author: daniel.sobral
'''
from django.core.management import BaseCommand
from django.contrib.auth.models import User
from django.conf import settings
from managing_files.models import MetaKey, UploadFiles
from constants.constants import Constants, TypeFile
from utils.parse_in_files import ParseInFiles
from utils.utils import Utils
import os, ntpath
import logging


class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Upload metadata file and fastq files to create samples."

	# logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")

	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)

	def add_arguments(self, parser):
		parser.add_argument('--metadata_file', nargs='?', type=str,
							required=True, help='Metadata file to upload')
		parser.add_argument('--user_login', nargs='?', type=str, required=True,
							help='User login to whom the samples will be associated')

	# A command must define handle()
	def handle(self, *args, **options):

		utils = Utils()

		metadata_file = options['metadata_file']
		account = options['user_login']
		self.stdout.write("Starting to upload: " + str(metadata_file))
		# This logger only makes sense if this script is ever integrated as part of the website functionality...
		# self.logger_production.info("Starting to upload: " + str(metadata_file))
		# self.logger_debug.info("Starting to upload: " + str(metadata_file))
		self.stdout.write("Warning: fastq files bigger than {} will be randomly reduced to this size".format(settings.MAX_FASTQ_FILE_UPLOAD))

		try:
			
			user = User.objects.get(username=account)

			self.stdout.write("Upload to user: " + str(user))

			metadata_full_path = os.path.join(getattr(settings, "MEDIA_ROOT", None), Constants.DIR_PROCESSED_FILES_UPLOADS, metadata_file)

			if(not os.path.exists(metadata_full_path)):
				self.stdout.write("Metadata file {} could not be found".format(metadata_full_path))
				return False

			# Process the metadata file to check if everything is ok
			parse_in_files = ParseInFiles()
			b_test_char_encoding = False
			parse_in_files.parse_sample_files(metadata_full_path, user, b_test_char_encoding, ParseInFiles.STATE_READ_all)

			if (parse_in_files.get_errors().has_errors()):
				self.stdout.write("Errors found processing the metadata file {}".format(metadata_file))
				self.stdout.write(str(parse_in_files.get_errors()))
				# self.logger_debug.error("Errors found processing the metadata table")
				# self.logger_debug.erro(str(parse_in_files.get_errors()))
				return False

			# Check if referenced fastQ files exist before creating the official upload of the sample metadata file
			fastq_files_to_upload = []
			missing_fastqs = False
			for sample in parse_in_files.get_vect_samples():
				fastq1 = sample[0].candidate_file_name_1
				
				fastq_full_path = os.path.join(getattr(settings, "MEDIA_ROOT", None), Constants.DIR_PROCESSED_FILES_UPLOADS, fastq1)
				if(not os.path.exists(fastq_full_path)): 
					self.stdout.write("Fastq file {} could not be found".format(fastq_full_path))
					missing_fastqs = True
				fastq_files_to_upload.append(fastq_full_path)

				fastq2 = ""
				# Carefull this may be sensitive to spaces in the table
				if(sample[0].candidate_file_name_2.strip() != ""):
					fastq2 = sample[0].candidate_file_name_2
					fastq_full_path = os.path.join(getattr(settings, "MEDIA_ROOT", None), Constants.DIR_PROCESSED_FILES_UPLOADS, fastq2)
					fastq_files_to_upload.append(fastq_full_path)
					if(not os.path.exists(fastq_full_path)): 
						self.stdout.write("Fastq file {} could not be found".format(fastq_full_path))
						missing_fastqs = True
				
				self.stdout.write("sample {} file(s) to be processed: {} {} ".format(sample[0].name, fastq1, fastq2))

			if(missing_fastqs):
				self.stdout.write("Fastq files are missing, cannot continue")
				return False

			self.stdout.write(" {} samples are going to be processed".format(len(parse_in_files.get_vect_samples())))
			
			# Make an atomic transaction for the whole process? 
			# so if it fails in the middle, the database will revert to the beginning...
			# Though potentially large files being copied will still be spread everywhere with no tracing...

			# Add the metadata file as an upload
			# may not be needed, but for consistency with website we do it
			sample_file_upload_files = UploadFiles()
			sample_file_upload_files.file_name = metadata_file
			sample_file_upload_files.is_valid = True
			sample_file_upload_files.is_processed = False
			sample_file_upload_files.is_deleted = False
			sample_file_upload_files.number_errors = 0
			sample_file_upload_files.number_files_processed = 0
			sample_file_upload_files.number_files_to_process = len(parse_in_files.get_vect_samples())
			try:
				type_file = MetaKey.objects.get(name=TypeFile.TYPE_FILE_sample_file)
			except MetaKey.DoesNotExist:
				type_file = MetaKey()
				type_file.name = TypeFile.TYPE_FILE_sample_file
				type_file.save()

			sample_file_upload_files.type_file = type_file
			sample_file_upload_files.owner = user

			sample_file_upload_files.description = ""

			# move the files to the right place
			sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), utils.get_path_upload_file(user.id,
									TypeFile.TYPE_FILE_sample_file), metadata_file)

			# get unique file name, as the user can upload files with same name...
			sz_file_to, path_added = utils.get_unique_file(sz_file_to)

			# Add this back in the end... to "consume" the file
			# utils.move_file(metadata_full_path, sz_file_to)
			utils.copy_file(metadata_full_path, sz_file_to)
			
			if path_added is None:
				sample_file_upload_files.path_name.name = os.path.join(self.utils.get_path_upload_file(self.request.user.id,\
									TypeFile.TYPE_FILE_sample_file), ntpath.basename(sz_file_to))
			else:
				sample_file_upload_files.path_name.name = os.path.join(self.utils.get_path_upload_file(self.request.user.id,\
									TypeFile.TYPE_FILE_sample_file), path_added, ntpath.basename(sz_file_to))
					
			self.stdout.write("{} file was processed".format(sample_file_upload_files.path_name.name))

			# Do not forget to add 
			uploads_to_save = [sample_file_upload_files]
			# sample_file_upload_files.save()

			# upload_files_by_djangoq = UploadFilesByDjangoQ()

			# Upload the sample fastq files

			for fastq_to_upload in fastq_files_to_upload:
				self.stdout.write("Fastq file to upload: {}".format(fastq_to_upload))
				
				fastq_upload_files = UploadFiles()
				fastq_upload_files.file_name = utils.clean_name(os.path.basename(fastq_to_upload))

				# move the files to the right place
				sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), utils.get_path_upload_file(user.id, TypeFile.TYPE_FILE_fastq_gz),
										fastq_upload_files.file_name)
				sz_file_to, path_added = utils.get_unique_file(sz_file_to)		## get unique file name, user can upload files with same name...
				#	utils.move_file(temp_file, sz_file_to)
				utils.copy_file(fastq_to_upload, sz_file_to)

				# test if file exists (may fail due to full disk or other error)
				if (not os.path.exists(sz_file_to) and os.path.getsize(sz_file_to) > 10):
					self.stdout.write(" Error copying file {} file to {}".format(fastq_to_upload, sz_file_to))					
					# If we do a return here then we need an atomic transaction or a way to prevent inconsistencies...
					return False
				
				if path_added is None:
					fastq_upload_files.path_name.name = os.path.join(self.utils.get_path_upload_file(self.request.user.id,\
										TypeFile.TYPE_FILE_fastq_gz), fastq_upload_files.file_name)
				else:
					fastq_upload_files.path_name.name = os.path.join(self.utils.get_path_upload_file(self.request.user.id,\
										TypeFile.TYPE_FILE_fastq_gz), path_added, fastq_upload_files.file_name)
					
				try:
					type_file = MetaKey.objects.get(name=TypeFile.TYPE_FILE_fastq_gz)
				except MetaKey.DoesNotExist:
					type_file = MetaKey()
					type_file.name = TypeFile.TYPE_FILE_fastq_gz
					type_file.save()
				
				fastq_upload_files.is_valid = True
				fastq_upload_files.is_processed = False			## True when all samples are set
				fastq_upload_files.owner = user
				fastq_upload_files.type_file = type_file
				fastq_upload_files.number_files_to_process = 1
				fastq_upload_files.number_files_processed = 0
				fastq_upload_files.description = ""

				uploads_to_save.append(fastq_upload_files)
				#fastq_upload_files.save()
			
			# Save all the UploadFiles
			for upload_file_to_save in uploads_to_save:
				upload_file_to_save.save()

			# Add the sample files and link the fastq files to the samples...					 
			#upload_files_by_djangoq = UploadFilesByDjangoQ()
			#upload_files_by_djangoq.read_sample_file(user, sample_file_upload_files, settings.RUN_TEST_IN_COMMAND_LINE)
			parse_in_files.create_samples(sample_file_upload_files, user)
			parse_in_files.link_files(user, False)

			self.stdout.write("End")

		except User.DoesNotExist as e:
			self.stdout.write("Error: User '{}' does not exist.".format(account))
			#self.logger_production.error(
			#	"Error: User '{}' does not exist.".format(account))
			#self.logger_debug.error(
			#	"Error: User '{}' does not exist.".format(account))
		except Exception as e:
			self.stdout.write("Error: {}.".format(e))
			#self.logger_production.error(
			#	"Error: {}.".format(e))
			#self.logger_debug.error(
			#	"Error: {}.".format(e))
