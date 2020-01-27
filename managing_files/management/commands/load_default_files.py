'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from constants.constants import Constants
from django.db import transaction
import logging

class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Reload default references."
	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	def upload_default_files(self):
		"""
		Upload default files
		"""
		## only runs once, wen start ans test if the file was uploaded with virus hypothesis
		from manage_virus.uploadFiles import UploadFiles
		from utils.software import Software
		uploadFiles = UploadFiles()
		## get version and pah
		b_test = False
		(version, path) = uploadFiles.get_file_to_upload(b_test)
		
		## uplaod
		uploadFile = uploadFiles.upload_file(version, path)

		# create the abricate database
		if (uploadFile != None):
			software= Software()
			if (not software.is_exist_database_abricate(uploadFile.abricate_name)):
				software.create_database_abricate(uploadFile.abricate_name, uploadFile.path)
	
	@transaction.atomic
	def upload_default_references(self):
		"""
		upload default files for reference
		"""
		from manage_virus.uploadFiles import UploadFiles
		from django.contrib.auth.models import User
		
		try:
			User.objects.get(username=Constants.DEFAULT_USER)
			### great, the default user exist 
		except User.DoesNotExist:
			self.stdout.write("Upload References failed because there's not default user.\n" +
							"Please, start the application first.")
			return
		
		uploadFiles = UploadFiles()
		b_test = False
		uploadFiles.upload_default_references(User.objects.get(username=Constants.DEFAULT_USER), b_test) 
			
	# A command must define handle()
	def handle(self, *args, **options):
		
		#### set default fields
		self.stdout.write("Upload abricate files")
		self.upload_default_files()
		
		#### set default references
		self.stdout.write("Upload References")
		self.upload_default_references()
		self.stdout.write("End")
		