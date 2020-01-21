'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from utils.parse_in_files import UpdateMetadataFileByDjangoQ
from managing_files.models import UploadFiles
from django.contrib.auth.models import User
import logging

class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Update samples in database."
	
	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")

	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	## https://docs.djangoproject.com/en/dev/howto/custom-management-commands/
	def add_arguments(self, parser):
		parser.add_argument('--upload_files_id', type=int, help='Upload files ID')
		parser.add_argument('--user_id', nargs='?', type=int, help='User id to join to this process')
	
	# A command must define handle()
	def handle(self, *args, **options):
		
		upload_files_by_djangoq = UpdateMetadataFileByDjangoQ()
		upload_files_id = options['upload_files_id']
		user_id = options['user_id']
		self.stdout.write("Starting for upload_files_id: " + str(upload_files_id))
		self.logger_production.info("Starting for upload_files_id: " + str(upload_files_id))
		self.logger_debug.info("Starting for upload_files_id: " + str(upload_files_id))
		try:
			upload_files = UploadFiles.objects.get(pk=upload_files_id)
			if (user_id == None): user = upload_files.owner
			else: user = User.objects.get(pk=user_id)
			upload_files_by_djangoq.update_sample_file(user, upload_files)
			self.stdout.write("End")
		except UploadFiles.DoesNotExist as e:
			self.stdout.write("Error: UploadFiles id '{}' does not exist.".format(upload_files_id))
			self.logger_production.info("Error: UploadFiles id '{}' does not exist.".format(upload_files_id))
			self.logger_debug.info("Error: UploadFiles id '{}' does not exist.".format(upload_files_id))
		except User.DoesNotExist as e:
			self.stdout.write("Error: User id '{}' does not exist.".format(user_id))
			self.logger_production.info("Error: User id '{}' does not exist.".format(user_id))
			self.logger_debug.info("Error: User id '{}' does not exist.".format(user_id))
		