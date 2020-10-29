'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from constants.constants import Constants, TypeFile
from django.contrib.auth.models import User
from managing_files.models import UploadFiles, MetaKey
from datetime import datetime
import logging

class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Unlock upload files for all accounts."
	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)


	def unlock_upload_file(self, upload_file):
		"""
		unlock upload file
		"""
		for sample in upload_file.samples.all():
			if (not sample.is_ready_for_projects and not sample.is_deleted):
				self.stdout.write("\t\tSample name/id: {}/{} is deleted...".format(sample.name, sample.id))
				sample.is_deleted = True
				sample.date_deleted = datetime.now()
				sample.save()
		upload_file.is_processed = True	
		upload_file.save()


	# A command must define handle()
	def handle(self, *args, **options):
		
		#### set default fields
		self.stdout.write("Unlock upload multiple files:")
		list_users = User.objects.all().order_by('date_joined');
		
		### try to find sample_file imported not processed yet
		try:
			metaKey = MetaKey.objects.get(name=TypeFile.TYPE_FILE_sample_file)
		except MetaKey.DoesNotExist:
			self.stdout.write("There's no files to unlock")
			return
		
		number_of_changes = 0
		for user in list_users:
			### default users, don't do it...
			if user.username in [Constants.DEFAULT_USER, Constants.USER_ANONYMOUS]: continue
			self.stdout.write("####################\nProcessing user: {}".format(user.username))
			
			lst_files = UploadFiles.objects.all().filter(is_deleted=False, is_processed=False, owner=user, type_file=metaKey)
			for uploadfile in lst_files:
				
				## to remove
				if uploadfile.number_files_processed == uploadfile.number_files_to_process: continue
				self.stdout.write("File to locked: '{}'   date: '{}'".format(uploadfile.file_name, uploadfile.creation_date))
				
				b_unlock = False
				while 1:
					input_anwser = input("\tDo you want to unlock [yes/no] (no): ")
					input_anwser = input_anwser.lower().strip()
					if len(input_anwser) == 0:
						self.stdout.write("\t\tAssumed 'no' answer...")
						break
					if input_anwser == 'yes':
						b_unlock = True
						break
					elif input_anwser == 'no': break
					else: self.stdout.write("\t\tOnly accept 'yes/no' answer, try again.")
				
				### going to unlock
				if b_unlock:
					self.unlock_upload_file(uploadfile)
					number_of_changes += 1
			
			self.stdout.write("\n")
			
		self.stdout.write("Number of unlocked files: {}".format(number_of_changes))




