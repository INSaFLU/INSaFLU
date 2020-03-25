'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from constants.constants import Constants
from constants.constants_mixed_infection import ConstantsMixedInfection
from managing_files.models import MixedInfectionsTag
from manage_virus.uploadFiles import UploadFiles
from django.contrib.auth.models import User
from managing_files.models import DataSet
from extend_user.models	import Profile
from utils.software import Software
from django.db import transaction
from django.conf import settings
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
		## only runs once, when it start and test if the file was uploaded with virus hypothesis
		uploadFiles = UploadFiles()
		## get version and path
		b_test = False
		(version, path) = uploadFiles.get_file_to_upload(b_test)
		
		## upload
		uploadFile = uploadFiles.upload_file(version, path)

		# create the Abricate database
		if (not uploadFile is None):
			software= Software()
			if (not software.is_exist_database_abricate(uploadFile.abricate_name)):
				self.stdout.write("Create new Abricate database: {}".format(uploadFile.abricate_name))
				software.create_database_abricate(uploadFile.abricate_name, uploadFile.path)
	
	@transaction.atomic
	def upload_default_references(self):
		"""
		upload default files for reference
		"""

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
	
	
	def default_database_fields(self):
		"""
		set default fields in database
		"""
		### MixedInfectionsTag
		constants_mixed_infection = ConstantsMixedInfection()
		for tag in constants_mixed_infection.vect_upload_to_database:
			try:
				mixed_infections_tag = MixedInfectionsTag.objects.get(name=tag)
			except MixedInfectionsTag.DoesNotExist as e:
				mixed_infections_tag = MixedInfectionsTag()
				mixed_infections_tag.name = tag
				mixed_infections_tag.save()

	def create_default_user(self):
		"""
		create a default user to link the default references...
		"""
		self.__create_account(Constants.DEFAULT_USER, Constants.DEFAULT_USER_PASS, settings.DEFAULT_USER_EMAIL, False)
		self.__create_account(Constants.USER_ANONYMOUS, Constants.USER_ANONYMOUS_PASS, settings.USER_ANONYMOUS_EMAIL, True)


	def __create_account(self, user_name, password, email, b_active):
		"""
		try to create a default accounts
		"""
		try:
			user = User.objects.get(username=user_name)
			### great, the default user exist
		except User.DoesNotExist:
			
			### need to create it
			self.stdout.write("Add user: {}".format(user_name))
			user = User()
			user.username = user_name
			user.set_password(password)
			user.first_name = user_name
			user.email = email
			user.is_active = b_active
			user.is_staff = False
			user.is_superuser = False
			user.save()
	
		### create generic dataset
		for user in User.objects.all():
			result = DataSet.objects.filter(owner__id=user.id)
			if (len(result) == 0):
				### need to create it
				dataSet = DataSet()
				dataSet.name = Constants.DATA_SET_GENERIC
				dataSet.owner = user
				dataSet.save()
		
		### for security reasons
		### set true for anonymous user always
		if (user_name == Constants.USER_ANONYMOUS):
			try:
				profile = Profile.objects.get(user__username=user_name)
				profile.only_view_project = True
				profile.save()
			except Profile.DoesNotExist:
				pass
			
	# A command must define handle()
	def handle(self, *args, **options):
		
		#### set default fields
		self.stdout.write("Upload abricate files")
		self.upload_default_files()
		
		self.stdout.write("Define default database fields")
		self.default_database_fields()
		
		self.stdout.write("Set default users...")
		self.create_default_user()
		
		#### set default references
		self.stdout.write("Upload References")
		self.upload_default_references()
		self.stdout.write("End")
		