from django.apps import AppConfig
# 
from constants.constants import Constants
from constants.constants_mixed_infection import ConstantsMixedInfection
from django.db import transaction


class ManagingFilesConfig(AppConfig):
	name = 'managing_files'
	verbose_name = "Managing Files"
	
	def ready(self):
		
		## create a default user
		self.create_default_user()
			
		#### Now upload the 
		self.upload_default_files()
		
		#### set default fields
		self.default_database_fields()
		
		#### set default references
		self.upload_default_references()
		
		pass

	def create_default_user(self):
		"""
		create a default user to link the default references...
		"""
		self.__create_account(Constants.DEFAULT_USER, Constants.DEFAULT_USER_PASS, Constants.DEFAULT_USER_EMAIL, False)
		self.__create_account(Constants.USER_ANONYMOUS, Constants.USER_ANONYMOUS_PASS, Constants.USER_ANONYMOUS_EMAIL, True)


	def __create_account(self, user_name, password, email, b_active):
		"""
		try to create a default accounts
		"""
		from django.contrib.auth.models import User
		from managing_files.models import DataSet
		from extend_user.models	import Profile
		try:
			user = User.objects.get(username=user_name)
			### great, the default user exist
		except User.DoesNotExist:
			
			### need to create it
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
				profile = Profile.objects.get(user=user)
				profile.only_view_project = True
				profile.save()
			except Profile.DoesNotExist:
				pass


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


	def default_database_fields(self):
		"""
		set default fields in database
		"""
		
		### MixedInfectionsTag
		from managing_files.models import MixedInfectionsTag
		constants_mixed_infection = ConstantsMixedInfection()
		for tag in constants_mixed_infection.vect_upload_to_database:
			try:
				mixed_infections_tag = MixedInfectionsTag.objects.get(name=tag)
			except MixedInfectionsTag.DoesNotExist as e:
				mixed_infections_tag = MixedInfectionsTag()
				mixed_infections_tag.name = tag
				mixed_infections_tag.save()
	
	
	@transaction.atomic
	def upload_default_references(self):
		"""
		upload default files for reference
		"""
		from manage_virus.uploadFiles import UploadFiles
		from django.contrib.auth.models import User
		
		uploadFiles = UploadFiles()
		self.create_default_user()
		b_test = False
		uploadFiles.upload_default_references(User.objects.get(username=Constants.DEFAULT_USER), b_test) 
		
