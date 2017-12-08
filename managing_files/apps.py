from django.apps import AppConfig
# 
from constants.constants import Constants

class ManagingFilesConfig(AppConfig):
	name = 'managing_files'
	verbose_name = "Managing Files"
	
	def ready(self):
		
		## create a default user
		self.create_default_user()
			
		#### Now upload the 
		self.upload_default_files()

	def create_default_user(self):
		"""
		create a default user to link the default references...
		"""
		
		from django.contrib.auth.models import User
		from .models import DataSet
		try:
			User.objects.get(username=Constants.DEFAULT_USER)
			### great, the default user exist 
		except User.DoesNotExist:
			
			### need to create it
			user = User()
			user.username = Constants.DEFAULT_USER
			user.password = Constants.DEFAULT_USER_PASS
			user.first_name = Constants.DEFAULT_USER
			user.is_active = False
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

