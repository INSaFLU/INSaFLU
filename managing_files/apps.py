from django.apps import AppConfig
# from django.contrib.auth.models import User
from constants.Constants import Constants

class ManagingFilesConfig(AppConfig):
	name = 'managing_files'
#	verbose_name = "Managing Files"
	
# 	constants = Constants()
# 	def ready(self):
# 		pass
		## create a default user
	#	self.create_default_user()
			
		#### Now upload the 
	#	self.upload_default_files()
		

# 	def create_default_user(self):
# 		"""
# 		create a default user to link the default references...
# 		"""
# 		try:
# 			User.objects.get(name=self.constants.DEFAULT_USER)
# 			### great, the default user exist 
# 		except User.DoesNotExist:
# 			
# 			### need to create it
# 			user = User()
# 			user.username = self.constants.DEFAULT_USER
# 			user.password = self.constants.DEFAULT_USER_PASS
# 			user.first_name = self.constants.DEFAULT_USER
# 			user.is_active = False
# 			user.is_staff = False
# 			user.is_superuser = False
# 			user.save()
# 	
# 	def upload_default_files(self):
# 		"""
# 		Upload default files
# 		"""