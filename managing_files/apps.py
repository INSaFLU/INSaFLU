from django.apps import AppConfig
# 
from utils.Constants import Constants

class ManagingFilesConfig(AppConfig):
	name = 'managing_files'
	verbose_name = "Managing Files"
	
# 	def ready(self):
# 		from django.contrib.auth.models import User
		## create a default user
	#	self.create_default_user()
			
		#### Now upload the 
	#	self.upload_default_files()
		

# 	def create_default_user(self):
# 		"""
# 		create a default user to link the default references...
# 		"""
# 		try:
# 			User.objects.get(name=self.utils.DEFAULT_USER)
# 			### great, the default user exist 
# 		except User.DoesNotExist:
# 			
# 			### need to create it
# 			user = User()
# 			user.username = self.utils.DEFAULT_USER
# 			user.password = self.utils.DEFAULT_USER_PASS
# 			user.first_name = self.utils.DEFAULT_USER
# 			user.is_active = False
# 			user.is_staff = False
# 			user.is_superuser = False
# 			user.save()
# 	
# 	def upload_default_files(self):
# 		"""
# 		Upload default files
# 		"""
# 		pass