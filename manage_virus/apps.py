from django.apps import AppConfig

class ManageVirusConfig(AppConfig):
	name = 'manage_virus'
	verbose_name = "Managing Virus Config"
	
	def ready(self):
		pass
		## only runs once, wen start ans test if the file was uploaded with virus hypothesis
#		from .uploadFiles import UploadFiles
#		from utils.Software import Software
# 		uploadFiles = UploadFiles()
# 		## get version and pah
# 		b_test = False
# 		(version, path) = uploadFiles.get_file_to_upload(b_test)
# 		
# 		## uplaod
#		uploadFile = uploadFiles.upload_file(version, path)

#		# create the abricate database
# 		if (uploadFile != None):
# 			software= Software()
# 			if (not software.is_exist_database_abricate(uploadFile.abricate_name)):
# 				software.create_database_abricate(uploadFile.abricate_name, uploadFile.path)


