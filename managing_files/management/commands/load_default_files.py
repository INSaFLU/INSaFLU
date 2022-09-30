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
from utils.software import Contigs2Sequences
from managing_files.models import Reference
from constants.constants import TypePath, FileExtensions
from utils.utils import Utils
from django.db import transaction
from django.conf import settings
import logging, os

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
		## Constants.DIR_TYPE_IDENTIFICATION = "db/type_identification/"
		b_test = False
		(version, path) = uploadFiles.get_file_to_upload(b_test)
		
		## upload
		uploadFile = uploadFiles.upload_file(version, path)

		# create the Abricate database
		b_testing = False
		contigs2sequences = Contigs2Sequences(b_testing)
		software= Software()
		if (not uploadFile is None):
			try:
				if (not software.is_exist_database_abricate(uploadFile.abricate_name)):
					self.stdout.write("Create new Abricate database: {}".format(uploadFile.abricate_name))
					software.create_database_abricate(uploadFile.abricate_name, uploadFile.path)
			except Exception as e:
				## try to create this database
				software.create_database_abricate(uploadFile.abricate_name, uploadFile.path)

		## tests other abricate database
		### get database file name, if it is not passed
		## Constants.DIR_TYPE_CONTIGS_2_SEQUENCES = "db/contigs2sequences/"
		(version, database_file_name) = contigs2sequences.get_most_recent_database()
		database_name = contigs2sequences.get_database_name()
		
		### first create database
		if (not software.is_exist_database_abricate(database_name)):
			software.create_database_abricate(database_name, database_file_name)
	
	### A command must define handle()
	def test_gff_and_bed_files(self):
		
		utils = Utils()
		software = Software()
		count = 0
		for reference in Reference.objects.all():
			if reference.is_deleted: continue
			count += 1
		
			### create bed and index for genbank
			test_file = reference.get_reference_bed(TypePath.MEDIA_ROOT)
			if not os.path.exists(reference.get_reference_gbk(TypePath.MEDIA_ROOT)):	## remove this reference
				
				self.stdout.write("Reference removed: '{}' owner: '{}'".format(reference.name,
									reference.owner.username))
				reference.is_deleted = True;
				reference.save()
				continue
			
			if not os.path.exists(test_file): utils.from_genbank_to_bed(reference.get_reference_gbk(TypePath.MEDIA_ROOT), test_file)
			
			test_file = reference.get_gff3(TypePath.MEDIA_ROOT)
			if not os.path.exists(test_file): software.run_genbank2gff3(reference.get_reference_gbk(TypePath.MEDIA_ROOT), test_file)
			
			test_file = reference.get_gff3_with_gene_annotation(TypePath.MEDIA_ROOT)
			if not os.path.exists(test_file): 
				with_gene_annotation = True
				software.run_genbank2gff3(reference.get_reference_gbk(TypePath.MEDIA_ROOT), test_file, with_gene_annotation)
			
			test_file = reference.get_gff3_comulative_positions(TypePath.MEDIA_ROOT)
			if not os.path.exists(test_file):
				software.run_genbank2gff3_positions_comulative(reference.get_reference_gbk(TypePath.MEDIA_ROOT), test_file)
			
			test_file = reference.get_reference_bed(TypePath.MEDIA_ROOT) + FileExtensions.FILE_IDX
			if not os.path.exists(test_file):
				software.create_index_files_from_igv_tools(reference.get_reference_bed(TypePath.MEDIA_ROOT))
		self.stdout.write("Number of references processed: {}".format(count))

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
	
	
	def set_species_tags_to_refences(self):
		"""
		set species tags to all references
		"""
		software = Software()
		count = 0
		for reference in Reference.objects.all():
			if reference.is_deleted: continue
			count += 1

			## set species tag		
			software.get_species_tag(reference)

		self.stdout.write("Number of references processed on species tag: {}".format(count))


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
		
		#### set all species tags to reference
		self.stdout.write("Set species tag to all references")
		self.set_species_tags_to_refences()
		
		#### set default references
		self.stdout.write("Test bed and gff files for all references")
		self.test_gff_and_bed_files()
		
		self.stdout.write("End")
		