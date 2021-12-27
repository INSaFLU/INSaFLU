'''
Created on 01/05/2020

@author: mmp
'''
from django.test import TransactionTestCase
from django.conf import settings 
from constants.constantsTestsCase import ConstantsTestsCase
from constants.constants_mixed_infection import ConstantsMixedInfection
from constants.constants import Constants, TypePath, FileType, FileExtensions, TypeFile
from constants.meta_key_and_values import MetaKeyAndValue
from constants.tag_names_constants import TagNamesConstants
from utils.software import Software, Contigs2Sequences
from constants.software_names import SoftwareNames
from utils.utils import Utils
from utils.parse_out_files import ParseOutFiles
from utils.result import DecodeObjects, Coverage, CountHits
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference, Statistics, MetaKey
from managing_files.manage_database import ManageDatabase
from manage_virus import uploadFiles
from django.test.utils import override_settings
from utils.tree import CreateTree
from utils.parse_in_files import ParseInFiles
from utils.result import DecodeObjects, MixedInfectionMainVector
from managing_files.models import CountVariations, MixedInfections
from utils.mixed_infections_management import MixedInfectionsManagement
from manage_virus.constants_virus import ConstantsVirus
from utils.process_SGE import ProcessSGE
from extend_user.models import Profile
import os, filecmp, csv, ntpath

class TestMultiProcess(TransactionTestCase):


	### static
	software = Software()
	utils = Utils()
	constants = Constants()
	constants_tests_case = ConstantsTestsCase()
	software_names = SoftwareNames()

	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		
		from manage_virus.uploadFiles import UploadFiles
		uploadFiles = UploadFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		self.assertEqual("2", version)
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
		uploadFiles.upload_file(version, file)	## upload file


	def tearDown(self):
		pass
	
	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_multi_files(self):
		"""
		test multi files
		"""
		from managing_files.models import UploadFiles
		manageDatabase = ManageDatabase()
		metaKeyAndValue = MetaKeyAndValue()
		
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))
		self.utils.remove_dir(settings.MEDIA_ROOT_TEST)
		self.utils.remove_dir(os.path.join(Constants.TEMP_DIRECTORY, Constants.COUNT_DNA_TEMP_DIRECTORY))
		
		csv_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_INPUT_FILES, ConstantsTestsCase.MANAGING_FILES_TEMPLATE_MULTIPLE_FILES_data_csv)
		self.assertTrue(os.path.exists(csv_file))
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.set_password(ConstantsTestsCase.TEST_USER_NAME)
			user.save()
			
		parse_in_files = ParseInFiles()
		b_test_char_encoding = False
		parse_in_files.parse_sample_files(csv_file, user, b_test_char_encoding, ParseInFiles.STATE_READ_all)
		self.assertEquals(0, parse_in_files.get_errors().get_len_vect_results())
		self.assertEquals(6, len(parse_in_files.get_vect_samples()))
		
		### upload multiple sample file...
		upload_files = UploadFiles()
		
		upload_files.is_valid = True
		upload_files.is_processed = False
		upload_files.is_deleted = False
		upload_files.number_errors = 0
		upload_files.number_files_processed = 0
		upload_files.number_files_to_process = len(parse_in_files.get_vect_samples())
			
		try:
			type_file_sample = MetaKey.objects.get(name=TypeFile.TYPE_FILE_sample_file)
		except MetaKey.DoesNotExist:
			type_file_sample = MetaKey()
			type_file_sample.name = TypeFile.TYPE_FILE_sample_file
			type_file_sample.save()
		
		upload_files.type_file = type_file_sample
		upload_files.file_name = ConstantsTestsCase.MANAGING_FILES_TEMPLATE_MULTIPLE_FILES_data_csv
		upload_files.owner = user
		upload_files.description = ""
		
		## move the files to the right place
		sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), self.utils.get_path_upload_file(user.id,\
												TypeFile.TYPE_FILE_sample_file), upload_files.file_name)
		sz_file_to = self.utils.get_unique_file(sz_file_to)		## get unique file name, user can upload files with same name...
		self.utils.copy_file(csv_file, sz_file_to)
		upload_files.path_name.name = os.path.join(self.utils.get_path_upload_file(user.id,\
									TypeFile.TYPE_FILE_sample_file), ntpath.basename(sz_file_to))
		upload_files.save()
			
		vect_wait_sge_ids = []
		try:
			process_SGE = ProcessSGE()
			b_test = True
			vect_wait_sge_ids.append(process_SGE.set_read_sample_file(upload_files, user, b_test))
		except:
			self.fail("Upload process description file in SGE")
		
		### wait for all sge ids			
		process_SGE.wait_until_finished(vect_wait_sge_ids)
		self.assertTrue(len(vect_wait_sge_ids) == 0)

		### check samples
		self.assertEqual(6, Sample.objects.all().count())
		self.assertEqual(6, Sample.objects.all().filter(is_deleted=False, is_valid_1=False, owner=user).count())

		### start upload files...
		vect_wait_sge_ids.append(self.upload_file(ConstantsTestsCase.FASTQ1_1, user))
		vect_wait_sge_ids.append(self.upload_file(ConstantsTestsCase.FASTQ1_2, user))
		vect_wait_sge_ids.append(self.upload_file(ConstantsTestsCase.FASTQ2_1, user))
		vect_wait_sge_ids.append(self.upload_file(ConstantsTestsCase.FASTQ2_2, user))
		vect_wait_sge_ids.append(self.upload_file(ConstantsTestsCase.FASTQ3_1, user))
		vect_wait_sge_ids.append(self.upload_file(ConstantsTestsCase.FASTQ3_2, user))
		vect_wait_sge_ids.append(self.upload_file(ConstantsTestsCase.FASTQ4_1, user))
		vect_wait_sge_ids.append(self.upload_file(ConstantsTestsCase.FASTQ4_2, user))
		vect_wait_sge_ids.append(self.upload_file(ConstantsTestsCase.FASTQ10_1, user))
		vect_wait_sge_ids.append(self.upload_file(ConstantsTestsCase.FASTQ10_2, user))
		vect_wait_sge_ids.append(self.upload_file(ConstantsTestsCase.FASTQ11_1, user))
		vect_wait_sge_ids.append(self.upload_file(ConstantsTestsCase.FASTQ11_2, user))

		### wait for all sge ids
		process_SGE.wait_until_finished(vect_wait_sge_ids)
		self.assertTrue(len(vect_wait_sge_ids) == 0)
		
		### need to wait for all process for link
		self.assertEqual(6, Sample.objects.all().filter(is_deleted=False, is_valid_1=True, is_valid_2=True, owner=user).count())
		self.assertNotEqual(6, Sample.objects.all().filter(is_deleted=False, is_valid_1=True, is_valid_2=True, is_ready_for_projects=True, owner=user).count())
		uploaf_file = UploadFiles.objects.get(type_file=type_file_sample, owner=user)
		self.assertEqual(6, uploaf_file.number_files_processed)
		
		### test result trimmomatic
		process_SGE.wait_until_finished([])
		
		### test if all of them processed
		self.assertEqual(6, Sample.objects.all().filter(is_deleted=False, is_valid_1=True, is_valid_2=True, is_ready_for_projects=True, owner=user).count())

		for sample in Sample.objects.all().filter(is_deleted=False, is_valid_1=True, is_valid_2=True, is_ready_for_projects=True, owner=user):
			if sample.pk == 6: self.assertEqual(sample.type_subtype, "A-H1N1")
			elif sample.pk == 5: self.assertEqual(sample.type_subtype, "A-H3N2")
			elif sample.pk == 4: self.assertEqual(sample.type_subtype, "A-H3N2")
			elif sample.pk == 3: self.assertEqual(sample.type_subtype, "A-H3N2")
			elif sample.pk == 2: self.assertEqual(sample.type_subtype, "A-H3N2")
			elif sample.pk == 1: self.assertEqual(sample.type_subtype, "A-H3N2")
			else: self.fail("Fail fot this sample.pk {}".format(sample.pk))

		### start making trees
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		ref_name = "second_stage_2_ test_create_tree"
		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.display_name = ref_name
			reference.reference_fasta.name = fasta_file
			reference.reference_fasta_name = os.path.basename(fasta_file)
			reference.reference_genbank.name = gb_file
			reference.reference_genbank_name = os.path.basename(gb_file)
			reference.owner = user
			reference.save()
		
		project_name = "several_names_test_create_tree"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = project_name
			project.reference = reference
			project.owner = user
			project.save()
		
		####
		(job_name_wait, job_name) = ("", "")
		for sample in Sample.objects.all().filter(is_deleted=False, is_valid_1=True, is_valid_2=True, is_ready_for_projects=True, owner=user):
			## create project_sample
			project_sample = ProjectSample()
			project_sample.project = project
			project_sample.sample = sample
			project_sample.save()
			
			print(job_name_wait, job_name)
			if len(job_name_wait) == 0: (job_name_wait, job_name) = user.profile.get_name_sge_seq(Profile.SGE_GLOBAL)
			taskID = process_SGE.set_second_stage_snippy(project_sample, user, job_name, job_name_wait)
				
			### set project sample queue ID
			manageDatabase.set_project_sample_metakey(project_sample, user,\
							metaKeyAndValue.get_meta_key_queue_by_project_sample_id(project_sample.id),\
							MetaKeyAndValue.META_VALUE_Queue, taskID)
		
		taskID = process_SGE.set_collect_global_files(project, user)
		manageDatabase.set_project_metakey(project, user, metaKeyAndValue.get_meta_key(\
							MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id), MetaKeyAndValue.META_VALUE_Queue, taskID)

		### test result trimmomatic
		process_SGE.wait_until_finished([])

	def upload_file(self, file_name, user):
		
		from managing_files.models import UploadFiles
		fastq_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, file_name)
		upload_files = UploadFiles()
		upload_files.file_name = file_name
		
		try:
			type_file = MetaKey.objects.get(name=TypeFile.TYPE_FILE_fastq_gz)
		except MetaKey.DoesNotExist:
			type_file = MetaKey()
			type_file.name = TypeFile.TYPE_FILE_fastq_gz
			type_file.save()

		upload_files.is_valid = True
		upload_files.is_processed = False			## True when all samples are set
		upload_files.owner = user
		upload_files.type_file = type_file
		upload_files.number_files_to_process = 1
		upload_files.number_files_processed = 0
		upload_files.description = ""
		
		sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), self.utils.get_path_upload_file(user.id,\
												TypeFile.TYPE_FILE_fastq_gz), upload_files.file_name)
		sz_file_to = self.utils.get_unique_file(sz_file_to)		## get unique file name, user can upload files with same name...
		self.utils.copy_file(fastq_file, sz_file_to)
		upload_files.path_name.name = os.path.join(self.utils.get_path_upload_file(user.id,\
									TypeFile.TYPE_FILE_fastq_gz), ntpath.basename(sz_file_to))
		
		upload_files.save()
		
		process_SGE = ProcessSGE()
		b_test = True
		taskID = process_SGE.set_link_files(user, b_test)
		return taskID





