'''
Created on Nov 27, 2017

@author: mmp
'''
import unittest, os, filecmp, datetime
from datetime import date
from django.conf import settings 
from utils.utils import Utils
from constants.constantsTestsCase import ConstantsTestsCase
from managing_files.manage_database import ManageDatabase
from datasets.manage_database import ManageDatabase as ManageDatabaseDataset
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference, TagName, TagNames, CountVariations
from datasets.models import Dataset, DatasetConsensus, UploadFiles, MetaKey
from constants.meta_key_and_values import MetaKeyAndValue
from utils.collect_dataset_data import CollectExtraDatasetData
from utils.parse_coverage_file import GetCoverage
from django.test.utils import override_settings
from constants.constants import TypePath, FileType, TypeFile
from constants.software_names import SoftwareNames
from utils.result import Result, SoftwareDesc

class Test(unittest.TestCase):

	constants_tests_case = ConstantsTestsCase()
	utils = Utils()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass

	def tearDown(self):
		pass

	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_create_tree_and_alignments(self):
		"""
 		test global method
 		"""
		software_names = SoftwareNames()
		manage_database = ManageDatabase()
		manage_database_datasets = ManageDatabaseDataset()
		
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT_TEST", None))
		path_destination = os.path.join(getattr(settings, "MEDIA_ROOT_TEST", None), ConstantsTestsCase.DIR_PROJECTS)
		self.utils.make_path(os.path.join(getattr(settings, "MEDIA_ROOT_TEST", None), ConstantsTestsCase.DIR_PROJECTS))
		cmd = 'cp -r {}/{}/* {}'.format(self.baseDirectory, ConstantsTestsCase.DIR_PROJECTS, path_destination)
		os.system(cmd)
		
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME + '5000')
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME + '5000'
			user.id = 5000
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

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
			project.id= 5000
			project.name = project_name
			project.reference = reference
			project.owner = user
			project.save()
		
		## get all fastq files
		vect_files = self.constants_tests_case.get_all_fastq_files(self.baseDirectory)
		
		tag_name_name = 'xpto'
		try:
			tag_name = TagName.objects.get(name='xpto')
		except TagName.DoesNotExist as e:
			tag_name = TagName()
			tag_name.owner = user
			tag_name.name = tag_name_name
			tag_name.is_meta_data = False
			tag_name.save()
		
		get_coverage = GetCoverage()
		ProjectSample.objects.all().delete()
		Sample.objects.all().delete()
		temp_dir = self.utils.get_temp_dir()
		count = 0
		
		for vect_file in vect_files:
			self.utils.copy_file(vect_file[0], os.path.join(temp_dir, os.path.basename(vect_file[0])))
			self.utils.copy_file(vect_file[1], os.path.join(temp_dir, os.path.basename(vect_file[1])))
			
			n_id = 5000 + count + 1
			sample_name = "_".join(os.path.basename(vect_file[0]).split('_')[0:2])
			try:
				sample = Sample.objects.get(pk = n_id)
			except Sample.DoesNotExist as e:
				sample = Sample()
				sample.id = n_id
				sample.name = sample_name
				sample.is_valid_1 = True
				sample.file_name_1 = os.path.basename(vect_file[0])
				sample.path_name_1.name = os.path.join(temp_dir, os.path.basename(vect_file[0]))
				sample.is_valid_2 = False
				sample.file_name_2 = os.path.basename(vect_file[1])
				sample.path_name_2.name = os.path.join(temp_dir, os.path.basename(vect_file[1]))
				sample.owner = user
				sample.is_ready_for_projects = True
				sample.is_obsolete = False
				sample.type_of_fastq = Sample.TYPE_OF_FASTQ_illumina
				sample.type_subtype = 'xpto, zpto'
				sample.save()

				## KEY META_KEY_Identify_Sample_Software and META_KEY_Fastq_Trimmomatic_Software
				parameters = "let's go again";
				result_all_2 = Result()
				result_all_2.add_software(SoftwareDesc(software_names.get_fastq_name(), 
						software_names.get_fastq_version(), parameters))
				manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Identify_Sample_Software,
					MetaKeyAndValue.META_VALUE_Success, result_all_2.to_json())
				
				result_all_2.add_software(SoftwareDesc(software_names.get_trimmomatic_name(), 
						software_names.get_trimmomatic_version(), parameters))
				result_all_2.add_software(SoftwareDesc(software_names.get_trimmomatic_name(), 
						software_names.get_trimmomatic_version(), parameters + "second version trimmomatic"))
				manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software,
					MetaKeyAndValue.META_VALUE_Success, result_all_2.to_json())
			## add tag names to sample
			tag_names = TagNames()
			tag_names.value = tag_name_name + " " + tag_name_name
			tag_names.tag_name = tag_name
			tag_names.sample = sample
			tag_names.save()
			
			result_all = Result()
			n_id = 5000 + count + 1
			try:
				project_sample = ProjectSample.objects.get(pk = n_id)
			except ProjectSample.DoesNotExist as e:
				## create project_sample
				project_sample = ProjectSample()
				project_sample.id = n_id
				project_sample.sample = sample
				project_sample.project = project
				project_sample.seq_name_all_consensus = project_sample.sample.name
				project_sample.is_finished = True
				project_sample.is_deleted = False
				project_sample.is_error = False
			
				count_variations = CountVariations()
				count_variations.var_bigger_50_90 = n_id + 1
				count_variations.var_bigger_90 = n_id + 12
				count_variations.var_less_50 = 12 + count
				count_variations.total = n_id + 1 + n_id + 12
				count_variations.save()
				project_sample.count_variations = count_variations
				project_sample.save()
				
				### KEY META_KEY_Snippy_Freebayes
				parameters = "let's go";
				result_all.add_software(SoftwareDesc(software_names.get_spades_name(), 
						software_names.get_spades_version(), parameters))
				parameters = "let's go";
				result_all.add_software(SoftwareDesc(software_names.get_snippy_name(), 
						software_names.get_snippy_version(), parameters + "_2222"))
				### don't save this one to stay empty in the file
				if (not vect_file[0].endswith("EVA011_S54_L001_R1_001.fastq.gz")):
					manage_database.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Snippy_Freebayes, 
						MetaKeyAndValue.META_VALUE_Success, result_all.to_json())
				
			coverage = get_coverage.get_coverage(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ,\
						software_names.get_snippy_name()), project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT))
			meta_value = manage_database.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage,
					MetaKeyAndValue.META_VALUE_Success, coverage.to_json())
			count += 1

		### test group set tag names
		self.assertTrue(project_sample.sample.get_tag_names().count() == 1) 

		### set dataset name
		dataset = Dataset()
		dataset.name = "xpto"
		dataset.owner = user
		dataset.save()
		
		consensus_add, reference_add = 0, 0
		for project_sample in ProjectSample.objects.filter(project=project, is_deleted= False):
			## don't add project samples not added
			if not project_sample.is_finished: continue
			try:
				dataset_consensus = DatasetConsensus.objects.get(project_sample=project_sample,
									dataset=dataset)
				
				if dataset_consensus.is_deleted or dataset_consensus.is_error:
					dataset_consensus.is_deleted = False
					dataset_consensus.is_error = False
					dataset_consensus.save()
					consensus_add += 1
				## now is finished and before is not finished yet
				elif (project_sample.is_finished and not dataset_consensus.is_project_sample_finished):
					dataset_consensus.is_project_sample_finished = True
					consensus_add += 1
				continue
			except DatasetConsensus.DoesNotExist:
				dataset_consensus = DatasetConsensus()
				dataset_consensus.name = project_sample.sample.name
				dataset_consensus.type_subtype = project_sample.sample.type_subtype
				dataset_consensus.dataset = dataset
				dataset_consensus.project_sample = project_sample
				dataset_consensus.is_project_sample_finished = project_sample.is_finished
				dataset_consensus.save() 
				consensus_add += 1
			
		### Add the reference of this project if not there yet
		try:
			dataset_consensus = DatasetConsensus.objects.get(reference=project.reference,
								dataset=dataset)
			if dataset_consensus.is_deleted or dataset_consensus.is_error:
				dataset_consensus.is_deleted = False
				dataset_consensus.is_error = False
				dataset_consensus.save()
				reference_add += 1
		except DatasetConsensus.DoesNotExist:
			dataset_consensus = DatasetConsensus()
			dataset_consensus.name = project.reference.name
			dataset_consensus.dataset = dataset
			dataset_consensus.reference = project.reference
			dataset_consensus.save() 
			reference_add += 1

		### necessary to calculate the global results again 
		if (reference_add > 0):
			dataset.last_change_date = datetime.datetime.now()
			dataset.number_of_sequences_from_projects += consensus_add
			dataset.number_of_sequences_from_references += reference_add
			dataset.totla_alerts = 1 if dataset.get_number_different_references() > 1 else 0
			dataset.save()

		## test reference name
		self.assertEqual(ref_name, dataset.get_first_reference_name())
		
		collect_extra_data = CollectExtraDatasetData()
		collect_extra_data.collect_extra_data_for_dataset(dataset, user)

		meta_data = manage_database_datasets.get_dataset_metakey(dataset, MetaKeyAndValue.META_KEY_Dataset_max_name_length, MetaKeyAndValue.META_VALUE_Success)
		self.assertEqual("32", meta_data.description)
		
		### test nextStrain metadata
		expected_file_nextstrain = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_DATASET_FILES, "out_nextstrain_metadata_expected.tsv")
		out_file = dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_RESULT_NEXTSTRAIN_TSV)
		self.assertTrue(os.path.exists(out_file))
		
		## update date
		temp_file_compare = self.utils.get_temp_file_from_dir(os.path.dirname(out_file), "nextstrain_compare", ".txt")
		vect_data = self.utils.read_text_file(out_file)
		with open(temp_file_compare, 'w') as handle_write:
			for line in vect_data:
				handle_write.write(line.strip().replace(date.today().strftime(settings.DATE_FORMAT_FOR_SHOW), "2022-10-24") + "\n")
		self.assertTrue(filecmp.cmp(temp_file_compare, expected_file_nextstrain))
		
		### get sample result file
		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))

	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_create_tree_and_alignments_with_metadata(self):
		"""
 		test global method
 		"""
		software_names = SoftwareNames()
		manage_database = ManageDatabase()
		manage_database_datasets = ManageDatabaseDataset()
		
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT_TEST", None))
		path_destination = os.path.join(getattr(settings, "MEDIA_ROOT_TEST", None), ConstantsTestsCase.DIR_PROJECTS)
		self.utils.make_path(os.path.join(getattr(settings, "MEDIA_ROOT_TEST", None), ConstantsTestsCase.DIR_PROJECTS))
		cmd = 'cp -r {}/{}/* {}'.format(self.baseDirectory, ConstantsTestsCase.DIR_PROJECTS, path_destination)
		os.system(cmd)
		
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME + '5000')
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME + '5000'
			user.id = 5000
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

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
			project.id= 5000
			project.name = project_name
			project.reference = reference
			project.owner = user
			project.save()
		
		## get all fastq files
		vect_files = self.constants_tests_case.get_all_fastq_files(self.baseDirectory)
		
		tag_name_name = 'xpto'
		try:
			tag_name = TagName.objects.get(name='xpto')
		except TagName.DoesNotExist as e:
			tag_name = TagName()
			tag_name.owner = user
			tag_name.name = tag_name_name
			tag_name.is_meta_data = False
			tag_name.save()
		
		get_coverage = GetCoverage()
		ProjectSample.objects.all().delete()
		Sample.objects.all().delete()
		temp_dir = self.utils.get_temp_dir()
		count = 0
		
		for vect_file in vect_files:
			self.utils.copy_file(vect_file[0], os.path.join(temp_dir, os.path.basename(vect_file[0])))
			self.utils.copy_file(vect_file[1], os.path.join(temp_dir, os.path.basename(vect_file[1])))
			
			n_id = 5000 + count + 1
			sample_name = "_".join(os.path.basename(vect_file[0]).split('_')[0:2])
			try:
				sample = Sample.objects.get(pk = n_id)
			except Sample.DoesNotExist as e:
				sample = Sample()
				sample.id = n_id
				sample.name = sample_name
				sample.is_valid_1 = True
				sample.file_name_1 = os.path.basename(vect_file[0])
				sample.path_name_1.name = os.path.join(temp_dir, os.path.basename(vect_file[0]))
				sample.is_valid_2 = False
				sample.file_name_2 = os.path.basename(vect_file[1])
				sample.path_name_2.name = os.path.join(temp_dir, os.path.basename(vect_file[1]))
				sample.owner = user
				sample.is_ready_for_projects = True
				sample.is_obsolete = False
				sample.type_of_fastq = Sample.TYPE_OF_FASTQ_illumina
				sample.type_subtype = 'xpto, zpto'
				sample.save()

				## KEY META_KEY_Identify_Sample_Software and META_KEY_Fastq_Trimmomatic_Software
				parameters = "let's go again";
				result_all_2 = Result()
				result_all_2.add_software(SoftwareDesc(software_names.get_fastq_name(), 
						software_names.get_fastq_version(), parameters))
				manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Identify_Sample_Software,
					MetaKeyAndValue.META_VALUE_Success, result_all_2.to_json())
				
				result_all_2.add_software(SoftwareDesc(software_names.get_trimmomatic_name(), 
						software_names.get_trimmomatic_version(), parameters))
				result_all_2.add_software(SoftwareDesc(software_names.get_trimmomatic_name(), 
						software_names.get_trimmomatic_version(), parameters + "second version trimmomatic"))
				manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software,
					MetaKeyAndValue.META_VALUE_Success, result_all_2.to_json())
			## add tag names to sample
			tag_names = TagNames()
			tag_names.value = tag_name_name + " " + tag_name_name
			tag_names.tag_name = tag_name
			tag_names.sample = sample
			tag_names.save()
			
			result_all = Result()
			n_id = 5000 + count + 1
			try:
				project_sample = ProjectSample.objects.get(pk = n_id)
			except ProjectSample.DoesNotExist as e:
				## create project_sample
				project_sample = ProjectSample()
				project_sample.id = n_id
				project_sample.sample = sample
				project_sample.project = project
				project_sample.seq_name_all_consensus = project_sample.sample.name
				project_sample.is_finished = True
				project_sample.is_deleted = False
				project_sample.is_error = False
			
				count_variations = CountVariations()
				count_variations.var_bigger_50_90 = n_id + 1
				count_variations.var_bigger_90 = n_id + 12
				count_variations.var_less_50 = 12 + count
				count_variations.total = n_id + 1 + n_id + 12
				count_variations.save()
				project_sample.count_variations = count_variations
				project_sample.save()
				
				### KEY META_KEY_Snippy_Freebayes
				parameters = "let's go";
				result_all.add_software(SoftwareDesc(software_names.get_spades_name(), 
						software_names.get_spades_version(), parameters))
				parameters = "let's go";
				result_all.add_software(SoftwareDesc(software_names.get_snippy_name(), 
						software_names.get_snippy_version(), parameters + "_2222"))
				### don't save this one to stay empty in the file
				if (not vect_file[0].endswith("EVA011_S54_L001_R1_001.fastq.gz")):
					manage_database.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Snippy_Freebayes, 
						MetaKeyAndValue.META_VALUE_Success, result_all.to_json())
				
			coverage = get_coverage.get_coverage(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ,\
						software_names.get_snippy_name()), project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT))
			meta_value = manage_database.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage,
					MetaKeyAndValue.META_VALUE_Success, coverage.to_json())
			count += 1

		### test group set tag names
		self.assertTrue(project_sample.sample.get_tag_names().count() == 1) 

		### set dataset name
		dataset = Dataset()
		dataset.name = "xpto"
		dataset.owner = user
		dataset.save()
		
		consensus_add, reference_add = 0, 0
		for project_sample in ProjectSample.objects.filter(project=project, is_deleted= False):
			## don't add project samples not added
			if not project_sample.is_finished: continue
			try:
				dataset_consensus = DatasetConsensus.objects.get(project_sample=project_sample,
									dataset=dataset)
				
				if dataset_consensus.is_deleted or dataset_consensus.is_error:
					dataset_consensus.is_deleted = False
					dataset_consensus.is_error = False
					dataset_consensus.save()
					consensus_add += 1
				## now is finished and before is not finished yet
				elif (project_sample.is_finished and not dataset_consensus.is_project_sample_finished):
					dataset_consensus.is_project_sample_finished = True
					consensus_add += 1
				continue
			except DatasetConsensus.DoesNotExist:
				dataset_consensus = DatasetConsensus()
				dataset_consensus.name = project_sample.sample.name
				dataset_consensus.type_subtype = project_sample.sample.type_subtype
				dataset_consensus.dataset = dataset
				dataset_consensus.project_sample = project_sample
				dataset_consensus.is_project_sample_finished = project_sample.is_finished
				dataset_consensus.save() 
				consensus_add += 1
			
		### Add the reference of this project if not there yet
		try:
			dataset_consensus = DatasetConsensus.objects.get(reference=project.reference,
								dataset=dataset)
			if dataset_consensus.is_deleted or dataset_consensus.is_error:
				dataset_consensus.is_deleted = False
				dataset_consensus.is_error = False
				dataset_consensus.save()
				reference_add += 1
		except DatasetConsensus.DoesNotExist:
			dataset_consensus = DatasetConsensus()
			dataset_consensus.name = project.reference.name
			dataset_consensus.dataset = dataset
			dataset_consensus.reference = project.reference
			dataset_consensus.save() 
			reference_add += 1

		### necessary to calculate the global results again 
		if (reference_add > 0):
			dataset.last_change_date = datetime.datetime.now()
			dataset.number_of_sequences_from_projects += consensus_add
			dataset.number_of_sequences_from_references += reference_add
			dataset.totla_alerts = 1 if dataset.get_number_different_references() > 1 else 0
			dataset.save()

		## add metadata file
		type_file_name = TypeFile.TYPE_FILE_dataset_file_metadata
		try:
			type_file = MetaKey.objects.get(name=type_file_name)
		except MetaKey.DoesNotExist as e:
			type_file = MetaKey()
			type_file.name = type_file_name
			type_file.save()
			
		nextstrain_metadata = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_DATASET_FILES, "nextstrain_metadata_input.tsv")
		self.assertTrue(os.path.exists(nextstrain_metadata))
		upload_files = UploadFiles()
		upload_files.path_name.name = nextstrain_metadata
		upload_files.file_name = os.path.basename(nextstrain_metadata)
		upload_files.dataset = dataset
		upload_files.owner = user
		upload_files.type_file = type_file
		upload_files.is_valid = True
		upload_files.save()
		
		### check if it was saved
		upload_metadata_file = UploadFiles.objects.filter(owner=user, is_deleted=False,\
            			type_file__name=TypeFile.TYPE_FILE_dataset_file_metadata, is_valid=True,
            	  		dataset=dataset).order_by('-creation_date').first()
		self.assertFalse(upload_metadata_file is None)
		
		## test reference name
		self.assertEqual(ref_name, dataset.get_first_reference_name())
		
		collect_extra_data = CollectExtraDatasetData()
		collect_extra_data.collect_extra_data_for_dataset(dataset, user)

		meta_data = manage_database_datasets.get_dataset_metakey(dataset, MetaKeyAndValue.META_KEY_Dataset_max_name_length, MetaKeyAndValue.META_VALUE_Success)
		self.assertEqual("32", meta_data.description)
		
		### test nextStrain metadata
		expected_file_nextstrain = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_DATASET_FILES, "out_nextstrain_metadata_expected_2.tsv")
		out_file = dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_RESULT_NEXTSTRAIN_TSV)
		self.assertTrue(os.path.exists(out_file))
		
		## update date
		temp_file_compare = self.utils.get_temp_file_from_dir(os.path.dirname(out_file), "nextstrain_compare", ".txt")
		vect_data = self.utils.read_text_file(out_file)
		with open(temp_file_compare, 'w') as handle_write:
			for line in vect_data:
				handle_write.write(line.strip().replace(date.today().strftime(settings.DATE_FORMAT_FOR_SHOW), "2022-10-24") + "\n")
		self.assertTrue(filecmp.cmp(temp_file_compare, expected_file_nextstrain))
		
		upload_metadata_file = UploadFiles.objects.filter(owner=user, is_deleted=False,\
            			type_file__name=TypeFile.TYPE_FILE_dataset_file_metadata, is_valid=True,
            	  		dataset=dataset).order_by('-creation_date').first()
		self.assertFalse(upload_metadata_file is None)
		self.assertEqual(True, upload_metadata_file.is_processed)
		self.assertEqual(True, upload_metadata_file.is_valid)
		
		### get sample result file
		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))
