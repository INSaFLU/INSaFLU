from django.test import TestCase
from constants.constantsTestsCase import ConstantsTestsCase
from django.conf import settings 
from utils.utils import Utils
from django.contrib.auth.models import User
from constants.constants import FileType, TypePath, FileExtensions
from django.core.exceptions import MultipleObjectsReturned
from managing_files.manage_database import ManageDatabase
from managing_files.models import Sample, MetaKey, ProjectSample, Project, Reference, Statistics, CountVariations, MetaKeyProject
from constants.meta_key_and_values import MetaKeyAndValue
from constants.tag_names_constants import TagNamesConstants
from utils.result import CountHits
import os, time

class testsReferenceFiles(TestCase):
	
	constantsTestsCase = ConstantsTestsCase()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), self.constantsTestsCase.MANAGING_TESTS)
		
	def test_fasta_and_gb_file(self):
		"""
		test fasta and genbank
		"""
		utils = Utils()
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		try:
			count_seqs = utils.is_fasta(fasta_file)
			self.assertEqual(8, count_seqs)
		except IOError as e:
			self.fail(e.args)
			
		try:
			self.assertFalse(utils.is_genbank(fasta_file))
			self.fail("Must throw exception")
		except IOError as e:
			pass
		
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FAIL_FASTA)
		try:
			utils.is_fasta(fasta_file)
			self.fail("Must throw exception")
		except IOError as e:
			pass
			
		try:
			self.assertFalse(utils.is_genbank(fasta_file))
			self.fail("Must throw exception")
		except IOError as e:
			pass
		
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		try:
			utils.is_fasta(gb_file)
			self.fail("Must throw exception")
		except IOError as e:
			pass
			
		try:
			count_seqs = utils.is_genbank(gb_file)
			self.assertEqual(8, count_seqs)
		except IOError as e:
			self.fail(e.args)
		
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK_MISS_ONE)
		try:
			utils.is_fasta(gb_file)
			self.fail("Must throw exception")
		except IOError as e:
			pass
			
		try:
			count_seqs = utils.is_genbank(gb_file)
			self.assertEqual(7, count_seqs)
		except IOError as e:
			self.fail(e.args)


	def test_compare_locus_fasta_gb(self):
		"""
		Compare if fasta file has the same name of locus
		"""
		utils = Utils()
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		
		## compare fasta and genbank
		try:
			count_seqs = utils.compare_locus_fasta_gb(fasta_file, gb_file)
			self.assertEqual(8, count_seqs)
		except ValueError as e:
			self.fail(e.message)
			
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK_MISS_ONE)
		try:
			utils.compare_locus_fasta_gb(fasta_file, gb_file)
			self.fail("Must throw exception")
		except ValueError as e:
			self.assertEqual("Number of locus are different from fasta to genbank.", e.args[0])
			
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK_DIFF_LENGTH)
		try:
			utils.compare_locus_fasta_gb(fasta_file, gb_file)
			self.fail("Must throw exception")
		except ValueError as e:
			self.assertEqual("Different length. Fasta seq: PB2 length: 2280; Fasta seq: PB2 length: 2277.", e.args[0])

		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA_2)
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK_2)
		try:
			utils.compare_locus_fasta_gb(fasta_file, gb_file)
			self.fail("Must throw exception")
		except ValueError as e:
			self.assertEqual("This locus 'locus_1' is not in fasta file.", e.args[0])
		except:
			self.fail("throw other exception")

	def test_meta_key(self):
		"""
		test the metakey system
		"""
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()
			
			
		sample_name = "xpto"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_rejected = False
			sample.path_name_1.name = "sdasa"
			sample.owner = user
			sample.save()
		
		manageDatabase = ManageDatabase()
		manageDatabase.set_metakey(sample, user, ConstantsTestsCase.META_KEY_TEST, ConstantsTestsCase.VALUE_TEST, "description")
		
		try:
			metaKey = MetaKey.objects.get(name=ConstantsTestsCase.META_KEY_TEST)
			self.assertEqual(ConstantsTestsCase.META_KEY_TEST, metaKey.name)
		except MetaKey.DoesNotExist:
			self.fail("must exist the meta_key")
			
		metaKey_sample = manageDatabase.get_metakey(sample, ConstantsTestsCase.META_KEY_TEST, ConstantsTestsCase.VALUE_TEST)
		self.assertFalse(metaKey_sample == None)
		self.assertEqual(sample_name, metaKey_sample.sample.name)
		self.assertEqual(ConstantsTestsCase.META_KEY_TEST, metaKey_sample.meta_tag.name)
		self.assertEqual(ConstantsTestsCase.VALUE_TEST, metaKey_sample.value)
		self.assertEqual("description", metaKey_sample.description)

		manageDatabase.set_metakey(sample, user, ConstantsTestsCase.META_KEY_TEST, ConstantsTestsCase.VALUE_TEST_2, "description")
		metaKey_sample = manageDatabase.get_metakey(sample, ConstantsTestsCase.META_KEY_TEST, ConstantsTestsCase.VALUE_TEST_2)
		self.assertFalse(metaKey_sample == None)
		self.assertEqual(sample_name, metaKey_sample.sample.name)
		self.assertEqual(ConstantsTestsCase.META_KEY_TEST, metaKey_sample.meta_tag.name)
		self.assertEqual(ConstantsTestsCase.VALUE_TEST_2, metaKey_sample.value)
		self.assertEqual("description", metaKey_sample.description)
		
		
		metaKey_sample_lst = manageDatabase.get_metakey(sample, ConstantsTestsCase.META_KEY_TEST, None)
		self.assertEquals(2, metaKey_sample_lst.count())
		self.assertEqual(sample_name, metaKey_sample_lst[0].sample.name)
		self.assertEqual(ConstantsTestsCase.VALUE_TEST_2, metaKey_sample_lst[0].value)
		self.assertEqual(ConstantsTestsCase.VALUE_TEST, metaKey_sample_lst[1].value)
		
	def test_meta_key_order(self):
		"""
		test the metakey system
		"""
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()
			
			
		sample_name = "xpto"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_rejected = False
			sample.path_name_1.name = "sdasa"
			sample.owner = user
			sample.save()
		
		manageDatabase = ManageDatabase()
		manageDatabase.set_metakey(sample, user, ConstantsTestsCase.META_KEY_TEST, ConstantsTestsCase.VALUE_TEST, "description")
		
		try:
			metaKey = MetaKey.objects.get(name=ConstantsTestsCase.META_KEY_TEST)
			self.assertEqual(ConstantsTestsCase.META_KEY_TEST, metaKey.name)
		except MetaKey.DoesNotExist:
			self.fail("must exist the meta_key")
			
		metaKey_sample = manageDatabase.get_metakey(sample, ConstantsTestsCase.META_KEY_TEST, ConstantsTestsCase.VALUE_TEST)
		self.assertFalse(metaKey_sample == None)
		self.assertEqual(sample_name, metaKey_sample.sample.name)
		self.assertEqual(ConstantsTestsCase.META_KEY_TEST, metaKey_sample.meta_tag.name)
		self.assertEqual(ConstantsTestsCase.VALUE_TEST, metaKey_sample.value)
		self.assertEqual("description", metaKey_sample.description)
		
		time.sleep(0.5)
		manageDatabase.set_metakey(sample, user, ConstantsTestsCase.META_KEY_TEST, ConstantsTestsCase.VALUE_TEST, "description_2")
		list_data = manageDatabase.get_metakey(sample, ConstantsTestsCase.META_KEY_TEST, None)
		self.assertEquals(2, list_data.count())
		self.assertEquals("description_2", list_data[0].description)
		self.assertEquals("description", list_data[1].description)


	def test_get_file_output(self):
		"""
		test the metakey system
		"""
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()
			
		sample_name = "file_name_2"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_rejected = False
			sample.path_name_1.name = "/tmp/zpt"
			sample.owner = user
			sample.save()
			
		project_name = "file_name_3"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = sample_name
			project.owner = user
			project.save()
		
		project_sample = ProjectSample()
		project_sample.sample = sample
		project_sample.project = project
		project_sample.save()
		path_to_project = 'projects/result/user_{0}/project_{1}/sample_{2}'.format(sample.owner.id, project.id, sample.id)
		self.assertEquals(os.path.join(settings.MEDIA_URL, path_to_project, "file_name_2.bam"), project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_BAM, ""))
		self.assertEquals(os.path.join(settings.MEDIA_URL, path_to_project, "file_name_2.bam.bai"), project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_BAM_BAI, ""))
		self.assertEquals(os.path.join(settings.MEDIA_URL, path_to_project, "file_name_2.consensus.fa"), project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_CONSENSUS_FA, ""))
		self.assertEquals(os.path.join(settings.MEDIA_URL, path_to_project, "file_name_2.consensus.fasta"), project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_CONSENSUS_FASTA, ""))
		self.assertEquals(os.path.join(settings.MEDIA_URL, path_to_project, "file_name_2.csv"), project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_CSV, ""))
		self.assertEquals(os.path.join(settings.MEDIA_URL, path_to_project, "file_name_2.depth.gz"), project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_DEPTH_GZ, ""))
		self.assertEquals(os.path.join(settings.MEDIA_URL, path_to_project, "file_name_2.depth.gz.tbi"), project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_DEPTH_GZ_TBI, ""))
		self.assertEquals(os.path.join(settings.MEDIA_URL, path_to_project, "file_name_2.tab"), project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_TAB, ""))
		self.assertEquals(os.path.join(settings.MEDIA_URL, path_to_project, "file_name_2.vcf"), project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_VCF, ""))
		self.assertEquals(os.path.join(settings.MEDIA_URL, path_to_project, "file_name_2.vcf.gz"), project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_VCF_GZ, ""))
		self.assertEquals(os.path.join(settings.MEDIA_URL, path_to_project, "file_name_2.vcf.gz.tbi"), project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_VCF_GZ_TBI, ""))



	def test_meta_key_project(self):
		"""
		test the metakey system
		"""
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		ref_name = "second_stage"
		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.reference_fasta.name = "fasta_file"
			reference.reference_fasta_name = "os.path.basename(fasta_file)"
			reference.reference_genbank.name = "gb_file"
			reference.reference_genbank_name = "os.path.basename(gb_file)"
			reference.owner = user
			reference.save()
			
		project_name = "xpto_2342"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = project_name
			project.reference = reference
			project.owner = user
			project.save()

		key_test = ConstantsTestsCase.META_KEY_TEST + "_sdcdd"
		manageDatabase = ManageDatabase()
		manageDatabase.set_project_metakey(project, user, key_test, ConstantsTestsCase.VALUE_TEST, "description")

		try:
			metaKey = MetaKey.objects.get(name=key_test)
			self.assertEqual(key_test, metaKey.name)
		except MetaKey.DoesNotExist:
			self.fail("must exist the meta_key")
			
		metaKey_sample = manageDatabase.get_project_metakey(project, key_test, ConstantsTestsCase.VALUE_TEST)
		self.assertFalse(metaKey_sample == None)
		self.assertEqual(project_name, metaKey_sample.project.name)
		self.assertEqual(key_test, metaKey_sample.meta_tag.name)
		self.assertEqual(ConstantsTestsCase.VALUE_TEST, metaKey_sample.value)
		self.assertEqual("description", metaKey_sample.description)
		
		time.sleep(0.5)
		manageDatabase.set_project_metakey(project, user, key_test, key_test, "description_2")
		list_data = manageDatabase.get_project_metakey(project, key_test, None)
		self.assertEquals(2, list_data.count())
		self.assertEquals("description_2", list_data[0].description)
		self.assertEquals("description", list_data[1].description)

		metaKey_sample = manageDatabase.get_project_metakey(project, key_test, ConstantsTestsCase.VALUE_TEST)
		self.assertFalse(metaKey_sample == None)
		manageDatabase.set_project_metakey(project, user, key_test, ConstantsTestsCase.VALUE_TEST, "description")
		try:
			metaKey_sample = manageDatabase.get_project_metakey(project, key_test, ConstantsTestsCase.VALUE_TEST)
			self.fail('must throw error...')
		except MultipleObjectsReturned as e:
			pass
		get_last_meta = manageDatabase.get_project_metakey_last(project, key_test, ConstantsTestsCase.VALUE_TEST)
		self.assertFalse(get_last_meta == None)
		self.assertTrue(get_last_meta.creation_date > metaKey_sample.creation_date)
		get_last_meta = manageDatabase.get_project_metakey_last(project, key_test, ConstantsTestsCase.VALUE_TEST + 'sddfsdfsdfs')
		self.assertTrue(get_last_meta == None)


	def test_meta_key_project_sample(self):
		"""
		test the metakey system
		"""
		
		meta_key_and_value = MetaKeyAndValue()
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()
		
		ref_name = "second_stage"
		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.reference_fasta.name = "fasta_file"
			reference.reference_fasta_name = "os.path.basename(fasta_file)"
			reference.reference_genbank.name = "gb_file"
			reference.reference_genbank_name = "os.path.basename(gb_file)"
			reference.owner = user
			reference.save()
	
		project_name = "xpto_2342"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = project_name
			project.reference = reference
			project.owner = user
			project.save()

		sample_name = "xpto"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_rejected = False
			sample.path_name_1.name = "sdasa"
			sample.owner = user
			sample.save()

		### project
		project_sample = ProjectSample()
		project_sample.sample = sample
		project_sample.project = project
		project_sample.save()
		
		manageDatabase = ManageDatabase()
		manageDatabase.set_project_sample_metakey(project_sample, user, ConstantsTestsCase.META_KEY_TEST + "xrpt", ConstantsTestsCase.VALUE_TEST, "description")
		
		try:
			metaKey = MetaKey.objects.get(name=ConstantsTestsCase.META_KEY_TEST + "xrpt")
			self.assertEqual(ConstantsTestsCase.META_KEY_TEST + "xrpt", metaKey.name)
		except MetaKey.DoesNotExist:
			self.fail("must exist the meta_key")

		metaKey_sample = manageDatabase.get_project_sample_metakey(project_sample, ConstantsTestsCase.META_KEY_TEST, ConstantsTestsCase.VALUE_TEST)
		self.assertTrue(metaKey_sample == None)

		metaKey_sample = manageDatabase.get_project_sample_metakey(project_sample, ConstantsTestsCase.META_KEY_TEST + "xrpt", ConstantsTestsCase.VALUE_TEST)
		self.assertFalse(metaKey_sample == None)
		self.assertEqual(sample_name, metaKey_sample.project_sample.sample.name)
		self.assertEqual(project_name, metaKey_sample.project_sample.project.name)
		self.assertEqual(ConstantsTestsCase.META_KEY_TEST + "xrpt", metaKey_sample.meta_tag.name)
		self.assertEqual(ConstantsTestsCase.VALUE_TEST, metaKey_sample.value)
		self.assertEqual("description", metaKey_sample.description)
		
		time.sleep(0.5)
		manageDatabase.set_project_sample_metakey(project_sample, user, ConstantsTestsCase.META_KEY_TEST + "xrpt", ConstantsTestsCase.VALUE_TEST, "description_2")
		list_data = manageDatabase.get_project_sample_metakey(project_sample, ConstantsTestsCase.META_KEY_TEST + "xrpt", None)
		self.assertEquals(2, list_data.count())
		self.assertEquals("description_2", list_data[0].description)
		self.assertEquals("description", list_data[1].description)

		metaKey_sample_2 = manageDatabase.get_project_sample_metakey_last(project_sample, ConstantsTestsCase.META_KEY_TEST + "xrpt", ConstantsTestsCase.VALUE_TEST)
		self.assertFalse(metaKey_sample_2 == None)
		self.assertTrue(metaKey_sample_2.creation_date > metaKey_sample.creation_date)

		### test alerts
		meta_key_value = meta_key_and_value.get_meta_key(MetaKeyAndValue.META_KEY_ALERT_COVERAGE_9, 'xpto')
		self.assertTrue(meta_key_value.startswith(MetaKeyAndValue.META_KEY_ALERT_COVERAGE) and meta_key_value.endswith('xpto'))
		manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_value, ConstantsTestsCase.VALUE_TEST + 'ddd', 'xxprorot')
		meta_key_value = meta_key_and_value.get_meta_key(MetaKeyAndValue.META_KEY_ALERT_COVERAGE_9, 'xpto1')
		manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_value, ConstantsTestsCase.VALUE_TEST, 'xxprorotw')
		meta_key_value = meta_key_and_value.get_meta_key(MetaKeyAndValue.META_KEY_ALERT_COVERAGE_0, 'xpto2')
		manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_value, ConstantsTestsCase.VALUE_TEST, 'xxprorotqwqwwq')
		result_list = manageDatabase.get_project_sample_starts_with_metakey(project_sample, MetaKeyAndValue.META_KEY_ALERT_COVERAGE)
		self.assertEquals(3, len(result_list))
		self.assertEquals('xxprorotqwqwwq', result_list[0].description)
		self.assertEquals(ConstantsTestsCase.VALUE_TEST, result_list[0].value)
		self.assertEquals(ConstantsTestsCase.VALUE_TEST + 'ddd', result_list[-1].value)
		self.assertEquals('xxprorot', result_list[-1].description)
		
	def test_project_sample(self):
		"""
		test the metakey system
		"""
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME + "eewewerw")
		except User.DoesNotExist:
			user = User()
			user.id=1000
			user.username = ConstantsTestsCase.TEST_USER_NAME + "eewewerw"
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()
		
		ref_name = "seco$ncds_d_stage_"
		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.reference_fasta.name = "fasta_file"
			reference.reference_fasta_name = "os.path.basename(fasta_file)"
			reference.reference_genbank.name = "gb_file"
			reference.reference_genbank_name = "os.path.basename(gb_file)"
			reference.owner = user
			reference.save()
	
		project_name = "xp_1_3_$to_2342"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.id=1000
			project.name = project_name
			project.reference = reference
			project.owner = user
			project.save()

		sample_name = "xpt_2_4_f_6o"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.id=1000
			sample.name = sample_name
			sample.is_rejected = False
			sample.path_name_1.name = "sdasa"
			sample.owner = user
			sample.save()

		### project
		project_sample = ProjectSample()
		project_sample.sample = sample
		project_sample.project = project
		project_sample.save()
		
		software = None
		self.assertTrue(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, software).\
						endswith("projects/result/user_1000/project_1000/sample_1000/" + sample_name + ".vcf.gz"))
		software = ""
		self.assertTrue(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, software).\
						endswith("projects/result/user_1000/project_1000/sample_1000/" + sample_name + ".vcf.gz"))
		software = "None"
		self.assertTrue(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, software).\
						endswith("projects/result/user_1000/project_1000/sample_1000/none/" + sample_name.lower() + ".vcf.gz"))
		
		self.assertTrue(project_sample.get_global_file_by_element(TypePath.MEDIA_ROOT, ProjectSample.PREFIX_FILE_COVERAGE, 'sequence_name', FileExtensions.FILE_PNG).\
						endswith(os.path.join("projects/result/user_1000/project_1000/sample_1000", ProjectSample.PATH_MAIN_RESULT,\
						'{}_{}{}'.format(ProjectSample.PREFIX_FILE_COVERAGE, 'sequence_name', FileExtensions.FILE_PNG))))
		
		### output global out files by element and type
		self.assertTrue(project.get_global_file_by_project(TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_MAFFT).\
						endswith("projects/result/user_1000/project_1000/" + Project.PATH_MAIN_RESULT + "/" +\
						Project.PROJECT_FILE_NAME_MAFFT))
		
		self.assertEquals(None, project.get_global_file_by_element(TypePath.MEDIA_ROOT, 'xpto', 'xpto'))

		self.assertTrue(project.get_global_file_by_element(TypePath.MEDIA_ROOT, 'xpto', Project.PROJECT_FILE_NAME_MAFFT).\
						endswith("projects/result/user_1000/project_1000/" + Project.PATH_MAIN_RESULT +\
						'/xpto/' + Project.PROJECT_FILE_NAME_MAFFT_element + '_xpto' + FileExtensions.FILE_FASTA))

		self.assertTrue(project.get_global_file_by_element(TypePath.MEDIA_ROOT, 'xpto', Project.PROJECT_FILE_NAME_FASTTREE_tree).\
						endswith("projects/result/user_1000/project_1000/" + Project.PATH_MAIN_RESULT +\
						'/xpto/' + Project.PROJECT_FILE_NAME_FASTTREE_element + '_xpto' + FileExtensions.FILE_TREE))



	def test_percentil(self):
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		ref_name = "second"
		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.reference_fasta.name = ""
			reference.reference_fasta_name = ""
			reference.reference_genbank.name = ""
			reference.reference_genbank_name = ""
			reference.owner = user
			reference.save()
		
		project_name = "several"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.id= 5000
			project.name = project_name
			project.reference = reference
			project.owner = user
			project.save()
		
		sample_name = "_fddsdffds"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_rejected = False
			sample.is_valid_1 = True
			sample.file_name_1 = 'vect_file[0]'
			sample.path_name_1.name = ''
			sample.is_valid_2 = False
			sample.file_name_2 = 'vect_file[1]'
			sample.path_name_2.name = ''
			sample.owner = user
			
			sample.is_ready_for_projects = True
			sample.is_obsolete = False
			sample.is_rejected = False
			sample.save()


		manage_database = ManageDatabase()
		tagNamesConstants = TagNamesConstants()
		self.assertEquals(None, manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_98, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL)))
		self.assertEquals(None, manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_98, TagNamesConstants.TAG_PERCENTIL_VAR_50)))
		self.assertEquals(None, manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_98, TagNamesConstants.TAG_PERCENTIL_VAR_50_90)))
		self.assertEquals(None, manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_95, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL)))
		self.assertEquals(None, manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_95, TagNamesConstants.TAG_PERCENTIL_VAR_50)))
		self.assertEquals(None, manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_95, TagNamesConstants.TAG_PERCENTIL_VAR_50_90)))
		self.assertEquals(None, manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_80, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL)))
		self.assertEquals(None, manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_80, TagNamesConstants.TAG_PERCENTIL_VAR_50)))
		self.assertEquals(None, manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_80, TagNamesConstants.TAG_PERCENTIL_VAR_50_90)))

		value_ = manage_database.get_percentil_value(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_98, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL))
		self.assertEquals(33, value_)
		
		total_numbers = 1000
		project_sample_last = None
		for i in range(0, total_numbers):
			count_hits = CountHits()
			count_hits.set_hits_less_50((i * 999) % 2)
			count_hits.set_hits_50_90((i * 2999) % 3)
			
			## create project_sample
			project_sample = ProjectSample()
			project_sample.sample = sample
			project_sample.project = project
			project_sample.is_finished = True
			project_sample.is_deleted = False
			project_sample.is_error = False
			project_sample.save()
			self.assertTrue(manage_database.add_variation_count(project_sample, user, count_hits))
			project_sample_last = project_sample

		### test if insert an None countHits
		self.assertFalse(manage_database.add_variation_count(project_sample_last, user, None))

		self.assertEquals(total_numbers, CountVariations.objects.all().count())
		raw_ = manage_database.get_percentil_value(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_98, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL))
		self.assertEquals(3, raw_)
		raw_ = manage_database.get_percentil_value(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_98, TagNamesConstants.TAG_PERCENTIL_VAR_50))
		self.assertEquals(1, raw_)
		raw_ = manage_database.get_percentil_value(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_98, TagNamesConstants.TAG_PERCENTIL_VAR_50_90))
		self.assertEquals(2, raw_)
		raw_ = manage_database.get_percentil_value(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_95, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL))
		self.assertEquals(3, raw_)
		raw_ = manage_database.get_percentil_value(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_95, TagNamesConstants.TAG_PERCENTIL_VAR_50))
		self.assertEquals(1, raw_)
		raw_ = manage_database.get_percentil_value(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_95, TagNamesConstants.TAG_PERCENTIL_VAR_50_90))
		self.assertEquals(2, raw_)
		raw_ = manage_database.get_percentil_value(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_80, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL))
		self.assertEquals(2, raw_)
		raw_ = manage_database.get_percentil_value(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_80, TagNamesConstants.TAG_PERCENTIL_VAR_50))
		self.assertEquals(1, raw_)
		raw_ = manage_database.get_percentil_value(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_80, TagNamesConstants.TAG_PERCENTIL_VAR_50_90))
		self.assertEquals(2, raw_)

		### get all statistics
		self.assertEquals(18, Statistics.objects.all().count())
		
		### get directly from database
		raw_ = manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_98, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL))
		self.assertEquals(3.0, raw_)
		raw_ = manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_98, TagNamesConstants.TAG_PERCENTIL_VAR_50))
		self.assertEquals(1.0, raw_)
		raw_ = manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_98, TagNamesConstants.TAG_PERCENTIL_VAR_50_90))
		self.assertEquals(2.0, raw_)
		raw_ = manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_95, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL))
		self.assertEquals(3.0, raw_)
		raw_ = manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_95, TagNamesConstants.TAG_PERCENTIL_VAR_50))
		self.assertEquals(1.0, raw_)
		raw_ = manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_95, TagNamesConstants.TAG_PERCENTIL_VAR_50_90))
		self.assertEquals(2.0, raw_)
		raw_ = manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_80, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL))
		self.assertEquals(2.0, raw_)
		raw_ = manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_80, TagNamesConstants.TAG_PERCENTIL_VAR_50))
		self.assertEquals(1.0, raw_)
		raw_ = manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_80, TagNamesConstants.TAG_PERCENTIL_VAR_50_90))
		self.assertEquals(2.0, raw_)
		
		raw_ = manage_database.get_percentil_value_from_db(tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_80, TagNamesConstants.TAG_PERCENTIL_VAR_50_90) + "dfsdd")
		self.assertEquals(None, raw_)
		
		count_hits = CountHits()
		count_hits.set_hits_50_90(1)
		count_hits.set_hits_less_50(0)
		###
		percentil_name = tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_80, TagNamesConstants.TAG_PERCENTIL_VAR_50_90)
		manage_database.set_percentis_alert(project_sample_last, user, count_hits, percentil_name)
		return_data = manage_database.get_project_sample_starts_with_metakey(project_sample_last, MetaKeyAndValue.META_KEY_ALERT_COUNT_VAR)
		self.assertEquals(0, len(return_data))

		count_hits.set_hits_50_90(2)
		count_hits.set_hits_less_50(0)
		manage_database.set_percentis_alert(project_sample_last, user, count_hits, percentil_name)
		return_data = manage_database.get_project_sample_starts_with_metakey(project_sample_last, MetaKeyAndValue.META_KEY_ALERT_COUNT_VAR)
		self.assertEquals(0, len(return_data))
		
		count_hits.set_hits_50_90(3)
		count_hits.set_hits_less_50(0)
		manage_database.set_percentis_alert(project_sample_last, user, count_hits, percentil_name)
		return_data = manage_database.get_project_sample_starts_with_metakey(project_sample_last, MetaKeyAndValue.META_KEY_ALERT_COUNT_VAR)
		self.assertEquals(1, len(return_data))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Error, return_data[0].value)
		self.assertEquals("Warning, this sample has a total of variations '3' bigger than the percentile 80 of all database", return_data[0].description)
		
		count_hits.set_hits_50_90(4)
		count_hits.set_hits_less_50(0)
		percentil_name = tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_98, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL)
		manage_database.set_percentis_alert(project_sample_last, user, count_hits, percentil_name)
		return_data = manage_database.get_project_sample_starts_with_metakey(project_sample_last, MetaKeyAndValue.META_KEY_ALERT_COUNT_VAR)
		self.assertEquals(2, len(return_data))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Error, return_data[0].value)
		self.assertEquals("Warning, this sample has a total of variations '4' bigger than the percentile 98 of all database", return_data[0].description)
		self.assertEquals("Warning, this sample has a total of variations '3' bigger than the percentile 80 of all database", return_data[1].description)


		### test the remove count variations
		self.assertTrue(manage_database.delete_variation_count(project_sample, user))
		self.assertFalse(manage_database.delete_variation_count(project_sample, user))


	def test_get_elements_and_genes_from_db(self):
		
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		ref_name = "second_stage_te_tree"
		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.reference_fasta.name = fasta_file
			reference.reference_fasta_name = os.path.basename(fasta_file)
			reference.reference_genbank.name = gb_file
			reference.reference_genbank_name = os.path.basename(gb_file)
			reference.owner = user
			reference.save()
		
		project_name = "several_tree"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = project_name
			project.reference = reference
			project.owner = user
			project.save()
			
		manage_database = ManageDatabase()
		metaKeyAndValue = MetaKeyAndValue()
		meta_key = metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_Elements_Project, project.id)
		meta_project = manage_database.get_project_metakey(project, meta_key, MetaKeyAndValue.META_VALUE_Success)
		self.assertEqual(None, meta_project)
	
		vect_data = manage_database.get_elements_and_genes_from_db(project, user)
		self.assertEqual(8, len(vect_data))
		self.assertEqual('HA', vect_data[0])
		self.assertEqual('PB2', vect_data[-1])

		meta_project = manage_database.get_project_metakey(project, meta_key, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_project != None)
		self.assertEqual(8, len(meta_project.description.split(',')))
		self.assertEqual('HA', meta_project.description.split(',')[0])
		self.assertEqual('PB2', meta_project.description.split(',')[-1])



