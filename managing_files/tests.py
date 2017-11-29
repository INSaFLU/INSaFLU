from django.test import TestCase
from utils.constantsTestsCase import ConstantsTestsCase
from django.conf import settings 
from utils.utils import Utils
from django.contrib.auth.models import User
from utils.constants import FileType, TypePath, FileExtensions
from managing_files.manage_database import ManageDatabase
from managing_files.models import Sample, MetaKey, ProjectSample, Project, Reference
import os, time

class testsReferenceFiles(TestCase):
	
	def test_fasta_and_gb_file(self):
		"""
		test fasta and genbank
		"""
		utils = Utils()
		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
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
		
		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FAIL_FASTA)
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
		
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
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
		
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK_MISS_ONE)
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
		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		
		## compare fasta and genbank
		try:
			count_seqs = utils.compare_locus_fasta_gb(fasta_file, gb_file)
			self.assertEqual(8, count_seqs)
		except ValueError as e:
			self.fail(e.message)
			
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK_MISS_ONE)
		try:
			utils.compare_locus_fasta_gb(fasta_file, gb_file)
			self.fail("Must throw exception")
		except ValueError as e:
			self.assertEqual("Number of locus are different from fasta to genbank.", e.args[0])
			
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK_DIFF_LENGTH)
		try:
			utils.compare_locus_fasta_gb(fasta_file, gb_file)
			self.fail("Must throw exception")
		except ValueError as e:
			self.assertEqual("Different length. Fasta seq: PB2 length: 2280; Fasta seq: PB2 length: 2277.", e.args[0])

		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA_2)
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK_2)
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


	def test_meta_key_project_sample(self):
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

