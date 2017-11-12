from django.test import TestCase
from utils.constantsTestsCase import ConstantsTestsCase
from django.conf import settings 
from utils.utils import Utils
from django.contrib.auth.models import User
from managing_files.manage_database import ManageDatabase
from managing_files.models import Sample, MetaKey
import os, time

class testsReferenceFiles(TestCase):
	
	def test_fasta_and_gb_file(self):
		"""
		test fasta and genbank
		"""
		utils = Utils()
		constantsTestsCase = ConstantsTestsCase()
		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, constantsTestsCase.MANAGING_DIR, constantsTestsCase.MANAGING_FILES_FASTA)
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
		
		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FAIL_FASTA)
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
		
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
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
		
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK_MISS_ONE)
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
		constantsTestsCase = ConstantsTestsCase()
		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, constantsTestsCase.MANAGING_DIR, constantsTestsCase.MANAGING_FILES_FASTA)
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		
		## compare fasta and genbank
		try:
			count_seqs = utils.compare_locus_fasta_gb(fasta_file, gb_file)
			self.assertEqual(8, count_seqs)
		except ValueError as e:
			self.fail(e.message)
			
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK_MISS_ONE)
		try:
			utils.compare_locus_fasta_gb(fasta_file, gb_file)
			self.fail("Must throw exception")
		except ValueError as e:
			self.assertEqual("Number of locus are different from fasta to genbank.", e.args[0])
			
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK_DIFF_LENGTH)
		try:
			utils.compare_locus_fasta_gb(fasta_file, gb_file)
			self.fail("Must throw exception")
		except ValueError as e:
			self.assertEqual("Different length. Fasta seq: PB2 length: 2280; Fasta seq: PB2 length: 2277.", e.args[0])

		fasta_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, constantsTestsCase.MANAGING_DIR, constantsTestsCase.MANAGING_FILES_FASTA_2)
		gb_file = os.path.join(getattr(settings, "STATIC_ROOT", None), constantsTestsCase.MANAGING_TESTS, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK_2)
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



