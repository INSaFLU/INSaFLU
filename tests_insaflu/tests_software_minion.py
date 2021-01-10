'''
Created on 01/01/2021

@author: mmp
'''
import os
from django.test import TestCase
from django.conf import settings
from constants.constantsTestsCase import ConstantsTestsCase
from django.contrib.auth.models import User
from managing_files.models import Sample, Reference, Project, ProjectSample
from constants.constants import Constants, TypePath, FileType, FileExtensions
from managing_files.manage_database import ManageDatabase
from utils.result import DecodeObjects, Coverage
from constants.meta_key_and_values import MetaKeyAndValue
from utils.software_minion import SoftwareMinion
from utils.utils import Utils
from constants.software_names import SoftwareNames
from django.test.utils import override_settings

class Test(TestCase):

	software_minion = SoftwareMinion()
	utils = Utils()
	constants = Constants()
	software_names = SoftwareNames()
	constants_tests_case = ConstantsTestsCase()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)

	def tearDown(self):
		pass


# 	def test_run_clean_minion(self):
# 		"""
# 		Test run fastq and trimmomatic all together
#  		"""
# 		file_name = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ_MINION_1)
# 		self.assertTrue(file_name)
# 		try:
# 			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
# 		except User.DoesNotExist:
# 			user = User()
# 			user.username = ConstantsTestsCase.TEST_USER_NAME
# 			user.is_active = False
# 			user.password = ConstantsTestsCase.TEST_USER_NAME
# 			user.save()
# 
# 		temp_dir = self.utils.get_temp_dir()
# 		self.utils.copy_file(file_name, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
# 			
# 		sample_name = "run_nanofilt"
# 		try:
# 			sample = Sample.objects.get(name=sample_name)
# 		except Sample.DoesNotExist:
# 			sample = Sample()
# 			sample.name = sample_name
# 			sample.is_valid_1 = True
# 			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
# 			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
# 			sample.is_valid_2 = False
# 			sample.owner = user
# 			sample.save()
# 		
# 		### run software
# 		self.assertTrue(self.software_minion.run_nanofilt_and_stat(sample, user))
# 		
# 		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
# 		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_nanofilt_file(TypePath.MEDIA_ROOT)))))
# 		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_rabbitQC_output(TypePath.MEDIA_ROOT)))))
# 		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_rabbitQC_nanofilt(TypePath.MEDIA_ROOT)))))
# 		
# 		manageDatabase = ManageDatabase()
# 		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None)
# 		self.assertTrue(len(list_meta) == 1)
# 		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
# 		self.assertEquals(MetaKeyAndValue.META_KEY_Number_And_Average_Reads, list_meta[0].meta_tag.name)
# 		
# 		### number of sequences
# 		manageDatabase = ManageDatabase()
# 		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None)
# 		self.assertTrue(len(list_meta) == 1)
# 		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
# 		self.assertEquals(MetaKeyAndValue.META_KEY_Number_And_Average_Reads, list_meta[0].meta_tag.name)
# 
# 		decodeResultAverageAndNumberReads = DecodeObjects()
# 		result_average = decodeResultAverageAndNumberReads.decode_result(list_meta[0].description)
# 		self.assertEqual('969', result_average.number_file_1)
# 		self.assertEqual('468.1', result_average.average_file_1)
# 		self.assertEqual(None, result_average.number_file_2)
# 		self.assertEqual(None, result_average.average_file_2)
# 
# 		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
# 		self.assertTrue(len(list_meta) == 0)
# 		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, None)
# 		self.assertTrue(len(list_meta) == 1)
# 		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
# 		self.assertEquals(MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, list_meta[0].meta_tag.name)
# 		self.assertEquals("Success, NanoStat(1.4.0), NanoFilt(2.6.0)", list_meta[0].description)
# 		
# 		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, None)
# 		self.assertTrue(len(list_meta) == 1)
# 		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
# 		self.assertEquals(MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, list_meta[0].meta_tag.name)
# 		result_data = decodeResultAverageAndNumberReads.decode_result(list_meta[0].description)
# 		vect_soft = result_data.get_list_software_instance(self.software_names.get_NanoStat_name())
# 		self.assertEqual(2, len(vect_soft))
# 		
# 		self.assertEqual("Mean read length", vect_soft[0].get_vect_key_values()[0].key)
# 		self.assertEqual("489.3", vect_soft[0].get_vect_key_values()[0].value)
# 		self.assertEqual("Total bases", vect_soft[0].get_vect_key_values()[-1].key)
# 		self.assertEqual("517,186.0", vect_soft[0].get_vect_key_values()[-1].value)
# 		self.assertEqual("Mean read length", vect_soft[1].get_vect_key_values()[0].key)
# 		self.assertEqual("468.1", vect_soft[1].get_vect_key_values()[0].value)
# 		self.assertEqual("Total bases", vect_soft[1].get_vect_key_values()[-1].key)
# 		self.assertEqual("453,569.0", vect_soft[1].get_vect_key_values()[-1].value)
# 		
# 		## remove all files
# 		cmd = "rm -r %s*" % (temp_dir); os.system(cmd)
		
		
	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_process_second_stage_medaka(self):
		"""
 		test global method
 		"""
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT_TEST", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))

		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_COVID_GBK)
		reference_fasta = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_COVID_FASTA)
		file_fasta = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_nanopore)

		self.assertTrue(os.path.exists(gb_file))
		self.assertTrue(os.path.exists(reference_fasta))
		self.assertTrue(os.path.exists(file_fasta))
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		ref_name = "second_stagis_single_covid"
		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.reference_fasta.name = reference_fasta
			reference.reference_fasta_name = os.path.basename(reference_fasta)
			reference.reference_genbank.name = gb_file
			reference.reference_genbank_name = os.path.basename(gb_file)
			reference.owner = user
			reference.save()
			
		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_fasta, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
			
		sample_name = "run_snippyis_single_covid"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = False
			sample.file_name_2 = None
			sample.path_name_2.name = None
			sample.owner = user
			sample.save()

		project_name = "file_naais_single_covid"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = project_name
			project.reference = reference
			project.owner = user
			project.save()
		
		## create project_sample
		project_sample = ProjectSample()
		project_sample.sample = sample
		project_sample.project = project
		project_sample.save()
		
		### copy to nanofilt file
		self.utils.copy_file(file_fasta, project_sample.sample.get_nanofilt_file(TypePath.MEDIA_ROOT))
		
		manageDatabase = ManageDatabase()
		metaKeyAndValue = MetaKeyAndValue()
		meta_key_project_sample = metaKeyAndValue.get_meta_key_queue_by_project_sample_id(project_sample.id)
		manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue, "meta_sample.description")
		self.assertTrue(self.software_minion.process_second_stage_medaka(project_sample, user))
		try:
			project_sample = ProjectSample.objects.get(pk=project_sample.id)
		except ProjectSample.DoesNotExist:
			self.fail("Must exist")
		self.assertTrue(project_sample.is_finished)
		self.assertTrue(project_sample.is_mask_consensus_sequences)
		
		### test the files
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertFalse(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FA, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, SoftwareNames.SOFTWARE_SNIPPY_name).index(SoftwareNames.SOFTWARE_SNIPPY_name.lower()) != -1)
		
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_BAM, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_BAM_BAI, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertFalse(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ_TBI, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CSV, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_REF_FASTA, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_REF_FASTA, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_REF_FASTA_FAI, SoftwareNames.SOFTWARE_SNIPPY_name)))
		
		## freebayes
		self.assertTrue(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_FREEBAYES_name).index(SoftwareNames.SOFTWARE_FREEBAYES_name.lower()) != -1)
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, SoftwareNames.SOFTWARE_FREEBAYES_name)))

		## test freebayes clean
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		self.assertTrue(os.path.getsize(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)) > 10)

		## test consensus file
		self.assertTrue(os.path.exists(project_sample.get_consensus_file(TypePath.MEDIA_ROOT)))
		self.assertTrue(os.path.getsize(project_sample.get_consensus_file(TypePath.MEDIA_ROOT)) > 10)

		### human file name, snippy tab
		self.assertTrue(os.path.exists(project_sample.get_file_output_human(TypePath.MEDIA_ROOT, FileType.FILE_TAB,  self.software_names.get_snippy_name())))
		self.assertTrue(os.path.getsize(project_sample.get_file_output_human(TypePath.MEDIA_ROOT, FileType.FILE_TAB,  self.software_names.get_snippy_name())) > 10)
		
		### set flag that is finished
		list_meta = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Count_Hits, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Count_Hits, list_meta[0].meta_tag.name)
		
		### get the hits value
		decode_coverage = DecodeObjects()
		count_hits = decode_coverage.decode_result(list_meta[0].description)
		self.assertEquals(0, count_hits.get_hits_50_90())
		self.assertEquals(3, count_hits.get_hits_less_50())
		self.assertEquals(3, count_hits.get_total_50_50_90())
		self.assertEquals(125, count_hits.get_total())
		
		self.assertEquals(3, project_sample.count_variations.var_less_50)
		self.assertEquals(0, project_sample.count_variations.var_bigger_50_90)
		self.assertEquals(122, project_sample.count_variations.var_bigger_90)
		self.assertEquals(0, project_sample.alert_first_level)
		self.assertEquals(0, project_sample.alert_second_level)
		
		### test mixed infections
		self.assertEquals('0.815100969586423', '{}'.format(project_sample.mixed_infections.average_value))
		self.assertEquals('No', project_sample.mixed_infections.tag.name)
		self.assertFalse(project_sample.mixed_infections.has_master_vector)
		
		list_meta = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Snippy_Freebayes, None)
		self.assertEquals(1, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Snippy_Freebayes, list_meta[0].meta_tag.name)
		
		decode_result = DecodeObjects()
		result = decode_result.decode_result(list_meta[0].description)
		self.assertTrue(result is not None)
		self.assertEquals("Snippy-3.2-dev; (--mapqual 20 --mincov 10 --minfrac 0.51)", result.get_software(self.software_names.get_snippy_name()))
		self.assertEquals("Freebayes-v1.1.0-54-g49413aa; (-p 2 -q 20 -m 20 --min-coverage 100 --min-alternate-fraction 0.01 --min-alternate-count 10 -V)",\
 						result.get_software(self.software_names.get_freebayes_name()))

		meta_value = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
		decode_coverage = DecodeObjects()
		coverage = decode_coverage.decode_result(meta_value.description)
		self.assertEqual(coverage.get_coverage('PA', Coverage.COVERAGE_ALL), "527.4")
		self.assertEqual(coverage.get_coverage('MP', Coverage.COVERAGE_ALL), "2198.8")
		self.assertEqual(coverage.get_coverage('HA', Coverage.COVERAGE_ALL), "1449.3")
		self.assertEqual(coverage.get_coverage('PA', Coverage.COVERAGE_MORE_9), "100.0")
		self.assertEqual(coverage.get_coverage('MP', Coverage.COVERAGE_MORE_9), "100.0")
		self.assertEqual(coverage.get_coverage('HA', Coverage.COVERAGE_MORE_9), "100.0")
		self.assertEqual(coverage.get_coverage('PA', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('MP', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('HA', Coverage.COVERAGE_MORE_0), "100.0")
		
		lst_meta_sample = manageDatabase.get_project_sample_metakey(project_sample, meta_key_project_sample, None)
		self.assertEquals(2, len(lst_meta_sample))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Queue, lst_meta_sample[1].value)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, lst_meta_sample[0].value)
		
		### check if the images exist
		for gene in self.utils.get_elements_from_db(reference, user):
			output_image = project_sample.get_global_file_by_element(TypePath.MEDIA_ROOT, ProjectSample.PREFIX_FILE_COVERAGE, gene, FileExtensions.FILE_PNG)
			self.assertTrue(os.path.exists(output_image))

		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))

