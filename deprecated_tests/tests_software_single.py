'''
Created on Oct 28, 2017

@author: mmp
'''

from django.test import TestCase
from django.conf import settings 
from constants.constantsTestsCase import ConstantsTestsCase
from constants.constants_mixed_infection import ConstantsMixedInfection
from constants.constants import Constants, TypePath, FileType, FileExtensions
from constants.meta_key_and_values import MetaKeyAndValue
from constants.tag_names_constants import TagNamesConstants
from utils.software import Software, Contigs2Sequences
from utils.software_pangolin import SoftwarePangolin
from constants.software_names import SoftwareNames
from utils.utils import Utils
from utils.parse_out_files import ParseOutFiles
from utils.result import DecodeObjects, Coverage, CountHits, GeneticElement, MaskingConsensus, Result, SoftwareDesc, KeyValue, KeyValues
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference, Statistics, TagNames
from manage_virus.uploadFiles import UploadFiles
from manage_virus.models import UploadFile
from managing_files.manage_database import ManageDatabase
from django.test.utils import override_settings
from utils.tree import CreateTree
import os, filecmp, csv
from utils.parse_in_files import ParseInFiles
from utils.result import DecodeObjects, MixedInfectionMainVector
from managing_files.models import CountVariations, MixedInfections
from utils.mixed_infections_management import MixedInfectionsManagement, MixedInfectionsTag
from manage_virus.constants_virus import ConstantsVirus
from manage_virus.models import Tags, SeqVirus, IdentifyVirus
from settings.default_software import DefaultSoftware
from settings.models import Software as SoftwareSettings, Parameter
from utils.parse_coverage_file import GetCoverage
from plotly.figure_factory._dendrogram import scs
from managing_files.models import Software as SoftwareModel, TagName
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from settings.constants_settings import ConstantsSettings
from settings.default_software_project_sample import DefaultProjectSoftware
from utils.software_minion import SoftwareMinion
from utils.collect_extra_data import CollectExtraData

class Test(TestCase):

	### static
	software = Software()
	software_minion = SoftwareMinion()
	software_names = SoftwareNames()
	software_pangolin = SoftwarePangolin()
	utils = Utils()
	constants = Constants()
	constants_tests_case = ConstantsTestsCase()

	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass


	def tearDown(self):
		pass
	
	def test_run_clean_minion_2(self):
		"""
		Test run fastq and trimmomatic all together
 		"""
		uploadFiles = UploadFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		self.assertEqual("2", version)
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
		uploadFiles.upload_file(version, file)	## upload file
		
		file_name = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, "covid/test_minion_seq.fastq.gz")
		self.assertTrue(file_name)
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_name, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
			
		sample_name = "run_nanofilt"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = False
			sample.type_of_fastq = Sample.TYPE_OF_FASTQ_minion
			sample.owner = user
			sample.save()
		
		### turn off abricate
		default_software = DefaultSoftware()
		default_software.test_all_defaults(user)
		self.assertTrue(default_software.is_software_to_run(SoftwareNames.SOFTWARE_ABRICATE_name,
							user, ConstantsSettings.TECHNOLOGY_minion))
		is_to_run = False
		default_software.set_software_to_run(SoftwareNames.SOFTWARE_ABRICATE_name, user,
							ConstantsSettings.TECHNOLOGY_minion, is_to_run)
		self.assertFalse(default_software.is_software_to_run(SoftwareNames.SOFTWARE_ABRICATE_name,
							user, ConstantsSettings.TECHNOLOGY_minion))

		b_make_identify_species = False
		self.assertTrue(self.software_minion.run_clean_minion(sample, user, b_make_identify_species))
		default_project_software = DefaultProjectSoftware()
		self.assertFalse(default_project_software.is_to_run_abricate(user, sample,
				ConstantsSettings.TECHNOLOGY_minion))
		
		try:
			software = SoftwareSettings.objects.get(name=SoftwareNames.SOFTWARE_NanoFilt_name, owner=user,\
							type_of_use = SoftwareSettings.TYPE_OF_USE_sample,
							technology__name = ConstantsSettings.TECHNOLOGY_minion,
							is_obsolete = False)
		except SoftwareSettings.DoesNotExist:
			self.fail("Must must exist")
		
		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(5, len(parameters))
		for parameter in parameters:
			if (parameter.name == "-q"):
				parameter.parameter = "5"
				parameter.save()
				break
		
		tag_name_mixed= "No"
		sample = Sample.objects.get(name=sample_name)
		if (b_make_identify_species):
			self.fail("Must have b_make_identify_species FALSE")
			self.assertEqual(tag_name_mixed, sample.mixed_infections_tag.name)
			self.assertEqual("A-H5N5", sample.type_subtype)
			
			sample.mixed_infections_tag = None
			sample.save()
			
			sample = Sample.objects.get(name=sample_name)
			self.assertEqual(None, sample.mixed_infections_tag)
			
			try:
				mixed_infections_tag = MixedInfectionsTag.objects.get(name=tag_name_mixed)
			except MixedInfectionsTag.DoesNotExist as e:
				self.fail("Must have tag mixed")
							
			if (not sample.mixed_infections_tag is None): sample.mixed_infections_tag.delete()
		else:
			self.assertEqual(Constants.EMPTY_VALUE_NA, sample.mixed_infections_tag.name)
			self.assertEqual(Constants.EMPTY_VALUE_NA, sample.type_subtype)
		
		self.assertEqual(0, sample.number_alerts)
		
		### run again
		default_project_software.set_abricate_to_run(user, sample,
				ConstantsSettings.TECHNOLOGY_minion, True)
		## global
		self.assertFalse(default_software.is_software_to_run(SoftwareNames.SOFTWARE_ABRICATE_name,
							user, ConstantsSettings.TECHNOLOGY_minion))
		### project sample
		self.assertTrue(default_project_software.is_to_run_abricate(user, 
					sample, ConstantsSettings.TECHNOLOGY_minion))
		self.assertTrue(self.software_minion.run_clean_minion(sample, user,
					default_project_software.is_to_run_abricate(user,
					sample, ConstantsSettings.TECHNOLOGY_minion) ))
		
		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(5, len(parameters))
		for parameter in parameters:
			if (parameter.name == "-q"):
				self.assertEqual("5", parameter.parameter)
				break;

		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_nanofilt_file(TypePath.MEDIA_ROOT)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_rabbitQC_output(TypePath.MEDIA_ROOT)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_rabbitQC_nanofilt(TypePath.MEDIA_ROOT)))))
		
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None)
		self.assertEqual(2, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Number_And_Average_Reads, list_meta[0].meta_tag.name)
		
		decodeResultAverageAndNumberReads = DecodeObjects()
		result_average = decodeResultAverageAndNumberReads.decode_result(list_meta[0].description)
		self.assertEqual('3000', result_average.number_file_1)
		self.assertEqual('1158.4', result_average.average_file_1)
		self.assertEqual(None, result_average.number_file_2)
		self.assertEqual(None, result_average.average_file_2)

		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
		self.assertTrue(len(list_meta) == 0)
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, None)
		self.assertEqual(2, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, list_meta[0].meta_tag.name)
		self.assertEquals("Success, NanoStat(1.4.0), NanoFilt(2.6.0)", list_meta[0].description)
		self.assertEquals("Success, NanoStat(1.4.0), NanoFilt(2.6.0)", list_meta[1].description)
		
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, None)
		self.assertEqual(2, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, list_meta[0].meta_tag.name)
		result_data = decodeResultAverageAndNumberReads.decode_result(list_meta[0].description)
		vect_soft = result_data.get_list_software_instance(self.software_names.get_NanoStat_name())
		self.assertEqual(2, len(vect_soft))

		self.assertEqual("Mean read length", vect_soft[0].get_vect_key_values()[0].key)
		self.assertEqual("1,298.4", vect_soft[0].get_vect_key_values()[0].value)
		self.assertEqual("Total bases", vect_soft[0].get_vect_key_values()[-1].key)
		self.assertEqual("3,895,144.0", vect_soft[0].get_vect_key_values()[-1].value)
		self.assertEqual("Mean read length", vect_soft[1].get_vect_key_values()[0].key)
		self.assertEqual("1,158.4", vect_soft[1].get_vect_key_values()[0].value)
		self.assertEqual("Total bases", vect_soft[1].get_vect_key_values()[-1].key)
		self.assertEqual("3,475,144.0", vect_soft[1].get_vect_key_values()[-1].value)
		
		stats_illumina, total_reads = self.software.get_stats_from_sample_reads(sample)
		self.assertEqual('Mean read length', stats_illumina[0][0])
		self.assertEqual('1,298.4', stats_illumina[0][1])
		self.assertEqual(3000, total_reads)
		
		sample = Sample.objects.get(name=sample_name)
		self.assertEqual("A-H5N5", sample.type_subtype)
	#	self.assertEqual("H5N5", sample.type_subtype)
		self.assertEqual("No", sample.mixed_infections_tag.name)
		self.assertEqual(0, sample.number_alerts)
		
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, None)
		self.assertEquals(2, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals("Success, NanoStat(1.4.0), NanoFilt(2.6.0)", list_meta[0].description)
		
		### software
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, None)
		self.assertEquals(2, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, list_meta[0].meta_tag.name)
		
		decode_result = DecodeObjects()
		result = decode_result.decode_result(list_meta[0].description)
		self.assertTrue(not result is None)
		self.assertEquals("NanoStat-1.4.0",
				result.get_software(self.software_names.get_NanoStat_name()))
		self.assertEquals("NanoFilt-2.6.0; (-q 5 -l 50 --headcrop 70 --tailcrop 70)",\
 						result.get_software(self.software_names.get_NanoFilt_name()))
		self.assertEquals("RabbitQC-0.0.1; (-w 3 -D)",\
 						result.get_software(self.software_names.get_rabbitQC_name()))
		
		## remove all files
		cmd = "rm -r %s*" % (temp_dir); os.system(cmd)