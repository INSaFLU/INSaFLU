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
	
	def test_run_pangolin(self):
		"""
		test run_pangolin method
		"""
		consensus_file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "SARSCoVDec200153.consensus.fasta")
		consensus_file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "SARSCoVDec200234.consensus.fasta")
		out_file_consensus = self.utils.get_temp_file("all_file_name", ".fasta")
		cmd = "cat {} {} > {}".format(consensus_file_1, consensus_file_2, out_file_consensus)
		os.system(cmd)
		
		try:
			software = SoftwareModel.objects.get(name=SoftwareNames.SOFTWARE_Pangolin_name)
			self.fail("Must not exist software name")
		except SoftwareModel.DoesNotExist:	## need to create with last version
			pass

		self.software_pangolin.run_pangolin_update()

		try:
			software = SoftwareModel.objects.get(name=SoftwareNames.SOFTWARE_Pangolin_name)
			software.is_updated_today()
			dt_software = software.get_version_long()
			self.assertTrue(len(dt_software) > 0)
			self.assertTrue(len(software.version) > 0)
		except SoftwareModel.DoesNotExist:	## need to create with last version
			self.fail("Must not exist software name")
		
		out_file = self.utils.get_temp_file("file_name", ".txt")
		self.software_pangolin.run_pangolin(out_file_consensus, out_file)
		
		vect_data = self.utils.read_text_file(out_file)
		self.assertEqual(3, len(vect_data))

		try:
			software = SoftwareModel.objects.get(name=SoftwareNames.SOFTWARE_Pangolin_name)
			dt_versions = software.get_version_long()
			self.assertTrue(len(dt_versions) > 0)
			self.assertTrue(len(software.version) > 0)
							
		except SoftwareModel.DoesNotExist:	## need to create with last version
			self.fail("Must exist software name")
		
		os.unlink(out_file)
		os.unlink(out_file_consensus)


