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
from utils.result import DecodeObjects, Coverage, CountHits
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference, Statistics
from manage_virus.uploadFiles import UploadFiles
from manage_virus.models import UploadFile
from managing_files.manage_database import ManageDatabase
from django.test.utils import override_settings
from utils.tree import CreateTree
import os, filecmp, csv
from utils.parse_in_files import ParseInFiles
from utils.result import DecodeObjects, MixedInfectionMainVector
from managing_files.models import CountVariations, MixedInfections
from utils.mixed_infections_management import MixedInfectionsManagement
from manage_virus.constants_virus import ConstantsVirus
from manage_virus.models import Tags, SeqVirus, IdentifyVirus
from settings.default_software import DefaultSoftware
from settings.models import Software as Software2, Parameter
from utils.parse_coverage_file import GetCoverage
from plotly.figure_factory._dendrogram import scs
from managing_files.models import Software as SoftwareModel

class Test(TestCase):

	### static
	software = Software()
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
	
	def test_run_fasttree(self):
		"""
 		test run mauve
 		create a VCF
 		"""
		in_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.FILE_OUT_MAFFT_RESULT)
		expect_file_nwk = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, ConstantsTestsCase.FILE_FASTTREE_RESULT_NWK)
		self.assertTrue(os.path.exists(expect_file_nwk))
		self.assertTrue(os.path.exists(in_file))
		
		out_file = self.utils.get_temp_file("fasttree", ".nwk")
		output_file = self.software.run_fasttree(in_file, out_file, self.software_names.get_fasttree_parameters())
		self.assertTrue(self.software_names.get_fasttree().endswith(self.software_names.get_fasttree_name()))
		self.assertEquals(output_file, out_file)
		self.assertTrue(filecmp.cmp(out_file, expect_file_nwk))
		
		### reroot the tree with fake leaf name		
		reroot_leaf = "EVA011_S54_"
		output_file = self.software.run_fasttree(in_file, out_file, self.software_names.get_fasttree_parameters(), reroot_leaf)
		self.assertTrue(filecmp.cmp(out_file, expect_file_nwk))
		
		### reroot the tree		
		reroot_leaf = "EVA011_S54"
		output_file = self.software.run_fasttree(in_file, out_file, self.software_names.get_fasttree_parameters(), reroot_leaf)
		expect_file_nwk = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "2_" + ConstantsTestsCase.FILE_FASTTREE_RESULT_NWK)
		self.assertTrue(filecmp.cmp(out_file, expect_file_nwk))
		os.unlink(out_file)
		
