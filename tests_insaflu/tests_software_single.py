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

class Test(TestCase):

	### static
	software = Software()
	utils = Utils()
	constants = Constants()
	constants_tests_case = ConstantsTestsCase()
	software_names = SoftwareNames()

	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass


	def tearDown(self):
		pass
	
	def test_run_snippy_vcf_to_tab(self):
		"""
  		test snippy_vcf_to_tab method
  		"""
		
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		vcf_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "temp_more_REF.vcf")
		result_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, "resutl_vcf_to_tab.tab")

		#### add the transform p.Val423Glu to p.V423G
		### this run a on SNP_EFF
		parse_out_files = ParseOutFiles()
		vcf_file_with_pp = parse_out_files.add_amino_single_letter_code(vcf_file)
		
		out_file = self.utils.get_temp_file("snippy_vcf_to_tab", ".tab")
		out_file_2 = self.software.run_snippy_vcf_to_tab(fasta_file, gb_file, vcf_file_with_pp, out_file)
		self.assertEquals(out_file, out_file_2)
		print(out_file_2, result_file)
		self.assertTrue(filecmp.cmp(out_file_2, result_file))
		os.unlink(out_file)
		os.unlink(vcf_file_with_pp)
		
		
		
		
		