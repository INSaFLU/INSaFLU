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

	def test_run_snpEff(self):
		"""
		test snpEff method
		"""
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		genbank_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		freebayes_vcf = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, ConstantsTestsCase.MANAGING_FILES_FREEBAYES_VCF)
		freebayes_expect_vcf = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_VCF, ConstantsTestsCase.MANAGING_FILES_FREEBAYES_ANNOTATED_VCF)
		
		out_file = self.utils.get_temp_file("file_name", ".vcf")
		out_file_2 = self.software.run_snpEff(fasta_file, genbank_file, freebayes_vcf, out_file)
		self.assertEquals(out_file, out_file_2)
		
		out_file_clean = self.utils.get_temp_file("file_name", ".vcf")
		cmd = "grep -v '{}' {} > {}".format(os.path.dirname(out_file), out_file, out_file_clean)
		os.system(cmd)
		print(out_file_clean)
		print(freebayes_expect_vcf)
		self.assertTrue(filecmp.cmp(out_file_clean, freebayes_expect_vcf))
		os.unlink(out_file_clean)
		os.unlink(out_file)