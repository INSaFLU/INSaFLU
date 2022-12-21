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
from managing_files.models import TagName
from managing_files.models import Software as SoftwareModel
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
	
	

		


