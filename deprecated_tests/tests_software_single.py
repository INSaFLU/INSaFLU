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
	#software = Software()
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
	
	
	def test_get_species_tag(self):
		"""
		test SARS cov 
		"""
		
		ref_name = "File 1"
		reference_fasta = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_COVID_FASTA)
		ref_name_2 = "File 2"
		consensus_file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "SARSCoVDec200153.consensus.fasta")
		ref_name_3 = "File 3"
		consensus_file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "A_H3N2_A_Hong_Kong_4801_2014.fasta")
		ref_name_4 = "File 4"
		consensus_file_3 = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "monkeypox_MT903344_MPXV_UK_P2_wo_NN.fasta")
		
		self.assertTrue(os.path.exists(reference_fasta))
		self.assertTrue(os.path.exists(consensus_file_1))
		self.assertTrue(os.path.exists(consensus_file_2))
		self.assertTrue(os.path.exists(consensus_file_3))
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.reference_fasta.name = reference_fasta
			reference.reference_fasta_name = os.path.basename(reference_fasta)
			reference.owner = user
			reference.save()
			
		try:
			reference1 = Reference.objects.get(name=ref_name_2)
		except Reference.DoesNotExist:
			reference1 = Reference()
			reference1.name = ref_name_2
			reference1.reference_fasta.name = consensus_file_1
			reference1.reference_fasta_name = os.path.basename(consensus_file_1)
			reference1.owner = user
			reference1.save()
		
		try:
			reference2 = Reference.objects.get(name=ref_name_3)
		except Reference.DoesNotExist:
			reference2 = Reference()
			reference2.name = ref_name_3
			reference2.reference_fasta.name = consensus_file_2
			reference2.reference_fasta_name = os.path.basename(consensus_file_2)
			reference2.owner = user
			reference2.save()
			
		try:
			reference3 = Reference.objects.get(name=ref_name_4)
		except Reference.DoesNotExist:
			reference3 = Reference()
			reference3.name = ref_name_3
			reference3.reference_fasta.name = consensus_file_3
			reference3.reference_fasta_name = os.path.basename(consensus_file_3)
			reference3.owner = user
			reference3.save()
		
		uploadFiles = UploadFiles()
		version = 1
		file = os.path.join(self.baseDirectory, Constants.DIR_TYPE_IDENTIFICATION, "test_covid_typing.fasta")
		self.assertTrue(os.path.exists(file))
		uploadFiles.upload_file(version, file)	## upload file
		
		os.unlink(out_file)
		os.unlink(out_file_consensus)

		self.assertEquals(Reference.SPECIES_INFLUENZA, self.software.get_species_tag(reference2))		
		self.assertEquals(Reference.SPECIES_MPXV, self.software.get_species_tag(reference3))
		self.assertEquals(Reference.SPECIES_SARS_COV_2, self.software.get_species_tag(reference))
		self.assertEquals(Reference.SPECIES_SARS_COV_2, self.software.get_species_tag(reference1))
		
		##
		try:
			reference_test = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			self.assertFail("Record must exist")
		
		self.assertEqual(Reference.SPECIES_SARS_COV_2, reference_test.specie_tag)