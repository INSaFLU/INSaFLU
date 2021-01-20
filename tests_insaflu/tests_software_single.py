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

	def test_run_snippy(self):
		"""
		run snippy
		"""

		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_1, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
		self.utils.copy_file(file_2, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2))
			
		sample_name = "run_snippy1_sdfs"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = False
			sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
			sample.path_name_2.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2)
			sample.owner = user
			sample.save()

		project_name = "file_name_3_dsf"
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
		
		software_names = SoftwareNames()
		out_put_path = self.software.run_snippy(sample.get_fastq(TypePath.MEDIA_ROOT, True), sample.get_fastq(TypePath.MEDIA_ROOT, False),
				fasta_file, gb_file, sample.name, software_names.get_snippy_parameters())
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_BAM, "")))))
		file_size = os.path.getsize(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_BAM, ""))))
		self.assertTrue(file_size > 1000 )
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_CONSENSUS_FA, "")))))
		self.assertFalse(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_CONSENSUS_FASTA, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_CSV, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_DEPTH_GZ, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_DEPTH_GZ_TBI, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_TAB, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_VCF, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_VCF_GZ, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, "reference", os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_REF_FASTA, "")))))
		remove_path = os.path.dirname(out_put_path)
		if (len(remove_path.split('/')) > 2): self.utils.remove_dir(remove_path)
		else: self.utils.remove_dir(out_put_path)

		out_put_path = self.software.run_snippy(sample.get_fastq(TypePath.MEDIA_ROOT, True), None, fasta_file, gb_file, sample.name, software_names.get_snippy_parameters())
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_BAM, "")))))
		file_size_2 = os.path.getsize(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_BAM, ""))))
		self.assertTrue(file_size_2 > 1000)
		self.assertTrue(file_size_2 < file_size)
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_CONSENSUS_FA, "")))))
		self.assertFalse(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_CONSENSUS_FASTA, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_CSV, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_DEPTH_GZ, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_DEPTH_GZ_TBI, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_TAB, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_VCF, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_VCF_GZ, "")))))
		self.assertTrue(os.path.exists(os.path.join(out_put_path, "reference", os.path.basename(project_sample.get_file_output(TypePath.MEDIA_URL, FileType.FILE_REF_FASTA, "")))))
		remove_path = os.path.dirname(out_put_path)
		if (len(remove_path.split('/')) > 2): self.utils.remove_dir(remove_path)
		else: self.utils.remove_dir(out_put_path)
		self.utils.remove_dir(temp_dir)
