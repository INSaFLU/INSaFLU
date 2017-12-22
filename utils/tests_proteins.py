'''
Created on Dec 16, 2017

@author: mmp
'''
import unittest, os, filecmp
from utils.utils import Utils 
from django.conf import settings 
from constants.constantsTestsCase import ConstantsTestsCase
from utils.proteins import Proteins
from utils.result import Gene
from utils.result import Coverage
from django.test.utils import override_settings
from managing_files.manage_database import ManageDatabase
from constants.meta_key_and_values import MetaKeyAndValue
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference
from utils.result import DecodeObjects
from constants.software_names import SoftwareNames
from constants.constants import TypePath

class Test(unittest.TestCase):


	utils = Utils()
	proteins = Proteins()
	constants_tests_case = ConstantsTestsCase()
	
	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass
	
	def test_save_protein_by_sequence_name_and_cds(self):

		genbank_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		consensus_fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_CONSENSUS_ALIGNMENT_PROTEIN)
		expect_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_OUT_PROTEIN)
		self.assertTrue(os.path.join(genbank_file))

		out_dir = self.utils.get_temp_dir()
		out_file = self.utils.get_temp_file_from_dir(out_dir, 'test_save_protein_by_sequence_name_and_cds', '.faa')
		sample_name = 'xpto'
		gene = Gene('PB2', 10, 20, 1)
		
		coverage = Coverage()
		coverage.add_coverage('PB2', Coverage.COVERAGE_MORE_9, 100.0)
		coverage.add_coverage('PB2', Coverage.COVERAGE_MORE_0, 100.0)
		coverage.add_coverage('PB2', Coverage.COVERAGE_ALL, 100.0)
		sequence_name = 'PB2'
		self.assertTrue(self.proteins.save_protein_by_sequence_name_and_cds(consensus_fasta_file, genbank_file,
							sample_name, sequence_name, gene, coverage, out_dir, out_file))
		self.assertTrue(filecmp.cmp(out_file, expect_file))
		self.utils.remove_dir(out_dir)
		

	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_create_alignement_for_element(self):
		"""
 		test global method
 		"""
		manageDatabase = ManageDatabase()
		metaKeyAndValue = MetaKeyAndValue()
		
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))
	
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME + '5100')
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME + '5100'
			user.id = 5100
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		ref_name = "second_stage_2_ test_create_tree"
		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.reference_fasta.name = fasta_file
			reference.reference_fasta_name = os.path.basename(fasta_file)
			reference.reference_genbank.name = gb_file
			reference.reference_genbank_name = os.path.basename(gb_file)
			reference.owner = user
			reference.save()
		
		project_name = "several_names_test_proteines"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.id= 5100
			project.name = project_name
			project.reference = reference
			project.owner = user
			project.save()
		
		## get all fastq files
		vect_files = self.constants_tests_case.get_all_fastq_files(self.baseDirectory)
		
		temp_dir = self.utils.get_temp_dir()
		count = 0
		b_at_least_one_less_than_100 = False
		for vect_file in vect_files:
			self.utils.copy_file(vect_file[0], os.path.join(temp_dir, os.path.basename(vect_file[0])))
			self.utils.copy_file(vect_file[1], os.path.join(temp_dir, os.path.basename(vect_file[1])))
				
			sample_name = "_".join(os.path.basename(vect_file[0]).split('_')[0:2])
			sample = Sample()
			sample.id = 5100 + count + 1
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = vect_file[0]
			sample.path_name_1.name = os.path.join(temp_dir, os.path.basename(vect_file[0]))
			sample.is_valid_2 = False
			sample.file_name_2 = vect_file[1]
			sample.path_name_2.name = os.path.join(temp_dir, os.path.basename(vect_file[1]))
			sample.owner = user
			sample.is_ready_for_projects = True
			sample.is_obsolete = False
			sample.save()

			## create project_sample
			project_sample = ProjectSample()
			project_sample.id = 5100 + count + 1
			project_sample.sample = sample
			project_sample.project = project
			project_sample.is_finished = True
			project_sample.is_deleted = False
			project_sample.is_error = False
			project_sample.save()
			count += 1

			meta_value = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
			self.assertTrue(meta_value == None)
			
			coverage = Coverage()
			coverage.add_coverage('MP', Coverage.COVERAGE_ALL, '2198.8')
			coverage.add_coverage('MP', Coverage.COVERAGE_MORE_0, '100.0')
			coverage.add_coverage('MP', Coverage.COVERAGE_MORE_9, '100.0')
			coverage.add_coverage('PA', Coverage.COVERAGE_ALL, '527.8')
			coverage.add_coverage('PA', Coverage.COVERAGE_MORE_0, '100.0')
			coverage.add_coverage('PA', Coverage.COVERAGE_MORE_9, '100.0')
			coverage.add_coverage('HA', Coverage.COVERAGE_ALL, '1449.8')
			coverage.add_coverage('HA', Coverage.COVERAGE_MORE_0, '100.0')
			coverage.add_coverage('HA', Coverage.COVERAGE_MORE_9, '100.0')
			coverage.add_coverage('PB1', Coverage.COVERAGE_ALL, '618.8')
			coverage.add_coverage('PB1', Coverage.COVERAGE_MORE_0, '100.0')
			if (sample_name == 'EVA003_S91'):
				coverage.add_coverage('PB1', Coverage.COVERAGE_MORE_9, '90.0')
				b_at_least_one_less_than_100 = True
			else: coverage.add_coverage('PB1', Coverage.COVERAGE_MORE_9, '100.0')
			coverage.add_coverage('NP', Coverage.COVERAGE_ALL, '439.8')
			coverage.add_coverage('NP', Coverage.COVERAGE_MORE_0, '100.0')
			coverage.add_coverage('NP', Coverage.COVERAGE_MORE_9, '100.0')
			coverage.add_coverage('NS', Coverage.COVERAGE_ALL, '1214.8')
			coverage.add_coverage('NS', Coverage.COVERAGE_MORE_0, '100.0')
			coverage.add_coverage('NS', Coverage.COVERAGE_MORE_9, '100.0')
			coverage.add_coverage('PB2', Coverage.COVERAGE_ALL, '690.8')
			coverage.add_coverage('PB2', Coverage.COVERAGE_MORE_0, '100.0')
			coverage.add_coverage('PB2', Coverage.COVERAGE_MORE_9, '100.0')
			coverage.add_coverage('NA', Coverage.COVERAGE_ALL, '1092.8')
			coverage.add_coverage('NA', Coverage.COVERAGE_MORE_0, '100.0')
			coverage.add_coverage('NA', Coverage.COVERAGE_MORE_9, '100.0')
			meta_value = manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success, coverage.to_json())

		### test if it has the less 100 in EVA003_S91 sample
		self.assertTrue(b_at_least_one_less_than_100)
		self.assertEquals(len(vect_files), count)
		
		## get generic element
		geneticElement = self.utils.get_elements_and_cds_from_db(reference, user)
		sequence_name = 'MP'
		
		## copy all the data for getattr(settings, "MEDIA_ROOT", None)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))
		path_destination = os.path.join(getattr(settings, "MEDIA_ROOT", None), ConstantsTestsCase.DIR_PROJECTS)
		self.utils.make_path(os.path.join(getattr(settings, "MEDIA_ROOT", None), ConstantsTestsCase.DIR_PROJECTS))
		cmd = 'cp -r {}/{}/* {}'.format(self.baseDirectory, ConstantsTestsCase.DIR_PROJECTS, path_destination)
		os.system(cmd)
				
		proteins = Proteins()
		self.assertTrue(proteins.create_alignement_for_element(project, user, geneticElement, sequence_name))
		
		decode_result = DecodeObjects()
		software_names = SoftwareNames()
		meta_key = metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_Run_Proteins_Alignment_By_Element, sequence_name)
		meta_sample = manageDatabase.get_project_metakey(project, meta_key, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_sample != None)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, meta_sample.value)
		result = decode_result.decode_result(meta_sample.description)
		self.assertEquals('', result.get_software(software_names.get_mauve_name()))
		self.assertEquals('FastTree-2.1.10 SSE3; (-gtr -boot 1000)', result.get_software(software_names.get_fasttree_name()))
		self.assertEquals('Mafft-7.313; (--maxiterate 1000 --localpair --preservecase --amino)', result.get_software(software_names.get_mafft_name()))
		self.assertEquals(3, result.get_number_softwares())
		
		meta_key = metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_Tree_Count_Protein_By_Element, sequence_name)
		meta_sample = manageDatabase.get_project_metakey(project, meta_key, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_sample != None)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, meta_sample.value)
		self.assertEquals('4', meta_sample.description)
		self.assertTrue(os.path.getsize(project.get_global_file_by_element_and_cds(TypePath.MEDIA_ROOT, sequence_name, 'M', project.PROJECT_FILE_NAME_FASTTREE)) == 83)
		
		meta_key = metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_Tree_Count_Protein_By_Element, 'sequence_name')
		meta_sample = manageDatabase.get_project_metakey(project, meta_key, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_sample == None)
		
		### get all elements and gene names
		for type_files in project.vect_clean_file:
			self.assertTrue(os.path.exists(project.get_global_file_by_element_and_cds(TypePath.MEDIA_ROOT, sequence_name, 'M', type_files)))
		self.assertTrue(os.path.getsize(project.get_global_file_by_element_and_cds(TypePath.MEDIA_ROOT, sequence_name, 'M', project.PROJECT_FILE_NAME_FASTTREE)) > 80)

		expect_file = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_TREE, ConstantsTestsCase.MANAGING_TREE_out_protein)
		self.assertTrue(os.path.join(expect_file))
		temp_file = self.utils.get_temp_file("file_name", ".txt")
		temp_file_1 = self.utils.get_temp_file("file_name", ".txt")
		cmd = "grep -E -v 'TITLE: Written by EMBOSS' {} > {}".format(\
				project.get_global_file_by_element_and_cds(TypePath.MEDIA_ROOT, sequence_name, 'M', project.PROJECT_FILE_NAME_nex), temp_file)
		os.system(cmd);
		cmd = "grep -E -v 'TITLE: Written by EMBOSS' {} > {}".format(expect_file, temp_file_1)
		os.system(cmd);
		self.assertTrue(filecmp.cmp(temp_file_1, temp_file))
		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))
		
		
		
		
		
		