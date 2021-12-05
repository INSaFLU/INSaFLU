'''
Created on Dec 16, 2017

@author: mmp
'''
import unittest, os, filecmp
from utils.utils import Utils 
from django.conf import settings 
from constants.constantsTestsCase import ConstantsTestsCase
from constants.constants import Constants, FileType
from utils.proteins import Proteins
from utils.result import Gene, GeneticElement
from Bio import SeqIO
from utils.result import Coverage
from django.test.utils import override_settings
from managing_files.manage_database import ManageDatabase
from constants.meta_key_and_values import MetaKeyAndValue
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference
from utils.result import DecodeObjects
from constants.software_names import SoftwareNames
from constants.constants import TypePath, FileExtensions

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
		gene = Gene('PB2', 0, 2269, 1)
		
		coverage = Coverage()
		coverage.add_coverage('PB2', Coverage.COVERAGE_MORE_9, 100.0)
		coverage.add_coverage('PB2', Coverage.COVERAGE_MORE_0, 100.0)
		coverage.add_coverage('PB2', Coverage.COVERAGE_ALL, 100.0)
		sequence_name = 'PB2'
		
		with open(consensus_fasta_file) as handle_consensus: 
			record_dict_consensus = SeqIO.to_dict(SeqIO.parse(handle_consensus, "fasta"))
		
		generic_element = GeneticElement()
		generic_element.add_gene("PB2", 10, gene)
		generic_element_consensus = GeneticElement()
		generic_element_consensus.add_gene("PB2", 10, gene)
	
		self.assertTrue(self.proteins.save_protein_by_sequence_name_and_cds(record_dict_consensus,
				generic_element_consensus, sample_name, sequence_name, gene, out_file))
				
		self.assertTrue(filecmp.cmp(out_file, expect_file))
		self.utils.remove_dir(out_dir)


	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_create_alignment_for_element(self):
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
			reference.display_name = ref_name
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
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina)
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
		self.assertEquals('FastTreeDbl-2.1.10 Double precision; (-gtr -boot 1000)', result.get_software(software_names.get_fasttree_name()))
#		self.assertEquals('Mafft-7.313; (--preservecase --amino)', result.get_software(software_names.get_mafft_name()))
		self.assertEquals(3, result.get_number_softwares())
		
		meta_key = metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_Tree_Count_Protein_By_Element, sequence_name)
		meta_sample = manageDatabase.get_project_metakey(project, meta_key, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_sample != None)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, meta_sample.value)
		self.assertEquals('4', meta_sample.description)
		self.assertEquals(99, os.path.getsize(project.get_global_file_by_element_and_cds(TypePath.MEDIA_ROOT, sequence_name, 'M', project.PROJECT_FILE_NAME_FASTTREE)))
		
		meta_key = metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_Tree_Count_Protein_By_Element, 'sequence_name')
		meta_sample = manageDatabase.get_project_metakey(project, meta_key, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_sample == None)
		
		### get all elements and gene names
		for type_files in project.vect_clean_file:
			if (type_files in project.vect_exclude_clean_file_from_proteins): continue
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
		

	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_create_alignement_for_element_2(self):
		"""
		Alignment_aa_SARS_CoV_2_ORF6.nex   SARSCoVDec200234
		Alignment_aa_SARS_CoV_2_ORF7a.nex  SARSCoVDec200234
		Alignment_aa_SARS_CoV_2_ORF8.nex   SARSCoVDec200153
 		"""
		temp_dir = self.utils.get_temp_dir()
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))
	
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_COVID_GBK)
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_COVID_FASTA)
		fasta_file_200153 = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "SARSCoVDec200153.consensus.fasta")
		file_expected_200153 = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "SARSCoVDec200234.protein_ORF8.fasta")
		
		self.assertTrue(os.path.exists(gb_file))
		self.assertTrue(os.path.exists(fasta_file))
		self.assertTrue(os.path.exists(fasta_file_200153))
		self.assertTrue(os.path.exists(file_expected_200153))
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME + '5110')
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME + '5110'
			user.id = 5110
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		ref_name = "second_stage_2_ test_create_tree_2"
		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.display_name = ref_name
			reference.reference_fasta.name = fasta_file
			reference.reference_fasta_name = os.path.basename(fasta_file)
			reference.reference_genbank.name = gb_file
			reference.reference_genbank_name = os.path.basename(gb_file)
			reference.owner = user
			reference.save()
		
		## get generic element
		sequence_name = 'MN908947'
		geneticElement = self.utils.get_elements_and_cds_from_db(reference, user)
		
		with open(fasta_file_200153) as handle_consensus: 
			record_dict_consensus = SeqIO.to_dict(SeqIO.parse(handle_consensus, "fasta"))
		
		coverage = Coverage()
		coverage.add_coverage(sequence_name, Coverage.COVERAGE_ALL, '2198.8')
		coverage.add_coverage(sequence_name, Coverage.COVERAGE_MORE_0, '100.0')
		coverage.add_coverage(sequence_name, Coverage.COVERAGE_MORE_9, '100.0')
			
		proteins_test = Proteins()
		limit_to_mask_consensus = -1
		genetic_element_from_sample = proteins_test.genetic_element_from_sample(fasta_file,
					record_dict_consensus, geneticElement, coverage, limit_to_mask_consensus, temp_dir)
		
		limit_to_mask_consensus = -1
		dt_out_files = {}
		dict_out_sample_name = {}
		n_count_samples_processed, n_with_sequences = 0, 0
		for gene in geneticElement.get_genes(sequence_name):	### can have more than one gene for each sequence
			### only do this one
			if gene.name != 'ORF8': continue
			
			## get file name
			if (gene.name not in dt_out_files): 
				dt_out_files[gene.name] = self.utils.get_temp_file_from_dir(temp_dir,\
						"{}_{}".format(sequence_name, self.utils.clean_name(gene.name)), FileExtensions.FILE_FAA)
			
			if (proteins_test.save_protein_by_sequence_name_and_cds(record_dict_consensus,
							genetic_element_from_sample,
							"xpto_sample", sequence_name, gene,
							dt_out_files[gene.name])):
				n_with_sequences += 1
			n_count_samples_processed += 1

			out_name = '{}_{}_{}'.format(reference.name.replace(' ', '_'), sequence_name, self.utils.clean_name(gene.name))
			dict_out_sample_name[out_name] = 1
			proteins_test.save_protein_reference_cds(gb_file,
							reference.display_name, sequence_name, gene,\
							dt_out_files[gene.name], dict_out_sample_name)
			### file with translated sequences
			print(dt_out_files[gene.name], file_expected_200153)
			self.assertTrue(filecmp.cmp(dt_out_files[gene.name], file_expected_200153))
			
		self.assertEqual(1, n_count_samples_processed)
		self.assertEqual(1, n_with_sequences)
		
		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))

	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_create_alignement_for_element_3(self):
		"""
		Alignment_aa_SARS_CoV_2_ORF6.nex   SARSCoVDec200234
		Alignment_aa_SARS_CoV_2_ORF7a.nex  SARSCoVDec200234
		Alignment_aa_SARS_CoV_2_ORF8.nex   SARSCoVDec200153
 		"""
		temp_dir = self.utils.get_temp_dir()
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))
	
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_COVID_GBK)
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_COVID_FASTA)
		fasta_file_200234 = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "SARSCoVDec200234.consensus.fasta")
		file_expected_200234_S = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "SARSCoVDec200153.protein_S.fasta")
		file_expected_200234_orf1ab = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, "SARSCoVDec200153.protein_orf1ab.fasta")
		
		self.assertTrue(os.path.exists(gb_file))
		self.assertTrue(os.path.exists(fasta_file))
		self.assertTrue(os.path.exists(fasta_file_200234))
		self.assertTrue(os.path.exists(file_expected_200234_S))
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME + '5110')
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME + '5110'
			user.id = 5110
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		ref_name = "second_stage_2_ test_create_tree_2"
		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.display_name = ref_name
			reference.reference_fasta.name = fasta_file
			reference.reference_fasta_name = os.path.basename(fasta_file)
			reference.reference_genbank.name = gb_file
			reference.reference_genbank_name = os.path.basename(gb_file)
			reference.owner = user
			reference.save()
		
		## get generic element
		sequence_name = 'MN908947'
		geneticElement = self.utils.get_elements_and_cds_from_db(reference, user)
		
		with open(fasta_file_200234) as handle_consensus: 
			record_dict_consensus = SeqIO.to_dict(SeqIO.parse(handle_consensus, "fasta"))
		
		coverage = Coverage()
		coverage.add_coverage(sequence_name, Coverage.COVERAGE_ALL, '2198.8')
		coverage.add_coverage(sequence_name, Coverage.COVERAGE_MORE_0, '100.0')
		coverage.add_coverage(sequence_name, Coverage.COVERAGE_MORE_9, '100.0')
		
		proteins_test = Proteins()
		limit_to_mask_consensus = -1
		genetic_element_from_sample = proteins_test.genetic_element_from_sample(fasta_file,
					record_dict_consensus, geneticElement, coverage, limit_to_mask_consensus, temp_dir)

		n_count_genes_compared = 0		
		dt_out_files = {}
		dict_out_sample_name = {}
		n_count_samples_processed, n_with_sequences = 0, 0
		for gene in geneticElement.get_genes(sequence_name):	### can have more than one gene for each sequence
			### only do this one
			if gene.name != 'orf1ab' and gene.name != 'ORF7a' and gene.name != 'ORF6' and gene.name != 'S': continue
			
			## get file name
			if (gene.name not in dt_out_files): 
				dt_out_files[gene.name] = self.utils.get_temp_file_from_dir(temp_dir,\
						"{}_{}".format(sequence_name, self.utils.clean_name(gene.name)), FileExtensions.FILE_FAA)
			
			n_count_samples_processed += 1
			if (proteins_test.save_protein_by_sequence_name_and_cds(record_dict_consensus,
							genetic_element_from_sample,
							"xpto_sample", sequence_name, gene,
							dt_out_files[gene.name])):
				n_with_sequences += 1
			else: continue

			out_name = '{}_{}_{}'.format(reference.name.replace(' ', '_'), sequence_name, self.utils.clean_name(gene.name))
			dict_out_sample_name[out_name] = 1
			proteins_test.save_protein_reference_cds(gb_file,
							reference.display_name, sequence_name, gene,\
							dt_out_files[gene.name], dict_out_sample_name)
			
			### file with translated sequences
			if (gene.name == 'S'):
				n_count_genes_compared += 1
				self.assertTrue(filecmp.cmp(dt_out_files[gene.name], file_expected_200234_S))
			if (gene.name == 'orf1ab'):
				self.assertTrue(filecmp.cmp(dt_out_files[gene.name], file_expected_200234_orf1ab))
				n_count_genes_compared += 1
		self.assertEqual(4, n_count_samples_processed)
		self.assertEqual(2, n_with_sequences)
		self.assertEqual(2, n_count_genes_compared)
		
		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))

		
	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_genetic_element_from_sample(self):
		"""
		get positions in the consensus file
		"""
		constantsTestsCase = ConstantsTestsCase()
		manageDatabase = ManageDatabase()
		proteins_test = Proteins()
		
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
			reference.display_name = ref_name
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
		vect_files = constantsTestsCase.get_all_fastq_files(self.baseDirectory)
		
		temp_dir = self.utils.get_temp_dir()
		count = 0
		b_at_least_one_less_than_100 = False
		for vect_file in vect_files:
			self.utils.copy_file(vect_file[0], os.path.join(temp_dir, os.path.basename(vect_file[0])))
			self.utils.copy_file(vect_file[1], os.path.join(temp_dir, os.path.basename(vect_file[1])))
				
			try:
				sample = Sample.objects.get(id=5100 + count + 1)
				sample_name = sample.name
			except Sample.DoesNotExist:
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
				sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina)
				sample.owner = user
				sample.is_ready_for_projects = True
				sample.is_obsolete = False
				sample.save()

			## create project_sample
			try:
				project_sample = ProjectSample.objects.get(id=51100 + count + 1)
			except ProjectSample.DoesNotExist:
				project_sample = ProjectSample()
				project_sample.id = 51100 + count + 1
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
		
		## copy all the data for getattr(settings, "MEDIA_ROOT", None)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))
		path_destination = os.path.join(getattr(settings, "MEDIA_ROOT", None), ConstantsTestsCase.DIR_PROJECTS)
		self.utils.make_path(os.path.join(getattr(settings, "MEDIA_ROOT", None), ConstantsTestsCase.DIR_PROJECTS))
		cmd = 'cp -r {}/{}/* {}'.format(self.baseDirectory, ConstantsTestsCase.DIR_PROJECTS, path_destination)
		os.system(cmd)

		limit_to_mask_consensus = 10
		for project_sample in project.project_samples.all():
			
			consensus_fasta = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA,
						SoftwareNames.SOFTWARE_SNIPPY_name \
						if project_sample.is_sample_illumina() else SoftwareNames.SOFTWARE_Medaka_name)
			
			### get dict consensus file
			with open(consensus_fasta) as handle_consensus: 
				record_dict_consensus = SeqIO.to_dict(SeqIO.parse(handle_consensus, "fasta"))
				if (record_dict_consensus is None): continue
				
			genetic_element_from_sample = proteins_test.genetic_element_from_sample(fasta_file,
					record_dict_consensus, geneticElement, coverage, limit_to_mask_consensus, temp_dir)
			
#			print(project_sample.sample.name)
			for sequence_name in geneticElement.get_sorted_elements():
				for _, gene in enumerate(geneticElement.get_genes(sequence_name)):
					print(project_sample.sample.name, gene.name)
					self.assertEquals(genetic_element_from_sample.get_genes(sequence_name)[_].start,
								geneticElement.get_genes(sequence_name)[_].start)
					self.assertEquals(genetic_element_from_sample.get_genes(sequence_name)[_].name,
								geneticElement.get_genes(sequence_name)[_].name)
					self.assertEquals(genetic_element_from_sample.get_genes(sequence_name)[_].end,
								geneticElement.get_genes(sequence_name)[_].end)
					self.assertEquals(genetic_element_from_sample.get_genes(sequence_name)[_].strand,
								geneticElement.get_genes(sequence_name)[_].strand)
# 			print(genetic_element_from_sample)
# 			print(geneticElement)

		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))



		