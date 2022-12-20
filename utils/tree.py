'''
Created on Nov 25, 2017

@author: mmp
'''
from utils.utils import Utils
from datasets.models import Dataset
from datasets.manage_database import ManageDatabase as ManageDatabaseDatasets
from managing_files.manage_database import ManageDatabase
from constants.meta_key_and_values import MetaKeyAndValue
from utils.result import DecodeObjects, Result
from constants.constants import TypePath, FileExtensions
from constants.software_names import SoftwareNames
from utils.result import SoftwareDesc
from utils.software import Software
from utils.proteins import Proteins
from constants.constants import Constants
from settings.default_software_project_sample import DefaultProjectSoftware
from settings.default_parameters import DefaultParameters
from settings.constants_settings import ConstantsSettings
from django.conf import settings
import os, logging, time


class CreateTree(object):
	'''
	classdocs
	'''

	utils = Utils()
	software_names = SoftwareNames()
	software = Software()
	
	if settings.DEBUG: logger = logging.getLogger("fluWebVirus.debug")
	else: logger = logging.getLogger("fluWebVirus.production")
	
	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	def create_tree_and_alignments(self, project, owner):
		"""
		create both trees and the alignments
		"""
		
		start= time.time()
		self.logger.info("START TREE and ALIGNEMTS:")
		
		### create tree and alignments for all genes
		self.create_tree_and_alignments_sample_by_sample(project, None, owner)
		self.logger.info("ENDE TRE and ALIGNEMTS: sequence_name {}  diff_time:{}".format("AllSequences", time.time() - start))
		start = time.time()
			
		proteins = Proteins()
		geneticElement = self.utils.get_elements_and_cds_from_db(project.reference, owner)
		### create for single sequences
		for sequence_name in geneticElement.get_sorted_elements():
			self.logger.info("MAKE TREE: sequence_name {}  diff_time:{}".format(sequence_name, time.time() - start))
			start = time.time()
			self.create_tree_and_alignments_sample_by_sample(project, sequence_name, owner)

			### create the protein alignments
			self.logger.info("MAKE ALIGNEMTS: sequence_name {}  diff_time:{}".format(sequence_name, time.time() - start))
			start = time.time()
			proteins.create_alignement_for_element(project, owner, geneticElement, sequence_name)
			
		self.logger.info("END TREE and ALIGNEMTS: diff_time:{}".format(time.time() - start))
		
	def create_tree_and_alignments_sample_by_sample(self, project, sequence_name, owner):
		"""
		create the tree and the alignments
		:param sequence_name Element name to process, None if they are all
		:out path to results, or None if some error
		"""
		metaKeyAndValue = MetaKeyAndValue()
		manageDatabase = ManageDatabase()
		temp_dir = self.utils.get_temp_dir()
		
		### get limit mask
		default_software = DefaultProjectSoftware()
		
		### get meta_key
		meta_key = MetaKeyAndValue.META_KEY_Run_Tree_All_Sequences\
				if sequence_name is None else metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_Run_Tree_By_Element, sequence_name)
		
#		n_max_elements_in_reference = 0 
		n_files_with_sequences = 0
		n_count_samples_processed = 0
		### create a fasta file with all consensus that pass in filters for each project_sample
		dict_out_sample_name = {}
		for project_sample in project.project_samples.all():
			if (not project_sample.get_is_ready_to_proccess()): continue
			if not os.path.exists(project_sample.get_consensus_file(TypePath.MEDIA_ROOT)): continue
			
			## test if it has to join in all consensus files
			if not default_software.include_consensus(project_sample): continue
			
			### get coverage
			meta_value = manageDatabase.get_project_sample_metakey_last(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
			if (meta_value is None): continue
			
			decode_coverage = DecodeObjects()
			coverage = decode_coverage.decode_result(meta_value.description)
			### get max number of elements
# 			if (len(coverage.dt_data) > n_max_elements_in_reference):
# 				n_max_elements_in_reference = len(coverage.dt_data)
			
			### get consensus
			if (project_sample.is_mask_consensus_sequences): 
				limit_to_mask_consensus = int(default_software.get_mask_consensus_single_parameter(project_sample,\
								DefaultParameters.MASK_CONSENSUS_threshold, ConstantsSettings.TECHNOLOGY_illumina \
								if project_sample.is_sample_illumina() else ConstantsSettings.TECHNOLOGY_minion))
			else: limit_to_mask_consensus = -1
			
			consensus_fasta = project_sample.get_consensus_file(TypePath.MEDIA_ROOT)
			if (sequence_name is None):		### join all elements
				if (self.utils.filter_fasta_all_sequences(consensus_fasta, project_sample.sample.name, coverage, limit_to_mask_consensus, temp_dir)):
					dict_out_sample_name[project_sample.sample.name] = 1
					n_files_with_sequences += 1
			elif (self.utils.filter_fasta_by_sequence_names(consensus_fasta, project_sample.sample.name, sequence_name, coverage, None,\
					limit_to_mask_consensus, temp_dir)):
				dict_out_sample_name[project_sample.sample.name + "_" + sequence_name] = 1
				n_files_with_sequences += 1
			n_count_samples_processed += 1
		
		### error, there's no enough files to create tree file
		if (n_files_with_sequences < Constants.MINIMUN_NUMER_SAMPLES_CACULATE_GLOBAL_FILES):
			manageDatabase.set_project_metakey(project, owner, meta_key,\
					MetaKeyAndValue.META_VALUE_Error, "Error: there's no enough valid sequences to create a tree.")
			self.utils.remove_dir(temp_dir)
			## remove files that are going to be created
			self.clean_file_by_vect(project, sequence_name, project.vect_clean_file)
			return False
		
		### copy the reference, get the Ref name
		reference_name = project.reference.display_name.replace(' ', '_')
		if (sequence_name is None):
			if (reference_name is None): reference_name = ""
			if (reference_name in dict_out_sample_name or len(reference_name) == 0): reference_name = 'Ref_' + reference_name 
			self.utils.copy_file(project.reference.get_reference_fasta(TypePath.MEDIA_ROOT), os.path.join(temp_dir, reference_name + FileExtensions.FILE_FASTA))
		else:
			## test if exist a sample that is equal
			if (reference_name in dict_out_sample_name or len(reference_name) == 0): reference_name = 'Ref_' + reference_name
			self.utils.filter_fasta_by_sequence_names(project.reference.get_reference_fasta(TypePath.MEDIA_ROOT),\
						reference_name, sequence_name, None, None, -1, temp_dir)
		n_files_with_sequences += 1

		##########################
		#####
		##### When we have only one element in the reference it's not necessary run progressive mauve
		### start processing the data
		result_all = Result()
		out_file_convert_mauve = self.utils.get_temp_file_from_dir(temp_dir, "convert_mauve", ".fasta")
		## need to join all sequences
		if (sequence_name is None):	## by element
			vect_elements = self.utils.get_elements_from_db(project.reference, project.owner) 
			self.utils.merge_fasta_and_join_sequences(temp_dir, vect_elements, out_file_convert_mauve)
		else: self.utils.merge_fasta_first_sequence(temp_dir, out_file_convert_mauve)	## altogether
			
# 		if (n_max_elements_in_reference == 1): ## only concatenate fasta files
# 			out_file_convert_mauve = self.utils.get_temp_file_from_dir(temp_dir, "convert_mauve", ".fasta")
# 			self.utils.merge_fasta_first_sequence(temp_dir, out_file_convert_mauve)
# 			
# 		else:	## because we have more than one element in reference we need to make a rough alignment first 
# 			### run progressive mauve in all fasta files
# 			try:
# 				out_file_mauve = self.utils.get_temp_file_from_dir(temp_dir, "progressive_mauve", ".xfma")
# 				self.software.run_mauve(temp_dir, out_file_mauve)
# 				result_all.add_software(SoftwareDesc(self.software_names.get_mauve_name(), self.software_names.get_mauve_version(), self.software_names.get_mauve_parameters()))
# 			except Exception:
# 				result = Result()
# 				result.set_error("Mauve (%s) fail to run" % (self.software_names.get_mauve_version()))
# 				result.add_software(SoftwareDesc(self.software_names.get_mauve_name(), self.software_names.get_mauve_version(), self.software_names.get_mauve_parameters()))
# 				manageDatabase.set_project_metakey(project, owner, meta_key, MetaKeyAndValue.META_VALUE_Error, result.to_json())
# 				self.utils.remove_dir(temp_dir)
# 				return False
# 			
# 			### run Convert.pl in mauve result
# 			try:
# 				out_file_convert_mauve = self.utils.get_temp_file_from_dir(temp_dir, "convert_mauve", ".fasta")
# 				out_file_convert_mauve = self.software.run_convert_mauve(out_file_mauve, out_file_convert_mauve)
# 				result_all.add_software(SoftwareDesc(self.software_names.get_convert_mauve_name(), self.software_names.get_convert_mauve_version(),\
# 								self.software_names.get_convert_mauve_parameters()))
# 			except Exception:
# 				result = Result()
# 				result.set_error("Convert mauve (%s) fail to run" % (self.software_names.get_convert_mauve_version()))
# 				result.add_software(SoftwareDesc(self.software_names.get_convert_mauve_name(), self.software_names.get_convert_mauve_version(),\
# 								self.software_names.get_convert_mauve_parameters()))
# 				manageDatabase.set_project_metakey(project, owner, meta_key, MetaKeyAndValue.META_VALUE_Error, result.to_json())
# 				self.utils.remove_dir(temp_dir)
# 				return False
		
		### run mafft
		try:
			out_file_mafft = self.utils.get_temp_file_from_dir(temp_dir, "mafft", ".fasta")
			self.software.run_mafft(out_file_convert_mauve, out_file_mafft, SoftwareNames.SOFTWARE_MAFFT_PARAMETERS)
			result_all.add_software(SoftwareDesc(self.software_names.get_mafft_name(), self.software_names.get_mafft_version(),\
							self.software_names.get_mafft_parameters()))
		except Exception as e:
			result = Result()
			result.set_error("Mafft (%s) fail to run" % (self.software_names.get_mafft_version()))
			result.add_software(SoftwareDesc(self.software_names.get_mafft_name(), self.software_names.get_mafft_version(),\
							self.software_names.get_mafft_parameters()))
			manageDatabase.set_project_metakey(project, owner, meta_key, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			self.utils.remove_dir(temp_dir)
			return False
		
		### run seqret to produce nex
		try:
			out_file_nex = self.utils.get_temp_file_from_dir(temp_dir, "seq_ret", ".nex")
			self.software.run_seqret_nex(out_file_mafft, out_file_nex)
			result_all.add_software(SoftwareDesc(self.software_names.get_seqret_name(), self.software_names.get_seqret_version(),\
							self.software_names.get_seqret_nex_parameters()))
		except Exception:
			result = Result()
			result.set_error("{} {} fail to run".format(self.software_names.get_seqret_name(), self.software_names.get_seqret_version()))
			result.add_software(SoftwareDesc(self.software_names.get_seqret_name(), self.software_names.get_seqret_version(),\
							self.software_names.get_seqret_nex_parameters()))
			manageDatabase.set_project_metakey(project, owner, meta_key, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			self.utils.remove_dir(temp_dir)
			return False

		### run fastTree
		try:
			## dvtditr
			out_file_fasttree = self.utils.get_temp_file_from_dir(temp_dir, "fasttree", FileExtensions.FILE_NWK)
			self.software.run_fasttree(out_file_mafft, out_file_fasttree, self.software_names.get_fasttree_parameters(),
								reference_name)
			result_all.add_software(SoftwareDesc(self.software_names.get_fasttree_name(), self.software_names.get_fasttree_version(),\
							self.software_names.get_fasttree_parameters()))
		except Exception:
			result = Result()
			result.set_error("FastTree (%s) fail to run" % (self.software_names.get_fasttree_version()))
			result.add_software(SoftwareDesc(self.software_names.get_fasttree_name(), self.software_names.get_fasttree_version(),\
							self.software_names.get_fasttree_parameters()))
			manageDatabase.set_project_metakey(project, owner, meta_key, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			self.utils.remove_dir(temp_dir)
			return False
		
		## set meta info
		manageDatabase.set_project_metakey(project, owner, meta_key, MetaKeyAndValue.META_VALUE_Success, result_all.to_json())
		meta_key = MetaKeyAndValue.META_KEY_Tree_Count_All_Sequences if sequence_name == None else\
			metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_Tree_Count_By_Element, sequence_name)
		manageDatabase.set_project_metakey(project, owner, meta_key, MetaKeyAndValue.META_VALUE_Success, str(n_files_with_sequences))
		
		## clean fasta file from alignment
		out_clean_fasta_file = self.utils.get_temp_file_from_dir(temp_dir, "clean", ".fasta")
		self.utils.clean_fasta_file(out_file_mafft, out_clean_fasta_file)
		
		### copy files for the project
		if (sequence_name == None):
			self.utils.copy_file(out_file_mafft, project.get_global_file_by_project(TypePath.MEDIA_ROOT, project.PROJECT_FILE_NAME_MAFFT))
			self.utils.copy_file(out_clean_fasta_file, project.get_global_file_by_project(TypePath.MEDIA_ROOT, project.PROJECT_FILE_NAME_FASTA))
			self.utils.copy_file(out_file_fasttree, project.get_global_file_by_project(TypePath.MEDIA_ROOT, project.PROJECT_FILE_NAME_FASTTREE))
			self.utils.copy_file(out_file_fasttree, project.get_global_file_by_project(TypePath.MEDIA_ROOT, project.PROJECT_FILE_NAME_FASTTREE_tree))
			self.utils.copy_file(out_file_nex, project.get_global_file_by_project(TypePath.MEDIA_ROOT, project.PROJECT_FILE_NAME_nex))
		else:
			self.utils.copy_file(out_file_mafft, project.get_global_file_by_element(TypePath.MEDIA_ROOT, sequence_name, project.PROJECT_FILE_NAME_MAFFT))
			self.utils.copy_file(out_clean_fasta_file, project.get_global_file_by_element(TypePath.MEDIA_ROOT, sequence_name, project.PROJECT_FILE_NAME_FASTA))
			self.utils.copy_file(out_file_fasttree, project.get_global_file_by_element(TypePath.MEDIA_ROOT, sequence_name, project.PROJECT_FILE_NAME_FASTTREE))
			self.utils.copy_file(out_file_fasttree, project.get_global_file_by_element(TypePath.MEDIA_ROOT, sequence_name, project.PROJECT_FILE_NAME_FASTTREE_tree))
			self.utils.copy_file(out_file_nex, project.get_global_file_by_element(TypePath.MEDIA_ROOT, sequence_name, project.PROJECT_FILE_NAME_nex))
		self.utils.remove_dir(temp_dir)
		return True


	def clean_file_by_vect(self, project, sequence_name, vect_type_file):
		"""
		type_file: project.PROJECT_FILE_NAME_MAFFT, project.PROJECT_FILE_NAME_FASTTREE, project.PROJECT_FILE_NAME_FASTTREE_tree
		"""
		for type_file in vect_type_file: 
			self.clean_file(project, sequence_name, type_file)


	def clean_file(self, project, sequence_name, type_file):
		"""
		type_file: project.PROJECT_FILE_NAME_MAFFT, project.PROJECT_FILE_NAME_FASTTREE, project.PROJECT_FILE_NAME_FASTTREE_tree
		"""
		if (sequence_name == None): path_file = project.get_global_file_by_project(TypePath.MEDIA_ROOT, type_file)
		else: path_file = project.get_global_file_by_element(TypePath.MEDIA_ROOT, sequence_name, type_file)
		if (os.path.exists(path_file)): os.unlink(path_file)
		

	def create_tree_and_alignments_dataset(self, dataset, owner):
		"""
		create both trees and the alignments
		"""
		
		start= time.time()
		self.logger.info("START TREE and ALIGNEMTS:")
		
		manageDatabase = ManageDatabaseDatasets()
		temp_dir = self.utils.get_temp_dir()
		
		### get meta_key
		meta_key = MetaKeyAndValue.META_KEY_Run_Tree_All_Sequences
		
		##########################
		#####
		##### When we have only one element in the reference it's not necessary run progressive mauve
		### start processing the data
		result_all = Result()
		out_file_convert_mauve = dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_RESULT_all_consensus)
		
		### create tree and alignments for all genes
		### run mafft
		try:
			out_file_mafft = self.utils.get_temp_file_from_dir(temp_dir, "mafft", ".fasta")
			self.software.run_mafft(out_file_convert_mauve, out_file_mafft, SoftwareNames.SOFTWARE_MAFFT_PARAMETERS)
			result_all.add_software(SoftwareDesc(self.software_names.get_mafft_name(), self.software_names.get_mafft_version(),\
							self.software_names.get_mafft_parameters()))
		except Exception as e:
			result = Result()
			result.set_error("Mafft (%s) fail to run" % (self.software_names.get_mafft_version()))
			result.add_software(SoftwareDesc(self.software_names.get_mafft_name(), self.software_names.get_mafft_version(),\
							self.software_names.get_mafft_parameters()))
			manageDatabase.set_dataset_metakey(dataset, owner, meta_key, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			self.utils.remove_dir(temp_dir)
			return False
		
		### run seqret to produce nex
		try:
			out_file_nex = self.utils.get_temp_file_from_dir(temp_dir, "seq_ret", ".nex")
			self.software.run_seqret_nex(out_file_mafft, out_file_nex)
			result_all.add_software(SoftwareDesc(self.software_names.get_seqret_name(), self.software_names.get_seqret_version(),\
							self.software_names.get_seqret_nex_parameters()))
		except Exception:
			result = Result()
			result.set_error("{} {} fail to run".format(self.software_names.get_seqret_name(), self.software_names.get_seqret_version()))
			result.add_software(SoftwareDesc(self.software_names.get_seqret_name(), self.software_names.get_seqret_version(),\
							self.software_names.get_seqret_nex_parameters()))
			manageDatabase.set_dataset_metakey(dataset, owner, meta_key, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			self.utils.remove_dir(temp_dir)
			return False

		### run fastTree
		try:
			## dvtditr
			out_file_fasttree = self.utils.get_temp_file_from_dir(temp_dir, "fasttree", FileExtensions.FILE_NWK)
			self.software.run_fasttree(out_file_mafft, out_file_fasttree, self.software_names.get_fasttree_parameters(),
								dataset.get_first_reference_name())
			result_all.add_software(SoftwareDesc(self.software_names.get_fasttree_name(), self.software_names.get_fasttree_version(),\
							self.software_names.get_fasttree_parameters()))
		except Exception:
			result = Result()
			result.set_error("FastTree (%s) fail to run" % (self.software_names.get_fasttree_version()))
			result.add_software(SoftwareDesc(self.software_names.get_fasttree_name(), self.software_names.get_fasttree_version(),\
							self.software_names.get_fasttree_parameters()))
			manageDatabase.set_dataset_metakey(dataset, owner, meta_key, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			self.utils.remove_dir(temp_dir)
			return False
		
		## set meta info
		manageDatabase.set_dataset_metakey(dataset, owner, meta_key, MetaKeyAndValue.META_VALUE_Success, result_all.to_json())
		
		### copy files for the dataset
		self.utils.copy_file(out_file_mafft, dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_MAFFT))
		self.utils.copy_file(out_file_nex, dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_nex))
		self.utils.copy_file(out_file_fasttree, dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_FASTTREE))
		self.utils.copy_file(out_file_fasttree, dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_FASTTREE_tree))
		self.utils.remove_dir(temp_dir)
		return True


		self.logger.info("ENDE TREE and ALIGNEMTS: sequence_name {}  diff_time:{}".format("AllSequences", time.time() - start))
		start = time.time()
		