'''
Created on Nov 25, 2017

@author: mmp
'''
from utils.utils import Utils
from managing_files.manage_database import ManageDatabase
from constants.meta_key_and_values import MetaKeyAndValue
from utils.result import DecodeCoverage, Result
from constants.constants import TypePath, FileType, FileExtensions
from utils.software_names import SoftwareNames
from utils.result import SoftwareDesc
from utils.software import Software
import os

class CreateTree(object):
	'''
	classdocs
	'''

	utils = Utils()
	software_names = SoftwareNames()
	software = Software()
	
	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	def create_tree_and_alignments(self, project, owner):
		"""
		create both trees and the alignments
		"""
		
		### create tree and alignments for all genes
		self.create_tree_and_alignments_sample_by_sample(project, None, owner)
		
		### get all elements and gene names
		dict_genes = self.utils.get_elements_and_genes(project.reference.reference_genbank.name)
		
		### create for single sequences
		for sequence_name in dict_genes.keys():
			self.create_tree_and_alignments_sample_by_sample(project, sequence_name, owner)

	def create_tree_and_alignments_sample_by_sample(self, project, sequence_name, owner):
		"""
		create the tree and the alignments
		return: path to results, or None if some error
		"""
		metaKeyAndValue = MetaKeyAndValue()
		manageDatabase = ManageDatabase()
		temp_dir = self.utils.get_temp_dir()
		
		### get meta_key
		meta_key = MetaKeyAndValue.META_KEY_Run_Tree_All_Sequences\
				if sequence_name == None else metaKeyAndValue.get_meta_key_by_element(MetaKeyAndValue.META_KEY_Run_Tree_By_Element, sequence_name)
				
		n_files_with_sequences = 0
		n_count_samples_processed = 0
		### create a fasta file with all consensus that pass in filters for each project_sample
		dict_out_sample_name = {}
		for project_sample in project.project_sample.all():
			if (not project_sample.get_is_ready_to_proccess()): continue
			### get coverage
			meta_value = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
			decode_coverage = DecodeCoverage()
			coverage = decode_coverage.decode_result(meta_value.description)
			
			### get consensus
			consensus_fasta = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, SoftwareNames.SOFTWARE_SNIPPY_name)
			if (not os.path.exists(consensus_fasta)):
				manageDatabase.set_project_metakey(project, owner, meta_key,\
						MetaKeyAndValue.META_VALUE_Error, "Error: fasta file doens't exist: " + consensus_fasta)
				self.utils.remove_dir(temp_dir)
				return False
			if (sequence_name == None):
				if (self.utils.filter_fasta_all_sequences(consensus_fasta, project_sample.sample.name, coverage, temp_dir)):
					dict_out_sample_name[project_sample.sample.name] = 1
					n_files_with_sequences += 1
			elif (self.utils.filter_fasta_by_sequence_names(consensus_fasta, project_sample.sample.name, sequence_name, coverage, temp_dir)):
					dict_out_sample_name[project_sample.sample.name + "_" + sequence_name] = 1
					n_files_with_sequences += 1
			n_count_samples_processed += 1
		
		
		### error, there's no enough files to create consensus file
		if (n_files_with_sequences < 2):
			manageDatabase.set_project_metakey(project, owner, meta_key,\
					MetaKeyAndValue.META_VALUE_Error, "Error: there's no enough valid sequences to create the tree.")
			self.utils.remove_dir(temp_dir)
			return False
		
		### copy the reference
		if (sequence_name != None):
			## test if exist a sample that is equal
			sample_name = self.utils.clean_extension(project.reference.reference_fasta_name) + '_' + sequence_name
			if (sample_name in dict_out_sample_name): sample_name = 'Ref_' + sample_name 
			self.utils.filter_fasta_by_sequence_names(project.reference.reference_fasta.name,\
						sample_name, sequence_name, None, temp_dir)
		else:
			sample_name = self.utils.clean_extension(os.path.basename(project.reference.reference_fasta.name))
			if (sample_name in dict_out_sample_name): sample_name = 'Ref_' + sample_name 
			self.utils.copy_file(project.reference.reference_fasta.name, os.path.join(temp_dir, sample_name + FileExtensions.FILE_FASTA))
		n_files_with_sequences += 1
		
		### start processing the data
		result_all = Result()
		### run progressive mauve in all fasta files
		try:
			temp_file_to_remove = self.utils.get_temp_file("progressive_mauve", ".xfma")
			os.unlink(temp_file_to_remove)
			out_file_mauve = os.path.join(temp_dir, os.path.basename(temp_file_to_remove))	
			self.software.run_mauve(temp_dir, out_file_mauve)
			result_all.add_software(SoftwareDesc(self.software_names.get_mauve_name(), self.software_names.get_mauve_version(), self.software_names.get_mauve_parameters()))
		except Exception:
			result = Result()
			result.set_error("Mauve (%s) fail to run" % (self.software_names.get_mauve_version()))
			result.add_software(SoftwareDesc(self.software_names.get_mauve_name(), self.software_names.get_mauve_version(), self.software_names.get_mauve_parameters()))
			manageDatabase.set_project_metakey(project, owner, meta_key, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			self.utils.remove_dir(temp_dir)
			return False
		
		### run Convert.pl in mauve result
		try:
			temp_file_to_remove = self.utils.get_temp_file("convert_mauve", ".fasta")
			os.unlink(temp_file_to_remove)
			out_file_convert_mauve = os.path.join(temp_dir, os.path.basename(temp_file_to_remove))	
			out_file_convert_mauve = self.software.run_convert_mauve(out_file_mauve, out_file_convert_mauve)
			result_all.add_software(SoftwareDesc(self.software_names.get_convert_mauve_name(), self.software_names.get_convert_mauve_version(),\
							self.software_names.get_convert_mauve_parameters()))
		except Exception:
			result = Result()
			result.set_error("Convert mauve (%s) fail to run" % (self.software_names.get_convert_mauve_version()))
			result.add_software(SoftwareDesc(self.software_names.get_convert_mauve_name(), self.software_names.get_convert_mauve_version(),\
							self.software_names.get_convert_mauve_parameters()))
			manageDatabase.set_project_metakey(project, owner, meta_key, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			self.utils.remove_dir(temp_dir)
			return False
		
		### run mafft
		try:
			temp_file_to_remove = self.utils.get_temp_file("mafft", ".fasta")
			os.unlink(temp_file_to_remove)
			out_file_mafft = os.path.join(temp_dir, os.path.basename(temp_file_to_remove))		### keep this one
			self.software.run_mafft(out_file_convert_mauve, out_file_mafft)
			result_all.add_software(SoftwareDesc(self.software_names.get_mafft_name(), self.software_names.get_mafft_version(),\
							self.software_names.get_mafft_parameters()))
		except Exception:
			result = Result()
			result.set_error("MAfft (%s) fail to run" % (self.software_names.get_mafft_version()))
			result.add_software(SoftwareDesc(self.software_names.get_mafft_name(), self.software_names.get_mafft_version(),\
							self.software_names.get_mafft_parameters()))
			manageDatabase.set_project_metakey(project, owner, meta_key, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			self.utils.remove_dir(temp_dir)
			return False

		### run fastTree
		try:
			temp_file_to_remove = self.utils.get_temp_file("fasttree", FileExtensions.FILE_NWK)
			os.unlink(temp_file_to_remove)
			out_file_fasttree = os.path.join(temp_dir, os.path.basename(temp_file_to_remove))	### keep this file
			self.software.run_fasttree(out_file_mafft, out_file_fasttree)
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
			metaKeyAndValue.get_meta_key_by_element(MetaKeyAndValue.META_KEY_Tree_Count_By_Element, sequence_name)
		manageDatabase.set_project_metakey(project, owner, meta_key, MetaKeyAndValue.META_VALUE_Success, str(n_files_with_sequences))
		
		### copy files for the project
		if (sequence_name == None):
			self.utils.copy_file(out_file_mafft, project.get_global_file_by_project(TypePath.MEDIA_ROOT, project.PROJECT_FILE_NAME_MAFFT))
			self.utils.copy_file(out_file_fasttree, project.get_global_file_by_project(TypePath.MEDIA_ROOT, project.PROJECT_FILE_NAME_FASTTREE))
			self.utils.copy_file(out_file_fasttree, project.get_global_file_by_project(TypePath.MEDIA_ROOT, project.PROJECT_FILE_NAME_FASTTREE_tree))
		else:
			self.utils.copy_file(out_file_mafft, project.get_global_file_by_element(TypePath.MEDIA_ROOT, sequence_name, project.PROJECT_FILE_NAME_MAFFT))
			self.utils.copy_file(out_file_fasttree, project.get_global_file_by_element(TypePath.MEDIA_ROOT, sequence_name, project.PROJECT_FILE_NAME_FASTTREE))
			self.utils.copy_file(out_file_fasttree, project.get_global_file_by_element(TypePath.MEDIA_ROOT, sequence_name, project.PROJECT_FILE_NAME_FASTTREE_tree))
		self.utils.remove_dir(temp_dir)
		return True


