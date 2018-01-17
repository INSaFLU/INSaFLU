'''
Created on Dec 15, 2017

@author: mmp
'''

import os
from constants.constants import TypePath, FileType, FileExtensions, Constants
from utils.utils import Utils
from utils.software import Software
from constants.software_names import SoftwareNames
from managing_files.manage_database import ManageDatabase
from constants.meta_key_and_values import MetaKeyAndValue
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from utils.result import DecodeObjects, Result
from utils.result import SoftwareDesc

class Proteins(object):
	'''
	classdocs
	'''
	utils = Utils()
	software = Software()
	constants = Constants()
	software_names = SoftwareNames()

	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	
	def create_alignement_for_element(self, project, user, geneticElement, sequence_name):
		"""
		create an alignment for an element
		"""
						
		metaKeyAndValue = MetaKeyAndValue()
		manageDatabase = ManageDatabase()
		temp_dir = self.utils.get_temp_dir()
		
		### get meta_key
		meta_key = metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_Run_Proteins_Alignment_By_Element, sequence_name)
		
		dt_out_files = {}	### dictonary with files
		dict_out_sample_name = {}
		n_files_with_sequences = 0
		n_count_samples_processed = 0
		for project_sample in project.project_samples.all():
			if (not project_sample.get_is_ready_to_proccess()): continue
			### get coverage
			meta_value = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
			decode_coverage = DecodeObjects()
			coverage = decode_coverage.decode_result(meta_value.description)
			
			## test the coverage
			if (not coverage.is_100_more_9(sequence_name)): continue
			
			### get consensus file name
			consensus_fasta = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, SoftwareNames.SOFTWARE_SNIPPY_name)
			if (not os.path.exists(consensus_fasta)):
				manageDatabase.set_project_metakey(project, user, meta_key,\
						MetaKeyAndValue.META_VALUE_Error, "Error: fasta file doens't exist: " + consensus_fasta)
				self.utils.remove_dir(temp_dir)
				return False
			
			b_first = True
			for gene in geneticElement.get_genes(sequence_name):
				## get file name
				if (gene.name not in dt_out_files): 
					dt_out_files[gene.name] = self.utils.get_temp_file_from_dir(temp_dir,\
							"{}_{}".format(sequence_name, gene.name), FileExtensions.FILE_FAA)
				
				if (self.save_protein_by_sequence_name_and_cds(consensus_fasta, project.reference.get_reference_gbk(TypePath.MEDIA_ROOT),
							project_sample.sample.name, sequence_name, gene, coverage, temp_dir, dt_out_files[gene.name])):
					if (b_first): n_files_with_sequences += 1
					dict_out_sample_name["{}_{}_{}".format(project_sample.sample.name, sequence_name, gene.name)] = 1
					b_first = False
			n_count_samples_processed += 1

		### error, there's no enough sequences to create tree file
		if (n_files_with_sequences < 2):
			manageDatabase.set_project_metakey(project, user, meta_key,\
					MetaKeyAndValue.META_VALUE_Error, "Error: there's no enough valid sequences to create a tree.")
			self.utils.remove_dir(temp_dir)
			## remove files that are going o be create
			self.clean_file_by_vect(project, sequence_name, geneticElement, project.vect_clean_file)
			return False
		
		### save the reference protein
		for gene in geneticElement.get_genes(sequence_name):
			self.save_protein_reference_cds(project.reference.get_reference_gbk(TypePath.MEDIA_ROOT),
							project.reference.display_name, sequence_name, gene,\
							dt_out_files[gene.name], dict_out_sample_name)
			
			
		### start processing the data
		result_all = Result()
		b_first = True	
		## make gene by gene
		for gene_name in dt_out_files:
			
			### run mafft
			try:
				out_file_mafft = self.utils.get_temp_file_from_dir(temp_dir, "mafft", ".fasta")
				self.software.run_mafft(dt_out_files[gene_name], out_file_mafft, SoftwareNames.SOFTWARE_MAFFT_PARAMETERS_PROTEIN)
				if (b_first): result_all.add_software(SoftwareDesc(self.software_names.get_mafft_name(), self.software_names.get_mafft_version(),\
								SoftwareNames.SOFTWARE_MAFFT_PARAMETERS_PROTEIN))
			except Exception:
				result = Result()
				result.set_error("Mafft (%s) fail to run" % (self.software_names.get_mafft_version()))
				result.add_software(SoftwareDesc(self.software_names.get_mafft_name(), self.software_names.get_mafft_version(),\
								SoftwareNames.SOFTWARE_MAFFT_PARAMETERS_PROTEIN))
				manageDatabase.set_project_metakey(project, user, meta_key, MetaKeyAndValue.META_VALUE_Error, result.to_json())
				self.utils.remove_dir(temp_dir)
				return False
			
			### run seqret to produce nex
			try:
				out_file_nex = self.utils.get_temp_file_from_dir(temp_dir, "seq_ret_proteins", ".nex")
				self.software.run_seqret_nex(out_file_mafft, out_file_nex)
				if (b_first): result_all.add_software(SoftwareDesc(self.software_names.get_seqret_name(), self.software_names.get_seqret_version(),\
								self.software_names.get_seqret_nex_parameters()))
			except Exception:
				result = Result()
				result.set_error("{} {} fail to run".format(self.software_names.get_seqret_name(), self.software_names.get_seqret_version()))
				result.add_software(SoftwareDesc(self.software_names.get_seqret_name(), self.software_names.get_seqret_version(),\
								self.software_names.get_seqret_nex_parameters()))
				manageDatabase.set_project_metakey(project, user, meta_key, MetaKeyAndValue.META_VALUE_Error, result.to_json())
				self.utils.remove_dir(temp_dir)
				return False
	
			### run fastTree
			try:
				out_file_fasttree = self.utils.get_temp_file_from_dir(temp_dir, "fasttree_proteins", FileExtensions.FILE_NWK)
				self.software.run_fasttree(out_file_mafft, out_file_fasttree, self.software_names.get_fasttree_parameters_protein())
				if (b_first): result_all.add_software(SoftwareDesc(self.software_names.get_fasttree_name(), self.software_names.get_fasttree_version(),\
								self.software_names.get_fasttree_parameters_protein()))
			except Exception:
				result = Result()
				result.set_error("FastTree (%s) fail to run" % (self.software_names.get_fasttree_version()))
				result.add_software(SoftwareDesc(self.software_names.get_fasttree_name(), self.software_names.get_fasttree_version(),\
								self.software_names.get_fasttree_parameters_protein()))
				manageDatabase.set_project_metakey(project, user, meta_key, MetaKeyAndValue.META_VALUE_Error, result.to_json())
				self.utils.remove_dir(temp_dir)
				return False

			### copy files for the project
			self.utils.copy_file(out_file_mafft, project.get_global_file_by_element_and_cds(TypePath.MEDIA_ROOT, sequence_name, gene_name, project.PROJECT_FILE_NAME_MAFFT))
			self.utils.copy_file(out_file_fasttree, project.get_global_file_by_element_and_cds(TypePath.MEDIA_ROOT, sequence_name, gene_name, project.PROJECT_FILE_NAME_FASTTREE))
			self.utils.copy_file(out_file_fasttree, project.get_global_file_by_element_and_cds(TypePath.MEDIA_ROOT, sequence_name, gene_name, project.PROJECT_FILE_NAME_FASTTREE_tree))
			self.utils.copy_file(out_file_nex, project.get_global_file_by_element_and_cds(TypePath.MEDIA_ROOT, sequence_name, gene_name, project.PROJECT_FILE_NAME_nex))
			b_first = False
		
		### remove dir
		self.utils.remove_dir(temp_dir)
		
		## set meta info
		manageDatabase.set_project_metakey(project, user, meta_key, MetaKeyAndValue.META_VALUE_Success, result_all.to_json())
		meta_key = metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_Tree_Count_Protein_By_Element, sequence_name)
		manageDatabase.set_project_metakey(project, user, meta_key, MetaKeyAndValue.META_VALUE_Success, str(n_files_with_sequences))
		return True


	def clean_file_by_vect(self, project, sequence_name, genetic_element, vect_type_file):
		"""
		type_file: project.PROJECT_FILE_NAME_MAFFT, project.PROJECT_FILE_NAME_FASTTREE, project.PROJECT_FILE_NAME_FASTTREE_tree
		"""
		for type_file in vect_type_file: 
			self.clean_file(project, sequence_name, genetic_element, type_file)


	def clean_file(self, project, sequence_name, genetic_element, type_file):
		"""
		type_file: project.PROJECT_FILE_NAME_MAFFT, project.PROJECT_FILE_NAME_FASTTREE, project.PROJECT_FILE_NAME_FASTTREE_tree
		"""
		for gene in genetic_element.get_genes(sequence_name):
			if (type_file in project.vect_exclude_clean_file_from_proteins): continue
			path_file = project.get_global_file_by_element_and_cds(TypePath.MEDIA_ROOT, sequence_name, gene.name, type_file)
			if (path_file != None and os.path.exists(path_file)): os.unlink(path_file)

	def save_protein_reference_cds(self, genbank_file, reference_name, sequence_name, gene,\
							out_file, dict_out_sample_name):
		"""
		add the reference protein to the file
		"""
		### get protein reference sequence
		seq_ref = self.utils.get_sequence_from_genbank(sequence_name, gene, genbank_file)
		## this is not necessary because if seq is reversed the Bio.Seq is returned in forward always  
		## if not gene.is_forward(): seq_ref = seq_ref.reverse_complement()
		coding_protein = seq_ref.translate(table=Constants.TRANSLATE_TABLE_NUMBER, to_stop=False)
		with open(out_file, 'a') as handle:
			out_name = '{}_{}_{}'.format(reference_name.replace(' ', '_'), sequence_name, gene.name)
			if (out_name in dict_out_sample_name): out_name += 'Ref_' + out_name
			handle.write('>{}\n{}\n'.format(out_name, str(coding_protein)))

	def save_protein_by_sequence_name_and_cds(self, consensus_fasta_file, genbank_file,
							sample_name, sequence_name, gene, coverage, out_dir, out_file):
		"""
		save the protein sequence in a file
		out_file, append file
		"""
		### get reference sequence
		seq_ref = self.utils.get_sequence_from_genbank(sequence_name, gene, genbank_file)

		## there's no sequences...
		if (seq_ref == None): return False
		
		### get consensus sequence
		file_name = self.utils.filter_fasta_by_sequence_names(consensus_fasta_file, sample_name, sequence_name, coverage, out_dir)
		if (file_name == None): return False
		
		ref_seq_name = 'ref_seq_name__'
		with open(file_name, 'a') as handle:
			handle.write('>{}\n{}\n'.format(ref_seq_name, str(seq_ref)))
		
		## make the alignment
		out_file_clustalo = self.utils.get_temp_file_from_dir(out_dir, 'sequence_name_and_cds', '.fna')
		try:
			## It's better mafft to make the alignment
			self.software.run_mafft(file_name, out_file_clustalo, SoftwareNames.SOFTWARE_MAFFT_PARAMETERS_TWO_SEQUENCES)
		except Exception as a:
			return False
		
		## test if the output is in fasta
		try:
			self.utils.is_fasta(out_file_clustalo)
		except IOError as e:
			return False
		
		## read file
		record_dict = SeqIO.index(out_file_clustalo, "fasta")
		if (len(record_dict) != 2): return False
		
		seq_ref = ""
		seq_other = ""
		for seq in record_dict:
			if (seq == ref_seq_name):	## ref seq
				seq_ref = record_dict[seq].seq.upper()
			else:
				seq_other = record_dict[seq].seq.upper()	
		
		### we have sequences
		sz_out = ""
		sz_out_temp = ""
		b_start = False
		if (len(seq_ref) > 0 and len(seq_ref) == len(seq_other)):
			for i in range(0, len(seq_ref)):
				if (seq_ref[i] != '-'):
					b_start = True
					sz_out += sz_out_temp + seq_other[i]
					sz_out_temp = ''
				elif (b_start):
					sz_out_temp += seq_other[i]

			sz_out = sz_out.replace('-', '')
			coding_dna = Seq(sz_out, generic_dna)
			## this is not necessary because if seq is reversed the Bio.Seq is returned in forward always  
			## if not gene.is_forward(): coding_dna = coding_dna.reverse_complement()
			coding_protein = coding_dna.translate(table=Constants.TRANSLATE_TABLE_NUMBER, to_stop=False)
			with open(out_file, 'a') as handle:
				handle.write('>{}\n{}\n'.format(sample_name.replace(' ', '_'), str(coding_protein)))
			return True
		return False
		