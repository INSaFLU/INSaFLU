'''
Created on Dec 15, 2017

@author: mmp
'''

import os
from constants.constants import TypePath, FileExtensions, Constants
from utils.utils import Utils
from utils.software import Software
from constants.software_names import SoftwareNames
from managing_files.manage_database import ManageDatabase
from constants.meta_key_and_values import MetaKeyAndValue
from settings.default_software_project_sample import DefaultProjectSoftware
from settings.default_parameters import DefaultParameters
from utils.result import GeneticElement, Gene, FeatureLocationSimple
from Bio import SeqIO
from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna
from utils.result import DecodeObjects, Result
from utils.result import SoftwareDesc
from settings.constants_settings import ConstantsSettings

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
		### has the alignments of the genes in the consensus
		self.dt_alignments_genes_consensus = {}
	
	
	def create_alignement_for_element(self, project, user, geneticElement, sequence_name):
		"""
		create an alignment for an element
		"""
						
		metaKeyAndValue = MetaKeyAndValue()
		manageDatabase = ManageDatabase()
		temp_dir = self.utils.get_temp_dir()
		
		### get limit mask
		default_software = DefaultProjectSoftware()
		
		### get meta_key
		meta_key = metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_Run_Proteins_Alignment_By_Element, sequence_name)
		
		dt_out_files = {}	### dictonary with files
		dict_out_sample_name = {}
		n_files_with_sequences = 0
		n_count_samples_processed = 0
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
			
			### get consensus
			if (project_sample.is_mask_consensus_sequences): 
				limit_to_mask_consensus = int(default_software.get_mask_consensus_single_parameter(project_sample,\
						DefaultParameters.MASK_CONSENSUS_threshold, ConstantsSettings.TECHNOLOGY_illumina \
						if project_sample.is_sample_illumina() else ConstantsSettings.TECHNOLOGY_minion))
			else: limit_to_mask_consensus = -1
			
			## test the coverage
			if (limit_to_mask_consensus == -1 and not coverage.is_100_more_9(sequence_name)) or\
				(limit_to_mask_consensus > 0 and not coverage.ratio_value_coverage_bigger_limit(sequence_name, limit_to_mask_consensus)):
				continue
			
			### get consensus file name
			consensus_fasta = project_sample.get_consensus_file(TypePath.MEDIA_ROOT)
			
			### get dict consensus file
			with open(consensus_fasta) as handle_consensus: 
				record_dict_consensus = SeqIO.to_dict(SeqIO.parse(handle_consensus, "fasta"))
				if (record_dict_consensus is None): continue
				
			### get positions of the genes in the consensus file
			if (not project_sample.id in self.dt_alignments_genes_consensus):
				self.dt_alignments_genes_consensus[project_sample.id] = \
					self.genetic_element_from_sample(project.reference.get_reference_fasta(TypePath.MEDIA_ROOT),
							record_dict_consensus, geneticElement, coverage, 
							limit_to_mask_consensus, temp_dir)
			
			b_first = True
			for gene in geneticElement.get_genes(sequence_name):	### can have more than one gene for each sequence
				## get file name
				if (gene.name not in dt_out_files): 
					dt_out_files[gene.name] = self.utils.get_temp_file_from_dir(temp_dir,\
							"{}_{}".format(sequence_name.replace('/', '_'), self.utils.clean_name(gene.name)), FileExtensions.FILE_FAA)
				
				if (self.save_protein_by_sequence_name_and_cds(record_dict_consensus,
							self.dt_alignments_genes_consensus[project_sample.id],
							project_sample.sample.name, sequence_name, gene,
							dt_out_files[gene.name])):
					if (b_first): n_files_with_sequences += 1
					dict_out_sample_name["{}_{}_{}".format(project_sample.sample.name, sequence_name, self.utils.clean_name(gene.name))] = 1
					b_first = False
			n_count_samples_processed += 1

		### error, there's no enough sequences to create tree file, at least one from sample + reference
		if (n_files_with_sequences < 1):
			manageDatabase.set_project_metakey(project, user, meta_key,\
					MetaKeyAndValue.META_VALUE_Error, "Error: there's no enough valid sequences to create a tree.")
			self.utils.remove_dir(temp_dir)
			## remove files that are going o be create
			self.clean_file_by_vect(project, sequence_name, geneticElement, project.vect_clean_file)
			return False
		
		### save the reference protein
		dt_reference_name_saved = {}
		for gene in geneticElement.get_genes(sequence_name):
			reference_name_saved = self.save_protein_reference_cds(project.reference.get_reference_gbk(TypePath.MEDIA_ROOT),
							project.reference.display_name, sequence_name, gene,\
							dt_out_files[gene.name], dict_out_sample_name)
			## keep the names of the reference saved in file
			dt_reference_name_saved[gene.name] = reference_name_saved
			
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
				self.software.run_fasttree(out_file_mafft, out_file_fasttree, self.software_names.get_fasttree_parameters_protein(),
										dt_reference_name_saved[gene_name])
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
			if (not path_file is None and os.path.exists(path_file)): os.unlink(path_file)

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
			reference_name_saved = '{}_{}_{}'.format(reference_name.replace(' ', '_'), sequence_name, self.utils.clean_name(gene.name))
			if (reference_name_saved in dict_out_sample_name): reference_name_saved += 'Ref_' + reference_name_saved
			handle.write('>{}\n{}\n'.format(reference_name_saved, str(coding_protein)))
		return reference_name_saved


	def save_protein_by_sequence_name_and_cds(self, record_dict_consensus, generic_element_consensus,
				sample_name, sequence_name, gene, out_file):
		"""
		save the protein sequence in a file
		out_file, append file
		"""
		
		### test if has genes
		if (generic_element_consensus.has_genes(sequence_name)):
		
			### get gene
			gene_to_translate_in_conensus = generic_element_consensus.get_gene(sequence_name, gene.name)
			if (gene_to_translate_in_conensus is None): return False
			
			gene_length = gene_to_translate_in_conensus.end - gene_to_translate_in_conensus.start
			seq_feature = gene_to_translate_in_conensus.get_seq_feature()
			if (not seq_feature is None):
				seq = Seq(str(record_dict_consensus[sequence_name].seq).replace('-', ''))
				sz_out = str(seq_feature.extract(seq))
			else:
				sz_out = str(record_dict_consensus[sequence_name].seq).replace('-', '')[gene_to_translate_in_conensus.start:
									gene_to_translate_in_conensus.end]
			
			### count N's if more than 20% discharge
			if (sz_out.count('N') / gene_length > 0.1): return False
			coding_dna = Seq(sz_out) ##, generic_dna)
			
			if (not gene_to_translate_in_conensus.is_forward()): coding_dna = coding_dna.reverse_complement()
			## this is not necessary because if seq is reversed the Bio.Seq is returned in forward always  
			## if not gene.is_forward(): coding_dna = coding_dna.reverse_complement()
			coding_protein = coding_dna.translate(table=Constants.TRANSLATE_TABLE_NUMBER, to_stop=False)
			count_stop = str(coding_protein)[:-1].count('*')
			
			### sometimes it's necessary to star in other frame to get the best translation
			count_stop_1 = 1000
			if (count_stop > 0):
				coding_dna_1 = Seq(sz_out[1:]) ##, generic_dna)
				coding_protein_1 = coding_dna_1.translate(table=Constants.TRANSLATE_TABLE_NUMBER, to_stop=False)
				count_stop_1 = str(coding_protein_1)[:-1].count('*')
				if (count_stop_1 < count_stop):
					coding_protein = coding_protein_1 

			if (count_stop > 0 and count_stop_1 > 0):
				coding_dna_2 = Seq(sz_out[2:]) ##, generic_dna)
				coding_protein_2 = coding_dna_2.translate(table=Constants.TRANSLATE_TABLE_NUMBER, to_stop=False)
				count_stop_2 = str(coding_protein_2)[:-1].count('*')
				if (count_stop > count_stop_2 and count_stop_1 > count_stop_2):
					coding_protein = coding_protein_2

			### can be zero
			if (len(str(coding_protein)) == 0): return False
			with open(out_file, 'a') as handle:
				handle.write('>{}\n{}\n'.format(sample_name.replace(' ', '_'), str(coding_protein)))
			return True
		return False
	
	def genetic_element_from_sample(self, reference_fasta_file, record_dict_consensus, genetic_element,
					coverage, limit_to_mask_consensus, temp_dir):
		"""
		get position where genes consensus from sample matches in the reference
		"""
		generic_consensus_element = GeneticElement()
		with open(reference_fasta_file) as handle_ref: 
			record_dict_ref = SeqIO.to_dict(SeqIO.parse(handle_ref, "fasta"))
			if (record_dict_ref is None): return generic_consensus_element

		for sequence_name in genetic_element.get_sorted_elements():
			
			if sequence_name in coverage.get_dict_data() and sequence_name in record_dict_ref and\
				sequence_name in record_dict_consensus and\
				( (limit_to_mask_consensus == -1 and coverage.is_100_more_9(sequence_name)) or\
				(limit_to_mask_consensus > 0 and coverage.ratio_value_coverage_bigger_limit(sequence_name, limit_to_mask_consensus)) ):
				
				### align two sequences
				seq_ref, seq_other = self.software.align_two_sequences(str(record_dict_ref[sequence_name].seq),
												str(record_dict_consensus[sequence_name].seq), temp_dir)
				
				#### get positions for genes				
				for gene in genetic_element.get_genes(sequence_name):
					### only has
					pos_ref = 0
					pos_con = 0
					cons_start = 0
					
					### can had some feature locations like join{[265:13468](+), [13467:21555](+)}
					vect_feature_location = []
					dt_positions = {}		## only for speed
					if (len(gene.get_feature_locations()) > 0):
						for feature in gene.get_feature_locations():
							dt_positions[feature.start] = 1
							dt_positions[feature.end] = 1
							vect_feature_location.append(FeatureLocationSimple(feature.start,
									feature.end, strand=feature.strand))
					else:
						dt_positions[gene.start] = 1
						dt_positions[gene.end] = 1
						vect_feature_location.append(FeatureLocationSimple(gene.start,
									gene.end, strand=gene.strand))
					
					###  let's start
					if (len(seq_ref) > 0 and len(seq_other) > 0):
						for i in range(0, len(seq_ref)):
							## has some position in seq features
							if (pos_ref in dt_positions):
								for _, feature in enumerate(gene.get_feature_locations()):
									if feature.start == pos_ref: vect_feature_location[_].start = pos_con
									if feature.end == pos_ref: vect_feature_location[_].end = pos_con

							### start and end position
							if (pos_ref == gene.start):
								cons_start = pos_con 
							if (pos_ref == gene.end):
								generic_consensus_element.add_gene(sequence_name,
										len(seq_other.replace('-', '')), Gene(
											gene.name, cons_start, pos_con,
											gene.strand,
											vect_feature_location))
								break
							
							if (seq_ref[i] != '-'): pos_ref += 1
							if (i < len(seq_other) and seq_other[i] != '-'): pos_con += 1
						
						### can have last pos
						if (pos_ref == gene.end):
							generic_consensus_element.add_gene(sequence_name,
									len(seq_other.replace('-', '')), Gene(
									gene.name, cons_start, pos_con,
									gene.strand, vect_feature_location))
		
		return generic_consensus_element

