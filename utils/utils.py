'''
Created on Oct 31, 2017

@author: mmp
'''
import glob
import gzip
import hashlib
import logging
import ntpath
import os
import random
import re
import stat

from Bio import SeqIO
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.Seq import MutableSeq, Seq
from Bio.SeqFeature import CompoundLocation
from Bio.SeqRecord import SeqRecord
from constants.constants import Constants, FileExtensions, TypeFile, TypePath
from constants.meta_key_and_values import MetaKeyAndValue
from constants.software_names import SoftwareNames
from datasets.models import DatasetConsensus
from managing_files.manage_database import ManageDatabase
from managing_files.models import ProjectSample

## from Bio.Alphabet import IUPAC    version 1.78 doesn't have Bio.Alphabet 
from utils.result import FeatureLocationSimple, Gene, GeneticElement

## Add 'Ter' to dictonary
## http://www.hgmd.cf.ac.uk/docs/cd_amino.html
protein_letters_3to1['Ter'] = 'X'
from datetime import datetime
from statistics import mean

from django.conf import settings
from django.core.mail import send_mail
from django.db import transaction
from django.utils.translation import ugettext_lazy as _
from pysam import pysam

from utils.result import CountHits, DecodeObjects, MaskingConsensus


class Utils(object):
	'''
	class docs
	'''

	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")

	def __init__(self):
		'''
		Constructor
		'''
		pass

	def get_path_to_reference_file(self, user_id, ref_id):
		"""
		get the path to reference
		"""
		return os.path.join(Constants.DIR_PROCESSED_FILES_REFERENCE, "userId_{0}".format(user_id), "refId_{0}".format(ref_id))

	def get_path_to_consensus_file(self, user_id, ref_id):
		"""
		get the path to reference
		"""
		return os.path.join(Constants.DIR_PROCESSED_FILES_CONSENSUS, "userId_{0}".format(user_id), "consensusId_{0}".format(ref_id))
	
	def get_path_to_fastq_file(self, user_id, sample_id):
		"""
		get the path to sample
		"""
		return os.path.join(Constants.DIR_PROCESSED_FILES_FASTQ, "userId_{0}".format(user_id), "sampleId_{0}".format(sample_id))

	def get_sample_list_by_user(self, user_id, type_path, extension):
		"""
		get the path to sample
		"""
		return os.path.join(getattr(settings, type_path, None),
					Constants.DIR_PROCESSED_FILES_FASTQ, "userId_{0}".format(user_id),
					Constants.SAMPLE_LIST_all_samples + extension)
		
	def get_project_list_by_user(self, user_id, type_path, extension):
		"""
		get the path to sample
		"""
		return os.path.join(getattr(settings, type_path, None),
					Constants.DIR_PROCESSED_FILES_PROJECT, "user_{0}".format(user_id),
					Constants.PROJECTS_LIST_all_samples + extension)
	
	def get_path_upload_file(self, user_id, type_file):
		"""
		user_id ->
		type_file -> TypeFile.TYPE_FILE_sample_file, TypeFile.TYPE_FILE_fastq_gz
		"""
		file_path = 'csv_sample_file'
		if type_file == TypeFile.TYPE_FILE_fastq_gz: file_path = 'fastq_files'
		if type_file == TypeFile.TYPE_FILE_dataset_file_metadata: file_path = 'tsv_dataset_file'
		return os.path.join(Constants.DIR_PROCESSED_FILES_MULTIPLE_SAMPLES, "userId_{0}".format(user_id),\
				"{}".format(file_path))
		
	def get_unique_file(self, file_name):
		"""
		get unique file name from a file_name
		return '<path file_name>/<file_name>'
		OR if exists
		return '<path file_name>/<random number>/<file_name>' path_added
		path_added Can be None is the file does not exist
		"""
		temp_file_name = ntpath.basename(file_name.replace(" ", "_"))
		main_path = os.path.dirname(file_name)
		if (not os.path.exists(main_path)): os.makedirs(main_path, exist_ok=True)
		path_added = None
		while 1:
			if (not os.path.exists(os.path.join(main_path, temp_file_name))): break
			path_added = str(random.randrange(10000000, 99999999, 10))
			temp_file_name = os.path.join(path_added, ntpath.basename(file_name))
		return os.path.join(main_path, temp_file_name.replace(" ", "_")), path_added

	def get_temp_file(self, file_name, sz_type):
		"""
		return a temp file name
		"""
		main_path = os.path.join(Constants.TEMP_DIRECTORY, Constants.COUNT_DNA_TEMP_DIRECTORY)
		if (not os.path.exists(main_path)): os.makedirs(main_path, exist_ok=True)
		self.touch_file(main_path)
		while 1:
			return_file = os.path.join(main_path, "insa_flu_file_" + file_name + "_" + str(random.randrange(10000000, 99999999, 10)) + "_file" + sz_type)
			if (os.path.exists(return_file)): continue
			try:
				os.close(os.open(return_file, os.O_CREAT | os.O_EXCL))
				return return_file
			except OSError:
				pass
			
	def get_temp_file_from_dir(self, dir_out, file_name, sz_type):
		"""
		return a temp file name
		"""
		if (not os.path.exists(dir_out)): os.makedirs(dir_out, exist_ok=True)
		self.touch_file(dir_out)		## up to date the path
		while 1:
			return_file = os.path.join(dir_out, "insa_flu_" + file_name + "_" + str(random.randrange(10000000, 99999999, 10)) + "_file" + sz_type)
			if (os.path.exists(return_file)): continue
			try:
				os.close(os.open(return_file, os.O_CREAT | os.O_EXCL))
				return return_file
			except OSError:
				pass

	def touch_file(self, file_name):
		"""
		Uptodate dir name/file to not be removed by file system 
		"""
		if (os.path.exists(file_name)):
			cmd = "touch {}".format(file_name)
			os.system(cmd)

	def get_temp_dir(self):
		"""
		return a temp directory
		"""
		main_path = os.path.join(Constants.TEMP_DIRECTORY, Constants.COUNT_DNA_TEMP_DIRECTORY)
		if (not os.path.exists(main_path)): os.makedirs(main_path, exist_ok=True)
		self.touch_file(main_path)		## up to date main path to not be removed by file system
		while 1:
			return_path = os.path.join(main_path, "influ_path_{}_{}".format(
				str(os.getpid()), str(random.randrange(10000000, 99999999, 10))) )
			if (not os.path.exists(return_path)):
				os.makedirs(return_path, exist_ok=True)
				return return_path

	def get_file_name_without_extension(self, file_name):
		"""
		return file name without extension
		"""
		return os.path.splitext(os.path.basename(file_name))[0]
		
		
	def remove_temp_file(self, sz_file_name):
		"""
		prevent to remove files outside of temp directory
		"""
		if (sz_file_name == None): return
		
		if os.path.exists(sz_file_name) and len(sz_file_name) > 0 and sz_file_name.startswith(Constants.TEMP_DIRECTORY):
			cmd = "rm " + sz_file_name
			exist_status = os.system(cmd)
			if (exist_status != 0):
				self.logger_production.error('Fail to run: ' + cmd)
				self.logger_debug.error('Fail to run: ' + cmd)
				raise Exception("Fail to remove a file") 

	def remove_file(self, sz_file_name):
		"""
		Remove files
		return True if the file exists and was removed
		"""
		if (sz_file_name == None): return False
		
		if os.path.exists(sz_file_name) and len(sz_file_name) > 0:
			cmd = "rm " + sz_file_name
			exist_status = os.system(cmd)
			if (exist_status != 0):
				self.logger_production.error('Fail to run: ' + cmd)
				self.logger_debug.error('Fail to run: ' + cmd)
				raise Exception("Fail to remove a file")
			return True
		return False

	def remove_dir(self, path_name):
		if (not path_name is None and os.path.isdir(path_name)):
			main_path = os.path.join(Constants.TEMP_DIRECTORY, Constants.COUNT_DNA_TEMP_DIRECTORY)
			## to prevent errors
			if path_name != main_path and path_name != (main_path + "/"):
				cmd = "rm -rf {}*".format(path_name)
				os.system(cmd)

	def move_file(self, sz_file_from, sz_file_to):
		if os.path.exists(sz_file_from):
			self.make_path(os.path.dirname(sz_file_to))
			cmd = "mv " + sz_file_from + " " + sz_file_to
			exist_status = os.system(cmd)
			if (exist_status != 0):
				self.logger_production.error('Fail to run: ' + cmd)
				self.logger_debug.error('Fail to run: ' + cmd)
				raise Exception("Fail to make a move a file") 
	
			### set attributes to file 664
			os.chmod(sz_file_to, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
			
	def link_file(self, sz_file_from, sz_file_to, test_destination_file = True):
		if os.path.exists(sz_file_from) and (not os.path.exists(sz_file_to) or not test_destination_file):
			self.make_path(os.path.dirname(sz_file_to))
			cmd = "ln -f -s " + sz_file_from + " " + sz_file_to
			exist_status = os.system(cmd)
			if (exist_status != 0):
				self.logger_production.error('Fail to run: ' + cmd)
				self.logger_debug.error('Fail to run: ' + cmd)
				raise Exception("Fail to link a file") 
			
	def copy_file(self, sz_file_from, sz_file_to):
		if os.path.exists(sz_file_from):
			self.make_path(os.path.dirname(sz_file_to))
			cmd = "cp " + sz_file_from + " " + sz_file_to
			exist_status = os.system(cmd)
			if (exist_status != 0):
				self.logger_production.error('Fail to run: ' + cmd)
				self.logger_debug.error('Fail to run: ' + cmd)
				raise Exception("Fail to make a copy a file") 
			
			### set attributes to file 664
			if os.path.isfile(sz_file_to):
				os.chmod(sz_file_to, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
			
	def make_path(self, path_name):
		if (not os.path.isdir(path_name) and not os.path.isfile(path_name)):
			cmd = "mkdir -p " + path_name
			exist_status = os.system(cmd)
			if (exist_status != 0):
				self.logger_production.error('Fail to run: ' + cmd)
				self.logger_debug.error('Fail to run: ' + cmd)
				raise Exception("Fail to make a path") 
		
	def is_integer(self, n_value):
		try:
			int(n_value)
			return True
		except ValueError: 
			return False


	def is_float(self, n_value):
		try:
			float(n_value)
			return True
		except ValueError: 
			return False

	def is_gzip(self, file_name):
		"""
		test if the file name ends in gzip
		""" 
		return file_name.endswith(".gz")
	
	def is_fastq_gz(self, file_name):
		"""
		test if the file name ends in gzip
		:out (True/False, type of file) OR raise Exception
		""" 
		if (not self.is_gzip(file_name)): raise Exception("File need to have suffix '.fastq.gz'/'.fq.gz'")
		
		try:
			sz_type = self.get_type_file(file_name)
			if (sz_type == Constants.FORMAT_FASTQ_illumina): return (True, Constants.FORMAT_FASTQ_illumina)
			if (sz_type == Constants.FORMAT_FASTQ_ont): return (True, Constants.FORMAT_FASTQ_ont)
		except OSError as e:
			self.logger_production.error("Fail to test '" + file_name + "' fastq.gz file: " + e.args[0])
			self.logger_debug.error("Fail to test '" + file_name + "' fastq.gz file: " + e.args[0])
		except Exception as e:
			self.logger_production.error("Fail to test '" + file_name + "' fastq.gz file: " + e.args[0])
			self.logger_debug.error("Fail to test '" + file_name + "' fastq.gz file: " + e.args[0])
			raise e
		raise Exception("File is not in fastq.gz format.")

	def get_type_file(self, file_name):
		"""
		return 'fasta' or 'fastq' 
		raise exception if can't detected
		"""
		
		vect_length = []
		with (gzip.open(file_name, mode='rt') if self.is_gzip(file_name) \
			else open(file_name, mode='r')) as handle_read:
			### read 100 lines
			try:
				count = 0
				for record in SeqIO.parse(handle_read, "fastq"):
					vect_length.append(len(str(record.seq)))
					count += 1
					if (count > 100): break
					#print("mean(vect_length): {} ".format(mean(vect_length)))
			except:
				pass
			
		### if read something in last SeqIO.parse
		if (len(vect_length) > 1):
			if (mean(sorted(vect_length, reverse=True)[:5]) <= Constants.MAX_LENGHT_ILLUMINA_FASQC_SEQ): return Constants.FORMAT_FASTQ_illumina
			if (mean(vect_length) > Constants.MIN_LENGHT_MINION_FASQC_SEQ): return Constants.FORMAT_FASTQ_ont
			raise Exception("Can not detect file format. Ensure Illumina fastq file.")
		
		## test fasta format
		with (gzip.open(file_name, mode='rt') if self.is_gzip(file_name) \
			else open(file_name, mode='r')) as handle_read:
		
			try:
				for record in SeqIO.parse(handle_read, Constants.FORMAT_FASTA):
					handle_read.close() 
					return Constants.FORMAT_FASTA
			except:
				pass
		
		raise Exception("File is not in fastq.gz format.")
	

	def is_fasta(self, sz_file_name):
		"""
		Test Fata file
		"""
		if (not os.path.exists(sz_file_name)): raise IOError(_("Error: File doens't exist: "  + sz_file_name))
		b_pass = False
		with open(sz_file_name) as handle:
			for line in handle:
				sz_temp = line.strip()
				if (len(sz_temp) == 0): continue
				if (sz_temp[0] == ">"): 
					b_pass = True
					break
				else: 
					raise IOError(_("Error: the file is not in FASTA format."))
		if (not b_pass): raise IOError(_("Error: file is not in FASTA format."))

		record_dict = SeqIO.index(sz_file_name, "fasta")
		if (len(record_dict) > 0): return len(record_dict)
		raise IOError(_("Error: file is not in FASTA format."))
	
	def test_sequences_same_length(self, sz_file_name):
		"""
		Test Fasta file if the sequences has the same length
		because of samtools faidx 
		"""
		if (not os.path.exists(sz_file_name)): raise IOError(_("Error: File doens't exist: "  + sz_file_name))
		handle = open(sz_file_name)
		vect_pos = []
		for line in handle:
			sz_temp = line.strip()
			if (len(sz_temp) == 0): continue
			if (sz_temp[0] == ">"):
				if (len(vect_pos) > 2):
					for i in range(1, len(vect_pos) - 1):
						if (vect_pos[0] != vect_pos[i]): return False
				vect_pos = [] 
			else:
				vect_pos.append(len(sz_temp))

		handle.close()
		if (len(vect_pos) > 2):
			for i in range(1, len(vect_pos) - 1):
				if (vect_pos[0] != vect_pos[i]): return False
		return True		


	def get_number_seqs_names_bigger_than(self, sz_file_name, size_limit_length, prefix_sum = 0):
		"""
		Test Fasta file
		"""
		if (not os.path.exists(sz_file_name)): raise IOError(_("Error: File doens't exist: "  + sz_file_name))
		record_dict = SeqIO.index(sz_file_name, "fasta")
		n_count = 0
		for key in record_dict:
			if ((len(key) + prefix_sum) > size_limit_length): n_count += 1
		return n_count
	
	def has_degenerated_bases(self, sz_file_name):
		"""
		Test Fasta file for degenerated bases
		"""
		if (not os.path.exists(sz_file_name)): raise IOError(_("Error: File doens't exist: "  + sz_file_name))
		handle = open(sz_file_name)
		sequence_name = ""
		going_to_test_sequence = False
		for line in handle:
			sz_temp = line.strip()
			if (len(sz_temp) == 0): continue
			if (sz_temp[0] == ">"):
				if (going_to_test_sequence):
					handle.close()
					raise Exception(_("Error: file is not in FASTA format."))
				sequence_name = sz_temp[1:].split(' ')[0]
				going_to_test_sequence = True 
			else:
				going_to_test_sequence = False
				if (len(sequence_name) == 0):
					handle.close()
					raise Exception(_("Error: file is not in FASTA format."))
				if (len(sz_temp) != len(re.findall("A|C|G|T", sz_temp.upper()))):
					handle.close()
					raise Exception(_("Error: sequence '{}' must have only 'A', 'C', 'T' and 'G' bases.".format(sequence_name)))
		handle.close()
		return False

	def get_max_length_fasta(self, sz_file_name):
		"""
		get max length fasta
		"""
		n_max = 0
		record_dict = SeqIO.index(sz_file_name, "fasta")
		if (len(record_dict) > 0): return len(record_dict)
		for seq in record_dict:
			if (len(record_dict[seq].seq) > n_max): n_max = len(record_dict[seq].seq)
		return n_max
	
	def get_total_length_fasta(self, sz_file_name):
		"""
		get max length fasta
		"""
		n_total = 0
		if os.path.exists(sz_file_name):
			record_dict = SeqIO.index(sz_file_name, "fasta")
			for seq in record_dict:
				n_total += len(record_dict[seq].seq)
		return n_total

							
	def is_genbank(self, sz_file_name):
		"""
		Test GenBank file
		"""
		with open(sz_file_name) as handle:
			b_pass = False
			for line in handle:
				sz_temp = line.strip()
				if (len(sz_temp) == 0): continue
				if (sz_temp.startswith("LOCUS", 0)): 
					b_pass = True
					break
				else: 
					raise IOError(_("Error: the file is not in GenBank format."))
			if (not b_pass): raise IOError(_("Error: the file is not in GenBank format."))

		n_number_locus = 0
		with open(sz_file_name) as handle:
			for record in SeqIO.parse(handle, "genbank"):
				n_number_locus += 1
			if (n_number_locus > 0): 
				return n_number_locus
		raise IOError(_("Error: the file is not in GenBank format."))


	def get_elements_and_genes(self, genbank_name):
		"""
		return a dictonary with elements and vect genes
		vect_genes = [[pos_start, pos_end, name, strand 1|-1], [...], ...]
		return: dt_data{ element_name : vect_genes, element_name_2 : vect_genes_2, ....} 
		"""
		geneticElement = GeneticElement()
		handle = open(genbank_name)
		for record in SeqIO.parse(handle, "genbank"):
			length = 0
			gene_add = False
			for features in record.features:
				if (features.type == 'source'):
					length = abs(features.location.end - features.location.start)
				elif (features.type == 'CDS'):
					for key_name in ['CDS', 'gene', 'locus_tag']:
						if (key_name in features.qualifiers):
							### has join or others... otherwise can be empty
							vect_feature_location = []
							if (isinstance(features.location, CompoundLocation)):
								for feature_location in features.location.parts:	## has more than one part
									vect_feature_location.append(FeatureLocationSimple(feature_location.start,
											feature_location.end, feature_location.strand))
								
							geneticElement.add_gene(record.name, length, Gene(features.qualifiers[key_name][0],
								int(features.location.start), int(features.location.end),
								features.location.strand, vect_feature_location))
							gene_add = True
							break
						
			### element without gene
			if (not gene_add and length > 0):
				geneticElement.add_gene(record.name, length, None)
				
		handle.close()
		return geneticElement
	
	def get_elements_from_db(self, reference, user):
		"""
		return vector with name of elements sorted
		"""
		manageDatabase = ManageDatabase()
		meta_key = MetaKeyAndValue.META_KEY_Elements_Reference
		meta_reference = manageDatabase.get_reference_metakey_last(reference, meta_key, MetaKeyAndValue.META_VALUE_Success)
		if (meta_reference == None):
			utils = Utils()
			geneticElement = utils.get_elements_and_genes(reference.get_reference_gbk(TypePath.MEDIA_ROOT))
			if (geneticElement == None): return None
			manageDatabase.set_reference_metakey(reference, user, meta_key,\
				MetaKeyAndValue.META_VALUE_Success, ','.join(geneticElement.get_sorted_elements()))
			return geneticElement.get_sorted_elements()
		else:
			return meta_reference.description.split(',')
		return None
	
	def get_elements_with_CDS_from_db(self, reference, user):
		"""
		return vector with name of elements sorted
		"""
		vect_elements = []
		for element in self.get_elements_from_db(reference, user):
			genetic_element = self.get_elements_and_cds_from_db(reference, user)
			if (not genetic_element is None and genetic_element.has_genes(element)):
				vect_elements.append(element)
		return vect_elements

	@transaction.atomic
	def get_elements_and_cds_from_db(self, reference, user):
		"""
		return geneticElement
		"""
		manageDatabase = ManageDatabase()
		meta_key = MetaKeyAndValue.META_KEY_Elements_And_CDS_Reference
		meta_reference = manageDatabase.get_reference_metakey_last(reference, meta_key, MetaKeyAndValue.META_VALUE_Success)
		if (meta_reference is None or meta_reference.creation_date < datetime.now()):
			utils = Utils()
			geneticElement = utils.get_elements_and_genes(reference.get_reference_gbk(TypePath.MEDIA_ROOT))
			if (geneticElement == None): return None
			
			if (not meta_reference is None):
				decodeCoverage = DecodeObjects()
				geneticElement_old = decodeCoverage.decode_result(meta_reference.description)
				if (geneticElement_old != geneticElement):
					manageDatabase.set_reference_metakey(reference, user, meta_key,\
						MetaKeyAndValue.META_VALUE_Success, geneticElement.to_json())
			else:
				manageDatabase.set_reference_metakey(reference, user, meta_key,\
					MetaKeyAndValue.META_VALUE_Success, geneticElement.to_json())
			return geneticElement
		else:
			decodeCoverage = DecodeObjects()
			geneticElement = decodeCoverage.decode_result(meta_reference.description)
			return geneticElement
		return None
	
	def get_vect_cds_from_element_from_db(self, element_name, reference, user):
		"""
		return geneticElement
		"""
		manageDatabase = ManageDatabase()
		meta_key = MetaKeyAndValue.META_KEY_Elements_And_CDS_Reference
		meta_reference = manageDatabase.get_reference_metakey_last(reference, meta_key, MetaKeyAndValue.META_VALUE_Success)
		if (meta_reference == None):
			utils = Utils()
			geneticElement = utils.get_elements_and_genes(reference.get_reference_gbk(TypePath.MEDIA_ROOT))
			if (geneticElement == None): return None
			manageDatabase.set_reference_metakey(reference, user, meta_key,\
				MetaKeyAndValue.META_VALUE_Success, geneticElement.to_json())
			return geneticElement.get_vect_gene_names(element_name)
		else:
			decodeCoverage = DecodeObjects()
			geneticElement = decodeCoverage.decode_result(meta_reference.description)
			return geneticElement.get_vect_gene_names(element_name)
		return None

	def read_text_file(self, file_name):
		"""
		read text file and put the result in an vector
		"""
		if (not os.path.exists(file_name)):
			self.logger_production.error("Fail to read '" + file_name)
			self.logger_debug.error("Fail to test '" + file_name)
			raise IOError(_("Error: file '" + file_name + "' doens't exist."))
		
		vect_out = []
		with open(file_name) as handle: 
			for line in handle:
				sz_temp = line.strip()
				if (len(sz_temp) == 0): continue
				vect_out.append(sz_temp)
		return vect_out
	
	def compare_locus_fasta_gb(self, fasta_file, gb_file):
		"""
		Test if all fasta locus are in gb file
		"""
		locus_fasta = self.is_fasta(fasta_file)
		locus_gb = self.is_genbank(gb_file)
		if (locus_gb != locus_fasta): raise ValueError(_("Number of locus are different from fasta to genbank."))
		
		record_dict = SeqIO.index(fasta_file, "fasta")
		try:
			handle_gb = open(gb_file)
			for record in SeqIO.parse(handle_gb, "genbank"):
				b_found = False
				for seq in record_dict:
					if (seq == record.name):
						if (len(record_dict[seq].seq) != len(record.seq)):
							handle_gb.close()
							raise ValueError(_("Different length. Fasta seq: %s length: %d; Fasta seq: %s length: %d." 
										% (seq, len(record_dict[seq].seq), record.name, len(record.seq))) )
						b_found = True
						break
				if (not b_found): 
					handle_gb.close()
					raise ValueError(_("This locus '" + record.name + "' is not in fasta file."))
		except AttributeError as e:
			handle_gb.close()
			raise ValueError(_(e.args[0]))
	
		handle_gb.close()
		return locus_fasta


	def get_all_files(self, directory):
		"""
		return all files from a directory
		"""
		return [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]


	def compress_files(self, software, file_name):
		"""
		compress files
		"""
		### get extension of output
		# extension = FileExtensions.FILE_BGZ if software == SoftwareNames.SOFTWARE_BGZIP_name else FileExtensions.FILE_GZ
		extension = FileExtensions.FILE_GZ
		
		## test if the file exists
		if (os.path.exists(file_name + extension)): return
		
		cmd = "{} -c {} > {}{}".format(software, file_name, file_name, extension)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to compress file") 
		
	def uncompress_files(self, software, file_name_in, file_name_out):
		"""
		compress files
		"""
		## test if the file exists
		if (not os.path.exists(file_name_in)): return
		
		cmd = "{} -cd {} > {}".format(software, file_name_in, file_name_out)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to compress file") 

	def str2bool(self, v):
		"""
		str to bool
		"""
		return v.lower() in ("yes", "true", "t", "1", "y")

# 	def is_all_tasks_finished(self, vect_tasks_id):
# 		"""
# 		return true if all tasks finished
# 		"""
# 		for task_id in vect_tasks_id:
# 			task = fetch(task_id)
# 			if (task == None): return False
# 		return True
	
	def is_all_tasks_finished_by_result(self, vect_tasks_id):
		"""
		return true if all tasks finished
		"""
		for task_id in vect_tasks_id:
			task_id.result(wait=-1)
		return True

# 	def count_tasks_finished_and_not(self, vect_tasks_id):
# 		"""
# 		return (count_finished, count_not_finished)
# 		"""
# 		(count_finished, count_not_finished) = (0, 0)
# 		for task_id in vect_tasks_id:
# 			task = fetch(task_id)
# 			if (task == None): count_not_finished += 1
# 			else: count_finished += 1
# 		return (count_finished, count_not_finished)

# 	def is_all_tasks_finished_success(self, vect_tasks_id):
# 		"""
# 		return true if all tasks finished
# 		"""
# 		for task_id in vect_tasks_id:
# 			task = fetch(task_id)
# 			if (task == None): continue
# 			if (not task.success): return False
# 		return True


	def add_freq_to_vcf(self, vcf_file, vcf_file_out):
		"""
		Used in snippy and freebayes
		add FREQ to VCF, FREQ=AO/DP
		vcffile must be gzip and tbi included
		"""
		FREQ = 'FREQ'
		
		#read the input file
		vcf_hanlder = pysam.VariantFile(vcf_file, "r")
		if (FREQ in vcf_hanlder.header.info): 
			vcf_hanlder.close()
			return
		
		vcf_hanlder_write = pysam.VariantFile(vcf_file_out, "w")
		vcf_hanlder.header.info.add(FREQ, number='A', type='Float', description='Ratio of AO/DP')

		## write the header
		for variant_header_records in vcf_hanlder.header.records:
			vcf_hanlder_write.header.add_record(variant_header_records)
		
		for variant_sample in vcf_hanlder.header.samples:
			vcf_hanlder_write.header.add_sample(variant_sample)
			
		for variant in vcf_hanlder:
			if ("DP" in variant.info and "AO" in variant.info):
				vect_ao_out = []
				for value_ in variant.info['AO']:
					vect_ao_out.append(float("{:.1f}".format((value_/float(variant.info['DP']) * 100))))
				variant.info[FREQ] = tuple(vect_ao_out)
			vcf_hanlder_write.write(variant)
			
		vcf_hanlder_write.close()
		vcf_hanlder.close()
		return vcf_file_out

	def add_freq_ao_ad_and_type_to_vcf(self, vcf_file, file_coverage, vcf_file_out, vcf_file_out_removed_by_filter,
					coverage_limit, freq_vcf_limit):
		""" add FREQ, AO, AF and TYPE to VCF, FREQ=AO/DP
		This case is used in MEDAKA only
		:param vcf_file_out_removed_by_filter -> can be None, keep the variants that are filter by  freq_vcf_limit
		vcffile must be gzip and tbi included
		:param coverage_limit -> filter by this coverage (this is necessary because medaka doesn't have)
		:param cut off for VCF freq
		returns: vcf file with freq, AO and AF 
		"""
		FREQ = 'FREQ'
		AO = 'AO'
		RO = 'RO'
		AF = 'AF'
		TYPE = "TYPE"
		DP_COMPOSED = "DP_COMPOSED"		### this is used to get 
		
		#read the input file
		vcf_hanlder = pysam.VariantFile(vcf_file, "r")
		if (FREQ in vcf_hanlder.header.info and AO in vcf_hanlder.header.info and\
			AF in vcf_hanlder.header.info): 
			vcf_hanlder.close()
			return
		
		vcf_hanlder_write = pysam.VariantFile(vcf_file_out, "w")
		if not vcf_file_out_removed_by_filter is None:
			vcf_hanlder_write_removed_by_filter = pysam.VariantFile(vcf_file_out_removed_by_filter, "w")
		if (not FREQ in vcf_hanlder.header.info): vcf_hanlder.header.info.add(FREQ, number='A', type='Float', description='Ratio of AO/(DPSP-AR)')
		if (not AO in vcf_hanlder.header.info): vcf_hanlder.header.info.add(AO, number='A', type='Integer', description='Alternate allele observation count, SR (alt1 fwd + alt1 rev, etc.)')
		if (not RO in vcf_hanlder.header.info): vcf_hanlder.header.info.add(RO, number='1', type='Integer', description='Reference allele observation count, SR (ref fwd + ref rev)')
		if (not AF in vcf_hanlder.header.info): vcf_hanlder.header.info.add(AF, number='R', type='Integer', description='Number of observation for each allele, SR (ref fwd + ref rev, alt1 fwd + alt1 rev, etc.)')
		if (not TYPE in vcf_hanlder.header.info): vcf_hanlder.header.info.add(TYPE, number='A', type='String', description='The type of allele, either snp, mnp, ins, del, or complex')
		if (not DP_COMPOSED in vcf_hanlder.header.info): vcf_hanlder.header.info.add(DP_COMPOSED, number='1', type='String',\
						description='Coverage at position (DPSP-AR)/(samtools -aa). First is collected by Medaka, Second is collected by samtools.')

		## write the header
		for variant_header_records in vcf_hanlder.header.records:
			vcf_hanlder_write.header.add_record(variant_header_records)
			if not vcf_file_out_removed_by_filter is None:
				vcf_hanlder_write_removed_by_filter.header.add_record(variant_header_records)
				
		for variant_sample in vcf_hanlder.header.samples:
			vcf_hanlder_write.header.add_sample(variant_sample)
			if not vcf_file_out_removed_by_filter is None:
				vcf_hanlder_write_removed_by_filter.header.add_sample(variant_sample)

		for variant in vcf_hanlder:
			### DP must be replaced by DPSP. DPSP is the sum of all reads Span and Ambiguous
			if ("SR" in variant.info and "DPSP" in variant.info and "AR" in variant.info):	## SR=0,0,15,6
				### don't process this VCF because has a low coverage
				total_deep = int(variant.info['DPSP']) - sum([int(_) for _ in variant.info['AR']])
				total_deep_samtools = self.get_coverage_by_pos(file_coverage,
									variant.chrom, variant.pos, variant.pos)
				if (coverage_limit > 0 and total_deep_samtools >= 0 and \
					total_deep_samtools < coverage_limit): continue
				if ( ((len(variant.info['SR']) // 2) - 1) != len(variant.alts)):
					#vcf_hanlder_write.write(variant) 
					continue		### different numbers of Alleles and References

				#### extra info				
				vect_out_ao = []	### AO
				out_ro = -1			### RO
				vect_out_af = []	### AF
				vect_out_freq = []	### FREQ
				vect_out_freq_filtered = []	### FREQ
				vect_out_type = []	### TYPE
				
				for value_ in range(0, len(variant.info['SR']), 2):
					if (value_ > 0):
						allele_count = int(variant.info['SR'][value_]) + int(variant.info['SR'][value_ + 1])
						
						if (total_deep > 0):
							### incongruences in Medaka, 
							### these values are collected in different stages of the Medaka workflow, (email from support@nanoporetech.com at 23 Dec 2020)
							if (total_deep <= allele_count): vect_out_freq.append(100)
							else:
								freq_value = allele_count/float(total_deep)
								if (freq_value >= freq_vcf_limit): 
									vect_out_freq.append(float("{:.1f}".format(freq_value * 100)))
								elif (not vcf_file_out_removed_by_filter is None): 
									vect_out_freq_filtered.append(float("{:.1f}".format(freq_value * 100)))
							#print(variant.pos, variant.ref, str(variant.alts), variant.info['DP'], vect_out_ao[-1], vect_out_freq[-1])
						
						vect_out_ao.append(allele_count)
						vect_out_type.append(self._get_type_variation(variant.ref, variant.alts[(value_ - 2) >> 1]))
					vect_out_af.append(int(variant.info['SR'][value_]) + int(variant.info['SR'][value_ + 1]))
					if (out_ro == -1): out_ro = int(variant.info['SR'][value_]) + int(variant.info['SR'][value_ + 1])
				
				### has some variant to save
				if (len(vect_out_freq) > 0):
					if (out_ro > -1): variant.info[RO] = tuple([out_ro])
					variant.info[AO] = tuple(vect_out_ao)
					variant.info[AF] = tuple(vect_out_af)
					variant.info[TYPE] = tuple(vect_out_type)	
					variant.info[DP_COMPOSED] = tuple(["{}/{}".format(total_deep, total_deep_samtools)])
					variant.info[FREQ] = tuple(vect_out_freq)
					
					### Only save the ones with FREQ
					vcf_hanlder_write.write(variant)
				
				### save the filtered
				if (len(vect_out_freq_filtered) > 0):
					if (out_ro > -1): variant.info[RO] = tuple([out_ro])
					variant.info[AO] = tuple(vect_out_ao)
					variant.info[AF] = tuple(vect_out_af)
					variant.info[TYPE] = tuple(vect_out_type)	
					variant.info[DP_COMPOSED] = tuple(["{}/{}".format(total_deep, total_deep_samtools)])
					variant.info[FREQ] = tuple(vect_out_freq_filtered)
					
					### Only save the ones with FREQ
					vcf_hanlder_write_removed_by_filter.write(variant)
			
		vcf_hanlder_write.close()
		vcf_hanlder.close()
		if not vcf_file_out_removed_by_filter is None: vcf_hanlder_write_removed_by_filter.close()
		return vcf_file_out

	def _get_type_variation(self, ref, alt):
		""" return type of variation based on change
		possible in variation type
				snp	Single Nucleotide Polymorphism	A => T
				mnp	Multiple Nuclotide Polymorphism	GC => AT
				ins	Insertion	ATT => AGTT
				del	Deletion	ACGG => ACG
				complex	Combination of snp/mnp 
		"""
		if (len(ref) == len(alt) and len(alt) == 1): return "snp"
		if (len(ref) > len(alt)): return "del"
		if (len(ref) < len(alt)): return "ins"
		if (len(ref) == len(alt) and len(alt) > 1): return "mnp"
		return "complex"

	def count_hits_from_tab(self, tab_file, vect_count_type):
		"""
		NP	864	snp	G	T	99.6855	CDS	+	864/1497	288/498	synonymous_variant c.864G>T p.Gly288Gly	locus_00005		Nucleoprotein
		count the hits by <50; 50-90
		possible in variation type
				snp	Single Nucleotide Polymorphism	A => T
				mnp	Multiple Nuclotide Polymorphism	GC => AT
				ins	Insertion	ATT => AGTT
				del	Deletion	ACGG => ACG
				complex	Combination of snp/mnp
		vect_count_type = ['snp', 'ins']
		"""
		count_hits = CountHits()
		with open(tab_file) as handle:
			for line in handle:
				sz_temp = line.strip()
				if (len(sz_temp) == 0 or sz_temp[0] == '#'): continue
				lst_data = sz_temp.split('\t')
				if (len(lst_data) > 5):
					lst_freq_data = lst_data[5].split(',')
					lst_type_var = lst_data[2].split(',')
					for i in range(0, len(lst_type_var)):
						if (lst_type_var[i] in vect_count_type):
							value_ = lst_freq_data[i]
							if (self.is_float(value_)):
								if (float(value_) < 50): count_hits.add_one_hits_less_50()
								elif (float(value_) < 91): count_hits.add_one_hits_50_90()
								else: count_hits.add_one_hits_more_90()
		return count_hits

	def get_variations_by_freq_from_tab(self, tab_file, vect_count_type):
		"""
		NP	864	snp	G	T	99.6855	CDS	+	864/1497	288/498	synonymous_variant c.864G>T p.Gly288Gly	locus_00005		Nucleoprotein
		count the hits by <50; 50-90
		possible in variation type
				snp	Single Nucleotide Polymorphism	A => T
				mnp	Multiple Nuclotide Polymorphism	GC => AT
				ins	Insertion	ATT => AGTT
				del	Deletion	ACGG => ACG
				complex	Combination of snp/mnp
		vect_count_type = ['snp', 'ins']
		out (dict_less_50, dict_more_50, dict_more_90)
		out: dict_less_50{ 'NP': [pos1, pos2, pos3, ...], 'BP1': [pos1, pos2, pos3, ...] ...}
		out: dict_more_50{ 'NP': [pos1, pos2, pos3, ...], 'BP1': [pos1, pos2, pos3, ...] ...}
		"""
		if (not os.path.exists(tab_file)): return ({}, {}, {})
		
		dict_less_50 = {}
		dict_more_50 = {}
		dict_more_90 = {}
		with open(tab_file) as handle:
			for line in handle:
				sz_temp = line.strip()
				if (len(sz_temp) == 0 or sz_temp[0] == '#'): continue
				lst_data = sz_temp.split('\t')
				if (len(lst_data) > 5):
					lst_freq_data = lst_data[5].split(',')
					lst_type_var = lst_data[2].split(',')
					for i in range(0, len(lst_type_var)):
						if (lst_type_var[i] in vect_count_type):
							value_ = lst_freq_data[i]
							if (self.is_float(value_)):
								if (float(value_) < 50):
									if (lst_data[0] in dict_less_50 and self.is_integer(lst_data[1])): 
										dict_less_50[lst_data[0]].append(int(lst_data[1]))
									elif (self.is_integer(lst_data[1])):
										dict_less_50[lst_data[0]] = [int(lst_data[1])]
								elif (float(value_) < 90):
									if (lst_data[0] in dict_more_50 and self.is_integer(lst_data[1])): 
										dict_more_50[lst_data[0]].append(int(lst_data[1]))
									elif (self.is_integer(lst_data[1])):
										dict_more_50[lst_data[0]] = [int(lst_data[1])]
								else:
									if (lst_data[0] in dict_more_90 and self.is_integer(lst_data[1])): 
										dict_more_90[lst_data[0]].append(int(lst_data[1]))
									elif (self.is_integer(lst_data[1])):
										dict_more_90[lst_data[0]] = [int(lst_data[1])]
		return (dict_less_50, dict_more_50, dict_more_90)


	def get_sequence_from_genbank(self, sequence_name, gene, genbank_file):
		"""
		get seq instance from genbank file
		"""
		handle_gb = open(genbank_file)
		for record in SeqIO.parse(handle_gb, "genbank"):
			if (record.name != sequence_name): continue
			for features in record.features:
				if (features.type == 'CDS'):
					for key_name in Constants.VECT_GENBANK_TAG_NAME:
						if (key_name in features.qualifiers and features.qualifiers[key_name][0] == gene.name):
							handle_gb.close()
							return features.location.extract(record).seq
		handle_gb.close()
		return None


	def filter_fasta_all_sequences(self, consensus_fasta, sample_name, coverage, limit_to_mask_consensus, out_dir):
		"""
		filter fasta file
		file name out: None if not saved, else output file name
		return True if has sequences, False doesn't have sequences
		"""
		file_name = os.path.join(out_dir, sample_name +  FileExtensions.FILE_FASTA)
		return self.filter_fasta_all_sequences_file(consensus_fasta, coverage, file_name, limit_to_mask_consensus, True)
	
	def filter_fasta_all_sequences_file(self, consensus_fasta, coverage, file_out, limit_to_mask_consensus, test_all_locus_good_coverage):
		"""
		filter fasta file
		:param limit_to_mask_consensus can be -1 if not defined for this project_sample
		file name out: None if not saved, else output file name
		return True if has sequences, False doesn't have sequences
		"""
		if (not os.path.exists(consensus_fasta)): return None
		locus_fasta = self.is_fasta(consensus_fasta)
		### if it has more than coverage, some problem exist
		### the number can be smaller because of consensus filter 
		if (locus_fasta > len(coverage.get_dict_data())): return None
		
		### need to have all with 100 more 9
		if (test_all_locus_good_coverage):
			for key in coverage.get_dict_data():
				if ( (limit_to_mask_consensus == -1 and not coverage.is_100_more_9(key)) or
					(limit_to_mask_consensus > 0 and not 
					coverage.ratio_value_coverage_bigger_limit(key, limit_to_mask_consensus)) ):
					return None
		
		### test if the out directory exists
		self.make_path(os.path.dirname(file_out))
		
		with open(consensus_fasta) as handle_consensus:
			record_dict = SeqIO.to_dict(SeqIO.parse(handle_consensus, "fasta"))
			with open(file_out, 'w') as handle_write:
				records = []
				for key in sorted(coverage.get_dict_data()):
					if ( (limit_to_mask_consensus == -1 and coverage.is_100_more_9(key)) or
						(limit_to_mask_consensus > 0 and coverage.ratio_value_coverage_bigger_limit(key, limit_to_mask_consensus)) ):
						records.append(record_dict[key])
				if (len(records) > 0): SeqIO.write(records, handle_write, "fasta")
		if (len(records) == 0 and os.path.exists(file_out)): os.unlink(file_out)
		return file_out if len(records) > 0 else None


	def filter_fasta_by_sequence_names(self, consensus_fasta, sample_name, sequence_name, coverage, gene, limit_to_mask_consensus, out_dir):
		"""
		:param limit_to_mask_consensus can be -1 if not defined for this project_sample
		:param coverage Can be None
		:param gene Can be None
		:param limit_to_mask_consensus if Coverage is None this value is not necessary
		Test if necessary reverse complement
		filter fasta file
		write a file for each element 
		file name out: None if not saved, else output file name
		coverage can be None, print all
		return True if has sequences, False doesn't have sequences
		"""
		if (not os.path.exists(consensus_fasta)): return None
		locus_fasta = self.is_fasta(consensus_fasta)
		### if it has more than coverage, some problem exist
		### the number can be smaller because of consensus filter
		if (not coverage is None and locus_fasta > len(coverage.get_dict_data())): return None
		
		file_name = os.path.join(out_dir, sample_name + FileExtensions.FILE_FASTA)
#		file_name = os.path.join(out_dir, sample_name + "_" + sequence_name + FileExtensions.FILE_FASTA)
		b_saved = False
		with open(consensus_fasta) as handle_consensus:
			record_dict = SeqIO.to_dict(SeqIO.parse(handle_consensus, "fasta"))
			with open(file_name, 'w') as handle:
				if (coverage is None):	### its the reference, there's no coverage and limit_to_mask
					if sequence_name in record_dict:
						seq_ref = record_dict[sequence_name]
						if not gene is None and not gene.is_forward(): seq_ref = seq_ref.reverse_complement()
						handle.write(">{}\n{}\n".format(sequence_name.replace(' ', '_'), str(seq_ref.seq).upper()))
						b_saved = True
				else:
					if sequence_name in coverage.get_dict_data() and sequence_name in record_dict and\
						( (limit_to_mask_consensus == -1 and coverage.is_100_more_9(sequence_name)) or\
						(limit_to_mask_consensus > 0 and coverage.ratio_value_coverage_bigger_limit(sequence_name, limit_to_mask_consensus)) ):
						
						seq_ref = record_dict[sequence_name]
						if not gene is None and not gene.is_forward(): seq_ref = seq_ref.reverse_complement()
						handle.write(">{}\n{}\n".format(sequence_name.replace(' ', '_'), str(seq_ref.seq).upper()))
						b_saved = True
		if (not b_saved and os.path.exists(file_name)): os.unlink(file_name)
		return file_name if b_saved else None
	
	def clean_fasta_names(self, vect_names_to_clean, in_file, out_file):
		"""
		clean name in fasta files
		Ex: eva0033_se.consensus.fasta/1..23 -> eva0033_se 
		Ex: vect_names_to_clean = [ FileExtensions.FILE_CONSENSUS_FASTA, FileExtensions.FILE_FASTA, FileExtensions.FILE_FA]
		"""
		if (not os.path.exists(in_file)): return None
		
		dict_out_name = {} 
		with open(in_file) as handle:
			with open(out_file, 'w') as handle_write:
				records = []
				for record in SeqIO.parse(handle, Constants.FORMAT_FASTA):
					name_clean = ""
					for names_to_clean in vect_names_to_clean:
						if (record.name.rfind(names_to_clean) != -1):
							name_clean = record.name[:record.name.rfind(names_to_clean)]
					if (len(name_clean) == 0): name_clean = record.name
					if (name_clean in dict_out_name): dict_out_name[name_clean] += 1
					else: dict_out_name[name_clean] = 0
					
					name_to_write = name_clean if dict_out_name[name_clean] == 0 else '{}_{}'.format(name_clean, dict_out_name[name_clean])
					records.append(SeqRecord( Seq(str(record.seq)), id = name_to_write, description=""))
				SeqIO.write(records, handle_write, "fasta")
		return out_file

	def clean_extension(self, file_name):
		"""
		remove extension
		"""
		if (file_name.rfind('.') != -1): return file_name[:file_name.rfind('.')]
		return file_name 


	def read_file_to_string(self, file_name):
		"""
		read a file to a string content
		"""
		if (not os.path.exists(file_name)): return None
		
		sz_return = ""
		with open(file_name) as handle:
			for line in handle:
				if (len(sz_return) > 0): sz_return += "\n"
				sz_return += line
		return sz_return
	
	
	def md5sum(self, filename):
		"""
		read file and transform to md5sum
		"""
		f = open(filename, mode='r')
		d = hashlib.md5()
		for buf in f.read(128):
			d.update(buf.encode())
		return d.hexdigest()
	
	
	def validate_date(self, date_text):
		"""
		The international format yyyy-mm-dd or yyyymmdd
		validate date time
		"""
		try:
			date_ = datetime.strptime(date_text, '%d-%m-%Y')
		except ValueError:
			try:
				date_ = datetime.strptime(date_text, '%d/%m/%Y')
			except ValueError:
				raise ValueError("Incorrect data format, should be dd/mm/YYYY")
		return date_
	
	def validate_date_format(self, date_text, format_, text_format):
		"""
		The international format yyyy-mm-dd or yyyy/mm/dd
		validate date time
		"""
		try:
			date_ = datetime.strptime(date_text, format_)
		except ValueError:
			raise ValueError("Incorrect data format, should be " + text_format)
		return date_
	
	def clean_fasta_file(self, in_file, out_file):
		"""
		clean fasta file from '-'
		"""
		if (not os.path.exists(in_file)): return
		
		with open(out_file, "w+") as output_file_handle:
			vect_sequences = []
			with open(in_file) as file_handle:
				for seq_record in SeqIO.parse(file_handle, "fasta"):
					# Take the current sequence
					#vect_sequences.append(SeqRecord(Seq(str(seq_record.seq).upper().replace('-', ''), IUPAC.ambiguous_dna), id=seq_record.id, description="", name=""))
					vect_sequences.append(SeqRecord(Seq(str(seq_record.seq).upper().replace('-', '')), id=seq_record.id, description="", name=""))
				SeqIO.write(vect_sequences, output_file_handle, "fasta")

	def from_genbank_to_bed(self, file_in, file_out):
		"""
		from genbank to bed
		"""
		header = """track name=Genes description="{} genes" itemRgb=On\n""".format(self.clean_extension(os.path.basename(file_in)) )
		with open(file_out, 'w+') as outh: 
			outh.write(header)
			with open(file_in) as file_handle:
				for record in SeqIO.parse(file_handle, "genbank") :
					dt_out = {}
					for feature in record.features:
						if feature.type == 'gene' or feature.type == 'CDS':
							start = feature.location.start.position
							stop = feature.location.end.position
							for key_name in Constants.VECT_GENBANK_TAG_NAME:
								name = self.__get_feature_from_seq_record(feature, key_name)
								if (name != None): break
							
							if (name == None or name in dt_out): continue
							dt_out[name] = 1
							if feature.strand < 0:
								strand = "-"
							else:
								strand = "+"
							bed_line = "{4}\t{0}\t{1}\t{2}\t1000\t{3}\t{0}\t{1}\t{5}\n".format(start, stop, name,\
													strand, record.id, '65,105,225' if strand == '+' else '105,180,65')
							outh.write(bed_line)
						
			outh.close()

	def __get_feature_from_seq_record(self, feature, tag_name):
		"""
		get feature from seq_record
		return None if not exist
		"""
		try:
			return feature.qualifiers[tag_name][0]
		except:
			return None


	def send_email(self, address, header, message):
		"""
		send an email
		"""
		send_mail(header, message, settings.DEFAULT_USER_EMAIL, [address])


	def grouped(self, l, n):
		"""
		group instances
		"""
		for i in range(0, len(l), n):
			yield l[i:i+n]

	def clean_name(self, name_to_clean, dict_to_clean = { ' ' : '_', '(' : '', ')' : '', '$' : '', '#' : '', '&' : '', '/' : '', '\\' : '' }):
		"""
		clean a name based on dictionary, dict_to_clean = { ' ' : '_', '(' : '' , ')' : '' }
		
		"""
		for key in dict_to_clean:
			name_to_clean = name_to_clean.replace(key, dict_to_clean[key])
		return name_to_clean
	
	def clean_genbank_version_name(self, file_name_in, file_name_out):
		"""
		clean version name and set it equal to ACESSION
		Equal: VERSION to ACCESSION
			ACCESSION   KX162693
			VERSION     KX162693
		THIS is important in snpEFF
		param	in: file_name_in
				out: file_name_out
		"""
		ACCESSION = "ACCESSION"
		VERSION = "VERSION"
		sub_accession = ""
		with open(file_name_in) as handle_in, open(file_name_out, 'w') as handle_out:
			for line in handle_in:
				sz_temp = line.strip()
				if (line.startswith(ACCESSION)):
					lst_data = sz_temp.split()
					if (len(lst_data) > 1): sub_accession = lst_data[1]
					handle_out.write(line)
					continue
				### now test version
				elif (line.startswith(VERSION) and len(sub_accession) > 0):
					handle_out.write("{}     {}\n".format(VERSION, sub_accession))
					sub_accession = ""
					continue
				else: sub_accession = ""
				handle_out.write(line)
	
	def get_number_sequeces_in_gff_file(self, sz_file_name):
		""" count the number of lines with sequences """
		if not os.path.exists(sz_file_name): return 0
		count = 0
		with open(sz_file_name) as handle_in:
			for line in handle_in:
				sz_temp = line.strip()
				if (len(sz_temp) == 0 or sz_temp[0] == '#'): continue
				count += 1
		return count
	
	
	def merge_fasta_first_sequence(self, path_where_files_are, out_file):
		"""
		Merge all fasta files into one
		"""
		vect_out_fasta = []
		dt_out_name = {}
		count = 1
		for file in glob.glob(path_where_files_are + "/*.fasta"):
			with open(file, "rU") as handle_fasta:
				for record in SeqIO.parse(handle_fasta, "fasta"):
					seq_name = os.path.splitext(os.path.basename(file))[0]
					record.id = seq_name
					while True:
						if record.id in dt_out_name:
							record.id = "{}_{}".format(seq_name, count)
							count += 1
						else: 
							dt_out_name[record.id] = 1
							break
						
					vect_out_fasta.append(record)
					break
					
		### write the output
		with open(out_file, "w") as handle_fasta_out:
			if (len(vect_out_fasta) > 0):
				SeqIO.write(vect_out_fasta, handle_fasta_out, "fasta")
		return len(vect_out_fasta)
	
	def merge_fasta_and_join_sequences(self, path_where_files_are, vect_elements, out_file):
		"""
		Merge all fasta files into one
		"""
		vect_out_fasta = []
		for file in glob.glob(path_where_files_are + "/*.fasta"):
			dict_seq_record = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
			sequence = ""
			b_not_all = False
			for element in vect_elements:
				if element in dict_seq_record:
					sequence += str(dict_seq_record[element].seq)
				else: b_not_all = True

			### 
			if (not b_not_all):
				seq_name = os.path.splitext(os.path.basename(file))[0]
				vect_out_fasta.append(SeqRecord(Seq(sequence), id=seq_name, description=""))
			
		### write the output
		with open(out_file, "w") as handle_fasta_out:
			if (len(vect_out_fasta) > 0):
				SeqIO.write(vect_out_fasta, handle_fasta_out, "fasta")
		return len(vect_out_fasta)
	
	def merge_fasta_files(self, vect_sample_path_and_name, out_file):
		"""
		:param vect_sample_path_and_name = [[path_file, sample_name],
					[path_file, sample_name], ... ]
		:param outfile file 
		"""
		vect_out_fasta = []
		dt_out_name = {}
		
		for data_file in vect_sample_path_and_name:
			count = 1
			if (not os.path.exists(data_file[0])): continue
			with open(data_file[0], "rU") as handle_fasta:
				for record in SeqIO.parse(handle_fasta, "fasta"):
					sample_name = data_file[1].replace(" ", "_")
					possible_name = "{}{}{}".format(sample_name,
						Constants.SEPARATOR_sample_record_id, record.id)
					while True:
						if possible_name in dt_out_name:
							possible_name = "{}{}{}_{}".format(sample_name,
								Constants.SEPARATOR_sample_record_id, record.id, count)
							count += 1
						else: 
							dt_out_name[possible_name] = 1
							break
					record.id = possible_name
					record.description = ""
					vect_out_fasta.append(record)
					
					### set all_consensus name in table project_sample
					if data_file[2] != -1:
						try:
							project_sample = ProjectSample.objects.get(id=data_file[2])
							## if the sample starts like this
							project_sample.seq_name_all_consensus = "{}{}".format(sample_name,
								Constants.SEPARATOR_sample_record_id)
							project_sample.save()
						except ProjectSample.DoesNotExist:	## need to create with last version
							continue
					
		### write the output
		with open(out_file, "w") as handle_fasta_out:
			if (len(vect_out_fasta) > 0):
				SeqIO.write(vect_out_fasta, handle_fasta_out, "fasta")
		return len(vect_out_fasta)

	def merge_fasta_files_and_join_multifasta(self, vect_sample_path_and_name, out_file, segment=None):
		"""
		:param vect_sample_path_and_name = [[path_file, name, ID],
					[path_file, name, ID], ... ]
		:param outfile file 
		"""
		vect_out_fasta_total = []
		dt_out_name = {}
		
		for data_file in vect_sample_path_and_name:
			
			if (not os.path.exists(data_file[0])): continue
			fasta_out = ""
			with open(data_file[0], "rU") as handle_fasta:
				for record in SeqIO.parse(handle_fasta, "fasta"):
					# In case of multiple segments we may concatenate all or retrieve a specific segment
					if( (segment is None) or (segment == record.id) ):
						fasta_out += str(record.seq)
			
			## none fasta sequence
			if (len(fasta_out) == 0): continue
			
			###
			count = 1
			sample_name = data_file[1].replace(" ", "_")
			possible_name = sample_name
			while True:
				if possible_name in dt_out_name:
					possible_name = "{}_{}".format(sample_name, count)
					count += 1
				else: 
					dt_out_name[possible_name] = 1
					break
						
			vect_out_fasta_total.append(SeqRecord(Seq(fasta_out), id=possible_name, description=""))
			
			### set all_consensus name in table dataset_consensus
			if data_file[2] != -1:
				try:
					dataset_consensus = DatasetConsensus.objects.get(id=data_file[2])
					## if the sample starts like this
					dataset_consensus.seq_name_all_consensus = possible_name
					dataset_consensus.save()
				except DatasetConsensus.DoesNotExist:	## need to create with last version
					continue
						
		### write the output
		with open(out_file, "w") as handle_fasta_out:
			if (len(vect_out_fasta_total) > 0):
				SeqIO.write(vect_out_fasta_total, handle_fasta_out, "fasta")
		return len(vect_out_fasta_total)

	def parse_amino_HGVS_code(self, amino_value):
		""" p.Asn292Asn -> p.Asn292Asn"""
		match = re.search("p.(?P<first_amino>[A-Za-z*]+)(?P<position>[0-9]+)(?P<second_amino>[A-Za-z*]+)", amino_value)
		
		if (not match is None):
			return "p.{}{}{}".format(self._get_amino_single_letter(match.group('first_amino')),
					match.group('position'),
					self._get_amino_single_letter(match.group('second_amino')))
		return ""

	def _get_amino_single_letter(self, string_amino):
		""" can have cases like this: SerVal """
		if (len(string_amino) % 3 != 0): return string_amino
		sz_out = ""
		for i in range(0, len(string_amino), 3):
			sz_out += protein_letters_3to1.get(string_amino[i:i+3], string_amino[i:i+3])
		return sz_out

	### get models from medaka
	def get_all_medaka_models(self):
		"""
		--model MODEL         Model to use. {r103_min_high_g345, r103_min_high_g360,
                        r103_prom_high_g360, r103_prom_snp_g3210,
                        r103_prom_variant_g3210, r10_min_high_g303,
                        r10_min_high_g340, r941_min_fast_g303,
                        r941_min_high_g303, r941_min_high_g330,
                        r941_min_high_g340_rle, r941_min_high_g344,
                        r941_min_high_g351, r941_min_high_g360,
                        r941_prom_fast_g303, r941_prom_high_g303,
                        r941_prom_high_g330, r941_prom_high_g344,
                        r941_prom_high_g360, r941_prom_high_g4011,
                        r941_prom_snp_g303, r941_prom_snp_g322,
                        r941_prom_snp_g360, r941_prom_variant_g303,
                        r941_prom_variant_g322, r941_prom_variant_g360}
                        (default: r941_min_high_g360)
        --threads THREADS     Number of threads used by inference. (default: 1)
        
        Medaka models are named to indicate: 
        i) the pore type, 
        ii) the sequencing device (min -> MinION, prom -> PromethION), 
        iii) the basecaller variant (only high and variant available in INSAFlu),
        iv) the Guppy basecaller version.
        Complete format:
        {pore}_{device}_{caller variant}_{caller version}

		:out return all models available
		"""
		software_names = SoftwareNames()
		temp_file = self.get_temp_file("medaka_models", ".txt")
 		
		cmd = "{} {} tools list_models > {}".format(
				software_names.get_medaka_env(),
				software_names.get_medaka(),
				temp_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			self.remove_file(temp_file)
			raise Exception("Fail to run medaka_consensus")
		
		vect_data = self.read_text_file(temp_file)
		vect_models = []
		for line in vect_data:
			if line.find('Available:') == 0:
				lst_data = line.replace('Available:', '').split(',')
				if len(lst_data) > 0:
					for model in lst_data:
						model_name = model.strip()
						if (len(model_name) > 0 and not self._exist_tag_name(model_name,
								software_names.get_medaka_remove_tags_model())): vect_models.append(model_name)
			elif (len(vect_models) > 0): break
		return vect_models

	def _exist_tag_name(self, model_name, vect_names_to_exclude):
		""" test if the model_name has some of the tags in vect_names_to_exclude"""
		for tag_to_test in vect_names_to_exclude:
			if (model_name.find(tag_to_test) != -1): return True
		return False
    
	def is_differente_fasta_size(self, file_name, percentage_diff):
		"""
		:out True if difference between first two sequences is more than percentage_diff
		"""
		vect_length = []
        
		### read 200 lines
		with open(file_name) as handle_in:
			for record in SeqIO.parse(handle_in, "fasta"):
				vect_length.append(len(str(record.seq)))
				if (len(vect_length) == 2):
					if (vect_length[0] > vect_length[1]):
						if ((100 - ((vect_length[1] / vect_length[0] * 100))) > percentage_diff): return True
						return False
					else:
						if ((100 - ((vect_length[0] / vect_length[1] * 100))) > percentage_diff): return True
						return False
		return False  
    
	def get_last_name_from_fasta(self, file_name):
		""" return last name in the fasta file """
		last_name = ""
		if (os.path.exists(file_name)):
			with open(file_name) as handle_in:
				for line in handle_in:
					sz_temp = line.strip()
					if (len(sz_temp) > 0 and sz_temp[0] == '>'): last_name = sz_temp.split()[0].replace('>', '')
		return last_name
    
	def get_number_sequences_fastq(self, file_name):
		""" return average and number of sequences
		:output (number seqs, average, std)
		"""
		temp_file =  self.get_temp_file("lines_and_average_", ".txt")
		cmd = "gzip -cd " + file_name + " | awk '{ s++; if ((s % 4) == 0) { count ++; size += length($0); " + \
			" sumsq += (length($0))^2; }  } END " + \
			" { print \"sequences: \", count,  \"average: \", size/count,  \"std: \", sqrt(sumsq/count - (size/count)^2) }' > " + temp_file
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			self.remove_temp_file(temp_file)
			raise Exception("Fail to run get_number of sequences in fastq file")
        
        ###
		vect_out = self.read_text_file(temp_file)
		if (len(vect_out) == 0):
			self.logger_production.error('can not read any data: ' + temp_file)
			self.logger_debug.error('can not read any data: ' + temp_file)
			self.remove_temp_file(temp_file)
			raise Exception("Can't read any data")
		vect_data = vect_out[0].split()
		if (len(vect_data) != 6):
			self.logger_production.error('can not parse this data: ' + vect_out[0])
			self.logger_debug.error('can not parse this data: ' + vect_out[0])
			self.remove_temp_file(temp_file)
			raise Exception("Can't read any data")
		average_value = "%.1f" % (float(vect_data[3]))
		std = "%.1f" % (float(vect_data[5]))
		
		self.remove_temp_file(temp_file)
		return (int(vect_data[1]), float(average_value), float(std))

	def get_coverage_by_pos(self, file_coverage, chr_name, position_start, position_end):
		"""
		:out coverage at a specific position
		"""
		if (not os.path.exists(file_coverage)): return -1

		software_names = SoftwareNames()
		temp_file = self.get_temp_file("get_number_fastq", ".txt")
		cmd = "{} {} {}:{}-{} > {}".format(software_names.get_tabix(), file_coverage,
								chr_name, position_start, position_end, temp_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			self.remove_file(temp_file)
			raise Exception("Fail to run gzip")

		### get number
		vect_lines = self.read_text_file(temp_file)
		self.remove_file(temp_file)
		if len(vect_lines) == 1 and len(vect_lines[0].split()) == 3 and self.is_integer(vect_lines[0].split()[2]):
			return int(vect_lines[0].split()[2])
        
		return -1
    
	def mask_sequence(self, sequence, mask_sites, mask_from_beginning, mask_from_end, mask_range):
		"""Mask characters at the given sites in a single sequence record, modifying the
		record in place.
		Parameters
		----------
		sequence : Bio.SeqIO.SeqRecord
			A sequence to be masked
		mask_sites: list[int]
			A list of site indexes to exclude from the FASTA.
		mask_from_beginning: int
			Number of sites to mask from the beginning of each sequence (default 0)
		mask_from_end: int
			Number of sites to mask from the end of each sequence (default 0)
		mask_invalid: bool
			Mask invalid nucleotides (default False)
		Returns
		-------
		Bio.SeqIO.SeqRecord
			Masked sequence in its original record object
		"""
		# Convert to a mutable sequence to enable masking with Ns.
		sequence_length = len(sequence.seq)
		beginning = int(mask_from_beginning) if not mask_from_beginning is None and len(mask_from_beginning) > 0 else 0 
		end = int(mask_from_end) if not mask_from_end is None and len(mask_from_end) > 0 else 0
		
		if beginning + end > sequence_length:
			beginning, end = sequence_length, 0
		
		seq = str(sequence.seq)[beginning:-end or None]
		masked_sequence = MutableSeq("N" * beginning + seq + "N" * end)
        
		# Replace all excluded sites with Ns.
		if not mask_sites is None and  len(mask_sites.split(',')[0]) > 0:
			for site in [int(_) - 1 for _ in mask_sites.split(',')]:
				if site < sequence_length:
					masked_sequence[site] = "N"
		
		if not mask_range is None:
			for data_ in mask_range.split(','):
				if (len(data_) > 0):
					for site in range(int(data_.split('-')[0]) - 1, int(data_.split('-')[1])):
						if site < sequence_length: masked_sequence[site] = "N"
		    
		sequence.seq = masked_sequence
		return sequence
	
	def mask_sequence_by_sites(self, consensus_from, consensus_to, genetic_elemets):
		vect_record_out = []
		## always work with the backup    
		with open(consensus_from, "rU") as handle_fasta:
			for record in SeqIO.parse(handle_fasta, "fasta"):
				masking_consensus = genetic_elemets.dt_elements_mask.get(record.id, MaskingConsensus())
				if masking_consensus.has_data():
					vect_record_out.append(self.mask_sequence(record,
						masking_consensus.mask_sites, masking_consensus.mask_from_beginning, 
						masking_consensus.mask_from_ends, masking_consensus.mask_regions))
				else: vect_record_out.append(record)
		if (len(vect_record_out) > 0):
			temp_file = self.get_temp_file("masked_seq_", ".fasta")
			with open(temp_file, "w") as handle_fasta_out:
				SeqIO.write(vect_record_out, handle_fasta_out, "fasta")

			### move temp consensus to original position, if has info
			if os.stat(temp_file).st_size > 0:
				self.move_file(temp_file, consensus_to)
			else: os.unlink(temp_file)


class ShowInfoMainPage(object):
	"""
	only a help class
	"""    
	def __init__(self):
		self.show_images_main_page = settings.SHOW_IMAGES_MAIN_PAGE
		if (len(settings.INSTITUTION_NAME) > 0):
			self.instutution_name = settings.INSTITUTION_NAME
			self.instutution_web_site = settings.INSTITUTION_WEB_SITE
