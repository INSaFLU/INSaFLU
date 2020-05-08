'''
Created on Oct 31, 2017

@author: mmp
'''
from constants.constants import Constants, FileExtensions, TypePath, TypeFile
from constants.meta_key_and_values import MetaKeyAndValue
from managing_files.manage_database import ManageDatabase
from Bio.Alphabet import IUPAC
from utils.result import GeneticElement, Gene
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from django.utils.translation import ugettext_lazy as _
from django.core.mail import send_mail
from utils.result import CountHits, DecodeObjects
from datetime import datetime
import os, random, gzip, hashlib, logging, ntpath, stat, re
from pysam import pysam
from django.conf import settings
from statistics import mean

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
	
	def get_path_to_fastq_file(self, user_id, sample_id):
		"""
		get the path to sample
		"""
		return os.path.join(Constants.DIR_PROCESSED_FILES_FASTQ, "userId_{0}".format(user_id), "sampleId_{0}".format(sample_id))

	def get_path_to_projec_file(self, user_id, project_id):
		"""
		get the path to project
		"""
		return os.path.join(Constants.DIR_PROCESSED_FILES_PROJECT, "userId_{0}".format(user_id), "projectId_{0}".format(project_id))
	
	def get_path_upload_file(self, user_id, type_file):
		"""
		user_id ->
		type_file -> TypeFile.TYPE_FILE_sample_file, TypeFile.TYPE_FILE_fastq_gz
		"""
		return os.path.join(Constants.DIR_PROCESSED_FILES_MULTIPLE_SAMPLES, "userId_{0}".format(user_id),\
				"{}".format('fastq_files' if type_file == TypeFile.TYPE_FILE_fastq_gz else 'csv_sample_file'))
		
	def get_unique_file(self, file_name):
		"""
		get unique file name from a file_name
		return '<path file_name>/<random number>_<file_name>'
		"""
		temp_file_name = "{}_{}".format(random.randrange(10000000, 99999999, 10), ntpath.basename(file_name))
		main_path = os.path.dirname(file_name)
		if (not os.path.exists(main_path)): os.makedirs(main_path)
		while 1:
			if (not os.path.exists(os.path.join(main_path, temp_file_name))): break
			temp_file_name = "{}_{}".format(random.randrange(10000000, 99999999, 10), ntpath.basename(file_name))
		return os.path.join(main_path, temp_file_name.replace(" ", "_"))

	def get_temp_file(self, file_name, sz_type):
		"""
		return a temp file name
		"""
		main_path = os.path.join(Constants.TEMP_DIRECTORY, Constants.COUNT_DNA_TEMP_DIRECTORY)
		if (not os.path.exists(main_path)): os.makedirs(main_path)
		while 1:
			return_file = os.path.join(main_path, "insa_flu_" + file_name + "_" + str(random.randrange(10000000, 99999999, 10)) + "_file" + sz_type)
			if (os.path.exists(return_file)): continue
			try:
				os.close(os.open(return_file, os.O_CREAT | os.O_EXCL))
				return return_file
			except FileExistsError:
				pass
			
	def get_temp_file_from_dir(self, dir_out, file_name, sz_type):
		"""
		return a temp file name
		"""
		if (not os.path.exists(dir_out)): os.makedirs(dir_out)
		while 1:
			return_file = os.path.join(dir_out, "insa_flu_" + file_name + "_" + str(random.randrange(10000000, 99999999, 10)) + "_file" + sz_type)
			if (os.path.exists(return_file)): continue
			try:
				os.close(os.open(return_file, os.O_CREAT | os.O_EXCL))
				return return_file
			except FileExistsError:
				pass


	def get_temp_dir(self):
		"""
		return a temp directory
		"""
		main_path = os.path.join(Constants.TEMP_DIRECTORY, Constants.COUNT_DNA_TEMP_DIRECTORY)
		if (not os.path.exists(main_path)): os.makedirs(main_path)
		while 1:
			return_path = os.path.join(main_path, "insa_flu_" + str(random.randrange(10000000, 99999999, 10)))
			if (not os.path.exists(return_path)):
				os.makedirs(return_path)
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
			if path_name == main_path or path_name == (main_path + "/"): cmd = "rm -r {}/*".format(path_name)
			else: cmd = "rm -r {}*".format(path_name)
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
			
	def link_file(self, sz_file_from, sz_file_to):
		if os.path.exists(sz_file_from):
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
				raise Exception("Fail to make a move a file") 
			
			### set attributes to file 664
			os.chmod(sz_file_to, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
			
	def make_path(self, path_name):
		if (not os.path.isdir(path_name) and not os.path.isfile(path_name)):
			cmd = "mkdir -p " + path_name
			os.system(cmd)
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
		raise Exception
		""" 
		if (not self.is_gzip(file_name)): raise Exception("File need to have suffix '.fastq.gz'")
		
		try:
			sz_type = self.get_type_file(file_name)
			if (sz_type == Constants.FORMAT_FASTQ): return True
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
		if (self.is_gzip(file_name)): handle = gzip.open(file_name, mode='rt')	## need to be opened in text mode, default it's in binary mode
		else: handle = open(file_name)
		vect_length = []
		try:
			count = 0
			for record in SeqIO.parse(handle, Constants.FORMAT_FASTQ):
				vect_length.append(len(str(record.seq)))
				if (count > 100): break
			handle.close()
#			print("mean(vect_length): {} ".format(mean(vect_length)))
		except:
			handle.close()
		
		### if read something in last SeqIO.parse
		if (len(vect_length) > 1):
			if (mean(vect_length) <= Constants.MAX_LENGHT_ILLUMINA_FASQC_SEQ): return Constants.FORMAT_FASTQ
			raise Exception("Can not detect file format. Ensure Illumina fastq file.")
		
		if (self.is_gzip(file_name)): handle = gzip.open(file_name, mode='rt')
		else: handle = open(file_name)
		try:
			for record in SeqIO.parse(handle, Constants.FORMAT_FASTA):
				handle.close() 
				return Constants.FORMAT_FASTA
		except:
			handle.close()
		
		raise Exception("File is not in fastq.gz format.")
	

	def is_fasta(self, sz_file_name):
		"""
		Test Fata file
		"""
		if (not os.path.exists(sz_file_name)): raise IOError(_("Error: File doens't exist: "  + sz_file_name))
		handle = open(sz_file_name)
		b_pass = False
		for line in handle:
			sz_temp = line.strip()
			if (len(sz_temp) == 0): continue
			if (sz_temp[0] == ">"): 
				b_pass = True
				break
			else: 
				handle.close()
				raise IOError(_("Error: the file is not in FASTA format."))
		handle.close()
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


	def get_number_seqs_names_bigger_than(self, sz_file_name, size_limit_length):
		"""
		Test Fasta file
		"""
		if (not os.path.exists(sz_file_name)): raise IOError(_("Error: File doens't exist: "  + sz_file_name))
		record_dict = SeqIO.index(sz_file_name, "fasta")
		n_count = 0
		for key in record_dict:
			if (len(key) > size_limit_length): n_count += 1
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
		record_dict = SeqIO.index(sz_file_name, "fasta")
		if (len(record_dict) > 0): return len(record_dict)
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
					for key_name in ['CDS', 'gene']:
						if (key_name in features.qualifiers):
							geneticElement.add_gene(record.name, length, Gene(features.qualifiers[key_name][0],
								int(features.location.start), int(features.location.end), features.location.strand))
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
		meta_reference = manageDatabase.get_reference_metakey(reference, meta_key, MetaKeyAndValue.META_VALUE_Success)
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


	def get_elements_and_cds_from_db(self, reference, user):
		"""
		return geneticElement
		"""
		manageDatabase = ManageDatabase()
		meta_key = MetaKeyAndValue.META_KEY_Elements_And_CDS_Reference
		meta_reference = manageDatabase.get_reference_metakey(reference, meta_key, MetaKeyAndValue.META_VALUE_Success)
		if (meta_reference == None):
			utils = Utils()
			geneticElement = utils.get_elements_and_genes(reference.get_reference_gbk(TypePath.MEDIA_ROOT))
			if (geneticElement == None): return None
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
		meta_reference = manageDatabase.get_reference_metakey(reference, meta_key, MetaKeyAndValue.META_VALUE_Success)
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
					vect_ao_out.append((value_/float(variant.info['DP']) * 100))
				variant.info[FREQ] = tuple(vect_ao_out)
			vcf_hanlder_write.write(variant)
			
		vcf_hanlder_write.close()
		vcf_hanlder.close()
		return vcf_file_out


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


	def filter_fasta_all_sequences(self, consensus_fasta, sample_name, coverage, out_dir):
		"""
		filter fasta file
		file name out: None if not saved, else output file name
		return True if has sequences, False doesn't have sequences
		"""
		file_name = os.path.join(out_dir, sample_name +  FileExtensions.FILE_FASTA)
		return self.filter_fasta_all_sequences_file(consensus_fasta, coverage, file_name, True)
	
	def filter_fasta_all_sequences_file(self, consensus_fasta, coverage, file_out, test_all_locus_good_coverage):
		"""
		filter fasta file
		file name out: None if not saved, else output file name
		return True if has sequences, False doesn't have sequences
		"""
		if (not os.path.exists(consensus_fasta)): return None
		locus_fasta = self.is_fasta(consensus_fasta)
		### doesn't have the same size, sequences in consensus/coverage
		if (locus_fasta != len(coverage.get_dict_data())): return None
		
		### need to have all with 100 more 9
		if (test_all_locus_good_coverage):
			for key in coverage.get_dict_data():
				if (not coverage.is_100_more_9(key)): return None
		
		### test if the out directory exists
		self.make_path(os.path.dirname(file_out))
		
		with open(consensus_fasta) as handle_consensus:
			record_dict = SeqIO.to_dict(SeqIO.parse(handle_consensus, "fasta"))
			with open(file_out, 'w') as handle_write:
				records = []
				for key in sorted(coverage.get_dict_data()):
					if (coverage.is_100_more_9(key)):
						records.append(SeqRecord(Seq(str(record_dict[key].seq)), id = key, description=""))
				if (len(records) > 0): SeqIO.write(records, handle_write, "fasta")
		if (len(records) == 0 and os.path.exists(file_out)): os.unlink(file_out)
		return file_out if len(records) > 0 else None


	def filter_fasta_by_sequence_names(self, consensus_fasta, sample_name, sequence_name, coverage, gene, out_dir):
		"""
		filter fasta file
		write a file for each element 
		file name out: None if not saved, else output file name
		coverage can be None, print all
		return True if has sequences, False doesn't have sequences
		"""
		if (not os.path.exists(consensus_fasta)): return None
		locus_fasta = self.is_fasta(consensus_fasta)
		### doesn't have the same size, sequences in consensus/coverage
		if (coverage != None and locus_fasta != len(coverage.get_dict_data())): return None
		
		file_name = os.path.join(out_dir, sample_name + FileExtensions.FILE_FASTA)
#		file_name = os.path.join(out_dir, sample_name + "_" + sequence_name + FileExtensions.FILE_FASTA)
		b_saved = False
		with open(consensus_fasta) as handle_consensus:
			record_dict = SeqIO.to_dict(SeqIO.parse(handle_consensus, "fasta"))
			with open(file_name, 'w') as handle:
				if (coverage == None):
					if sequence_name in record_dict:
						seq_ref = record_dict[sequence_name]
						if not gene is None and not gene.is_forward(): seq_ref = seq_ref.reverse_complement()
						handle.write(">{}\n{}\n".format(sequence_name.replace(' ', '_'), str(seq_ref.seq).upper()))
						b_saved = True
				else:
					if sequence_name in coverage.get_dict_data() and sequence_name in record_dict and coverage.is_100_more_9(sequence_name):
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
					vect_sequences.append(SeqRecord(Seq(str(seq_record.seq).upper().replace('-', ''), IUPAC.ambiguous_dna), id=seq_record.id, description="", name=""))
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
		send_mail(header, message, 'insaflu@insa.min-saude.pt', [address])


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
		
		
class ShowInfoMainPage(object):
	"""
	only a help class
	"""	
	def __init__(self):
		self.show_images_main_page = settings.SHOW_IMAGES_MAIN_PAGE
		if (len(settings.INSTITUTION_NAME) > 0):
			self.instutution_name = settings.INSTITUTION_NAME
			self.instutution_web_site = settings.INSTITUTION_WEB_SITE

		
		