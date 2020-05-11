'''
Created on Oct 28, 2017

@author: mmp
'''
import os, gzip, logging, cmd, re, humanfriendly
from utils.coverage import DrawAllCoverage
from utils.utils import Utils
from utils.parse_out_files import ParseOutFiles
from constants.constants import Constants, TypePath, FileType, FileExtensions
from constants.meta_key_and_values import MetaKeyAndValue
from manage_virus.models import UploadFile
from managing_files.models import Sample, ProjectSample, MixedInfectionsTag, ProcessControler
from manage_virus.uploadFiles import UploadFiles
from managing_files.manage_database import ManageDatabase
from utils.result import Result, SoftwareDesc, ResultAverageAndNumberReads, CountHits
from utils.parse_coverage_file import GetCoverage
from utils.mixed_infections_management import MixedInfectionsManagement
from settings.default_software import DefaultSoftware
from constants.software_names import SoftwareNames
from Bio import SeqIO
from BCBio import GFF
from django.conf import settings
from utils.process_SGE import ProcessSGE
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class Software(object):
	'''
	classdocs
	'''
	utils = Utils()
	software_names = SoftwareNames()
	
	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	def test_bgzip_and_tbi_in_vcf(self, vcf_file):
		"""
		test if a a vcf file has a gzip file, if not create it
		"""
		### create the gzip file
		self.utils.compress_files(self.software_names.get_bgzip(), vcf_file)
		### create the tabix
		if (vcf_file.endswith('.gz')):
			self.create_index_files(vcf_file)
		else:
			self.create_index_files(vcf_file + ".gz")


	def get_vect_type_files_to_copy(self, software):
		"""
		get type of files to copy
		"""
		if (software == SoftwareNames.SOFTWARE_SNIPPY_name):
			return [FileType.FILE_BAM, FileType.FILE_BAM_BAI, FileType.FILE_CONSENSUS_FA, FileType.FILE_DEPTH_GZ, FileType.FILE_DEPTH_GZ_TBI,\
				FileType.FILE_TAB, FileType.FILE_VCF_GZ, FileType.FILE_VCF, FileType.FILE_VCF_GZ_TBI, FileType.FILE_CSV, FileType.FILE_REF_FASTA]
		elif (software == SoftwareNames.SOFTWARE_FREEBAYES_name):
			return [FileType.FILE_VCF, FileType.FILE_TAB]


	def copy_files_to_project(self, project_sample, software, path_from):
		"""
		copy files to the project
		software : SOFTWARE_SNIPPY_name, SOFTWARE_FREEBAYES_name
		"""
		for type_file in self.get_vect_type_files_to_copy(software):
			if (type_file == FileType.FILE_CONSENSUS_FA):	## if .fa file pass to .fasta
				self.utils.copy_file(os.path.join(path_from, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, type_file, software))),\
					project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, software))
			elif (type_file == FileType.FILE_VCF and software == SoftwareNames.SOFTWARE_FREEBAYES_name):	## vcf file
				### create the bgzip file
				self.utils.compress_files(self.software_names.get_bgzip(), os.path.join(path_from, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, type_file, software))))
				### create the tabix
				self.create_index_files(os.path.join(path_from, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, software))))
	
				### copy both
				self.utils.copy_file(os.path.join(path_from, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, software))),\
					project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, software))
				self.utils.copy_file(os.path.join(path_from, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ_TBI, software))),\
					project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ_TBI, software))
				
				### if snippy copy also the vcf
				self.utils.copy_file(os.path.join(path_from, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, type_file, software))),\
								project_sample.get_file_output(TypePath.MEDIA_ROOT, type_file, software))
			elif (type_file == FileType.FILE_REF_FASTA):
				self.utils.copy_file(os.path.join(path_from, 'reference', os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, type_file, software))),\
					project_sample.get_file_output(TypePath.MEDIA_ROOT, type_file, software))
				## create the FAI index
				self.create_fai_fasta(project_sample.get_file_output(TypePath.MEDIA_ROOT, type_file, software))
			else:
				self.utils.copy_file(os.path.join(path_from, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, type_file, software))),\
					project_sample.get_file_output(TypePath.MEDIA_ROOT, type_file, software))

	def create_fai_fasta(self, fileFastaName):
		"""
		Create fai for a fasta file
		"""
		cmd = "%s faidx %s" % (self.software_names.get_samtools(), fileFastaName)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run samtools")
		return cmd

	def creat_new_reference_to_snippy(self, project_sample):
		
		### get temp file
		temp_file = self.utils.get_temp_file("new_reference", FileExtensions.FILE_FASTA)
		cmd = "perl %s %s %s" % (self.software_names.get_create_new_reference_to_snippy(), \
				project_sample.project.reference.get_reference_gbk(TypePath.MEDIA_ROOT), temp_file)
#		print(cmd)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			if (os.path.exists(temp_file)): os.unlink(temp_file)
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run create_new_reference_to_snippy")
		path_dest = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_REF_FASTA, self.software_names.get_snippy_name())
#		print("Copy {} to {}".format(temp_file, path_dest))
		self.utils.copy_file(temp_file, path_dest)
		self.create_fai_fasta(path_dest)
		if (os.path.exists(temp_file)): os.unlink(temp_file)

	def run_spades(self, fastq_1, fastq_2, out_dir):
		"""
		Run spades
		IF you have problems running spades.py change the spades.py file from:
		#!/usr/bin/env python
		to
		#!/usr/bin/env python3
		"""
		if (not os.path.exists(fastq_1)):
			self.logger_production.error('Fastq 1 not found: ' + fastq_1)
			self.logger_debug.error('Fastq 1 not found: ' + fastq_1)
			raise Exception('Fastq 1 not found: ' + fastq_1)
		
		if (fastq_2 is None or len(fastq_2) == 0 or not os.path.exists(fastq_2)): 
			cmd = "%s -s %s %s -t %d -o %s" % (self.software_names.get_spades(), fastq_1,
					self.software_names.get_spades_parameters_single(), settings.THREADS_TO_RUN_FAST, out_dir)
		else: cmd = "%s --pe1-1 %s --pe1-2 %s %s -t %d -o %s" % (self.software_names.get_spades(), fastq_1, fastq_2,\
					self.software_names.get_spades_parameters(), settings.THREADS_TO_RUN_FAST, out_dir)
		
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run spades")
		return cmd
	
	def is_exist_database_abricate(self, database):
		"""
		DATABASE	SEQUENCES	DATE
		argannot	1749	2017-Oct-30
		card	2158	2017-Oct-30
		"""
		temp_file = self.utils.get_temp_file('list_abricate', FileExtensions.FILE_TXT)
		cmd = "{} --list > {}".format(self.software_names.get_abricate(), temp_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to collect list from abricate")
		vect_out = self.utils.read_text_file(temp_file)
		if (os.path.exists(temp_file)): os.unlink(temp_file)
		for line in vect_out:
			for field in line.split('\t'):
				if (field.lower() == database.lower()): return True
		return False


	def create_database_abricate(self, database, file_name):
		"""
		create a database on abricate
		"""
		if (not os.path.isfile(file_name)): raise IOError("File not found: " + file_name) 
		cmd = "mkdir -p %s/%s" % (self.software_names.SOFTWARE_ABRICATE_DB, database)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run make directory")
		
		## copy the file
		self.utils.copy_file(file_name, os.path.join(self.software_names.SOFTWARE_ABRICATE_DB, database, "sequences"))
		
		cmd = "%s --setupdb" % (self.software_names.get_abricate())
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run abricate --setupdb")
		

	def run_abricate(self, database, file_name, parameters, out_file):
		"""
		Run abricator
		"""
		cmd = "%s --db %s %s --quiet %s > %s" % (self.software_names.get_abricate(), database,\
				parameters, file_name, out_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run abricate")
		return cmd


	def get_lines_and_average_reads(self, file_name):
		"""
		return (lines, average)
		
		raise Exception
		"""
		temp_file =  self.utils.get_temp_file("lines_and_average_", ".txt")
		cmd = "gzip -cd " + file_name + " | wc -l > " + temp_file
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run get_lines_and_average_reads")
		
		###
		vect_out = self.utils.read_text_file(temp_file)
		if (len(vect_out) == 0): return (0, 0)
		if (int(vect_out[0]) == 0): return (0, 0)
		
		cmd = "gzip -cd " + file_name + " | awk '{ s++; if ((s % 4) == 0) { count ++; size += length($0); }  } END { print \"sequences: \", count,  \"average: \", size/count }' > " + temp_file
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run get_lines_and_average_reads")
		
		###
		vect_out = self.utils.read_text_file(temp_file)
		if (len(vect_out) == 0):
			self.logger_production.error('can not read any data: ' + temp_file)
			self.logger_debug.error('can not read any data: ' + temp_file)
			raise Exception("Can't read any data")
		vect_data = vect_out[0].split()
		if (len(vect_data) != 4):
			self.logger_production.error('can not parse this data: ' + vect_out[0])
			self.logger_debug.error('can not parse this data: ' + vect_out[0])
			raise Exception("Can't read any data")
		if (self.utils.is_float(vect_data[3])): average_value = "%.1f" % (float(vect_data[3]))
		else: average_value = vect_data[3]
		
		self.utils.remove_temp_file(temp_file)
		return (vect_data[1], average_value)
	
	"""
	Global processing
	"""
#	@transaction.atomic
	def identify_type_and_sub_type(self, sample, fastq1_1, fastq1_2, owner, b_run_tests = False):
		"""
		Identify type and sub_type
		Because of the tests need to pass the files also as parameters
		"""
		
		manageDatabase = ManageDatabase()
		### temp dir out spades		
		out_dir_spades = self.utils.get_temp_dir()
		result_all = Result()
		try:
			cmd = self.run_spades(fastq1_1, fastq1_2, out_dir_spades)
			parameters = self.software_names.get_spades_parameters()
			if (fastq1_2 == None or len(fastq1_2) == 0): self.software_names.get_spades_parameters_single()
			result_all.add_software(SoftwareDesc(self.software_names.get_spades_name(), self.software_names.get_spades_version(), parameters))
		except Exception:
			result = Result()
			result.set_error("Spades (%s) fail to run" % (self.software_names.get_spades_version()))
			result.add_software(SoftwareDesc(self.software_names.get_spades_name(), self.software_names.get_spades_version(), self.software_names.get_spades_parameters()))
			manageDatabase.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False
		
		file_out = os.path.join(out_dir_spades, "contigs.fasta")
		if (not os.path.exists(file_out) or os.path.getsize(file_out) < 50):
			## save error in MetaKeySample
			result = Result()
			result.set_error("Spades (%s) fail to run, empty contigs.fasta file." % (self.software_names.get_spades_version()))
			result.add_software(SoftwareDesc(self.software_names.get_spades_name(), self.software_names.get_spades_version(), self.software_names.get_spades_parameters()))
			manageDatabase.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False
		
		### test id abricate has the database
		try:
			uploadFile = UploadFile.objects.order_by('-version')[0]
		except UploadFile.DoesNotExist:
			## save error in MetaKeySample
			result = Result()
			result.set_error("Abricate (%s) fail to run" % (self.software_names.get_abricate_version()))
			result.add_software(SoftwareDesc(self.software_names.get_abricate_name(), self.software_names.get_abricate_version(), self.software_names.get_abricate_parameters()))
			manageDatabase.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False

		if (not self.is_exist_database_abricate(uploadFile.abricate_name)):
			try:
				self.create_database_abricate(uploadFile.abricate_name, uploadFile.path)
			except Exception:
				result = Result()
				result.set_error("Abricate (%s) fail to run --setupdb" % (self.software_names.get_abricate_version()))
				result.add_software(SoftwareDesc(self.software_names.get_abricate_name(), self.software_names.get_abricate_version(), self.software_names.get_abricate_parameters()))
				manageDatabase.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
				cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
				return False
		
		## run abricate
		out_file_abricate = self.utils.get_temp_file("temp_abricate", ".txt")
		try:
			cmd = self.run_abricate(uploadFile.abricate_name, file_out, SoftwareNames.SOFTWARE_ABRICATE_PARAMETERS, out_file_abricate)
			result_all.add_software(SoftwareDesc(self.software_names.get_abricate_name(), self.software_names.get_abricate_version(),\
						self.software_names.get_abricate_parameters() + " for type/subtype identification"))
		except Exception:
			result = Result()
			result.set_error("Abricate (%s) fail to run" % (self.software_names.get_abricate_version()))
			result.add_software(SoftwareDesc(self.software_names.get_abricate_name(), self.software_names.get_abricate_version(),\
							self.software_names.get_abricate_parameters() + " for type/subtype identification"))
			manageDatabase.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False
		
		if (not os.path.exists(out_file_abricate)):
			## save error in MetaKeySample
			result = Result()
			result.set_error("Abricate (%s) fail to run" % (self.software_names.get_abricate_version()))
			result.add_software(SoftwareDesc(self.software_names.get_abricate(), self.software_names.get_abricate_version(), self.software_names.get_abricate_parameters()))
			manageDatabase.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False

		parseOutFiles = ParseOutFiles()
		(vect_data, clean_abricate_file) = parseOutFiles.parse_abricate_file(out_file_abricate, os.path.basename(sample.get_abricate_output(TypePath.MEDIA_ROOT)),\
											SoftwareNames.SOFTWARE_SPAdes_CLEAN_HITS_BELLOW_VALUE)
		## copy the abricate output 
		self.utils.copy_file(clean_abricate_file, sample.get_abricate_output(TypePath.MEDIA_ROOT))

		try:
			contigs_2_sequences = Contigs2Sequences(b_run_tests)
			(out_file_clean, clean_abricate_file) = contigs_2_sequences.identify_contigs(file_out,\
					os.path.basename(sample.get_draft_contigs_abricate_output(TypePath.MEDIA_ROOT)))
			## copy the contigs from spades
			if (os.path.exists(out_file_clean)): self.utils.copy_file(out_file_clean, sample.get_draft_contigs_output(TypePath.MEDIA_ROOT))
			if (os.path.exists(clean_abricate_file)): self.utils.copy_file(clean_abricate_file, sample.get_draft_contigs_abricate_output(TypePath.MEDIA_ROOT))
			result_all.add_software(SoftwareDesc(self.software_names.get_abricate_name(), self.software_names.get_abricate_version(),\
						self.software_names.get_abricate_parameters_mincov_30() + " for segments/references assignment"))
		except Exception:
			result = Result()
			result.set_error("Abricate (%s) fail to run" % (self.software_names.get_abricate_version()))
			result.add_software(SoftwareDesc(self.software_names.get_abricate_name(), self.software_names.get_abricate_version(),\
							self.software_names.get_abricate_parameters_mincov_30() + " for segments/references assignment"))
			manageDatabase.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			return False
		
		uploadFiles = UploadFiles()
		vect_data = uploadFiles.uploadIdentifyVirus(vect_data, uploadFile.abricate_name)
		if (len(vect_data) == 0):
			## save error in MetaKeySample
			result = Result()
			result.set_error("Fail to identify type and sub type")
			result.add_software(SoftwareDesc(self.software_names.get_abricate(), self.software_names.get_abricate_version(), self.software_names.get_abricate_parameters()))
			manageDatabase.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm %s" % (out_file_abricate); os.system(cmd)
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False
		
		
		for identify_virus in vect_data:
			sample.identify_virus.add(identify_virus)
		sample.save()
		
		## save everything OK
		manageDatabase.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Success, "Success, Spades(%s), Abricate(%s)" % (self.software_names.get_spades_version(), self.software_names.get_abricate_version()))
		manageDatabase.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample_Software, MetaKeyAndValue.META_VALUE_Success, result_all.to_json())
		cmd = "rm %s" % (out_file_abricate); os.system(cmd)
		self.utils.remove_dir(out_dir_spades)
		self.utils.remove_file(out_file_clean)
		self.utils.remove_file(clean_abricate_file)
		return True


	def run_fastq(self, file_name_1, file_name_2):
		"""
		run fastQ, return output directory
		-o OUT_FOLDER --nogroup --format fastq --threads 10 --dir OUT_FOLDER FILE1 FILE2
		"""
		temp_dir = self.utils.get_temp_dir()
		if (not file_name_2 is None and len(file_name_2) > 0):
			cmd = "%s -o %s --nogroup --format fastq --threads %d --dir %s %s %s" % (self.software_names.get_fastq(), temp_dir, settings.THREADS_TO_RUN_FASTQC, 
										temp_dir, file_name_1, file_name_2)
		else: 
			cmd = "%s -o %s --nogroup --format fastq --threads %d --dir %s %s" % (self.software_names.get_fastq(), temp_dir, settings.THREADS_TO_RUN_FASTQC, temp_dir, file_name_1)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run run_fastq")
		return temp_dir


	def run_prokka(self, fasta_file_name, original_file_name):
		"""
		run prokka, in fasta file
		out: genbank
		{bpipe_prokka} FILE1 --kingdom Viruses --locustag locus --kingdom Viruses --locustag locus --genus Influenzavirus 
			--species Influenzavirus --strain ref_PREFIX_FILES_OUT --outdir OUT_FOLDER/PREFIX_FILES_OUT --prefix PREFIX_FILES_OUT
		"""
		if (not os.path.exists(fasta_file_name) or os.path.getsize(fasta_file_name) == 0): raise Exception("File doesn't exist")
		temp_dir = self.utils.get_temp_dir()
		name_strain = self.utils.clean_extension(os.path.basename(original_file_name))
		cmd = "{} {} {} --strain {} --force --outdir {} --prefix {}".format(\
					self.software_names.get_prokka(), fasta_file_name, self.software_names.get_prokka_parameters(), name_strain, temp_dir, name_strain)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run prokka")
		
		## clean /mol_type=
		gbk_file = os.path.join(temp_dir, self.utils.clean_extension(original_file_name) + FileExtensions.FILE_GBK)
		temp_file = self.utils.get_temp_file_from_dir(temp_dir, 'new_file', FileExtensions.FILE_GBK)
		cmd = "grep -E -v '/mol_type=' {} > {}".format(gbk_file, temp_file)
		os.system(cmd);
		self.utils.move_file(temp_file, gbk_file)
		return temp_dir
	
	def run_mauve(self, dir_to_process, out_file):
		"""
		run mauve
		out: out_file
		--output=alignment_all_96_samples_H3.xmfa *fasta
		"""
		dir_present = os.getcwd()
		os.chdir(dir_to_process)
		cmd = "{} --output={} *fasta".format(self.software_names.get_mauve(), out_file)
		exist_status = os.system(cmd)
		os.chdir(dir_present)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run progressive mauve")
		return out_file

	def run_convert_mauve(self, input_file, out_file):
		"""
		run convert mauve
		out: out_file
		"""
		temp_file = self.utils.get_temp_file('convert_mauve_software', FileExtensions.FILE_FASTA)
		cmd = "perl {} -c -i {} -o {} -f fasta".format(self.software_names.get_convert_mauve(), input_file, temp_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run convert mauve")
		
		#### clean names
		### important, must have this order
		vect_names_to_clean = [ FileExtensions.FILE_CONSENSUS_FASTA, FileExtensions.FILE_FASTA, FileExtensions.FILE_FA] 
		self.utils.clean_fasta_names(vect_names_to_clean, temp_file, out_file)
		os.unlink(temp_file)
		return out_file

	def run_mafft(self, input_file, out_file, parameters):
		"""
		run mafft
		out: out_file
		"""
		cmd = "{}; {} {} --thread {} {} > {}".format(self.software_names.get_mafft_set_env_variable(), self.software_names.get_mafft(),\
							parameters, settings.THREADS_TO_RUN_SLOW, input_file, out_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run mafft")
		return out_file
	
	def run_seqret_nex(self, input_file, out_file):
		"""
		run mafft
		out: out_file
		"""
		cmd = "{} -sequence {} {} -outseq {}".format(self.software_names.get_seqret(), input_file,\
							self.software_names.get_seqret_nex_parameters(), out_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run seqret")
		return out_file

	def run_fasttree(self, input_file, out_file, parameters):
		"""
		run fasttree
		out: out_file
		"""
		cmd = "{} {} {} > {}".format(self.software_names.get_fasttree(), parameters, input_file, out_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run fasttree")
		return out_file
	
	def run_trimmomatic(self, file_name_1, file_name_2, sample_name, user = None):
		"""
		run trimmomatic
		return output directory
		
		#${bpipe_trimmomatic} PE -threads 3 -basein FILE1 -baseout OUT_FOLDER/PREFIX_FILES_OUT.fastq.gz SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:55 TOPHRED33
		
		PE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
   			or: 
		SE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] <inputFile> <outputFile> <trimmer1>
		"""

		## get dynamic parameters
		if (user is None):
			parameters = self.software_names.get_trimmomatic_parameters()
		else:
			default_software = DefaultSoftware()
			parameters = default_software.get_parameters(self.software_names.get_trimmomatic_name(), user)

		### run software
		temp_dir = self.utils.get_temp_dir()
		if (file_name_2 is None or len(file_name_2) == 0):
			cmd = "java -jar %s SE -threads %d %s %s_1P.fastq.gz %s" % (self.software_names.get_trimmomatic(), settings.THREADS_TO_RUN_FAST, file_name_1, 
					os.path.join(temp_dir, sample_name), parameters)
		else:
			### need to make links the files to trimmomatic identify the _R1_ and _R2_ 
			new_file_name = os.path.join(temp_dir, 'name_R1_001.fastq.gz')
			cmd = "ln -s {} {}".format(file_name_1, new_file_name)
			os.system(cmd)
			cmd = "ln -s {} {}".format(file_name_2, os.path.join(temp_dir, 'name_R2_001.fastq.gz'))
			os.system(cmd)
			cmd = "java -jar %s PE -threads %d -basein %s -baseout %s.fastq.gz %s" % (self.software_names.get_trimmomatic(), settings.THREADS_TO_RUN_FAST, 
										new_file_name, os.path.join(temp_dir, sample_name), parameters)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run trimmomatic")
		return (temp_dir, parameters)

	def run_snippy(self, file_name_1, file_name_2, path_reference, sample_name):
		"""
		run snippy
		return output directory
		
		## ./snippy --cpus 3 --outdir /tmp/insafli/xpto --prefix xpto --ref ~/fluWeb/ref_H3.fasta --mapqual 20 --mincov 10 --minfrac 0.51 --R1 ~/fluWeb/EVA001_S66_L001_R1_001.fastq.gz --R2 ~/fluWeb/EVA001_S66_L001_R2_001.fastq.gz SNIPPY

		Out file		
		[06:29:16] * /tmp/insafli/xpto/xpto.aligned.fa
		[06:29:16] * /tmp/insafli/xpto/xpto.bam
		[06:29:16] * /tmp/insafli/xpto/xpto.bam.bai
		[06:29:16] * /tmp/insafli/xpto/xpto.bed
		[06:29:16] * /tmp/insafli/xpto/xpto.consensus.fa
		[06:29:16] * /tmp/insafli/xpto/xpto.consensus.subs.fa
		[06:29:16] * /tmp/insafli/xpto/xpto.csv
		[06:29:16] * /tmp/insafli/xpto/xpto.depth.gz
		[06:29:16] * /tmp/insafli/xpto/xpto.depth.gz.tbi
		[06:29:16] * /tmp/insafli/xpto/xpto.filt.subs.vcf
		[06:29:16] * /tmp/insafli/xpto/xpto.filt.subs.vcf.gz
		[06:29:16] * /tmp/insafli/xpto/xpto.filt.subs.vcf.gz.tbi
		[06:29:16] * /tmp/insafli/xpto/xpto.filt.vcf
		[06:29:16] * /tmp/insafli/xpto/xpto.gff
		[06:29:16] * /tmp/insafli/xpto/xpto.html
		[06:29:16] * /tmp/insafli/xpto/xpto.log
		[06:29:16] * /tmp/insafli/xpto/xpto.raw.vcf
		[06:29:16] * /tmp/insafli/xpto/xpto.tab
		[06:29:16] * /tmp/insafli/xpto/xpto.txt
		[06:29:16] * /tmp/insafli/xpto/xpto.vcf
		[06:29:16] * /tmp/insafli/xpto/xpto.vcf.gz
		[06:29:16] * /tmp/insafli/xpto/xpto.vcf.gz.tbi
		"""
		temp_dir = os.path.join(self.utils.get_temp_dir(), sample_name)
		if (file_name_2 is None or len(file_name_2) == 0):
			cmd = "%s --cpus %d --outdir %s --prefix %s --ref %s %s --se %s" %\
				(self.software_names.get_snippy(), settings.THREADS_TO_RUN_SLOW, temp_dir, sample_name,
				path_reference, self.software_names.get_snippy_parameters(), file_name_1)
		else:
			cmd = "%s --cpus %d --outdir %s --prefix %s --ref %s %s --R1 %s --R2 %s" %\
				(self.software_names.get_snippy(), settings.THREADS_TO_RUN_SLOW, temp_dir, sample_name,
				path_reference, self.software_names.get_snippy_parameters(), file_name_1, file_name_2)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			cmd = "rm -r %s*" % (temp_dir); os.system(cmd)
			raise Exception("Fail to run snippy")
		return temp_dir

	def run_genbank2gff3(self, genbank, out_file):
		"""
		
		"""
		temp_file = self.utils.get_temp_file("gbk_to_gff3", ".txt")
		
		## set VERSION ID equal to ACCESSION
		out_file_gb = self.utils.get_temp_file("file_name", ".gb")
		self.utils.clean_genbank_version_name(genbank, out_file_gb)
		
		with open(temp_file, "w") as out_handle:
			with open(out_file_gb) as in_handle:
				GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)

		### remove clean version gb file
		os.unlink(out_file_gb)
		
		### filter file
	#	vect_filter = ['remark', 'source', 'gene']
		vect_pass = ['CDS']
		with open(temp_file) as handle:
			with open(out_file, "w") as handle_write:
				for line in handle:
					sz_temp = line.strip()
					if (len(sz_temp) == 0 or sz_temp.find('# Input') == 0 or sz_temp.find('# GFF3 saved') == 0): continue
					if (sz_temp.find('##FASTA') == 0): break
					if (sz_temp[0] == '#'): handle_write.write(sz_temp + "\n")
					elif (len(sz_temp.split('\t')) > 3 and sz_temp.split('\t')[2] in vect_pass): handle_write.write(sz_temp + "\n")
		os.unlink(temp_file)
		return out_file

	def get_snpeff_config(self, fasta_file):
		"""
		return a file for snpEff config for the fasta
		genome_name: this name is going to set in properties in the file
		"""
		temp_file = self.utils.get_temp_file("snpEff_config", ".config")
		if (not os.path.exists(self.software_names.get_snp_eff_config())):
			raise IOError("Error: file not found {}".format(self.software_names.get_snp_eff_config()))
		
		self.utils.copy_file(self.software_names.get_snp_eff_config(), temp_file)
		handle = open(temp_file, 'a')
		base_file_name = os.path.basename(fasta_file)
		base_file_name = base_file_name[0: base_file_name.rfind('.')]
		handle.write("{}.genome : {} reference".format(base_file_name, Constants.INSAFLU_NAME))
		
		### open fasta
		if (self.utils.is_gzip(fasta_file)): handle_fasta = gzip.open(fasta_file, mode='rt')
		else: handle_fasta = open(fasta_file)
		vect_names = []
		for record in SeqIO.parse(handle_fasta, Constants.FORMAT_FASTA):
			vect_names.append(record.name)
		handle_fasta.close()
		
		handle.write("\n{}.chromosome : {}".format(base_file_name, ", ".join(vect_names)))
		for chromosome in vect_names:
			handle.write("\n{}.{}.codonTable : Bacterial_and_Plant_Plastid".format(base_file_name, chromosome))
		handle.close()
		return (base_file_name, temp_file)

	def run_snpEff(self, fasta_file, genbank, vcf_file, out_file):
		"""
		./snpEff ann -no-downstream -no-upstream -no-intergenic -no-utr -c ../path_to/reference/snpeff.config -dataDir . -noStats ref sample.vcf > sample_annot.vcf
		"""
		temp_dir = self.utils.get_temp_dir()
		(fasta_file_name, snpeff_config) = self.get_snpeff_config(fasta_file)
		
		temp_vcf_file = os.path.join(temp_dir, os.path.basename(vcf_file))
		self.utils.copy_file(vcf_file, temp_vcf_file)
		
		## create the database
		out_gff_file = self.utils.get_temp_file('temp_gff', '.gff')
		self.run_genbank2gff3(genbank, out_gff_file)
		
		### count sequences, if none return None
		if (self.utils.get_number_sequeces_in_gff_file(out_gff_file) == 0):
			os.unlink(out_gff_file)
			os.unlink(snpeff_config)
			self.utils.remove_dir(temp_dir)
			return None

		datase_dir = "{}".format(fasta_file_name)
		cmd = "mkdir -p {}".format(os.path.join(temp_dir, datase_dir)); os.system(cmd)
		cmd = "mkdir -p {}".format(os.path.join(temp_dir, 'genomes')); os.system(cmd)
		self.utils.copy_file(out_gff_file, os.path.join(temp_dir, datase_dir, 'genes.gff'))
		temp_file = os.path.join(temp_dir, 'genomes', fasta_file_name + '.fa')
		self.utils.copy_file(fasta_file, temp_file)
		os.unlink(out_gff_file)
		
		## indexing database
		## snpEff build -c reference/snpeff.config -dataDir . -gff3 ref 2>> run_snippy2_1.log
		cmd = "%s build -c %s -dataDir %s -gff3 %s" % (self.software_names.get_snp_eff(),\
						snpeff_config, temp_dir, fasta_file_name)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			os.unlink(snpeff_config)
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to create snpEff database")
		
		### create the annotation
		cmd = "%s ann %s -c %s -dataDir %s %s %s > %s" % (self.software_names.get_snp_eff(), self.software_names.get_snp_eff_parameters(),\
						snpeff_config, temp_dir, fasta_file_name, temp_vcf_file, out_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			os.unlink(snpeff_config)
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run snpEff")
		os.unlink(snpeff_config)
		self.utils.remove_dir(temp_dir)
		return out_file
	
	def run_snippy_vcf_to_tab(self, fasta, genbank, vcf_file, out_file):
		"""
		./snippy-vcf_to_tab [options] --ref ref.fa [--gff ref.gff] --vcf snps.vcf > snp.tab
		"""
		
		temp_file = self.utils.get_temp_file("snippy_vcf_to_tab", ".tab") 
		self.run_genbank2gff3(genbank, temp_file)
		
		cmd = "%s --ref %s --gff %s --vcf %s > %s" % (self.software_names.get_snippy_vcf_to_tab(), fasta, temp_file, vcf_file, out_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			os.unlink(temp_file)
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run snippy-vcf-to-tab")
		os.unlink(temp_file)
		return out_file

	def run_freebayes(self, bam_file, reference_fasta, genbank_file, sample_name):
		"""
		run freebayes
		return output directory
		
		## freebayes -p 1 -q 20 -m 20 --min-coverage 100 --min-alternate-fraction 0.01 --min-alternate-count 10 -V -f ../$input2 -b {} > {.}.vcf'
		"""
		temp_dir = os.path.join(self.utils.get_temp_dir())
		
		file_to_process = os.path.join(temp_dir, sample_name + ".bam" )
		cmd = "ln -s {} {}".format(bam_file, os.path.join(temp_dir, sample_name + ".bam" ))
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run freebayes")

		cmd = "ln -s {} {}.bai".format(bam_file + ".bai", file_to_process)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run freebayes")

		reference_fasta_temp = os.path.join(temp_dir, os.path.basename(reference_fasta))
		cmd = "ln -s {} {}".format(reference_fasta, reference_fasta_temp)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run freebayes")

		reference_fasta_fai = reference_fasta + FileExtensions.FILE_FAI
		if (not os.path.exists(reference_fasta_fai)):
			self.logger_production.error('Files doesnt exist: ' + reference_fasta_fai)
			self.logger_debug.error('Files doesnt exist: ' + reference_fasta_fai)
			raise Exception("Fail to run freebayes")
		
		reference_fasta_temp_fai = os.path.join(temp_dir, os.path.basename(reference_fasta) + FileExtensions.FILE_FAI)
		cmd = "ln -s {} {}".format(reference_fasta_fai, reference_fasta_temp_fai)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run freebayes")

		temp_file = self.utils.get_temp_file('freebayes_temp', '.vcf')
		cmd = "%s %s -f %s -b %s > %s" %\
 				(self.software_names.get_freebayes(), self.software_names.get_freebayes_parameters(),
 				reference_fasta_temp, file_to_process, temp_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run freebayes")
		
		if (os.path.getsize(temp_file) > 0):
			### run snpEff
			temp_file_2 = self.utils.get_temp_file("vcf_file", ".vcf")
			output_file = self.run_snpEff(reference_fasta, genbank_file, temp_file, os.path.join(temp_dir, os.path.basename(temp_file_2)))
			
			if output_file is None:	## sometimes the gff does not have amino sequences
				self.utils.copy_file(temp_file, os.path.join(temp_dir, os.path.basename(temp_file_2)))
				
			self.test_bgzip_and_tbi_in_vcf(os.path.join(temp_dir, os.path.basename(temp_file_2)))
			
			### add FREQ to vcf file
			vcf_file_out_temp = self.utils.add_freq_to_vcf(os.path.join(temp_dir, os.path.basename(temp_file_2)), os.path.join(temp_dir, sample_name + '.vcf'))
			os.unlink(temp_file)
			if (os.path.exists(temp_file_2)): os.unlink(temp_file_2)
			
			### pass vcf to tab
			self.run_snippy_vcf_to_tab(reference_fasta, genbank_file, vcf_file_out_temp, "{}.tab".format(os.path.join(temp_dir, sample_name)))
		return temp_dir

	def run_freebayes_parallel(self, bam_file, reference_fasta, genbank_file, sample_name):
		"""
		run freebayes
		return output directory
		
		## freebayes -p 1 -q 20 -m 20 --min-coverage 100 --min-alternate-fraction 0.01 --min-alternate-count 10 -V -f ../$input2 -b {} > {.}.vcf'
		"""
		temp_dir = os.path.join(self.utils.get_temp_dir())
		
		file_to_process = os.path.join(temp_dir, sample_name + ".bam" )
		cmd = "ln -s {} {}".format(bam_file, os.path.join(temp_dir, sample_name + ".bam" ))
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run freebayes parallel")

		cmd = "ln -s {} {}.bai".format(bam_file + ".bai", file_to_process)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run freebayes parallel")

		reference_fasta_temp = os.path.join(temp_dir, os.path.basename(reference_fasta))
		cmd = "ln -s {} {}".format(reference_fasta, reference_fasta_temp)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run freebayes parallel")

		reference_fasta_fai = reference_fasta + FileExtensions.FILE_FAI
		if (not os.path.exists(reference_fasta_fai)):
			self.logger_production.error('Files doesnt exist: ' + reference_fasta_fai)
			self.logger_debug.error('Files doesnt exist: ' + reference_fasta_fai)
			raise Exception("Fail to run freebayes parallel")
		
		reference_fasta_temp_fai = os.path.join(temp_dir, os.path.basename(reference_fasta) + FileExtensions.FILE_FAI)
		cmd = "ln -s {} {}".format(reference_fasta_fai, reference_fasta_temp_fai)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run freebayes parallel")
		
		### create regions
		temp_file_bam_coverage = self.utils.get_temp_file('bamtools_coverage', '.txt')
		
		# bamtools coverage -in aln.bam | coverage_to_regions.py ref.fa 500 >ref.fa.500.regions
		cmd = "%s coverage -in %s > %s" % (self.software_names.get_bamtools(), bam_file, temp_file_bam_coverage)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			os.unlink(temp_file_bam_coverage)
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run bamtools coverage")
	
		#### test if it has some data, if all zero does not have data
		if (os.path.getsize(temp_file_bam_coverage)) == 0:
			if (os.path.exists(temp_file_bam_coverage)): os.unlink(temp_file_bam_coverage)
			self.utils.remove_dir(temp_dir)
			return None
	
		temp_file_regions = self.utils.get_temp_file('freebayes_regions', '.txt')
		cmd = "cat %s | %s %s 500 > %s" %\
				(temp_file_bam_coverage, self.software_names.get_coverage_to_regions(),
				reference_fasta_temp_fai, temp_file_regions)

		exist_status = os.system(cmd)
		if (exist_status != 0):
			os.unlink(temp_file_bam_coverage)
			os.unlink(temp_file_regions)
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run freebayes parallel")
		
		temp_file = self.utils.get_temp_file('freebayes_temp', '.vcf')
		cmd = "%s %s %s %s -f %s -b %s > %s" %\
				(self.software_names.get_freebayes_parallel(), temp_file_regions, settings.THREADS_TO_RUN_SLOW,
				self.software_names.get_freebayes_parameters(), reference_fasta_temp, file_to_process, temp_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			os.unlink(temp_file_regions)
			os.unlink(temp_file_bam_coverage)
			os.unlink(temp_file)
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run freebayes parallel")
		
		### run snpEff
		if (os.path.exists(temp_file)):
			temp_file_2 = self.utils.get_temp_file("vcf_file", ".vcf")
			output_file = self.run_snpEff(reference_fasta, genbank_file, temp_file, os.path.join(temp_dir, os.path.basename(temp_file_2)))
			if output_file is None:	## sometimes the gff does not have amino sequences
				self.utils.copy_file(temp_file, os.path.join(temp_dir, os.path.basename(temp_file_2)))
				
			self.test_bgzip_and_tbi_in_vcf(os.path.join(temp_dir, os.path.basename(temp_file_2)))
		
			### add FREQ to vcf file
			vcf_file_out_temp = self.utils.add_freq_to_vcf(os.path.join(temp_dir, os.path.basename(temp_file_2)), os.path.join(temp_dir, sample_name + '.vcf'))
			### pass vcf to tab
			self.run_snippy_vcf_to_tab(reference_fasta, genbank_file, vcf_file_out_temp, "{}.tab".format(os.path.join(temp_dir, sample_name)))
		
		if (os.path.exists(temp_file)): os.unlink(temp_file)
		if (os.path.exists(temp_file_2)): os.unlink(temp_file_2)
		if (os.path.exists(temp_file_regions)): os.unlink(temp_file_regions)
		if (os.path.exists(temp_file_bam_coverage)): os.unlink(temp_file_bam_coverage)
		return temp_dir


	"""
	Global processing
	"""
#	@transaction.atomic
	def run_fastq_and_trimmomatic(self, sample, owner):
		"""
		run fastq and trimmomatic
		Upload average and sequence numbers
		"""
		manage_database = ManageDatabase()
		result_all = Result()
		
		### first try run downsize if necessary
		if (settings.DOWN_SIZE_FASTQ_FILES):
			(is_downsized, file_name_1, file_name_2) = self.make_downsize(sample.get_fastq(TypePath.MEDIA_ROOT, True),\
						sample.get_fastq(TypePath.MEDIA_ROOT, False), settings.MAX_FASTQ_FILE_UPLOAD)
			if (is_downsized):
				if (os.path.exists(file_name_1) and os.path.getsize(file_name_1) > 100):
					self.utils.move_file(file_name_1, sample.get_fastq(TypePath.MEDIA_ROOT, True))
				if (file_name_2 != None and len(file_name_2) > 0 and os.path.exists(file_name_2) and os.path.getsize(file_name_2) > 100): 
					self.utils.move_file(file_name_2, sample.get_fastq(TypePath.MEDIA_ROOT, False))
				self.utils.remove_dir(os.path.dirname(file_name_1))
				
				### set the downsize message
				manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_ALERT_DOWNSIZE_OF_FASTQ_FILES,\
											MetaKeyAndValue.META_VALUE_Success,\
											"Fastq files were down sized to values ~{}.".format( humanfriendly.format_size(int(settings.MAX_FASTQ_FILE_UPLOAD)) ))
		
		### first run fastqc
		try:
			temp_dir = self.run_fastq(sample.get_fastq(TypePath.MEDIA_ROOT, True), sample.get_fastq(TypePath.MEDIA_ROOT, False))
			result_all.add_software(SoftwareDesc(self.software_names.get_fastq_name(), self.software_names.get_fastq_version(), self.software_names.get_fastq_parameters()))
			
			### need to copy the files to samples/user path
			self.utils.copy_file(os.path.join(temp_dir, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, True))), sample.get_fastq_output(TypePath.MEDIA_ROOT, True))
			if (sample.exist_file_2()): self.utils.copy_file(os.path.join(temp_dir, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, False))), sample.get_fastq_output(TypePath.MEDIA_ROOT, False))
		except Exception as e:
			result = Result()
			result.set_error("Fail to run fastq software: " + e.args[0])
			result.add_software(SoftwareDesc(self.software_names.get_fastq_name(), self.software_names.get_fastq(), self.software_names.get_fastq_parameters()))
			manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			self.utils.remove_dir(temp_dir)
			return False
		self.utils.remove_dir(temp_dir)
		
		### run trimmomatic
		try:
			(temp_dir, parameters) = self.run_trimmomatic(sample.get_fastq(TypePath.MEDIA_ROOT, True), sample.get_fastq(TypePath.MEDIA_ROOT, False), sample.name, owner)
			result_all.add_software(SoftwareDesc(self.software_names.get_trimmomatic_name(), self.software_names.get_trimmomatic_version(), parameters))
			### need to copy the files to samples/user path
			self.utils.copy_file(os.path.join(temp_dir, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True))), sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True))
			if (sample.exist_file_2()): self.utils.copy_file(os.path.join(temp_dir, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False))), sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False))
		except Exception as e:
			result = Result()
			result.set_error("Fail to run trimmomatic software: " + e.args[0])
			result.add_software(SoftwareDesc(self.software_names.get_trimmomatic_name(), self.software_names.get_trimmomatic_version(), self.software_names.get_trimmomatic_parameters()))
			manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			self.utils.remove_dir(temp_dir)
			return False
		self.utils.remove_dir(temp_dir)
		
		
		### run fastq again
		try:
			temp_dir = self.run_fastq(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True), sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False))
											
			### need to copy the files to samples/user path
			self.utils.copy_file(os.path.join(temp_dir, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True))), sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True))
			if (sample.exist_file_2()): self.utils.copy_file(os.path.join(temp_dir, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False))), sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False))
		except Exception as e:
			result = Result()
			result.set_error("Fail to run fastq software: " + e.args[0])
			result.add_software(SoftwareDesc(self.software_names.get_fastq_name(), self.software_names.get_fastq(), self.software_names.get_fastq_parameters()))
			manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			self.utils.remove_dir(temp_dir)
			return False
		self.utils.remove_dir(temp_dir)

		### collect numbers
		(lines_1, average_1) = self.get_lines_and_average_reads(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True))
		if (sample.exist_file_2()): (lines_2, average_2) = self.get_lines_and_average_reads(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False))
		else: (lines_2, average_2) = (None, None)
		result_average = ResultAverageAndNumberReads(lines_1, average_1, lines_2, average_2)
		manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, MetaKeyAndValue.META_VALUE_Success, result_average.to_json())

		## save everything OK
		manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, MetaKeyAndValue.META_VALUE_Success, "Success, Fastq(%s), Trimmomatic(%s)" %\
							(self.software_names.get_fastq_version(), self.software_names.get_trimmomatic_version()))
		manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software, MetaKeyAndValue.META_VALUE_Success, result_all.to_json())

		### set the flag of the end of the task		
		meta_sample = manage_database.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue)
		if (meta_sample != None):
			manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Success, meta_sample.description)
		return result_average.has_reads()

	"""
	Global processing, fastQ, trimmomatic and GetSpecies
	"""
	def run_fastq_and_trimmomatic_and_identify_species(self, sample, user):

		print("Start ProcessControler")
		process_controler = ProcessControler()
		process_SGE = ProcessSGE()
		process_SGE.set_process_controler(user, process_controler.get_name_sample(sample), ProcessControler.FLAG_RUNNING)
		
		### it can be deleted
		if (sample.is_deleted or not sample.is_valid_1):
			process_SGE.set_process_controler(user, process_controler.get_name_sample(sample), ProcessControler.FLAG_FINISHED)
			return True

		try:
			print("Start run_fastq_and_trimmomatic")
			### run trimmomatics
			b_return = self.run_fastq_and_trimmomatic(sample, user)
			
			print("Result run_fastq_and_trimmomatic: " + str(b_return))
			
			### queue the quality check and
			if (b_return):	## don't run for single file because spades doesn't work for one single file
				self.identify_type_and_sub_type(sample, sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True),\
					sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False), user)
	
			## set the flag that is ready for process
			if (b_return):
				sample_to_update = Sample.objects.get(pk=sample.id)
				sample_to_update.is_ready_for_projects = True
				sample_to_update.type_subtype = sample_to_update.get_type_sub_type()
				
				(tag_mixed_infection, alert, message) = sample_to_update.get_mixed_infection()
				if (sample_to_update.number_alerts == None): sample_to_update.number_alerts = alert
				else: sample_to_update.number_alerts += alert
					
				manage_database = ManageDatabase()
				if (message != None and len(message) > 0):
					manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE,\
										MetaKeyAndValue.META_VALUE_Success, message)
	
				### save tag mixed_infecion
				manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_TAG_MIXED_INFECTION_TYPE_SUBTYPE,\
							MetaKeyAndValue.META_VALUE_Success, tag_mixed_infection)

				try:
					mixed_infections_tag = MixedInfectionsTag.objects.get(name=tag_mixed_infection)
				except MixedInfectionsTag.DoesNotExist as e:
					mixed_infections_tag = MixedInfectionsTag()
					mixed_infections_tag.name = tag_mixed_infection
					mixed_infections_tag.save()
				
				sample_to_update.mixed_infections_tag = mixed_infections_tag
				sample_to_update.save()
			else:
				sample_to_update = Sample.objects.get(pk=sample.id)
				manage_database = ManageDatabase()
				manage_database.set_sample_metakey(sample_to_update, user, MetaKeyAndValue.META_KEY_ALERT_NO_READS_AFTER_FILTERING,\
										MetaKeyAndValue.META_VALUE_Success, "Warning: no reads left after filtering.")
				
				if (sample_to_update.number_alerts == None): sample_to_update.number_alerts = 1
				else: sample_to_update.number_alerts += 1
				sample_to_update.type_subtype = Constants.EMPTY_VALUE_TYPE_SUBTYPE
				sample_to_update.save()
				
		except:
			process_SGE.set_process_controler(user, process_controler.get_name_sample(sample), ProcessControler.FLAG_ERROR)
			return False
		
		### finished
		process_SGE.set_process_controler(user, process_controler.get_name_sample(sample), ProcessControler.FLAG_FINISHED)
		return b_return


	def process_second_stage_snippy_coverage_freebayes(self, project_sample, user):
		"""
		"""
		### make it running 
		process_controler = ProcessControler()
		process_SGE = ProcessSGE()
		process_SGE.set_process_controler(user, process_controler.get_name_project_sample(project_sample), ProcessControler.FLAG_RUNNING)
		
		## run collect data
		return self.__process_second_stage_snippy_coverage_freebayes(project_sample, user)
		
	"""
	Global processing, Snippy, Coverage, Freebayes and MixedInfections
	"""
#	@transaction.atomic
	def __process_second_stage_snippy_coverage_freebayes(self, project_sample, user):
		"""
		Global processing, snippy, coverage, 
		"""
		
		process_controler = ProcessControler()
		process_SGE = ProcessSGE()
		manageDatabase = ManageDatabase()
		result_all = Result()
		### metakey for this process
		metaKeyAndValue = MetaKeyAndValue()
		try:	
			meta_key_project_sample = metaKeyAndValue.get_meta_key_queue_by_project_sample_id(project_sample.id)
	
			### Test if this sample already run		
			meta_sample = manageDatabase.get_project_sample_metakey_last(project_sample, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue)
			if (meta_sample != None and meta_sample.value == MetaKeyAndValue.META_VALUE_Success): return 
	
			## process snippy
			try:
				out_put_path = self.run_snippy(project_sample.sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True),\
						project_sample.sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False),\
						project_sample.project.reference.get_reference_gbk(TypePath.MEDIA_ROOT),\
						project_sample.sample.name)
				result_all.add_software(SoftwareDesc(self.software_names.get_snippy_name(), self.software_names.get_snippy_version(), self.software_names.get_snippy_parameters()))
			except Exception as e:
				result = Result()
				result.set_error(e.args[0])
				result.add_software(SoftwareDesc(self.software_names.get_snippy_name(), self.software_names.get_snippy_version(), self.software_names.get_snippy_parameters()))
				manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Snippy, MetaKeyAndValue.META_VALUE_Error, result.to_json())
				
				### get again and set error
				project_sample = ProjectSample.objects.get(pk=project_sample.id)
				project_sample.is_error = True
				project_sample.save()
				
				meta_sample = manageDatabase.get_project_sample_metakey(project_sample, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue)
				if (meta_sample != None):
					manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Error, meta_sample.description)
				process_SGE.set_process_controler(user, process_controler.get_name_project_sample(project_sample), ProcessControler.FLAG_ERROR)
				return False
	
			## copy the files to the project sample directories
			self.copy_files_to_project(project_sample, self.software_names.get_snippy_name(), out_put_path)
			remove_path = os.path.dirname(out_put_path)
			if (len(remove_path.split('/')) > 2): self.utils.remove_dir(remove_path)
			else: self.utils.remove_dir(out_put_path)
	
			### make the link for the new tab file name
			path_snippy_tab = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB,  self.software_names.get_snippy_name())
			if (os.path.exists(path_snippy_tab)):
				sz_file_to = project_sample.get_file_output_human(TypePath.MEDIA_ROOT, FileType.FILE_TAB,  self.software_names.get_snippy_name())
				self.utils.link_file(path_snippy_tab, sz_file_to)
			
			## get coverage from deep file
			get_coverage = GetCoverage()
			try:
				coverage = get_coverage.get_coverage(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ,\
							self.software_names.get_snippy_name()), project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT))

				################################
				##################################
				### set the alerts in the coverage
				project_sample = ProjectSample.objects.get(pk=project_sample.id)
				for element in coverage.get_dict_data():
					if (not coverage.is_100_more_9(element)):
						project_sample.alert_second_level += 1
						meta_key = metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_ALERT_COVERAGE_9, element)
						manageDatabase.set_project_sample_metakey(project_sample, user, meta_key, MetaKeyAndValue.META_VALUE_Success, coverage.get_fault_message_9(element))
					elif (not coverage.is_100_more_0(element)):
						project_sample.alert_first_level += 1
						meta_key = metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_ALERT_COVERAGE_0, element)
						manageDatabase.set_project_sample_metakey(project_sample, user, meta_key, MetaKeyAndValue.META_VALUE_Success, coverage.get_fault_message_0(element))
				project_sample.save()
			
				## set the coverage in database
				meta_sample = manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage,\
									MetaKeyAndValue.META_VALUE_Success, coverage.to_json())
			except Exception as e:
				result = Result()
				result.set_error("Fail to get coverage: " + e.args[0])
				result.add_software(SoftwareDesc(self.software_names.get_coverage_name(), self.software_names.get_coverage_version(), self.software_names.get_coverage_parameters()))
				manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Error, result.to_json())
				
				### get again and set error
				project_sample = ProjectSample.objects.get(pk=project_sample.id)
				project_sample.is_error = True
				project_sample.save()
				
				meta_sample = manageDatabase.get_project_sample_metakey(project_sample, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue)
				if (meta_sample != None):
					manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Error, meta_sample.description)
				process_SGE.set_process_controler(user, process_controler.get_name_project_sample(project_sample), ProcessControler.FLAG_ERROR)
				return False
			
			## identify VARIANTS IN INCOMPLETE LOCUS in all locus, set yes in variants if are in ares with coverage problems
			parse_out_files = ParseOutFiles()
			parse_out_files.add_variants_in_incomplete_locus(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, 
										self.software_names.get_snippy_name()), coverage)
			
			## run freebayes if at least one segment has some coverage
			try:
				out_put_path = self.run_freebayes_parallel(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_BAM, self.software_names.get_snippy_name()),\
							project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT), project_sample.project.reference.get_reference_gbk(TypePath.MEDIA_ROOT),\
							project_sample.sample.name)
				result_all.add_software(SoftwareDesc(self.software_names.get_freebayes_name(), self.software_names.get_freebayes_version(), self.software_names.get_freebayes_parameters()))
			except Exception as e:
				
				### can fail the freebayes parallel and try the regular one
				try:
					out_put_path = self.run_freebayes(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_BAM, self.software_names.get_snippy_name()),\
								project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT), project_sample.project.reference.get_reference_gbk(TypePath.MEDIA_ROOT),\
								project_sample.sample.name)
					result_all.add_software(SoftwareDesc(self.software_names.get_freebayes_name(), self.software_names.get_freebayes_version(), self.software_names.get_freebayes_parameters()))
				except Exception as e:
					result = Result()
					result.set_error(e.args[0])
					result.add_software(SoftwareDesc(self.software_names.get_freebayes_name(), self.software_names.get_freebayes_version(), self.software_names.get_freebayes_parameters()))
					manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Freebayes, MetaKeyAndValue.META_VALUE_Error, result.to_json())
					
					### get again and set error
					project_sample = ProjectSample.objects.get(pk=project_sample.id)
					project_sample.is_error = True
					project_sample.save()
					
					meta_sample = manageDatabase.get_project_sample_metakey(project_sample, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue)
					if (meta_sample != None):
						manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Error, meta_sample.description)
					process_SGE.set_process_controler(user, process_controler.get_name_project_sample(project_sample), ProcessControler.FLAG_ERROR)
					return False
			
			## count hits from tab file
			count_hits = CountHits()
			if (not out_put_path is None):
				file_tab = os.path.join(out_put_path, project_sample.sample.name + ".tab")
				if (os.path.exists(file_tab)):
					vect_count_type = ['snp']	## only detects snp
					count_hits = self.utils.count_hits_from_tab(file_tab, vect_count_type)
					### set flag that is finished
					manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Count_Hits, MetaKeyAndValue.META_VALUE_Success, count_hits.to_json())
				else:
					result = Result()
					result.set_error("Fail to collect tab file from freebayes")
					result.add_software(SoftwareDesc(self.software_names.get_freebayes_name(), self.software_names.get_freebayes_version(), self.software_names.get_freebayes_parameters()))
					manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Freebayes, MetaKeyAndValue.META_VALUE_Error, result.to_json())
					
					### get again and set error
					project_sample = ProjectSample.objects.get(pk=project_sample.id)
					project_sample.is_error = True
					project_sample.save()
					
					meta_sample = manageDatabase.get_project_sample_metakey(project_sample, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue)
					if (meta_sample != None):
						manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Error, meta_sample.description)
					process_SGE.set_process_controler(user, process_controler.get_name_project_sample(project_sample), ProcessControler.FLAG_ERROR)
					return False
				
				self.copy_files_to_project(project_sample, self.software_names.get_freebayes_name(), out_put_path)
				## remove path dir if exist
				self.utils.remove_dir(out_put_path)
			else:
				### set count hits to zero
				manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Count_Hits, MetaKeyAndValue.META_VALUE_Success, count_hits.to_json())
			
			### draw coverage
			try:
				### make the coverage images
				draw_all_coverage = DrawAllCoverage()
				draw_all_coverage.draw_all_coverages(project_sample)
			except:
				result = Result()
				result.set_error("Fail to draw coverage images")
				result.add_software(SoftwareDesc('In house software', '1.0', ''))
				manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Error, result.to_json())
				
				### get again and set error
				project_sample = ProjectSample.objects.get(pk=project_sample.id)
				project_sample.is_error = True
				project_sample.save()
				
				meta_sample = manageDatabase.get_project_sample_metakey(project_sample, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue)
				if (meta_sample != None):
					manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Error, meta_sample.description)
				process_SGE.set_process_controler(user, process_controler.get_name_project_sample(project_sample), ProcessControler.FLAG_ERROR)
				return False

			### mixed infection
			try:
				## get instances
				mixed_infections_management = MixedInfectionsManagement()
				
				## set the alert also
				mixed_infection = mixed_infections_management.get_mixed_infections(project_sample, user, count_hits)
			except:
				result = Result()
				result.set_error("Fail to calculate mixed infextion")
				result.add_software(SoftwareDesc('In house software', '1.0', ''))
				manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Mixed_Infection, MetaKeyAndValue.META_VALUE_Error, result.to_json())
				
				### get again and set error
				project_sample = ProjectSample.objects.get(pk=project_sample.id)
				project_sample.is_error = True
				project_sample.save()
				
				meta_sample = manageDatabase.get_project_sample_metakey(project_sample, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue)
				if (meta_sample != None):
					manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Error, meta_sample.description)
				process_SGE.set_process_controler(user, process_controler.get_name_project_sample(project_sample), ProcessControler.FLAG_ERROR)
				return False
			
			### get again
			manage_database = ManageDatabase()
			project_sample = ProjectSample.objects.get(pk=project_sample.id)
			project_sample.is_finished = True
			project_sample.is_deleted = False
			project_sample.is_deleted_in_file_system = False
			project_sample.date_deleted = None
			project_sample.is_error = False
			project_sample.count_variations = manage_database.get_variation_count(count_hits)
			project_sample.mixed_infections = mixed_infection
			project_sample.save()
			
			from utils.collect_extra_data import CollectExtraData
			collect_extra_data = CollectExtraData()
			
			### get a clean freebayes file
			tab_freebayes_file = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)
			if (os.path.exists(tab_freebayes_file)):
				file_out = project_sample.get_file_output_human(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)
				collect_extra_data.collect_variations_freebayes_only_one_file(tab_freebayes_file, file_out)
			
			### get clean consensus file
			consensus_fasta = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, SoftwareNames.SOFTWARE_SNIPPY_name)
			if (os.path.exists(consensus_fasta)):
				file_out = project_sample.get_consensus_file(TypePath.MEDIA_ROOT)
				self.utils.filter_fasta_all_sequences_file(consensus_fasta, coverage, file_out, False)
			
			### set the tag of result OK 
			manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Snippy_Freebayes, MetaKeyAndValue.META_VALUE_Success, result_all.to_json())
			
			### set the flag of the end of the task		
			meta_sample = manageDatabase.get_project_sample_metakey_last(project_sample, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue)
			if (meta_sample != None):
				manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Success, meta_sample.description)
		except:
			## finished with error
			process_SGE.set_process_controler(user, process_controler.get_name_project_sample(project_sample), ProcessControler.FLAG_ERROR)
			return False
		
		### finished
		process_SGE.set_process_controler(user, process_controler.get_name_project_sample(project_sample), ProcessControler.FLAG_FINISHED)
		return True

	def create_index_files(self, file_name):
		"""
		create index, need to be .bz
		"""
		file_to_index = file_name
		if (not file_to_index.endswith(FileExtensions.FILE_VCF_GZ)): file_to_index += FileExtensions.FILE_VCF_GZ
		if (not os.path.exists(file_to_index)):
			self.logger_production.error("File doesn't exist: " + file_to_index)
			self.logger_debug.error("Fail doesn't exist: " + file_to_index)
			raise Exception("File doesn't exist")
		
		## test if tbi exists
		if (os.path.exists(file_to_index + FileExtensions.FILE_TBI)): return

		cmd = "{} {}".format(self.software_names.get_tabix(), file_name)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to create index") 
	
	
	def create_index_files_from_igv_tools(self, file_name):
		"""
		Create index from igvtools
		"""
		cmd = "java -jar {} index {}".format(self.software_names.get_igvtools(), file_name)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to create index") 


	def dos_2_unix(self, file_name):
		"""
		convert dos 2 unix
		"""
		if (not os.path.exists(file_name)): return
		cmd = "dos2unix {}".format(file_name)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to create index") 
	

	def make_downsize(self, path_1, path_2, max_fastq_file):
		"""
		Reduce file size from X to Constants.MAX_FASTQ_FILE
		return (TRUE|FALSE if was reduce or not, new_path_1_reduced, new_path_2_reduced)
		"""
		file_size_max = 0
		if (os.path.exists(path_1)):
			file_size_max = os.path.getsize(path_1)
		if (path_2 != None and os.path.exists(path_2) and file_size_max < os.path.getsize(path_2)): 
			file_size_max = os.path.getsize(path_2)
			
		### need to make the down size
		if (file_size_max > max_fastq_file):
			ratio = max_fastq_file / file_size_max
		
			## smmall correction
			if (ratio < 0.8): ratio += 0.1

			path_to_work = self.utils.get_temp_dir()
			path_1_temp = self.utils.get_temp_file_from_dir(path_to_work, 'fastq_1', '.fastq')
			path_2_temp = self.utils.get_temp_file_from_dir(path_to_work, 'fastq_2', '.fastq')
			
			self.utils.uncompress_files(self.software_names.get_bgzip(), path_1, path_1_temp)
			file_names = path_1_temp
			if (path_2 != None and os.path.exists(path_2)):
				self.utils.uncompress_files(self.software_names.get_bgzip(), path_2, path_2_temp)
				file_names += " " + path_2_temp

			cmd = "{} -p {:.2f} -o {}/sample {}".format(self.software_names.get_fastqtools_sample(), ratio, path_to_work, file_names)
#			print(cmd)
			exist_status = os.system(cmd)
			if (exist_status != 0):
				self.logger_production.error('Fail to run: ' + cmd)
				self.logger_debug.error('Fail to run: ' + cmd)
				raise Exception("Fail to downsize") 
		
			if (path_2 != None and os.path.exists(path_2)):
				self.utils.compress_files(self.software_names.get_bgzip(), os.path.join(path_to_work, 'sample.1.fastq'))
				path_1_temp = os.path.join(path_to_work, 'sample.1.fastq' + FileExtensions.FILE_GZ)
				self.utils.compress_files(self.software_names.get_bgzip(), os.path.join(path_to_work, 'sample.2.fastq'))
				path_2_temp = os.path.join(path_to_work, 'sample.2.fastq' + FileExtensions.FILE_GZ)
			else:
				self.utils.compress_files(self.software_names.get_bgzip(), os.path.join(path_to_work, 'sample.fastq'))
				path_1_temp = os.path.join(path_to_work, 'sample.fastq' + FileExtensions.FILE_GZ)
				
			return(True, path_1_temp, path_2_temp if (path_2 != None and os.path.exists(path_2)) else None)
		return(False, path_1, path_2)


class Contigs2Sequences(object):
	'''
	classdocs
	'''
	DATABASE_NAME = "influenza_assign_segments2contigs"
	
	utils = Utils()

	def __init__(self, b_testing):
		'''
		Constructor
		'''
		self.b_testing = b_testing
		self.root_path = os.path.join(getattr(settings, "STATIC_ROOT", None), "tests" if b_testing else "")

	def get_most_recent_database(self):
		"""
		
		"""
		version = 0
		path_to_return = None
		path_to_find = os.path.join(self.root_path, Constants.DIR_TYPE_CONTIGS_2_SEQUENCES)
		for file in self.utils.get_all_files(path_to_find):
			base_name = os.path.basename(file)
			match = re.search('\w+(_[v|V]\d+)\.\w+', base_name)
			if (match == None): continue
			if (len(match.regs) == 2):
				try:
					self.utils.is_fasta(os.path.join(path_to_find, file))
				except IOError as e:
					continue
				
				(temp, path) = (int(match.group(1).lower().replace("_v", "")), os.path.join(path_to_find, file))
				if (temp > version):
					version = temp
					path_to_return = path
		return (str(version), path_to_return)

	def get_database_name(self):
		### get database file name
		(version, database_file_name) = self.get_most_recent_database()
		return self.utils.clean_extension(os.path.basename(database_file_name))


	def identify_contigs(self, file_name, file_name_out):
		"""
		identify contigs
		params in: fasta file from spades
		out: fasta file with low coverage removed and Elements ID in description
		"""
		software = Software()
		
		### get database file name
		(version, database_file_name) = self.get_most_recent_database()
		dabasename = self.get_database_name()
		
		### first create database
		if (not software.is_exist_database_abricate(dabasename)):
			software.create_database_abricate(dabasename, database_file_name)
		
		out_file = self.utils.get_temp_file('abricate_contig2seq', FileExtensions.FILE_TXT)
		### run abricate
		software.run_abricate(dabasename, file_name, SoftwareNames.SOFTWARE_ABRICATE_PARAMETERS_mincov_30, out_file)

		parseOutFiles = ParseOutFiles()
		(vect_data, clean_abricate_file) = parseOutFiles.parse_abricate_file(out_file, file_name_out, SoftwareNames.SOFTWARE_SPAdes_CLEAN_HITS_BELLOW_VALUE)
		
		vect_out_fasta = []
		vect_out_fasta_without_id = []
		out_file_fasta = self.utils.get_temp_file('spades_out_identified', FileExtensions.FILE_FASTA)
		with open(file_name) as handle_in, open(out_file_fasta, 'w') as handle_out:
			for record in SeqIO.parse(handle_in, Constants.FORMAT_FASTA):
				vect_possible_id = []
				for dict_data in vect_data:
					if (dict_data['Seq_Name'] == record.name): vect_possible_id.append(dict_data['Gene'])
				if (len(vect_possible_id)):
					vect_out_fasta.append(SeqRecord(Seq(str(record.seq)), id = "_".join(record.id.split('.')[0].split('_')[:4]),\
												description=";".join(vect_possible_id)))
				elif (float(record.id.split('_')[-1]) > SoftwareNames.SOFTWARE_SPAdes_CLEAN_HITS_BELLOW_VALUE):
					vect_out_fasta_without_id.append(SeqRecord(Seq(str(record.seq)), id = record.id, description=""))

			if (len(vect_out_fasta) > 0 or len(vect_out_fasta_without_id) > 0):
				vect_out_fasta.extend(vect_out_fasta_without_id)
				SeqIO.write(vect_out_fasta, handle_out, "fasta")
		
		if (os.path.exists(out_file)): os.unlink(out_file)
		return (out_file_fasta, clean_abricate_file)

	
			
