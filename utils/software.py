'''
Created on Oct 28, 2017

@author: mmp
'''
import os, gzip
import logging
import cmd
import subprocess
from utils.coverage import DrawAllCoverage
from utils.utils import Utils
from utils.parse_out_files import ParseOutFiles
from constants.constants import Constants, TypePath, FileType, FileExtensions
from constants.meta_key_and_values import MetaKeyAndValue
from manage_virus.models import UploadFile
from managing_files.models import Sample, ProjectSample
from manage_virus.uploadFiles import UploadFiles
from managing_files.manage_database import ManageDatabase
from utils.result import Result, SoftwareDesc, ResultAverageAndNumberReads
from utils.parse_coverage_file import GetCoverage
from django.db import transaction
from constants.software_names import SoftwareNames
from constants.tag_names_constants import TagNamesConstants
from Bio import SeqIO

class Software(object):
	'''
	classdocs
	'''
	utils = Utils()
	software_names = SoftwareNames()
	CORES_TO_USE = 3
	
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
			self.utils.create_index_files(self.software_names.get_tabix(), vcf_file)
		else:
			self.utils.create_index_files(self.software_names.get_tabix(), vcf_file + ".gz")


	def get_vect_type_files_to_copy(self, software):
		"""
		get type of files to copy
		"""
		if (software == SoftwareNames.SOFTWARE_SNIPPY_name):
			return [FileType.FILE_BAM, FileType.FILE_BAM_BAI, FileType.FILE_CONSENSUS_FA, FileType.FILE_DEPTH_GZ, FileType.FILE_DEPTH_GZ_TBI,\
				FileType.FILE_TAB, FileType.FILE_VCF_GZ, FileType.FILE_VCF_GZ_TBI, FileType.FILE_CSV]
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
			elif (type_file == FileType.FILE_VCF):	## vcf file
				### create the gzip file
				self.utils.compress_files(self.software_names.get_bgzip(), os.path.join(path_from, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, type_file, software))))
				### create the tabix
				self.utils.create_index_files(self.software_names.get_tabix(), os.path.join(path_from, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, software))))
				
				### copy both
				self.utils.copy_file(os.path.join(path_from, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, software))),\
					project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, software))
				self.utils.copy_file(os.path.join(path_from, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ_TBI, software))),\
					project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ_TBI, software))
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


	def run_spades(self, fastq_1, fastq_2, out_dir):
		"""
		Run spades
		"""
		if (fastq_2 is None or len(fastq_2) == 0): cmd = "%s -s %s --meta --only-assembler -t %d -o %s" % (self.software_names.get_spades(), fastq_1, self.CORES_TO_USE, out_dir)
		else: cmd = "%s --pe1-1 %s --pe1-2 %s --meta --only-assembler -t %d -o %s" % (self.software_names.get_spades(), fastq_1, fastq_2, self.CORES_TO_USE, out_dir)
		
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
		cmd = "%s --list" % (self.software_names.get_abricate())
		proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
		(out, err) = proc.communicate()
		if (err != None):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run abricate")
		out_str = out.decode("utf-8")
		for line in out_str.split('\n'):
			for field in line.split('\t'):
				if (field.lower() == database.lower()): return True
		return False


	def create_database_abricate(self, database, file_name):
		"""
		create a database
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
		

	def run_abricate(self, database, file_name, out_file):
		"""
		Run abricator
		"""
		cmd = "%s --db %s --quiet %s > %s" % (self.software_names.get_abricate(), database, file_name, out_file)
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
	@transaction.atomic
	def identify_type_and_sub_type(self, sample, fastq1_1, fastq1_2, owner):
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
			result_all.add_software(SoftwareDesc(self.software_names.get_spades_name(), self.software_names.get_spades_version(), self.software_names.get_spades_parameters()))
		except Exception:
			result = Result()
			result.set_error("Spades (%s) fail to run" % (self.software_names.get_spades_version()))
			result.add_software(SoftwareDesc(self.software_names.get_spades_name(), self.software_names.get_spades_version(), self.software_names.get_spades_parameters()))
			manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False
		
		file_out = os.path.join(out_dir_spades, "contigs.fasta")
		if (not os.path.exists(file_out) or os.path.getsize(file_out) < 100):
			## save error in MetaKeySample
			result = Result()
			result.set_error("Spades (%s) fail to run" % (self.software_names.get_spades_version()))
			result.add_software(SoftwareDesc(self.software_names.get_spades_name(), self.software_names.get_spades_version(), self.software_names.get_spades_parameters()))
			manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False
		
		try:
			uploadFile = UploadFile.objects.order_by('-version')[0]
		except UploadFile.DoesNotExist:
			## save error in MetaKeySample
			result = Result()
			result.set_error("Abricate (%s) fail to run" % (self.software_names.get_abricate_version()))
			result.add_software(SoftwareDesc(self.software_names.get_abricate_name(), self.software_names.get_abricate_version(), self.software_names.get_abricate_parameters()))
			manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False

		if (not self.is_exist_database_abricate(uploadFile.abricate_name)):
			try:
				self.create_database_abricate(uploadFile.abricate_name, uploadFile.path)
			except Exception:
				result = Result()
				result.set_error("Abricate (%s) fail to run --setupdb" % (self.software_names.get_abricate_version()))
				result.add_software(SoftwareDesc(self.software_names.get_abricate_name(), self.software_names.get_abricate_version(), self.software_names.get_abricate_parameters()))
				manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
				cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
				return False
		
		## run abricate
		out_file_abricate = self.utils.get_temp_file("temp_abricate", ".txt")
		try:
			cmd = self.run_abricate(uploadFile.abricate_name, file_out, out_file_abricate)
			result_all.add_software(SoftwareDesc(self.software_names.get_abricate_name(), self.software_names.get_abricate_version(), self.software_names.get_abricate_parameters()))
		except Exception:
			result = Result()
			result.set_error("Abricate (%s) fail to run" % (self.software_names.get_abricate_version()))
			result.add_software(SoftwareDesc(self.software_names.get_abricate_name(), self.software_names.get_abricate_version(), self.software_names.get_abricate_parameters()))
			manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False
		
		if (not os.path.exists(out_file_abricate)):
			## save error in MetaKeySample
			result = Result()
			result.set_error("Abricate (%s) fail to run" % (self.software_names.get_abricate_version()))
			result.add_software(SoftwareDesc(self.software_names.get_abricate(), self.software_names.get_abricate_version(), self.software_names.get_abricate_parameters()))
			manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False

		parseOutFiles = ParseOutFiles()
		vect_data = parseOutFiles.parse_abricate_file(out_file_abricate)
		## copy the abricate output 
		self.utils.copy_file(out_file_abricate, sample.get_abricate_output(TypePath.MEDIA_ROOT))
		
		uploadFiles = UploadFiles()
		vect_data = uploadFiles.uploadIdentifyVirus(vect_data, uploadFile.abricate_name)
		if (len(vect_data) == 0):
			## save error in MetaKeySample
			result = Result()
			result.set_error("Fail to identify type and sub type")
			result.add_software(SoftwareDesc(self.software_names.get_abricate(), self.software_names.get_abricate_version(), self.software_names.get_abricate_parameters()))
			manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm %s" % (out_file_abricate); os.system(cmd)
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False
		
		
		for identify_virus in vect_data:
			sample.identify_virus.add(identify_virus)
		sample.save()
		
		## save everything OK
		manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Success, "Success, Spades(%s), Abricate(%s)" % (self.software_names.get_spades_version(), self.software_names.get_abricate_version()))
		manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample_Software, MetaKeyAndValue.META_VALUE_Success, result_all.to_json())
		cmd = "rm %s" % (out_file_abricate); os.system(cmd)
		self.utils.remove_dir(out_dir_spades)
		return True


	def run_fastq(self, file_name_1, file_name_2):
		"""
		run fastQ, return output directory
		-o OUT_FOLDER --nogroup --format fastq --threads 10 --dir OUT_FOLDER FILE1 FILE2
		"""
		temp_dir = self.utils.get_temp_dir()
		if (not file_name_2 is None and len(file_name_2) > 0):
			cmd = "%s -o %s --nogroup --format fastq --threads %d --dir %s %s %s" % (self.software_names.get_fastq(), temp_dir, Software.CORES_TO_USE, 
										temp_dir, file_name_1, file_name_2)
		else: 
			cmd = "%s -o %s --nogroup --format fastq --threads %d --dir %s %s" % (self.software_names.get_fastq(), temp_dir, Software.CORES_TO_USE, temp_dir, file_name_1)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run run_fastq")
		return temp_dir


	def run_prokka(self, fasta_file_name):
		"""
		run prokka, in fasta file
		out: genbank
		{bpipe_prokka} FILE1 --kingdom Viruses --locustag locus --kingdom Viruses --locustag locus --genus Influenzavirus 
			--species Influenzavirus --strain ref_PREFIX_FILES_OUT --outdir OUT_FOLDER/PREFIX_FILES_OUT --prefix PREFIX_FILES_OUT
		"""
		if (not os.path.exists(fasta_file_name)): raise Exception("File doesn't exist")
		temp_dir = self.utils.get_temp_dir()
		name_strain = os.path.basename(fasta_file_name)
		name_strain = name_strain[:name_strain.rfind('.')]
		cmd = "{} {} {} --strain {} --force --outdir {} --prefix {}".format(\
					self.software_names.get_prokka(), fasta_file_name, self.software_names.get_prokka_parameters(), name_strain, temp_dir, name_strain)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run prokka")
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
		temp_file = self.utils.get_temp_file('clean_fasta_names', FileExtensions.FILE_FASTA)
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

	def run_mafft(self, input_file, out_file):
		"""
		run mafft
		out: out_file
		"""
		cmd = "{}; {} {} > {}".format(self.software_names.get_mafft_set_env_variable(), self.software_names.get_mafft(),\
							input_file, out_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run mafft")
		return out_file
	
	def run_fasttree(self, input_file, out_file):
		"""
		run fasttree
		out: out_file
		"""
		cmd = "{} {} {} > {}".format(self.software_names.get_fasttree(), self.software_names.get_fasttree_parameters(), input_file, out_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run progressive mauve")
		return out_file
	
	def run_trimmomatic(self, file_name_1, file_name_2, sample_name):
		"""
		run trimmomatic
		return output directory
		
		#${bpipe_trimmomatic} PE -threads 3 -basein FILE1 -baseout OUT_FOLDER/PREFIX_FILES_OUT.fastq.gz SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:55 TOPHRED33
		
		PE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
   			or: 
		SE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] <inputFile> <outputFile> <trimmer1>
		"""

		temp_dir = self.utils.get_temp_dir()
		if (file_name_2 is None or len(file_name_2) == 0):
			cmd = "java -jar %s SE -threads %d %s %s_1P.fastq.gz %s" % (self.software_names.get_trimmomatic(), Software.CORES_TO_USE, file_name_1, 
					os.path.join(temp_dir, sample_name), self.software_names.get_trimmomatic_parameters())
		else:
			cmd = "java -jar %s PE -threads %d -basein %s -baseout %s.fastq.gz %s" % (self.software_names.get_trimmomatic(), Software.CORES_TO_USE, 
										file_name_1, os.path.join(temp_dir, sample_name), self.software_names.get_trimmomatic_parameters())
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run trimmomatic")
		return temp_dir

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
				(self.software_names.get_snippy(), Software.CORES_TO_USE, temp_dir, sample_name,
				path_reference, self.software_names.get_snippy_parameters(), file_name_1)
		else:
			cmd = "%s --cpus %d --outdir %s --prefix %s --ref %s %s --R1 %s --R2 %s" %\
				(self.software_names.get_snippy(), Software.CORES_TO_USE, temp_dir, sample_name,
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
		bp_genbank2gff3 --filter gene --filter region --outdir stdout --quiet static/tests/managing_files/A_H3N2_A_Hong_Kong_4801_2014.gbk
		"""
		vect_filter = ['gene', 'region']
		temp_file = self.utils.get_temp_file("gbk_to_gff3", ".txt") 
		cmd = "%s --filter %s --outdir stdout %s > %s" %\
				(self.software_names.get_genbank2gff3(), " --filter ".join(vect_filter), genbank, temp_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run snippy-vcf-to-tab")
		
		### filter file
		handle = open(temp_file)
		handle_write = open(out_file, 'w')
		for line in handle:
			sz_temp = line.strip()
			if (len(sz_temp) == 0 or sz_temp.find('# Input') == 0 or sz_temp.find('# GFF3 saved') == 0): continue
			if (sz_temp.find('##FASTA') == 0): break
			if (sz_temp[0] == '#'): handle_write.write(sz_temp + "\n")
			elif (len(sz_temp.split('\t')) > 3 and sz_temp.split('\t')[2] in vect_filter): continue
			else: handle_write.write(sz_temp + "\n")
		handle.close()
		handle_write.close()
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
			#os.unlink(snpeff_config)
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to create snpEff database")
		
		### create the annotation
		cmd = "%s ann %s -c %s -dataDir %s %s %s > %s" % (self.software_names.get_snp_eff(), self.software_names.get_snp_eff_parameters(),\
						snpeff_config, temp_dir, fasta_file_name, temp_vcf_file, out_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			#os.unlink(snpeff_config)
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
		
		temp_file = self.utils.get_temp_file("snippy_vcf_to_tab", ".gff") 
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

		### create the FAI index
		self.create_fai_fasta(reference_fasta_temp)
		
		temp_file = self.utils.get_temp_file('freebayes_temp', '.vcf')
		cmd = "%s %s -f %s -b %s > %s" %\
				(self.software_names.get_freebayes(), self.software_names.get_freebayes_parameters(),
				reference_fasta_temp, file_to_process, temp_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run freebayes")
		
		### run snpEff
		temp_file_2 = self.utils.get_temp_file("vcf_file", ".vcf")
		self.run_snpEff(reference_fasta, genbank_file, temp_file, os.path.join(temp_dir, os.path.basename(temp_file_2)))
		
		self.test_bgzip_and_tbi_in_vcf(os.path.join(temp_dir, os.path.basename(temp_file_2)))
		
		### add FREQ to vcf file
		vcf_file_out_temp = self.utils.add_freq_to_vcf(os.path.join(temp_dir, os.path.basename(temp_file_2)), os.path.join(temp_dir, sample_name + '.vcf'))
		os.unlink(temp_file)
		os.unlink(temp_file_2)
		
		### pass vcf to tab
		self.run_snippy_vcf_to_tab(reference_fasta, genbank_file, vcf_file_out_temp, "{}.tab".format(os.path.join(temp_dir, sample_name)))
		return temp_dir


	"""
	Global processing
	"""
	@transaction.atomic
	def run_fastq_and_trimmomatic(self, sample, owner):
		"""
		run fastq and trimmomatic
		Upload average and sequence numbers
		"""
		manageDatabase = ManageDatabase()
		result_all = Result()
		### first run fastq
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
			manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (temp_dir); os.system(cmd)
			return False
		cmd = "rm -r %s*" % (temp_dir); os.system(cmd)
		
		### run trimmomatic
		try:
			temp_dir = self.run_trimmomatic(sample.get_fastq(TypePath.MEDIA_ROOT, True), sample.get_fastq(TypePath.MEDIA_ROOT, False), sample.name)
			result_all.add_software(SoftwareDesc(self.software_names.get_trimmomatic_name(), self.software_names.get_trimmomatic_version(), self.software_names.get_trimmomatic_parameters()))
			### need to copy the files to samples/user path
			self.utils.copy_file(os.path.join(temp_dir, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True))), sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True))
			if (sample.exist_file_2()): self.utils.copy_file(os.path.join(temp_dir, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False))), sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False))
		except Exception as e:
			result = Result()
			result.set_error("Fail to run trimmomatic software: " + e.args[0])
			result.add_software(SoftwareDesc(self.software_names.get_trimmomatic_name(), self.software_names.get_trimmomatic_version(), self.software_names.get_trimmomatic_parameters()))
			manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (temp_dir); os.system(cmd)
			return False
		cmd = "rm -r %s*" % (temp_dir); os.system(cmd)
		
		
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
			manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (temp_dir); os.system(cmd)
			return False
		cmd = "rm -r %s*" % (temp_dir); os.system(cmd)

		### collect numbers
		(lines_1, average_1) = self.get_lines_and_average_reads(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True))
		if (sample.exist_file_2()): (lines_2, average_2) = self.get_lines_and_average_reads(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False))
		else: (lines_2, average_2) = (None, None)
		result_average = ResultAverageAndNumberReads(lines_1, average_1, lines_2, average_2)
		manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, MetaKeyAndValue.META_VALUE_Success, result_average.to_json())

		## save everything OK
		manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, MetaKeyAndValue.META_VALUE_Success, "Success, Fastq(%s), Trimmomatic(%s)" %\
							(self.software_names.get_fastq_version(), self.software_names.get_trimmomatic_version()))
		manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software, MetaKeyAndValue.META_VALUE_Success, result_all.to_json())

		### set the flag of the end of the task		
		meta_sample = manageDatabase.get_metakey(sample, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue)
		if (meta_sample != None):
			manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Success, meta_sample.description)
		return True

	"""
	Global processing, fastQ, trimmomatic and GetSpecies
	"""
	def run_fastq_and_trimmomatic_and_identify_species(self, sample, user):

		### run trimmomatics
		b_return = self.run_fastq_and_trimmomatic(sample, user)
		
		### queue the quality check and
		if (b_return and sample.exist_file_2()):	## don't run for single file because spades doesn't work for one single file
			self.identify_type_and_sub_type(sample, sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True),\
				sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False), user)

		## set the flag that is ready for process
		if (b_return):
			sample_to_update = Sample.objects.get(pk=sample.id)
			sample_to_update.is_ready_for_projects = True
			sample_to_update.type_subtype = sample_to_update.get_type_sub_type()
			sample_to_update.save()
		return b_return


	"""
	Global processing, fastQ, trimmomatic and GetSpecies
	"""
	@transaction.atomic
	def process_second_stage_snippy_coverage_freebayes(self, project_sample, user):
		"""
		Global processing, snippy, coverage, 
		"""
		
		manageDatabase = ManageDatabase()
		result_all = Result()
		### metakey for this process
		metaKeyAndValue = MetaKeyAndValue()
		meta_key_project_sample = metaKeyAndValue.get_meta_key_queue_by_project_sample_id(project_sample.id)

		### Test if this sample already run		
		meta_sample = manageDatabase.get_project_sample_metakey_last(project_sample, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue)
		if (meta_sample != None and meta_sample.value == MetaKeyAndValue.META_VALUE_Success): return 

		## process snippy
		try:
			out_put_path = self.run_snippy(project_sample.sample.get_fastq(TypePath.MEDIA_ROOT, True),\
					project_sample.sample.get_fastq(TypePath.MEDIA_ROOT, False),\
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
			return False

		## copy the files to the project sample directories
		self.copy_files_to_project(project_sample, self.software_names.get_snippy_name(), out_put_path)
		remove_path = os.path.dirname(out_put_path)
		if (len(remove_path.split('/')) > 2): self.utils.remove_dir(remove_path)
		else: self.utils.remove_dir(out_put_path)

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
					meta_key = metaKeyAndValue.get_meta_key_alert_coverage(MetaKeyAndValue.META_KEY_ALERT_COVERAGE_9, element)
					manageDatabase.set_project_sample_metakey(project_sample, user, meta_key, MetaKeyAndValue.META_VALUE_Error, coverage.get_fault_message_9(element))
				elif (not coverage.is_100_more_0(element)):
					project_sample.alert_first_level += 1
					meta_key = metaKeyAndValue.get_meta_key_alert_coverage(MetaKeyAndValue.META_KEY_ALERT_COVERAGE_0, element)
					manageDatabase.set_project_sample_metakey(project_sample, user, meta_key, MetaKeyAndValue.META_VALUE_Error, coverage.get_fault_message_0(element))
			project_sample.save()
			
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
			return False
		
		meta_sample = manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage,\
								MetaKeyAndValue.META_VALUE_Success, coverage.to_json())

		## run freebayes
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
			return False
		
		## count hits from tab file
		file_tab = os.path.join(out_put_path, project_sample.sample.name + ".tab")
		if (os.path.exists(file_tab)):
			vect_count_type = ['snp']
			count_hits = self.utils.count_hits_from_tab(file_tab, vect_count_type)
			### set flag that is finished
			manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Count_Hits, MetaKeyAndValue.META_VALUE_Success, count_hits.to_json())
		
		self.copy_files_to_project(project_sample, self.software_names.get_freebayes_name(), out_put_path)
		self.utils.remove_dir(out_put_path)
		
		### make the coverage images
		draw_all_coverage = DrawAllCoverage()
		draw_all_coverage.draw_all_coverages(project_sample)
		
		### get again
		project_sample = ProjectSample.objects.get(pk=project_sample.id)
		project_sample.is_finished = True
		project_sample.is_deleted = False
		project_sample.is_error = False
		project_sample.save()
		
		#### set the alerts of count_hits based on percentiles
		manage_database = ManageDatabase()
		tagNamesConstants = TagNamesConstants()
		manage_database.add_variation_count(project_sample, user, count_hits)
		percentil_name = tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_INSAFLU, TagNamesConstants.TAG_PERCENTIL_VAR_INSAFLU)
		manage_database.set_percentis_alert(project_sample, user, count_hits, percentil_name)
		
		### set the flag of the end of the task		
		meta_sample = manageDatabase.get_project_sample_metakey(project_sample, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue)
		if (meta_sample != None):
			manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Success, meta_sample.description)
		return True


