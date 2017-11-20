'''
Created on Oct 28, 2017

@author: mmp
'''
import os
import logging
import cmd
import subprocess
from utils.utils import Utils
from utils.parseOutFiles import ParseOutFiles
from utils.constants import Constants, TypePath, FileType
from utils.meta_key_and_values import MetaKeyAndValue
from manage_virus.models import UploadFile
from managing_files.models import Sample, ProjectSample
from manage_virus.uploadFiles import UploadFiles
from managing_files.manage_database import ManageDatabase
from utils.result import Result, SoftwareDesc, ResultAverageAndNumberReads
from utils.parse_coverage_file import GetCoverage
from django.db import transaction

class Software(object):
	'''
	classdocs
	'''

	CORES_TO_USE = 3
	
	## dir with software
	DIR_SOFTWARE = "/usr/local/software/insaflu"
	SOFTWARE_SAMTOOLS = os.path.join(DIR_SOFTWARE, "snippy/bin/samtools")
	SOFTWARE_SAMTOOLS_name = "Samtools"
	SOFTWARE_SAMTOOLS_VERSION = "1.3"
	SOFTWARE_SAMTOOLS_PARAMETERS = ""
	SOFTWARE_BGZIP = os.path.join(DIR_SOFTWARE, "snippy/bin/bgzip")
	SOFTWARE_BGZIP_name = "bgzip"
	SOFTWARE_BGZIP_VERSION = "1.3"
	SOFTWARE_BGZIP_PARAMETERS = ""
	SOFTWARE_TABIX = os.path.join(DIR_SOFTWARE, "snippy/bin/tabix")
	SOFTWARE_TABIX_name = "tabix"
	SOFTWARE_TABIX_VERSION = "1.3"
	SOFTWARE_TABIX_PARAMETERS = ""

	SOFTWARE_SPAdes = os.path.join(DIR_SOFTWARE, "SPAdes-3.11.1-Linux/bin/spades.py") 
	SOFTWARE_SPAdes_name = "SPAdes" 
	SOFTWARE_SPAdes_VERSION = "3.11.1"
	SOFTWARE_SPAdes_PARAMETERS = ""
	SOFTWARE_ABRICATE = os.path.join(DIR_SOFTWARE, "abricate/bin/abricate")
	SOFTWARE_ABRICATE_name = "Abricate"
	SOFTWARE_ABRICATE_DB = os.path.join(DIR_SOFTWARE, "abricate/db")
	SOFTWARE_ABRICATE_VERSION = "0.8-dev"
	SOFTWARE_ABRICATE_PARAMETERS = ""
	SOFTWARE_FASTQ = os.path.join(DIR_SOFTWARE, "FastQC/fastqc")
	SOFTWARE_FASTQ_name = "FastQC"
	SOFTWARE_FASTQ_VERSION = "0.11.5"
	SOFTWARE_FASTQ_PARAMETERS = ""
	SOFTWARE_TRIMMOMATIC = os.path.join(DIR_SOFTWARE, "trimmomatic/classes/trimmomatic.jar")
	SOFTWARE_TRIMMOMATIC_name = "Trimmomatic"
	SOFTWARE_TRIMMOMATIC_VERSION = "0.27"
	SOFTWARE_TRIMMOMATIC_PARAMETERS = "SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:55 TOPHRED33"
	SOFTWARE_SNIPPY = os.path.join(DIR_SOFTWARE, "snippy/bin/snippy")
	SOFTWARE_SNIPPY_name = "Snippy"
	SOFTWARE_SNIPPY_VERSION = "3.2-dev"
	SOFTWARE_SNIPPY_PARAMETERS = "--mapqual 20 --mincov 10 --minfrac 0.51"
	SOFTWARE_FREEBAYES = os.path.join(DIR_SOFTWARE, "snippy/bin/freebayes")
	SOFTWARE_FREEBAYES_name = "Freebayes"
	SOFTWARE_FREEBAYES_VERSION = "v1.1.0-54-g49413aa"
	SOFTWARE_FREEBAYES_PARAMETERS = "-p 1 -q 20 -m 20 --min-coverage 100 --min-alternate-fraction 0.01 --min-alternate-count 10 -V"
	SOFTWARE_COVERAGE = "Coverage, inhouse script"
	SOFTWARE_COVERAGE_name = "Coverage"
	SOFTWARE_COVERAGE_VERSION = "v1.1"
	SOFTWARE_COVERAGE_PARAMETERS = ""


	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	utils = Utils()
	constants = Constants()

	def __init__(self):
		'''
		Constructor
		'''
		pass

	"""
	return samtools software
	"""
	def get_samtools(self): return self.SOFTWARE_SAMTOOLS
	def get_samtools_version(self): return self.SOFTWARE_SAMTOOLS_VERSION

	"""
	return spades software
	"""
	def get_spades(self): return self.SOFTWARE_SPAdes
	def get_spades_name(self): return self.SOFTWARE_SPAdes_name
	def get_spades_version(self): return self.SOFTWARE_SPAdes_VERSION
	def get_spades_parameters(self): return self.SOFTWARE_SPAdes_PARAMETERS

	"""
	return abricate software
	"""
	def get_abricate(self): return self.SOFTWARE_ABRICATE
	def get_abricate_name(self): return self.SOFTWARE_ABRICATE_name
	def get_abricate_version(self): return self.SOFTWARE_ABRICATE_VERSION
	def get_abricate_parameters(self): return self.SOFTWARE_ABRICATE_PARAMETERS

	"""
	return FASTq software
	"""
	def get_fastq(self): return self.SOFTWARE_FASTQ
	def get_fastq_name(self): return self.SOFTWARE_FASTQ_name
	def get_fastq_version(self): return self.SOFTWARE_FASTQ_VERSION
	def get_fastq_parameters(self): return self.SOFTWARE_FASTQ_PARAMETERS
	
	"""
	return trimmomatic software
	"""
	def get_trimmomatic(self): return self.SOFTWARE_TRIMMOMATIC
	def get_trimmomatic_name(self): return self.SOFTWARE_TRIMMOMATIC_name
	def get_trimmomatic_version(self): return self.SOFTWARE_TRIMMOMATIC_VERSION
	def get_trimmomatic_parameters(self): return self.SOFTWARE_TRIMMOMATIC_PARAMETERS
	
	"""
	return snippy software
	"""
	def get_snippy(self): return self.SOFTWARE_SNIPPY
	def get_snippy_name(self): return self.SOFTWARE_SNIPPY_name
	def get_snippy_version(self): return self.SOFTWARE_SNIPPY_VERSION
	def get_snippy_parameters(self): return self.SOFTWARE_SNIPPY_PARAMETERS

	"""
	return freebayes software
	"""
	def get_freebayes(self): return self.SOFTWARE_FREEBAYES
	def get_freebayes_name(self): return self.SOFTWARE_FREEBAYES_name
	def get_freebayes_version(self): return self.SOFTWARE_FREEBAYES_VERSION
	def get_freebayes_parameters(self): return self.SOFTWARE_FREEBAYES_PARAMETERS

	"""
	return bgzip software
	"""
	def get_bgzip(self): return self.SOFTWARE_BGZIP
	def get_bgzip_name(self): return self.SOFTWARE_BGZIP_name
	def get_bgzip_version(self): return self.SOFTWARE_BGZIP_VERSION
	def get_bgzip_parameters(self): return self.SOFTWARE_BGZIP_PARAMETERS

	"""
	return tabix software
	"""
	def get_tabix(self): return self.SOFTWARE_TABIX
	def get_tabix_name(self): return self.SOFTWARE_TABIX_name
	def get_tabix_version(self): return self.SOFTWARE_TABIX_VERSION
	def get_tabix_parameters(self): return self.SOFTWARE_TABIX_PARAMETERS
	
	"""
	return Coverage software
	"""
	def get_coverage(self): return self.SOFTWARE_COVERAGE
	def get_coverage_name(self): return self.SOFTWARE_COVERAGE_name
	def get_coverage_version(self): return self.SOFTWARE_COVERAGEVERSION
	def get_coverage_parameters(self): return self.SOFTWARE_COVERAGE_PARAMETERS


	def get_vect_type_files_to_copy(self, software):
		"""
		get type of files to copy
		"""
		if (software == Software.SOFTWARE_SNIPPY_name):
			return [FileType.FILE_BAM, FileType.FILE_BAM_BAI, FileType.FILE_CONSENSUS_FA, FileType.FILE_DEPTH_GZ, FileType.FILE_DEPTH_GZ_TBI,\
				FileType.FILE_TAB, FileType.FILE_VCF_GZ, FileType.FILE_VCF_GZ_TBI, FileType.FILE_CSV]
		elif (software == Software.SOFTWARE_FREEBAYES_name):
			return [FileType.FILE_VCF]


	def copy_files_to_project(self, project_sample, software, path_from):
		"""
		copy files to the project
		software : SOFTWARE_SNIPPY_name, SOFTWARE_FREEBAYES_name
		"""
		utils = Utils()
		for type_file in self.get_vect_type_files_to_copy(software):
			if (type_file == FileType.FILE_CONSENSUS_FA):	## if .fa file pass to .fasta
				utils.copy_file(os.path.join(path_from, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, type_file, software))),\
					project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, software))
			elif (type_file == FileType.FILE_VCF):	## if .fa file pass to .fasta
				### create the gzip file
				utils.compress_files(self.get_bgzip(), os.path.join(path_from, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, type_file, software))))
				### create the tabix
				utils.create_index_files(self.get_tabix(), os.path.join(path_from, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, software))))
				
				### copy both
				utils.copy_file(os.path.join(path_from, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, software))),\
					project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, software))
				utils.copy_file(os.path.join(path_from, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ_TBI, software))),\
					project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ_TBI, software))
			else:
				utils.copy_file(os.path.join(path_from, os.path.basename(project_sample.get_file_output(TypePath.MEDIA_ROOT, type_file, software))),\
					project_sample.get_file_output(TypePath.MEDIA_ROOT, type_file, software))

	def create_fai_fasta(self, fileFastaName):
		"""
		Create fai for a fasta file
		"""
		cmd = "%s faidx %s" % (self.get_samtools(), fileFastaName)
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
		if (fastq_2 is None or len(fastq_2) == 0): cmd = "%s -s %s --meta --only-assembler -t %d -o %s" % (self.get_spades(), fastq_1, self.CORES_TO_USE, out_dir)
		else: cmd = "%s --pe1-1 %s --pe1-2 %s --meta --only-assembler -t %d -o %s" % (self.get_spades(), fastq_1, fastq_2, self.CORES_TO_USE, out_dir)
		print(cmd)
		
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
		cmd = "%s --list" % (self.get_abricate())
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
		cmd = "mkdir -p %s/%s" % (self.SOFTWARE_ABRICATE_DB, database)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run make directory")
		
		## copy the file
		self.utils.copy_file(file_name, os.path.join(self.SOFTWARE_ABRICATE_DB, database, "sequences"))
		
		cmd = "%s --setupdb" % (self.get_abricate())
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run abricate --setupdb")
		

	def run_abricate(self, database, file_name, out_file):
		"""
		Run abricator
		"""
		cmd = "%s --db %s --quiet %s > %s" % (self.get_abricate(), database, file_name, out_file)
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
		utils = Utils()
		temp_file =  utils.get_temp_file("lines_and_average_", ".txt")
		cmd = "gzip -cd " + file_name + " | awk '{ s++; if ((s % 4) == 0) { count ++; size += length($0); }  } END { print \"sequences: \", count,  \"average: \", size/count }' > " + temp_file
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run get_lines_and_average_reads")
		
		###
		vect_out = utils.read_text_file(temp_file)
		if (len(vect_out) == 0):
			self.logger_production.error('can not read any data: ' + temp_file)
			self.logger_debug.error('can not read any data: ' + temp_file)
			raise Exception("Can't read any data")
		vect_data = vect_out[0].split()
		if (len(vect_data) != 4):
			self.logger_production.error('can not parse this data: ' + vect_out[0])
			self.logger_debug.error('can not parse this data: ' + vect_out[0])
			raise Exception("Can't read any data")
		if (utils.is_float(vect_data[3])): average_value = "%.1f" % (float(vect_data[3]))
		else: average_value = vect_data[3]
		
		utils.remove_temp_file(temp_file)
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
			result_all.add_software(SoftwareDesc(self.get_spades_name(), self.get_spades_version(), self.get_spades_parameters()))
		except Exception:
			result = Result()
			result.set_error("Spades (%s) fail to run" % (self.get_spades_version()))
			result.add_software(SoftwareDesc(self.get_spades_name(), self.get_spades_version(), self.get_spades_parameters()))
			manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False
		
		file_out = os.path.join(out_dir_spades, "contigs.fasta")
		if (not os.path.exists(file_out) or os.path.getsize(file_out) < 100):
			## save error in MetaKeySample
			result = Result()
			result.set_error("Spades (%s) fail to run" % (self.get_spades_version()))
			result.add_software(SoftwareDesc(self.get_spades_name(), self.get_spades_version(), self.get_spades_parameters()))
			manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False
		
		try:
			uploadFile = UploadFile.objects.order_by('-version')[0]
		except UploadFile.DoesNotExist:
			## save error in MetaKeySample
			result = Result()
			result.set_error("Abricate (%s) fail to run" % (self.get_abricate_version()))
			result.add_software(SoftwareDesc(self.get_abricate_name(), self.get_abricate_version(), self.get_abricate_parameters()))
			manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False

		if (not self.is_exist_database_abricate(uploadFile.abricate_name)):
			try:
				self.create_database_abricate(uploadFile.abricate_name, uploadFile.path)
			except Exception:
				result = Result()
				result.set_error("Abricate (%s) fail to run --setupdb" % (self.get_abricate_version()))
				result.add_software(SoftwareDesc(self.get_abricate_name(), self.get_abricate_version(), self.get_abricate_parameters()))
				manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
				cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
				return False
		
		## run abricate
		out_file_abricate = self.utils.get_temp_file("temp_abricate", ".txt")
		try:
			cmd = self.run_abricate(uploadFile.abricate_name, file_out, out_file_abricate)
			result_all.add_software(SoftwareDesc(self.get_abricate_name(), self.get_abricate_version(), self.get_abricate_parameters()))
		except Exception:
			result = Result()
			result.set_error("Abricate (%s) fail to run" % (self.get_abricate_version()))
			result.add_software(SoftwareDesc(self.get_abricate_name(), self.get_abricate_version(), self.get_abricate_parameters()))
			manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False
		
		if (not os.path.exists(out_file_abricate)):
			## save error in MetaKeySample
			result = Result()
			result.set_error("Abricate (%s) fail to run" % (self.get_abricate_version()))
			result.add_software(SoftwareDesc(self.get_abricate(), self.get_abricate_version(), self.get_abricate_parameters()))
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
			result.add_software(SoftwareDesc(self.get_abricate(), self.get_abricate_version(), self.get_abricate_parameters()))
			manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm %s" % (out_file_abricate); os.system(cmd)
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False
		
		
		for identify_virus in vect_data:
			sample.identify_virus.add(identify_virus)
		sample.save()
		
		## save everything OK
		manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Success, "Success, Spades(%s), Abricate(%s)" % (self.get_spades_version(), self.get_abricate_version()))
		manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Identify_Sample_Software, MetaKeyAndValue.META_VALUE_Success, result_all.to_json())
		cmd = "rm %s" % (out_file_abricate); os.system(cmd)
		self.utils.remove_dir(out_dir_spades)
		return True


	def run_fastq(self, file_name_1, file_name_2):
		"""
		run fastQ, return output directory
		-o OUT_FOLDER --nogroup --format fastq --threads 10 --dir OUT_FOLDER FILE1 FILE2
		"""
		utils = Utils()
		temp_dir = utils.get_temp_dir()
		if (not file_name_2 is None and len(file_name_2) > 0):
			cmd = "%s -o %s --nogroup --format fastq --threads %d --dir %s %s %s" % (self.get_fastq(), temp_dir, Software.CORES_TO_USE, 
										temp_dir, file_name_1, file_name_2)
		else: 
			cmd = "%s -o %s --nogroup --format fastq --threads %d --dir %s %s" % (self.get_fastq(), temp_dir, Software.CORES_TO_USE, temp_dir, file_name_1)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run run_fastq")
		return temp_dir


	def run_trimmomatic(self, file_name_1, file_name_2, sample_name):
		"""
		run trimmomatic
		return output directory
		
		#${bpipe_trimmomatic} PE -threads 3 -basein FILE1 -baseout OUT_FOLDER/PREFIX_FILES_OUT.fastq.gz SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:55 TOPHRED33
		
		PE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
   			or: 
		SE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] <inputFile> <outputFile> <trimmer1>
		"""

		utils = Utils()
		temp_dir = utils.get_temp_dir()
		if (file_name_2 is None or len(file_name_2) == 0):
			cmd = "java -jar %s SE -threads %d %s %s_1P.fastq.gz %s" % (self.get_trimmomatic(), Software.CORES_TO_USE, file_name_1, 
					os.path.join(temp_dir, sample_name), self.get_trimmomatic_parameters())
		else:
			cmd = "java -jar %s PE -threads %d -basein %s -baseout %s.fastq.gz %s" % (self.get_trimmomatic(), Software.CORES_TO_USE, 
										file_name_1, os.path.join(temp_dir, sample_name), self.get_trimmomatic_parameters())
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
		utils = Utils()
		temp_dir = os.path.join(utils.get_temp_dir(), sample_name)
		if (file_name_2 is None or len(file_name_2) == 0):
			cmd = "%s --cpus %d --outdir %s --prefix %s --ref %s %s --se %s" %\
				(self.get_snippy(), Software.CORES_TO_USE, temp_dir, sample_name,
				path_reference, self.get_snippy_parameters(), file_name_1)
		else:
			cmd = "%s --cpus %d --outdir %s --prefix %s --ref %s %s --R1 %s --R2 %s" %\
				(self.get_snippy(), Software.CORES_TO_USE, temp_dir, sample_name,
				path_reference, self.get_snippy_parameters(), file_name_1, file_name_2)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run snippy")
		return temp_dir


	def run_freebayes(self, bam_file, reference_fasta, sample_name):
		"""
		run freebayes
		return output directory
		
		## freebayes -p 1 -q 20 -m 20 --min-coverage 100 --min-alternate-fraction 0.01 --min-alternate-count 10 -V -f ../$input2 -b {} > {.}.vcf'
		"""
		utils = Utils()
		temp_dir = os.path.join(utils.get_temp_dir())
		
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
		
		cmd = "%s %s -f %s -b %s > %s.vcf" %\
				(self.get_freebayes(), self.get_freebayes_parameters(),
				reference_fasta_temp, file_to_process, os.path.join(temp_dir, sample_name))
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run freebayes")
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
			result_all.add_software(SoftwareDesc(self.get_fastq_name(), self.get_fastq_version(), self.get_fastq_parameters()))
			
			### need to copy the files to samples/user path
			self.utils.copy_file(os.path.join(temp_dir, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, True))), sample.get_fastq_output(TypePath.MEDIA_ROOT, True))
			if (sample.exist_file_2()): self.utils.copy_file(os.path.join(temp_dir, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, False))), sample.get_fastq_output(TypePath.MEDIA_ROOT, False))
		except Exception as e:
			result = Result()
			result.set_error("Fail to run fastq software: " + e.args[0])
			result.add_software(SoftwareDesc(self.get_fastq_name(), self.get_fastq(), self.get_fastq_parameters()))
			manageDatabase.set_metakey(sample, owner, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (temp_dir); os.system(cmd)
			return False
		cmd = "rm -r %s*" % (temp_dir); os.system(cmd)
		
		### run trimmomatic
		try:
			temp_dir = self.run_trimmomatic(sample.get_fastq(TypePath.MEDIA_ROOT, True), sample.get_fastq(TypePath.MEDIA_ROOT, False), sample.name)
			result_all.add_software(SoftwareDesc(self.get_trimmomatic_name(), self.get_trimmomatic_version(), self.get_trimmomatic_parameters()))
			### need to copy the files to samples/user path
			self.utils.copy_file(os.path.join(temp_dir, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True))), sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True))
			if (sample.exist_file_2()): self.utils.copy_file(os.path.join(temp_dir, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False))), sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False))
		except Exception as e:
			result = Result()
			result.set_error("Fail to run trimmomatic software: " + e.args[0])
			result.add_software(SoftwareDesc(self.get_trimmomatic_name(), self.get_trimmomatic_version(), self.get_trimmomatic_parameters()))
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
			result.add_software(SoftwareDesc(self.get_fastq_name(), self.get_fastq(), self.get_fastq_parameters()))
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
							(self.get_fastq_version(), self.get_trimmomatic_version()))
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
		## process snippy
		try:
			out_put_path = self.run_snippy(project_sample.sample.get_fastq(TypePath.MEDIA_ROOT, True),\
					project_sample.sample.get_fastq(TypePath.MEDIA_ROOT, False), project_sample.project.reference.reference_genbank.name,\
					project_sample.sample.name)
			result_all.add_software(SoftwareDesc(self.get_snippy_name(), self.get_snippy_version(), self.get_snippy_parameters()))
		except Exception as e:
			result = Result()
			result.set_error("Fail to run fastq software: " + e.args[0])
			result.add_software(SoftwareDesc(self.get_snippy_name(), self.get_snippy_version(), self.get_snippy_parameters()))
			manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Snippy, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_put_path); os.system(cmd)
			
			### get again and set error
			project_sample = ProjectSample.objects.get(pk=project_sample.id)
			project_sample.is_error = True
			project_sample.save()
			return False

		## copy the files to the project sample directories
		self.copy_files_to_project(project_sample, Software.SOFTWARE_SNIPPY_name, out_put_path)
		self.utils.remove_dir(out_put_path)

		## get coverage from deep file
		get_coverage = GetCoverage()
		try:
			coverage = get_coverage.get_coverage(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ,\
						Software.SOFTWARE_SNIPPY_name), project_sample.project.reference.reference_fasta.name)
		except Exception as e:
			result = Result()
			result.set_error("Fail to get coverage: " + e.args[0])
			result.add_software(SoftwareDesc(self.get_coverage_name(), self.get_coverage_version(), self.get_coverage_parameters()))
			manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_put_path); os.system(cmd)
			
			### get again and set error
			project_sample = ProjectSample.objects.get(pk=project_sample.id)
			project_sample.is_error = True
			project_sample.save()
			return False
		
		meta_sample = manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage,\
								MetaKeyAndValue.META_VALUE_Success, coverage.to_json())

		## run freebayes
		try:
			out_put_path = self.run_freebayes(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_BAM, Software.SOFTWARE_SNIPPY_name),\
						project_sample.project.reference.reference_fasta.name, project_sample.sample.name)
			result_all.add_software(SoftwareDesc(self.get_freebayes_name(), self.get_freebayes_version(), self.get_freebayes_parameters()))
		except Exception as e:
			result = Result()
			result.set_error("Fail to run freebayes software: " + e.args[0])
			result.add_software(SoftwareDesc(self.get_freebayes_name(), self.get_freebayes_version(), self.get_freebayes_parameters()))
			manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Freebayes, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_put_path); os.system(cmd)
			
			### get again and set error
			project_sample = ProjectSample.objects.get(pk=project_sample.id)
			project_sample.is_error = True
			project_sample.save()
			return False
		
		self.copy_files_to_project(project_sample, Software.SOFTWARE_FREEBAYES_name, out_put_path)
		self.utils.remove_dir(out_put_path)
		
		### set flag that is finished
		manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Snippy_Freebayes, MetaKeyAndValue.META_VALUE_Success, result_all.to_json())
		
		### get again
		project_sample = ProjectSample.objects.get(pk=project_sample.id)
		project_sample.is_finished = True
		project_sample.save()
		
		### set the flag of the end of the task		
		meta_sample = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue)
		if (meta_sample != None):
			manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Success, meta_sample.description)
		return True


