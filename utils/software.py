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
from utils.constants import Constants
from manage_virus.models import UploadFile
from manage_virus.uploadFiles import UploadFiles
from managing_files.manage_database import ManageDatabase
from utils.result import Result, SoftwareDesc

class Software(object):
	'''
	classdocs
	'''

	CORES_TO_USE = 3
	
	## dir with software
	DIR_SOFTWARE = "/usr/local/software/insaflu"
	SOFTWARE_SAMTOOLS = "/usr/bin/samtools"
	SOFTWARE_SAMTOOLS_VERSION = ""
	SOFTWARE_SPAdes = os.path.join(DIR_SOFTWARE, "SPAdes-3.11.1-Linux/bin/spades.py") 
	SOFTWARE_SPAdes_VERSION = "3.11.1"
	SOFTWARE_SPAdes_PARAMETERS = ""
	SOFTWARE_ABRICATE = os.path.join(DIR_SOFTWARE, "abricate/bin/abricate")
	SOFTWARE_ABRICATE_DB = os.path.join(DIR_SOFTWARE, "abricate/db")
	SOFTWARE_ABRICATE_VERSION = "0.8-dev"
	SOFTWARE_ABRICATE_PARAMETERS = ""
	SOFTWARE_FASTQ = os.path.join(DIR_SOFTWARE, "FastQC/fastqc")
	SOFTWARE_FASTQ_VERSION = "0.11.5"
	SOFTWARE_FASTQ_PARAMETERS = ""
	SOFTWARE_TIMMOMATIC = os.path.join(DIR_SOFTWARE, "trimmomatic/classes/trimmomatic.jar")
	SOFTWARE_TIMMOMATIC_VERSION = "0.27"
	SOFTWARE_TIMMOMATIC_PARAMETERS = "SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:55 TOPHRED33"

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
	def get_spades_version(self): return self.SOFTWARE_SPAdes_VERSION
	def get_spades_parameters(self): return self.SOFTWARE_SPAdes_PARAMETERS
	
	"""
	return abricate software
	"""
	def get_abricate(self): return self.SOFTWARE_ABRICATE
	def get_abricate_version(self): return self.SOFTWARE_ABRICATE_VERSION
	def get_abricate_parameters(self): return self.SOFTWARE_ABRICATE_PARAMETERS
	
	"""
	return FASTq software
	"""
	def get_fastq(self): return self.SOFTWARE_FASTQ
	def get_fastq_version(self): return self.SOFTWARE_FASTQ_VERSION
	def get_fastq_parameters(self): return self.SOFTWARE_FASTQ_PARAMETERS
	
	"""
	return trimmomatic software
	"""
	def get_trimmomantic(self): return self.SOFTWARE_TIMMOMATIC
	def get_trimmomatic_version(self): return self.SOFTWARE_TIMMOMATIC_VERSION
	def get_trimmomatic_parameters(self): return self.SOFTWARE_TIMMOMATIC_PARAMETERS
	
	def createFaiToFastaFile(self, fileFastaName):
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
		if (fastq_2 is None or len(fastq_2) == 0): cmd = "%s --pe1-1 %s --meta --only-assembler -t %d -o %s" % (self.get_spades(), fastq_1, self.CORES_TO_USE, out_dir)
		else: cmd = "%s --pe1-1 %s --pe1-2 %s --meta --only-assembler -t %d -o %s" % (self.get_spades(), fastq_1, fastq_2, self.CORES_TO_USE, out_dir)
		
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
	def identify_type_and_sub_type(self, sample, owner):
		"""
		Identify type and sub_type
		"""
		fastq1_1 = sample.path_name_1.name
		if (sample.exist_file_2()): fastq1_2 = sample.path_name_2.name
		
		manageDatabase = ManageDatabase()
		### temp dir out spades		
		out_dir_spades = self.utils.get_temp_dir()
		try:
			cmd = self.run_spades(fastq1_1, fastq1_2, out_dir_spades)
		except Exception:
			result = Result()
			result.set_error("Spades (%s) fail to run" % (self.get_spades_version()))
			result.add_software(Software(self.get_spades(), self.get_spades_version(), self.get_spades_parameters()))
			manageDatabase.set_metakey(sample, owner, Constants.META_KEY_Identify_Sample, Constants.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False
		
		file_out = os.path.join(out_dir_spades, "contigs.fasta")
		if (not os.path.exists(file_out) or os.path.getsize(file_out) < 100):
			## save error in MetaKeySample
			result = Result()
			result.set_error("Spades (%s) fail to run" % (self.get_spades_version()))
			result.add_software(SoftwareDesc(self.get_spades(), self.get_spades_version(), self.get_spades_parameters()))
			manageDatabase.set_metakey(sample, owner, Constants.META_KEY_Identify_Sample, Constants.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False
		
		try:
			uploadFile = UploadFile.objects.order_by('-version')[0]
		except UploadFile.DoesNotExist:
			## save error in MetaKeySample
			result = Result()
			result.set_error("Abricate (%s) fail to run" % (self.get_abricate_version()))
			result.add_software(SoftwareDesc(self.get_abricate(), self.get_abricate_version(), self.get_abricate_parameters()))
			manageDatabase.set_metakey(sample, owner, Constants.META_KEY_Identify_Sample, Constants.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False

		if (not self.is_exist_database_abricate(uploadFile.abricate_name)):
			try:
				self.create_database_abricate(uploadFile.abricate_name, uploadFile.path)
			except Exception:
				result = Result()
				result.set_error("Abricate (%s) fail to run --setupdb" % (self.get_abricate_version()))
				result.add_software(SoftwareDesc(self.get_abricate(), self.get_abricate_version(), self.get_abricate_parameters()))
				manageDatabase.set_metakey(sample, owner, Constants.META_KEY_Identify_Sample, Constants.META_VALUE_Error, result.to_json())
				cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
				return False
		
		## run abricate
		out_file_abricate = self.utils.get_temp_file("temp_abricate", ".txt")
		try:
			cmd = self.run_abricate(uploadFile.abricate_name, file_out, out_file_abricate)
		except Exception:
			result = Result()
			result.set_error("Abricate (%s) fail to run" % (self.get_abricate_version()))
			result.add_software(SoftwareDesc(self.get_abricate(), self.get_abricate_version(), self.get_abricate_parameters()))
			manageDatabase.set_metakey(sample, owner, Constants.META_KEY_Identify_Sample, Constants.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False
		
		if (not os.path.exists(out_file_abricate)):
			## save error in MetaKeySample
			result = Result()
			result.set_error("Abricate (%s) fail to run" % (self.get_abricate_version()))
			result.add_software(SoftwareDesc(self.get_abricate(), self.get_abricate_version(), self.get_abricate_parameters()))
			manageDatabase.set_metakey(sample, owner, Constants.META_KEY_Identify_Sample, Constants.META_VALUE_Error, result.to_json())
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False

		parseOutFiles = ParseOutFiles()
		vect_data = parseOutFiles.parse_abricate_file(out_file_abricate)
		## copy the abricate output 
		self.utils.copy_file(out_file_abricate, self.constants.get_abricate_output(sample.path_name_1))
		
		uploadFiles = UploadFiles()
		vect_data = uploadFiles.uploadIdentifyVirus(vect_data, uploadFile.abricate_name)
		if (len(vect_data) == 0):
			## save error in MetaKeySample
			result = Result()
			result.set_error("Fail to identify type and sub type")
			result.add_software(SoftwareDesc(self.get_abricate(), self.get_abricate_version(), self.get_abricate_parameters()))
			manageDatabase.set_metakey(sample, owner, Constants.META_KEY_Identify_Sample, Constants.META_VALUE_Error, result.to_json())
			cmd = "rm %s" % (out_file_abricate); os.system(cmd)
			cmd = "rm -r %s*" % (out_dir_spades); os.system(cmd)
			return False
		
		
		for identify_virus in vect_data:
			sample.identify_virus.add(identify_virus)
		sample.save()
		
		## save everything OK
		manageDatabase.set_metakey(sample, owner, Constants.META_KEY_Identify_Sample, Constants.META_VALUE_Success, "Success, Spades(%s), Abricate(%s)" % (self.get_spades_version(), self.get_abricate_version()))
		cmd = "rm %s" % (out_file_abricate); os.system(cmd)
		self.utils.remove_dir(out_dir_spades)
		return True

	def run_fastq(self, file_name_1, file_name_2):
		"""
		run fastQ, return output direcory
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
		return output direcotry
		
		#${bpipe_trimmomatic} PE -threads 3 -basein FILE1 -baseout OUT_FOLDER/PREFIX_FILES_OUT.fastq.gz SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:55 TOPHRED33
		
		PE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
   			or: 
		SE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] <inputFile> <outputFile> <trimmer1>
	   
		"""
		utils = Utils()
		temp_dir = utils.get_temp_dir()
		if (file_name_2 is None or len(file_name_2) == 0):
			b_first_file = True
			cmd = "java -jar %s SE -threads %d %s %s %s" % (self.get_trimmomantic(), Software.CORES_TO_USE, file_name_1, 
					os.path.join(temp_dir, os.path.basename(self.constants.get_trimmomatic_output(temp_dir, sample_name, b_first_file))), 
					self.get_trimmomatic_parameters())
			print(cmd)
		else:
			cmd = "java -jar %s PE -threads %d -basein %s -baseout %s.fastq.gz %s" % (self.get_trimmomantic(), Software.CORES_TO_USE, 
										file_name_1, os.path.join(temp_dir, sample_name), self.get_trimmomatic_parameters())
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run trimmomatic")
		return temp_dir
		
	"""
	Global processing
	"""
	def run_fastq_and_trimmomatic(self, sample, owner):
		"""
		run fastq and trimmomatic
		Upload several tags about the quality
		"""
		pass
