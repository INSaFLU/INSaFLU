'''
Created on Oct 28, 2017

@author: mmp
'''
import os
import logging
import cmd
import subprocess
from utils.Utils import Utils
from utils.ParseOutFiles import ParseOutFiles
from utils.Constants import Constants
from manage_virus.models import UploadFile
from manage_virus.uploadFiles import UploadFiles
from managing_files.manage_database import ManageDatabase

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
	SOFTWARE_ABRICATE = os.path.join(DIR_SOFTWARE, "abricate/bin/abricate")
	SOFTWARE_ABRICATE_DB = os.path.join(DIR_SOFTWARE, "abricate/db")
	SOFTWARE_ABRICATE_VERSION = "0.8-dev"
	
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	utils = Utils()
	
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
	
	"""
	return spades software
	"""
	def get_abricate(self): return self.SOFTWARE_ABRICATE
	def get_abricate_version(self): return self.SOFTWARE_ABRICATE_VERSION
	
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
			raise Exception("Fail to run samtools")
		
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


	def identify_type_and_sub_type(self, sample, owner):
		"""
		Identify the type and sub_type
		"""
		fastq1_1 = sample.files.path_name_1.name
		fastq1_2 = None if sample.files.path_name_2 == None else sample.files.path_name_2.name
		
		manageDatabase = ManageDatabase()
		### temp dir out spades		
		out_dir = self.utils.get_temp_dir()
		try:
			cmd = self.software.run_spades(fastq1_1, fastq1_2, out_dir)
		except Exception:
			manageDatabase.set_metakey(sample, owner, Constants.META_KEY_Identify_Sample, Constants.META_VALUE_Error, "Spades (%s) fail to run" % (self.get_spades_version()))
			cmd = "rm -r %s*" % (out_dir)
			os.system(cmd)
			return
		
		file_out = os.path.join(out_dir, "contigs.fasta")

		if (not os.path.exists(file_out) or os.path.getsize(file_out) > 100):
			## save error in MetaKeySample
			manageDatabase.set_metakey(sample, owner, Constants.META_KEY_Identify_Sample, Constants.META_VALUE_Error, "Spades (%s) fail to run" % (self.get_spades_version()))
			cmd = "rm -r %s*" % (out_dir)
			os.system(cmd)
			return
		try:
			uploadFile = UploadFile.objects.order_by('-version')[0]
		except UploadFile.DoesNotExist:
			## save error in MetaKeySample
			manageDatabase.set_metakey(sample, owner, Constants.META_KEY_Identify_Sample, Constants.META_VALUE_Error, "There's no 'abricate' files to create the database")
			cmd = "rm -r %s*" % (out_dir)
			os.system(cmd)
			return

		if (not self.software.is_exist_database_abricate(uploadFile.abricate_name)):
			self.software.create_database_abricate(uploadFile.abricate_name, uploadFile.path)
		
		## run abricate
		out_file = self.utils.get_temp_file("temp_abricate", ".txt")
		try:
			cmd = self.software.run_abricate(uploadFile.abricate_name, file_out, out_file)
		except Exception:
			manageDatabase.set_metakey(sample, owner, Constants.META_KEY_Identify_Sample, Constants.META_VALUE_Error, "Abricate (%s) fail to run" % (self.get_abricate_version()))
			cmd = "rm -r %s*" % (out_dir)
			os.system(cmd)
			return
		
		if (not os.path.exists(out_file)):
			## save error in MetaKeySample
			manageDatabase.set_metakey(sample, owner, Constants.META_KEY_Identify_Sample, Constants.META_VALUE_Error, "Abricate (%s) fail to run" % (self.get_abricate_version()))
			cmd = "rm -r %s*" % (out_dir); os.system(cmd)
			return

		parseOutFiles = ParseOutFiles()
		vect_data = parseOutFiles.parse_abricate_file(out_file)
		
		uploadFiles = UploadFiles()
		vect_data = uploadFiles.uploadIdentifyVirus(vect_data, uploadFile.abricate_name)
		if (len(vect_data) == 0):
			## save error in MetaKeySample
			manageDatabase.set_metakey(sample, owner, Constants.META_KEY_Identify_Sample, Constants.META_VALUE_Error, "Fail to upload the Abricate (%s) to database" % (self.get_abricate_version()))
			cmd = "rm %s" % (out_file); os.system(cmd)
			cmd = "rm -r %s*" % (out_dir); os.system(cmd)
			return
		
		for identify_virus in vect_data:
			sample.identify_virus.add(identify_virus)
		sample.save()
		
		## save everything OK
		manageDatabase.set_metakey(sample, owner, Constants.META_KEY_Identify_Sample, Constants.META_VALUE_Success, "Success, Spades(%s), Abricate(%s)" % (self.get_spades_version(), self.get_abricate_version()))
		cmd = "rm %s" % (out_file); os.system(cmd)
		cmd = "rm -r %s*" % (out_dir); os.system(cmd)

	