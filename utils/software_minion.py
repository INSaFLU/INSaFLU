'''
Created on 01/01/2021

@author: mmp
'''
import os, logging, humanfriendly
from constants.constants import Constants, TypePath
from constants.meta_key_and_values import MetaKeyAndValue
from django.conf import settings
from utils.process_SGE import ProcessSGE
from managing_files.models import Sample, MixedInfectionsTag, ProcessControler
from managing_files.manage_database import ManageDatabase
from constants.software_names import SoftwareNames
from utils.utils import Utils
from utils.result import Result, SoftwareDesc, ResultAverageAndNumberReads
from utils.software import Software
from settings.default_software import DefaultSoftware

######################################
####   Minion methods
class SoftwareMinion(object):
	'''
	classdocs
	'''
	utils = Utils()
	software_names = SoftwareNames()
	software = Software()
	
	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")

	def __init__(self):
		'''
		Constructor
		'''
		pass
	"""
	Global processing, RabbitQC, NanoStat, NanoFilt and GetSpecies
	"""
	def run_clean_minion(self, sample, user):

		print("Start ProcessControler")
		process_controler = ProcessControler()
		process_SGE = ProcessSGE()
		process_SGE.set_process_controler(user, process_controler.get_name_sample(sample), ProcessControler.FLAG_RUNNING)
		
		### it can be deleted
		if (sample.is_deleted or not sample.is_valid_1):
			process_SGE.set_process_controler(user, process_controler.get_name_sample(sample), ProcessControler.FLAG_FINISHED)
			return True

		try:
			print("Start run_clean_minion")
			### run stat and rabbit for Images
			b_return = self.run_nanofilt_and_stat(sample, user)
			
			print("Result run_clean_minion: " + str(b_return))
			
			### queue the quality check and
			if (b_return):	## don't run for single file because spades doesn't work for one single file
				self.software.identify_type_and_sub_type(sample, sample.get_nanofilt_file(TypePath.MEDIA_ROOT),\
					None, user)
	
			## set the flag that is ready for process
			if (b_return):
				sample_to_update = Sample.objects.get(pk=sample.id)
				sample_to_update.is_ready_for_projects = True
				sample_to_update.type_subtype = Constants.EMPTY_VALUE_TYPE_SUBTYPE
				
				manage_database = ManageDatabase()
				message = "Warning: INSaFLU don't identify species for Nanopore sequences."
				manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE,\
							MetaKeyAndValue.META_VALUE_Success, message)
				
				if (sample_to_update.number_alerts == None): sample_to_update.number_alerts = 1
				else: sample_to_update.number_alerts += 1
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
	

	def run_nanofilt_and_stat(self, sample, owner):
		"""
		run clean and stat before and after
		"""
		manage_database = ManageDatabase()
		result_all = Result()
		
		### first try run down size if necessary
		if (settings.DOWN_SIZE_FASTQ_FILES):
			(is_downsized, file_name_1, file_name_2) = self.software.make_downsize(sample.get_fastq(TypePath.MEDIA_ROOT, True),\
						sample.get_fastq(TypePath.MEDIA_ROOT, False), settings.MAX_FASTQ_FILE_UPLOAD)
			if (is_downsized):
				if (os.path.exists(file_name_1) and os.path.getsize(file_name_1) > 100):
					self.utils.move_file(file_name_1, sample.get_fastq(TypePath.MEDIA_ROOT, True))
				if (file_name_2 != None and len(file_name_2) > 0 and os.path.exists(file_name_2) and os.path.getsize(file_name_2) > 100): 
					self.utils.move_file(file_name_2, sample.get_fastq(TypePath.MEDIA_ROOT, False))
				
				### set the downsize message
				manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_ALERT_DOWNSIZE_OF_FASTQ_FILES,\
											MetaKeyAndValue.META_VALUE_Success,\
											"Fastq files were down sized to values ~{}.".format( humanfriendly.format_size(int(settings.MAX_FASTQ_FILE_UPLOAD)) ))
		
		### first run stat
		try:
			result_nano_stat = self.run_nanostat(sample.get_fastq(TypePath.MEDIA_ROOT, True))
			result_all.add_software(SoftwareDesc(self.software_names.get_NanoStat_name(), self.software_names.get_NanoStat_version(),
					self.software_names.get_NanoStat_parameters()))
			
		except Exception as e:
			result = Result()
			result.set_error("Fail to run nanostat software: " + e.args[0])
			result.add_software(SoftwareDesc(self.software_names.get_NanoStat_name(), self.software_names.get_NanoStat_version(), self.software_names.get_NanoStat_parameters()))
			manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			return False
		
		### run rabbitQC, only to have a image of the data
		try:
			out_file_html = self.run_rabbitQC(sample.get_fastq(TypePath.MEDIA_ROOT, True))
			result_all.add_software(SoftwareDesc(self.software_names.get_rabbitQC_name(), self.software_names.get_rabbitQC_version(),
					self.software_names.get_rabbitQC_parameters()))
			self.utils.copy_file(out_file_html, sample.get_rabbitQC_output(TypePath.MEDIA_ROOT))
			self.utils.remove_file(out_file_html)
		except Exception as e:
			result = Result()
			result.set_error("Fail to run rabbitQC software: " + e.args[0])
			result.add_software(SoftwareDesc(self.software_names.get_rabbitQC_name(), self.software_names.get_rabbitQC_version(),
					self.software_names.get_rabbitQC_parameters()))
			manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			return False
		
		### run nanoflit
		try:
			(result_file, parameters) = self.run_nanofilt(sample.get_fastq(TypePath.MEDIA_ROOT, True), owner)
			result_all.add_software(SoftwareDesc(self.software_names.get_NanoFilt_name(), self.software_names.get_NanoFilt_version(), parameters))
			
			### need to copy the files to samples/user path
			self.utils.copy_file(result_file, sample.get_nanofilt_file(TypePath.MEDIA_ROOT))
			self.utils.remove_file(result_file)
		except Exception as e:
			result = Result()
			result.set_error("Fail to run NanoFilt software: " + e.args[0])
			result.add_software(SoftwareDesc(self.software_names.get_NanoFilt_name(), self.software_names.get_NanoFilt_version(),
						self.software_names.get_NanoFilt_parameters()))
			manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			return False
		
		### run stat again
		try:
			result_nano_stat = self.run_nanostat(sample.get_nanofilt_file(TypePath.MEDIA_ROOT))
			result_all.add_software(SoftwareDesc(self.software_names.get_NanoStat_name(), self.software_names.get_NanoStat_version(),
					self.software_names.get_NanoStat_parameters()))
			
		except Exception as e:
			result = Result()
			result.set_error("Fail to run nanostat software: " + e.args[0])
			result.add_software(SoftwareDesc(self.software_names.get_NanoStat_name(), self.software_names.get_NanoStat_version(), self.software_names.get_NanoStat_parameters()))
			manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			return False
		
		### run rabbitQC, only to have a image of the data
		try:
			out_file_html = self.run_rabbitQC(sample.get_nanofilt_file(TypePath.MEDIA_ROOT))
			result_all.add_software(SoftwareDesc(self.software_names.get_rabbitQC_name(), self.software_names.get_rabbitQC_version(),
					self.software_names.get_rabbitQC_parameters()))
			self.utils.copy_file(out_file_html, sample.get_rabbitQC_nanofilt(TypePath.MEDIA_ROOT))
			self.utils.remove_file(out_file_html)
		except Exception as e:
			result = Result()
			result.set_error("Fail to run rabbitQC software: " + e.args[0])
			result.add_software(SoftwareDesc(self.software_names.get_rabbitQC_name(), self.software_names.get_rabbitQC_version(),
					self.software_names.get_rabbitQC_parameters()))
			manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, MetaKeyAndValue.META_VALUE_Error, result.to_json())
			return False
		
		### collect numbers
		(lines_1, average_1) = self.software.get_lines_and_average_reads(sample.get_nanofilt_file(TypePath.MEDIA_ROOT))
		result_average = ResultAverageAndNumberReads(lines_1, average_1, None, None)
		manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, MetaKeyAndValue.META_VALUE_Success, result_average.to_json())

		## save everything OK
		manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt, MetaKeyAndValue.META_VALUE_Success, "Success, NanoStat(%s), NanoFilt(%s)" %\
							(self.software_names.get_NanoStat_version(), self.software_names.get_NanoFilt_version()))
		manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software, MetaKeyAndValue.META_VALUE_Success, result_all.to_json())


		### set the flag of the end of the task		
		meta_sample = manage_database.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue)
		if (meta_sample != None):
			manage_database.set_sample_metakey(sample, owner, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Success, meta_sample.description)
		return result_average.has_reads()

	
	def run_nanostat(self, file_name):
		"""
		run nanoStat, return result with keys
		crate a text file with output
			General summary:         
			Mean read length:                  490.4
			Mean read quality:                  12.6
			Median read length:                488.0
			Median read quality:                12.7
			Number of reads:               211,388.0
			Read length N50:                   489.0
			STDEV read length:                  17.2
			Total bases:               103,672,824.0
			Number, percentage and megabases of reads above quality cutoffs
			>Q5:	211388 (100.0%) 103.7Mb
			>Q7:	211388 (100.0%) 103.7Mb
			>Q10:	191968 (90.8%) 93.9Mb
			>Q12:	136331 (64.5%) 66.5Mb
			>Q15:	16830 (8.0%) 8.2Mb
		"""
		temp_dir = self.utils.get_temp_dir()
		out_file_name = "temp.txt"
		out_path_file_name = os.path.join(temp_dir, out_file_name)
		cmd = "{} -o {} -n {} --fastq {}".format(self.software_names.get_NanoStat(), temp_dir, out_file_name, file_name)
		exist_status = os.system(cmd)
		if (exist_status != 0 or not os.path.exists(out_path_file_name)):
			self.utils.remove_dir(temp_dir)
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run run_fastq")
		
		### process out file name TOFO
		result = Result()
	
		self.utils.remove_dir(temp_dir)
		return result
	
	
	def run_rabbitQC(self, file_name):
		"""
		run rabbitQC, return output directory
		only for visual data
		return a html file
		"""
		temp_file = self.utils.get_temp_file("rabbit", ".html")
		cmd = "{} {} -i {} -h {}".format(self.software_names.get_rabbitQC(), self.software_names.get_rabbitQC_parameters(),\
				file_name, temp_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.utils.remove_file(temp_file)
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run rabbitQC")
		return temp_file
	
	
	def run_nanofilt(self, file_name, user = None):
		"""
		run nanofilt
		:out clean file and parameters
		"""

		## get dynamic parameters
		if (user is None):
			parameters = self.software_names.get_NanoFilt_parameters()
		else:
			default_software = DefaultSoftware()
			parameters = default_software.get_parameters(self.software_names.get_NanoFilt_name(), user)

		### run software
		temp_file_name = self.utils.get_temp_file("nanofilt", "fastq.fz")
		cmd = "gzip -cd {} | {} {} | gzip > {}".format(
				file_name, self.software_names.get_NanoFilt(), parameters,
				temp_file_name)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.utils.remove_file(temp_file_name)
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run NanoFilt")
		return (temp_file_name, parameters)

