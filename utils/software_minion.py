'''
Created on 01/01/2021

@author: mmp
'''
import os, logging, humanfriendly, datetime
from constants.constants import Constants, TypePath, FileType
from constants.meta_key_and_values import MetaKeyAndValue
from settings.default_software_project_sample import DefaultProjectSoftware
from utils.coverage import DrawAllCoverage
from django.conf import settings
from utils.process_SGE import ProcessSGE
from utils.parse_out_files import ParseOutFiles
from managing_files.models import Sample, ProcessControler, ProjectSample
from managing_files.manage_database import ManageDatabase
from constants.software_names import SoftwareNames
from utils.utils import Utils
from utils.result import Result, SoftwareDesc, ResultAverageAndNumberReads, CountHits
from utils.software import Software
from settings.default_software import DefaultSoftware
from utils.parse_coverage_file import GetCoverage
from utils.mixed_infections_management import MixedInfectionsManagement
from utils.result import KeyValue

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
					self.software_names.get_NanoStat_parameters(), result_nano_stat.key_values))
			
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
					self.software_names.get_NanoStat_parameters(), result_nano_stat.key_values))
			
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
		out_path_file_name = os.path.join(temp_dir, "temp.txt")
		cmd = "{} -o {} -n {} --fastq {}".format(self.software_names.get_NanoStat(), temp_dir, out_path_file_name, file_name)
		exist_status = os.system(cmd)
		if (exist_status != 0 or not os.path.exists(out_path_file_name)):
			self.utils.remove_dir(temp_dir)
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run run_fastq")
		
		### process out STATs 
		result = Result()
		with open(out_path_file_name) as handle_in:
			for line in handle_in:
				sz_temp = line.strip()
				if (len(sz_temp) == 0 or sz_temp[0] == '#'): continue
				
				### add tags
				lst_data = sz_temp.split(':')
				if len(lst_data) == 2 and lst_data[0].strip() in SoftwareNames.SOFTWARE_NANOSTAT_vect_info_to_collect:
					result.add_key_value(KeyValue(lst_data[0].strip(), lst_data[1].strip()))
		
		### remove dir
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
		
	"""
	Global processing, Medaka, Coverage, and MixedInfections
	"""
	def process_second_stage_medaka(self, project_sample, user):
		"""
		Global processing, medaka and coverage
		"""
		### make it running 
		process_controler = ProcessControler()
		process_SGE = ProcessSGE()
		manageDatabase = ManageDatabase()
		result_all = Result()

		process_SGE.set_process_controler(user, process_controler.get_name_project_sample(project_sample), ProcessControler.FLAG_RUNNING)
		
		### metakey for this process
		metaKeyAndValue = MetaKeyAndValue()
		try:	
			meta_key_project_sample = metaKeyAndValue.get_meta_key_queue_by_project_sample_id(project_sample.id)
	
			### Test if this sample already run		
			meta_sample = manageDatabase.get_project_sample_metakey_last(project_sample, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue)
			if (meta_sample != None and meta_sample.value == MetaKeyAndValue.META_VALUE_Success): return 
	
			## process medaka
			try:
				#### need to check the models
				parameters = "-m r941_min_high_g360"
				### get snippy parameters
				out_put_path = self.run_medaka(project_sample.sample.get_nanofilt_file(TypePath.MEDIA_ROOT),\
						project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT),\
						project_sample.project.reference.get_reference_gbk(TypePath.MEDIA_ROOT),\
						project_sample.sample.name, parameters)
				result_all.add_software(SoftwareDesc(self.software_names.get_medaka_name(),\
								self.software_names.get_medaka_version(), "consensus " + parameters))
				result_all.add_software(SoftwareDesc(self.software_names.get_medaka_name(),\
								self.software_names.get_medaka_version(), "depth -aa -q 0"))
			except Exception as e:
				result = Result()
				result.set_error(e.args[0])
				result.add_software(SoftwareDesc(self.software_names.get_medaka_name(), self.software_names.get_medaka_version(),\
						self.software_names.get_medaka_parameters()))
				manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Medaka, MetaKeyAndValue.META_VALUE_Error, result.to_json())
				
				### get again and set error
				project_sample = ProjectSample.objects.get(pk=project_sample.id)
				project_sample.is_error = True
				project_sample.save()
				
				meta_sample = manageDatabase.get_project_sample_metakey_last(project_sample, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue)
				if (not meta_sample is None):
					manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Error, meta_sample.description)
				process_SGE.set_process_controler(user, process_controler.get_name_project_sample(project_sample), ProcessControler.FLAG_ERROR)
				return False
	
			## copy the files to the project sample directories
			self.software.copy_files_to_project(project_sample, self.software_names.get_medaka_name(), out_put_path)
			self.utils.remove_dir(out_put_path)
	
			### make the link for the new tab file name
			path_medaka_tab = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, self.software_names.get_medaka_name())
			if (os.path.exists(path_medaka_tab)):
				sz_file_to = project_sample.get_file_output_human(TypePath.MEDIA_ROOT, FileType.FILE_TAB, self.software_names.get_medaka_name())
				self.utils.link_file(path_medaka_tab, sz_file_to)
			
			## get coverage from deep file
			default_software = DefaultProjectSoftware()
			get_coverage = GetCoverage()
			try:
				
				### limit of the coverage for a project, can be None, if not exist
				### TODO
				b_coverage_default = True
				coverage_for_project = 10
				default_coverage_value = 10	### not used
				coverage = get_coverage.get_coverage(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ,\
 							self.software_names.get_medaka_name()),\
 							project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT),
 							None, coverage_for_project)

# 				coverage_for_project = default_software.get_snippy_single_parameter_for_project(project_sample.project, DefaultProjectSoftware.SNIPPY_COVERAGE_NAME)
# 				if (not coverage_for_project is None): coverage_for_project = int(coverage_for_project)
# 				
# 				b_coverage_default = True
# 				if (default_software.is_snippy_single_parameter_default(project_sample, DefaultProjectSoftware.SNIPPY_COVERAGE_NAME)):
# 					coverage = get_coverage.get_coverage(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ,\
# 							self.software_names.get_snippy_name()),\
# 							project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT),
# 							None, coverage_for_project)
# 				else:
# 					b_coverage_default = False
# 					default_coverage_value = default_software.get_snippy_single_parameter(project_sample, DefaultProjectSoftware.SNIPPY_COVERAGE_NAME)
# 					coverage = get_coverage.get_coverage(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ,\
# 							self.software_names.get_snippy_name()),\
# 							project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT),\
# 							int(default_coverage_value), coverage_for_project)
				################################
				##################################
				### set the alerts in the coverage
				### remove possible previous alerts from others run
				for keys_to_remove in MetaKeyAndValue.VECT_TO_REMOVE_RUN_PROJECT_SAMPLE:
					manageDatabase.remove_project_sample_start_metakey(project_sample, keys_to_remove)
				
				project_sample = ProjectSample.objects.get(pk=project_sample.id)
				project_sample.alert_second_level = 0
				project_sample.alert_first_level = 0
				for element in coverage.get_dict_data():
					if (not coverage.is_100_more_9(element) and b_coverage_default):
						project_sample.alert_second_level += 1
						meta_key = metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_ALERT_COVERAGE_9, element)
						manageDatabase.set_project_sample_metakey(project_sample, user, meta_key, MetaKeyAndValue.META_VALUE_Success, coverage.get_fault_message_9(element))
					elif (not coverage.is_100_more_defined_by_user(element) and not b_coverage_default):
						project_sample.alert_second_level += 1
						meta_key = metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_ALERT_COVERAGE_value_defined_by_user, element)
						manageDatabase.set_project_sample_metakey(project_sample, user, meta_key, MetaKeyAndValue.META_VALUE_Success,
										coverage.get_fault_message_defined_by_user(element, default_coverage_value))
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
				
				meta_sample = manageDatabase.get_project_sample_metakey_last(project_sample, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue)
				if (meta_sample != None):
					manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Error, meta_sample.description)
				process_SGE.set_process_controler(user, process_controler.get_name_project_sample(project_sample), ProcessControler.FLAG_ERROR)
				return False
			
			#####################
			###
			### make mask the consensus SoftwareNames.SOFTWARE_MSA_MASKER
			limit_to_mask_consensus = int(default_software.get_mask_consensus_single_parameter(project_sample,\
							DefaultProjectSoftware.MASK_CONSENSUS_threshold, SoftwareNames.TECHNOLOGY_minion))
			msa_parameters = self.software.make_mask_consensus( 
				project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, self.software_names.get_medaka_name()), 
				project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT),
				project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ, self.software_names.get_medaka_name()),
				coverage, project_sample.sample.name, limit_to_mask_consensus)
			result_all.add_software(SoftwareDesc(self.software_names.get_msa_masker_name(), self.software_names.get_msa_masker_version(),\
					"{}; for coverages less than {} in {}% of the regions.".format(msa_parameters,\
					10,									
					100 - limit_to_mask_consensus) ))
			
			## identify VARIANTS IN INCOMPLETE LOCUS in all locus, set yes in variants if are in areas with coverage problems
			parse_out_files = ParseOutFiles()
			parse_out_files.add_variants_in_incomplete_locus(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, 
										self.software_names.get_medaka_name()), coverage)
			
			### draw coverage
			try:
				### make the coverage images
				draw_all_coverage = DrawAllCoverage()
				draw_all_coverage.draw_all_coverages(project_sample, SoftwareNames.SOFTWARE_Medaka_name)
			except:
				result = Result()
				result.set_error("Fail to draw coverage images")
				result.add_software(SoftwareDesc('In house software', '1.0', ''))
				manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Error, result.to_json())
				
				### get again and set error
				project_sample = ProjectSample.objects.get(pk=project_sample.id)
				project_sample.is_error = True
				project_sample.save()
				
				meta_sample = manageDatabase.get_project_sample_metakey_last(project_sample, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue)
				if (meta_sample != None):
					manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Error, meta_sample.description)
				process_SGE.set_process_controler(user, process_controler.get_name_project_sample(project_sample), ProcessControler.FLAG_ERROR)
				return False

			## count hits from tab file
			count_hits = CountHits()
			if (not out_put_path is None):
				file_tab = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, self.software_names.get_medaka_name())
				if (os.path.exists(file_tab)):
					vect_count_type = ['snp']	## only detects snp
					count_hits = self.utils.count_hits_from_tab(file_tab, vect_count_type)
					### set flag that is finished
					manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Count_Hits, MetaKeyAndValue.META_VALUE_Success, count_hits.to_json())

			### mixed infection
			try:
				## get instances
				mixed_infections_management = MixedInfectionsManagement()
				
				## set the alert also
				mixed_infection = mixed_infections_management.get_mixed_infections(project_sample, user, count_hits)
			except:
				result = Result()
				result.set_error("Fail to calculate mixed infection")
				result.add_software(SoftwareDesc('In house software', '1.0', ''))
				manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Mixed_Infection, MetaKeyAndValue.META_VALUE_Error, result.to_json())
				
				### get again and set error
				project_sample = ProjectSample.objects.get(pk=project_sample.id)
				project_sample.is_error = True
				project_sample.save()
				
				meta_sample = manageDatabase.get_project_sample_metakey_last(project_sample, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue)
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
			project_sample.is_mask_consensus_sequences = True
			project_sample.count_variations = manage_database.get_variation_count(count_hits)
			project_sample.mixed_infections = mixed_infection
			project_sample.save()
			
			### add today date, last change
			project = project_sample.project
			project.last_change_date = datetime.datetime.now()
			project.save()
			
			### get clean consensus file
			consensus_fasta = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, SoftwareNames.SOFTWARE_Medaka_name)
			if (os.path.exists(consensus_fasta)):
				file_out = project_sample.get_consensus_file(TypePath.MEDIA_ROOT)		### this is going to main path
				self.utils.filter_fasta_all_sequences_file(consensus_fasta, coverage, file_out, limit_to_mask_consensus, False)
			
			### set the tag of result OK 
			manageDatabase.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_Medaka, MetaKeyAndValue.META_VALUE_Success, result_all.to_json())
			
			### set the flag of the end of the task		
			meta_sample = manageDatabase.get_project_sample_metakey_last(project_sample, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue)
			if (meta_sample != None):
				manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Success, meta_sample.description)
		except Exception as e:
			## finished with error
			process_SGE.set_process_controler(user, process_controler.get_name_project_sample(project_sample), ProcessControler.FLAG_ERROR)
			return False
		
		### finished
		process_SGE.set_process_controler(user, process_controler.get_name_project_sample(project_sample), ProcessControler.FLAG_FINISHED)
		return True


	def run_medaka(self, file_fastq, reference_fasta, reference_gbk, sample_name, parameters):
		"""
		run medaka
		return output directory of snippy, try to do something most close possible
		
		
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
		
		:param
		(default: r941_min_high_g360).
        Available: r103_min_high_g345, r103_min_high_g360, r103_prom_high_g360, r103_prom_snp_g3210, r103_prom_variant_g3210,
        	r10_min_high_g303, r10_min_high_g340, r941_min_fast_g303, r941_min_high_g303, r941_min_high_g330,
        	r941_min_high_g340_rle, r941_min_high_g344, r941_min_high_g351, r941_min_high_g360, r941_prom_fast_g303,
        	r941_prom_high_g303, r941_prom_high_g330, r941_prom_high_g344, r941_prom_high_g360, r941_prom_high_g4011,
        	r941_prom_snp_g303, r941_prom_snp_g322, r941_prom_snp_g360, r941_prom_variant_g303, r941_prom_variant_g322,
        	r941_prom_variant_g360.
		### output medaka files
		
		bcftools   1.9    
		bgzip      1.9    
		minimap2   2.11     
		samtools   1.9    
		tabix      1.9
		
		consensus.fasta
		calls_to_draft.bam
		
		"""
		### make medaka consensus and bam
		temp_dir = os.path.join(self.utils.get_temp_dir(), sample_name)
		cmd = "{} {}_consensus -i {} -d {} -o {} -t {} {}".format(
				self.software_names.get_medaka_env(),
				self.software_names.get_medaka(), file_fastq, reference_fasta,
				temp_dir, settings.THREADS_TO_RUN_SLOW, parameters)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			self.utils.remove_dir(temp_dir)
			raise Exception("Fail to run medaka_consensus")

		### test output files
		hdf_file = os.path.join(temp_dir, "consensus_probs.hdf")
		bam_file = os.path.join(temp_dir, "calls_to_draft.bam")
		consensus_file = os.path.join(temp_dir, "consensus.fasta")
		for file_to_test in [hdf_file, bam_file, consensus_file]:
			if (not os.path.exists(file_to_test)):
				message = "File '{}' not found after medaka_consensus process: ".format(file_to_test) + cmd
				self.logger_production.error(message)
				self.logger_debug.error(message)
				self.utils.remove_dir(temp_dir)
				raise Exception(message)
		
		### change bam file names
		self.utils.move_file(bam_file, os.path.join(os.path.dirname(bam_file), sample_name + ".bam"))
		self.utils.move_file(bam_file + ".bai", os.path.join(os.path.dirname(bam_file), sample_name + ".bam.bai"))
		self.utils.move_file(consensus_file, os.path.join(os.path.dirname(bam_file), sample_name + ".consensus.fa"))
		bam_file = os.path.join(os.path.dirname(bam_file), sample_name + ".bam")
				
		### create depth
		depth_file = os.path.join(temp_dir, sample_name + ".depth.gz")
##		cmd =  "{} depth -aa -q 10 {} | {} -c > {}".format(   ### with quality
		cmd =  "{} depth -aa {} | {} -c > {}".format(
			self.software_names.get_samtools(),
			bam_file,
			self.software_names.get_bgzip(), 
			depth_file) 
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			self.utils.remove_dir(temp_dir)
			raise Exception("Fail to run samtools depth in nanopore")
		
		### create tbi
		self.software.create_index_with_tabix(depth_file, "bed")
		
		### vcf
		vcf_before_file = os.path.join(temp_dir, sample_name + "_before_annotation.vcf")
		cmd =  "{} {} variant --verbose {} {} {};".format(
#		cmd =  "{} {} snp --verbose {} {} {}".format(
			self.software_names.get_medaka_env(),
			self.software_names.get_medaka(),
			reference_fasta, hdf_file, vcf_before_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			self.utils.remove_dir(temp_dir)
			raise Exception("Fail to run medaka variant")
		
		
		### annotate vcf
		vcf_file = os.path.join(temp_dir, sample_name + ".vcf")
		cmd =  "{} {} tools annotate {} {} {} {}".format(
			self.software_names.get_medaka_env(),
			self.software_names.get_medaka(),
			vcf_before_file, reference_fasta, bam_file, vcf_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			self.utils.remove_dir(temp_dir)
			raise Exception("Fail to run medaka variant")
		
		### Add
		##INFO=<ID=SR,Number=.,Type=Integer,Description="Depth of spanning reads by strand which best align to each allele (ref fwd, ref rev, alt1 fwd, alt1 rev, etc.)">
		##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
		##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
		##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
		##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Number of observation for each allele">
		
		if (os.path.getsize(vcf_file) > 0):
			### Add ANN anotation with snpEff
			temp_file_2 = os.path.join(temp_dir, sample_name + '_ann.vcf')
			output_file = self.software.run_snpEff(reference_fasta, reference_gbk, vcf_file, temp_file_2)
			
			if not output_file is None:	## sometimes the gff does not have amino sequences
				self.utils.copy_file(output_file, vcf_file)
				
			### add FREQ to VCF file
			final_vcf = os.path.join(temp_dir, sample_name + '_2.vcf')
			self.utils.add_freq_ao_ad_and_type_to_vcf(vcf_file, final_vcf)
			self.utils.move_file(final_vcf, vcf_file)
			self.software.test_bgzip_and_tbi_in_vcf(vcf_file)
			
			### create TAB file
			self.software.run_snippy_vcf_to_tab_freq_and_evidence(reference_fasta,\
								reference_gbk, vcf_file,\
								os.path.join(temp_dir, sample_name + '.tab'))
		return temp_dir
		
		
		
		
		
		