'''
Created on 19/04/2021

@author: mmp
'''
import os, logging
from django.conf import settings
from managing_files.models import Software as SoftwareModel
from utils.lock_atomic_transaction import LockedAtomicTransaction
from constants.software_names import SoftwareNames
from manage_virus.models import UploadFile
from manage_virus.uploadFiles import UploadFiles
from utils.parse_out_files import ParseOutFiles
from utils.software import Software
from managing_files.manage_database import ManageDatabase
from constants.meta_key_and_values import MetaKeyAndValue
from utils.result import DecodeObjects
from utils.utils import Utils

class SoftwarePangolin(object):
	'''
	classdocs
	'''
	utils = Utils()
	software_names = SoftwareNames()

	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	def run_pangolin_update(self):
		"""
		Identifies pangolin, locks the table and 
		"""
		#### test if there's a new update
		# Lock the table, and with this approach only one job can run at the time
		with LockedAtomicTransaction(SoftwareModel):
			try:
				software = SoftwareModel.objects.get(name=SoftwareNames.SOFTWARE_Pangolin_name)
			except SoftwareModel.DoesNotExist:	## need to create with last version
				
				### run pangolin update first time
				dt_result_version = self._run_pangolin_update()
				software = SoftwareModel()
				software.name = SoftwareNames.SOFTWARE_Pangolin_name
				software.set_version(dt_result_version.get(SoftwareNames.SOFTWARE_Pangolin_name.lower(), ""))
				software.set_version_long(dt_result_version)
				software.save()
				return
			
			## return if the software was updated today
			if (software.is_updated_today()): 
				return software.get_version_long()
			
			dt_result_version = self._run_pangolin_update()
			if len(dt_result_version) > 0:
				software.set_version_long(dt_result_version)
				software.set_version(dt_result_version.get(SoftwareNames.SOFTWARE_Pangolin_name.lower(), ""))
				software.set_last_update_today()
				software.save()
			else:
				dt_result_version = software.get_version_long()
		return dt_result_version

	def _run_pangolin_update(self):
		"""
		Run update
		:return (pangolin version, pangolinLearn Version)
		"""
		dt_result_version = {}
		if (not settings.DEBUG):
			## default version
			temp_file = self.utils.get_temp_file("pangolin_verion", ".txt")
			cmd = "{} {} --update > {} 2>&1".format(self.software_names.get_pangolin_env(),
				self.software_names.get_pangolin(), temp_file)
			exist_status = os.system(cmd)
			if (exist_status != 0):
				self.logger_production.error('Fail to run: ' + cmd)
				self.logger_debug.error('Fail to run: ' + cmd)
				raise Exception("Fail to run pangolin --update")
			
			### check all versions installed 
# cmd = "{} {} --all-versions > {} 2>&1".format(self.software_names.get_pangolin_env(),
# 	self.software_names.get_pangolin(), temp_file)
# exist_status = os.system(cmd)
# if (exist_status != 0):
# 	self.logger_production.error('Fail to run: ' + cmd)
# 	self.logger_debug.error('Fail to run: ' + cmd)
# 	raise Exception("Fail to run pangolin --all-versions")

			vect_lines = self.utils.read_text_file(temp_file)
			## Important, go through pangolin --all-versions
			#$ pangolin --all-versions	## now aldo uses conda
			# pangolin: 4.0.5
			# pangolin-data: 1.3
			# constellations: v0.1.6
			# scorpio: 0.3.16
			
			for line in vect_lines:
				## must be lower case
				lst_data = [_.lower().strip() for _ in line.strip().split()]
				if len(lst_data) > 1:
					dt_result_version[lst_data[0]] = self._get_verion_tag(lst_data[-1].replace('(', '').replace(')', ''))
			self.utils.remove_file(temp_file)
		return dt_result_version


	def _get_verion_tag(self, description):
		"""
		For "pango-designation (vv1.2.13) is newer than latest stable release (v1.2.86), not updating." return (vv1.2.13)
		For "pangolearn updated to 2021-09-28" return (2021-09-28)
		For "pangolin already latest release (v3.1.14)" return (v3.1.14)
		For "pango-designation used by pangoLEARN/Usher: v1.2.101" return (v1.2.101)
		For "pango-designation aliases: 1.2.106" return (v1.2.101)
		"""
		for word in description.split():
			count_word = 0
			for char_ in word:
				if self.utils.is_integer(char_): count_word += 1
				if (count_word > 1):
					return word
		return ""
	
	def run_pangolin(self, file_in_fasta, file_out_csv):
		"""
		Identifies pangolin
		"""
		
		### update pangolin, if necessary
		self.run_pangolin_update()
		
		tem_dir = self.utils.get_temp_dir()
		
		### run pangolin
		cmd = "{} {} {} -o {} --analysis-mode {} -t 2".format(self.software_names.get_pangolin_env(),
			self.software_names.get_pangolin(), file_in_fasta, tem_dir, settings.RUN_PANGOLIN_MODEL)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run pangolin identifcation")
		
		output_file_expected = os.path.join(tem_dir, "lineage_report.csv")
		if (os.path.exists(output_file_expected)):
			self.utils.move_file(output_file_expected, file_out_csv)
			self.utils.remove_dir(tem_dir)
			return file_out_csv
		self.utils.remove_dir(tem_dir)
		return None

	def is_ref_sars_cov_2(self, fasta_ref_name):
		"""
		test if it's SARS_COV
		"""

		#self.logger_debug.info("is_ref_sars_cov_2 {}".format(fasta_ref_name))

		software = Software()
		### test id abricate has the database
		try:
			uploadFile = UploadFile.objects.order_by('-version')[0]
		except UploadFile.DoesNotExist as e:
			#self.logger_debug.error("is_ref_sars_cov_2 getting db {}".format(e))
			return False

		#self.logger_debug.info("is_ref_sars_cov_2 {}".format(uploadFile.abricate_name))

		## if not exist try to create it
		if (not software.is_exist_database_abricate(uploadFile.abricate_name)):
			#self.logger_debug.info("is_ref_sars_cov_2 {} db does not exist".format(uploadFile.abricate_name))
			try:
				#self.logger_debug.info("is_ref_sars_cov_2 creating db for {}".format(uploadFile.abricate_name))
				software.create_database_abricate(uploadFile.abricate_name, uploadFile.path)
				#self.logger_debug.info("is_ref_sars_cov_2 created db for {}".format(uploadFile.abricate_name))				
			except Exception as e:
				#self.logger_debug.error("is_ref_sars_cov_2 failed creating db {}".format(e))
				return False
		
		## run abricate
		out_file_abricate = self.utils.get_temp_file("temp_abricate", ".txt")
		try:
			#self.logger_debug.info("is_ref_sars_cov_2 going to run abricate")				
			software.run_abricate(uploadFile.abricate_name, fasta_ref_name, 
				SoftwareNames.SOFTWARE_ABRICATE_PARAMETERS, out_file_abricate)
			#self.logger_debug.info("is_ref_sars_cov_2 ran abricate")
		except Exception as e:
			#self.logger_debug.error("is_ref_sars_cov_2 failed running abricate {}".format(e))
			return False
		
		if (not os.path.exists(out_file_abricate)):
			#self.logger_debug.info("is_ref_sars_cov_2 could not find {}".format(out_file_abricate))
			return False

		parseOutFiles = ParseOutFiles()
		(vect_data, clean_abricate_file) = parseOutFiles.parse_abricate_file(out_file_abricate,
				"doesnt_mather.txt",\
				SoftwareNames.SOFTWARE_SPAdes_CLEAN_HITS_BELLOW_VALUE)
		## don't need this file, only need vect_data_file
		self.utils.remove_file(clean_abricate_file)
		self.utils.remove_file(out_file_abricate)
		
		uploadFiles = UploadFiles()
		try:
			#self.logger_debug.info("is_ref_sars_cov_2 trying to identify")
			## could fail some identification, so, return False if it occurs
			vect_data = uploadFiles.uploadIdentifyVirus(vect_data, uploadFile.abricate_name)
			if (len(vect_data) == 0):
				#self.logger_debug.info("is_ref_sars_cov_2 failed to identify")
				return False
		except Exception as e:
			#self.logger_debug.info("is_ref_sars_cov_2 error trying to identify {}".format(e))
			return False
		
		### test number of right segments
		number_right_bet_cov =0 
		for identify_virus in vect_data:
			if (identify_virus.seq_virus.name == "BetaCoV" and
				identify_virus.seq_virus.kind_type.name == "Genus"): number_right_bet_cov += 1
			### need to read from the file 
			elif (identify_virus.seq_virus.name in ("SARS_CoV_2", "SARS_CoV", "SCoV2_potential_Omicron", "HCoV_OC43", "HCoV_HKU1", "MERS_CoV") and
				identify_virus.seq_virus.kind_type.name == "Human"): number_right_bet_cov += 1

		## if right at least two		
		if (number_right_bet_cov > 1): return True
		return False


	def pangolin_results_out_date(self, project):
		"""
		"""
		### test if there any collect data running
		decode_result = DecodeObjects()
		manage_database = ManageDatabase()
		metaKeyAndValue = MetaKeyAndValue()
		meta_project_queue = manage_database.get_project_metakey_last(project, metaKeyAndValue.get_meta_key(\
						MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id), MetaKeyAndValue.META_VALUE_Queue)
		
		if (not meta_project_queue is None):
			meta_project_success = manage_database.get_project_metakey_last(project, metaKeyAndValue.get_meta_key(\
						MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id), MetaKeyAndValue.META_VALUE_Success)
			if (not meta_project_success is None and meta_project_success.creation_date < meta_project_queue.creation_date):
				return False	## can't run because there are at least one in queue
		
		
		## run always this to get last version
		dt_pangolin_versions = self.run_pangolin_update()
		
		## get version
		meta_sample = manage_database.get_project_metakey_last(project,
							MetaKeyAndValue.META_KEY_Identify_pangolin,\
							MetaKeyAndValue.META_VALUE_Success)
		if (not meta_sample is None):
			result_pangolin = decode_result.decode_result(meta_sample.description)
			for soft_name in dt_pangolin_versions:
				## if actual version does not have this software it is because out of date
				version_run_before = result_pangolin.get_software_version(soft_name)
				if version_run_before is None: return True
				if version_run_before != dt_pangolin_versions[soft_name]: return True
			return False
		return True

	def get_update_message(self, project):
		"""
		"""
		decode_result = DecodeObjects()
		manage_database = ManageDatabase()
		
		## run always this to get last version
		dt_pangolin_versions = self.run_pangolin_update()
		
		## get version
		meta_sample = manage_database.get_project_metakey_last(project,
							MetaKeyAndValue.META_KEY_Identify_pangolin,\
							MetaKeyAndValue.META_VALUE_Success)
		
		sz_out = ""
		if (not meta_sample is None):
			result_pangolin = decode_result.decode_result(meta_sample.description)
			for soft_name in self.__get_order_pangolin_names(dt_pangolin_versions):
				## if actual version does not have this software it is because out of date
				version_run_before = result_pangolin.get_software_version(soft_name)
				if version_run_before is None:
					if (len(sz_out) > 0): sz_out += "<br /><br />"
					sz_out += "Do you want to run {} version ({})?".format(soft_name, dt_pangolin_versions.get(soft_name))
				elif (version_run_before != dt_pangolin_versions.get(soft_name)):
					if (len(sz_out) > 0): sz_out += "<br /><br />"
					sz_out += "Do you want to update {} version from ({}) to ({})?".format(soft_name, version_run_before,
									dt_pangolin_versions.get(soft_name))
		else:
			for soft_name in dt_pangolin_versions:
				if (len(sz_out) > 0): sz_out += "<br /><br />"
				sz_out += "Do you want to run {} version ({})?".format(soft_name, dt_pangolin_versions.get(soft_name))

		return sz_out

	def __get_order_pangolin_names(self, dt_pangolin_versions):
		"""
		Return Pangolin first, then pangolinXXX, at last others
		"""
		vect_out = []
		if SoftwareNames.SOFTWARE_Pangolin_name.lower() in dt_pangolin_versions:
			vect_out = [SoftwareNames.SOFTWARE_Pangolin_name.lower()]
		for soft_name in dt_pangolin_versions:
			if soft_name.startswith(SoftwareNames.SOFTWARE_Pangolin_name.lower()) and \
				not soft_name in vect_out:
				vect_out.append(soft_name)
		for soft_name in dt_pangolin_versions:
			if not soft_name in vect_out: vect_out.append(soft_name)
		return vect_out
		
		
		
