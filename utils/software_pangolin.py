'''
Created on 19/04/2021

@author: mmp
'''
import os, logging
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
				(pangolin_version, pangolin_learn_version) = self._run_pangolin_update()
				software = SoftwareModel()
				software.name = SoftwareNames.SOFTWARE_Pangolin_name
				software.set_dual_version(pangolin_version, pangolin_learn_version)
				software.save()
				return
			
			## return if the software was updated today
			if (software.is_updated_today()): 
				return software.get_dual_version()
			
			(pangolin_version, pangolin_learn_version) = self._run_pangolin_update()
			software.set_dual_version(pangolin_version, pangolin_learn_version)
			software.set_last_update_today()
			software.save()
		
		return (pangolin_version, pangolin_learn_version)

	def _run_pangolin_update(self):
		"""
		Run update
		:return (pangolin version, pangolinLearn Version)
		"""
		## default version
		(pangolin_version, pangolin_learn_version) = self.software_names.get_pangolin_version(),\
				self.software_names.get_pangolin_learn_version()
		
		temp_file = self.utils.get_temp_file("pangolin_verion", ".txt")
		cmd = "{} {} --update > {} 2>&1".format(self.software_names.get_pangolin_env(),
			self.software_names.get_pangolin(), temp_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run pangolin --update")
		
		vect_lines = self.utils.read_text_file(temp_file)
		#$ pangolin --update
		#pangolin updated to v2.3.8
		#pangoLEARN updated to 2021-04-01
		#$ pangolin --update
		#pangolin already latest release (v2.3.8)
		#pangoLEARN already latest release (2021-04-01)
		for line in vect_lines:
			lst_data = line.strip().split()
			if (lst_data[0] == "pangolin"): pangolin_version = lst_data[-1].replace('(', '').replace(')', '')
			if (lst_data[0] == "pangoLEARN"): pangolin_learn_version = lst_data[-1].replace('(', '').replace(')', '')
		self.utils.remove_file(temp_file)
		return (pangolin_version, pangolin_learn_version)


	def run_pangolin(self, file_in_fasta, file_out_csv):
		"""
		Identifies pangoiln
		"""
		
		### update pangolin, if necessary
		self.run_pangolin_update()
		
		tem_dir = self.utils.get_temp_dir()
		
		### run pangolin
		cmd = "{} {} {} -o {}".format(self.software_names.get_pangolin_env(),
			self.software_names.get_pangolin(), file_in_fasta, tem_dir)
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
		software = Software()
		### test id abricate has the database
		try:
			uploadFile = UploadFile.objects.order_by('-version')[0]
		except UploadFile.DoesNotExist:
			return False

		## if not exist try to create it
		if (not software.is_exist_database_abricate(uploadFile.abricate_name)):
			try:
				self.create_database_abricate(uploadFile.abricate_name, uploadFile.path)
			except Exception:
				return False
		
		## run abricate
		out_file_abricate = self.utils.get_temp_file("temp_abricate", ".txt")
		try:
			cmd = software.run_abricate(uploadFile.abricate_name, fasta_ref_name, 
				SoftwareNames.SOFTWARE_ABRICATE_PARAMETERS, out_file_abricate)
		except Exception:
			return False
		
		if (not os.path.exists(out_file_abricate)):
			return False

		parseOutFiles = ParseOutFiles()
		(vect_data, clean_abricate_file) = parseOutFiles.parse_abricate_file(out_file_abricate,
				"doesnt_mather.txt",\
				SoftwareNames.SOFTWARE_SPAdes_CLEAN_HITS_BELLOW_VALUE)
		
		uploadFiles = UploadFiles()
		vect_data = uploadFiles.uploadIdentifyVirus(vect_data, uploadFile.abricate_name)
		if (len(vect_data) == 0):
			return False
		
		### test number of right segments
		number_right =0 
		for identify_virus in vect_data:
			if (identify_virus.seq_virus.name == "BetaCoV" and
				identify_virus.seq_virus.kind_type.name == "Genus"): number_right += 1
			elif (identify_virus.seq_virus.name == "SARS_CoV_2" and
				identify_virus.seq_virus.kind_type.name == "Human"): number_right += 1

		## if right at least two		
		if (number_right > 1): return True
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
		(pangolin_version, pangolin_learn_version) = self.run_pangolin_update()
		
		## get version
		meta_sample = manage_database.get_project_metakey_last(project,
							MetaKeyAndValue.META_KEY_Identify_pangolin,\
							MetaKeyAndValue.META_VALUE_Success)
		if (not meta_sample is None):
			software_names = SoftwareNames()
			result_pangolin = decode_result.decode_result(meta_sample.description)
			pangolin_version_run = result_pangolin.get_software_version(software_names.get_pangolin_name())
			pangolin_learn_version_run = result_pangolin.get_software_version(software_names.get_pangolin_learn_name())
			if (not pangolin_version_run is None and not pangolin_learn_version_run is None and \
				pangolin_version == pangolin_version_run and pangolin_learn_version == pangolin_learn_version_run):
				return False
		return True

	def get_update_message(self, project):
		"""
		"""
		decode_result = DecodeObjects()
		manage_database = ManageDatabase()
		
		## run always this to get last version
		(pangolin_version, pangolin_learn_version) = self.run_pangolin_update()
		
		## get version
		meta_sample = manage_database.get_project_metakey_last(project,
							MetaKeyAndValue.META_KEY_Identify_pangolin,\
							MetaKeyAndValue.META_VALUE_Success)
		pangolin_version_run = None
		pangolin_learn_version_run = None
		if (not meta_sample is None):
			software_names = SoftwareNames()
			result_pangolin = decode_result.decode_result(meta_sample.description)
			pangolin_version_run = result_pangolin.get_software(software_names.get_pangolin_name())
			pangolin_learn_version_run = result_pangolin.get_software(software_names.get_pangolin_learn_name())
			if (not pangolin_version_run is None and not pangolin_learn_version_run is None and \
				pangolin_version == pangolin_version_run and pangolin_learn_version == pangolin_learn_version_run):
				return False
		sz_out = ""
		if (pangolin_version_run is None):
			sz_out += "Do you want to run pangolin version ({})?".format(pangolin_version)
		elif (pangolin_version != pangolin_version_run):
			sz_out += "Do you want to update pangolin version from ({}) to ({})?".format(pangolin_version_run, pangolin_version)
		if (pangolin_learn_version_run is None):
			if (len(sz_out) > 0): sz_out += "<br /><br />"
			sz_out += "Do you want to run pangolinLearn version ({})?".format(pangolin_learn_version)
		elif (pangolin_learn_version != pangolin_learn_version_run):
			if (len(sz_out) > 0): sz_out += "<br /><br />"
			sz_out += "Do you want to update pangolinLearn version from ({}) to ({})?".format(pangolin_learn_version_run, pangolin_learn_version)
		return sz_out

