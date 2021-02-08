'''
Created on 03/05/2020

@author: mmp
'''
from settings.models import Software, Parameter
from constants.software_names import SoftwareNames
from utils.lock_atomic_transaction import LockedAtomicTransaction
from settings.default_software import DefaultSoftware
from managing_files.models import Project, ProjectSample, Sample
		
class DefaultProjectSoftware(object):
	'''
	classdocs
	'''
	software_names = SoftwareNames()

	### used in snippy
	SNIPPY_COVERAGE_NAME = "--mincov"
	SNIPPY_MAPQUAL_NAME = "--mapqual"

	### used in NANOfilt
	NANOfilt_quality_read = "-q"
	
	### used in Medaka
	MEDAKA_model = "-m"
	
	### used in mask consensus
	MASK_CONSENSUS_threshold = "Threshold"

	def __init__(self):
		""" change values """
		self.change_values_software = {}	### the key is the name of the software
		
	def test_all_defaults(self, user, type_of_use, project, project_sample, sample):
		"""
		test all defaults for all software available
		"""
		## only for project and for all technology
		if (not project is None):
			self.test_default_db(SoftwareNames.SOFTWARE_SNIPPY_name, user, type_of_use, project,
						project_sample, None, SoftwareNames.TECHNOLOGY_illumina)
	## Not in used yet
	##		self.test_default_db(SoftwareNames.SOFTWARE_FREEBAYES_name, user, project)
	
			self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,\
							user, type_of_use, project, project_sample, None, SoftwareNames.TECHNOLOGY_illumina)
			self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,\
							user, type_of_use, project, project_sample, None, SoftwareNames.TECHNOLOGY_minion)
			self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,\
							user, type_of_use, project, project_sample, None, SoftwareNames.TECHNOLOGY_minion)
			self.test_default_db(SoftwareNames.SOFTWARE_Medaka_name_consensus,\
							user, type_of_use, project, project_sample, None, SoftwareNames.TECHNOLOGY_minion)
			self.test_default_db(SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT,\
							user, type_of_use, project, project_sample, None, SoftwareNames.TECHNOLOGY_minion)
		## only for project sample and by technology
		elif (not project_sample is None):
			
			if (project_sample.sample.is_type_fastq_gz_sequencing()): ### illumina
				self.test_default_db(SoftwareNames.SOFTWARE_SNIPPY_name, user, type_of_use, project,
						project_sample, None, SoftwareNames.TECHNOLOGY_illumina)
	## Not in used yet
	##		self.test_default_db(SoftwareNames.SOFTWARE_FREEBAYES_name, user, project)
	
				self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,\
							user, type_of_use, project, project_sample, None, SoftwareNames.TECHNOLOGY_illumina)
			else:
				self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,\
								user, type_of_use, project, project_sample, None, SoftwareNames.TECHNOLOGY_minion)
				self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,\
								user, type_of_use, project, project_sample, None, SoftwareNames.TECHNOLOGY_minion)
				self.test_default_db(SoftwareNames.SOFTWARE_Medaka_name_consensus,\
								user, type_of_use, project, project_sample, None, SoftwareNames.TECHNOLOGY_minion)
				self.test_default_db(SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT,\
								user, type_of_use, project, project_sample, None, SoftwareNames.TECHNOLOGY_minion)
		
		### only for sample and ONT technology
		if (not sample is None and sample.is_type_fastq_gz_sequencing(Sample.TYPE_OF_FASTQ_minion)):
			self.test_default_db(SoftwareNames.SOFTWARE_NanoFilt_name,\
						user, Software.TYPE_OF_USE_sample, None, None, sample, SoftwareNames.TECHNOLOGY_minion)


	def test_default_db(self, software_name, user, type_of_use, project, project_sample, 
					sample, technology_name):
		"""
		test if exist, if not persist in database
		"""
		default_software = DefaultSoftware()
		list_software = Software.objects.filter(name=software_name, owner=user, type_of_use=type_of_use,
					parameter__project=project, parameter__project_sample=project_sample,
					parameter__sample=sample,
					technology__name = technology_name).distinct("name")
		if len(list_software) == 0:
			### if it is Minion is because that does not exist at all. 
			### Previous versions didn't have TechnologyName
			vect_parameters = self._get_default_parameters(software_name, user, type_of_use, project,
					project_sample, sample, technology_name)
			if (technology_name == SoftwareNames.TECHNOLOGY_illumina):
				list_software = Software.objects.filter(name=software_name, owner=user,
					type_of_use = type_of_use,
					parameter__project=project,
					parameter__project_sample=project_sample).distinct("name")

				### if exist set illumina in technology					
				if (len(list_software) == 1):
					list_software[0].technology = default_software.get_technology_instance(technology_name)
					list_software[0].save()
				else:
					with LockedAtomicTransaction(Software), LockedAtomicTransaction(Parameter):
						self._persist_parameters(vect_parameters, type_of_use, technology_name)
			else:			
				### create a default one for this user
				with LockedAtomicTransaction(Software), LockedAtomicTransaction(Parameter):
					### persist 
					if (len(vect_parameters) > 0):
						self._persist_parameters(vect_parameters, type_of_use, technology_name)


	def _get_default_parameters(self, software_name, user, type_of_use, project, project_sample,
					sample, technology_name):
		if (software_name == SoftwareNames.SOFTWARE_SNIPPY_name):
			vect_parameters = self._get_snippy_default(user, type_of_use, project, project_sample)		### base values
			if (not project is None): vect_parameters = self._get_default_project(user,\
					SoftwareNames.SOFTWARE_SNIPPY_name, None, vect_parameters,
					technology_name)		### base values
			if (not project_sample is None): vect_parameters = self._get_default_project(user,\
					SoftwareNames.SOFTWARE_SNIPPY_name, project_sample.project, vect_parameters,
					technology_name)		### base values
			return vect_parameters
		elif (software_name == SoftwareNames.SOFTWARE_FREEBAYES_name):
			return self._get_freebayes_default(user, type_of_use, project, project_sample)
		elif (software_name == SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name):
			vect_parameters = self._get_mask_consensus_threshold_default(user, type_of_use, project, project_sample)
			if (not project is None): vect_parameters = self._get_default_project(user,\
				software_name, None, vect_parameters,
				technology_name)		### base values
			if (not project_sample is None): vect_parameters = self._get_default_project(user,\
				software_name, project_sample.project,
				vect_parameters, technology_name)		### base values
			return vect_parameters
		elif (software_name == SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name):
			vect_parameters = self._get_limit_coverage_ONT_threshold_default(user, type_of_use, project, project_sample)
			if (not project is None): vect_parameters = self._get_default_project(user,\
				software_name, None, vect_parameters,
				technology_name)		### base values
			if (not project_sample is None): vect_parameters = self._get_default_project(user,\
				software_name, project_sample.project,
				vect_parameters, technology_name)		### base values
			return vect_parameters
		elif (software_name == SoftwareNames.SOFTWARE_Medaka_name_consensus):
			vect_parameters = self._get_medaka_model_default(user, type_of_use, project, project_sample)
			if (not project is None): vect_parameters = self._get_default_project(user,\
				software_name, None, vect_parameters,
				technology_name)		### base values
			if (not project_sample is None): vect_parameters = self._get_default_project(user,\
				software_name, project_sample.project,
				vect_parameters, technology_name)		### base values
			return vect_parameters
		elif (software_name == SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT):
			vect_parameters = self._get_samtools_depth_default_ONT(user, type_of_use, project, project_sample)
			if (not project is None): vect_parameters = self._get_default_project(user,\
				software_name, None, vect_parameters,
				technology_name)		### base values
			if (not project_sample is None): vect_parameters = self._get_default_project(user,\
				software_name, project_sample.project,
				vect_parameters, technology_name)		### base values
			return vect_parameters
		elif (software_name == SoftwareNames.SOFTWARE_NanoFilt_name):
			vect_parameters = self._get_nanofilt_default(user, type_of_use, sample)
			vect_parameters = self._get_default_project(user,\
				software_name, sample,
				vect_parameters, technology_name)		### base values
			return vect_parameters		
		return []

	#####################################################
	#####
	#####		snippy
	#####
	
	def get_snippy_parameters(self, user, type_of_use, project, project_sample):
		"""
		get snippy parameters
		"""
		return self._get_parameters(SoftwareNames.SOFTWARE_SNIPPY_name, user, type_of_use,
			project, project_sample, None, SoftwareNames.TECHNOLOGY_illumina)
	
	def get_snippy_parameters_all_possibilities(self, user, project_sample):
		"""
		get snippy parameters for project_sample, project and default
		"""
		
		### Test project_sample first
		parameters = self._get_parameters(SoftwareNames.SOFTWARE_SNIPPY_name, user,\
				Software.TYPE_OF_USE_project_sample, None, project_sample, None,
				SoftwareNames.TECHNOLOGY_illumina)
		if (not parameters is None): return parameters
		
		### Test project
		parameters = self._get_parameters(SoftwareNames.SOFTWARE_SNIPPY_name, user,\
				Software.TYPE_OF_USE_project, project_sample.project, None, None,
				SoftwareNames.TECHNOLOGY_illumina)
		if (not parameters is None): return parameters
		
		### can be a default one
		default_software = DefaultSoftware()
		parameters = default_software.get_snippy_parameters(user)
		if (len(parameters) > 0): return parameters
		
		software_names = SoftwareNames()
		return software_names.get_snippy_parameters()
	
	
	def get_snippy_parameters_for_project(self, user, project):
		"""
		get snippy parameters only for project or default
		"""
		
		### Test project
		parameters = self._get_parameters(SoftwareNames.SOFTWARE_SNIPPY_name, user,\
				Software.TYPE_OF_USE_project, project, None, None, SoftwareNames.TECHNOLOGY_illumina)
		if (not parameters is None): return parameters
		
		### can be a default one
		default_software = DefaultSoftware()
		parameters = default_software.get_snippy_parameters(user)
		if (len(parameters) > 0): return parameters
		return None
	
	def get_snippy_single_parameter_default(self, parameter_name):
		"""
		:param parameter_name -> Only these two possibilities available SNIPPY_COVERAGE_NAME; SNIPPY_MAPQUAL_NAME
		:return value of parameter
		"""
		vect_parameters = self._get_snippy_default(None, None, None, None)
		for parameters in vect_parameters:
			if parameters.name == parameter_name:
				return parameters.parameter
		return None

	def is_snippy_single_parameter_default(self, project_sample, parameter_name):
		"""
		test if a specific parameter is default SNIPPY_COVERAGE_NAME; SNIPPY_MAPQUAL_NAME
		"""
		
		value_default_parameter = self.get_snippy_single_parameter_default(parameter_name)
		if (value_default_parameter is None): return False
		
		parameter_defined = self.get_snippy_single_parameter(project_sample, parameter_name)
		if not parameter_defined is None and parameter_defined == value_default_parameter: return True
		return False
		
	def get_snippy_single_parameter(self, project_sample, parameter_name):
		"""
		get snippy single parameters
		:param parameter_name -> Only these two possibilities available SNIPPY_COVERAGE_NAME; SNIPPY_MAPQUAL_NAME
		"""
		
		parameters_string = self.get_snippy_parameters_all_possibilities(project_sample.project.owner, project_sample)
		if (parameters_string is None): return None
		lst_data = parameters_string.split(parameter_name)
		if len(lst_data) == 2: return lst_data[1].split()[0]
		return None
	
	def is_snippy_single_parameter_default_for_project(self, project, parameter_name):
		"""
		test if a specific parameter is default SNIPPY_COVERAGE_NAME; SNIPPY_MAPQUAL_NAME
		"""
		
		value_default_parameter = self.get_snippy_single_parameter_default(parameter_name)
		if (value_default_parameter is None): return False
		
		parameter_defined = self.get_snippy_single_parameter_for_project(project, parameter_name)
		if not parameter_defined is None and parameter_defined == value_default_parameter: return True
		return False
	
	def get_snippy_single_parameter_for_project(self, project, parameter_name):
		"""
		get snippy single parameters
		:param parameter_name -> Only these two possibilities available SNIPPY_COVERAGE_NAME; SNIPPY_MAPQUAL_NAME
		"""
		
		parameters_string = self.get_snippy_parameters_for_project(project.owner, project)
		if (parameters_string is None): return None
		lst_data = parameters_string.split(parameter_name)
		if len(lst_data) == 2: return lst_data[1].split()[0]
		return None

	#####
	#####		END snippy
	#####
	#####################################################
	
	#####################################################
	#####
	#####		nanofilt
	#####
	
	def get_nanofilt_parameters(self, user, type_of_use, sample):
		"""
		get nanofilt parameters
		"""
		return self._get_parameters(SoftwareNames.SOFTWARE_NanoFilt_name, user, type_of_use,
			None, None, sample, SoftwareNames.TECHNOLOGY_minion)
	
	def get_nanofilt_parameters_all_possibilities(self, user, sample):
		"""
		get nanofilt parameters for project_sample, project and default
		"""
		
		### Test project_sample first
		parameters = self._get_parameters(SoftwareNames.SOFTWARE_NanoFilt_name, user,\
				Software.TYPE_OF_USE_sample, None, None, sample,
				SoftwareNames.TECHNOLOGY_minion)
		if (not parameters is None): return parameters
		
		### can be a default one
		default_software = DefaultSoftware()
		parameters = default_software.get_nanofilt_parameters(user)
		if (len(parameters) > 0): return parameters
		
		software_names = SoftwareNames()
		return software_names.get_nanofilt_parameters()
	
	
	def get_nanofilt_parameters_for_sample(self, user, sample):
		"""
		get nanofilt parameters only for project or default
		"""
		
		### Test project
		parameters = self._get_parameters(SoftwareNames.SOFTWARE_NanoFilt_name, user,\
				Software.TYPE_OF_USE_sample, None, None, sample, SoftwareNames.TECHNOLOGY_minion)
		if (not parameters is None): return parameters
		
		### can be a default one
		default_software = DefaultSoftware()
		parameters = default_software.get_nanofilt_parameters(user)
		if (len(parameters) > 0): return parameters
		return None
	
	def get_nanofilt_single_parameter_default(self, parameter_name):
		"""
		:param parameter_name -> NANOfilt_quality_read
		:return value of parameter
		"""
		vect_parameters = self._get_nanofilt_default(None, None, None)
		for parameters in vect_parameters:
			if parameters.name == parameter_name:
				return parameters.parameter
		return None

	def is_nanofilt_single_parameter_default(self, sample, parameter_name):
		"""
		test if a specific parameter is default NANOfilt_quality_read
		"""
		
		value_default_parameter = self.get_nanofilt_single_parameter_default(parameter_name)
		if (value_default_parameter is None): return False
		
		parameter_defined = self.get_nanofilt_single_parameter(sample, parameter_name)
		if not parameter_defined is None and parameter_defined == value_default_parameter: return True
		return False
		
	def get_nanofilt_single_parameter(self, sample, parameter_name):
		"""
		get nanofilt single parameters
		:param parameter_name -> Only these two possibilities available NANOfilt_quality_read
		"""
		
		parameters_string = self.get_nanofilt_parameters_all_possibilities(sample.owner, sample)
		if (parameters_string is None): return None
		lst_data = parameters_string.split(parameter_name)
		if len(lst_data) == 2: return lst_data[1].split()[0]
		return None
	
	def is_nanofilt_single_parameter_default_for_project(self, sample, parameter_name):
		"""
		test if a specific parameter is default NANOfilt_quality_read
		"""
		
		value_default_parameter = self.get_nanofilt_single_parameter_default(parameter_name)
		if (value_default_parameter is None): return False
		
		parameter_defined = self.get_nanofilt_single_parameter_for_sample(sample, parameter_name)
		if not parameter_defined is None and parameter_defined == value_default_parameter: return True
		return False
	
	def get_nanofilt_single_parameter_for_sample(self, sample, parameter_name):
		"""
		get nanofilt single parameters
		:param parameter_name -> Only these two possibilities available NANOfilt_quality_read
		"""
		
		parameters_string = self.get_nanofilt_parameters_for_sample(sample.owner, sample)
		if (parameters_string is None): return None
		lst_data = parameters_string.split(parameter_name)
		if len(lst_data) == 2: return lst_data[1].split()[0]
		return None

	#####
	#####		END nanofilt
	#####
	#####################################################
	
	#####################################################
	#####
	#####		Medaka
	#####
	
	def get_medaka_parameters(self, user, type_of_use, project, project_sample):
		"""
		get medaka parameters
		"""
		return self._get_parameters(SoftwareNames.SOFTWARE_Medaka_name_consensus, user, type_of_use,
			project, project_sample, None, SoftwareNames.TECHNOLOGY_minion)
	
	def get_medaka_parameters_all_possibilities(self, user, project_sample):
		"""
		get medaka parameters for project_sample, project and default
		get the model to run, this is the only possibility
		"""
		
		### Test project_sample first
		parameters = self._get_parameters(SoftwareNames.SOFTWARE_Medaka_name_consensus, user,\
				Software.TYPE_OF_USE_project_sample, None, project_sample, None,
				SoftwareNames.TECHNOLOGY_minion)
		if (not parameters is None): return parameters
		
		### Test project
		parameters = self._get_parameters(SoftwareNames.SOFTWARE_Medaka_name_consensus, user,\
				Software.TYPE_OF_USE_project, project_sample.project, None, None,
				SoftwareNames.TECHNOLOGY_minion)
		if (not parameters is None): return parameters
		
		### can be a default one
		default_software = DefaultSoftware()
		parameters = default_software.get_medaka_parameters_consensus(user)
		if (len(parameters) > 0): return parameters
		
		software_names = SoftwareNames()
		return software_names.get_medaka_parameters_consensus()
	
	
	def get_medaka_parameters_for_project(self, user, project):
		"""
		get medaka parameters only for project or default
		"""
		
		### Test project
		parameters = self._get_parameters(SoftwareNames.SOFTWARE_Medaka_name_consensus, user,\
				Software.TYPE_OF_USE_project, project, None, None, SoftwareNames.TECHNOLOGY_minion)
		if (not parameters is None): return parameters
		
		### can be a default one
		default_software = DefaultSoftware()
		parameters = default_software.get_medaka_parameters_consensus(user)
		if (len(parameters) > 0): return parameters
		return None
	
	def get_medaka_single_parameter_default(self, parameter_name):
		"""
		:param parameter_name -> Only these two possibilities available MEDAKA_model
		:return value of parameter
		"""
		vect_parameters = self._get_medaka_model_default(None, None, None, None)
		for parameters in vect_parameters:
			if parameters.name == parameter_name:
				return parameters.parameter
		return None

	def is_medaka_single_parameter_default(self, project_sample, parameter_name):
		"""
		test if a specific parameter is default MEDAKA_model
		"""
		
		value_default_parameter = self.get_medaka_single_parameter_default(parameter_name)
		if (value_default_parameter is None): return False
		
		parameter_defined = self.get_medaka_single_parameter(project_sample, parameter_name)
		if not parameter_defined is None and parameter_defined == value_default_parameter: return True
		return False
		
	def get_medaka_single_parameter(self, project_sample, parameter_name):
		"""
		get medaka single parameters
		:param parameter_name -> Only these two possibilities available MEDAKA_model
		"""
		
		parameters_string = self.get_medaka_parameters_all_possibilities(project_sample.project.owner, project_sample)
		if (parameters_string is None): return None
		lst_data = parameters_string.split(parameter_name)
		if len(lst_data) == 2: return lst_data[1].split()[0]
		return None
	
	def is_medaka_single_parameter_default_for_project(self, project, parameter_name):
		"""
		test if a specific parameter is default MEDAKA_model
		"""
		
		value_default_parameter = self.get_medaka_single_parameter_default(parameter_name)
		if (value_default_parameter is None): return False
		
		parameter_defined = self.get_medaka_single_parameter_for_project(project, parameter_name)
		if not parameter_defined is None and parameter_defined == value_default_parameter: return True
		return False
	
	def get_medaka_single_parameter_for_project(self, project, parameter_name):
		"""
		get medaka single parameters
		:param parameter_name -> Only these two possibilities available MEDAKA_model
		"""
		
		parameters_string = self.get_medaka_parameters_for_project(project.owner, project)
		if (parameters_string is None): return None
		lst_data = parameters_string.split(parameter_name)
		if len(lst_data) == 2: return lst_data[1].split()[0]
		return None

	#####
	#####		END medaka
	#####
	#####################################################
	
	#####################################################
	#####
	#####		Samtools ONT
	#####
	
	def get_samtools_parameters_ONT(self, user, type_of_use, project, project_sample):
		"""
		get samtools parameters
		"""
		return self._get_parameters(SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT, user, type_of_use,
			project, project_sample, None, SoftwareNames.TECHNOLOGY_minion)
	
	def get_samtools_parameters_all_possibilities_ONT(self, user, project_sample):
		"""
		get samtools parameters ONT for project_sample, project and default
		get the model to run, this is the only possibility
		"""
		
		### Test project_sample first
		parameters = self._get_parameters(SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT, user,\
				Software.TYPE_OF_USE_project_sample, None, project_sample, None,
				SoftwareNames.TECHNOLOGY_minion)
		if (not parameters is None): return parameters
		
		### Test project
		parameters = self._get_parameters(SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT, user,\
				Software.TYPE_OF_USE_project, project_sample.project, None, None,
				SoftwareNames.TECHNOLOGY_minion)
		if (not parameters is None): return parameters
		
		### can be a default one
		default_software = DefaultSoftware()
		parameters = default_software.get_samtools_parameters_depth_ONT(user)
		if (len(parameters) > 0): return parameters
		
		software_names = SoftwareNames()
		return software_names.get_samtools_parameters_depth_ONT()
	
	
	def get_samtools_parameters_for_project_ONT(self, user, project):
		"""
		get samtools parameters only for project or default
		"""
		
		### Test project
		parameters = self._get_parameters(SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT, user,\
				Software.TYPE_OF_USE_project, project, None, None, SoftwareNames.TECHNOLOGY_minion)
		if (not parameters is None): return parameters
		
		### can be a default one
		default_software = DefaultSoftware()
		parameters = default_software.get_samtools_parameters_depth_ONT(user)
		if (len(parameters) > 0): return parameters
		return None
	
	def get_samtools_single_parameter_default_ONT(self, parameter_name):
		"""
		:param parameter_name -> Only these two possibilities available Samtools ONT
		:return value of parameter
		"""
		vect_parameters = self._get_samtools_depth_default_ONT(None, None, None, None)
		for parameters in vect_parameters:
			if parameters.name == parameter_name:
				return parameters.parameter
		return None

	def is_samtools_single_parameter_default_ONT(self, project_sample, parameter_name):
		"""
		test if a specific parameter is default Samtools ONT
		"""
		
		value_default_parameter = self.get_samtools_single_parameter_default_ONT(parameter_name)
		if (value_default_parameter is None): return False
		
		parameter_defined = self.get_samtools_single_parameter_ONT(project_sample, parameter_name)
		if not parameter_defined is None and parameter_defined == value_default_parameter: return True
		return False
		
	def get_samtools_single_parameter_ONT(self, project_sample, parameter_name):
		"""
		get samtools single parameters
		:param parameter_name -> Only these two possibilities available Samtools ONT
		"""
		
		parameters_string = self.get_samtools_parameters_all_possibilities_ONT(
					project_sample.project.owner, project_sample)
		if (parameters_string is None): return None
		lst_data = parameters_string.split(parameter_name)
		if len(lst_data) == 2: return lst_data[1].split()[0]
		return None
	
	def is_samtools_single_parameter_default_for_project_ONT(self, project, parameter_name):
		"""
		test if a specific parameter is default Samtools ONT
		"""
		
		value_default_parameter = self.get_samtools_single_parameter_default_ONT(parameter_name)
		if (value_default_parameter is None): return False
		
		parameter_defined = self.get_samtools_single_parameter_for_project_ONT(project, parameter_name)
		if not parameter_defined is None and parameter_defined == value_default_parameter: return True
		return False
	
	def get_samtools_single_parameter_for_project_ONT(self, project, parameter_name):
		"""
		get samtools single parameters
		:param parameter_name -> Only these two possibilities available MEDAKA_model
		"""
		
		parameters_string = self.get_samtools_parameters_for_project_ONT(project.owner, project)
		if (parameters_string is None): return None
		lst_data = parameters_string.split(parameter_name)
		if len(lst_data) == 2: return lst_data[1].split()[0]
		return None

	#####
	#####		END samtools ONT
	#####
	#####################################################
	
	
	#####################################################
	#####
	#####		Mask consensus
	#####
	
	def get_mask_consensus_parameters(self, user, type_of_use, project, project_sample, technology_name):
		"""
		get mask_consensus parameters
		"""
		return self._get_parameters(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, user, type_of_use,
					project, project_sample, None, technology_name)
	
	def get_mask_consensus_parameters_all_possibilities(self, user, project_sample, technology_name):
		"""
		get mask_consensus parameters for project and default
		"""
		
		### Test project_sample first
		parameters = self._get_parameters(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, user,\
				Software.TYPE_OF_USE_project_sample, None, project_sample, None, technology_name)
		if (not parameters is None): return parameters
		
		### Test project
		parameters = self._get_parameters(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, user,\
				Software.TYPE_OF_USE_project, project_sample.project, None, None, technology_name)
		if (not parameters is None): return parameters
		
		### can be a default one
		default_software = DefaultSoftware()
		parameters = default_software.get_mask_consensus_threshold_parameters(user, technology_name)
		if (len(parameters) > 0): return parameters
		
		software_names = SoftwareNames()
		return software_names.get_insaflu_parameter_mask_consensus_parameters()
	
	
	def get_mask_consensus_parameters_for_project(self, user, project, technology_name):
		"""
		get mask_consensus parameters only for project or default
		"""
		
		### Test project
		parameters = self._get_parameters(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, user,\
				Software.TYPE_OF_USE_project, project, None, None, technology_name)
		if (not parameters is None): return parameters
		
		### can be a default one
		default_software = DefaultSoftware()
		parameters = default_software.get_mask_consensus_threshold_parameters(user, technology_name)
		if (len(parameters) > 0): return parameters
		return None
	
	def get_mask_consensus_single_parameter_default(self, parameter_name):
		"""
		:param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
		:return value of parameter
		"""
		vect_parameters = self._get_mask_consensus_threshold_default(None, None, None, None)
		for parameters in vect_parameters:
			if parameters.name == parameter_name:
				return parameters.parameter
		return None
	
	def is_mask_consensus_single_parameter_default(self, project_sample, parameter_name, technology_name):
		"""
		 one possibility available MASK_CONSENSUS_threshold
		"""
		
		value_default_parameter = self.get_mask_consensus_single_parameter_default(parameter_name)
		if (value_default_parameter is None): return False
		
		parameter_defined = self.get_mask_consensus_single_parameter(project_sample, parameter_name, technology_name)
		if not parameter_defined is None and parameter_defined == value_default_parameter: return True
		return False
	
	def get_mask_consensus_single_parameter(self, project_sample, parameter_name, technology_name):
		"""
		get mask_consensus single parameters
		:param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
		"""
		
		parameters_string = self.get_mask_consensus_parameters_all_possibilities(project_sample.project.owner,
								project_sample, technology_name)
		if (parameters_string is None): return None
		lst_data = parameters_string.split(parameter_name)
		if len(lst_data) == 2: return lst_data[1][1:]
		return None
	
	def is_mask_consensus_single_parameter_default_for_project(self, project, parameter_name, technology_name):
		"""
		:param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
		"""
		
		value_default_parameter = self.get_mask_consensus_single_parameter_default(parameter_name)
		if (value_default_parameter is None): return False
		
		parameter_defined = self.get_mask_consensus_single_parameter_for_project(project, parameter_name, technology_name)
		if not parameter_defined is None and parameter_defined == value_default_parameter: return True
		return False
	

	def get_mask_consensus_single_parameter_for_project(self, project, parameter_name, technology_name):
		"""
		get mask_consensus single parameters
		:param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
		"""
		
		parameters_string = self.get_mask_consensus_parameters_for_project(project.owner, project, technology_name)
		if (parameters_string is None): return None
		lst_data = parameters_string.split(parameter_name)
		if len(lst_data) == 2: return lst_data[1][1:]
		return None

	#####
	#####		END Mask consensus
	#####
	#####################################################

	
	
	#####################################################
	#####
	#####		Coverage limit ONT
	#####
	
	def get_limit_coverage_ONT_parameters(self, user, type_of_use, project, project_sample):
		"""
		get limit_coverage_ONT parameters
		"""
		return self._get_parameters(SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name, user, type_of_use,
					project, project_sample, None, SoftwareNames.TECHNOLOGY_minion)
	
	def get_limit_coverage_ONT_parameters_all_possibilities(self, user, project_sample):
		"""
		get limit_coverage parameters for project and default
		"""
		
		### Test project_sample first
		parameters = self._get_parameters(SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name, user,\
				Software.TYPE_OF_USE_project_sample, None, project_sample, None, SoftwareNames.TECHNOLOGY_minion)
		if (not parameters is None): return parameters
		
		### Test project
		parameters = self._get_parameters(SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name, user,\
				Software.TYPE_OF_USE_project, project_sample.project, None, None, SoftwareNames.TECHNOLOGY_minion)
		if (not parameters is None): return parameters
		
		### can be a default one
		default_software = DefaultSoftware()
		parameters = default_software.get_limit_coverage_ONT_parameters(user)
		if (len(parameters) > 0): return parameters
		
		software_names = SoftwareNames()
		return software_names.get_insaflu_parameter_limit_coverage_parameters()
	
	
	def get_limit_coverage_ONT_parameters_for_project(self, user, project):
		"""
		get limit_coverage_ONT parameters only for project or default
		"""
		
		### Test project
		parameters = self._get_parameters(SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name, user,\
				Software.TYPE_OF_USE_project, project, None, None, SoftwareNames.TECHNOLOGY_minion)
		if (not parameters is None): return parameters
		
		### can be a default one
		default_software = DefaultSoftware()
		parameters = default_software.get_limit_coverage_ONT_parameters(user)
		if (len(parameters) > 0): return parameters
		return None
	
	def get_limit_coverage_ONT_single_parameter_default(self, parameter_name):
		"""
		:param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
		:return value of parameter
		"""
		vect_parameters = self._get_limit_coverage_ONT_threshold_default(None, None, None, None)
		for parameters in vect_parameters:
			if parameters.name == parameter_name:
				return parameters.parameter
		return None
	
	def is_limit_coverage_ONT_single_parameter_default(self, project_sample, parameter_name):
		"""
		 one possibility available MASK_CONSENSUS_threshold
		"""
		
		value_default_parameter = self.get_limit_coverage_ONT_single_parameter_default(parameter_name)
		if (value_default_parameter is None): return False
		
		parameter_defined = self.get_limit_coverage_ONT_single_parameter(project_sample, parameter_name)
		if not parameter_defined is None and parameter_defined == value_default_parameter: return True
		return False
	
	def get_limit_coverage_ONT_single_parameter(self, project_sample, parameter_name):
		"""
		get limit_coverage_ONT single parameters
		:param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
		"""
		
		parameters_string = self.get_limit_coverage_ONT_parameters_all_possibilities(project_sample.project.owner,
								project_sample)
		if (parameters_string is None): return None
		lst_data = parameters_string.split(parameter_name)
		if len(lst_data) == 2: return lst_data[1][1:]
		return None
	
	def is_limit_coverage_ONT_single_parameter_default_for_project(self, project, parameter_name):
		"""
		:param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
		"""
		
		value_default_parameter = self.get_limit_coverage_ONT_single_parameter_default(parameter_name)
		if (value_default_parameter is None): return False
		
		parameter_defined = self.get_limit_coverage_ONT_single_parameter_for_project(project, parameter_name)
		if not parameter_defined is None and parameter_defined == value_default_parameter: return True
		return False
	

	def get_limit_coverage_ONT_single_parameter_for_project(self, project, parameter_name):
		"""
		get limit_coverage_ONT single parameters
		:param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
		"""
		
		parameters_string = self.get_limit_coverage_ONT_parameters_for_project(project.owner, project)
		if (parameters_string is None): return None
		lst_data = parameters_string.split(parameter_name)
		if len(lst_data) == 2: return lst_data[1][1:]
		return None
	
	#####
	#####		END Coverage limit ONT
	#####
	#####################################################
	
	
	def get_freebayes_parameters(self, user, type_of_use, project, project_sample):
		"""
		get freebayes parameters
		Add extra -V to the end
		"""
		return self._get_parameters(SoftwareNames.SOFTWARE_FREEBAYES_name, user, type_of_use, project,
				project_sample, None, SoftwareNames.TECHNOLOGY_illumina) + " -V"
	
	
	def _get_parameters(self, software_name, user, type_of_use, project, project_sample, sample, technology_name):
		"""
		:out return parameters for a specific software
		"""
		try:
			software = Software.objects.get(name=software_name, owner=user, type_of_use=type_of_use,
					technology__name = technology_name)
		except Software.DoesNotExist:
			return None

		## get parameters for a specific user
		parameters = Parameter.objects.filter(software=software, project=project,
						project_sample=project_sample, sample=sample)
		
		### parse them
		dict_out = {}
		vect_order_ouput = []
		for parameter in parameters:
			### don't set the not set parameters
			if (not parameter.not_set_value is None and parameter.parameter == parameter.not_set_value): continue
			
			### create a dict with parameters
			if (parameter.name in dict_out): 
				dict_out[parameter.name][1].append(parameter.parameter)
				dict_out[parameter.name][0].append(parameter.union_char)
			else:
				vect_order_ouput.append(parameter.name) 
				dict_out[parameter.name] = [[parameter.union_char], [parameter.parameter]]
			
		return_parameter = ""
		for par_name in vect_order_ouput:
			if (len(dict_out[par_name][1]) == 1 and len(dict_out[par_name][1][0]) == 0):
				return_parameter += " {}".format(par_name)
			else:
				return_parameter += " {}".format(par_name)
				for _ in range(len(dict_out[par_name][0])):
					return_parameter += "{}{}".format(dict_out[par_name][0][_], dict_out[par_name][1][_])
					
		#### This is the case where all the options can be "not to set"
		if (len(return_parameter.strip()) == 0 and len(parameters) == 0): return None
		return return_parameter.strip()

# 	def get_vect_parameters_higher_priority(self, software_name, user, type_of_use, project, project_sample):
# 		"""
# 		:out return parameters for a specific software
# 		"""
# 		try:
# 			software = Software.objects.get(name=software_name, owner=user, type_of_use=type_of_use)
# 		except Software.DoesNotExist:
# 			return ""
# 
# 		## get parameters for a specific user
# 		parameters = Parameter.objects.filter(software=software, project=project, project_sample=project_sample)
# 		
# 		### parse them
# 		dict_out = {}
# 		for parameter in parameters:
# 			### don't set the not set parameters
# 			if (not parameter.not_set_value is None and parameter.parameter == parameter.not_set_value): continue
# 			
# 			### create a dict with parameters
# 			if (parameter.name in dict_out): 
# 				dict_out[parameter.name][1].append(parameter.parameter)
# 				dict_out[parameter.name][0].append(parameter.union_char)
# 			else:
# 				dict_out[parameter.name] = [[parameter.union_char], [parameter.parameter]]
# 		return dict_out


	def set_default_software(self, software, user, type_of_use, project, project_sample, sample):
		"""
		Set default parameters for a software
		"""
		vect_parameters = self._get_default_parameters(software.name, user, type_of_use, project,
					project_sample, sample, software.technology.name)
		parameters = Parameter.objects.filter(software=software, project=project,
					project_sample=project_sample, sample=sample)
		key_value = "{}_{}".format(software.name, SoftwareNames.TECHNOLOGY_illumina if\
					software.technology is None else software.technology.name)
		self.change_values_software[key_value] = False
		for parameter in parameters:
			if parameter.can_change:
				for parameter_to_set_default in vect_parameters:
					if (parameter_to_set_default.sequence_out == parameter.sequence_out):
						###   if change software name
						if (parameter.parameter != parameter_to_set_default.parameter):
							self.change_values_software[key_value] = True
						
						parameter.parameter = parameter_to_set_default.parameter
						parameter.save()
						break

	def is_change_values_for_software(self, software_name, technology_name):
		""" Return if the software has a value changed"""
		key_value = "{}_{}".format(software_name, technology_name)
		return self.change_values_software.get(key_value, False)

	def _persist_parameters(self, vect_parameters, type_of_use, technology_name):
		"""
		presist a specific software by default
		param: type_of_use Can by Software.TYPE_OF_USE_project; Software.TYPE_OF_USE_project_sample
		"""
		default_software = DefaultSoftware()
		software = None
		dt_out_sequential = {}
		for parameter in vect_parameters:
			assert parameter.sequence_out not in dt_out_sequential
			if software is None:
				software = parameter.software
				software.technology = default_software.get_technology_instance(technology_name)
				try:
					software = Software.objects.get(name=parameter.software.name, owner=parameter.software.owner,\
						type_of_use=type_of_use, technology__name=technology_name)
				except Software.DoesNotExist:
					software.save()
			parameter.software = software
			parameter.save()
			
			## set sequential number
			dt_out_sequential[parameter.sequence_out] = 1

	def get_parameters(self, software_name, user, type_of_use, project, project_sample, sample,
				technology_name = SoftwareNames.TECHNOLOGY_illumina):
		"""
		"""
		self.test_default_db(software_name, user, type_of_use, project, project_sample, sample, technology_name)
		return self._get_parameters(software_name, user, type_of_use, project, project_sample, sample, technology_name)

	def get_all_software(self):
		"""
		get all softwares available by this class
		"""
		vect_software = []
		vect_software.append(self.software_names.get_snippy_name())
		vect_software.append(self.software_names.get_medaka_name_consensus())
		vect_software.append(self.software_names.get_samtools_name_depth_ONT())
		vect_software.append(self.software_names.get_NanoFilt_name())
		vect_software.append(self.software_names.get_insaflu_parameter_mask_consensus_name())
		vect_software.append(self.software_names.get_insaflu_parameter_limit_coverage_name())
#		vect_software.append(self.software_names.get_freebayes_name())
		return vect_software

	def _get_default_project(self, user, software_name, project, vect_parameters, technology_name):
		"""
		:param software_name name of the software
		:param project is None pass to global
		try to get project parameters
		"""
		type_of_use = Software.TYPE_OF_USE_global
		if project is None: type_of_use = Software.TYPE_OF_USE_global
		elif type(project) is Project: type_of_use = Software.TYPE_OF_USE_project
		elif type(project) is ProjectSample: type_of_use = Software.TYPE_OF_USE_project
		elif type(project) is Sample: type_of_use = Software.TYPE_OF_USE_sample
		
		try:
			software = Software.objects.get(name=software_name, owner=user,\
				type_of_use=type_of_use,
				technology__name = technology_name)
		except Software.DoesNotExist:
			return vect_parameters

		## get parameters for a specific user and software
		if project is None: parameters = Parameter.objects.filter(software=software)
		elif type(project) is Project: parameters = Parameter.objects.filter(software=software, project=project)
		elif type(project) is ProjectSample: parameters = Parameter.objects.filter(software=software, project_sample=project)
		elif type(project) is Sample: parameters = Parameter.objects.filter(software=software, sample=project)
		
		### parse them
		for parameter in parameters:
			### don't set the not set parameters
			if (not parameter.not_set_value is None and parameter.parameter == parameter.not_set_value): continue
			
			for previous_parameter in vect_parameters:
				if previous_parameter.sequence_out == parameter.sequence_out:
					previous_parameter.parameter = parameter.parameter
					break
		return vect_parameters

	def _get_snippy_default(self, user, type_of_use, project, project_sample):
		"""
		–mapqual: minimum mapping quality to allow (–mapqual 20)
		—mincov: minimum coverage of variant site (–mincov 10)
		–minfrac: minumum proportion for variant evidence (–minfrac 0.51)
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_SNIPPY_name
		software.name_extended = SoftwareNames.SOFTWARE_SNIPPY_name_extended
		software.version = SoftwareNames.SOFTWARE_SNIPPY_VERSION
		software.type_of_use = type_of_use
		software.type_of_software = Software.TYPE_SOFTWARE
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = DefaultProjectSoftware.SNIPPY_MAPQUAL_NAME
		parameter.parameter = "20"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[5:100]"
		parameter.range_max = "50"
		parameter.range_min = "10"
		parameter.description = "MAPQUAL: is the minimum mapping quality to accept in variant calling (–mapqual 20)."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = DefaultProjectSoftware.SNIPPY_COVERAGE_NAME
		parameter.parameter = "10"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 2
		parameter.range_available = "[4:100]"
		parameter.range_max = "100"
		parameter.range_min = "4"
		parameter.description = "MINCOV: the minimum number of reads covering a site to be considered (–mincov 10)."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--minfrac"
		parameter.parameter = "0.51"
		parameter.type_data = Parameter.PARAMETER_float
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 3
		parameter.range_available = "[0.5:1.0]"
		parameter.range_max = "1.0"
		parameter.range_min = "0.5"
		parameter.description = "MINFRAC: minumum proportion for variant evidence (–minfrac 0.51)"
		vect_parameters.append(parameter)
		
		return vect_parameters

	def _get_freebayes_default(self, user, type_of_use, project, project_sample):
		"""
		–min-mapping-quality: excludes read alignments from analysis if they have a mapping quality less than Q (–min-mapping-quality 20)
		—min-base-quality: excludes alleles from iSNV analysis if their supporting base quality is less than Q (–min-base-quality 20)
		–min-coverage: requires at least 100-fold of coverage to process a site (–min-coverage 100)
		—min-alternate-count: require at least 10 reads supporting an alternate allele within a single individual in order to evaluate the position (–min-alternate-count 10)
		–min-alternate-fraction: defines the minimum intra-host frequency of the alternate allele to be assumed (–min-alternate-fraction 0.01). This frequency is contingent on the depth of coverage of each processed site since min-alternate-count is set to 10, i.e., the identification of iSNV sites at frequencies of 10%, 2% and 1% is only allowed for sites with depth of coverage of at least 100-fold, 500-fold and 1000-fold, respectively.
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_FREEBAYES_name
		software.name_extended = SoftwareNames.SOFTWARE_FREEBAYES_name_extended
		software.version = SoftwareNames.SOFTWARE_FREEBAYES_VERSION
		software.type_of_use = type_of_use
		software.type_of_software = Software.TYPE_SOFTWARE
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "--min-mapping-quality"
		parameter.parameter = "20"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[10:50]"
		parameter.range_max = "50"
		parameter.range_min = "10"
		parameter.not_set_value = "0"
		parameter.description = "min-mapping-quality: excludes read alignments from analysis if they have a mapping quality less than Q (–min-mapping-quality 20)."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--min-base-quality"
		parameter.parameter = "20"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 2
		parameter.range_available = "[10:50]"
		parameter.range_max = "50"
		parameter.range_min = "10"
		parameter.not_set_value = "0"
		parameter.description = "min-base-quality: excludes alleles from iSNV analysis if their supporting base quality is less than Q (–min-base-quality 20)."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--min-coverage"
		parameter.parameter = "100"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 3
		parameter.range_available = "[20:500]"
		parameter.range_max = "500"
		parameter.range_min = "20"
		parameter.description = "min-coverage: requires at least 100-fold of coverage to process a site (–min-coverage 100)."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--min-alternate-count"
		parameter.parameter = "10"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 4
		parameter.range_available = "[5:100]"
		parameter.range_max = "100"
		parameter.range_min = "5"
		parameter.description = "min-alternate-count: require at least 10 reads supporting an alternate allele within a single " +\
			"individual in order to evaluate the position (–min-alternate-count 10)."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--min-alternate-fraction"
		parameter.parameter = "100"
		parameter.type_data = Parameter.PARAMETER_float
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 5
		parameter.range_available = "[0.01:0.5]"
		parameter.range_max = "0.5"
		parameter.range_min = "0.01"
		parameter.description = "min-alternate-fraction: defines the minimum intra-host frequency of the alternate allele to be assumed (–min-alternate-fraction 0.01)."
		vect_parameters.append(parameter)

		parameter = Parameter()
		parameter.name = "--ploidy"
		parameter.parameter = "2"
		parameter.type_data = Parameter.PARAMETER_float
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 6
		parameter.range_available = "[1:5]"
		parameter.range_max = "5"
		parameter.range_min = "1"
		parameter.description = "ploidy:  Sets the default ploidy for the analysis to N.  default: 2"
		vect_parameters.append(parameter)
		
		return vect_parameters

	def _get_mask_consensus_threshold_default(self, user, type_of_use, project, project_sample):
		"""
		Threshold of mask not consensus coverage
		"""
		software = Software()
		software.name = SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name
		software.name_extended = SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name_extended
		software.type_of_use = type_of_use
		software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
		software.version = "1.0"
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = DefaultProjectSoftware.MASK_CONSENSUS_threshold
		parameter.parameter = "70"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[50:100]"
		parameter.range_max = "100"
		parameter.range_min = "50"
		parameter.description = "Minimum percentage of locus horizontal coverage with depth of coverage equal or above –mincov (see Snippy) to generate consensus sequence."
		vect_parameters.append(parameter)
		return vect_parameters

	def _get_limit_coverage_ONT_threshold_default(self, user, type_of_use, project, project_sample):
		"""
		Minimum depth of coverage per site to validate the sequence (default: –mincov 30)
		Where to use this cut-off:
			This cut-off is used to exclude from vcf files sites with DEPTH <=30
			This cut-off is used to mask consensus sequences with DEPTH <=30 in msa_masker (-c: 30-1 = 29)
		"""
		software = Software()
		software.name = SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name
		software.name_extended = SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name_extended
		software.type_of_use = type_of_use
		software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
		software.version = "1.0"
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "Threshold"
		parameter.parameter = "30"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[4:100]"
		parameter.range_max = "100"
		parameter.range_min = "4"
		parameter.description = "This cut-off is used to exclude from vcf files sites and to mask consensus sequences. Value in percentage"
		vect_parameters.append(parameter)
		return vect_parameters

	def _get_medaka_model_default(self, user, type_of_use, project, project_sample):
		"""
		Minimum depth of coverage per site to validate the sequence (default: –mincov 30)
		Where to use this cut-off:
			This cut-off is used to exclude from vcf files sites with DEPTH <=30
			This cut-off is used to mask consensus sequences with DEPTH <=30 in msa_masker (-c: 30-1 = 29)
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_Medaka_name_consensus
		software.name_extended = SoftwareNames.SOFTWARE_Medaka_name_extended_consensus
		software.type_of_use = type_of_use
		software.type_of_software = Software.TYPE_SOFTWARE
		software.version = SoftwareNames.SOFTWARE_Medaka_VERSION
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "-m"
		parameter.parameter = SoftwareNames.SOFTWARE_Medaka_default_model	## default model
		parameter.type_data = Parameter.PARAMETER_char_list
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.description = "Medaka models are named to indicate: " +\
			"i) the pore type; " +\
			"ii) the sequencing device (min -> MinION, prom -> PromethION); " +\
			"iii) the basecaller variant (only high and variant available in INSAFlu); " +\
			"iv) the Guppy basecaller version. " +\
			"Complete format: " +\
			"{pore}_{device}_{caller variant}_{caller version}"
		vect_parameters.append(parameter)
		return vect_parameters

	def _get_nanofilt_default(self, user, type_of_use, sample):
		"""
		-l <LENGTH>, Filter on a minimum read length. Range: [50:1000].
		--maxlength Filter on a maximum read length
		dá para colocar outro parametro, por exemplo: -ml <MAX_LENGTH>, Set the maximum read length.
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_NanoFilt_name
		software.name_extended = SoftwareNames.SOFTWARE_NanoFilt_name_extended
		software.version = SoftwareNames.SOFTWARE_NanoFilt_VERSION
		software.type_of_use = type_of_use
		software.type_of_software = Software.TYPE_SOFTWARE
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "-q"
		parameter.parameter = "10"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[5:30]"
		parameter.range_max = "30"
		parameter.range_min = "5"
		parameter.description = "-q <QUALITY>, Filter on a minimum average read quality score."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "-l"
		parameter.parameter = "50"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 2
		parameter.range_available = "[50:1000]"
		parameter.range_max = "1000"
		parameter.range_min = "50"
		parameter.description = "-l <LENGTH>, Filter on a minimum read length."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--headcrop"
		parameter.parameter = "10"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 3
		parameter.range_available = "[1:1000]"
		parameter.range_max = "1000"
		parameter.range_min = "1"
		parameter.not_set_value = "0"
		parameter.description = "--headcrop <HEADCROP>, Trim n nucleotides from start of read."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--tailcrop"
		parameter.parameter = "10"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 4
		parameter.range_available = "[1:1000]"
		parameter.range_max = "1000"
		parameter.range_min = "1"
		parameter.not_set_value = "0"
		parameter.description = "--tailcrop <TAILCROP>, Trim n nucleotides from end of read."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--maxlength"
		parameter.parameter = "0"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 5
		parameter.range_available = "[{}:50000]".format(DefaultSoftware.NANOFILT_MINIMUN_MAX)
		parameter.range_max = "50000"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "--maxlength <LENGTH>, Set a maximum read length."
		vect_parameters.append(parameter)
		
		return vect_parameters

	def _get_samtools_depth_default_ONT(self, user, type_of_use, project, project_sample):
		"""
		samtools depth for ONT
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT
		software.name_extended = SoftwareNames.SOFTWARE_SAMTOOLS_name_extended_depth_ONT
		software.type_of_use = type_of_use
		software.type_of_software = Software.TYPE_SOFTWARE
		software.version = SoftwareNames.SOFTWARE_SAMTOOLS_VERSION
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "-q"
		parameter.parameter = "0"	## default model
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[0:100]"
		parameter.range_max = "100"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "-q <Quality> base quality threshold."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "-Q"
		parameter.parameter = "0"	## default model
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 2
		parameter.range_available = "[0:100]"
		parameter.range_max = "100"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "-Q <Quality> mapping quality threshold."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "-aa"
		parameter.parameter = ""	## default model
		parameter.type_data = Parameter.PARAMETER_null
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = ""
		parameter.can_change = False
		parameter.sequence_out = 3
		parameter.description = "Output absolutely all positions, including unused ref. sequences."
		vect_parameters.append(parameter)
		
		return vect_parameters