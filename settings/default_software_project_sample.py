'''
Created on 03/05/2020

@author: mmp
'''
from settings.models import Software, Parameter
from constants.software_names import SoftwareNames
from utils.lock_atomic_transaction import LockedAtomicTransaction
from django.core.exceptions import MultipleObjectsReturned
from settings.default_software import DefaultSoftware

class DefaultProjectSoftware(object):
	'''
	classdocs
	'''
	software_names = SoftwareNames()

	### used in snippy
	SNIPPY_COVERAGE_NAME = "--mincov"
	SNIPPY_MAPQUAL_NAME = "--mapqual"

	### used in mask consensus
	MASK_CONSENSUS_threshold = "Threshold"

	def __init__(self):
		""" change values """
		self.change_values_software = {}	### the key is the name of the software
		
	def test_all_defaults(self, user, type_of_use, project, project_sample):
		"""
		test all defaults for all software available
		"""
		self.test_default_db(SoftwareNames.SOFTWARE_SNIPPY_name, user, type_of_use, project, project_sample)
## Not in used yet
##		self.test_default_db(SoftwareNames.SOFTWARE_FREEBAYES_name, user, project)

		
		self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,\
						user, type_of_use, project, project_sample)

	def test_default_db(self, software_name, user, type_of_use, project, project_sample):
		"""
		test if exist, if not persist in database
		"""
		try:
			Software.objects.get(name=software_name, owner=user, type_of_use=type_of_use,
					parameter__project=project, parameter__project_sample=project_sample)
		except Software.DoesNotExist:
			### create a default one for this user
			with LockedAtomicTransaction(Software), LockedAtomicTransaction(Parameter):
				vect_parameters = self._get_default_parameters(software_name, user, type_of_use, project, project_sample)
				### persist 
				if (len(vect_parameters) > 0):
					self._persist_parameters(vect_parameters, type_of_use)
		except MultipleObjectsReturned:
			pass

	def _get_default_parameters(self, software_name, user, type_of_use, project, project_sample):
		if (software_name == SoftwareNames.SOFTWARE_SNIPPY_name):
			vect_parameters = self._get_snippy_default(user, type_of_use, project, project_sample)		### base values
			if (not project is None): vect_parameters = self._get_default_project(user,\
					SoftwareNames.SOFTWARE_SNIPPY_name, None, vect_parameters)		### base values
			if (not project_sample is None): vect_parameters = self._get_default_project(user,\
					SoftwareNames.SOFTWARE_SNIPPY_name, project_sample.project, vect_parameters)		### base values
			return vect_parameters
		elif (software_name == SoftwareNames.SOFTWARE_FREEBAYES_name):
			return self._get_freebayes_default(user, type_of_use, project, project_sample)
		elif (software_name == SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name):
			vect_parameters = self._get_mask_consensus_threshold_default(user, type_of_use, project, project_sample)
			if (not project is None): vect_parameters = self._get_default_project(user,\
				SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, None, vect_parameters)		### base values
			if (not project_sample is None): vect_parameters = self._get_default_project(user,\
				SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, project_sample.project, vect_parameters)		### base values
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
		return self._get_parameters(SoftwareNames.SOFTWARE_SNIPPY_name, user, type_of_use, project, project_sample)
	
	def get_snippy_parameters_all_possibilities(self, user, project_sample):
		"""
		get snippy parameters for project_sample, project and default
		"""
		
		### Test project_sample first
		parameters = self._get_parameters(SoftwareNames.SOFTWARE_SNIPPY_name, user,\
				Software.TYPE_OF_USE_project_sample, None, project_sample)
		if (len(parameters) > 0): return parameters
		
		### Test project
		parameters = self._get_parameters(SoftwareNames.SOFTWARE_SNIPPY_name, user,\
				Software.TYPE_OF_USE_project, project_sample.project, None)
		if (len(parameters) > 0): return parameters
		
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
				Software.TYPE_OF_USE_project, project, None)
		if (len(parameters) > 0): return parameters
		
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
	#####		Mask consensus
	#####
	
	def get_mask_consensus_parameters(self, user, type_of_use, project, project_sample):
		"""
		get snippy parameters
		"""
		return self._get_parameters(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, user, type_of_use, project, project_sample)
	
	def get_mask_consensus_parameters_all_possibilities(self, user, project_sample):
		"""
		get mask_consensus parameters for project and default
		"""
		
		### Test project_sample first
		parameters = self._get_parameters(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, user,\
				Software.TYPE_OF_USE_project_sample, None, project_sample)
		if (len(parameters) > 0): return parameters
		
		### Test project
		parameters = self._get_parameters(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, user,\
				Software.TYPE_OF_USE_project, project_sample.project, None)
		if (len(parameters) > 0): return parameters
		
		### can be a default one
		default_software = DefaultSoftware()
		parameters = default_software.get_mask_consensus_threshold_parameters(user)
		if (len(parameters) > 0): return parameters
		
		software_names = SoftwareNames()
		return software_names.get_insaflu_parameter_mask_consensus_parameters()
	
	
	def get_mask_consensus_parameters_for_project(self, user, project):
		"""
		get snippy parameters only for project or default
		"""
		
		### Test project
		parameters = self._get_parameters(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, user,\
				Software.TYPE_OF_USE_project, project, None)
		if (len(parameters) > 0): return parameters
		
		### can be a default one
		default_software = DefaultSoftware()
		parameters = default_software.get_mask_consensus_threshold_parameters(user)
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

	def is_mask_consensus_single_parameter_default(self, project_sample, parameter_name):
		"""
		 one possibility available MASK_CONSENSUS_threshold
		"""
		
		value_default_parameter = self.get_mask_consensus_single_parameter_default(parameter_name)
		if (value_default_parameter is None): return False
		
		parameter_defined = self.get_mask_consensus_single_parameter(project_sample, parameter_name)
		if not parameter_defined is None and parameter_defined == value_default_parameter: return True
		return False
		
	def get_mask_consensus_single_parameter(self, project_sample, parameter_name):
		"""
		get snippy single parameters
		:param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
		"""
		
		parameters_string = self.get_mask_consensus_parameters_all_possibilities(project_sample.project.owner, project_sample)
		if (parameters_string is None): return None
		lst_data = parameters_string.split(parameter_name)
		if len(lst_data) == 2: return lst_data[1][1:]
		return None
	
	def is_mask_consensus_single_parameter_default_for_project(self, project, parameter_name):
		"""
		:param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
		"""
		
		value_default_parameter = self.get_mask_consensus_single_parameter_default(parameter_name)
		if (value_default_parameter is None): return False
		
		parameter_defined = self.get_mask_consensus_single_parameter_for_project(project, parameter_name)
		if not parameter_defined is None and parameter_defined == value_default_parameter: return True
		return False
	

	def get_mask_consensus_single_parameter_for_project(self, project, parameter_name):
		"""
		get snippy single parameters
		:param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
		"""
		
		parameters_string = self.get_mask_consensus_parameters_for_project(project.owner, project)
		if (parameters_string is None): return None
		lst_data = parameters_string.split(parameter_name)
		if len(lst_data) == 2: return lst_data[1][1:]
		return None


	#####
	#####		END Mask consensus
	#####
	#####################################################
	
	def get_freebayes_parameters(self, user, type_of_use, project, project_sample):
		"""
		get freebayes parameters
		Add extra -V to the end
		"""
		return self._get_parameters(SoftwareNames.SOFTWARE_FREEBAYES_name, user, type_of_use, project, project_sample) + " -V"
	
	
	def _get_parameters(self, software_name, user, type_of_use, project, project_sample):
		"""
		:out return parameters for a specific software
		"""
		try:
			software = Software.objects.get(name=software_name, owner=user, type_of_use=type_of_use)
		except Software.DoesNotExist:
			return ""

		## get parameters for a specific user
		parameters = Parameter.objects.filter(software=software, project=project, project_sample=project_sample)
		
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


	def set_default_software(self, software, user, type_of_use, project, project_sample):
		"""
		Set default software
		"""
		vect_parameters = self._get_default_parameters(software.name, user, type_of_use, project, project_sample)
		parameters = Parameter.objects.filter(software=software, project=project, project_sample=project_sample)
		self.change_values_software[software.name] = False
		for parameter in parameters:
			if parameter.can_change:
				for parameter_to_set_default in vect_parameters:
					if (parameter_to_set_default.sequence_out == parameter.sequence_out):
						###   if change software name
						if (parameter.parameter != parameter_to_set_default.parameter):
							self.change_values_software[software.name] = True
						
						parameter.parameter = parameter_to_set_default.parameter
						parameter.save()
						break

	def is_change_values_for_software(self, software_name):
		""" Return if the software has a value changed"""
		return self.change_values_software.get(software_name, False)

	def _persist_parameters(self, vect_parameters, type_of_use):
		"""
		presist a specific software by default
		param: type_of_use Can by Software.TYPE_OF_USE_project; Software.TYPE_OF_USE_project_sample
		"""
		software = None
		dt_out_sequential = {}
		for parameter in vect_parameters:
			assert parameter.sequence_out not in dt_out_sequential
			if software is None:
				software = parameter.software
				try:
					software = Software.objects.get(name=parameter.software.name, owner=parameter.software.owner,\
						type_of_use=type_of_use)
				except Software.DoesNotExist:
					software.save()
			parameter.software = software
			parameter.save()
			
			## set sequential number
			dt_out_sequential[parameter.sequence_out] = 1

	def get_parameters(self, software_name, user, type_of_use, project, project_sample):
		"""
		"""
		self.test_default_db(software_name, user, type_of_use, project, project_sample)
		return self._get_parameters(software_name, user, type_of_use, project, project_sample)

	def get_all_software(self):
		"""
		get all softwares available by this class
		"""
		vect_software = []
		vect_software.append(self.software_names.get_snippy_name())
#		vect_software.append(self.software_names.get_freebayes_name())
		return vect_software

	def _get_default_project(self, user, software_name, project, vect_parameters):
		"""
		:param software_name name of the software
		:param project is None pass to global
		try to get project parameters
		"""
		try:
			software = Software.objects.get(name=software_name, owner=user,\
				type_of_use=Software.TYPE_OF_USE_global if project is None else Software.TYPE_OF_USE_project)
		except Software.DoesNotExist:
			return vect_parameters

		## get parameters for a specific user and software
		if project is None: parameters = Parameter.objects.filter(software=software)
		else: parameters = Parameter.objects.filter(software=software, project=project)
		
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


	def _get_mask_consensus_default_project(self, user, project, vect_parameters):
		"""
		try to get project parameters
		"""
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, owner=user,\
						type_of_use=Software.TYPE_OF_USE_project)
		except Software.DoesNotExist:
			return vect_parameters

		## get parameters for a specific user and software
		parameters = Parameter.objects.filter(software=software, project=project)
		
		### parse them
		for parameter in parameters:
			### don't set the not set parameters
			if (not parameter.not_set_value is None and parameter.parameter == parameter.not_set_value): continue
			
			for previous_parameter in vect_parameters:
				if previous_parameter.sequence_out == parameter.sequence_out:
					previous_parameter.parameter = parameter.parameter
					break
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

