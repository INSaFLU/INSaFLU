'''
Created on 03/05/2020

@author: mmp
'''
from settings.models import Software, Parameter
from constants.software_names import SoftwareNames
from utils.lock_atomic_transaction import LockedAtomicTransaction

class DefaultSoftware(object):
	'''
	classdocs
	'''
	software_names = SoftwareNames()
	
	def __init__(self):
		""" change values """
		self.change_values_software = {}	### the key is the name of the software

	def test_all_defaults(self, user):
		### test all defaults
		self.test_default_trimmomatic_db(user)
		self.test_default_snippy_db(user)

	def test_default_trimmomatic_db(self, user):
		"""
		test if exist, if not persist in database
		"""
		try:
			Software.objects.get(name=SoftwareNames.SOFTWARE_TRIMMOMATIC_name, owner=user,\
						type_of_use = Software.TYPE_OF_USE_global)
		except Software.DoesNotExist:
			### create a default one for this user
			with LockedAtomicTransaction(Software), LockedAtomicTransaction(Parameter):
				vect_parameters = self._get_trimmomatic_default(user)
				self._persist_parameters(vect_parameters)

	def test_default_snippy_db(self, user):
		"""
		test if exist, if not persist in database
		"""
		try:
			Software.objects.get(name=SoftwareNames.SOFTWARE_SNIPPY_name, owner=user,\
						type_of_use = Software.TYPE_OF_USE_global)
		except Software.DoesNotExist:
			### create a default one for this user
			with LockedAtomicTransaction(Software), LockedAtomicTransaction(Parameter):
				vect_parameters = self._get_snippy_default(user)
				self._persist_parameters(vect_parameters)

	def get_trimmomatic_parameters(self, user):
		return self._get_parameters(user, SoftwareNames.SOFTWARE_TRIMMOMATIC_name)
	def get_snippy_parameters(self, user):
		return self._get_parameters(user, SoftwareNames.SOFTWARE_SNIPPY_name)
	
	def _get_parameters(self, user, software_name):
		"""
		get trimmomatic parameters
		"""
		try:
			software = Software.objects.get(name=software_name, owner=user,\
						type_of_use=Software.TYPE_OF_USE_global)
		except Software.DoesNotExist:
			return ""

		## get parameters for a specific user
		parameters = Parameter.objects.filter(software=software)
		
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

	def set_default_software(self, software, user):
		if (software.name == SoftwareNames.SOFTWARE_SNIPPY_name): vect_parameters = self._get_snippy_default(user)
		elif (software.name == SoftwareNames.SOFTWARE_TRIMMOMATIC_name): vect_parameters = self._get_trimmomatic_default(user)
		else: return
		
		parameters = Parameter.objects.filter(software=software)
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

	def is_change_values_for_software(self, software):
		""" Return if the software has a value changed"""
		return self.change_values_software.get(software.name, False)

	def _persist_parameters(self, vect_parameters):
		"""
		presist a specific software by default 
		"""
		software = None
		dt_out_sequential = {}
		for parameter in vect_parameters:
			assert parameter.sequence_out not in dt_out_sequential
			if software is None:
				software = parameter.software 
				software.save()
			parameter.software = software
			parameter.save()
			
			## set sequential number
			dt_out_sequential[parameter.sequence_out] = 1

	def get_parameters(self, software_name, user):
		"""
		"""
		if (software_name == SoftwareNames.SOFTWARE_TRIMMOMATIC_name):
			self.test_default_trimmomatic_db(user) 
			return self.get_trimmomatic_parameters(user)
		if (software_name == SoftwareNames.SOFTWARE_SNIPPY_name):
			self.test_default_snippy_db(user) 
			return self.get_snippy_parameters(user)
		return ""


	def get_all_software(self):
		"""
		get all softwares available by this class
		"""
		vect_software = []
		vect_software.append(self.software_names.get_trimmomatic_name())
		vect_software.append(self.software_names.get_snippy_name())
		return vect_software

	def _get_snippy_default(self, user):
		"""
		–mapqual: minimum mapping quality to allow (–mapqual 20)
		—mincov: minimum coverage of variant site (–mincov 10)
		–minfrac: minumum proportion for variant evidence (–minfrac 0.51)
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_SNIPPY_name
		software.version = SoftwareNames.SOFTWARE_SNIPPY_VERSION
		software.type_of_use = Software.TYPE_OF_USE_global
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "--mapqual"
		parameter.parameter = "20"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[5:100]"
		parameter.range_max = "50"
		parameter.range_min = "10"
		parameter.description = "MAPQUAL: is the minimum mapping quality to accept in variant calling (–mapqual 20)."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--mincov"
		parameter.parameter = "10"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
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
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 3
		parameter.range_available = "[0.5:1.0]"
		parameter.range_max = "1.0"
		parameter.range_min = "0.5"
		parameter.description = "MINFRAC: minumum proportion for variant evidence (–minfrac 0.51)"
		vect_parameters.append(parameter)
		
		return vect_parameters


	def _get_trimmomatic_default(self, user):
		"""
		SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_TRIMMOMATIC_name
		software.version = SoftwareNames.SOFTWARE_TRIMMOMATIC_VERSION
		software.type_of_use = Software.TYPE_OF_USE_global
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "HEADCROP"
		parameter.parameter = "0"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[0:100]"
		parameter.range_max = "100"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "HEADCROP:<length> Cut the specified number of bases from the start of the read."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "CROP"
		parameter.parameter = "0"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 2
		parameter.range_available = "[0:400]"
		parameter.range_max = "400"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "CROP:<length> Cut the read to a specified length."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "SLIDINGWINDOW"
		parameter.parameter = "5"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 3
		parameter.range_available = "[3:50]"
		parameter.range_max = "50"
		parameter.range_min = "3"
		parameter.description = "SLIDINGWINDOW:<windowSize> specifies the number of bases to average across"
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "SLIDINGWINDOW"
		parameter.parameter = "20"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 4
		parameter.range_available = "[10:100]"
		parameter.range_max = "100"
		parameter.range_min = "10"
		parameter.description = "SLIDINGWINDOW:<requiredQuality> specifies the average quality required"
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "LEADING"
		parameter.parameter = "3"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 5
		parameter.range_available = "[0:100]"
		parameter.range_max = "100"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "LEADING:<quality> Remove low quality bases from the beginning."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "TRAILING"
		parameter.parameter = "3"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 6
		parameter.range_available = "[0:100]"
		parameter.range_max = "100"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "TRAILING:<quality> Remove low quality bases from the end."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "MINLEN"
		parameter.parameter = "35"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 7
		parameter.range_available = "[5:500]"
		parameter.range_max = "500"
		parameter.range_min = "5"
		parameter.description = "SMINLEN:<length> This module removes reads that fall below the specified minimal length."
		vect_parameters.append(parameter)

##		Only available in 0.30 version		
#
# 		parameter = Parameter()
# 		parameter.name = "AVGQUAL"
# 		parameter.parameter = "0"
# 		parameter.type_data = Parameter.PARAMETER_int
# 		parameter.software = software
# 		parameter.union_char = ":"
# 		parameter.can_change = True
# 		parameter.sequence_out = 8
# 		parameter.range_available = "[0:100]"
# 		parameter.range_max = "100"
# 		parameter.range_min = "0"
# 		parameter.not_set_value = "0"
# 		parameter.description = "AVGQUAL:<quality> Drop the read if the average quality is below the specified level."
# 		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "TOPHRED33"
		parameter.parameter = ""
		parameter.type_data = Parameter.PARAMETER_null
		parameter.software = software
		parameter.union_char = ""
		parameter.can_change = False
		parameter.sequence_out = 8
		parameter.description = "This (re)encodes the quality part of the FASTQ file to base 33."
		vect_parameters.append(parameter)
		return vect_parameters






