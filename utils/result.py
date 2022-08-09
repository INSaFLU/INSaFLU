'''
Created on Nov 2, 2017

@author: mmp
'''
import json, os
from constants.constants import Constants
from constants.constants_mixed_infection import ConstantsMixedInfection
from Bio.SeqFeature import SeqFeature, CompoundLocation, FeatureLocation

def is_integer(n_value):
	try:
		int(n_value)
		return True
	except ValueError: 
		return False


class DecodeObjects(object):

	def __init__(self):
		pass
	
	def decode_result(self, sz_temp):
		return json.loads(sz_temp, object_hook=self.decode_object)
		
	def decode_object(self, o):
		if '__Coverage__' in o:
			a = Coverage()
			a.__dict__.update(o['__Coverage__'])
			return a
		elif '__MaskingConsensus__' in o:
			a = MaskingConsensus()
			a.__dict__.update(o['__MaskingConsensus__'])
			return a
		elif '__CountHits__' in o:
			a = CountHits()
			a.__dict__.update(o['__CountHits__'])
			return a
		elif '__TasksToProcess__' in o:
			a = TasksToProcess()
			a.__dict__.update(o['__TasksToProcess__'])
			return a
		elif '__FeatureLocationSimple__' in o:
			a = FeatureLocationSimple(o['__FeatureLocationSimple__']['start'], o['__FeatureLocationSimple__']['end'],\
				o['__FeatureLocationSimple__']['strand'])
			a.__dict__.update(o['__FeatureLocationSimple__'])
			return a
		elif '__Gene__' in o:
			a = Gene(o['__Gene__']['name'], o['__Gene__']['start'], o['__Gene__']['end'],\
				o['__Gene__']['strand'])
			a.__dict__.update(o['__Gene__'])
			return a
		elif '__GeneticElement__' in o:
			a = GeneticElement()
			a.__dict__.update(o['__GeneticElement__'])
			return a
		elif '__CoverageElement__' in o:
			a = CoverageElement(o['__CoverageElement__']['element'])
			a.__dict__.update(o['__CoverageElement__'])
			return a
		elif '__Outputs__' in o:
			a = Outputs()
			a.__dict__.update(o['__Outputs__'])
			return a
		elif '__KeyValues__' in o:
			a = KeyValues()
			a.__dict__.update(o['__KeyValues__'])
			return a
		elif '__Softwares__' in o:
			a = Softwares()
			a.__dict__.update(o['__Softwares__'])
			return a
		elif '__Result__' in o:
			a = Result()
			a.__dict__.update(o['__Result__'])
			return a
		elif '__KeyValue__' in o:
			a = KeyValue(o['__KeyValue__']['key'], o['__KeyValue__']['value'])
			a.__dict__.update(o['__KeyValue__'])
			return a
		elif '__software_desc__' in o:
			a = SoftwareDesc(o['__software_desc__']['name'], o['__software_desc__']['version'], o['__software_desc__']['parameters'])
			a.__dict__.update(o['__software_desc__'])
			return a
		elif '__output__' in o:
			a = Output(o['__output__']['file_name'], o['__output__']['path'])
			a.__dict__.update(o['__output__'])
			return a
		elif '__MixedInfectionMainVector__' in o:
			a = MixedInfectionMainVector()
			a.__dict__.update(o['__MixedInfectionMainVector__'])
			return a
		elif '__SingleResult__' in o:
			a = SingleResult(o['__SingleResult__']['result'], o['__SingleResult__']['message'])
			a.__dict__.update(o['__SingleResult__'])
			return a
		elif '__ProcessResults__' in o:
			a = ProcessResults()
			a.__dict__.update(o['__ProcessResults__'])
			return a
		elif '__ResultAverageAndNumberReads__' in o:
			a = ResultAverageAndNumberReads(o['__ResultAverageAndNumberReads__']['number_file_1'],\
						o['__ResultAverageAndNumberReads__']['average_file_1'], o['__ResultAverageAndNumberReads__']['number_file_2'],\
						o['__ResultAverageAndNumberReads__']['average_file_2'], )
			a.__dict__.update(o['__ResultAverageAndNumberReads__'])
			return a
		return o

class ObjectEncoder(json.JSONEncoder):

	def default(self, o):
		return {'__{}__'.format(o.__class__.__name__): o.__dict__}

	
class SoftwareDesc(object):
	
	def __init__(self, name, version, parameters, key_values = None):
		"""
		:param key_values has the values from the analysis of the software
		"""
		self.name = name
		self.version = version
		self.parameters = parameters
		self.key_values = None
		if (not key_values is None):
			self.key_values = KeyValues(key_values)
		
	def __eq__(self, other):
		return other.name == self.name and other.version == self.version and other.parameters == self.parameters

	def __str__(self):
		return self.name

	def get_vect_key_values(self):
		return None if self.key_values is None else self.key_values.get_key_value() 


class Softwares(object):
	
	def __init__(self):
		self.list_software = []
	
	def add_software(self, software):
		self.list_software.append(software)

	def get_software(self, sz_name, b_not_add_software_name = False):
		"""   return software name, version and parameters"""
		list_return = []
		for software_desc in self.list_software:
			if (software_desc.name.lower().startswith(sz_name.lower())):
				version_data = ""
				## version
				if (len(software_desc.version) > 0): version_data = "-{}".format(software_desc.version)
				###
				if (not software_desc.parameters is None and len(software_desc.parameters) > 0):
					if b_not_add_software_name: description = "{}".format(software_desc.parameters)
					else: description = "{}{}; ({})".format(software_desc.name, version_data, software_desc.parameters)
					if (not description in list_return): list_return.append(description)
				else:
					description = "{}{}".format(software_desc.name, version_data)
					if (not description in list_return): list_return.append(description)
		if (len(list_return) > 0): return "/".join(list_return)
		return ""
	
	def is_software_present(self, sz_name):
		"""   test if a specific software is present in the list.
			This is useful to test if the software it runs or Not. 
		"""
		for software_desc in self.list_software:
			if (software_desc.name.startswith(sz_name)): return True
		return False
	
	def get_software_instance(self, sz_name):
		""" return software instance """
		for software_desc in self.list_software:
			if (software_desc.name == sz_name): return software_desc
		return None
	
	def get_software_version(self, sz_name):
		""" return software instance """
		for software_desc in self.list_software:
			if (software_desc.name == sz_name): return software_desc.version
		return None
	
	def get_list_software_instance(self, sz_name):
		""" return software instance """
		vect_software = []
		for software_desc in self.list_software:
			if (software_desc.name == sz_name): vect_software.append(software_desc)
		return vect_software
			
	def get_number_softwares(self):
		return len(self.list_software)

	def get_all_software_names(self):
		"""
		out: return all software names in an arrays
		"""
		vect_out = []
		for software in self.list_software:
			if (not software.name in vect_out): vect_out.append(software.name)
		return vect_out

class Output(object):
	def __init__(self, file_name, path):
		self.file_name = file_name
		self.path = path
		self.software = None

	def set_software(self, software):
		self.software = software

	def __eq__(self, other):
		return other.file_name == self.file_name and other.path == self.path
		
	def __str__(self):
		return self.file_name
	
class Outputs(object):
	def __init__(self):
		self.list_output = []
	
	def add_output(self, output):
		self.list_output.append(output)
		
class KeyValue(object):
	def __init__(self, key, value):
		self.key = key
		self.value = value
	
	def __eq__(self, other):
		return other.key == self.key and other.value == self.value
		
	def __str__(self):
		return "{} {}".format(self.key, self.value)

class KeyValues(object):
	def __init__(self, key_values = None):
		self.list_key_value = []
		
		## copy constructor
		if (not key_values is None):
			for key_value in key_values.get_key_value():
				self.list_key_value.append(KeyValue(key_value.key, key_value.value))
	
	def add_key_value(self, key_value):
		self.list_key_value.append(key_value)

	def get_key_value(self):
		return self.list_key_value
	
	def get_value_by_key(self, key):
		for data_ in self.list_key_value:
			if data_[0] == key: return data_[1]
		return None
	
class Result(object):
	'''
	classdocs
	'''
	SUCCESS = "Success"
	ERROR = "Error"
	def __init__(self):
		'''
		Constructor
		'''
		self.result = ""
		self.message = ""
		self.key_values = KeyValues()
		self.outputs = Outputs()
		self.softwares = Softwares()
	
	def set_error(self, message):
		self.result = self.ERROR
		self.message = message
	
	def set_success(self, message):
		self.result = self.SUCCESS
		self.message = message
	
	def add_output(self, output):
		self.outputs.add_output(output)
		
	def add_key_value(self, key_value):
		self.key_values.add_key_value(key_value)
		
	def add_software(self, software):
		self.softwares.add_software(software)
		
	def to_json(self):
		return json.dumps(self, indent=4, cls=ResultEncoder)

	def get_software(self, sz_name, b_not_add_software_name = False):
		return self.softwares.get_software(sz_name, b_not_add_software_name)
	
	def get_software_instance(self, sz_name):
		return self.softwares.get_software_instance(sz_name)
	
	def get_software_version(self, sz_name):
		return self.softwares.get_software_version(sz_name)
	
	def get_list_software_instance(self, sz_name):
		""" return a list of softwares """
		return self.softwares.get_list_software_instance(sz_name)
	
	def get_number_softwares(self):
		return self.softwares.get_number_softwares()
	
	def get_all_software_names(self):
		return self.softwares.get_all_software_names()
	
	def is_software_present(self, sz_name):
		return self.softwares.is_software_present(sz_name)
	
	def is_success(self):
		return self.result == Result.SUCCESS

	def get_key_value(self):
		return self.key_values.get_key_value()

	def get_value_by_key(self, key):
		return self.key_values.get_value_by_key(key)
	
class ResultEncoder(json.JSONEncoder):

	def default(self, o):
		if isinstance(o, SoftwareDesc):
			return {'__software_desc__': o.__dict__}
		elif isinstance(o, Output):
			return {'__output__': o.__dict__}
		return {'__{}__'.format(o.__class__.__name__): o.__dict__}
	

class ResultAverageAndNumberReads(object):
	"""
	Only have the number of reads and average
	"""
	def __init__(self, number_file_1, average_file_1, number_file_2, average_file_2):
		self.number_file_1 = None if number_file_1 is None else str(number_file_1)
		self.average_file_1 = None if average_file_1 is None else str(average_file_1)
		self.number_file_2 = None if number_file_2 is None else str(number_file_2)
		self.average_file_2 = None if average_file_2 is None else str(average_file_2)
	
	def to_json(self):
		return json.dumps(self, indent=4, cls=ObjectEncoder)

	def get_result_number(self):
		return ""		
	
	def __eq__(self, other):
		return other != None and other.number_file_1 == self.number_file_1 and other.average_file_1 == self.average_file_1 and\
			other.number_file_2 == self.number_file_2 and other.average_file_2 == self.average_file_2

	def has_reads(self):
		"""
		Test if has reads
		"""
		if (self.number_file_1 == '0' and
			(self.number_file_2 == '0' or self.number_file_2 is None)): return False
		return True

class CoverageElement(object):
	"""
	Only have the number of reads and average
	"""
	
	def __init__(self, element):
		self.element = element
		self.dt_data = {}
		
	def add_coverage(self, type_coverage, coverage):
		self.dt_data[type_coverage] = coverage
		
	def get_coverage(self, type_coverage):
		return self.dt_data.get(type_coverage, None)
	
	def __str__(self):
		return "Element: {}  {}".format(self.element, self.dt_data)


class Coverage(object):
	"""
	Only have the number of reads and average
	"""
	COVERAGE_ALL = "CoverageAll"
	COVERAGE_MORE_0 = "CoverageMore0"
	COVERAGE_MORE_9 = "CoverageMore9"
	COVERAGE_MORE_DEFINED_BY_USER = "CoverageMoreDefinedByUser"
	COVERAGE_PROJECT = "CoverageProject"

	def __init__(self, limit_defined_by_user = None, limit_defined_to_project = None):
		self.limit_defined_by_user = limit_defined_by_user			### if there is limit defined by user put it here
																	## limit_defined_by_user is called first
		self.limit_defined_to_project = limit_defined_to_project	### if there is limit defined to the project
		self.dt_data = {}
		self.ratio_value_defined_by_user = 0
		self.ratio_value_project = 0
		self.ratio_value_9 = 0
		self.ratio_value_0 = 0

	def is_exist_limit_defined_by_user(self):
		""" test if limit defined by user exists """
		return not self.limit_defined_by_user is None
	
	def is_exist_limit_defined_to_project(self):
		""" test if limit defined by user exists """
		return not self.limit_defined_to_project is None
	
	def get_middle_limit(self):
		"""
		:return threshold of middle limit 
		"""
		if (not self.limit_defined_by_user is None): return self.limit_defined_by_user
		if (not self.limit_defined_to_project is None): return self.limit_defined_to_project
		return 10
	
	def get_type_coverage_middle_limit(self):
		"""
		:return threshold of middle limit 
		"""
		if (not self.limit_defined_by_user is None): return self.COVERAGE_MORE_DEFINED_BY_USER
		if (not self.limit_defined_to_project is None): return self.COVERAGE_PROJECT
		return self.COVERAGE_MORE_9

	def get_dict_data(self): return self.dt_data
	
	def get_sorted_elements_name(self): return sorted(self.dt_data.keys())
	
	def add_coverage(self, element, type_coverage, coverage):
		if (element in self.dt_data): self.dt_data[element].add_coverage(type_coverage, coverage)
		else:
			self.dt_data[element] = CoverageElement(element)
			self.dt_data[element].add_coverage(type_coverage, coverage)
	
	def get_coverage_by_middle_tag(self, element):
		type_coverage = self.get_type_coverage_middle_limit()
		if (element in self.dt_data): return self.dt_data[element].get_coverage(type_coverage)
		raise Exception("Error: there's no key like this: " + element)
	
	def get_coverage(self, element, type_coverage):
		if (element in self.dt_data): return self.dt_data[element].get_coverage(type_coverage)
		raise Exception("Error: there's no key like this: " + element)

	def exist_this_element(self, element):
		return (element in self.dt_data)

	def to_json(self):
		return json.dumps(self, indent=4, cls=ObjectEncoder)

	def get_result_number(self):
		return ""
	
	def ratio_value_coverage_bigger_limit(self, element, limit_to_mask_consensus):
		"""
		return the True if ratio bigger than limit_to_mask_consensus
		"""
		### set ratio value variable
		self.is_100_more_9(element)
		
		if (self.ratio_value_project != 0 and self.ratio_value_project > limit_to_mask_consensus): return True
		if (self.ratio_value_defined_by_user != 0 and self.ratio_value_defined_by_user > limit_to_mask_consensus): return True
		if (self.ratio_value_9 != 0 and self.ratio_value_9 > limit_to_mask_consensus): return True
		return False

	def is_100_total(self, element):
		value_coverage = self.get_coverage(element, self.COVERAGE_ALL)
		return self.__is_100__(value_coverage)
	
	def is_100_more_9(self, element, b_only_project = False):
		"""
		:param b_only_project -> Test only project project limits 
		"""
		if (b_only_project):		### test only for a value in a project
			#### if exist value defined by the user call next method
			if (self.is_exist_limit_defined_to_project()):
				return self.is_100_more_defined_to_project(element)
		#### if exist value defined by the user call next method
		elif (self.is_exist_limit_defined_by_user()):
			return self.is_100_more_defined_by_user(element)
		
		### regular method
		value_coverage = self.get_coverage(element, self.COVERAGE_MORE_9)
		self.ratio_value_9 = divmod(float(value_coverage), 1)[0]
		return self.__is_100__(value_coverage)

	def is_100_more_defined_by_user(self, element):
		value_coverage = self.get_coverage(element, self.COVERAGE_MORE_DEFINED_BY_USER)
		if (value_coverage == None): return False
		self.ratio_value_defined_by_user = divmod(float(value_coverage), 1)[0]
		return self.__is_100__(value_coverage)
	
	def is_100_more_defined_to_project(self, element):
		"""
		this is only used for testing coverage when this value is defined to the a project
		"""
		value_coverage = self.get_coverage(element, self.COVERAGE_PROJECT)
		if (value_coverage == None): return False
		self.ratio_value_project = divmod(float(value_coverage), 1)[0]
		return self.__is_100__(value_coverage)
	
	def is_100_more_0(self, element):
		value_coverage = self.get_coverage(element, self.COVERAGE_MORE_0)
		self.ratio_value_0 = divmod(float(value_coverage), 1)[0]
		return self.__is_100__(value_coverage)
		
	def __is_100__(self, value_coverage):
		if (value_coverage == None): return False
		try:
			(i, d) = divmod(float(value_coverage), 1)
			if (int(i) == 100): return True
		except ValueError: 
			pass
		return False

	def __str__(self):
		sz_return = ""
		for key in self.dt_data:
			sz_return += "{} - All {}; 9- {}; 0- {}\n".format(key, self.get_coverage(key, Coverage.COVERAGE_ALL),\
				self.get_coverage(key, Coverage.COVERAGE_MORE_9), self.get_coverage(key, Coverage.COVERAGE_MORE_0))
		return sz_return

	def get_fault_message_9(self, element_name):
		return "Fail, locus '{}': the % of locus size covered by at least 10-fold is '{}%' (below 100%)".format(element_name, self.ratio_value_9)
	def get_fault_message_defined_by_user(self, element_name, limit_defined_by_user):
		return "Fail, locus '{}': the % of locus size covered by at least {}-fold is '{}%' (below 100%)".format(\
				element_name, limit_defined_by_user, self.ratio_value_defined_by_user)
	def get_fault_message_0(self, element_name):
		return "Fail, locus '{}': the % of locus size covered by at least 1-fold is '{}%' (below 100%)".format(element_name, self.ratio_value_0)

	def get_message_to_show_in_web_site(self, sample_name, element):
		"""
		get message for web site about coverage in More than 9
		"""
		if (not self.limit_defined_by_user is None):
			return "Sample name: {}\nLocus: {}\n\nMean depth of coverage: {}\n% of size covered by at least 1-fold: {}%\n% of size covered by at least {}-fold: {}%".format(
					sample_name, element,
					self.get_coverage(element, self.COVERAGE_ALL), self.get_coverage(element, self.COVERAGE_MORE_0),\
					self.limit_defined_by_user,\
					self.get_coverage(element, self.COVERAGE_MORE_DEFINED_BY_USER))
		elif (not self.limit_defined_to_project is None):
			return "Sample name: {}\nLocus: {}\n\nMean depth of coverage: {}\n% of size covered by at least 1-fold: {}%\n% of size covered by at least {}-fold: {}%".format(
					sample_name, element,
					self.get_coverage(element, self.COVERAGE_ALL), self.get_coverage(element, self.COVERAGE_MORE_0),\
					self.limit_defined_to_project,\
					self.get_coverage(element, self.COVERAGE_PROJECT))
		else:
			return "Sample name: {}\nLocus: {}\n\nMean depth of coverage: {}\n% of size covered by at least 1-fold: {}%\n% of size covered by at least 10-fold: {}%".format(
					sample_name, element,
					self.get_coverage(element, self.COVERAGE_ALL), self.get_coverage(element, self.COVERAGE_MORE_0),\
					self.get_coverage(element, self.COVERAGE_MORE_9))
	
	def get_icon(self, element, limit_to_mask_consensus):
		
		"""
		get coverage for COVERAGE_MORE_9
		:param limit_to_mask_consensus can be -1 or integer
		"""
		### GREEN
		if (not self.limit_defined_by_user is None):	### not defined by user
			if (self.is_100_more_defined_by_user(element)): return os.path.join("/" + Constants.DIR_STATIC, Constants.DIR_ICONS, Constants.ICON_GREEN_16_16)
		elif (not self.limit_defined_to_project is None):
			if (self.is_100_more_defined_to_project(element)): return os.path.join("/" + Constants.DIR_STATIC, Constants.DIR_ICONS, Constants.ICON_GREEN_16_16)
		elif (self.is_100_more_9(element)): return os.path.join("/" + Constants.DIR_STATIC, Constants.DIR_ICONS, Constants.ICON_GREEN_16_16)

		if (limit_to_mask_consensus > 0 and self.ratio_value_coverage_bigger_limit(element, limit_to_mask_consensus)):
			return os.path.join("/" + Constants.DIR_STATIC, Constants.DIR_ICONS, Constants.ICON_YELLOW_16_16)
		return os.path.join("/" + Constants.DIR_STATIC, Constants.DIR_ICONS, Constants.ICON_RED_16_16)

	def get_color(self, element, limit_to_mask_consensus):
		
		"""
		get coverage for COVERAGE_MORE_9
		:param limit_to_mask_consensus can be -1 or integer
		:out return three color "red", "yellow", "green"
		"""
		### GREEN
		if (not self.limit_defined_by_user is None):	### not defined by user
			if (self.is_100_more_defined_by_user(element)): return "green"
		elif (not self.limit_defined_to_project is None):
			if (self.is_100_more_defined_to_project(element)): return "green"
		elif (self.is_100_more_9(element)): return "green"

		if (limit_to_mask_consensus > 0 and self.ratio_value_coverage_bigger_limit(element, limit_to_mask_consensus)):
			return "yellow"
		return "red"
	
class CountHits(object):
	"""
	Count the hits in the variations
	"""
	TOTAL_GRATHER_THAN_MIXED = 100	## used to test mixed infection
	
	def __init__(self):
		self.hits_more_90 = 0
		self.hits_50_90 = 0
		self.hits_less_50 = 0
		
	def set_hits_50_90(self, hits_50_90):
		self.hits_50_90 = hits_50_90
		
	def set_hits_more_90(self, hits_more_90):
		self.hits_more_90 = hits_more_90
		
	def set_hits_less_50(self, hits_less_50):
		self.hits_less_50 = hits_less_50
	
	def add_one_hits_less_50(self):
		self.hits_less_50 += 1
	
	def add_one_hits_50_90(self):
		self.hits_50_90 += 1
	
	def add_one_hits_more_90(self):
		self.hits_more_90 += 1
	
	def get_hits_50_90(self):
		return self.hits_50_90
	
	def get_hits_more_90(self):
		return self.hits_more_90
	
	def get_hits_less_50(self):
		return self.hits_less_50
	
	def get_total_50_50_90(self):
		return self.hits_50_90 + self.hits_less_50
	
	def get_mixed_infection_ratio(self):
		if (self.hits_less_50 == 0): return 0.0;
		return self.hits_50_90 / float(self.hits_less_50)
	
	def get_mixed_infection_ratio_str(self):
		return '%.1f' % (self.get_mixed_infection_ratio())
	
	def get_total(self):
		return self.hits_50_90 + self.hits_less_50 + self.hits_more_90

	def get_vect_mixed_infections(self):
		return [self.get_hits_less_50(), self.get_hits_50_90()]

	def is_mixed_infection_ratio_test(self):
		"""
		mixed infection positive
		0.5 < ratio < 1.5 and total > 20
		"""
		if (self.get_mixed_infection_ratio() > 0.5 and self.get_mixed_infection_ratio() < 2.0 and self.get_total_50_50_90() > 20): return True
		return False

	def total_grather_than_mixed_infection(self):
		"""
		if its bigger than 200 is also mixed infection
		"""
		if (self.get_total_50_50_90() > CountHits.TOTAL_GRATHER_THAN_MIXED): return True
		return False
	
	def to_json(self):
		return json.dumps(self, indent=4, cls=ObjectEncoder)

	def __eq__(self, other):
		return other != None and other.hits_50_90 == self.hits_50_90 and other.hits_more_90 == self.hits_more_90\
			and other.hits_less_50 == self.hits_less_50

	def __str__(self):
		return "hits_more_90: {}  hits_50_90: {}   hits_less_50: {}".format(self.hits_more_90, self.hits_50_90, self.hits_less_50)


class TasksToProcess(object):
	"""
	task to process for a project
	"""
	def __init__(self):
		self.vect_tasks_id = []

	def add_taskd_id(self, task_id):
		if (task_id not in self.vect_tasks_id): self.vect_tasks_id.append(task_id)

	def get_tasks_id(self):
		return self.vect_tasks_id

	def to_json(self):
		return json.dumps(self, indent=4, cls=ObjectEncoder)
	
	def __eq__(self, other):
		if (other == None or len(other.get_tasks_id()) != len(self.vect_tasks_id)): return False
		for value_ in self.vect_tasks_id:
			if (value_ not in other.vect_tasks_id): return False
		return True
	
class FeatureLocationSimple(object):
	
	def __init__(self, start, end, strand):
		self.start = start
		self.end = end
		self.strand = strand	## 1 forward, -1 reverse
		
	def is_forward(self):
		return self.strand > 0
		
	def __eq__(self, other):
		if (other == None): return False 
		return self.start == other.start and self.end == other.end and self.strand == other.strand
		
	def __str__(self):
		return "[{}-{}] ({})".format(self.start, self.end, self.strand)

	def get_feature_location(self):
		return FeatureLocation(self.start, self.end, strand=self.strand)

class Gene(object):
	
	def __init__(self, name, start, end, strand, vect_feature_locations = []):
		self.name = name
		self.start = start
		self.end = end
		self.strand = strand	## 1 forward, -1 reverse
		self.vect_feature_locations = vect_feature_locations
		
	def is_forward(self):
		return self.strand == 1
		
	def __eq__(self, other):
		if (other == None): return False 
		return self.name == other.name and self.start == other.start and \
			self.end == other.end and self.strand == other.strand and \
			self.vect_feature_locations == other.vect_feature_locations
		
	def __str__(self):
		return "{} {}/{} {}".format(self.name, self.start, self.end, self.strand)

	def get_feature_locations(self):
		""" return feature locations"""
		return self.vect_feature_locations

	def has_feature_locations(self):
		""" return feature locations"""
		return len(self.vect_feature_locations) > 0
	
	def get_seq_feature(self):
		if (self.has_feature_locations()):
			vect_compound_location = []
			for feature_location in self.vect_feature_locations:
				vect_compound_location.append(feature_location.get_feature_location())
			try:	## some files doesn't have this well formated
				return SeqFeature(CompoundLocation(vect_compound_location), type="CDS")
			except:
				pass
		return SeqFeature(FeatureLocation(self.start, self.end, strand=self.strand), type="CDS")
		
class GeneticElement(object):
	
	def __init__(self):
		self.name = ""
		self.dt_elements = {}
		self.dt_elements_size = {}
		self.dt_elements_mask = {}
		
	def add_gene(self, element_name, length, gene):
		if (element_name in self.dt_elements):
			if (gene == None): return False
			if (gene not in self.dt_elements[element_name]): self.dt_elements[element_name].append(gene)
			else: return False
		else: 
			if (gene == None): self.dt_elements[element_name] = []
			else: self.dt_elements[element_name] = [gene]
			self.dt_elements_size[element_name] = length
		return True

	def get_genes(self, element_name):
		if (element_name in self.dt_elements): return self.dt_elements[element_name]
		return None
	
	def get_gene(self, element_name, gene_name):
		if (element_name in self.dt_elements):
			for gene in self.dt_elements[element_name]:
				if (gene.name == gene_name): return gene
		return None
	
	def has_genes(self, element_name):
		if (element_name in self.dt_elements): return len(self.dt_elements[element_name]) > 0
		return None
	
	def get_size_element(self, element_name):
		if (element_name in self.dt_elements_size): return self.dt_elements_size[element_name]
		return None
	
	def get_mask_consensus_element(self, element_name):
		if (element_name in self.dt_elements_mask): return self.dt_elements_mask[element_name]
		return None
	
	def set_mask_consensus_element(self, element_name, mask_consensus):
		""" MaskingConsensus """
		self.dt_elements_mask[element_name] = mask_consensus
	
	def get_vect_gene_names(self, element_name):
		if (element_name in self.dt_elements): 
			return [gene.name for gene in self.dt_elements[element_name]]
		return None
	
	def get_sorted_elements(self):
		return sorted(self.dt_elements.keys())
	
	def has_masking_data(self):
		""" testing if has masking data """
		for element in self.dt_elements_mask:
			if self.dt_elements_mask[element].has_data(): return True
		return False
	
	def get_message_mask_to_show_in_web_site(self):
		sz_return = ""
		for element in self.get_sorted_elements():
			if (len(sz_return) > 0): sz_return += "\n#####################\n"
			else: sz_return += "#####################\n"
			sz_return += self.dt_elements_mask[element].get_message_to_show_in_web_site(element)
		return sz_return

	def get_message_to_show_in_csv_file(self):
		sz_return = ""
		for element in self.get_sorted_elements():
			if not element in self.dt_elements_mask: continue
			sz_out = self.dt_elements_mask[element].get_message_to_show_in_csv_file(element)
			if (len(sz_out) > 0):
				if (len(sz_return) > 0): sz_return += " "
				sz_return += sz_out
		return sz_return
	
	def cleaning_mask_results(self):
		""" cleaning masking values """			
		for element in self.get_sorted_elements():
			self.dt_elements_mask[element].cleaning_mask_results()
			
	def to_json(self):
		return json.dumps(self, indent=4, cls=ObjectEncoder)

	def __eq__(self, other):
		if (other == None or len(other.dt_elements) != len(self.dt_elements)): return False
		for value_ in self.dt_elements:
			if (value_ not in other.dt_elements or other.dt_elements[value_] != self.dt_elements[value_]): return False
		for value_ in self.dt_elements_size:
			if (value_ not in other.dt_elements_size or other.dt_elements_size[value_] != self.dt_elements_size[value_]): return False
		for value_ in self.dt_elements_mask:
			if (value_ not in other.dt_elements_mask or other.dt_elements_mask[value_] != self.dt_elements_mask[value_]): return False
		return True
		
	def __str__(self):
		sz_out = ""
		for element_name in self.dt_elements:
			if (len(self.dt_elements[element_name]) == 0): sz_out += "-- E.Name: {} ".format(element_name)
			else: sz_out += ";; Name: {}; Gene: {}; ".format(element_name, "" if len(self.dt_elements[element_name]) == 0 \
									else "-".join([str(gene) for gene in self.dt_elements[element_name]]))
		return sz_out

class MixedInfectionMainVector(object):
	
	def __init__(self):
		'''
		Constructor
		'''
		self.vect_start_compare = [[b for b in a] for a in ConstantsMixedInfection.vect_start_compare]
	
	def add_vector(self, vector):
		if (vector not in self.vect_start_compare): self.vect_start_compare.append(vector)
		
	def to_json(self):
		return json.dumps(self, indent=4, cls=ObjectEncoder)
	
	def get_vector(self):
		return self.vect_start_compare

	def __eq__(self, other):
		if (other == None): return False
		return self.vect_start_compare == other.vect_start_compare
	
	def __str__(self):
		return '; '.join(['-'.join([str(b) for b in a]) for a in self.vect_start_compare])


class SingleResult(object):
	'''
	single status result
	'''
	SUCCESS = "Success"
	ERROR = "Error"
	def __init__(self, result, message):
		'''
		Constructor
		'''
		self.result = result
		self.message = message

	def set_error(self, message):
		self.result = self.ERROR
		self.message = message
	
	def set_success(self, message):
		self.result = self.SUCCESS
		self.message = message
		
	def is_success(self):
		return self.result == Result.SUCCESS
	
	def __eq__(self, other):
		if (other == None): return False
		return self.result == other.result and self.message == other.message
	
	def __str__(self):
		return '{} - {}'.format(self.result, self.message)


class ProcessResults(object):
	"""
	global process result
	"""
	def __init__(self):
		self.vect_results = []
		
	def add_single_result(self, single_result):
		self.vect_results.append(single_result)

	def get_vect_results(self):
		return self.vect_results

	def get_len_vect_results(self):
		return len(self.vect_results)
	
	def get_error(self, index):
		if (index < len(self.vect_results)): return self.vect_results[index]
		return None
	
	def has_errors(self):
		return self.get_len_vect_results() > 0
	
	def to_json(self):
		return json.dumps(self, indent=4, cls=ObjectEncoder)

	def __eq__(self, other):
		if (other == None): return False
		return self.vect_results == other.vect_results

	def __str__(self):
		return '\n'.join([str(a) for a in self.vect_results])


class MaskingConsensus(object):
	"""
	Mask Consensus sequences...
	"""
	def __init__(self):
		self.mask_sites = "" 			### <number>,<number>,...
		self.mask_from_beginning = ""	### <number>
		self.mask_from_ends = ""		### <number>
		self.mask_regions = ""			### [<number>-<number>],[<number>-<number>],...
		
	def set_mask_sites(self, mask_sites):
		lst_data = mask_sites.split(',')
		vect_data = []
		for data_ in lst_data:
			if (is_integer(data_) and not abs(int(data_)) in vect_data): vect_data.append(abs(int(data_)))
		vect_data = sorted(vect_data)
		self.mask_sites = ",".join([str(_) for _ in vect_data])
	
	def set_mask_from_beginning(self, mask_from_beginning):
		if (is_integer(mask_from_beginning)):
			self.mask_from_beginning = "{}".format(abs(int(mask_from_beginning)))
		else: self.mask_from_beginning = ""
	
	def set_mask_from_ends(self, mask_from_ends):
		if (is_integer(mask_from_ends)):
			self.mask_from_ends = "{}".format(abs(int(mask_from_ends)))
		else: self.mask_from_ends = ""
	
	def set_mask_regions(self, mask_regions):
		self._clean_mask_regions(mask_regions)
	
	def _clean_mask_regions(self, mask_regions):
		lst_data = mask_regions.split(',')
		vect_data = []
		for data_ in lst_data:
			lst_positions = data_.strip().replace('[','').replace(']','').split('-')
			if len(lst_positions) == 2 and is_integer(lst_positions[0]) and \
				is_integer(lst_positions[1]):
				pos_1 = abs(int(lst_positions[0]))
				pos_2 = abs(int(lst_positions[1]))
				if pos_1 > pos_2: vect_data.append([pos_2, pos_1])
				elif pos_1 < pos_2: vect_data.append([pos_1, pos_2])
		
		## merge data
		if len(vect_data) > 0:
			vect_data = self._merge(vect_data)
			self.mask_regions = ",".join(["{}-{}".format(data_[0], data_[1]) for data_ in vect_data])
		else: self.mask_regions = ""
	
	def cleaning_mask_results(self):
		""" cleaning all results """
		self.set_mask_sites(self.mask_sites)
		self.set_mask_from_beginning(self.mask_from_beginning)
		self.set_mask_from_ends(self.mask_from_ends)
		self._clean_mask_regions(self.mask_regions)
		
	# if the two intervals overlaps
	def _is_overlaping(self, a, b):
		""" is overlapping """
		if b[0] >= a[0] and b[0] <= a[1]: return True
		else: return False

	# merge the intervals
	def _merge(self, vect_data):
		""" merge intervals """
		#sort the intervals by its first value
		vect_data = sorted(vect_data, key = lambda x: x[0])
		merged_list= []
		merged_list.append(vect_data[0])
		for i in range(1, len(vect_data)):
			pop_element = merged_list.pop()
			if self._is_overlaping(pop_element, vect_data[i]):
				new_element = pop_element[0], max(pop_element[1], vect_data[i][1])
				merged_list.append(new_element)
			else:
				merged_list.append(pop_element)
				merged_list.append(vect_data[i])
		return merged_list
	
	def get_header(self, separator = ","):
		""" return header """	
		return "Mask sites{}Mask from beginning{}Mask from end{}Mask regions".format(separator,
				separator, separator)
	def get_vect_header(self):
		""" return header """	
		return ["Mask sites", "Mask from beginning", "Mask from end", "Mask regions"]
		
	def get_mask_sites(self):
		return self.mask_sites
	
	def get_mask_from_beginning(self):
		return self.mask_from_beginning
	
	def get_mask_from_ends(self):
		return self.mask_from_ends
	
	def get_mask_regions(self):
		return self.mask_regions

	def get_message_to_show_in_web_site(self, element_name):
		"""  get message for web site about regions to be masked  """
		if not self.has_data(): return "Element:{} -> No mask".format(element_name)
		return "Element:{}\n\nMask sites:{}  Mask regions:{}\nMask from beginning:{}  Mask from ends:{}".format(
			element_name, self.mask_sites, self.mask_regions,
			self.mask_from_beginning, self.mask_from_ends)
		
	def get_message_to_show_in_csv_file(self, element_name):
		"""  get message for web site about regions to be masked  """
		if not self.has_data(): return ""
		sz_out = "Element:{} -> ".format(element_name)
		b_add = False
		if len(self.mask_sites) > 0:
			sz_out += "{}Mask sites:{}".format(" " if b_add else "(", self.mask_sites)
			b_add = True
		if len(self.mask_regions) > 0:
			sz_out += "{}Mask regions:{}".format(" " if b_add else "(", self.mask_regions)
			b_add = True
		if len(self.mask_from_beginning) > 0:
			sz_out += "{}Mask from beginning:{}".format(" " if b_add else "(", self.mask_from_beginning)
			b_add = True
		if len(self.mask_from_ends) > 0:
			sz_out += "{}Mask from end:{}".format(" " if b_add else "(", self.mask_from_ends)
		return sz_out + ")"
	
	def has_data(self):
		return len(self.mask_sites) > 0 or len(self.mask_from_beginning) > 0 or \
			len(self.mask_from_ends) > 0 or len(self.mask_regions) > 0
		
	def to_json(self):
		return json.dumps(self, indent=4, cls=ObjectEncoder)

	def __eq__(self, other):
		return other != None and other.mask_regions == self.mask_regions and other.mask_from_ends == self.mask_from_ends\
			and other.mask_from_beginning == self.mask_from_beginning and other.mask_sites == self.mask_sites

	def __str__(self):
		return "mask_sites: {}  mask_from_beginning: {}   mask_from_ends: {}".format(self.mask_sites, self.mask_from_beginning, self.mask_from_ends) + \
			"  mask_regions: {}".format(self.mask_regions)



