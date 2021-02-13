'''
Created on Nov 2, 2017

@author: mmp
'''
import json, os
from constants.constants import Constants
from constants.constants_mixed_infection import ConstantsMixedInfection

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
		elif '__CountHits__' in o:
			a = CountHits()
			a.__dict__.update(o['__CountHits__'])
			return a
		elif '__TasksToProcess__' in o:
			a = TasksToProcess()
			a.__dict__.update(o['__TasksToProcess__'])
			return a
		elif '__Gene__' in o:
			a = Gene(o['__Gene__']['name'], o['__Gene__']['start'], o['__Gene__']['end'], o['__Gene__']['strand'])
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

	def get_software(self, sz_name):
		"""   return software name, version and parameters"""
		list_return = []
		for software_desc in self.list_software:
			if (software_desc.name.startswith(sz_name)):
				version_data = ""
				if (len(software_desc.version) > 0): version_data = "-{}".format(software_desc.version)
				if (not software_desc.parameters is None and len(software_desc.parameters) > 0):
					description = "{}{}; ({})".format(software_desc.name, version_data, software_desc.parameters)
					if (not description in list_return): list_return.append(description)
				else:
					description = "{}{}".format(software_desc.name, version_data)
					if (not description in list_return): list_return.append(description)
		if (len(list_return) > 0): return "/".join(list_return)
		return ""
	
	def get_software_instance(self, sz_name):
		""" return software instance """
		for software_desc in self.list_software:
			if (software_desc.name == sz_name): return software_desc
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

	def get_software(self, sz_name):
		return self.softwares.get_software(sz_name)
	
	def get_software_instance(self, sz_name):
		return self.softwares.get_software_instance(sz_name)
	
	def get_list_software_instance(self, sz_name):
		""" return a list of softwares """
		return self.softwares.get_list_software_instance(sz_name)
	
	def get_number_softwares(self):
		return self.softwares.get_number_softwares()
	
	def get_all_software_names(self):
		return self.softwares.get_all_software_names()
	
	def is_success(self):
		return self.result == Result.SUCCESS

	def get_key_value(self):
		return self.key_values.get_key_value()

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
		self.number_file_1 = number_file_1
		self.average_file_1 = average_file_1
		self.number_file_2 = number_file_2
		self.average_file_2 = average_file_2
	
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
		if ((self.number_file_1 == 0 and self.average_file_1 == 0) and
			((self.number_file_2 == 0 and self.average_file_2 == 0) or 
			(self.number_file_2 is None and self.average_file_2 is None))): return False
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
		self.limit_defined_by_user = limit_defined_by_user		### if there is limit defined by user put it here
		self.limit_defined_to_project = limit_defined_to_project		### if there is limit defined to the project
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

	def get_dict_data(self): return self.dt_data
	
	def get_sorted_elements_name(self): return sorted(self.dt_data.keys())
	
	def add_coverage(self, element, type_coverage, coverage):
		if (element in self.dt_data): self.dt_data[element].add_coverage(type_coverage, coverage)
		else:
			self.dt_data[element] = CoverageElement(element)
			self.dt_data[element].add_coverage(type_coverage, coverage)

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
		else:#### if exist value defined by the user call next method
			if (self.is_exist_limit_defined_by_user()):
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

	def get_message_to_show_in_web_site(self, element):
		"""
		get message for web site about coverage in More than 9
		"""
		if (not self.limit_defined_by_user is None):
			return "Locus: {}\n\nMean depth of coverage: {}\n% of size covered by at least 1-fold: {}%\n% of size covered by at least {}-fold: {}%".format(element,\
					self.get_coverage(element, self.COVERAGE_ALL), self.get_coverage(element, self.COVERAGE_MORE_0),\
					self.limit_defined_by_user,\
					self.get_coverage(element, self.COVERAGE_MORE_DEFINED_BY_USER))
		elif (not self.limit_defined_to_project is None):
			return "Locus: {}\n\nMean depth of coverage: {}\n% of size covered by at least 1-fold: {}%\n% of size covered by at least {}-fold: {}%".format(element,\
					self.get_coverage(element, self.COVERAGE_ALL), self.get_coverage(element, self.COVERAGE_MORE_0),\
					self.limit_defined_to_project,\
					self.get_coverage(element, self.COVERAGE_PROJECT))
		else:
			return "Locus: {}\n\nMean depth of coverage: {}\n% of size covered by at least 1-fold: {}%\n% of size covered by at least 10-fold: {}%".format(element,\
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

class Gene(object):
	
	def __init__(self, name, start, end, strand):
		self.name = name
		self.start = start
		self.end = end
		self.strand = strand	## 1 forward, -1 reverse
		
	def is_forward(self):
		return self.strand == 1
		
	def __eq__(self, other):
		if (other == None): return False 
		return self.name == other.name
		
	def __str__(self):
		return "{} {}/{} {}".format(self.name, self.start, self.end, self.strand)

class GeneticElement(object):
	
	def __init__(self):
		self.name = ""
		self.dt_elements = {}
		self.dt_elements_size = {}
		
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
	
	def has_genes(self, element_name):
		if (element_name in self.dt_elements): return len(self.dt_elements[element_name]) > 0
		return None
	
	def get_size_element(self, element_name):
		if (element_name in self.dt_elements_size): return self.dt_elements_size[element_name]
		return None
	
	def get_vect_gene_names(self, element_name):
		if (element_name in self.dt_elements): 
			return [gene.name for gene in self.dt_elements[element_name]]
		return None
	
	def get_sorted_elements(self):
		return sorted(self.dt_elements.keys())
	
	def to_json(self):
		return json.dumps(self, indent=4, cls=ObjectEncoder)

	def __eq__(self, other):
		if (other == None or len(other.dt_elements) != len(self.dt_elements)): return False
		for value_ in self.dt_elements:
			if (value_ not in other.dt_elements or other.dt_elements[value_] != self.dt_elements[value_]): return False
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


