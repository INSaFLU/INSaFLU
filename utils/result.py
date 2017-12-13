'''
Created on Nov 2, 2017

@author: mmp
'''
import json, os
from constants.constants import Constants

class SoftwareDesc(object):
	
	def __init__(self, name, version, parameters):
		self.name = name
		self.version = version
		self.parameters = parameters
		
	def __eq__(self, other):
		return other.name == self.name and other.version == self.version and other.parameters == self.parameters

	def __str__(self):
		return self.name


class Softwares(object):
	
	def __init__(self):
		self.list_software = []
	
	def add_software(self, software):
		self.list_software.append(software)

	def get_software(self, sz_name):
		for software_desc in self.list_software:
			if (software_desc.name == sz_name):
				if (software_desc.parameters != None and len(software_desc.parameters) > 0):
					return "{}-{}; ({})".format(software_desc.name, software_desc.version, software_desc.parameters)
				return "{}-{}".format(software_desc.name, software_desc.version)
		return ""
	
	def get_number_softwares(self):
		return len(self.list_software)

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
		self.outputs = Outputs()
		self.softwares = Softwares()
	
	def set_error(self, message):
		self.result = self.SUCCESS
		self.message = message
	
	def set_success(self, message):
		self.result = self.ERROR
		self.message = message
	
	def add_output(self, output):
		self.outputs.add_output(output)
		
	def add_software(self, software):
		self.softwares.add_software(software)
		
	def to_json(self):
		return json.dumps(self, indent=4, cls=ResultEncoder)

	def get_software(self, sz_name):
		return self.softwares.get_software(sz_name)
	
	def get_number_softwares(self):
		return self.softwares.get_number_softwares()
	
class ResultEncoder(json.JSONEncoder):

	def default(self, o):
		if isinstance(o, SoftwareDesc):
			return {'__software_desc__': o.__dict__}
		elif isinstance(o, Output):
			return {'__output__': o.__dict__}
		return {'__{}__'.format(o.__class__.__name__): o.__dict__}
	
class DecodeResult(object):
	
	def __init__(self):
		pass
	
	def decode_result(self, sz_temp):
		return json.loads(sz_temp, object_hook=self.decode_object)
	
	def decode_object(self, o):
		
		if '__Outputs__' in o:
			a = Outputs()
			a.__dict__.update(o['__Outputs__'])
			return a
		elif '__Softwares__' in o:
			a = Softwares()
			a.__dict__.update(o['__Softwares__'])
			return a
		elif '__Result__' in o:
			a = Result()
			a.__dict__.update(o['__Result__'])
			return a
		elif '__software_desc__' in o:
			a = SoftwareDesc(o['__software_desc__']['name'], o['__software_desc__']['version'], o['__software_desc__']['parameters'])
			a.__dict__.update(o['__software_desc__'])
			return a
		elif '__output__' in o:
			a = Output(o['__output__']['file_name'], o['__output__']['path'])
			a.__dict__.update(o['__output__'])
			return a
		
		return o


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
		return json.dumps(self, indent=4, cls=ResultAverageEncoder)

	def get_result_number(self):
		return ""		
	
	def __eq__(self, other):
		return other != None and other.number_file_1 == self.number_file_1 and other.average_file_1 == self.average_file_1 and\
			other.number_file_2 == self.number_file_2 and other.average_file_2 == self.average_file_2


class ResultAverageEncoder(json.JSONEncoder):

	def default(self, o):
		return {'__{}__'.format(o.__class__.__name__): o.__dict__}


class DecodeResultAverageAndNumberReads(object):
	
	def __init__(self):
		pass
	
	def decode_result(self, sz_temp):
		return json.loads(sz_temp, object_hook=self.decode_object)
		
	def decode_object(self, o):
		if '__ResultAverageAndNumberReads__' in o:
			a = ResultAverageAndNumberReads(o['__ResultAverageAndNumberReads__']['number_file_1'],\
						o['__ResultAverageAndNumberReads__']['average_file_1'], o['__ResultAverageAndNumberReads__']['number_file_2'],\
						o['__ResultAverageAndNumberReads__']['average_file_2'], )
			a.__dict__.update(o['__ResultAverageAndNumberReads__'])
			return a
		return o


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
		if (type_coverage in self.dt_data): return self.dt_data[type_coverage]
		raise Exception("Error: there's no key like this: " + type_coverage)


class Coverage(object):
	"""
	Only have the number of reads and average
	"""
	COVERAGE_ALL = "CoverageAll"
	COVERAGE_MORE_0 = "CoverageMore0"
	COVERAGE_MORE_9 = "CoverageMore9"

	def __init__(self):
		self.dt_data = {}

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

	def to_json(self):
		return json.dumps(self, indent=4, cls=CoverageEncoder)

	def get_result_number(self):
		return ""
	
	def is_100_total(self, element):
		value_coverage = self.get_coverage(element, self.COVERAGE_ALL)
		return self.__is_100__(value_coverage)
	
	def is_100_more_9(self, element):
		value_coverage = self.get_coverage(element, self.COVERAGE_MORE_9)
		return self.__is_100__(value_coverage)
	
	def is_100_more_0(self, element):
		value_coverage = self.get_coverage(element, self.COVERAGE_MORE_0)
		return self.__is_100__(value_coverage)
		
	def __is_100__(self, value_coverage):
		try:
			(i, d) = divmod(float(value_coverage), 1)
			if (int(i) == 100): return True
		except ValueError: 
			pass
		return False

	def __str__(self):
		sz_return = ""
		for key in self.dt_data:
			sz_return += "{} - All {}; 0 {}; 9 {}\n".format(key, self.get_coverage(key, Coverage.COVERAGE_ALL),\
				self.get_coverage(key, Coverage.COVERAGE_MORE_0), self.get_coverage(key, Coverage.COVERAGE_MORE_9))
		return sz_return

	def get_fault_message_9(self, element_name):
		return "Fail, element '{}' has a ratio under of 100 for a coverage of 10 or bigger".format(element_name)
	def get_fault_message_0(self, element_name):
		return "Fail, element '{}' has a ratio under of 100 for a coverage of 1 or bigger".format(element_name)

	def get_message_to_show_in_web_site(self, element):
		"""
		get message for web site about coverage in More than 9
		"""
		return "Element: {}\n\nCoverage: {}\nRatio coverage >9: {}%\nRatio coverage >0: {}%".format(element,\
				self.get_coverage(element, self.COVERAGE_ALL), self.get_coverage(element, self.COVERAGE_MORE_9),\
				self.get_coverage(element, self.COVERAGE_MORE_0))
	
	def get_icon(self, element):
		"""
		get coverage for COVERAGE_MORE_9
		"""
		if (self.is_100_more_9(element)): return os.path.join("/" + Constants.DIR_STATIC, Constants.DIR_ICONS, Constants.ICON_GREEN_16_16)
		if (self.is_100_more_0(element)): return os.path.join("/" + Constants.DIR_STATIC, Constants.DIR_ICONS, Constants.ICON_YELLOW_16_16)
		return os.path.join("/" + Constants.DIR_STATIC, Constants.DIR_ICONS, Constants.ICON_RED_16_16)

class CoverageEncoder(json.JSONEncoder):

	def default(self, o):
		return {'__{}__'.format(o.__class__.__name__): o.__dict__}


class DecodeCoverage(object):
	
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
		return o


class CountHits(object):
	"""
	Count the hits in the variations
	"""
	def __init__(self):
		self.hits_50_90 = 0
		self.hits_less_50 = 0
		
	def set_hits_50_90(self, hits_50_90):
		self.hits_50_90 = hits_50_90
		
	def set_hits_less_50(self, hits_less_50):
		self.hits_less_50 = hits_less_50
	
	def add_one_hits_less_50(self):
		self.hits_less_50 += 1
	
	def add_one_hits_50_90(self):
		self.hits_50_90 += 1
	
	def get_hits_50_90(self):
		return self.hits_50_90
	
	def get_hits_less_50(self):
		return self.hits_less_50
	
	def get_total(self):
		return self.hits_50_90 + self.hits_less_50

	def to_json(self):
		return json.dumps(self, indent=4, cls=CoverageEncoder)

	def __eq__(self, other):
		return other != None and other.hits_50_90 == self.hits_50_90 and other.hits_less_50 == self.hits_less_50

	def __str__(self):
		return "hits_50_90: {}   hits_less_50: {}".format(self.hits_50_90, self.hits_less_50)


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
		return json.dumps(self, indent=4, cls=CoverageEncoder)
	
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
		return self.name == other.name
		
class GeneticElement(object):
	
	def __init__(self):
		self.name = ""
		self.dt_elements = {}
		
	def add_gene(self, element_name, gene):
		if (element_name in self.dt_elements): 
			if (gene not in self.dt_elements[element_name]): self.dt_elements[element_name].append(gene)
			else: return False
		else: self.dt_elements[element_name] = [gene]
		return True

	def get_genes(self, element_name):
		if (element_name in self.dt_elements): return self.dt_elements[element_name]
		return None
	
	def get_sorted_elements(self):
		return sorted(self.dt_elements.keys())
	
	def to_json(self):
		return json.dumps(self, indent=4, cls=CoverageEncoder)

	def __eq__(self, other):
		if (other == None or len(other.dt_elements) != len(self.dt_elements)): return False
		for value_ in self.dt_elements:
			if (value_ not in other.dt_elements or other.dt_elements[value_] != self.dt_elements[value_]): return False
		return True
		
		