'''
Created on Nov 2, 2017

@author: mmp
'''
import json
import humanfriendly

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
		return json.dumps(self, indent=4, cls=CustomEncoder)

class CustomEncoder(json.JSONEncoder):

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
		return json.dumps(self, indent=4, cls=CustomEncoderResult)

	def get_result_number(self):
		return ""		
	
	def __eq__(self, other):
		return other.number_file_1 == self.number_file_1 and other.average_file_1 == self.average_file_1 and\
			other.number_file_2 == self.number_file_2 and other.average_file_2 == self.average_file_2



class CustomEncoderResult(json.JSONEncoder):

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
