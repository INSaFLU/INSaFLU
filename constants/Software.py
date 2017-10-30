'''
Created on Oct 28, 2017

@author: mmp
'''
import os
import logging

class Software(object):
	'''
	classdocs
	'''

	## dir with software
	DIR_SOFTWARE = "/usr/local/software"
	SOFTWARE_SAMTOOLS = "/usr/bin/samtools"
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
		
	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	"""
	return samtools software
	"""
	def getSamtools(self): return self.SOFTWARE_SAMTOOLS
	
	
	def createFaiToFastaFile(self, fileFastaName):
		"""
		Create fai for a fasta file
		"""
		cmd = "%s faidx %s" % (self.getSamtools(), fileFastaName)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to run samtools") 
		