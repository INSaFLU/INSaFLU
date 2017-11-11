'''
Created on Nov 10, 2017

@author: mmp
'''

class MetaKeyAndValue(object):
	'''
	has all the constants used in MetaKeySample and MetaTag
	'''

	### MetaKeys
	META_KEY_Identify_Sample = "IdentifySamples"					## Used to identify the status on run software.Software.identify_type_and_sub_type 
																	## Value: META_VALUE_Error|META_VALUE_Success; Description: result.Result
	META_KEY_Fastq_Trimmomatic = "FastqTrimmomatic"					## Used to run fastq and trimmomatic software.Software.run_fastq_and_trimmomatic 
																	## Value: META_VALUE_Error|META_VALUE_Success; Description: result.Result
	META_KEY_Number_And_Average_Reads = "NumberAndAverageReads"  	## Used to count the number of the reads on run software.Software.get_lines_and_average_reads 
																	## Value: META_VALUE_Error|META_VALUE_Success; Description: result.Result
	
	## meta value
	META_VALUE_Error = "Error"
	META_VALUE_Success = "Success"

	def __init__(self):
		'''
		Constructor
		'''
		pass