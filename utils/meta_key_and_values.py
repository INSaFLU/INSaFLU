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
	META_KEY_Identify_Sample_Software = "IdentifySamplesSoftware"			## Used to identify the status on run software.Software.identify_type_and_sub_type 
																	## Value: META_VALUE_Error|META_VALUE_Success; Description: result.Result
	META_KEY_Fastq_Trimmomatic = "FastqTrimmomatic"					## Used to run fastq and trimmomatic software.Software.run_fastq_and_trimmomatic 
	META_KEY_Fastq_Trimmomatic_Software = "FastqTrimmomaticSoftware"	## Used to run fastq and trimmomatic software.Software.run_fastq_and_trimmomatic 
																	## Value: META_VALUE_Error|META_VALUE_Success; Description: result.Result
	META_KEY_Number_And_Average_Reads = "NumberAndAverageReads"  	## Used to count the number of the reads on run software.Software.get_lines_and_average_reads 
																	## Value: META_VALUE_Error|META_VALUE_Success; Description: result.Result
																	
	### Task ID by Job
	META_KEY_Import_Sample_Import_Queue_TaskID = "SampleImportQueueTaskID"
	
	## meta value
	META_VALUE_Error = "Error"
	META_VALUE_Success = "Success"
	META_VALUE_Queue = "Queue"

	def __init__(self):
		'''
		Constructor
		'''
		pass