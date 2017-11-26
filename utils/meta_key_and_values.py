'''
Created on Nov 10, 2017

@author: mmp
'''

class MetaKeyAndValue(object):
	'''
	has all the constants used in MetaKeySample and MetaTag
	'''

	### MetaKeys
	META_KEY_Identify_Sample = "IdentifySamples"						## Used to identify the status on run software.Software.identify_type_and_sub_type 
	META_KEY_Identify_Sample_Software = "IdentifySamplesSoftware"		## Used to identify the status on run software.Software.identify_type_and_sub_type 
																		## Value: META_VALUE_Error|META_VALUE_Success; Description: result.Result
	META_KEY_Fastq_Trimmomatic = "FastqTrimmomatic"						## Used to run fastq and trimmomatic software.Software.run_fastq_and_trimmomatic 
	META_KEY_Fastq_Trimmomatic_Software = "FastqTrimmomaticSoftware"	## Used to run fastq and trimmomatic software.Software.run_fastq_and_trimmomatic 
																		## Value: META_VALUE_Error|META_VALUE_Success; Description: result.Result
	META_KEY_Number_And_Average_Reads = "NumberAndAverageReads"  		## Used to count the number of the reads on run software.Software.get_lines_and_average_reads 
																		## Value: META_VALUE_Error|META_VALUE_Success; Description: result.Result
	META_KEY_Snippy_Freebayes = "SnippyFreeBayes"						## Used to run snippy and freebayes 
																		## Value: META_VALUE_Error|META_VALUE_Success; Description: result.Result
	META_KEY_Count_Hits = "Count Hits"									## Has the hits (50-90) (<50)
	META_KEY_Tree_All_Sequences = "Tree_All_Sequences"					## has the number of samples processed to build the tree
	META_KEY_Tree_By_Element = "Tree_By_Element"						## has the number of samples processed to build the tree by element
	META_KEY_Run_Tree_All_Sequences = "Tree_All_Sequences"				## Used to identify the status on run tree.CreateTree.create_tree_and_alignments_all 
	META_KEY_Run_Tree_By_Element = "Tree_By_Element"					## Used to identify the status on run tree.CreateTree.create_tree_and_alignments_by_sequence 
	
	META_KEY_Snippy = "Snippy"
	META_KEY_Freebayes = "Freebayes"
	META_KEY_Coverage = "Coverage"
	
	## coverage about bam file
	META_KEY_Coverage = "Coverage"
																	
	### Task ID by Job
	META_KEY_Queue_TaskID = "QueueTaskID"
	
	## meta value
	META_VALUE_Error = "Error"
	META_VALUE_Success = "Success"
	META_VALUE_Queue = "Queue"

	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	def get_meta_key_by_element(self, meta_key_element, element):
		"""
		in: MetaKeyAndValue.META_KEY_Tree_All_Sequences, MetaKeyAndValue.META_KEY_Run_Tree_By_Element
		return metakey by element
		
		"""
		return '{}_{}'.format(meta_key_element, element)