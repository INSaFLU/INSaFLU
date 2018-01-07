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
	META_KEY_Tree_Count_All_Sequences = "Tree_Count_All_Sequences"		## has the number of samples processed to build the tree
	META_KEY_Tree_Count_By_Element = "Tree_Count_By_Element"			## has the number of samples processed to build the tree by element
	META_KEY_Tree_Count_Protein_By_Element = "Tree_Count_Protein By_Element"			## has the number of samples processed to build the tree by protein element
	META_KEY_Run_Tree_All_Sequences = "Tree_All_Sequences"				## Used to identify the status on run tree.CreateTree.create_tree_and_alignments_all 
	META_KEY_Run_Tree_By_Element = "Tree_By_Element"					## Used to identify the status on run tree.CreateTree.create_tree_and_alignments_by_sequence 
	META_KEY_Run_Proteins_Alignment_By_Element = "Proteins_Alignment_By_Element"		## Used to identify the status on run tree.CreateTree.create_alignement_for_element 
	META_KEY_Count_Samples_Var_Graph = "Count_samples_var_graph"		## has the number of samples in var graph
	META_KEY_Snippy = "Snippy"
	META_KEY_Freebayes = "Freebayes"
	META_KEY_Coverage = "Coverage"
	META_KEY_Mixed_Infection = "Mixed Infection"
	
	## metakey about the alerts
	META_KEY_ALERT_COVERAGE = "Alert Coverage"
	META_KEY_ALERT_COVERAGE_9 = META_KEY_ALERT_COVERAGE + " >9"
	META_KEY_ALERT_COVERAGE_0 = META_KEY_ALERT_COVERAGE + " >0"
	META_KEY_ALERT_MIXED_INFECTION_COSINE_DISTANCE = "AlertMixedInfectionCosineDistance"
	META_KEY_ALERT_MIXED_INFECTION_RATIO_TEST = "AlertMixedInfectionRatioTest"
	META_KEY_ALERT_COUNT_VAR = "Alert count variation percentile"	## obsolete
	
	## coverage about bam file
	META_KEY_Coverage = "Coverage"
																	
	### Task ID by Job
	META_KEY_Queue_TaskID = "QueueTaskID"		### Global queueTaskID
	
	### Queue Task ID by Job plus project_sample_id
	## to check if this project_sample_ID
	META_KEY_Queue_TaskID_ProjectSample = "TaskID_ProjectSample"	### has the metaKey plus project_sample_id
	META_KEY_Queue_TaskID_Project = "TaskID_Project"				### has the metaKey plus project_sample_id
	META_KEY_Elements_Reference = "Elements"						### has the elements to a specific project, separated by comma and sorted
	META_KEY_Elements_And_CDS_Reference = "Elements And CDS"		### has the elements and the CDS
	
	### has the max sample length for a specific project
	META_KEY_Project_max_sample_length = "Project_max_sample_length"### has the max sample length for a specific project

	### MetaKey for the alerts
	META_KEY_Alert_First_level = "AlertFirstLevel"			### has the alert description for highest level alert
	META_KEY_Alert_Second_level = "AlertSecondLevel"		### has the alert description for lowest level alert
	
	## meta value
	META_VALUE_Error = "Error"
	META_VALUE_Success = "Success"
	META_VALUE_Queue = "Queue"

	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	def get_meta_key(self, meta_key, entity):
		"""
		in: MetaKeyAndValue.META_KEY_Tree_All_Sequences, MetaKeyAndValue.META_KEY_Run_Tree_By_Element,
				META_KEY_ALERT_COVERAGE_9, META_KEY_ALERT_COVERAGE_0
				META_KEY_Queue_TaskID_Project, META_KEY_Elements_Project
		return metakey by element
		"""
		return '{} {}'.format(meta_key, entity)
	
	def get_meta_key_queue_by_project_sample_id(self, project_sample_id):
		"""
		in: project_sample_id
		return META_KEY_Queue_TaskID_ProjectSample + project_sample_id
		Can assume the value: META_VALUE_Queue, META_VALUE_Success, META_VALUE_Error
		Each project_sample_id as a metakey, that is the final value of META_VALUE_Success or META_VALUE_Error
		"""
		return self.get_meta_key(self.META_KEY_Queue_TaskID_ProjectSample, project_sample_id)
	
	def get_meta_key_queue_by_project_id(self, project_id):
		"""
		in: project_id
		return META_KEY_Queue_TaskID_Project + project_sample_id
		Can assume the value: META_VALUE_Queue, META_VALUE_Success, META_VALUE_Error
		Each project_id as a metakey, that is the final value of META_VALUE_Success or META_VALUE_Error
		"""
		return self.get_meta_key(self.META_KEY_Queue_TaskID_Project, project_id)
	


