'''
Created on Nov 10, 2017

@author: mmp
'''
from settings.constants_settings import ConstantsSettings

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
	META_KEY_NanoStat_NanoFilt = "NanoStatNanoFilt"						## Used to run NanoStat and NanoFilt errors 
	META_KEY_NanoStat_NanoFilt_Software = "NanoStatNanoFiltSoftware"	## Used to run NanoStat and NanoFilt 
	META_KEY_Number_And_Average_Reads = "NumberAndAverageReads"  		## Used to count the number of the reads on run software.Software.get_lines_and_average_reads 
																		## Value: META_VALUE_Error|META_VALUE_Success; Description: result.Result
	META_KEY_Snippy_Freebayes = "SnippyFreeBayes"						## Used to run snippy and freebayes 
																		## Value: META_VALUE_Error|META_VALUE_Success; Description: result.Result
	META_KEY_Identify_pangolin = "IdentifyPangolin"						## Tag to use in pangolin
																		## Value: META_VALUE_Error|META_VALUE_Success; Description: result.Result
	META_KEY_Count_Hits = "Count Hits"									## Has the hits (50-90) (<50)
	META_KEY_bam_stats = "Bam statistics"								## Has the bam statistics, reads available, mapped
	META_KEY_Tree_Count_All_Sequences = "Tree_Count_All_Sequences"		## has the number of samples processed to build the tree
	META_KEY_Tree_Count_By_Element = "Tree_Count_By_Element"			## has the number of samples processed to build the tree by element
	META_KEY_Tree_Count_Protein_By_Element = "Tree_Count_Protein By_Element"			## has the number of samples processed to build the tree by protein element
	META_KEY_Run_Tree_All_Sequences = "Tree_All_Sequences"				## Used to identify the status on run tree.CreateTree.create_tree_and_alignments_all 
	META_KEY_Run_Tree_By_Element = "Tree_By_Element"					## Used to identify the status on run tree.CreateTree.create_tree_and_alignments_by_sequence 
	META_KEY_Run_Proteins_Alignment_By_Element = "Proteins_Alignment_By_Element"		## Used to identify the status on run tree.CreateTree.create_alignement_for_element 
	META_KEY_Count_Samples_Var_Graph = "Count_samples_var_graph"		## has the number of samples in var graph
	META_KEY_Snippy = "Snippy"
	META_KEY_Medaka = "Medaka"
	META_KEY_Freebayes = "Freebayes"
	META_KEY_Coverage = "Coverage"
	META_KEY_Mixed_Infection = "Mixed Infection"
	
	## metakey about the alerts
	META_KEY_ALERT_COVERAGE = "Alert Coverage"
	META_KEY_ALERT_COVERAGE_value_defined_by_user = META_KEY_ALERT_COVERAGE + " >X"
	META_KEY_ALERT_COVERAGE_9 = META_KEY_ALERT_COVERAGE + " >9"
	META_KEY_ALERT_COVERAGE_0 = META_KEY_ALERT_COVERAGE + " >0"
	META_KEY_ALERT_MIXED_INFECTION_COSINE_DISTANCE = "AlertMixedInfectionCosineDistance"
	META_KEY_ALERT_MIXED_INFECTION_RATIO_TEST = "AlertMixedInfectionRatioTest"
	META_KEY_ALERT_MIXED_INFECTION_SUM_TEST = "AlertMixedInfectionSumTest"
	META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE = "AlertMixedInfectionTypeSubtype"
	META_KEY_TAG_MIXED_INFECTION_TYPE_SUBTYPE = "TagMixedInfectionTypeSubtype"
	META_KEY_ALERT_NOT_ASSIGNED_TYPE_SUBTYPE = "AlertNotAssignedTypeSubtype"
	META_KEY_ALERT_COUNT_VAR = "Alert count variation percentile"				## obsolete
	META_KEY_ALERT_DOWNSIZE_OF_FASTQ_FILES = "Downsize fastq files"
	META_KEY_ALERT_REMOVE_ORIGINAL_FASTQ_FILES = "Remove Original Fastq Files"	### remove original fastq files
	META_KEY_ALERT_NO_READS_AFTER_FILTERING = "No reads after filtering."
	
	#### KEYS to remove if we run Snippy And FreeBayes Again
	VECT_TO_REMOVE_RUN_PROJECT_SAMPLE = [META_KEY_ALERT_COVERAGE,\
						META_KEY_ALERT_MIXED_INFECTION_RATIO_TEST,\
						META_KEY_ALERT_MIXED_INFECTION_SUM_TEST,\
						META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE]

	#### KEYS to remove if we run trimmomatic or NanoFilt again
	VECT_TO_REMOVE_RUN_SAMPLE = [
						META_KEY_TAG_MIXED_INFECTION_TYPE_SUBTYPE,\
						META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE,\
						META_KEY_ALERT_NO_READS_AFTER_FILTERING,\
						META_KEY_Identify_Sample]

	### VECT message for show in web page
	VECT_MESSAGE_ALERT_COVERAGE = [META_KEY_ALERT_COVERAGE_value_defined_by_user,\
						META_KEY_ALERT_COVERAGE_9,\
						META_KEY_ALERT_COVERAGE_0]

	### Key softwares to show in description of tab/csv file
	### { SAMPLE/PROJECT_SAMPLE : [KEY_TO_SHOW, KEY_TO_SHOW_1, KEY_TO_SHOW_2, ...],
	###   SAMPLE/PROJECT_SAMPLE : [KEY_TO_SHOW_2, KEY_TO_SHOW_4, KEY_TO_SHOW_2, ...]
	### }
	NAME_sample = "Sample"						### if the software was used in sample
	NAME_project_sample = "ProjectSample"		### if the software was used in project sample
	
	### order out
	VECT_TECHNOLOGIES_OUT_REPORT = [ConstantsSettings.TECHNOLOGY_illumina,
								ConstantsSettings.TECHNOLOGY_minion]
	### key to show in report
	DICT_SOFTWARE_SHOW_IN_RESULTS = {
		ConstantsSettings.TECHNOLOGY_illumina : {
			NAME_sample : [META_KEY_Fastq_Trimmomatic_Software, META_KEY_Identify_Sample_Software, ],
			NAME_project_sample : [META_KEY_Snippy_Freebayes, ],
		},
		ConstantsSettings.TECHNOLOGY_minion : {
			NAME_sample : [	META_KEY_NanoStat_NanoFilt_Software, ],
			NAME_project_sample : [ META_KEY_Medaka, ],
		}
	}
# 	DICT_SOFTWARE_SHOW_IN_SAMPLE_RESULTS = { 
# 		NAME_sample : [	META_KEY_NanoStat_NanoFilt_Software, ],
# 		NAME_project_sample : [ META_KEY_Medaka, META_KEY_limit_ONT_coverage ],
# 	}
	## coverage about bam file
	META_KEY_Coverage = "Coverage"
								
	## tag for masking consensus values
	META_KEY_Masking_consensus = "Masking consensus"
	META_KEY_Masking_consensus_by_minfrac_VCF_medaka = "Masking consensus by Minfrac VCF medaka"
										
	### Task ID by Job
	META_KEY_Queue_TaskID = "QueueTaskID"		### Global queueTaskID
	
	### Queue Task ID by Job plus project_sample_id
	## to check if this project_sample_ID
	META_KEY_Queue_TaskID_ProjectSample = "TaskID_ProjectSample"	### has the metaKey plus project_sample_id
	META_KEY_Queue_TaskID_Project = "TaskID_Project"				### has the metaKey plus project_sample_id
	META_KEY_Elements_Reference = "Elements"						### has the elements to a specific project, separated by comma and sorted
	META_KEY_Elements_And_CDS_Reference = "Elements And CDS"		### has the elements and the CDS
	
	### has the max sample length for a specific project
	META_KEY_Project_max_sample_length = "Project_max_sample_length"	### has the max sample length for a specific project
	META_KEY_Dataset_max_name_length = "Dataset_max_name_length"			### has the max sample length for a specific dataset

	### MetaKey for the alerts
	META_KEY_Alert_First_level = "AlertFirstLevel"			### has the alert description for highest level alert
	META_KEY_Alert_Second_level = "AlertSecondLevel"		### has the alert description for lowest level alert
	
	## meta value
	META_VALUE_Error = "Error"
	META_VALUE_Success = "Success"
	META_VALUE_Queue = "Queue"

	## Keys for statistics samtools mapped reads
	SAMTOOLS_flagstat_total_reads = "Total reads"
	SAMTOOLS_flagstat_mapped_reads = "Mapped reads"
	
	
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
	
	def get_keys_show_alerts_in_sample_projects_details_view(self):
		"""
		this show the alerts in SampleProjectsDetailsView
		"""
		return [MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_RATIO_TEST, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_SUM_TEST,\
			MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE]
		
	def get_keys_show_alerts_in_sample_details_view(self):
		"""
		this show the alerts in SampleDetailsView
		"""
		return [MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE, MetaKeyAndValue.META_KEY_ALERT_NOT_ASSIGNED_TYPE_SUBTYPE,\
				MetaKeyAndValue.META_KEY_ALERT_DOWNSIZE_OF_FASTQ_FILES, MetaKeyAndValue.META_KEY_ALERT_REMOVE_ORIGINAL_FASTQ_FILES,\
				MetaKeyAndValue.META_KEY_ALERT_NO_READS_AFTER_FILTERING]
