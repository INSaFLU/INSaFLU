'''
Created on Dec 14, 2017

@author: mmp
'''

class ConstantsSettings(object):
	'''
	classdocs
	'''
	PIPELINE_NAME_read_quality_analysis = 'Read quality analysis and improvement'
	PIPELINE_NAME_type_and_subtype_analysis = 'Classification'
	PIPELINE_NAME_variant_detection = 'Mutation detection and consensus generation'
	PIPELINE_NAME_coverage_analysis = 'Coverage analysis'
	PIPELINE_NAME_alignment = 'Alignment/Phylogeny'
	PIPELINE_NAME_intra_host_minor_variant_detection = 'Intra-host minor variant detection'
	
	## values to upload to database
	vect_pipeline_names = [
						PIPELINE_NAME_read_quality_analysis,
						PIPELINE_NAME_type_and_subtype_analysis,
						PIPELINE_NAME_variant_detection,
						PIPELINE_NAME_coverage_analysis,
						PIPELINE_NAME_alignment,
						PIPELINE_NAME_intra_host_minor_variant_detection,
					]

	###############################
	### technology available
	TECHNOLOGY_illumina = "Illumina"
	TECHNOLOGY_minion = "ONT"
	TECHNOLOGY_all = "All"
	TECHNOLOGY_Undefined = "Undefined"

	### technology names
	vect_technology = [
					TECHNOLOGY_illumina,
					TECHNOLOGY_minion,
					TECHNOLOGY_Undefined,
					TECHNOLOGY_all,
				]
	
	
	###### Relation between software and technology
	def get_list_software_names_by_technology(self, technology_name):
		"""
		:param technology_name TECHNOLOGY_illumina, TECHNOLOGY_minion
		"""
		return self.DICT_SOFTWARE_RELATION.get(technology_name, [])
	
	def is_software_in_technology(self, technology_name, software_name):
		"""
		:param technology_name TECHNOLOGY_illumina, TECHNOLOGY_minion
		"""
		return software_name in self.get_list_software_names_by_technology(technology_name)
