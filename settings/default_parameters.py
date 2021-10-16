'''
Created on 10/04/2021

@author: mmp
'''
import os
from django.conf import settings
from settings.models import Software, Parameter, Technology
from constants.software_names import SoftwareNames
from utils.lock_atomic_transaction import LockedAtomicTransaction

class DefaultParameters(object):
	'''
	Define default values for all softwares
	'''
	software_names = SoftwareNames()

	### used in snippy
	SNIPPY_COVERAGE_NAME = "--mincov"
	SNIPPY_MAPQUAL_NAME = "--mapqual"

	### used in NANOfilt
	NANOfilt_quality_read = "-q"
	
	### used in Medaka
	MEDAKA_model = "-m"
	
	### used in mask consensus
	MASK_CONSENSUS_threshold = "Threshold"

	### clean human reads
	MASK_CLEAN_HUMAN_READS = software_names.SOFTWARE_CLEAN_HUMAN_READS_name
	
	### MINIMUN of MAX of NanoFilt
	NANOFILT_MINIMUN_MAX = 100
	
	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	def get_software_parameters_version(self, software_name):
		"""
		With time it is possible to add or remove parameters available in softwares
		Add one every time that is necessary to add or remove parameters in software 
		"""
		if (software_name == SoftwareNames.SOFTWARE_TRIMMOMATIC_name): return 1
		## all others remain zero, didn't change anything
		return 0
	
	def _get_software_obsolete(self):
		"""
		Return all softwares/parameters version obsolete 
		"""
		vect_return = []
		vect_return.append([SoftwareNames.SOFTWARE_TRIMMOMATIC_name, 0])
		return vect_return
	
	def set_software_obsolete(self):
		"""
		set software obsolete
		"""
		vect_to_obsolete = self._get_software_obsolete()
		with LockedAtomicTransaction(Software):
			for data_obsolete in vect_to_obsolete:
				query_set = Software.objects.filter(name=data_obsolete[0],
					version_parameters = data_obsolete[1], is_obsolete=False)
				for software in query_set:
					software.is_obsolete = True
					software.save()
		
			
	def get_technology_instance(self, technology_name):
		"""
		:out return an instance with technology name
		"""
		with LockedAtomicTransaction(Technology):
			try:
				technology = Technology.objects.get(name=technology_name)
			except Technology.DoesNotExist:
				technology = Technology()
				technology.name = technology_name
				technology.save()
			return technology

	def persist_parameters(self, vect_parameters, type_of_use, technology_name):
		"""
		presist a specific software by default
		param: type_of_use Can by Software.TYPE_OF_USE_project; Software.TYPE_OF_USE_project_sample
		"""
		software = None
		dt_out_sequential = {}
		for parameter in vect_parameters:
			assert parameter.sequence_out not in dt_out_sequential
			if software is None:
				try:
					software = Software.objects.get(name=parameter.software.name, owner=parameter.software.owner,
						type_of_use=type_of_use, technology__name=technology_name,
						version_parameters = parameter.software.version_parameters)
				except Software.DoesNotExist:
					software = parameter.software
					software.technology = self.get_technology_instance(technology_name)
					software.save()
			parameter.software = software
			parameter.save()
			
			## set sequential number
			dt_out_sequential[parameter.sequence_out] = 1
	
	#####################################
	#### TODO
	# update settings_parameter set software_id = 229 where software_id = 50
	# update settings_parameter set software_id = 228 where software_id = 49
	# delete from settings_software where id = 49
	# delete from settings_software where id = 50
	
	# select * from settings_software where owner_id=38 and technology_id is not null order by name, type_of_use
	#
	# select * from settings_parameter where software_id = 50 order by id
	# select * from settings_parameter where software_id = 229  order by id
	# select * from settings_parameter where software_id = 49  order by id
	# select * from settings_parameter where software_id = 228  order by id

	def get_parameters(self, software_name, user, type_of_use, project, project_sample,
				sample, technology_name = SoftwareNames.TECHNOLOGY_illumina):
		"""
		get software_name parameters
		"""
		try:
			software = Software.objects.get(name=software_name, owner=user,\
						type_of_use = type_of_use,
						technology__name = technology_name,
						version_parameters = self.get_software_parameters_version(software_name))
		except Software.DoesNotExist:
			if (type_of_use == Software.TYPE_OF_USE_global):
				try:
					software = Software.objects.get(name=software_name, owner=user,\
							type_of_use=type_of_use,
							version_parameters = self.get_software_parameters_version(software_name))
				except Software.DoesNotExist:
					return None
			else: return None

		## get parameters for a specific user
		parameters = Parameter.objects.filter(software=software, project=project,
						project_sample=project_sample, sample=sample)
		
		### parse them
		dict_out = {}
		vect_order_ouput = []
		for parameter in parameters:
			### don't set the not set parameters
			if (not parameter.not_set_value is None and parameter.parameter == parameter.not_set_value): continue
			
			### create a dict with parameters
			if (parameter.name in dict_out): 
				dict_out[parameter.name][1].append(parameter.parameter)
				dict_out[parameter.name][0].append(parameter.union_char)
			else:
				vect_order_ouput.append(parameter.name) 
				dict_out[parameter.name] = [[parameter.union_char], [parameter.parameter]]
			
		return_parameter = ""
		for par_name in vect_order_ouput:
			if (len(dict_out[par_name][1]) == 1 and len(dict_out[par_name][1][0]) == 0):
				return_parameter += " {}".format(par_name)
			else:
				return_parameter += " {}".format(par_name)
				### exception SOFTWARE_TRIMMOMATIC_illuminaclip SOFTWARE_TRIMMOMATIC_name
				if (software_name == SoftwareNames.SOFTWARE_TRIMMOMATIC_name and \
					par_name == SoftwareNames.SOFTWARE_TRIMMOMATIC_illuminaclip):
					return_parameter += "{}{}{}".format(dict_out[par_name][0][0],
									os.path.join(settings.DIR_SOFTWARE, "trimmomatic/adapters",
									dict_out[par_name][1][0]),
									SoftwareNames.SOFTWARE_TRIMMOMATIC_addapter_trim_used_to_assemble)
				else:
					for _ in range(len(dict_out[par_name][0])):
						return_parameter += "{}{}".format(dict_out[par_name][0][_], dict_out[par_name][1][_])
					
		#### This is the case where all the options can be "not to set"
		if (len(return_parameter.strip()) == 0 and len(parameters) == 0): return None
		return return_parameter.strip()

	def get_snippy_default(self, user, type_of_use, project = None, project_sample = None):
		"""
		–mapqual: minimum mapping quality to allow (–mapqual 20)
		—mincov: minimum coverage of variant site (–mincov 10)
		–minfrac: minumum proportion for variant evidence (–minfrac 0.51)
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_SNIPPY_name
		software.name_extended = SoftwareNames.SOFTWARE_SNIPPY_name_extended
		software.version = SoftwareNames.SOFTWARE_SNIPPY_VERSION
		software.type_of_use = type_of_use
		software.type_of_software = Software.TYPE_SOFTWARE
		software.version_parameters = self.get_software_parameters_version(software.name)
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = DefaultParameters.SNIPPY_MAPQUAL_NAME
		parameter.parameter = "20"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[5:100]"
		parameter.range_max = "50"
		parameter.range_min = "10"
		parameter.description = "MAPQUAL: is the minimum mapping quality to accept in variant calling (–mapqual 20)."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = DefaultParameters.SNIPPY_COVERAGE_NAME
		parameter.parameter = "10"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 2
		parameter.range_available = "[4:100]"
		parameter.range_max = "100"
		parameter.range_min = "4"
		parameter.description = "MINCOV: the minimum number of reads covering a site to be considered (–mincov 10)."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--minfrac"
		parameter.parameter = "0.51"
		parameter.type_data = Parameter.PARAMETER_float
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 3
		parameter.range_available = "[0.5:1.0]"
		parameter.range_max = "1.0"
		parameter.range_min = "0.5"
		parameter.description = "MINFRAC: minimum proportion for variant evidence (–minfrac 0.51)"
		vect_parameters.append(parameter)
		
		return vect_parameters

	def get_freebayes_default(self, user, type_of_use, project = None, project_sample = None):
		"""
		–min-mapping-quality: excludes read alignments from analysis if they have a mapping quality less than Q (–min-mapping-quality 20)
		—min-base-quality: excludes alleles from iSNV analysis if their supporting base quality is less than Q (–min-base-quality 20)
		–min-coverage: requires at least 100-fold of coverage to process a site (–min-coverage 100)
		—min-alternate-count: require at least 10 reads supporting an alternate allele within a single individual in order to evaluate the position (–min-alternate-count 10)
		–min-alternate-fraction: defines the minimum intra-host frequency of the alternate allele to be assumed (–min-alternate-fraction 0.01). This frequency is contingent on the depth of coverage of each processed site since min-alternate-count is set to 10, i.e., the identification of iSNV sites at frequencies of 10%, 2% and 1% is only allowed for sites with depth of coverage of at least 100-fold, 500-fold and 1000-fold, respectively.
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_FREEBAYES_name
		software.name_extended = SoftwareNames.SOFTWARE_FREEBAYES_name_extended
		software.version = SoftwareNames.SOFTWARE_FREEBAYES_VERSION
		software.type_of_use = type_of_use
		software.type_of_software = Software.TYPE_SOFTWARE
		software.version_parameters = self.get_software_parameters_version(software.name)
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "--min-mapping-quality"
		parameter.parameter = "20"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[10:50]"
		parameter.range_max = "50"
		parameter.range_min = "10"
		parameter.not_set_value = "0"
		parameter.description = "min-mapping-quality: excludes read alignments from analysis if they have a mapping quality less than Q (–min-mapping-quality 20)."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--min-base-quality"
		parameter.parameter = "20"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 2
		parameter.range_available = "[10:50]"
		parameter.range_max = "50"
		parameter.range_min = "10"
		parameter.not_set_value = "0"
		parameter.description = "min-base-quality: excludes alleles from iSNV analysis if their supporting base quality is less than Q (–min-base-quality 20)."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--min-coverage"
		parameter.parameter = "100"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 3
		parameter.range_available = "[20:500]"
		parameter.range_max = "500"
		parameter.range_min = "20"
		parameter.description = "min-coverage: requires at least 100-fold of coverage to process a site (–min-coverage 100)."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--min-alternate-count"
		parameter.parameter = "10"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 4
		parameter.range_available = "[5:100]"
		parameter.range_max = "100"
		parameter.range_min = "5"
		parameter.description = "min-alternate-count: require at least 10 reads supporting an alternate allele within a single " +\
			"individual in order to evaluate the position (–min-alternate-count 10)."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--min-alternate-fraction"
		parameter.parameter = "100"
		parameter.type_data = Parameter.PARAMETER_float
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 5
		parameter.range_available = "[0.01:0.5]"
		parameter.range_max = "0.5"
		parameter.range_min = "0.01"
		parameter.description = "min-alternate-fraction: defines the minimum intra-host frequency of the alternate allele to be assumed (–min-alternate-fraction 0.01)."
		vect_parameters.append(parameter)

		parameter = Parameter()
		parameter.name = "--ploidy"
		parameter.parameter = "2"
		parameter.type_data = Parameter.PARAMETER_float
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 6
		parameter.range_available = "[1:5]"
		parameter.range_max = "5"
		parameter.range_min = "1"
		parameter.description = "ploidy:  Sets the default ploidy for the analysis to N.  default: 2"
		vect_parameters.append(parameter)
		
		return vect_parameters

	def get_mask_consensus_threshold_default(self, user, type_of_use, project = None, project_sample = None):
		"""
		Threshold of mask not consensus coverage
		"""
		software = Software()
		software.name = SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name
		software.name_extended = SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name_extended
		software.type_of_use = type_of_use
		software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
		software.version = "1.0"
		software.version_parameters = self.get_software_parameters_version(software.name)
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = DefaultParameters.MASK_CONSENSUS_threshold
		parameter.parameter = "70"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[50:100]"
		parameter.range_max = "100"
		parameter.range_min = "50"
		parameter.description = "Minimum percentage of locus horizontal coverage with depth of coverage equal or above –mincov (see Snippy) to generate consensus sequence."
		vect_parameters.append(parameter)
		return vect_parameters

	def get_clean_human_reads_default(self, user, type_of_use):
		"""
		Threshold of mask not consensus coverage
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_CLEAN_HUMAN_READS_name
		software.name_extended = SoftwareNames.SOFTWARE_CLEAN_HUMAN_READS_extended
		software.type_of_use = type_of_use
		software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
		software.version = "1.0"
		software.version_parameters = self.get_software_parameters_version(software.name)
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = DefaultParameters.MASK_CLEAN_HUMAN_READS
		parameter.parameter = SoftwareNames.SOFTWARE_TAG_no
		parameter.type_data = Parameter.PARAMETER_char_list
		parameter.software = software
		parameter.union_char = ": "
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.description = "Clean human reads from fastq files."
		vect_parameters.append(parameter)
		return vect_parameters
	
	def get_limit_coverage_ONT_threshold_default(self, user, type_of_use, project = None, project_sample = None):
		"""
		Minimum depth of coverage per site to validate the sequence (default: –mincov 30)
		Where to use this cut-off:
			This cut-off is used to exclude from vcf files sites with DEPTH <=30
			This cut-off is used to mask consensus sequences with DEPTH <=30 in msa_masker (-c: 30-1 = 29)
		"""
		software = Software()
		software.name = SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name
		software.name_extended = SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name_extended
		software.type_of_use = type_of_use
		software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
		software.version = "1.0"
		software.version_parameters = self.get_software_parameters_version(software.name)
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "Threshold"
		parameter.parameter = "30"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[4:100]"
		parameter.range_max = "100"
		parameter.range_min = "4"
		parameter.description = "This cut-off is used to exclude from vcf files sites and to mask consensus sequences. Value in percentage"
		vect_parameters.append(parameter)
		return vect_parameters

	def get_vcf_freq_ONT_threshold_default(self, user, type_of_use, project = None, project_sample = None):
		"""
		MINFRAC: minumum proportion for variant evidence (–minfrac 51) Range: [10:100]
		
		"""
		software = Software()
		software.name = SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name
		software.name_extended = SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name_extended
		software.type_of_use = type_of_use
		software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
		software.version = "1.0"
		software.version_parameters = self.get_software_parameters_version(software.name)
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "Threshold"
		parameter.parameter = "0.80"
		parameter.type_data = Parameter.PARAMETER_float
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[0.10:1.0]"
		parameter.range_max = "1.0"
		parameter.range_min = "0.10"
		parameter.description = "MINFRAC: minimum proportion for variant evidence (–minfrac) Range: [0.1:1.0]"
		vect_parameters.append(parameter)
		return vect_parameters

	def get_medaka_model_default(self, user, type_of_use, project = None, project_sample = None):
		"""
		Minimum depth of coverage per site to validate the sequence (default: –mincov 30)
		Where to use this cut-off:
			This cut-off is used to exclude from vcf files sites with DEPTH <=30
			This cut-off is used to mask consensus sequences with DEPTH <=30 in msa_masker (-c: 30-1 = 29)
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_Medaka_name_consensus
		software.name_extended = SoftwareNames.SOFTWARE_Medaka_name_extended_consensus
		software.type_of_use = type_of_use
		software.type_of_software = Software.TYPE_SOFTWARE
		software.version = SoftwareNames.SOFTWARE_Medaka_VERSION
		software.version_parameters = self.get_software_parameters_version(software.name)
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "-m"
		parameter.parameter = SoftwareNames.SOFTWARE_Medaka_default_model	## default model
		parameter.type_data = Parameter.PARAMETER_char_list
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.description = "Medaka models are named to indicate: " +\
			"i) the pore type; " +\
			"ii) the sequencing device (min -> MinION, prom -> PromethION); " +\
			"iii) the basecaller variant (only high and variant available in INSAFlu); " +\
			"iv) the Guppy basecaller version. " +\
			"Complete format: " +\
			"{pore}_{device}_{caller variant}_{caller version}"
		vect_parameters.append(parameter)
		return vect_parameters

	def get_nanofilt_default(self, user, type_of_use, sample = None):
		"""
		-l <LENGTH>, Filter on a minimum read length. Range: [50:1000].
		--maxlength Filter on a maximum read length
		dá para colocar outro parametro, por exemplo: -ml <MAX_LENGTH>, Set the maximum read length.
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_NanoFilt_name
		software.name_extended = SoftwareNames.SOFTWARE_NanoFilt_name_extended
		software.version = SoftwareNames.SOFTWARE_NanoFilt_VERSION
		software.type_of_use = type_of_use
		software.type_of_software = Software.TYPE_SOFTWARE
		software.version_parameters = self.get_software_parameters_version(software.name)
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "-q"
		parameter.parameter = "10"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[5:30]"
		parameter.range_max = "30"
		parameter.range_min = "5"
		parameter.description = "-q <QUALITY>, Filter on a minimum average read quality score."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "-l"
		parameter.parameter = "50"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 2
		parameter.range_available = "[50:1000]"
		parameter.range_max = "1000"
		parameter.range_min = "50"
		parameter.description = "-l <LENGTH>, Filter on a minimum read length."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--headcrop"
		parameter.parameter = "70"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 3
		parameter.range_available = "[1:1000]"
		parameter.range_max = "1000"
		parameter.range_min = "1"
		parameter.not_set_value = "0"
		parameter.description = "--headcrop <HEADCROP>, Trim n nucleotides from start of read."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--tailcrop"
		parameter.parameter = "70"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 4
		parameter.range_available = "[1:1000]"
		parameter.range_max = "1000"
		parameter.range_min = "1"
		parameter.not_set_value = "0"
		parameter.description = "--tailcrop <TAILCROP>, Trim n nucleotides from end of read."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--maxlength"
		parameter.parameter = "0"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 5
		parameter.range_available = "[{}:50000]".format(DefaultParameters.NANOFILT_MINIMUN_MAX)
		parameter.range_max = "50000"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "--maxlength <LENGTH>, Set a maximum read length."
		vect_parameters.append(parameter)
		
		return vect_parameters

	def get_samtools_depth_default_ONT(self, user, type_of_use, project = None, project_sample = None):
		"""
		samtools depth for ONT
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT
		software.name_extended = SoftwareNames.SOFTWARE_SAMTOOLS_name_extended_depth_ONT
		software.type_of_use = type_of_use
		software.type_of_software = Software.TYPE_SOFTWARE
		software.version = SoftwareNames.SOFTWARE_SAMTOOLS_VERSION
		software.version_parameters = self.get_software_parameters_version(software.name)
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "-q"
		parameter.parameter = "0"	## default model
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[0:100]"
		parameter.range_max = "100"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "-q <Quality> base quality threshold."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "-Q"
		parameter.parameter = "0"	## default model
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 2
		parameter.range_available = "[0:100]"
		parameter.range_max = "100"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "-Q <Quality> mapping quality threshold."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "-aa"
		parameter.parameter = ""	## default model
		parameter.type_data = Parameter.PARAMETER_null
		parameter.software = software
		parameter.project = project
		parameter.project_sample = project_sample
		parameter.union_char = ""
		parameter.can_change = False
		parameter.sequence_out = 3
		parameter.description = "Output absolutely all positions, including unused ref. sequences."
		vect_parameters.append(parameter)
		
		return vect_parameters
	
	def get_trimmomatic_default(self, user, type_of_use, sample = None):
		
		software = Software()
		software.name = SoftwareNames.SOFTWARE_TRIMMOMATIC_name
		software.name_extended = SoftwareNames.SOFTWARE_TRIMMOMATIC_name_extended
		software.version = SoftwareNames.SOFTWARE_TRIMMOMATIC_VERSION
		software.type_of_use = type_of_use
		software.type_of_software = Software.TYPE_SOFTWARE
		software.version_parameters = self.get_software_parameters_version(software.name)
		software.owner = user
		
		vect_parameters =  []
		
		# ILLUMINACLIP:/usr/local/software/bioinformatics/trimmomatic/adapters/adapters.fa:3:30:10:6:true
		parameter = Parameter()
		parameter.name = SoftwareNames.SOFTWARE_TRIMMOMATIC_illuminaclip
		parameter.parameter = SoftwareNames.SOFTWARE_TRIMMOMATIC_addapter_not_apply
		parameter.type_data = Parameter.PARAMETER_char_list
		parameter.software = software
		parameter.sample = sample
		parameter.not_set_value = SoftwareNames.SOFTWARE_TRIMMOMATIC_addapter_not_apply
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.description = "To clip the Illumina adapters from the input file using the adapter sequences.\n" +\
			"ILLUMINACLIP:<ADAPTER_FILE>:3:30:10:6:true"
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "HEADCROP"
		parameter.parameter = "0"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 2
		parameter.range_available = "[0:100]"
		parameter.range_max = "100"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "HEADCROP:<length> Cut the specified number of bases from the start of the read."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "CROP"
		parameter.parameter = "0"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 3
		parameter.range_available = "[0:400]"
		parameter.range_max = "400"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "CROP:<length> Cut the read to a specified length."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "SLIDINGWINDOW"
		parameter.parameter = "5"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 4
		parameter.range_available = "[3:50]"
		parameter.range_max = "50"
		parameter.range_min = "3"
		parameter.description = "SLIDINGWINDOW:<windowSize> specifies the number of bases to average across"
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "SLIDINGWINDOW"
		parameter.parameter = "20"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 5
		parameter.range_available = "[10:100]"
		parameter.range_max = "100"
		parameter.range_min = "10"
		parameter.description = "SLIDINGWINDOW:<requiredQuality> specifies the average quality required"
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "LEADING"
		parameter.parameter = "3"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 6
		parameter.range_available = "[0:100]"
		parameter.range_max = "100"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "LEADING:<quality> Remove low quality bases from the beginning."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "TRAILING"
		parameter.parameter = "3"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 7
		parameter.range_available = "[0:100]"
		parameter.range_max = "100"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "TRAILING:<quality> Remove low quality bases from the end."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "MINLEN"
		parameter.parameter = "35"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 8
		parameter.range_available = "[5:500]"
		parameter.range_max = "500"
		parameter.range_min = "5"
		parameter.description = "MINLEN:<length> This module removes reads that fall below the specified minimal length."
		vect_parameters.append(parameter)

##		Only available in 0.30 version		
#
# 		parameter = Parameter()
# 		parameter.name = "AVGQUAL"
# 		parameter.parameter = "0"
# 		parameter.type_data = Parameter.PARAMETER_int
# 		parameter.software = software
#		parameter.sample = sample
# 		parameter.union_char = ":"
# 		parameter.can_change = True
# 		parameter.sequence_out = 8
# 		parameter.range_available = "[0:100]"
# 		parameter.range_max = "100"
# 		parameter.range_min = "0"
# 		parameter.not_set_value = "0"
# 		parameter.description = "AVGQUAL:<quality> Drop the read if the average quality is below the specified level."
# 		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "TOPHRED33"
		parameter.parameter = ""
		parameter.type_data = Parameter.PARAMETER_null
		parameter.software = software
		parameter.sample = sample
		parameter.union_char = ""
		parameter.can_change = False
		parameter.sequence_out = 9
		parameter.description = "This (re)encodes the quality part of the FASTQ file to base 33."
		vect_parameters.append(parameter)
		return vect_parameters