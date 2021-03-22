'''
Created on 03/05/2020

@author: mmp
'''
from settings.models import Software, Parameter, Technology
from constants.software_names import SoftwareNames
from utils.lock_atomic_transaction import LockedAtomicTransaction

""""
--mapqual is the minimum mapping quality to accept in variant calling. BWA MEM using 60 to mean a read is "uniquely mapped".

--basequal is minimum quality a nucleotide needs to be used in variant calling. We use 13 which corresponds to error probability of ~5%. It is a traditional SAMtools value.

--maxsoft is how many bases of an alignment to allow to be soft-clipped before discarding the alignment. This is to encourage global over local alignment, and is passed to the samclip tool.

--mincov and --minfrac are used to apply hard thresholds to the variant calling beyond the existing statistical measure.. The optimal values depend on your sequencing depth and contamination rate. Values of 10 and 0.9 are commonly used.
"""

class DefaultSoftware(object):
	'''
	classdocs
	'''
	### MINIMUN of MAX of NanoFilt
	NANOFILT_MINIMUN_MAX = 100
	
	software_names = SoftwareNames()
	
	def __init__(self):
		""" change values """
		self.change_values_software = {}	### the key is the name of the software

	def test_all_defaults(self, user):
		### test all defaults
		self.test_default_db(SoftwareNames.SOFTWARE_TRIMMOMATIC_name, self._get_trimmomatic_default(user), user)
		self.test_default_db(SoftwareNames.SOFTWARE_SNIPPY_name, self._get_snippy_default(user), user)
		self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
				self._get_mask_consensus_threshold_default(user), user)
		
		## ONT software
		self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
				self._get_mask_consensus_threshold_default(user), user,
				SoftwareNames.TECHNOLOGY_minion)
		self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,
				self._get_limit_coverage_ONT_threshold_default(user), user,
				SoftwareNames.TECHNOLOGY_minion)
		self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name,
				self._get_vcf_freq_ONT_threshold_default(user), user,
				SoftwareNames.TECHNOLOGY_minion)
		self.test_default_db(SoftwareNames.SOFTWARE_Medaka_name_consensus,
				self._get_medaka_model_default(user), user,
				SoftwareNames.TECHNOLOGY_minion)
# 		self.test_default_db(SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT,
# 				self._get_samtools_depth_default_ONT(user), user,
# 				SoftwareNames.TECHNOLOGY_minion)
		self.test_default_db(SoftwareNames.SOFTWARE_NanoFilt_name, self._get_nanofilt_default(user), user,
				SoftwareNames.TECHNOLOGY_minion)


	def test_default_db(self, software_name, vect_parameters, user, 
					technology_name = SoftwareNames.TECHNOLOGY_illumina):
		"""
		test if exist, if not persist in database
		"""
		try:
			Software.objects.get(name=software_name, owner=user,\
						type_of_use = Software.TYPE_OF_USE_global,
						technology__name = technology_name)
		except Software.DoesNotExist:
			
			### if it is Minion is because that does not exist at all. 
			### Previous versions didn't have TechnologyName
			if (technology_name == SoftwareNames.TECHNOLOGY_illumina):
				try:
					software = Software.objects.get(name=software_name, owner=user,\
						type_of_use = Software.TYPE_OF_USE_global, 
						technology__name=technology_name)
					### if exist set illumina in technology
					software.technology = self.get_technology_instance(technology_name)
					software.save()
				except Software.DoesNotExist:
					with LockedAtomicTransaction(Software), LockedAtomicTransaction(Parameter):
						self._persist_parameters(vect_parameters, technology_name)
			else:	### do this all the time for 
				### create a default one for this user TECHNOLOGY_minion
				with LockedAtomicTransaction(Software), LockedAtomicTransaction(Parameter):
					self._persist_parameters(vect_parameters, technology_name)

	def get_trimmomatic_parameters(self, user):
		result = self._get_parameters(user, SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
					SoftwareNames.TECHNOLOGY_illumina)
		return "" if result is None else result
	def get_snippy_parameters(self, user):
		result = self._get_parameters(user, SoftwareNames.SOFTWARE_SNIPPY_name,
					SoftwareNames.TECHNOLOGY_illumina)
		return "" if result is None else result
	def get_nanofilt_parameters(self, user):
		result = self._get_parameters(user, SoftwareNames.SOFTWARE_NanoFilt_name,
					SoftwareNames.TECHNOLOGY_minion)
		return "" if result is None else result
	def get_mask_consensus_threshold_parameters(self, user, technology_name):
		result = self._get_parameters(user, SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
								technology_name)
		return "" if result is None else result
	def get_limit_coverage_ONT_parameters(self, user):
		result = self._get_parameters(user, SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,
								SoftwareNames.TECHNOLOGY_minion)
		return "" if result is None else result
	def get_vcf_freq_ONT_parameters(self, user):
		result = self._get_parameters(user, SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name,
								SoftwareNames.TECHNOLOGY_minion)
		return "" if result is None else result
	def get_medaka_parameters_consensus(self, user):
		result = self._get_parameters(user, SoftwareNames.SOFTWARE_Medaka_name_consensus,
								SoftwareNames.TECHNOLOGY_minion)
		return "" if result is None else result
	def get_samtools_parameters_depth_ONT(self, user):
		result = self._get_parameters(user, SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT,
								SoftwareNames.TECHNOLOGY_minion)
		return "" if result is None else result
	
	def _get_parameters(self, user, software_name, technology_name = SoftwareNames.TECHNOLOGY_illumina):
		"""
		get software_name parameters
		"""
		try:
			software = Software.objects.get(name=software_name, owner=user,\
						type_of_use=Software.TYPE_OF_USE_global,
						technology__name=technology_name)
		except Software.DoesNotExist:
			try:
				software = Software.objects.get(name=software_name, owner=user,\
						type_of_use=Software.TYPE_OF_USE_global)
			except Software.DoesNotExist:
				return None

		## get parameters for a specific user
		parameters = Parameter.objects.filter(software=software)
		
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
				for _ in range(len(dict_out[par_name][0])):
					return_parameter += "{}{}".format(dict_out[par_name][0][_], dict_out[par_name][1][_])
					
		#### This is the case where all the options can be "not to set"
		if (len(return_parameter.strip()) == 0 and len(parameters) == 0): return None
		return return_parameter.strip()

	def set_default_software(self, software, user):
		if (software.name == SoftwareNames.SOFTWARE_SNIPPY_name):
			vect_parameters = self._get_snippy_default(user)
		elif (software.name == SoftwareNames.SOFTWARE_TRIMMOMATIC_name):
			vect_parameters = self._get_trimmomatic_default(user)
		elif (software.name == SoftwareNames.SOFTWARE_NanoFilt_name):
			vect_parameters = self._get_nanofilt_default(user)
		elif (software.name == SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name):
			vect_parameters = self._get_mask_consensus_threshold_default(user)
		elif (software.name == SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name):
			vect_parameters = self._get_limit_coverage_ONT_threshold_default(user)
		elif (software.name == SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name):
			vect_parameters = self._get_vcf_freq_ONT_threshold_default(user)
		elif (software.name == SoftwareNames.SOFTWARE_Medaka_name_consensus):
			vect_parameters = self._get_medaka_model_default(user)
		elif (software.name == SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT):
			vect_parameters = self._get_samtools_depth_default_ONT(user)
		else: return
		
		key_value = "{}_{}".format(software.name, SoftwareNames.TECHNOLOGY_illumina if\
					software.technology is None else software.technology.name)
		parameters = Parameter.objects.filter(software=software)
		self.change_values_software[key_value] = False
		for parameter in parameters:
			if parameter.can_change:
				for parameter_to_set_default in vect_parameters:
					if (parameter_to_set_default.sequence_out == parameter.sequence_out):
						###   if change software name
						if (parameter.parameter != parameter_to_set_default.parameter):
							self.change_values_software[key_value] = True
							parameter.parameter = parameter_to_set_default.parameter
							parameter.save()
						break

	def is_change_values_for_software(self, software):
		""" Return if the software has a value changed"""
		key_value = "{}_{}".format(software.name, SoftwareNames.TECHNOLOGY_illumina if\
					software.technology is None else software.technology.name)
		return self.change_values_software.get(key_value, False)

	def _persist_parameters(self, vect_parameters, technology_name):
		"""
		presist a specific software by default 
		"""
		software = None
		technology = self.get_technology_instance(technology_name)
		dt_out_sequential = {}
		for parameter in vect_parameters:
			assert parameter.sequence_out not in dt_out_sequential
			if software is None:
				software = parameter.software
				software.technology = technology 
				software.save()
			parameter.software = software
			parameter.save()
			
			## set sequential number
			dt_out_sequential[parameter.sequence_out] = 1

	def get_parameters(self, software_name, user, technology_name = SoftwareNames.TECHNOLOGY_illumina):
		"""
		"""
		if (software_name == SoftwareNames.SOFTWARE_TRIMMOMATIC_name):
			self.test_default_db(SoftwareNames.SOFTWARE_TRIMMOMATIC_name, self._get_trimmomatic_default(user),
								user, SoftwareNames.TECHNOLOGY_illumina)
			return self.get_trimmomatic_parameters(user)
		if (software_name == SoftwareNames.SOFTWARE_SNIPPY_name):
			self.test_default_db(SoftwareNames.SOFTWARE_SNIPPY_name, self._get_snippy_default(user), user,
								SoftwareNames.TECHNOLOGY_illumina)
			return self.get_snippy_parameters(user)
		if (software_name == SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name):
			self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, self._get_mask_consensus_threshold_default(user),
								user, technology_name)
			return self.get_mask_consensus_threshold_parameters(user, technology_name)
		if (software_name == SoftwareNames.SOFTWARE_NanoFilt_name):
			self.test_default_db(SoftwareNames.SOFTWARE_NanoFilt_name, self._get_nanofilt_default(user), user,
								SoftwareNames.TECHNOLOGY_minion)
			return self.get_nanofilt_parameters(user)
		if (software_name == SoftwareNames.SOFTWARE_Medaka_name_consensus):
			self.test_default_db(SoftwareNames.SOFTWARE_Medaka_name_consensus,
								self._get_medaka_model_default(user), user,
								SoftwareNames.TECHNOLOGY_minion)
			return self.get_medaka_parameters_consensus(user)
		if (software_name == SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT):
			self.test_default_db(software_name, self._get_samtools_depth_default_ONT(user), user,
								SoftwareNames.TECHNOLOGY_minion)
			return self.get_samtools_parameters_depth_ONT(user)
		if (software_name == SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name):
			self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,
								self._get_limit_coverage_ONT_threshold_default(user), user,
								SoftwareNames.TECHNOLOGY_minion)
			return self.get_limit_coverage_ONT_parameters(user)
		if (software_name == SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name):
			self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name,
								self._get_vcf_freq_ONT_threshold_default(user), user,
								SoftwareNames.TECHNOLOGY_minion)
			return self.get_vcf_freq_ONT_parameters(user)
		return ""
		
	def get_all_software(self):
		"""
		get all softwares available by this class
		"""
		vect_software = []
		vect_software.append(self.software_names.get_trimmomatic_name())
		vect_software.append(self.software_names.get_snippy_name())
		vect_software.append(self.software_names.get_NanoFilt_name())
		vect_software.append(self.software_names.get_insaflu_parameter_mask_consensus_name())
		vect_software.append(self.software_names.get_medaka_name_extended_consensus())
		vect_software.append(self.software_names.get_samtools_name_depth_ONT())
		vect_software.append(self.software_names.get_insaflu_parameter_limit_coverage_name())
		vect_software.append(self.software_names.get_insaflu_parameter_freq_vcf_name())
		return vect_software

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

	def _get_snippy_default(self, user):
		"""
		–mapqual: minimum mapping quality to allow (–mapqual 20)
		—mincov: minimum coverage of variant site (–mincov 10)
		–minfrac: minimum proportion for variant evidence (–minfrac 0.51)
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_SNIPPY_name
		software.name_extended = SoftwareNames.SOFTWARE_SNIPPY_name_extended
		software.version = SoftwareNames.SOFTWARE_SNIPPY_VERSION
		software.type_of_use = Software.TYPE_OF_USE_global
		software.type_of_software = Software.TYPE_SOFTWARE
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "--mapqual"
		parameter.parameter = "20"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[5:100]"
		parameter.range_max = "50"
		parameter.range_min = "10"
		parameter.description = "MAPQUAL: is the minimum mapping quality to accept in variant calling (–mapqual 20)."
		vect_parameters.append(parameter)
		
		parameter = Parameter()
		parameter.name = "--mincov"
		parameter.parameter = "10"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
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
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 3
		parameter.range_available = "[0.5:1.0]"
		parameter.range_max = "1.0"
		parameter.range_min = "0.5"
		parameter.description = "MINFRAC: minimum proportion for variant evidence (–minfrac 0.51)"
		vect_parameters.append(parameter)
		
		return vect_parameters

	

	def _get_trimmomatic_default(self, user):
		"""
		SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_TRIMMOMATIC_name
		software.name_extended = SoftwareNames.SOFTWARE_TRIMMOMATIC_name_extended
		software.version = SoftwareNames.SOFTWARE_TRIMMOMATIC_VERSION
		software.type_of_use = Software.TYPE_OF_USE_global
		software.type_of_software = Software.TYPE_SOFTWARE
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "HEADCROP"
		parameter.parameter = "0"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 1
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
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 2
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
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 3
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
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 4
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
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 5
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
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 6
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
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 7
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
		parameter.union_char = ""
		parameter.can_change = False
		parameter.sequence_out = 8
		parameter.description = "This (re)encodes the quality part of the FASTQ file to base 33."
		vect_parameters.append(parameter)
		return vect_parameters

	def _get_nanofilt_default(self, user):
		"""
		-l <LENGTH>, Filter on a minimum read length. Range: [50:1000].
		--maxlength Filter on a maximum read length
		dá para colocar outro parametro, por exemplo: -ml <MAX_LENGTH>, Set the maximum read length.
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_NanoFilt_name
		software.name_extended = SoftwareNames.SOFTWARE_NanoFilt_name_extended
		software.version = SoftwareNames.SOFTWARE_NanoFilt_VERSION
		software.type_of_use = Software.TYPE_OF_USE_global
		software.type_of_software = Software.TYPE_SOFTWARE
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "-q"
		parameter.parameter = "10"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
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
		parameter.union_char = " "
		parameter.can_change = True
		parameter.sequence_out = 5
		parameter.range_available = "[{}:50000]".format(self.NANOFILT_MINIMUN_MAX)
		parameter.range_max = "50000"
		parameter.range_min = "0"
		parameter.not_set_value = "0"
		parameter.description = "--maxlength <LENGTH>, Set a maximum read length."
		vect_parameters.append(parameter)
		
		return vect_parameters

	def _get_mask_consensus_threshold_default(self, user):
		"""
		Threshold of mask not consensus coverage
		"""
		software = Software()
		software.name = SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name
		software.name_extended = SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name_extended
		software.type_of_use = Software.TYPE_OF_USE_global
		software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
		software.version = "1.0"
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "Threshold"
		parameter.parameter = "70"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[50:100]"
		parameter.range_max = "100"
		parameter.range_min = "50"
		parameter.description = "Threshold to mask failed consensus coverage. Value in percentage"
		vect_parameters.append(parameter)
		return vect_parameters
	
	def _get_limit_coverage_ONT_threshold_default(self, user):
		"""
		Minimum depth of coverage per site to validate the sequence (default: –mincov 30)
		Where to use this cut-off:
			This cut-off is used to exclude from vcf files sites with DEPTH <=30
			This cut-off is used to mask consensus sequences with DEPTH <=30 in msa_masker (-c: 30-1 = 29)
		"""
		software = Software()
		software.name = SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name
		software.name_extended = SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name_extended
		software.type_of_use = Software.TYPE_OF_USE_global
		software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
		software.version = "1.0"
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "Threshold"
		parameter.parameter = "30"
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[4:100]"
		parameter.range_max = "100"
		parameter.range_min = "4"
		parameter.description = "This cut-off is used to exclude from vcf files sites and to mask consensus sequences. Value in percentage"
		vect_parameters.append(parameter)
		return vect_parameters

	def _get_vcf_freq_ONT_threshold_default(self, user):
		"""
		MINFRAC: minimum proportion for variant evidence (–minfrac 0.80) Range: [10:100]
		
		"""
		software = Software()
		software.name = SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name
		software.name_extended = SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name_extended
		software.type_of_use = Software.TYPE_OF_USE_global
		software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
		software.version = "1.0"
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "Threshold"
		parameter.parameter = "0.80"
		parameter.type_data = Parameter.PARAMETER_float
		parameter.software = software
		parameter.union_char = ":"
		parameter.can_change = True
		parameter.sequence_out = 1
		parameter.range_available = "[0.10:1.0]"
		parameter.range_max = "1.0"
		parameter.range_min = "0.10"
		parameter.description = "MINFRAC: minimum proportion for variant evidence (–minfrac) Range: [0.1:1.0]"
		vect_parameters.append(parameter)
		return vect_parameters

	def _get_medaka_model_default(self, user):
		"""
		Model for medaka
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_Medaka_name_consensus
		software.name_extended = SoftwareNames.SOFTWARE_Medaka_name_extended_consensus
		software.type_of_use = Software.TYPE_OF_USE_global
		software.type_of_software = Software.TYPE_SOFTWARE
		software.version = SoftwareNames.SOFTWARE_Medaka_VERSION
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "-m"
		parameter.parameter = SoftwareNames.SOFTWARE_Medaka_default_model	## default model
		parameter.type_data = Parameter.PARAMETER_char_list
		parameter.software = software
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
	
	def _get_samtools_depth_default_ONT(self, user):
		"""
		samtools depth for ONT
		"""
		software = Software()
		software.name = SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT
		software.name_extended = SoftwareNames.SOFTWARE_SAMTOOLS_name_extended_depth_ONT
		software.type_of_use = Software.TYPE_OF_USE_global
		software.type_of_software = Software.TYPE_SOFTWARE
		software.version = SoftwareNames.SOFTWARE_SAMTOOLS_VERSION
		software.owner = user
		
		vect_parameters =  []
		
		parameter = Parameter()
		parameter.name = "-q"
		parameter.parameter = "0"	## default model
		parameter.type_data = Parameter.PARAMETER_int
		parameter.software = software
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
		parameter.union_char = ""
		parameter.can_change = False
		parameter.sequence_out = 3
		parameter.description = "Output absolutely all positions, including unused ref. sequences."
		vect_parameters.append(parameter)
		
		return vect_parameters
	
