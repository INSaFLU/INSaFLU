'''
Created on 03/05/2020

@author: mmp
'''
from settings.models import Software, Parameter
from constants.software_names import SoftwareNames
from utils.lock_atomic_transaction import LockedAtomicTransaction
from settings.default_parameters import DefaultParameters
from settings.constants_settings import ConstantsSettings

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
	
	software_names = SoftwareNames()
	
	def __init__(self):
		""" change values """
		self.default_parameters = DefaultParameters()
		self.change_values_software = {}	### the key is the name of the software

	def test_all_defaults(self, user):
		
		### test all defaults
		self.test_default_db(SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
			self.default_parameters.get_trimmomatic_default(user, Software.TYPE_OF_USE_global,
			ConstantsSettings.TECHNOLOGY_illumina), user)
		self.test_default_db(SoftwareNames.SOFTWARE_SNIPPY_name,
				self.default_parameters.get_snippy_default(user, Software.TYPE_OF_USE_global,
				ConstantsSettings.TECHNOLOGY_illumina), user)
		self.test_default_db(SoftwareNames.SOFTWARE_FREEBAYES_name,
				self.default_parameters.get_freebayes_default(user, Software.TYPE_OF_USE_global,
				ConstantsSettings.TECHNOLOGY_illumina), user)
		self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
				self.default_parameters.get_mask_consensus_threshold_default(user,
				Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_illumina), user)
		self.test_default_db(SoftwareNames.SOFTWARE_ABRICATE_name,
				self.default_parameters.get_abricate_default(user,
				Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_illumina), user)
# 		self.test_default_db(SoftwareNames.SOFTWARE_CLEAN_HUMAN_READS_name,
# 				self.default_parameters.get_clean_human_reads_default(user,
#				Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_illumina), user)

		## ONT software
		self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
				self.default_parameters.get_mask_consensus_threshold_default(user,
				Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_minion), user)
		self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,
				self.default_parameters.get_limit_coverage_ONT_threshold_default(user,
				Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_minion), user)
		self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name,
				self.default_parameters.get_vcf_freq_ONT_threshold_default(user,
				Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_minion), user)
		self.test_default_db(SoftwareNames.SOFTWARE_Medaka_name_consensus,
				self.default_parameters.get_medaka_model_default(user,
				Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_minion), user)
		
		self.test_default_db(SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT,
 				self.default_parameters.get_samtools_depth_default_ONT(user,
				Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_minion), user)
		
		self.test_default_db(SoftwareNames.SOFTWARE_NanoFilt_name,
				self.default_parameters.get_nanofilt_default(user,
				Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_minion), user)
		
		self.test_default_db(SoftwareNames.SOFTWARE_ABRICATE_name,
				self.default_parameters.get_abricate_default(user,
				Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_minion), user)
# 		self.test_default_db(SoftwareNames.SOFTWARE_CLEAN_HUMAN_READS_name,
# 				self.default_parameters.get_clean_human_reads_default(user,
#				Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_minion), user)

	def test_default_db(self, software_name, vect_parameters, user):
		"""
		test if exist, if not persist in database
		"""
		type_of_use = Software.TYPE_OF_USE_global
		## lock because more than one process can duplicate software names
		with LockedAtomicTransaction(Software), LockedAtomicTransaction(Parameter):
			try:
				Software.objects.get(name=software_name, owner=user,\
							type_of_use = type_of_use,
							technology__name = vect_parameters[0].software.technology.name,
							version_parameters = \
							self.default_parameters.get_software_parameters_version(software_name))
			except Software.DoesNotExist: ### if not exist save it
				self.default_parameters.persist_parameters(vect_parameters, type_of_use)

	def get_trimmomatic_parameters(self, user):
		result = self.default_parameters.get_parameters(SoftwareNames.SOFTWARE_TRIMMOMATIC_name, user,
					Software.TYPE_OF_USE_global, None, None, None, ConstantsSettings.TECHNOLOGY_illumina)
		return "" if result is None else result
	def get_snippy_parameters(self, user):
		result = self.default_parameters.get_parameters(SoftwareNames.SOFTWARE_SNIPPY_name, user,
					Software.TYPE_OF_USE_global, None, None, None, ConstantsSettings.TECHNOLOGY_illumina)
		return "" if result is None else result
	def get_freebayes_parameters(self, user):
		result = self.default_parameters.get_parameters(SoftwareNames.SOFTWARE_FREEBAYES_name, user,
					Software.TYPE_OF_USE_global, None, None, None, ConstantsSettings.TECHNOLOGY_illumina)
		return "" if result is None else result
	def get_nanofilt_parameters(self, user):
		result = self.default_parameters.get_parameters(SoftwareNames.SOFTWARE_NanoFilt_name, user,
					Software.TYPE_OF_USE_global, None, None, None, ConstantsSettings.TECHNOLOGY_minion)
		return "" if result is None else result
	def get_mask_consensus_threshold_parameters(self, user, technology_name):
		result = self.default_parameters.get_parameters(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, user,
					Software.TYPE_OF_USE_global, None, None, None, technology_name)
		return "" if result is None else result
	def get_clean_human_reads_parameters(self, user, technology_name):
		result = self.default_parameters.get_parameters(SoftwareNames.SOFTWARE_CLEAN_HUMAN_READS_name, user,
					Software.TYPE_OF_USE_global, None, None, None, technology_name)
		return "" if result is None else result
	def get_limit_coverage_ONT_parameters(self, user):
		result = self.default_parameters.get_parameters(SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name, user,
					Software.TYPE_OF_USE_global, None, None, None, ConstantsSettings.TECHNOLOGY_minion)
		return "" if result is None else result
	def get_vcf_freq_ONT_parameters(self, user):
		result = self.default_parameters.get_parameters(SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name, user,
					Software.TYPE_OF_USE_global, None, None, None, ConstantsSettings.TECHNOLOGY_minion)
		return "" if result is None else result
	def get_medaka_parameters_consensus(self, user):
		result = self.default_parameters.get_parameters(SoftwareNames.SOFTWARE_Medaka_name_consensus, user,
					Software.TYPE_OF_USE_global, None, None, None, ConstantsSettings.TECHNOLOGY_minion)
		return "" if result is None else result
	def get_samtools_parameters_depth_ONT(self, user):
		result = self.default_parameters.get_parameters(SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT, user,
					Software.TYPE_OF_USE_global, None, None, None, ConstantsSettings.TECHNOLOGY_minion)
		return "" if result is None else result
	def get_abricate_parameters(self, user, technology_name):
		result = self.default_parameters.get_parameters(SoftwareNames.SOFTWARE_ABRICATE_name, user,
					Software.TYPE_OF_USE_global, None, None, None, technology_name)
		return "" if result is None else result
	
	def set_default_software(self, software):
		""" Set a default  
		"""
		vect_parameters = self.default_parameters.get_vect_parameters(software)
		if vect_parameters is None: return
		
		key_value = "{}_{}".format(software.name, software.technology.name)
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
		""" Return if the software has a value changed """
		key_value = "{}_{}".format(software.name, software.technology.name)
		return self.change_values_software.get(key_value, False)

	def get_parameters(self, software_name, user, technology_name = ConstantsSettings.TECHNOLOGY_illumina):
		"""
		"""
		if (software_name == SoftwareNames.SOFTWARE_TRIMMOMATIC_name):
			self.test_default_db(SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
				self.default_parameters.get_trimmomatic_default(user, Software.TYPE_OF_USE_global,
				ConstantsSettings.TECHNOLOGY_illumina), user)
			return self.get_trimmomatic_parameters(user)
		if (software_name == SoftwareNames.SOFTWARE_SNIPPY_name):
			self.test_default_db(SoftwareNames.SOFTWARE_SNIPPY_name,
				self.default_parameters.get_snippy_default(user, Software.TYPE_OF_USE_global,
				ConstantsSettings.TECHNOLOGY_illumina), user)
			return self.get_snippy_parameters(user)
		if (software_name == SoftwareNames.SOFTWARE_FREEBAYES_name):
			self.test_default_db(SoftwareNames.SOFTWARE_FREEBAYES_name,
				self.default_parameters.get_freebayes_default(user, Software.TYPE_OF_USE_global,
				ConstantsSettings.TECHNOLOGY_illumina), user)
			return self.get_freebayes_parameters(user)
		if (software_name == SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name):
			self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
				self.default_parameters.get_mask_consensus_threshold_default(user, Software.TYPE_OF_USE_global,
				technology_name), user)
			return self.get_mask_consensus_threshold_parameters(user, technology_name)
		if (software_name == SoftwareNames.SOFTWARE_CLEAN_HUMAN_READS_name):
			self.test_default_db(SoftwareNames.SOFTWARE_CLEAN_HUMAN_READS_name,
				self.default_parameters.get_clean_human_reads_default(user, Software.TYPE_OF_USE_global,
				technology_name), user)
			return self.get_clean_human_reads_parameters(user, technology_name)
		if (software_name == SoftwareNames.SOFTWARE_NanoFilt_name):
			self.test_default_db(SoftwareNames.SOFTWARE_NanoFilt_name,
				self.default_parameters.get_nanofilt_default(user, Software.TYPE_OF_USE_global,
				ConstantsSettings.TECHNOLOGY_minion), user)
			return self.get_nanofilt_parameters(user)
		if (software_name == SoftwareNames.SOFTWARE_Medaka_name_consensus):
			self.test_default_db(SoftwareNames.SOFTWARE_Medaka_name_consensus,
				self.default_parameters.get_medaka_model_default(user, Software.TYPE_OF_USE_global,
				ConstantsSettings.TECHNOLOGY_minion), user)
			return self.get_medaka_parameters_consensus(user)
		if (software_name == SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT):
			self.test_default_db(software_name, self.default_parameters.get_samtools_depth_default_ONT(user, Software.TYPE_OF_USE_global,
				ConstantsSettings.TECHNOLOGY_minion), user)
			return self.get_samtools_parameters_depth_ONT(user)
		if (software_name == SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name):
			self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,
				self.default_parameters.get_limit_coverage_ONT_threshold_default(user, Software.TYPE_OF_USE_global,
				ConstantsSettings.TECHNOLOGY_minion), user)
			return self.get_limit_coverage_ONT_parameters(user)
		if (software_name == SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name):
			self.test_default_db(SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name,
				self.default_parameters.get_vcf_freq_ONT_threshold_default(user, Software.TYPE_OF_USE_global,
				ConstantsSettings.TECHNOLOGY_minion), user)
			return self.get_vcf_freq_ONT_parameters(user)
		if (software_name == SoftwareNames.SOFTWARE_ABRICATE_name):
			self.test_default_db(SoftwareNames.SOFTWARE_ABRICATE_name,
				self.default_parameters.get_abricate_default(user, Software.TYPE_OF_USE_global, technology_name), user)
			return self.get_abricate_parameters(user, technology_name)
		return ""
		
	def get_all_software(self):
		"""
		get all softwares available by this class
		"""
		vect_software = []
		vect_software.append(self.software_names.get_trimmomatic_name())
		vect_software.append(self.software_names.get_snippy_name())
		vect_software.append(self.software_names.get_freebayes_name())
		vect_software.append(self.software_names.get_NanoFilt_name())
		vect_software.append(self.software_names.get_insaflu_parameter_mask_consensus_name())
		vect_software.append(self.software_names.get_medaka_name_extended_consensus())
		vect_software.append(self.software_names.get_samtools_name_depth_ONT())
		vect_software.append(self.software_names.get_insaflu_parameter_limit_coverage_name())
		vect_software.append(self.software_names.get_insaflu_parameter_freq_vcf_name())
		vect_software.append(self.software_names.get_abricate_name())
		return vect_software

