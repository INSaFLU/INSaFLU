from django.test import TestCase
from settings.models import Software, Parameter
from constants.software_names import SoftwareNames, Constants
from settings.default_software import DefaultSoftware
from settings.default_software_project_sample import DefaultProjectSoftware
from settings.default_parameters import DefaultParameters
from settings.constants_settings import ConstantsSettings
from constants.constantsTestsCase import ConstantsTestsCase
from django.contrib.auth.models import User
from managing_files.models import Project, Sample, ProjectSample

# Create your tests here.

class testsDefaultSoftwares(TestCase):
	
	def setUp(self):
		pass
	
	def test_software_obsolete(self):
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME + " set obsolete")
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME + " set obsolete"
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()
		
		default_software = DefaultSoftware()
		vect_software = default_software.get_all_software()
		self.assertEqual(10, len(vect_software))
		self.assertEqual(SoftwareNames.SOFTWARE_TRIMMOMATIC_name, vect_software[0])
		self.assertEqual(SoftwareNames.SOFTWARE_SNIPPY_name, vect_software[1])
		self.assertEqual(SoftwareNames.SOFTWARE_FREEBAYES_name, vect_software[2])
		self.assertEqual(SoftwareNames.SOFTWARE_NanoFilt_name, vect_software[3])
		self.assertEqual(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, vect_software[4])
		
		### test all defaults
		default_software.test_all_defaults(user)
		
		self.assertEqual("SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33",
					default_software.get_trimmomatic_parameters(user))
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51", default_software.get_snippy_parameters(user))
		self.assertEqual("-q 10 -l 50 --headcrop 70 --tailcrop 70", default_software.get_nanofilt_parameters(user))
		self.assertEqual("--minid 70 --mincov 30", default_software.get_abricate_parameters(user, ConstantsSettings.TECHNOLOGY_illumina))
		self.assertEqual("--minid 70 --mincov 30", default_software.get_abricate_parameters(user, ConstantsSettings.TECHNOLOGY_minion))
		
		### set new software with version 0
		software = Software()
		software.name = SoftwareNames.SOFTWARE_TRIMMOMATIC_name
		software.name_extended = SoftwareNames.SOFTWARE_TRIMMOMATIC_name_extended
		software.version = SoftwareNames.SOFTWARE_TRIMMOMATIC_VERSION
		software.type_of_use = Software.TYPE_OF_USE_global
		software.type_of_software = Software.TYPE_SOFTWARE
		software.version_parameters = 0
		software.owner = user
		software.save()
		
		###
		query_set = Software.objects.filter(name=SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
					type_of_use = Software.TYPE_OF_USE_global,
					is_obsolete=False, owner=user)
		self.assertEqual(2, len(query_set))
		
		default_parameters = DefaultParameters()
		default_parameters.set_software_obsolete()
		
		###
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_TRIMMOMATIC_name, owner=user,\
						type_of_use = Software.TYPE_OF_USE_global,
						version_parameters = default_parameters.get_software_parameters_version(SoftwareNames.SOFTWARE_TRIMMOMATIC_name))
			self.assertEqual(default_parameters.get_software_parameters_version(SoftwareNames.SOFTWARE_TRIMMOMATIC_name),
							software.version_parameters)
			self.assertFalse(software.is_obsolete)
		except Software.DoesNotExist:
			self.fail("Must exist")
			
		###
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_TRIMMOMATIC_name, owner=user,\
						type_of_use = Software.TYPE_OF_USE_global,
						version_parameters = 0)
			self.assertEqual(0, software.version_parameters)
			self.assertTrue(software.is_obsolete)
		except Software.DoesNotExist:
			self.fail("Must exist")
						
	def test_default_software(self):
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()
		
		default_software = DefaultSoftware()
		vect_software = default_software.get_all_software()
		self.assertEqual(10, len(vect_software))
		self.assertEqual(SoftwareNames.SOFTWARE_TRIMMOMATIC_name, vect_software[0])
		self.assertEqual(SoftwareNames.SOFTWARE_SNIPPY_name, vect_software[1])
		self.assertEqual(SoftwareNames.SOFTWARE_FREEBAYES_name, vect_software[2])
		self.assertEqual(SoftwareNames.SOFTWARE_NanoFilt_name, vect_software[3])
		self.assertEqual(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, vect_software[4])
		
		### test all defaults
		default_software.test_all_defaults(user)
		
		self.assertEqual("SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33",
					default_software.get_trimmomatic_parameters(user))
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51", default_software.get_snippy_parameters(user))
		self.assertEqual("-q 10 -l 50 --headcrop 70 --tailcrop 70", default_software.get_nanofilt_parameters(user))
		self.assertEqual("-q 10 -l 50 --headcrop 70 --tailcrop 70", 
				default_software.get_parameters(SoftwareNames.SOFTWARE_NanoFilt_name, user))
		self.assertEqual("Threshold:70", default_software.get_mask_consensus_threshold_parameters(user, ConstantsSettings.TECHNOLOGY_illumina))
		self.assertEqual("Threshold:70", default_software.get_mask_consensus_threshold_parameters(user, ConstantsSettings.TECHNOLOGY_minion))
		self.assertTrue(default_software.is_software_to_run(SoftwareNames.SOFTWARE_NanoFilt_name,
							user, ConstantsSettings.TECHNOLOGY_minion))
		is_to_run = False
		default_software.set_software_to_run(SoftwareNames.SOFTWARE_NanoFilt_name,
							user, ConstantsSettings.TECHNOLOGY_minion, is_to_run)
		self.assertFalse(default_software.is_software_to_run(SoftwareNames.SOFTWARE_NanoFilt_name,
							user, ConstantsSettings.TECHNOLOGY_minion))
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_TRIMMOMATIC_name, owner=user,
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = ConstantsSettings.TECHNOLOGY_illumina)
			self.assertTrue(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(9, len(parameters))
		
		### test set default
		self.assertEqual("5", parameters[3].parameter)
		self.assertEqual("0", parameters[2].parameter)
		parameter = parameters[2]
		parameter.parameter = "42334"
		parameter.save()
		
		parameters = Parameter.objects.filter(software=software)
		self.assertEqual("42334", parameters[2].parameter)
		
		parameter = parameters[0]
		parameter.parameter = SoftwareNames.SOFTWARE_TRIMMOMATIC_addapter_vect_available[1]
		parameter.save()
		self.assertEqual("ILLUMINACLIP:/usr/local/software/insaflu/trimmomatic/adapters/All-adapters.fa:3:30:10:6:true CROP:42334 " +
					"SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33",
					default_software.get_trimmomatic_parameters(user))
			
		default_software.set_default_software(software)
		parameters = Parameter.objects.filter(software=software)
		self.assertEqual("5", parameters[3].parameter)
		self.assertNotEqual("42334", parameters[2].parameter)
		
		#####
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_SNIPPY_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = ConstantsSettings.TECHNOLOGY_illumina)
			self.assertTrue(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(3, len(parameters))
		
		### test set default
		self.assertEqual("20", parameters[0].parameter)
		self.assertEqual("10", parameters[1].parameter)
		self.assertEqual("0.51", parameters[2].parameter)
		parameter = parameters[2]
		parameter.parameter = "42334"
		parameter.save()
		
		parameters = Parameter.objects.filter(software=software)
		self.assertEqual("20", parameters[0].parameter)
		self.assertEqual("10", parameters[1].parameter)
		self.assertEqual("42334", parameters[2].parameter)
		
		default_software.set_default_software(software)
		parameters = Parameter.objects.filter(software=software)
		self.assertEqual("20", parameters[0].parameter)
		self.assertEqual("10", parameters[1].parameter)
		self.assertEqual("0.51", parameters[2].parameter)
		self.assertNotEqual("42334", parameters[2].parameter)

		self.assertTrue(default_software.is_software_to_run(SoftwareNames.SOFTWARE_SNIPPY_name,
							user, ConstantsSettings.TECHNOLOGY_illumina))
		## this software can not be ON/OFF
		is_to_run = False
		self.assertFalse(default_software.set_software_to_run(SoftwareNames.SOFTWARE_SNIPPY_name,
							user, ConstantsSettings.TECHNOLOGY_illumina, is_to_run))
		self.assertTrue(default_software.is_software_to_run(SoftwareNames.SOFTWARE_SNIPPY_name,
							user, ConstantsSettings.TECHNOLOGY_illumina))
		
		##########################
		### test nanofilt
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_NanoFilt_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = ConstantsSettings.TECHNOLOGY_illumina)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_NanoFilt_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = ConstantsSettings.TECHNOLOGY_minion)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertFalse(software.is_used_in_project())
			self.assertFalse(software.is_used_in_sample())
			self.assertTrue(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		
		### test set default
		self.assertEqual("10", parameters[0].parameter)
		self.assertEqual("50", parameters[1].parameter)
		self.assertEqual("-l", parameters[1].name)
		parameter = parameters[0]
		parameter.parameter = "20"
		parameter.save()
		self.assertEqual("-q 20 -l 50 --headcrop 70 --tailcrop 70", default_software.get_nanofilt_parameters(user))
		
		self.assertEqual("--maxlength", parameters[4].name)
		parameter = parameters[4]
		parameter.parameter = "20"
		parameter.save()
		self.assertEqual("-q 20 -l 50 --headcrop 70 --tailcrop 70 --maxlength 20", default_software.get_nanofilt_parameters(user))
		self.assertEqual("-q 20 -l 50 --headcrop 70 --tailcrop 70 --maxlength 20", 
					default_software.get_parameters(SoftwareNames.SOFTWARE_NanoFilt_name,
							user, ConstantsSettings.TECHNOLOGY_minion))
		
		##########################
		### test medaka
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_Medaka_name_consensus, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = ConstantsSettings.TECHNOLOGY_illumina)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_Medaka_name_consensus, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = ConstantsSettings.TECHNOLOGY_minion)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertFalse(software.is_used_in_project())
			self.assertFalse(software.is_used_in_sample())
			self.assertTrue(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		
		### test set default
		self.assertEqual("r941_min_high_g360", parameters[0].parameter)
		self.assertEqual("-m", parameters[0].name)
		self.assertEqual("-m r941_min_high_g360", default_software.get_medaka_parameters_consensus(user))
		parameter = parameters[0]
		parameter.parameter = "20"
		parameter.save()
		self.assertEqual("-m 20", default_software.get_medaka_parameters_consensus(user))
		self.assertEqual("-m 20", default_software.get_parameters(
			SoftwareNames.SOFTWARE_Medaka_name_consensus,
			user, ConstantsSettings.TECHNOLOGY_minion))
		
		##########################
		### test limit_coverage_ONT_parameters
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = ConstantsSettings.TECHNOLOGY_illumina)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = ConstantsSettings.TECHNOLOGY_minion)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertFalse(software.is_used_in_project())
			self.assertFalse(software.is_used_in_sample())
			self.assertTrue(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))

		### test set default
		self.assertEqual("30", parameters[0].parameter)
		self.assertEqual("Threshold", parameters[0].name)
		self.assertEqual("Threshold:30", default_software.get_limit_coverage_ONT_parameters(user))
		parameter = parameters[0]
		parameter.parameter = "20"
		parameter.save()
		self.assertEqual("Threshold:20", default_software.get_limit_coverage_ONT_parameters(user))
		self.assertEqual("Threshold:20", default_software.get_parameters(SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,
													user, ConstantsSettings.TECHNOLOGY_minion))
		
		
		##########################
		### test freq_vcf ONT_parameters
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = ConstantsSettings.TECHNOLOGY_illumina)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = ConstantsSettings.TECHNOLOGY_minion)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertFalse(software.is_used_in_project())
			self.assertFalse(software.is_used_in_sample())
			self.assertTrue(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))

		### test set default
		self.assertEqual("0.80", parameters[0].parameter)
		self.assertEqual("Threshold", parameters[0].name)
		self.assertEqual("Threshold:0.80", default_software.get_vcf_freq_ONT_parameters(user))
		parameter = parameters[0]
		parameter.parameter = "0.1"
		parameter.save()
		self.assertEqual("Threshold:0.1", default_software.get_vcf_freq_ONT_parameters(user))
		self.assertEqual("Threshold:0.1", default_software.get_parameters(SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name,
													user, ConstantsSettings.TECHNOLOGY_minion))
		
		##########################
		### test Samtools  ONT
# try:
# 	software = Software.objects.get(name=SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT, owner=user,\
# 					type_of_use = Software.TYPE_OF_USE_global,
# 					technology__name = ConstantsSettings.TECHNOLOGY_illumina)
# 	self.fail("must fail")
# except Software.DoesNotExist:
# 	pass
#
# try:
# 	software = Software.objects.get(name=SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT, owner=user,\
# 					type_of_use = Software.TYPE_OF_USE_global,
# 					technology__name = ConstantsSettings.TECHNOLOGY_minion)
# 	self.assertFalse(software.is_used_in_project_sample())
# 	self.assertFalse(software.is_used_in_project())
# 	self.assertFalse(software.is_used_in_sample())
# 	self.assertTrue(software.is_used_in_global())
# except Software.DoesNotExist:
# 	self.fail("must exist this software name")
#
# parameters = Parameter.objects.filter(software=software)
# self.assertTrue(1, len(parameters))
#
# ### test set default
# self.assertEqual("0", parameters[0].parameter)
# self.assertEqual("-q", parameters[0].name)
# self.assertEqual("-aa", default_software.get_samtools_parameters_depth_ONT(user))
# parameter = parameters[0]
# parameter.parameter = "20"
# parameter.save()
# parameter = parameters[1]
# parameter.parameter = "40"
# parameter.save()
# self.assertEqual("-q 20 -Q 40 -aa", default_software.get_samtools_parameters_depth_ONT(user))
#
# parameter = parameters[1]
# parameter.parameter = "0"
# parameter.save()
# self.assertEqual("-q 20 -aa", default_software.get_samtools_parameters_depth_ONT(user))
		
				### Test Abricate
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_ABRICATE_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = ConstantsSettings.TECHNOLOGY_illumina)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertFalse(software.is_used_in_project())
			self.assertFalse(software.is_used_in_sample())
			self.assertTrue(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")
			
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_ABRICATE_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = ConstantsSettings.TECHNOLOGY_minion)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertFalse(software.is_used_in_project())
			self.assertFalse(software.is_used_in_sample())
			self.assertTrue(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(2, len(parameters))
		
		### test set default
		self.assertEqual("70", parameters[0].parameter)
		self.assertEqual("--minid", parameters[0].name)
		self.assertEqual("30", parameters[1].parameter)
		self.assertEqual("--mincov", parameters[1].name)

		parameter = parameters[0]
		parameter.parameter = "60"
		parameter.save()
		parameter = parameters[1]
		parameter.parameter = "20"
		parameter.save()
		self.assertEqual("--minid 60 --mincov 20", default_software.get_abricate_parameters(user, ConstantsSettings.TECHNOLOGY_minion))
		self.assertEqual("--minid 70 --mincov 30", default_software.get_abricate_parameters(user, ConstantsSettings.TECHNOLOGY_illumina))


	def test_default_project_sample_software_1(self):
		""" Only test project_sample_software """
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()
			
		project_name = "file_name_3_4"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = project_name
			project.owner = user
			project.save()
		
		### define a project
		sample_name = "file_name_2_222"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.path_name_1.name = "/tmp/zpt"
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_ont)
			sample.owner = user
			sample.save()
			
		###  save project sample
		project_sample = ProjectSample()
		project_sample.project = project
		project_sample.sample = sample
		project_sample.save()
			
		#### save a project_sample
		default_software = DefaultProjectSoftware()
		default_software.test_all_defaults(user, project, None, None)
		default_software.test_all_defaults(user, None, None, sample)
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51",\
				default_software.get_snippy_parameters_all_possibilities(user, project_sample))
		self.assertEqual("Threshold:70",\
				default_software.get_mask_consensus_parameters_all_possibilities(user, project_sample, ConstantsSettings.TECHNOLOGY_illumina))
		self.assertEqual("Threshold:70",\
				default_software.get_mask_consensus_parameters_all_possibilities(user, project_sample, ConstantsSettings.TECHNOLOGY_minion))

		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project_sample,
							technology__name = ConstantsSettings.TECHNOLOGY_minion)
			self.fail("Must not exist this software name")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project,
							technology__name = ConstantsSettings.TECHNOLOGY_minion)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertTrue(software.is_used_in_project())
			self.assertFalse(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")
			
		### change for project maskConsensus
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project,
							technology__name = ConstantsSettings.TECHNOLOGY_illumina)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertTrue(software.is_used_in_project())
			self.assertFalse(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		
		### test set default
		self.assertEqual("70", parameters[0].parameter)
		parameter = parameters[0]
		parameter.parameter = "20"
		parameter.save()
		self.assertEqual("Threshold:70",\
				default_software.get_mask_consensus_parameters_all_possibilities(user, project_sample, ConstantsSettings.TECHNOLOGY_minion))
		self.assertEqual("Threshold:20",\
				default_software.get_mask_consensus_parameters_all_possibilities(user, project_sample, ConstantsSettings.TECHNOLOGY_illumina))
		self.assertTrue(default_software.can_change_values_for_this_software(software, project, None, None))
		
		##########################
		### test nanofilt
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_NanoFilt_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_sample,
							technology__name = ConstantsSettings.TECHNOLOGY_illumina)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_NanoFilt_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_sample,
							technology__name = ConstantsSettings.TECHNOLOGY_minion)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertFalse(software.is_used_in_project())
			self.assertFalse(software.is_used_in_global())
			self.assertTrue(software.is_used_in_sample())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		
		### test set default
		self.assertEqual("10", parameters[0].parameter)
		self.assertEqual("50", parameters[1].parameter)
		self.assertEqual("-q", parameters[0].name)
		parameter = parameters[0]
		parameter.parameter = "20"
		parameter.save()
		self.assertEqual("-q 20 -l 50 --headcrop 70 --tailcrop 70", default_software.get_nanofilt_parameters(user,
									Software.TYPE_OF_USE_sample, sample))
		self.assertEqual("-q 20 -l 50 --headcrop 70 --tailcrop 70", default_software.get_nanofilt_parameters_all_possibilities(user,
									sample))
		self.assertEqual("20", default_software.get_nanofilt_single_parameter(sample,
									DefaultParameters.NANOfilt_quality_read))
		
		self.assertFalse(default_software.is_nanofilt_single_parameter_default_for_project(sample,
									DefaultParameters.NANOfilt_quality_read))
		self.assertTrue(default_software.can_change_values_for_this_software(software, None, None, sample))
		
		##########################
		### test medaka
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_Medaka_name_consensus, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = ConstantsSettings.TECHNOLOGY_illumina)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_Medaka_name_consensus, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project,
							technology__name = ConstantsSettings.TECHNOLOGY_minion)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertTrue(software.is_used_in_project())
			self.assertFalse(software.is_used_in_sample())
			self.assertFalse(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		
		### test set default
		self.assertEqual("r941_min_high_g360", parameters[0].parameter)
		self.assertEqual("-m", parameters[0].name)
		self.assertEqual("-m r941_min_high_g360", default_software.get_medaka_parameters(user,
											Software.TYPE_OF_USE_project, project, None))
		parameter = parameters[0]
		parameter.parameter = "20"
		parameter.save()
		self.assertEqual("-m 20", default_software.get_medaka_parameters(user,
									Software.TYPE_OF_USE_project, project, None))
		self.assertFalse(default_software.is_medaka_single_parameter_default_for_project(project,
											DefaultParameters.MEDAKA_model))

		##########################
		### test limit_coverage_ONT_parameters
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = ConstantsSettings.TECHNOLOGY_illumina)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project,
							technology__name = ConstantsSettings.TECHNOLOGY_minion)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertTrue(software.is_used_in_project())
			self.assertFalse(software.is_used_in_sample())
			self.assertFalse(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		
		### test set default
		self.assertEqual("30", parameters[0].parameter)
		self.assertEqual("Threshold", parameters[0].name)
		self.assertEqual("Threshold:30", default_software.get_limit_coverage_ONT_parameters(user,
											Software.TYPE_OF_USE_project, project, None))
		parameter = parameters[0]
		parameter.parameter = "20"
		parameter.save()
		self.assertEqual("Threshold:20", default_software.get_limit_coverage_ONT_parameters(user,
											Software.TYPE_OF_USE_project, project, None))
		self.assertFalse(default_software.is_limit_coverage_ONT_single_parameter_default_for_project(project,
											DefaultParameters.MASK_CONSENSUS_threshold))
		self.assertEqual("20", default_software.get_limit_coverage_ONT_single_parameter(project_sample,
											DefaultParameters.MASK_CONSENSUS_threshold))
		
		
		##########################
		### test vcf Freq ONT_parameters
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = ConstantsSettings.TECHNOLOGY_illumina)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project,
							technology__name = ConstantsSettings.TECHNOLOGY_minion)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertTrue(software.is_used_in_project())
			self.assertFalse(software.is_used_in_sample())
			self.assertFalse(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		
		### test set default
		self.assertEqual("0.80", parameters[0].parameter)
		self.assertEqual("Threshold", parameters[0].name)
		self.assertEqual("Threshold:0.80", default_software.get_freq_vcf_ONT_parameters(user,
										Software.TYPE_OF_USE_project, project, None))
		self.assertTrue(default_software.is_freq_vcf_ONT_single_parameter_default_for_project(project,
										DefaultParameters.MASK_CONSENSUS_threshold))

		parameter = parameters[0]
		parameter.parameter = "0.12"
		parameter.save()
		self.assertEqual("Threshold:0.12", default_software.get_freq_vcf_ONT_parameters(user,
											Software.TYPE_OF_USE_project, project, None))
		self.assertFalse(default_software.is_freq_vcf_ONT_single_parameter_default_for_project(project,
											DefaultParameters.MASK_CONSENSUS_threshold))
		self.assertEqual("0.12", default_software.get_freq_vcf_ONT_single_parameter(project_sample,
											DefaultParameters.MASK_CONSENSUS_threshold))
		self.assertEqual(ConstantsSettings.PIPELINE_NAME_variant_detection, software.pipeline_step.name)
		
		##########################
		### test samtools ONT
# try:
# 	software = Software.objects.get(name=SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT, owner=user,\
# 					type_of_use = Software.TYPE_OF_USE_global,
# 					technology__name = ConstantsSettings.TECHNOLOGY_illumina)
# 	self.fail("Must fail")
# except Software.DoesNotExist:
# 	pass
#
# try:
# 	software = Software.objects.get(name=SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT, owner=user,\
# 					type_of_use = Software.TYPE_OF_USE_project,
# 					technology__name = ConstantsSettings.TECHNOLOGY_minion)
# 	self.assertFalse(software.is_used_in_project_sample())
# 	self.assertTrue(software.is_used_in_project())
# 	self.assertFalse(software.is_used_in_sample())
# 	self.assertFalse(software.is_used_in_global())
# except Software.DoesNotExist:
# 	self.fail("Must exist this software name")
#
# parameters = Parameter.objects.filter(software=software)
# self.assertTrue(1, len(parameters))
# self.assertFalse(default_software.can_change_values_for_this_software(software, None, None, None))
#
# ### test set default
# self.assertEqual("0", parameters[0].parameter)
# self.assertEqual("-q", parameters[0].name)
# self.assertEqual(False, parameters[0].can_change)
# self.assertEqual(None, default_software.get_samtools_single_parameter_ONT(
# 									project_sample, "-q"))
# self.assertEqual(None, default_software.get_samtools_single_parameter_ONT(
# 									project_sample, "-Q"))
# parameter = parameters[0]
# parameter.parameter = "20"
# parameter.save()
# parameter = parameters[1]
# parameter.parameter = "30"
# parameter.save()
#
# self.assertEqual("-q 20 -Q 30 -aa", default_software.get_samtools_parameters_ONT(user,
# 								Software.TYPE_OF_USE_project, project, None))
# self.assertFalse(default_software.is_samtools_single_parameter_default_for_project_ONT(
# 			project, "-q"))
# self.assertEqual("20", default_software.get_samtools_single_parameter_ONT(
# 			project_sample, "-q"))


	def test_default_software_2(self):
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		###################		
		default_software = DefaultSoftware()
		default_software.test_all_defaults(user)
		
		#####
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_SNIPPY_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global)
			self.assertTrue(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(3, len(parameters))
		
		### test set default
		self.assertEqual("20", parameters[0].parameter)
		self.assertEqual("10", parameters[1].parameter)
		self.assertEqual("0.51", parameters[2].parameter)
		parameter = parameters[2]
		parameter.parameter = "42334"
		parameter.save()
		
		parameters = Parameter.objects.filter(software=software)
		self.assertEqual("20", parameters[0].parameter)
		self.assertEqual("10", parameters[1].parameter)
		self.assertEqual("42334", parameters[2].parameter)
		
		project_name = "file_name_3_4"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = project_name
			project.owner = user
			project.save()
		
		### define a project
		sample_name = "file_name_2_222"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.path_name_1.name = "/tmp/zpt"
			sample.owner = user
			sample.save()
			
		###  save project sample
		project_sample = ProjectSample()
		project_sample.project = project
		project_sample.sample = sample
		project_sample.save()
			
		#### save a project_sample
		default_software = DefaultProjectSoftware()
		default_software.test_all_defaults(user, project, None, None)
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 42334", default_software.get_snippy_parameters_all_possibilities(user, project_sample))
		self.assertEqual(ConstantsSettings.PIPELINE_NAME_variant_detection, software.pipeline_step.name)

	def test_default_software_project(self):
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()
		
		### define a project
		sample_name = "file_name_2_232"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.path_name_1.name = "/tmp/zpt"
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina)
			sample.owner = user
			sample.save()
			
		project_name = "file_name_3"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = project_name
			project.owner = user
			project.save()
			
		default_software = DefaultProjectSoftware()
		vect_software = default_software.get_all_software()
		self.assertEqual(10, len(vect_software))
		self.assertEqual(SoftwareNames.SOFTWARE_SNIPPY_name, vect_software[0])
	#	self.assertEqual(SoftwareNames.SOFTWARE_FREEBAYES_name, vect_software[1])
		
		### test all defaults
		default_software.test_all_defaults(user, project, None, None)
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51", default_software.get_snippy_parameters(user, Software.TYPE_OF_USE_project, project, None))
		self.assertEqual("--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10 " +\
				"--min-alternate-fraction 0.01 --ploidy 2 -V", default_software.get_freebayes_parameters(user, Software.TYPE_OF_USE_project, project, None))
		self.assertEqual("--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10 " +\
				"--min-alternate-fraction 0.01 --ploidy 2 -V", SoftwareNames.SOFTWARE_FREEBAYES_PARAMETERS)
		self.assertEqual("Threshold:70", default_software.get_mask_consensus_parameters_for_project(user, project, ConstantsSettings.TECHNOLOGY_illumina))
		self.assertEqual("Threshold:70", default_software.get_mask_consensus_parameters_for_project(user, project, ConstantsSettings.TECHNOLOGY_minion))

		
		try:
			software = Software.objects.filter(name=SoftwareNames.SOFTWARE_SNIPPY_name, owner=user, type_of_use=Software.TYPE_OF_USE_project,\
							parameter__project=project).distinct()
			self.assertEqual(1, software.count())
			self.assertFalse(software[0].is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist Snippy software for this project")

		parameters = Parameter.objects.filter(software=software[0], project=project)
		self.assertTrue(3, len(parameters))
		
		### test set default
		self.assertEqual("0.51", parameters[2].parameter)
		parameter = parameters[2]
		parameter.parameter = "0.222"
		parameter.save()
		
		parameters = Parameter.objects.filter(software=software, project=project)
		self.assertEqual("0.222", parameters[2].parameter)
		
		default_software.set_default_software(software[0], user, Software.TYPE_OF_USE_project, project, None, None)
		parameters = Parameter.objects.filter(software=software[0], project=project)
		self.assertEqual("0.51", parameters[2].parameter)
		self.assertNotEqual("0.222", parameters[2].parameter)

		parameter = parameters[2]
		parameter.parameter = "0.222"
		parameter.save()
		parameters = Parameter.objects.filter(software=software, project=project)
		self.assertEqual("0.222", parameters[2].parameter)
		
		#####################################
		####################################
		##### new project
		project_name = "file_name_4erww"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = project_name
			project.owner = user
			project.save()
		
		###  save project sample
		project_sample = ProjectSample()
		project_sample.project = project
		project_sample.sample = sample
		project_sample.save()
		
		## default
		default_software = DefaultProjectSoftware()
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51", default_software.get_snippy_parameters_all_possibilities(user, project_sample))
		self.assertEqual("Threshold:70", default_software.get_mask_consensus_parameters_all_possibilities(user, project_sample, ConstantsSettings.TECHNOLOGY_illumina))
		self.assertEqual(True, default_software.is_mask_consensus_single_parameter_default(project_sample,
				DefaultParameters.MASK_CONSENSUS_threshold, ConstantsSettings.TECHNOLOGY_illumina))

		default_software = DefaultProjectSoftware()
		
		### test all defaults
		default_software.test_all_defaults(user, project, None, None)
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51", default_software.get_snippy_parameters(user, Software.TYPE_OF_USE_project, project, None))
		self.assertEqual("--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10 " +\
				"--min-alternate-fraction 0.01 --ploidy 2 -V", default_software.get_freebayes_parameters(user, Software.TYPE_OF_USE_project, project, None))
		
		try:
			software = Software.objects.filter(name=SoftwareNames.SOFTWARE_SNIPPY_name, owner=user, type_of_use=Software.TYPE_OF_USE_project,\
							parameter__project=project, technology__name=ConstantsSettings.TECHNOLOGY_illumina).distinct()
			self.assertEqual(1, software.count())
			self.assertFalse(software[0].is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist Snippy software for this project")

		parameters = Parameter.objects.filter(software=software[0], project=project)
		self.assertTrue(3, len(parameters))
		
		### test set default
		self.assertEqual("0.51", parameters[2].parameter)
		parameter = parameters[1]
		parameter.parameter = "30"
		parameter.save()
		parameter = parameters[2]
		parameter.parameter = "0.222"
		parameter.save()
		
		#### 
		self.assertEqual("--mapqual 20 --mincov 30 --minfrac 0.222", default_software.get_snippy_parameters_all_possibilities(user, project_sample))
		
		self.assertFalse(default_software.is_snippy_single_parameter_default(project_sample, DefaultParameters.SNIPPY_COVERAGE_NAME))
		self.assertTrue(default_software.is_snippy_single_parameter_default(project_sample, DefaultParameters.SNIPPY_MAPQUAL_NAME))
		self.assertEqual('30', default_software.get_snippy_single_parameter(project_sample, DefaultParameters.SNIPPY_COVERAGE_NAME)) 
		self.assertEqual('20', default_software.get_snippy_single_parameter(project_sample, DefaultParameters.SNIPPY_MAPQUAL_NAME)) 
		self.assertEqual('20', default_software.get_snippy_single_parameter_for_project(project, DefaultParameters.SNIPPY_MAPQUAL_NAME)) 

		
		### must pass
		default_software.test_all_defaults(user, project, None, None)
		try:
			software = Software.objects.filter(name=SoftwareNames.SOFTWARE_SNIPPY_name, owner=user, type_of_use=Software.TYPE_OF_USE_project,\
							parameter__project=project).distinct()
			self.assertEqual(1, software.count())
			self.assertFalse(software[0].is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist Snippy software for this project")
			
		parameters = Parameter.objects.filter(software=software, project=project)
		self.assertTrue(3, len(parameters))
		self.assertEqual("0.222", parameters[2].parameter)
		
		#### save a project_sample
		default_software.test_all_defaults(user, None, project_sample, None)
		self.assertEqual("--mapqual 20 --mincov 30 --minfrac 0.222", default_software.get_snippy_parameters_all_possibilities(user, project_sample))
		
		default_software.set_default_software(software[0], user, Software.TYPE_OF_USE_project, project, None, None)
		parameters = Parameter.objects.filter(software=software[0], project=project)
		self.assertTrue(3, len(parameters))
		self.assertEqual("0.51", parameters[2].parameter)
		self.assertNotEqual("0.222", parameters[2].parameter)
		self.assertTrue(default_software.is_change_values_for_software(software[0].name, ConstantsSettings.TECHNOLOGY_illumina,
												ConstantsTestsCase.TEST_USER_NAME))
		
		default_software.set_default_software(software[0], user, Software.TYPE_OF_USE_project, project, None, None)
		self.assertFalse(default_software.is_change_values_for_software(software[0].name, ConstantsSettings.TECHNOLOGY_illumina,
												ConstantsTestsCase.TEST_USER_NAME))

		#### save a project_sample
		self.assertEqual("--mapqual 20 --mincov 30 --minfrac 0.222", default_software.get_snippy_parameters_all_possibilities(user, project_sample))
		
		### test extra mask maker
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project, technology__name = ConstantsSettings.TECHNOLOGY_illumina)
			self.assertTrue(software.is_used_in_project())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		self.assertEqual(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name_extended, software.name_extended)
		self.assertEqual(Software.TYPE_INSAFLU_PARAMETER, software.type_of_software)
		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		
		self.assertEqual("70", parameters[0].parameter)
		parameter = parameters[0]
		parameter.parameter = "80"
		parameter.save()
		
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			self.fail("Must exist project name")
		self.assertEqual("Threshold:80", default_software.get_mask_consensus_parameters_for_project(user, project,
					ConstantsSettings.TECHNOLOGY_illumina))
		self.assertEqual("80", default_software.get_mask_consensus_single_parameter_for_project(project,\
				DefaultParameters.MASK_CONSENSUS_threshold,
				ConstantsSettings.TECHNOLOGY_illumina))
		self.assertEqual("Threshold:70", default_software.get_mask_consensus_parameters_for_project(user, project,
					ConstantsSettings.TECHNOLOGY_minion))
		self.assertEqual("70", default_software.get_mask_consensus_single_parameter_for_project(project,\
				DefaultParameters.MASK_CONSENSUS_threshold,
				ConstantsSettings.TECHNOLOGY_minion))
		
		project_name = "file_name_3_4ww"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = project_name
			project.owner = user
			project.save()
		
		##############################
		### test mask site consensus
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_MASK_CONSENSUS_BY_SITE_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project,
							technology__name = ConstantsSettings.TECHNOLOGY_generic)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertTrue(software.is_used_in_project())
			self.assertFalse(software.is_used_in_global())
			self.assertFalse(software.is_used_in_sample())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")
		self.assertEqual(ConstantsSettings.PIPELINE_NAME_variant_detection, software.pipeline_step.name)
		self.assertEqual(False, software.can_be_on_off_in_pipeline)
		
		######################################################################
		######################################################################
		###  save project sample
		project_sample = ProjectSample()
		project_sample.project = project
		project_sample.sample = sample
		project_sample.save()

		default_software = DefaultProjectSoftware()
		default_software.test_all_defaults(user, None, project_sample, None)
		self.assertEqual("Threshold:70", default_software.get_mask_consensus_parameters_all_possibilities(user,
				project_sample, ConstantsSettings.TECHNOLOGY_illumina))
		self.assertEqual("70", default_software.get_mask_consensus_single_parameter(project_sample,\
				DefaultParameters.MASK_CONSENSUS_threshold, ConstantsSettings.TECHNOLOGY_illumina))
		
		### for project sample
		default_software.test_all_defaults(user, None, project_sample, None)
		### test extra mask maker
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project_sample,
							technology__name=ConstantsSettings.TECHNOLOGY_illumina)
			self.assertTrue(software.is_used_in_project_sample())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")
		
		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		self.assertEqual("70", parameters[0].parameter)
		parameter = parameters[0]
		parameter.parameter = "90"
		parameter.save()
		
		self.assertEqual("90", default_software.get_mask_consensus_single_parameter(project_sample,\
				DefaultParameters.MASK_CONSENSUS_threshold, ConstantsSettings.TECHNOLOGY_illumina))
		self.assertEqual("Threshold:90", default_software.get_mask_consensus_parameters_all_possibilities(user,
				project_sample, ConstantsSettings.TECHNOLOGY_illumina))
		
		##############################
		### test mask site consensus
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_MASK_CONSENSUS_BY_SITE_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project_sample,
							technology__name = ConstantsSettings.TECHNOLOGY_generic)
			self.assertTrue(software.is_used_in_project_sample())
			self.assertFalse(software.is_used_in_project())
			self.assertFalse(software.is_used_in_global())
			self.assertFalse(software.is_used_in_sample())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")
		self.assertEqual(ConstantsSettings.PIPELINE_NAME_variant_detection, software.pipeline_step.name)
		self.assertEqual(False, software.can_be_on_off_in_pipeline)
		
		##########################
		### test Trimmomatic
		default_software.test_all_defaults(user, None, None, sample)
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_TRIMMOMATIC_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_sample,
							technology__name = ConstantsSettings.TECHNOLOGY_minion)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_TRIMMOMATIC_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_sample,
							technology__name = ConstantsSettings.TECHNOLOGY_illumina)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertFalse(software.is_used_in_project())
			self.assertFalse(software.is_used_in_global())
			self.assertTrue(software.is_used_in_sample())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		self.assertEqual(ConstantsSettings.PIPELINE_NAME_read_quality_analysis, software.pipeline_step.name)
		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		
		### test set default
		self.assertEqual("Not apply", parameters[0].parameter)
		self.assertEqual("0", parameters[1].parameter)
		self.assertEqual("5", parameters[3].parameter)
		self.assertEqual("HEADCROP", parameters[1].name)
		self.assertEqual("SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33",
						default_software.get_trimmomatic_parameters(user,
						Software.TYPE_OF_USE_sample, sample))
		parameter = parameters[1]
		parameter.parameter = "20"
		parameter.save()
		self.assertEqual("HEADCROP:20 SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33",
						default_software.get_trimmomatic_parameters(user,
						Software.TYPE_OF_USE_sample, sample))
		
		self.assertFalse(default_software.is_trimmomatic_single_parameter_default_for_project(sample,
									"HEADCROP"))
		
		##############################
		### test abricate
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_ABRICATE_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_sample,
							technology__name = ConstantsSettings.TECHNOLOGY_illumina)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertFalse(software.is_used_in_project())
			self.assertFalse(software.is_used_in_global())
			self.assertTrue(software.is_used_in_sample())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")
			
		self.assertEqual(None, default_software.get_abricate_parameters(user,
				Software.TYPE_OF_USE_sample, sample, ConstantsSettings.TECHNOLOGY_minion))
		self.assertEqual("--minid 70 --mincov 30", default_software.get_abricate_parameters(user,
				Software.TYPE_OF_USE_sample, sample, ConstantsSettings.TECHNOLOGY_illumina))
		self.assertTrue(default_software.is_abricate_single_parameter_default_for_project(sample,
									"--minid", ConstantsSettings.TECHNOLOGY_illumina))
		
		self.assertEqual(ConstantsSettings.PIPELINE_NAME_type_and_subtype_analysis, software.pipeline_step.name)
		

	def test_default_software_project_sample(self):
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()
		
		### define a project
		sample_name = "file_name_2"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.path_name_1.name = "/tmp/zpt"
			sample.set_type_of_fastq_sequencing(Constants.FORMAT_FASTQ_illumina)
			sample.owner = user
			sample.save()
			
		project_name = "file_name_3"
		try:
			project_ = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project_ = Project()
			project_.name = project_name
			project_.owner = user
			project_.save()

		###  save project sample
		project_sample = ProjectSample()
		project_sample.project = project_
		project_sample.sample = sample
		project_sample.save()
		
		## default
		default_software = DefaultProjectSoftware()
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51", default_software.get_snippy_parameters_all_possibilities(user, project_sample))
		
		vect_software = default_software.get_all_software()
		self.assertEqual(10, len(vect_software))
		self.assertEqual(SoftwareNames.SOFTWARE_SNIPPY_name, vect_software[0])
		self.assertEqual(SoftwareNames.SOFTWARE_Medaka_name_consensus, vect_software[1])
		
		### test all defaults
		default_software.test_all_defaults(user, None, project_sample, None)
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51", default_software.get_snippy_parameters(user, Software.TYPE_OF_USE_project_sample, None, project_sample))
		self.assertEqual("--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10 " +\
				"--min-alternate-fraction 0.01 --ploidy 2 -V", default_software.get_freebayes_parameters(user, Software.TYPE_OF_USE_project_sample, None, project_sample))
		
		try:
			software = Software.objects.filter(name=SoftwareNames.SOFTWARE_SNIPPY_name, owner=user, type_of_use=Software.TYPE_OF_USE_project_sample,\
							parameter__project_sample=project_sample).distinct()
			self.assertEqual(1, software.count())
			self.assertFalse(software[0].is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist Snippy software for this project")

		parameters = Parameter.objects.filter(software=software[0], project_sample=project_sample)
		self.assertTrue(3, len(parameters))
		
		### test set default
		self.assertEqual("0.51", parameters[2].parameter)
		parameter = parameters[2]
		parameter.parameter = "0.222"
		parameter.save()
		
		parameters = Parameter.objects.filter(software=software, project_sample=project_sample)
		self.assertEqual("0.222", parameters[2].parameter)
		
		default_software.set_default_software(software[0], user, Software.TYPE_OF_USE_project_sample, None, project_sample, None)
		parameters = Parameter.objects.filter(software=software[0], project_sample=project_sample)
		self.assertEqual("0.51", parameters[2].parameter)
		self.assertNotEqual("0.222", parameters[2].parameter)

		parameter = parameters[2]
		parameter.parameter = "0.222"
		parameter.save()
		parameters = Parameter.objects.filter(software=software, project_sample=project_sample)
		self.assertEqual("0.222", parameters[2].parameter)
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.222", default_software.get_snippy_parameters_all_possibilities(user, project_sample))
		
		self.assertTrue(default_software.is_snippy_single_parameter_default(project_sample, DefaultParameters.SNIPPY_COVERAGE_NAME))
		self.assertTrue(default_software.is_snippy_single_parameter_default(project_sample, DefaultParameters.SNIPPY_MAPQUAL_NAME))
		self.assertEqual('10', default_software.get_snippy_single_parameter(project_sample, DefaultParameters.SNIPPY_COVERAGE_NAME)) 
		self.assertEqual('20', default_software.get_snippy_single_parameter(project_sample, DefaultParameters.SNIPPY_MAPQUAL_NAME)) 
		
		#####################################
		####################################
		##### new project
		project_name = "file_name_4erww"
		try:
			project_ = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project_ = Project()
			project_.name = project_name
			project_.owner = user
			project_.save()

		###  save project sample
		project_sample = ProjectSample()
		project_sample.project = project_
		project_sample.sample = sample
		project_sample.save()
		
		default_software = DefaultProjectSoftware()
		vect_software = default_software.get_all_software()
		self.assertEqual(10, len(vect_software))
		self.assertEqual(SoftwareNames.SOFTWARE_SNIPPY_name, vect_software[0])
		self.assertEqual(SoftwareNames.SOFTWARE_Medaka_name_consensus, vect_software[1])
		
		### test all defaults
		default_software.test_all_defaults(user, None, project_sample, None)
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51", default_software.get_snippy_parameters(user, Software.TYPE_OF_USE_project_sample, None, project_sample))
		self.assertEqual("--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10 " +\
				"--min-alternate-fraction 0.01 --ploidy 2 -V", default_software.get_freebayes_parameters(user, Software.TYPE_OF_USE_project_sample, None, project_sample))
		
		try:
			software = Software.objects.filter(name=SoftwareNames.SOFTWARE_SNIPPY_name, owner=user, type_of_use=Software.TYPE_OF_USE_project_sample,\
							parameter__project_sample=project_sample).distinct()
			self.assertEqual(1, software.count())
			self.assertFalse(software[0].is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist Snippy software for this project")

		parameters = Parameter.objects.filter(software=software[0], project_sample=project_sample)
		self.assertTrue(3, len(parameters))
		
		### test set default
		self.assertEqual("0.51", parameters[2].parameter)
		parameter = parameters[2]
		parameter.parameter = "0.222"
		parameter.save()
		
		### must pass
		default_software.test_all_defaults(user, None, project_sample, None)
		try:
			software = Software.objects.filter(name=SoftwareNames.SOFTWARE_SNIPPY_name, owner=user, type_of_use=Software.TYPE_OF_USE_project_sample,\
							parameter__project_sample=project_sample).distinct()
			self.assertEqual(1, software.count())
			self.assertFalse(software[0].is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist Snippy software for this project")
			
		parameters = Parameter.objects.filter(software=software, project_sample=project_sample)
		self.assertTrue(3, len(parameters))
		self.assertEqual("0.222", parameters[2].parameter)
		
		default_software.set_default_software(software[0], user, Software.TYPE_OF_USE_project_sample, None, project_sample, None)
		parameters = Parameter.objects.filter(software=software[0], project_sample=project_sample)
		self.assertTrue(3, len(parameters))
		self.assertEqual("0.51", parameters[2].parameter)
		self.assertNotEqual("0.222", parameters[2].parameter)
		
		self.assertTrue(default_software.is_change_values_for_software(software[0].name, ConstantsSettings.TECHNOLOGY_illumina,
														ConstantsTestsCase.TEST_USER_NAME))
		default_software.set_default_software(software[0], user, Software.TYPE_OF_USE_project, None, project_sample, None)
		self.assertFalse(default_software.is_change_values_for_software(software[0].name, ConstantsSettings.TECHNOLOGY_illumina,
														ConstantsTestsCase.TEST_USER_NAME))
		
		software_names = SoftwareNames()
		self.assertFalse(default_software.is_change_values_for_software(software_names.get_freebayes_name(),
												ConstantsSettings.TECHNOLOGY_illumina, ConstantsTestsCase.TEST_USER_NAME))


