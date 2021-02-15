from django.test import TestCase
from settings.models import Software, Parameter
from constants.software_names import SoftwareNames, Constants
from settings.default_software import DefaultSoftware
from settings.default_software_project_sample import DefaultProjectSoftware
from constants.constantsTestsCase import ConstantsTestsCase
from django.contrib.auth.models import User
from managing_files.models import Project, Sample, ProjectSample
# Create your tests here.

class testsDefaultSoftwares(TestCase):
	
	def setUp(self):
		pass
		
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
		self.assertEqual(7, len(vect_software))
		self.assertEqual(SoftwareNames.SOFTWARE_TRIMMOMATIC_name, vect_software[0])
		self.assertEqual(SoftwareNames.SOFTWARE_SNIPPY_name, vect_software[1])
		self.assertEqual(SoftwareNames.SOFTWARE_NanoFilt_name, vect_software[2])
		self.assertEqual(SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, vect_software[3])
		
		### test all defaults
		default_software.test_all_defaults(user)
		
		self.assertEqual("SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33", default_software.get_trimmomatic_parameters(user))
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51", default_software.get_snippy_parameters(user))
		self.assertEqual("-q 10 -l 50 --headcrop 10 --tailcrop 10", default_software.get_nanofilt_parameters(user))
		self.assertEqual("-q 10 -l 50 --headcrop 10 --tailcrop 10", 
				default_software.get_parameters(SoftwareNames.SOFTWARE_NanoFilt_name, user))
		self.assertEqual("Threshold:70", default_software.get_mask_consensus_threshold_parameters(user, SoftwareNames.TECHNOLOGY_illumina))
		self.assertEqual("Threshold:70", default_software.get_mask_consensus_threshold_parameters(user, SoftwareNames.TECHNOLOGY_minion))
		
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_TRIMMOMATIC_name, owner=user,
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = SoftwareNames.TECHNOLOGY_illumina)
			self.assertTrue(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(9, len(parameters))
		
		### test set default
		self.assertEqual("5", parameters[2].parameter)
		parameter = parameters[2]
		parameter.parameter = "42334"
		parameter.save()
		
		parameters = Parameter.objects.filter(software=software)
		self.assertEqual("42334", parameters[2].parameter)
		
		default_software.set_default_software(software, user)
		parameters = Parameter.objects.filter(software=software)
		self.assertEqual("5", parameters[2].parameter)
		self.assertNotEqual("42334", parameters[2].parameter)
		
		#####
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_SNIPPY_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = SoftwareNames.TECHNOLOGY_illumina)
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
		
		default_software.set_default_software(software, user)
		parameters = Parameter.objects.filter(software=software)
		self.assertEqual("20", parameters[0].parameter)
		self.assertEqual("10", parameters[1].parameter)
		self.assertEqual("0.51", parameters[2].parameter)
		self.assertNotEqual("42334", parameters[2].parameter)

		##########################
		### test nanofilt
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_NanoFilt_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = SoftwareNames.TECHNOLOGY_illumina)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_NanoFilt_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = SoftwareNames.TECHNOLOGY_minion)
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
		self.assertEqual("-q 20 -l 50 --headcrop 10 --tailcrop 10", default_software.get_nanofilt_parameters(user))
		
		self.assertEqual("--maxlength", parameters[4].name)
		parameter = parameters[4]
		parameter.parameter = "20"
		parameter.save()
		self.assertEqual("-q 20 -l 50 --headcrop 10 --tailcrop 10 --maxlength 20", default_software.get_nanofilt_parameters(user))
		self.assertEqual("-q 20 -l 50 --headcrop 10 --tailcrop 10 --maxlength 20", 
					default_software.get_parameters(SoftwareNames.SOFTWARE_NanoFilt_name,
							user, SoftwareNames.TECHNOLOGY_minion))
		
		##########################
		### test medaka
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_Medaka_name_consensus, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = SoftwareNames.TECHNOLOGY_illumina)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_Medaka_name_consensus, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = SoftwareNames.TECHNOLOGY_minion)
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
			user, SoftwareNames.TECHNOLOGY_minion))
		
		##########################
		### test limit_coverage_ONT_parameters
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = SoftwareNames.TECHNOLOGY_illumina)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = SoftwareNames.TECHNOLOGY_minion)
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
													user, SoftwareNames.TECHNOLOGY_minion))
		
		##########################
		### test Samtools 
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = SoftwareNames.TECHNOLOGY_illumina)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = SoftwareNames.TECHNOLOGY_minion)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertFalse(software.is_used_in_project())
			self.assertFalse(software.is_used_in_sample())
			self.assertTrue(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		
		### test set default
		self.assertEqual("0", parameters[0].parameter)
		self.assertEqual("-q", parameters[0].name)
		self.assertEqual("-aa", default_software.get_samtools_parameters_depth_ONT(user))
		parameter = parameters[0]
		parameter.parameter = "20"
		parameter.save()
		parameter = parameters[1]
		parameter.parameter = "40"
		parameter.save()
		self.assertEqual("-q 20 -Q 40 -aa", default_software.get_samtools_parameters_depth_ONT(user))
		
		parameter = parameters[1]
		parameter.parameter = "0"
		parameter.save()
		self.assertEqual("-q 20 -aa", default_software.get_samtools_parameters_depth_ONT(user))
		
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
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project, project, None, None)
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_sample, None, None, sample)
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51",\
				default_software.get_snippy_parameters_all_possibilities(user, project_sample))
		self.assertEqual("Threshold:70",\
				default_software.get_mask_consensus_parameters_all_possibilities(user, project_sample, SoftwareNames.TECHNOLOGY_illumina))
		self.assertEqual("Threshold:70",\
				default_software.get_mask_consensus_parameters_all_possibilities(user, project_sample, SoftwareNames.TECHNOLOGY_minion))

		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project_sample,
							technology__name = SoftwareNames.TECHNOLOGY_minion)
			self.fail("Must not exist this software name")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project,
							technology__name = SoftwareNames.TECHNOLOGY_minion)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertTrue(software.is_used_in_project())
			self.assertFalse(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")
			
		### change for project maskConsensus
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project,
							technology__name = SoftwareNames.TECHNOLOGY_illumina)
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
				default_software.get_mask_consensus_parameters_all_possibilities(user, project_sample, SoftwareNames.TECHNOLOGY_minion))
		self.assertEqual("Threshold:20",\
				default_software.get_mask_consensus_parameters_all_possibilities(user, project_sample, SoftwareNames.TECHNOLOGY_illumina))

		##########################
		### test nanofilt
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_NanoFilt_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_sample,
							technology__name = SoftwareNames.TECHNOLOGY_illumina)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_NanoFilt_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_sample,
							technology__name = SoftwareNames.TECHNOLOGY_minion)
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
		self.assertEqual("-q 20 -l 50 --headcrop 10 --tailcrop 10", default_software.get_nanofilt_parameters(user,
									Software.TYPE_OF_USE_sample, sample))
		self.assertEqual("-q 20 -l 50 --headcrop 10 --tailcrop 10", default_software.get_nanofilt_parameters_all_possibilities(user,
									sample))
		self.assertEqual("20", default_software.get_nanofilt_single_parameter(sample,
									DefaultProjectSoftware.NANOfilt_quality_read))
		
		self.assertFalse(default_software.is_nanofilt_single_parameter_default_for_project(sample,
									DefaultProjectSoftware.NANOfilt_quality_read))
		
		##########################
		### test medaka
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_Medaka_name_consensus, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = SoftwareNames.TECHNOLOGY_illumina)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_Medaka_name_consensus, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project,
							technology__name = SoftwareNames.TECHNOLOGY_minion)
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
											DefaultProjectSoftware.MEDAKA_model))

		##########################
		### test limit_coverage_ONT_parameters
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = SoftwareNames.TECHNOLOGY_illumina)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project,
							technology__name = SoftwareNames.TECHNOLOGY_minion)
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
											DefaultProjectSoftware.MASK_CONSENSUS_threshold))
		self.assertEqual("20", default_software.get_limit_coverage_ONT_single_parameter(project_sample,
											DefaultProjectSoftware.MASK_CONSENSUS_threshold))
		
		
		##########################
		### test samtools ONT
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT, owner=user,\
							type_of_use = Software.TYPE_OF_USE_global,
							technology__name = SoftwareNames.TECHNOLOGY_illumina)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project,
							technology__name = SoftwareNames.TECHNOLOGY_minion)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertTrue(software.is_used_in_project())
			self.assertFalse(software.is_used_in_sample())
			self.assertFalse(software.is_used_in_global())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		
		### test set default
		self.assertEqual("0", parameters[0].parameter)
		self.assertEqual("-q", parameters[0].name)
		self.assertEqual(None, default_software.get_samtools_single_parameter_ONT(
											project_sample, "-q"))
		self.assertEqual(None, default_software.get_samtools_single_parameter_ONT(
											project_sample, "-Q"))
		parameter = parameters[0]
		parameter.parameter = "20"
		parameter.save()
		parameter = parameters[1]
		parameter.parameter = "30"
		parameter.save()
		
		self.assertEqual("-q 20 -Q 30 -aa", default_software.get_samtools_parameters_ONT(user,
										Software.TYPE_OF_USE_project, project, None))
		self.assertFalse(default_software.is_samtools_single_parameter_default_for_project_ONT(
					project, "-q"))
		self.assertEqual("20", default_software.get_samtools_single_parameter_ONT(
					project_sample, "-q"))



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
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project_sample, project, None, None)
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 42334", default_software.get_snippy_parameters_all_possibilities(user, project_sample))


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
		self.assertEqual(7, len(vect_software))
		self.assertEqual(SoftwareNames.SOFTWARE_SNIPPY_name, vect_software[0])
	#	self.assertEqual(SoftwareNames.SOFTWARE_FREEBAYES_name, vect_software[1])
		
		### test all defaults
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project, project, None, None)
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51", default_software.get_snippy_parameters(user, Software.TYPE_OF_USE_project, project, None))
#		self.assertEqual("--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10 " +\
#				"--min-alternate-fraction 100 --ploidy 2 -V", default_software.get_freebayes_parameters(user, project))
		self.assertEqual("Threshold:70", default_software.get_mask_consensus_parameters_for_project(user, project, SoftwareNames.TECHNOLOGY_illumina))
		self.assertEqual("Threshold:70", default_software.get_mask_consensus_parameters_for_project(user, project, SoftwareNames.TECHNOLOGY_minion))

		
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
		self.assertEqual("Threshold:70", default_software.get_mask_consensus_parameters_all_possibilities(user, project_sample, SoftwareNames.TECHNOLOGY_illumina))
		self.assertEqual(True, default_software.is_mask_consensus_single_parameter_default(project_sample,
				DefaultProjectSoftware.MASK_CONSENSUS_threshold, SoftwareNames.TECHNOLOGY_illumina))

		default_software = DefaultProjectSoftware()
		
		### test all defaults
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project, project, None, None)
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51", default_software.get_snippy_parameters(user, Software.TYPE_OF_USE_project, project, None))
# 		self.assertEqual("--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10 " +\
# 				"--min-alternate-fraction 100 --ploidy 2 -V", default_software.get_freebayes_parameters(user, project))
		
		try:
			software = Software.objects.filter(name=SoftwareNames.SOFTWARE_SNIPPY_name, owner=user, type_of_use=Software.TYPE_OF_USE_project,\
							parameter__project=project, technology__name=SoftwareNames.TECHNOLOGY_illumina).distinct()
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
		
		self.assertFalse(default_software.is_snippy_single_parameter_default(project_sample, DefaultProjectSoftware.SNIPPY_COVERAGE_NAME))
		self.assertTrue(default_software.is_snippy_single_parameter_default(project_sample, DefaultProjectSoftware.SNIPPY_MAPQUAL_NAME))
		self.assertEqual('30', default_software.get_snippy_single_parameter(project_sample, DefaultProjectSoftware.SNIPPY_COVERAGE_NAME)) 
		self.assertEqual('20', default_software.get_snippy_single_parameter(project_sample, DefaultProjectSoftware.SNIPPY_MAPQUAL_NAME)) 
		self.assertEqual('20', default_software.get_snippy_single_parameter_for_project(project, DefaultProjectSoftware.SNIPPY_MAPQUAL_NAME)) 

		
		### must pass
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project, project, None, None)
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
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project_sample, None, project_sample, None)
		self.assertEqual("--mapqual 20 --mincov 30 --minfrac 0.222", default_software.get_snippy_parameters_all_possibilities(user, project_sample))
		
		default_software.set_default_software(software[0], user, Software.TYPE_OF_USE_project, project, None, None)
		parameters = Parameter.objects.filter(software=software[0], project=project)
		self.assertTrue(3, len(parameters))
		self.assertEqual("0.51", parameters[2].parameter)
		self.assertNotEqual("0.222", parameters[2].parameter)
		self.assertTrue(default_software.is_change_values_for_software(software[0].name, SoftwareNames.TECHNOLOGY_illumina))
		
		default_software.set_default_software(software[0], user, Software.TYPE_OF_USE_project, project, None, None)
		self.assertFalse(default_software.is_change_values_for_software(software[0].name, SoftwareNames.TECHNOLOGY_illumina))

		#### save a project_sample
		self.assertEqual("--mapqual 20 --mincov 30 --minfrac 0.222", default_software.get_snippy_parameters_all_possibilities(user, project_sample))
		
		### test extra mask maker
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project, technology__name = SoftwareNames.TECHNOLOGY_illumina)
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
					SoftwareNames.TECHNOLOGY_illumina))
		self.assertEqual("80", default_software.get_mask_consensus_single_parameter_for_project(project,\
				DefaultProjectSoftware.MASK_CONSENSUS_threshold,
				SoftwareNames.TECHNOLOGY_illumina))
		self.assertEqual("Threshold:70", default_software.get_mask_consensus_parameters_for_project(user, project,
					SoftwareNames.TECHNOLOGY_minion))
		self.assertEqual("70", default_software.get_mask_consensus_single_parameter_for_project(project,\
				DefaultProjectSoftware.MASK_CONSENSUS_threshold,
				SoftwareNames.TECHNOLOGY_minion))
		
		project_name = "file_name_3_4ww"
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

		default_software = DefaultProjectSoftware()
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project, project, None, None)
		self.assertEqual("Threshold:70", default_software.get_mask_consensus_parameters_all_possibilities(user,
				project_sample, SoftwareNames.TECHNOLOGY_illumina))
		self.assertEqual("70", default_software.get_mask_consensus_single_parameter(project_sample,\
				DefaultProjectSoftware.MASK_CONSENSUS_threshold, SoftwareNames.TECHNOLOGY_illumina))
		
		### for project sample
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project_sample, None, project_sample, None)
		### test extra mask maker
		try:
			software = Software.objects.get(name=SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_project_sample,
							technology__name=SoftwareNames.TECHNOLOGY_illumina)
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
				DefaultProjectSoftware.MASK_CONSENSUS_threshold, SoftwareNames.TECHNOLOGY_illumina))
		self.assertEqual("Threshold:90", default_software.get_mask_consensus_parameters_all_possibilities(user,
				project_sample, SoftwareNames.TECHNOLOGY_illumina))
		
		
		##########################
		### test Trimmomatic
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_sample, None, None, sample)
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_TRIMMOMATIC_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_sample,
							technology__name = SoftwareNames.TECHNOLOGY_minion)
			self.fail("Must fail")
		except Software.DoesNotExist:
			pass
			
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_TRIMMOMATIC_name, owner=user,\
							type_of_use = Software.TYPE_OF_USE_sample,
							technology__name = SoftwareNames.TECHNOLOGY_illumina)
			self.assertFalse(software.is_used_in_project_sample())
			self.assertFalse(software.is_used_in_project())
			self.assertFalse(software.is_used_in_global())
			self.assertTrue(software.is_used_in_sample())
		except Software.DoesNotExist:
			self.fail("Must exist this software name")

		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(1, len(parameters))
		
		### test set default
		self.assertEqual("0", parameters[0].parameter)
		self.assertEqual("5", parameters[2].parameter)
		self.assertEqual("HEADCROP", parameters[0].name)
		self.assertEqual("SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33",
						default_software.get_trimmomatic_parameters(user,
						Software.TYPE_OF_USE_sample, sample))
		parameter = parameters[0]
		parameter.parameter = "20"
		parameter.save()
		self.assertEqual("HEADCROP:20 SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33",
						default_software.get_trimmomatic_parameters(user,
						Software.TYPE_OF_USE_sample, sample))
		
		self.assertFalse(default_software.is_trimmomatic_single_parameter_default_for_project(sample,
									"HEADCROP"))
		
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
		self.assertEqual(7, len(vect_software))
		self.assertEqual(SoftwareNames.SOFTWARE_SNIPPY_name, vect_software[0])
	#	self.assertEqual(SoftwareNames.SOFTWARE_FREEBAYES_name, vect_software[1])
		
		### test all defaults
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project_sample, None, project_sample, None)
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51", default_software.get_snippy_parameters(user, Software.TYPE_OF_USE_project_sample, None, project_sample))
#		self.assertEqual("--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10 " +\
#				"--min-alternate-fraction 100 --ploidy 2 -V", default_software.get_freebayes_parameters(user, project))
		
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
		
		self.assertTrue(default_software.is_snippy_single_parameter_default(project_sample, DefaultProjectSoftware.SNIPPY_COVERAGE_NAME))
		self.assertTrue(default_software.is_snippy_single_parameter_default(project_sample, DefaultProjectSoftware.SNIPPY_MAPQUAL_NAME))
		self.assertEqual('10', default_software.get_snippy_single_parameter(project_sample, DefaultProjectSoftware.SNIPPY_COVERAGE_NAME)) 
		self.assertEqual('20', default_software.get_snippy_single_parameter(project_sample, DefaultProjectSoftware.SNIPPY_MAPQUAL_NAME)) 
		
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
		self.assertEqual(7, len(vect_software))
		self.assertEqual(SoftwareNames.SOFTWARE_SNIPPY_name, vect_software[0])
#		self.assertEqual(SoftwareNames.SOFTWARE_FREEBAYES_name, vect_software[1])
		
		### test all defaults
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project_sample, None, project_sample, None)
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51", default_software.get_snippy_parameters(user, Software.TYPE_OF_USE_project_sample, None, project_sample))
# 		self.assertEqual("--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10 " +\
# 				"--min-alternate-fraction 100 --ploidy 2 -V", default_software.get_freebayes_parameters(user, project))
		
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
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project_sample, None, project_sample, None)
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
		
		self.assertTrue(default_software.is_change_values_for_software(software[0].name, SoftwareNames.TECHNOLOGY_illumina))
		default_software.set_default_software(software[0], user, Software.TYPE_OF_USE_project, None, project_sample, None)
		self.assertFalse(default_software.is_change_values_for_software(software[0].name, SoftwareNames.TECHNOLOGY_illumina))
		
		software_names = SoftwareNames()
		self.assertFalse(default_software.is_change_values_for_software(software_names.get_freebayes_name(), SoftwareNames.TECHNOLOGY_illumina))

