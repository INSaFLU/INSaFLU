from django.test import TestCase
from settings.models import Software, Parameter
from constants.software_names import SoftwareNames
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
		self.assertEqual(2, len(vect_software))
		self.assertEqual(SoftwareNames.SOFTWARE_TRIMMOMATIC_name, vect_software[0])
		self.assertEqual(SoftwareNames.SOFTWARE_SNIPPY_name, vect_software[1])
		
		### test all defaults
		default_software.test_all_defaults(user)
		
		self.assertEqual("SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33", default_software.get_trimmomatic_parameters(user))
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51", default_software.get_snippy_parameters(user))
		
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_TRIMMOMATIC_name, owner=user,
							type_of_use = Software.TYPE_OF_USE_global)
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
		
		default_software.set_default_software(software, user)
		parameters = Parameter.objects.filter(software=software)
		self.assertEqual("20", parameters[0].parameter)
		self.assertEqual("10", parameters[1].parameter)
		self.assertEqual("0.51", parameters[2].parameter)
		self.assertNotEqual("42334", parameters[2].parameter)

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
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project_sample, project, None)
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51", default_software.get_snippy_parameters_all_possibilities(user, project_sample))


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
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project_sample, project, None)
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
		sample_name = "file_name_2"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.path_name_1.name = "/tmp/zpt"
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
		self.assertEqual(1, len(vect_software))
		self.assertEqual(SoftwareNames.SOFTWARE_SNIPPY_name, vect_software[0])
	#	self.assertEqual(SoftwareNames.SOFTWARE_FREEBAYES_name, vect_software[1])
		
		### test all defaults
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project, project, None)
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51", default_software.get_snippy_parameters(user, Software.TYPE_OF_USE_project, project, None))
#		self.assertEqual("--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10 " +\
#				"--min-alternate-fraction 100 --ploidy 2 -V", default_software.get_freebayes_parameters(user, project))
		
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
		
		default_software.set_default_software(software[0], user, Software.TYPE_OF_USE_project, project, None)
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

		default_software = DefaultProjectSoftware()
		vect_software = default_software.get_all_software()
		self.assertEqual(1, len(vect_software))
		self.assertEqual(SoftwareNames.SOFTWARE_SNIPPY_name, vect_software[0])
#		self.assertEqual(SoftwareNames.SOFTWARE_FREEBAYES_name, vect_software[1])
		
		### test all defaults
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project, project, None)
		self.assertEqual("--mapqual 20 --mincov 10 --minfrac 0.51", default_software.get_snippy_parameters(user, Software.TYPE_OF_USE_project, project, None))
# 		self.assertEqual("--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10 " +\
# 				"--min-alternate-fraction 100 --ploidy 2 -V", default_software.get_freebayes_parameters(user, project))
		
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
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project, project, None)
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
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project_sample, None, project_sample)
		self.assertEqual("--mapqual 20 --mincov 30 --minfrac 0.222", default_software.get_snippy_parameters_all_possibilities(user, project_sample))
		
		default_software.set_default_software(software[0], user, Software.TYPE_OF_USE_project, project, None)
		parameters = Parameter.objects.filter(software=software[0], project=project)
		self.assertTrue(3, len(parameters))
		self.assertEqual("0.51", parameters[2].parameter)
		self.assertNotEqual("0.222", parameters[2].parameter)
		self.assertTrue(default_software.is_change_values_for_software(software[0].name))
		
		default_software.set_default_software(software[0], user, Software.TYPE_OF_USE_project, project, None)
		self.assertFalse(default_software.is_change_values_for_software(software[0].name))

		#### save a project_sample
		self.assertEqual("--mapqual 20 --mincov 30 --minfrac 0.222", default_software.get_snippy_parameters_all_possibilities(user, project_sample))
		
		
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
		self.assertEqual(1, len(vect_software))
		self.assertEqual(SoftwareNames.SOFTWARE_SNIPPY_name, vect_software[0])
	#	self.assertEqual(SoftwareNames.SOFTWARE_FREEBAYES_name, vect_software[1])
		
		### test all defaults
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project_sample, None, project_sample)
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
		
		default_software.set_default_software(software[0], user, Software.TYPE_OF_USE_project_sample, None, project_sample)
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
		self.assertEqual(1, len(vect_software))
		self.assertEqual(SoftwareNames.SOFTWARE_SNIPPY_name, vect_software[0])
#		self.assertEqual(SoftwareNames.SOFTWARE_FREEBAYES_name, vect_software[1])
		
		### test all defaults
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project_sample, None, project_sample)
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
		default_software.test_all_defaults(user, Software.TYPE_OF_USE_project_sample, None, project_sample)
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
		
		default_software.set_default_software(software[0], user, Software.TYPE_OF_USE_project_sample, None, project_sample)
		parameters = Parameter.objects.filter(software=software[0], project_sample=project_sample)
		self.assertTrue(3, len(parameters))
		self.assertEqual("0.51", parameters[2].parameter)
		self.assertNotEqual("0.222", parameters[2].parameter)
		
		self.assertTrue(default_software.is_change_values_for_software(software[0].name))
		default_software.set_default_software(software[0], user, Software.TYPE_OF_USE_project, None, project_sample)
		self.assertFalse(default_software.is_change_values_for_software(software[0].name))
		
		software_names = SoftwareNames()
		self.assertFalse(default_software.is_change_values_for_software(software_names.get_freebayes_name()))

