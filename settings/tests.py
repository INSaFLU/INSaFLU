from django.test import TestCase
from settings.models import Software, Parameter
from constants.software_names import SoftwareNames
from settings.default_software import DefaultSoftware
from constants.constantsTestsCase import ConstantsTestsCase
from django.contrib.auth.models import User

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
		self.assertEqual(1, len(vect_software))
		self.assertEqual(SoftwareNames.SOFTWARE_TRIMMOMATIC_name, vect_software[0])
		
		### test all defaults
		default_software.test_all_defaults(user)
		
		self.assertEqual("SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33", default_software.get_trimmomatic_parameters(user))
		
		try:
			software = Software.objects.get(name=SoftwareNames.SOFTWARE_TRIMMOMATIC_name, owner=user)
		except Software.DoesNotExist:
			self.assertFail("Must exist this software name")
			
		parameters = Parameter.objects.filter(software=software)
		self.assertTrue(9, len(parameters))
		
		### test set default
		self.assertEqual("5", parameters[0].parameter)
		parameter = parameters[0]
		parameter.parameter = "42334"
		parameter.save()
		
		parameters = Parameter.objects.filter(software=software)
		self.assertEqual("42334", parameters[0].parameter)
		
		default_software.set_default_software(software, user)
		parameters = Parameter.objects.filter(software=software)
		self.assertEqual("5", parameters[0].parameter)
		self.assertNotEqual("42334", parameters[0].parameter)
		


