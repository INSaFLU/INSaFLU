'''
Created on Jan 2, 2018

@author: mmp
'''
import unittest, requests, os
from django.conf import settings 
from constants.constantsTestsCase import ConstantsTestsCase
from django.utils.safestring import mark_safe

class Test(unittest.TestCase):


	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		self.baseDirectoryUrl = os.path.join(getattr(settings, "MEDIA_URL", None), ConstantsTestsCase.MANAGING_TESTS)


	def tearDown(self):
		pass


# 	def test_ranges(self):
# 		url = mark_safe('http://127.0.0.1:8000/media/uploads/references/userId_1/refId_1/ref_H3.fasta')
# 		headers = {"Range": "bytes=0-100"}
# 		r = requests.get(url, headers=headers)
# 		self.assertEquals(101, len(str(r.text)))

