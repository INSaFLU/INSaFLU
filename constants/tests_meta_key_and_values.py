'''
Created on Oct 28, 2017

@author: mmp
'''
import unittest
from constants.meta_key_and_values import MetaKeyAndValue

class Test(unittest.TestCase):


	### static
	metaKeyAndValue = MetaKeyAndValue()

	def test_meta_keys(self):
		self.assertTrue("", self.metaKeyAndValue.get_meta_key_by_project_id(MetaKeyAndValue.META_KEY_Queue_TaskID_Project, 20))
		self.assertTrue("", self.metaKeyAndValue.get_meta_key_by_project_id(MetaKeyAndValue.META_KEY_Elements_Project, 120))
		self.assertTrue("", self.metaKeyAndValue.get_meta_key_by_element(MetaKeyAndValue.META_KEY_ALERT_COVERAGE_9, 230))


