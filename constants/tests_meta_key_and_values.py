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
		self.assertEqual("TaskID_Project 20", self.metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_Queue_TaskID_Project, 20))
		self.assertEqual("Elements 120", self.metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_Elements_Reference, 120))
		self.assertEqual("Elements And CDS 120", self.metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_Elements_And_CDS_Reference, 120))
		self.assertEqual("Alert Coverage >9 230", self.metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_ALERT_COVERAGE_9, 230))


