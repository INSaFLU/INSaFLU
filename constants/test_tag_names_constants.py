'''
Created on Nov 30, 2017

@author: mmp
'''
import unittest
from constants.tag_names_constants import TagNamesConstants

class Test(unittest.TestCase):


	tagNamesConstants = TagNamesConstants()
	def test_is_meta_tag_name(self):
		
		self.assertTrue(self.tagNamesConstants.is_meta_tag_name(\
			self.tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_99, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL)))
		self.assertTrue(self.tagNamesConstants.is_meta_tag_name(\
			self.tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_80, TagNamesConstants.TAG_PERCENTIL_VAR_50_90)))

		self.assertEquals('{} {} {}'.format(TagNamesConstants.TAG_PERCENTIL_PREFIX, TagNamesConstants.TAG_PERCENTIL_98, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL),\
			self.tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_98, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL))


	def test_get_number_percentil_from_tag(self):
		
		percentil_name = self.tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_98, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL)
		self.assertEquals(98, self.tagNamesConstants.get_number_percentil_from_tag(percentil_name))
		self.assertEquals(None, self.tagNamesConstants.get_number_percentil_from_tag('percentil_name'))
		percentil_name = self.tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_95, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL)
		self.assertEquals(95, self.tagNamesConstants.get_number_percentil_from_tag(percentil_name))

	
	def test_is_which_var(self):
		
		percentil_name = self.tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_95, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL)
		self.assertEquals(True, self.tagNamesConstants.is_which_var(percentil_name, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL))
		self.assertEquals(False, self.tagNamesConstants.is_which_var(percentil_name, TagNamesConstants.TAG_PERCENTIL_VAR_50))
		percentil_name = self.tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_80, TagNamesConstants.TAG_PERCENTIL_VAR_50_90)
		self.assertEquals(True, self.tagNamesConstants.is_which_var(percentil_name, TagNamesConstants.TAG_PERCENTIL_VAR_50_90))
		self.assertEquals(False, self.tagNamesConstants.is_which_var(percentil_name, TagNamesConstants.TAG_PERCENTIL_VAR_50))


	def test_number_percentils_names(self):
		self.assertEquals(18, len(self.tagNamesConstants.get_all_tags_percentil()))
		self.assertEquals('Percentile 99 Total Variation', self.tagNamesConstants.get_all_tags_percentil()[0])
		self.assertEquals('Percentile 80 50<Variation<90', self.tagNamesConstants.get_all_tags_percentil()[-1])



