'''
Created on Dec 14, 2017

@author: mmp
'''

class ConstantsMixedInfection(object):
	'''
	classdocs
	'''
	## [<50, 50<value<90], 
	vect_start_compare = [[78, 56],[77, 56],[76, 57],[76, 53],[55, 51],[50, 36],[49, 35],[26,11]]

	TAGS_MIXED_INFECTION_YES = 'Yes'
	TAGS_MIXED_INFECTION_MAY_BE = 'May be'
	TAGS_MIXED_INFECTION_NO = 'No'
	
	## values to upload to database
	vect_upload_to_database = [TAGS_MIXED_INFECTION_YES, TAGS_MIXED_INFECTION_MAY_BE, TAGS_MIXED_INFECTION_NO]
	
	### all other values are NO
	threshold_yes = 0.98	##	>=
	threshold_may_be = 0.97	##	>=
	
	def get_tag_by_value(self, value):
		"""
		get tab by value
		"""
		if (value >= self.threshold_yes): return self.TAGS_MIXED_INFECTION_YES
		if (value >= self.threshold_may_be): return self.TAGS_MIXED_INFECTION_MAY_BE
		return self.TAGS_MIXED_INFECTION_NO
	
	def is_alert(self, tag):
		"""
		return true if is necessary generate an alert
		"""
		if (tag == self.TAGS_MIXED_INFECTION_YES or tag == self.TAGS_MIXED_INFECTION_MAY_BE): return True
		return False