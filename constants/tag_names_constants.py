'''
Created on Nov 30, 2017

@author: mmp
'''

class TagNamesConstants(object):
	'''
	OBSOLETE, not used for now
	'''
	
	MINIMUM_NUMBER_TO_COUNT_STATISTICS = 110	## has the minimum number of instances in table 'CountVariations' to calculate the statistics
	TAG_PERCENTIL_VAR_TOTAL = 'Total Variation'
	TAG_PERCENTIL_VAR_50 = 'Variation <50'
	TAG_PERCENTIL_VAR_50_90 = '50<Variation<90'
	TAG_PERCENTIL_PREFIX = 'Percentile' 
	TAG_PERCENTIL_99 = '99'
	TAG_PERCENTIL_98 = '98'
	TAG_PERCENTIL_95 = '95'
	TAG_PERCENTIL_90 = '90'
	TAG_PERCENTIL_85 = '85'
	TAG_PERCENTIL_80 = '80'
	
	###### This is the values that are used to get the alert...
	TAG_PERCENTIL_INSAFLU = TAG_PERCENTIL_99
	TAG_PERCENTIL_VAR_INSAFLU = TAG_PERCENTIL_VAR_TOTAL
	
	### set the empirical values for percentils, is the number in database are less than 100
	dt_empirical = { TAG_PERCENTIL_99 : 35, TAG_PERCENTIL_98 : 33, TAG_PERCENTIL_95 : 30,\
				TAG_PERCENTIL_90 : 20, TAG_PERCENTIL_85 : 17, TAG_PERCENTIL_80 : 15, }
	
	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	def is_meta_tag_name(self, tag_name):
		"""
		return true if it's a meta tag name
		"""
		if (tag_name.startswith(self.TAG_PERCENTIL_PREFIX)): return True
		return False

	
	def get_all_tags_percentil(self):
		"""
		type_percentil : TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL, TagNamesConstants.TAG_PERCENTIL_VAR_50, TagNamesConstants.TAG_PERCENTIL_VAR_50_90 
		return all tags percentil
		"""
		vect_return = []
		vect_return.extend(self.__get_tags_percentile__(TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL))
		vect_return.extend(self.__get_tags_percentile__(TagNamesConstants.TAG_PERCENTIL_VAR_50))
		vect_return.extend(self.__get_tags_percentile__(TagNamesConstants.TAG_PERCENTIL_VAR_50_90))
		return vect_return

	def get_percentil_tag_name(self, percentile, type_percentile):
		"""
		in percentile: TagNamesConstants.TAG_PERCENTIL_98, TagNamesConstants.TAG_PERCENTIL_95, ..
		in type_percentil: TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL, TagNamesConstants.TAG_PERCENTIL_VAR_50, TagNamesConstants.TAG_PERCENTIL_VAR_50_90
		
		return, name of percentile
		"""
		return '{} {} {}'.format(TagNamesConstants.TAG_PERCENTIL_PREFIX, percentile, type_percentile)


	def get_number_percentil_from_tag(self, percentile_name):
		"""
		in: 'percentile 95 Variation <50', and ohers variations
		return: float(95)
		"""
		lst_data = percentile_name.split(' ')
		for data_ in lst_data:
			if (self.is_float(data_)): return float(data_)
		return None

	def is_which_var(self, percentile_name, type_var):
		"""
		return if total var
		"""
		return percentile_name.endswith(type_var)


	def get_emperical_percentil_values(self, value_percentile):
		"""
		float value_percentile, 80, 95, 99, etc
		"""
		value_percentile_ = value_percentile
		if (isinstance(value_percentile_, float)): value_percentile_ = int(value_percentile_)
		if (isinstance(value_percentile_, int)): value_percentile_ = str(value_percentile_)
		if (value_percentile_ in self.dt_empirical): return self.dt_empirical[value_percentile_]
		return None


	def __get_tags_percentile__(self, type_percentil):
		"""
		in type_percentile: TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL, TagNamesConstants.TAG_PERCENTIL_VAR_50, TagNamesConstants.TAG_PERCENTIL_VAR_50_90
		"""
		vect_return = []
		vect_return.append(self.get_percentil_tag_name(self.TAG_PERCENTIL_99, type_percentil))
		vect_return.append(self.get_percentil_tag_name(self.TAG_PERCENTIL_98, type_percentil))
		vect_return.append(self.get_percentil_tag_name(self.TAG_PERCENTIL_95, type_percentil))
		vect_return.append(self.get_percentil_tag_name(self.TAG_PERCENTIL_90, type_percentil))
		vect_return.append(self.get_percentil_tag_name(self.TAG_PERCENTIL_85, type_percentil))
		vect_return.append(self.get_percentil_tag_name(self.TAG_PERCENTIL_80, type_percentil))
		return vect_return

	
	def is_float(self, n_value):
		"""
		it also in utils but we can call utils here 
		"""
		try:
			float(n_value)
			return True
		except ValueError: 
			return False

