'''
Created on Nov 1, 2017

@author: mmp
'''
from .models import MetaKeySample, MetaKey

class ManageDatabase(object):
	'''
	classdocs
	'''


	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	def set_metakey(self, sample, owner, meta_key_name, value, description):
		"""
		save a meta key
		"""
		try:
			metaKey = MetaKey.objects.get(name=meta_key_name)
		except MetaKey.DoesNotExist:
			metaKey = MetaKey()
			metaKey.name = meta_key_name
			metaKey.save()
		
		metaKeySample = MetaKeySample()
		metaKeySample.sample = sample
		metaKeySample.meta_tag = metaKey
		metaKeySample.owner = owner
		metaKeySample.value = value
		metaKeySample.description = description
		metaKeySample.save()
		return metaKeySample
	
	def get_metakey(self, sample, meta_key_name, value):
		"""
		value = None, return a list
		"""
		try:
			if (value == None): return MetaKeySample.objects.filter(sample__id=sample.id, meta_tag__name=meta_key_name)
			return MetaKeySample.objects.get(sample__id=sample.id, meta_tag__name=meta_key_name, value=value)
		except MetaKeySample.DoesNotExist:
			return None