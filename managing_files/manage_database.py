'''
Created on Nov 1, 2017

@author: mmp
'''
from managing_files.models import MetaKeySample, MetaKey, MetaKeyProject, MetaKeyProjectSample

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
			if (value == None): return MetaKeySample.objects.filter(sample__id=sample.id, meta_tag__name=meta_key_name).order_by('-creation_date')
			return MetaKeySample.objects.get(sample__id=sample.id, meta_tag__name=meta_key_name, value=value)
		except MetaKeySample.DoesNotExist:
			return None
	
	def set_project_metakey(self, project, owner, meta_key_name, value, description):
		"""
		save a meta key
		"""
		try:
			metaKey = MetaKey.objects.get(name=meta_key_name)
		except MetaKey.DoesNotExist:
			metaKey = MetaKey()
			metaKey.name = meta_key_name
			metaKey.save()
		
		metaKeyProject = MetaKeyProject()
		metaKeyProject.project = project
		metaKeyProject.meta_tag = metaKey
		metaKeyProject.owner = owner
		metaKeyProject.value = value
		metaKeyProject.description = description
		metaKeyProject.save()
		return MetaKeyProject
	
	def get_project_metakey(self, project, meta_key_name, value):
		"""
		value = None, return a list
		"""
		try:
			if (value == None): return MetaKeyProject.objects.filter(project__id=project.id, meta_tag__name=meta_key_name).order_by('-creation_date')
			return MetaKeyProject.objects.get(project__id=project.id, meta_tag__name=meta_key_name, value=value)
		except MetaKeyProject.DoesNotExist:
			return None
		
	def get_project_metakey_last(self, project, meta_key_name, value):
		"""
		value = None, return a list
		"""
		try:
			if (value == None): query_set = MetaKeyProject.objects.filter(project__id=project.id, meta_tag__name=meta_key_name).order_by('-creation_date')
			else: query_set = MetaKeyProject.objects.filter(project__id=project.id, meta_tag__name=meta_key_name, value=value).order_by('-creation_date')
			if (query_set.count() > 0 ): return query_set[0]
			return None
		except MetaKeyProject.DoesNotExist:
			return None	

	def set_project_sample_metakey(self, project_sample, owner, meta_key_name, value, description):
		"""
		save a meta key
		"""
		try:
			metaKey = MetaKey.objects.get(name=meta_key_name)
		except MetaKey.DoesNotExist:
			metaKey = MetaKey()
			metaKey.name = meta_key_name
			metaKey.save()
		
		metaKeyProjectSample = MetaKeyProjectSample()
		metaKeyProjectSample.project_sample = project_sample
		metaKeyProjectSample.meta_tag = metaKey
		metaKeyProjectSample.owner = owner
		metaKeyProjectSample.value = value
		metaKeyProjectSample.description = description
		metaKeyProjectSample.save()
		return metaKeyProjectSample
	
	def get_project_sample_metakey(self, project_sample, meta_key_name, value):
		"""
		value = None, return a list
		"""
		try:
			if (value == None): return MetaKeyProjectSample.objects.filter(project_sample__id=project_sample.id, meta_tag__name=meta_key_name).order_by('-creation_date')
			return MetaKeyProjectSample.objects.get(project_sample__id=project_sample.id, meta_tag__name=meta_key_name, value=value)
		except MetaKeyProjectSample.DoesNotExist:
			return None

	def get_project_sample_metakey_last(self, project_sample, meta_key_name, value):
		"""
		value = None, return a list
		"""
		try:
			if (value == None): query_set = MetaKeyProjectSample.objects.filter(project_sample__id=project_sample.id, meta_tag__name=meta_key_name).order_by('-creation_date')
			else: query_set = MetaKeyProjectSample.objects.filter(project_sample__id=project_sample.id, meta_tag__name=meta_key_name, value=value).order_by('-creation_date')
			if (query_set.count() > 0 ): return query_set[0]
			return None
		except MetaKeyProject.DoesNotExist:
			return None	


