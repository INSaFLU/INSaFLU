'''
Created on Nov 1, 2017

@author: mmp
'''
from managing_files.models import MetaKeySample, MetaKey, MetaKeyProject, MetaKeyProjectSample, MetaKeyReference
from managing_files.models import CountVariations, Statistics, TagName, ProjectSample
from constants.tag_names_constants import TagNamesConstants
from constants.meta_key_and_values import MetaKeyAndValue
from utils.lock_atomic_transaction import LockedAtomicTransaction
from django.db import connection

class ManageDatabase(object):
	'''
	classdocs
	'''
	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	def _get_metakey(self, meta_key_name):
		"""
		get metakey with lobk table
		"""
		with LockedAtomicTransaction(MetaKey):
			try:
				metaKey = MetaKey.objects.get(name=meta_key_name)
			except MetaKey.DoesNotExist:
				metaKey = MetaKey()
				metaKey.name = meta_key_name
				metaKey.save()
		return metaKey
	
	
	def set_reference_metakey(self, reference, owner, meta_key_name, value, description):
		"""
		save a meta key
		"""
		metaKey = self._get_metakey(meta_key_name)
		
		metaKeyReference = MetaKeyReference()
		metaKeyReference.reference = reference
		metaKeyReference.meta_tag = metaKey
		metaKeyReference.owner = owner
		metaKeyReference.value = value
		metaKeyReference.description = description
		metaKeyReference.save()
		return metaKeyReference

	def get_reference_metakey(self, reference, meta_key_name, value):
		"""
		value = None, return a list
		"""
		try:
			if (value == None): return MetaKeyReference.objects.filter(reference__id=reference.id, meta_tag__name=meta_key_name).order_by('-creation_date')
			return MetaKeyReference.objects.get(reference__id=reference.id, meta_tag__name=meta_key_name, value=value)
		except MetaKeyReference.DoesNotExist:
			return None
	
	def get_reference_metakey_last(self, reference, meta_key_name, value):
		"""
		value = None, return a list
		"""
		try:
			if (value == None): query_set = MetaKeyReference.objects.filter(reference__id=reference.id, meta_tag__name=meta_key_name).order_by('-creation_date')
			else: query_set = MetaKeyReference.objects.filter(reference__id=reference.id, meta_tag__name=meta_key_name, value=value).order_by('-creation_date')
			if (query_set.count() > 0 ): return query_set[0]
			return None
		except MetaKeyReference.DoesNotExist:
			return None	
		
	
	def set_sample_metakey(self, sample, owner, meta_key_name, value, description):
		"""
		save a meta key
		"""
		metaKey = self._get_metakey(meta_key_name)
		
		metaKeySample = MetaKeySample()
		metaKeySample.sample = sample
		metaKeySample.meta_tag = metaKey
		metaKeySample.owner = owner
		metaKeySample.value = value
		metaKeySample.description = description
		metaKeySample.save()
		return metaKeySample
	
	def get_sample_metakey(self, sample, meta_key_name, value):
		"""
		value = None, return a list
		"""
		try:
			if (value == None): return MetaKeySample.objects.filter(sample__id=sample.id, meta_tag__name=meta_key_name).order_by('-creation_date')
			return MetaKeySample.objects.get(sample__id=sample.id, meta_tag__name=meta_key_name, value=value)
		except MetaKeySample.DoesNotExist:
			return None
	
	def is_sample_wating_fastq_file(self, sample):
		"""
			META_VALUE_Error = "Error"
			META_VALUE_Success = "Success"
			META_VALUE_Queue = "Queue"
			the end of a process is 
		""" 
	
		meta_sample_queue = self.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue)
		if (meta_sample_queue is None): return True
		return False
	
	def is_sample_processing_step(self, sample):
		"""
			META_VALUE_Error = "Error"
			META_VALUE_Success = "Success"
			META_VALUE_Queue = "Queue"
			the end of a process is 
		""" 
	
		meta_sample_queue = self.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue)
		if (meta_sample_queue is None): return True
		
		### try to find errors
		number_errors = MetaKeySample.objects.filter(sample=sample, value=MetaKeyAndValue.META_VALUE_Error,
					creation_date__gt=meta_sample_queue.creation_date).order_by('-creation_date').count()
		### has a error
		if (number_errors > 0): return False

		### try to find queue_success 
		meta_sample_queue_success = self.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Success)		
		if (meta_sample_queue_success is None): return True
		if (meta_sample_queue_success.creation_date > meta_sample_queue.creation_date and \
			meta_sample_queue_success.description == meta_sample_queue.description): return False
		return True
	
	def get_sample_metakey_last(self, sample, meta_key_name, value):
		"""
		value = None, return a list
		"""
		try:
			if (value == None): query_set = MetaKeySample.objects.filter(sample__id=sample.id, meta_tag__name=meta_key_name).order_by('-creation_date')
			else: query_set = MetaKeySample.objects.filter(sample__id=sample.id, meta_tag__name=meta_key_name, value=value).order_by('-creation_date')
			if (query_set.count() > 0 ): return query_set[0]
			return None
		except MetaKeySample.DoesNotExist:
			return None	
		
	def remove_sample_start_metakey(self, sample, meta_key_name):
		"""
		in: meta_key_name + '%'
		return a list with match
		"""
		return MetaKeySample.objects.filter(sample__id=sample.id,\
					meta_tag__name__startswith=meta_key_name).delete()
					
	def set_project_metakey(self, project, owner, meta_key_name, value, description):
		"""
		save a meta key
		"""
		metaKey = self._get_metakey(meta_key_name)
		
		metaKeyProject = MetaKeyProject()
		metaKeyProject.project = project
		metaKeyProject.meta_tag = metaKey
		metaKeyProject.owner = owner
		metaKeyProject.value = value
		metaKeyProject.description = description
		metaKeyProject.save()
		return MetaKeyProject
	
	def update_project_metakey(self, project, owner, meta_key_name, value, description):
		"""
		update the meta_key, if not exist create other
		"""
		metaKey = self._get_metakey(meta_key_name)
		
		metaKeyProject_list = MetaKeyProject.objects.filter(project=project, meta_tag=metaKey, owner=owner)
		if (metaKeyProject_list.count() == 0):
			return self.set_project_metakey(project, owner, meta_key_name, value, description)
		else:
			meta_key_project = metaKeyProject_list[0]
			meta_key_project.value = value
			meta_key_project.description = description
			meta_key_project.save()
		return meta_key_project

	
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


	def is_project_processing_step(self, project, tag_key):
		"""
			META_VALUE_Error = "Error"
			META_VALUE_Success = "Success"
			META_VALUE_Queue = "Queue"
			the end of a process is 
		""" 
	
		meta_project_queue = self.get_project_metakey_last(project, tag_key, MetaKeyAndValue.META_VALUE_Queue)
		if (meta_project_queue is None): return True
		
		### try to find errors
		number_errors = MetaKeySample.objects.filter(project=project, value=MetaKeyAndValue.META_VALUE_Error,
					creation_date__gt=meta_project_queue.creation_date).order_by('-creation_date').count()
		### has a error
		if (number_errors > 0): return False

		### try to find queue_success 
		meta_project_queue_success = self.get_project_metakey_last(project, tag_key, MetaKeyAndValue.META_VALUE_Success)		
		if (meta_project_queue_success is None): return True
		if (meta_project_queue_success.creation_date > meta_project_queue.creation_date and \
			meta_project_queue_success.description == meta_project_queue.description): return False
		return True
	
	def set_project_sample_metakey(self, project_sample, owner, meta_key_name, value, description):
		"""
		save a meta key
		"""
		metaKey = self._get_metakey(meta_key_name)
		
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

	def get_project_sample_starts_with_metakey(self, project_sample, meta_key_name):
		"""
		in: meta_key_name + '%'
		return a list with match
		"""
		try:
			return list(MetaKeyProjectSample.objects.filter(project_sample__id=project_sample.id, meta_tag__name__startswith=meta_key_name).order_by('-creation_date'))
		except MetaKeyProject.DoesNotExist:
			return []

	def remove_project_sample_start_metakey(self, project_sample, meta_key_name):
		"""
		in: meta_key_name + '%'
		return a list with match
		"""
		return MetaKeyProjectSample.objects.filter(project_sample__id=project_sample.id,\
					meta_tag__name__startswith=meta_key_name).delete()
	


	def get_variation_count(self, count_hits):
		"""
		Add values in project_sample
		"""
		count_variations = CountVariations()
		if (not count_hits is None):
			count_variations.var_bigger_50_90 = count_hits.get_hits_50_90()
			count_variations.var_bigger_90 = count_hits.get_hits_more_90()
			count_variations.var_less_50 = count_hits.get_hits_less_50()
			count_variations.total = count_hits.get_total()
		count_variations.save()
		return count_variations

	def is_sample_downsized(self, sample):
		"""
		:out True  
		""" 
	
		meta_sample_queue = self.get_sample_metakey_last(sample,
			MetaKeyAndValue.META_KEY_ALERT_DOWNSIZE_OF_FASTQ_FILES,
			MetaKeyAndValue.META_VALUE_Success)
		if (not meta_sample_queue is None): return True
		return False
											
	#######################################
	###
	###		Other methods
	###
	
	def get_max_length_label(self, project, user, b_calculate_again):
		"""
		Get the max length of the samples in a specific project in chars
		b_calculate_again = True calculate again
		"""
		
		if (not b_calculate_again):
			meta_data = self.get_project_metakey_last(project, MetaKeyAndValue.META_KEY_Project_max_sample_length,\
										MetaKeyAndValue.META_VALUE_Success)
			if (meta_data != None):
				return int(meta_data.description)
			
		### calculate and save
		n_max_value = 0
		for project_sample in project.project_samples.all():
			if (not project_sample.is_finished): continue
			if (project_sample.is_deleted): continue
			if (project_sample.is_error): continue
			if (len(project_sample.sample.name) > n_max_value): n_max_value = len(project_sample.sample.name)
		## calculate for the reference too
		if (len(project.reference.display_name) > n_max_value): n_max_value = len(project.reference.display_name)
			
		meta_data = self.update_project_metakey(project, user, MetaKeyAndValue.META_KEY_Project_max_sample_length,\
							MetaKeyAndValue.META_VALUE_Success, str(n_max_value))
		return n_max_value
	###
	###
	###
	#######################################
	
	#######################################
	###
	###		deal with percentils, not used now
	###

	
	def percentils_refresh(self, user):
		"""
		refresh the percentils
		"""
		## don't calc the percentils
		if (CountVariations.objects.count() < TagNamesConstants.MINIMUM_NUMBER_TO_COUNT_STATISTICS): return
		
		tagNamesConstants = TagNamesConstants()
		for percentil_tag in tagNamesConstants.get_all_tags_percentil():
			try: 
				statistics = Statistics.objects.get(tag__name=percentil_tag)
			except Statistics.DoesNotExist:
				
				with LockedAtomicTransaction(TagName):
					try:
						tag_name = TagName.objects.get(name=percentil_tag, owner__id=user.id)
					except TagName.DoesNotExist:
						tag_name = TagName()
						tag_name.name = percentil_tag
						tag_name.owner = user
						tag_name.is_meta_data = tagNamesConstants.is_meta_tag_name(percentil_tag)
						tag_name.save()
				
				statistics = Statistics()
				statistics.tag = tag_name
			
			### set the values
			percentil_value = self.get_percentil_value(percentil_tag)
			if (percentil_value != None):
				statistics.value = percentil_value
				statistics.save()
	
	
	def set_percentis_alert(self, project_sample, user, count_hits, percentil_name):
		"""
		in: percentil_name from tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_95, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL), or other
		return True if has an alert
		"""
		tagNamesConstants = TagNamesConstants()
		
		### get the percentil for total
		percentil_value = self.get_percentil_value_from_db(percentil_name)
		if (percentil_value != None):
			b_percentil_fault = False
			if (percentil_name.endswith(TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL) and count_hits.get_total() > percentil_value): b_percentil_fault = True
			elif (percentil_name.endswith(TagNamesConstants.TAG_PERCENTIL_VAR_50) and count_hits.get_hits_less_50() > percentil_value): b_percentil_fault = True
			elif (percentil_name.endswith(TagNamesConstants.TAG_PERCENTIL_VAR_50_90) and count_hits.get_hits_50_90() > percentil_value): b_percentil_fault = True
			
			if (b_percentil_fault):
				### set the alert
				project_sample = ProjectSample.objects.get(pk=project_sample.id)
				project_sample.alert_first_level += 1
				project_sample.save()
				self.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_ALERT_COUNT_VAR,\
								MetaKeyAndValue.META_VALUE_Error,\
								"Warning, this sample has a total of variations '{}' bigger than the percentile {} of all database".\
								format(count_hits.get_total(), int(tagNamesConstants.get_number_percentil_from_tag(percentil_name))))
			return True
		return False


	def get_percentil_value(self, percentil_name):
		"""
		in: percentil_name from tagNamesConstants.get_percentil_tag_name(TagNamesConstants.TAG_PERCENTIL_95, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL)
		return  the float percentil value from database
		
		select total_matched, ntile(100) over (order by total_matched), cume_dist() OVER (ORDER BY total_matched) as percentile from runner;

		select avg(total_matched) from (select total_matched, ntile(100) over (order by total_matched), cume_dist() OVER (ORDER BY total_matched) as percentile from market_book) where ntile = 99;
		select avg(total_matched) from (select total_matched, ntile(100) over (order by total_matched), cume_dist() OVER (ORDER BY total_matched) as percentile from market_book) ss where ntile = 99;
		select avg(total) from (select total, ntile(100) over (order by total), cume_dist() OVER (ORDER BY  total) as percentile from managing_files_countvariations) ss where ntile = 99;
		"""
		## get the percentil float to get the value, 97.5, 95, 90, or other
		tagNamesConstants = TagNamesConstants()
		percentil_float = tagNamesConstants.get_number_percentil_from_tag(percentil_name)
		
		if (CountVariations.objects.count() < TagNamesConstants.MINIMUM_NUMBER_TO_COUNT_STATISTICS):
			return tagNamesConstants.get_emperical_percentil_values(percentil_float)
		
		## calculate for total
		cursor = connection.cursor()
		if (tagNamesConstants.is_which_var(percentil_name, TagNamesConstants.TAG_PERCENTIL_VAR_TOTAL)):
			cursor.execute('select avg(total) from (select total, ntile(100) over (order by total), cume_dist() OVER (ORDER BY total)' +\
					' as percentile from managing_files_countvariations) ss where ntile = {};'.format(percentil_float))
		elif (tagNamesConstants.is_which_var(percentil_name, TagNamesConstants.TAG_PERCENTIL_VAR_50)):
			cursor.execute('select avg(var_less_50) from (select var_less_50, ntile(100) over (order by var_less_50), cume_dist() OVER (ORDER BY var_less_50)' +\
					' as percentile from managing_files_countvariations) ss where ntile = {};'.format(percentil_float))
		elif (tagNamesConstants.is_which_var(percentil_name, TagNamesConstants.TAG_PERCENTIL_VAR_50_90)):
			cursor.execute('select avg(var_bigger_50) from (select var_bigger_50, ntile(100) over (order by var_bigger_50), cume_dist() OVER (ORDER BY var_bigger_50)' +\
					' as percentile from managing_files_countvariations) ss where ntile = {};'.format(percentil_float))
		row = cursor.fetchall()
		if (row != None and len(row) > 0 and len(row[0]) > 0): return int(row[0][0])
		return None


	def get_percentil_value_from_db(self, percentil_tag):
		"""
		get the statistic that is in database
		return None if not exist 
		"""
		try: 
			statistics = Statistics.objects.get(tag__name=percentil_tag)
			return statistics.value
		except Statistics.DoesNotExist:
			pass
		return None
		
	###
	###		END deal with percentils
	###
	#######################################
	
	
		