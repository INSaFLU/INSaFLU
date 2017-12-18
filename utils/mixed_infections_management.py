'''
Created on Dec 14, 2017

@author: mmp
'''
from constants.meta_key_and_values import MetaKeyAndValue
from constants.constants_mixed_infection import ConstantsMixedInfection
from managing_files.manage_database import ManageDatabase
from managing_files.models import ProjectSample, MixedInfections, MixedInfectionsTag
from utils.result import DecodeObjects, MixedInfectionMainVector
from scipy import spatial

class MixedInfectionsManagement(object):
	'''
	classdocs
	'''

	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	def get_mixed_infections(self, project_sample, user, count_hits):
		"""
		return a MixedInfections instance and set an alert if necessary
		"""
		
		### calculate mixed infection value
		value = self.get_value_mixed_infection(count_hits)
		
		constants_mixed_infection = ConstantsMixedInfection()
		tag = constants_mixed_infection.get_tag_by_value(value)
		
		### get tag
		try:
			mixed_infections_tag = MixedInfectionsTag.objects.get(name=tag)
		except MixedInfectionsTag.DoesNotExist as e:
			mixed_infections_tag = MixedInfectionsTag()
			mixed_infections_tag.name = tag
			mixed_infections_tag.save()
		
		mixed_infections = MixedInfections()
		mixed_infections.tag = mixed_infections_tag
		mixed_infections.average_value = value
		mixed_infections.description = self.get_mixed_infection_main_vector().to_json()
		mixed_infections.save()
		
		##  set the alert
		if (constants_mixed_infection.is_alert(tag)):
			project_sample_ = ProjectSample.objects.get(pk=project_sample.id)
			project_sample_.alert_first_level += 1
			project_sample_.save()
			
			manage_database = ManageDatabase()
			manage_database.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION,\
									MetaKeyAndValue.META_VALUE_Success, "Warning, this sample has an average cosine distance " +\
									"of '{}'\nSuggest mixed infection".format(value))
		return mixed_infections
		
	
	def get_value_mixed_infection(self, count_hits):
		"""
		in: count_hits
		return float with cosine distance
		"""
		
		## return main vector
		mixed_infections_main_vector = self.get_mixed_infection_main_vector()

		### get the average of cosine distance
		f_total = 0.0
		vect_data_to_test = count_hits.get_vect_mixed_infections()
		for vect_data in mixed_infections_main_vector.get_vector():
	#		print("[{}-{}] - [{}-{}] {}".format(vect_data_to_test[0], vect_data_to_test[1], vect_data[0], vect_data[1],\
	#							1 - spatial.distance.cosine(vect_data, vect_data_to_test)))
			f_total += 1 - spatial.distance.cosine(vect_data, vect_data_to_test)
		
		if (len(mixed_infections_main_vector.get_vector()) == 0): return 0.0
		return f_total / float(len(mixed_infections_main_vector.get_vector()))
	
	
	def get_mixed_infection_main_vector(self):
		"""
		only has the positive main samples with counts
		get mixed infection main vector
		return: instance of MixedInfectionMainVector
		"""
		try:
			mixed_infections_main_vector = MixedInfections.objects.get(has_master_vector=True)
			decodeResult = DecodeObjects()
			return decodeResult.decode_result(mixed_infections_main_vector.description)
		except MixedInfections.DoesNotExist as e:
			mixed_infection_main_vector = MixedInfectionMainVector()
			mixed_infections_main_vector = MixedInfections()
			mixed_infections_main_vector.has_master_vector = True
			mixed_infections_main_vector.description = mixed_infection_main_vector.to_json()
			mixed_infections_main_vector.save()
			return mixed_infection_main_vector
		
		
		
		
		