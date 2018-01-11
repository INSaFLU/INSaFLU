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

		## Old version of cosine, not used anymore
		#constants_mixed_infection = ConstantsMixedInfection()
		#tag = constants_mixed_infection.get_tag_by_value(value)
		
		tag = ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO
		## test ratio method
		if (count_hits.is_mixed_infection_ratio_test() or count_hits.total_grather_than_mixed_infection()):	## doesn't matter the other
			tag = ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES
			
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
		manage_database = ManageDatabase()
#		if (constants_mixed_infection.is_alert(tag) or count_hits.is_mixed_infection_ratio_test()):
		if (count_hits.is_mixed_infection_ratio_test() or count_hits.total_grather_than_mixed_infection()):
			project_sample_ = ProjectSample.objects.get(pk=project_sample.id)
			project_sample_.alert_first_level += 1
			project_sample_.save()
			
# 			manage_database.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_COSINE_DISTANCE,\
# 									MetaKeyAndValue.META_VALUE_Success, "Warning, this sample has an average cosine distance " +\
# 									"of '{}'.\nSuggest mixed infection.".format(value))
			
			## test mixed infection by empirical values
			if (count_hits.is_mixed_infection_ratio_test()):
				manage_database.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_RATIO_TEST,\
								MetaKeyAndValue.META_VALUE_Success, "Warning, this sample has a ratio of '{}' and a total of '{}' variations.\nSuggest mixed infection.".\
								format(count_hits.get_mixed_infection_ratio_str(), count_hits.get_total_50_50_90()))
			elif (count_hits.total_grather_than_mixed_infection()):
				manage_database.set_project_sample_metakey(project_sample, user, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_SUM_TEST,\
								MetaKeyAndValue.META_VALUE_Success, "Warning, this sample has a sum '{}' bigger than '{}'.\nSuggest mixed infection.".\
								format(count_hits.get_total_50_50_90(), count_hits.TOTAL_GRATHER_THAN_MIXED))
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
		if (sum(vect_data_to_test) == 0): return 0.0
		
		for vect_data in mixed_infections_main_vector.get_vector():
#			print("[{}-{}] - [{}-{}] {}".format(vect_data[0], vect_data[1], vect_data_to_test[0], vect_data_to_test[1],\
#								1 - spatial.distance.cosine(vect_data, vect_data_to_test)))
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
		
		
		
		
		