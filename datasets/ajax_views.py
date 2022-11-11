'''
Created on Dec 6, 2017
@author: mmp
'''

import logging, os, csv, json
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_protect
from extend_user.models import Profile
from django.conf import settings
from datetime import datetime
from django.db import transaction
from django.utils.safestring import mark_safe
from datasets.models import Dataset, Consensus, DatasetConsensus
from settings.default_parameters import DefaultParameters
from settings.models import Software
from utils.process_SGE import ProcessSGE
from utils.utils import Utils
from constants.constants import TypePath
from datasets.manage_database import ManageDatabase

### Logger
if settings.DEBUG: logger = logging.getLogger("fluWebVirus.debug")
else: logger = logging.getLogger("fluWebVirus.production")
	
######################################
###
###		AJAX methods for check box in session
###

@transaction.atomic
@csrf_protect
def remove_dataset(request):
	"""
	remove a dataset.
	"""
	if request.is_ajax():
		data = { 'is_ok' : False }
		dataset_id_a = 'dataset_id'
		
		if (dataset_id_a in request.GET):
			
			## some pre-requisites
			if (not request.user.is_active or not request.user.is_authenticated): return JsonResponse(data)
			try:
				profile = Profile.objects.get(user__pk=request.user.pk)
			except Profile.DoesNotExist:
				return JsonResponse(data)
			if (profile.only_view_project): return JsonResponse(data)
			
			dataset_id = request.GET[dataset_id_a]
			try:
				dataset = Dataset.objects.get(pk=dataset_id)
			except Dataset.DoesNotExist:
				return JsonResponse(data)
			
			## different owner or belong to a project not deleted
			if (dataset.owner.pk != request.user.pk): return JsonResponse(data)
			
			### now you can remove
			dataset.is_deleted = True
			dataset.is_deleted_in_file_system = False
			dataset.date_deleted = datetime.now()
			dataset.save()
			
			### delete all Sequences Datasets samples
			### this is only necessary for consistency
			for dataset_consensus in dataset.dataset_consensus.all():
				dataset_consensus.is_deleted = True
				dataset_consensus.save()
			data = { 'is_ok' : True }
		return JsonResponse(data)


@transaction.atomic
@csrf_protect
def add_dataset_name(request):
	"""
	add new dataset
	"""
	if request.is_ajax():
		data = { 'is_ok' : False }
		dataset_name = 'dataset_name'
		
		if (dataset_name in request.GET):
			
			## some pre-requisites
			if (not request.user.is_active or not request.user.is_authenticated): return JsonResponse(data)
			try:
				profile = Profile.objects.get(user__pk=request.user.pk)
			except Profile.DoesNotExist:
				return JsonResponse(data)
			if (profile.only_view_project): return JsonResponse(data)
			
			dataset_name_str = request.GET[dataset_name].strip()
			if len(dataset_name_str) == 0: return JsonResponse(data)
			
			try:
				dataset = Dataset.objects.get(name=dataset_name_str, owner__pk=request.user.pk,
										is_deleted=False)
				
				### name exist
				data['message'] = "This name '{}' already exist.".format(dataset_name_str)
				return JsonResponse(data)
			except Dataset.DoesNotExist:
				pass
	
			dataset = Dataset()
			### now you can remove
			dataset.name = dataset_name_str
			dataset.owner = profile.user
			dataset.save()

			message = "The Dataset '{}' was created".format(dataset_name_str)
			# TODO If there are specified parameters, change them before saving...
			try:
				default_parameters =  DefaultParameters()
				vect_parameters = default_parameters.get_nextstrain_default(user=request.user, dataset=dataset)
				# This will repeat the same software over and over for each dataset, but well...
				# TODO reuse existing software entry (only save new parameter)
				default_parameters.persist_parameters(vect_parameters, type_of_use=Software.TYPE_OF_USE_dataset)
			except Exception as e:
				message = "The Dataset '{}' was created, but with problems: {}".format(dataset_name_str, str(e))				

			data = { 
				'is_ok' : True,
				'dataset_name' : dataset_name_str,
				'id' : dataset.pk,
				'date_created' : dataset.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE),
				'message' : message
			}
		return JsonResponse(data)



@csrf_protect
def test_dataset_name(request):
	"""
	test dataset name, Return True if exits
	"""
	if request.is_ajax():
		data = { 'is_taken' : False }
		dataset_name = 'dataset_name'
		
		if (dataset_name in request.GET):
			
			## some pre-requisites
			if (not request.user.is_active or not request.user.is_authenticated): return JsonResponse(data)
			try:
				profile = Profile.objects.get(user__pk=request.user.pk)
			except Profile.DoesNotExist:
				return JsonResponse(data)
			if (profile.only_view_project): return JsonResponse(data)
			
			dataset_name_str = request.GET[dataset_name].strip()
			if (len(dataset_name_str) == 0): JsonResponse(data)
				
			try:
				dataset = Dataset.objects.get(name=dataset_name_str, owner__pk=request.user.pk,
											is_deleted=False)
				### name exist
				data = { 
					'is_taken' : True,
					'error_message' : "This name '{}' already exists".format(dataset_name_str)
				}
			except Dataset.DoesNotExist:
				return JsonResponse(data)
		return JsonResponse(data)


@csrf_protect
def test_consensus_name(request):
	"""
	test dataset name, Return True if exits
	"""
	if request.is_ajax():
		data = { 'is_taken' : False }
		consensus_name = 'consensus_name'
		
		if (consensus_name in request.GET):
			
			## some pre-requisites
			if (not request.user.is_active or not request.user.is_authenticated): return JsonResponse(data)
			try:
				profile = Profile.objects.get(user__pk=request.user.pk)
			except Profile.DoesNotExist:
				return JsonResponse(data)
			if (profile.only_view_project): return JsonResponse(data)
			
			dataset_name_str = request.GET[consensus_name].strip()
			if (len(dataset_name_str) == 0): JsonResponse(data)
				
			try:
				dataset = Consensus.objects.get(name=dataset_name_str, owner__pk=request.user.pk,
											is_deleted=False)
				### name exist
				data = { 
					'is_taken' : True,
					'error_message' : "This name '{}' already exists".format(dataset_name_str)
				}
			except Consensus.DoesNotExist:
				return JsonResponse(data)
		return JsonResponse(data)

@csrf_protect
def add_consensus_name(request):
	"""
	add new dataset
	"""
	if request.is_ajax():
		data = { 'is_ok' : False }
		consensus_name = 'consensus_name'
		
		if (consensus_name in request.POST):
			
			## some pre-requisites
			if (not request.user.is_active or not request.user.is_authenticated): return JsonResponse(data)
			try:
				profile = Profile.objects.get(user__pk=request.user.pk)
			except Profile.DoesNotExist:
				return JsonResponse(data)
			if (profile.only_view_project): return JsonResponse(data)
			
			consensus_name_str = request.GET[consensus_name]
			try:
				dataset = Dataset.objects.get(name=consensus_name_str, owner__pk=request.user.pk,
										is_deleted=False)
				
				### name exist
				data['message'] = "This name '{}' already exist.".format(consensus_name_str)
				return JsonResponse(data)
			except Dataset.DoesNotExist:
				pass
	
			consensus = Consensus()
			### now you can remove
			consensus.name = consensus_name_str
			consensus.owner = profile.user
			consensus.save()
			
			data = { 
				'is_ok' : True,
				'consensus_name' : consensus_name_str,
				'id' : dataset.pk,
				'date_created' : dataset.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE),
				'message' : "The Consensus '{}' was created".format(consensus_name_str)
			}
		return JsonResponse(data)

@transaction.atomic
@csrf_protect
def remove_consensus(request):
	"""
	remove a dataset.
	"""
	if request.is_ajax():
		data = { 'is_ok' : False }
		consensus_id_a = 'consensus_id'
		
		if (consensus_id_a in request.GET):
			
			## some pre-requisites
			if (not request.user.is_active or not request.user.is_authenticated): return JsonResponse(data)
			try:
				profile = Profile.objects.get(user__pk=request.user.pk)
			except Profile.DoesNotExist:
				return JsonResponse(data)
			if (profile.only_view_project): return JsonResponse(data)
			
			consensus_id = request.GET[consensus_id_a]
			try:
				consensus = Consensus.objects.get(pk=consensus_id)
			except Dataset.DoesNotExist:
				return JsonResponse(data)
			
			## different owner or belong to a project not deleted
			if (consensus.owner.pk != request.user.pk): return JsonResponse(data)
			
			### now you can remove
			consensus.is_deleted = True
			consensus.is_deleted_in_file_system = False
			consensus.date_deleted = datetime.now()
			consensus.save()
			
			data = { 'is_ok' : True }
		return JsonResponse(data)


@transaction.atomic
@csrf_protect
def remove_consensus_in_dataset(request):
	"""
	remove a dataset.
	"""
	if request.is_ajax():
		data = { 
			'is_ok' : False,
			'message' : "Something went wrong. Fail to remove."
		}
		consensus_id_a = 'consensus_id'
		
		if (consensus_id_a in request.GET):
			
			## some pre-requisites
			if (not request.user.is_active or not request.user.is_authenticated): return JsonResponse(data)
			try:
				profile = Profile.objects.get(user__pk=request.user.pk)
			except Profile.DoesNotExist:
				
				return JsonResponse(data)
			if (profile.only_view_project): return JsonResponse(data)
			
			consensus_id = request.GET[consensus_id_a]
			try:
				dataset_consensus = DatasetConsensus.objects.get(pk=consensus_id)
			except Dataset.DoesNotExist:
				data['message'] = "Something went wrong. Fail to remove."
				return JsonResponse(data)
			
			## different owner or belong to a project not deleted
			if (dataset_consensus.dataset.owner.pk != request.user.pk): return JsonResponse(data)
			
			### test how many references exist in this dataset
			if not dataset_consensus.reference is None and \
				DatasetConsensus.objects.filter(is_deleted=False, is_error=False,
							reference__isnull=False).count() < 2:
				data['message'] = "You can not remove the remain reference sequence. At least one must be present."
				return JsonResponse(data)
			
			### now you can remove
			dataset_consensus.is_deleted = True
			dataset_consensus.date_deleted = datetime.now()
			dataset_consensus.save()
			
			if not dataset_consensus.project_sample is None:
				dataset_consensus.dataset.number_of_sequences_from_projects -= 1
				data['message'] = "The consensus sequence {} from project {} was removed".format(
					dataset_consensus.project_sample.sample.name,
					dataset_consensus.project_sample.project.name)
				dataset_consensus.dataset.save() 
			elif not dataset_consensus.consensus is None:
				dataset_consensus.dataset.number_of_sequences_from_consensus -= 1
				data['message'] = "The consensus {} was removed".format(
					dataset_consensus.consensus.name)
				dataset_consensus.dataset.save()
			elif not dataset_consensus.reference is None:
				dataset_consensus.dataset.number_of_sequences_from_references -= 1
				data['message'] = "The reference {} was removed".format(
					dataset_consensus.reference.name)
				dataset_consensus.dataset.save()
			
			## need to run processing
			try:
				process_SGE = ProcessSGE()
				taskID =  process_SGE.set_collect_dataset_global_files(dataset_consensus.dataset,
										request.user)
			except:
				pass
			data['is_ok'] = True
		return JsonResponse(data)


@csrf_protect
def validate_consensus_name(request):
	"""
	test if exist this reference name
	"""
	if request.is_ajax():
		consensus_name = request.GET.get('consensus_name')
		
		data = {
			'is_taken': Consensus.objects.filter(name__iexact=consensus_name,
				is_deleted=False, owner=request.user).exists()
		}
		if (data['is_taken']): data['error_message'] = 'Exists a consensus with this name.'
		return JsonResponse(data)


@csrf_protect
def dataset_rebuild(request):
	"""
	Rebuild results
	"""
	if request.is_ajax():

		data = { 
			'is_ok' : False,
			'message' : "Something went wrong."
		}
		key_with_dataset_id = 'dataset_id'
		if (key_with_dataset_id in request.GET):
			dataset_id = int(request.GET.get(key_with_dataset_id))
			try:
				dataset = Dataset.objects.get(id=dataset_id)
				## need to run processing
				try:
					process_SGE = ProcessSGE()
					process_SGE.set_collect_dataset_global_files(dataset, request.user)
					data['is_ok'] = True
					data['message'] = "alls well that ends well."					
				except:
					data['message'] = "Something went wrong. Could not Rebuild Dataset."
			except Dataset.DoesNotExist:
				data['message'] = "Something went wrong. Dataset could not be found."

		return JsonResponse(data)


@csrf_protect
def show_msa_nucleotide(request):
	"""
	manage msa nucleotide alignments
	"""
	if request.is_ajax():
		data = { 'is_ok' : False }
		key_with_dataset_id = 'dataset_id'
		if (key_with_dataset_id in request.GET):
			dataset_id = int(request.GET.get(key_with_dataset_id))
			try:
				manage_database = ManageDatabase()
				dataset = Dataset.objects.get(id=dataset_id)
				file_name_fasta = dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_MAFFT)
				
				if (os.path.exists(file_name_fasta)):
					file_name_fasta = dataset.get_global_file_by_dataset(TypePath.MEDIA_URL, Dataset.DATASET_FILE_NAME_MAFFT)
					data['alignment_fasta_show_id'] = mark_safe(request.build_absolute_uri(file_name_fasta))
					url_file_name_fasta = '<a href="{}" download="{}"> {}</a>'.format(file_name_fasta,
						os.path.basename(file_name_fasta), os.path.basename(file_name_fasta))
				else: 
					url_file_name_fasta = 'File not available'
					data['alignment_fasta_show_id'] = '#'
				
				file_name_nex = dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_nex)
				if (os.path.exists(file_name_nex)):
					file_name_nex = dataset.get_global_file_by_dataset(TypePath.MEDIA_URL, Dataset.DATASET_FILE_NAME_nex)
					url_file_name_nex = '<a href="{}" download="{}"> {}</a>'.format(file_name_nex,
						os.path.basename(file_name_fasta), os.path.basename(file_name_nex))
				else: 
					url_file_name_nex = 'File not available'
						
				data['is_ok'] = True
				data['alignment_fasta_id'] = mark_safe("<strong>Alignment (.fasta):</strong> {}".format(url_file_name_fasta))
				data['alignment_nex_id'] = mark_safe("<strong>Alignment (.nex):</strong> {}".format(url_file_name_nex))
				b_calculate_again = False
				data['max_length_label'] = manage_database.get_max_length_label(dataset, request.user, b_calculate_again)
			except Dataset.DoesNotExist:
				pass
		return JsonResponse(data)

@csrf_protect
def show_phylo_canvas(request):
	"""
	manage check boxes through ajax
	"""
	
	if request.is_ajax():
		data = { 'is_ok' : False }
		utils = Utils()
		key_with_dataset_id = 'dataset_id'
		if (key_with_dataset_id in request.GET):
			dataset_id = int(request.GET.get(key_with_dataset_id))
			try:
				vect_remove_keys = ['strain', 'fastq1', 'fastq2', 'data set', 'latitude', 'longitude']
				dataset = Dataset.objects.get(id=dataset_id)
				file_name_root_json = dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_RESULT_json)
				file_name_url_json = dataset.get_global_file_by_dataset(TypePath.MEDIA_URL, Dataset.DATASET_FILE_NAME_RESULT_json)
				## the input it is DATASET_FILE_NAME_RESULT_NEXTSTRAIN_CSV
				file_name_root_sample = dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_RESULT_NEXTSTRAIN_CSV)
					
				file_name_root_nwk = dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_FASTTREE_tree)
				file_name_nwk = dataset.get_global_file_by_dataset(TypePath.MEDIA_URL, Dataset.DATASET_FILE_NAME_FASTTREE)
				file_name_tree = dataset.get_global_file_by_dataset(TypePath.MEDIA_URL, Dataset.DATASET_FILE_NAME_FASTTREE_tree)

				if (os.path.exists(file_name_root_nwk) and os.path.exists(file_name_root_sample)):
					string_file_content = utils.read_file_to_string(file_name_root_nwk).strip()
					
					if not os.path.exists(file_name_root_json) or os.path.getsize(file_name_root_json) == 0:
						with open(file_name_root_json, 'w', encoding='utf-8') as handle_write, open(file_name_root_sample) as handle_in_csv:
							reader = csv.DictReader(handle_in_csv)
							all_data = json.loads(json.dumps(list(reader)))
							dt_result = {}
							for dict_data in all_data:
								if ('strain' in dict_data):
									dt_out = dict_data.copy()
									for key_to_remove in vect_remove_keys:
										try:
											del dt_out[key_to_remove]
										except KeyError:
											pass
									dt_result[dict_data['strain']] = {key: '' if dt_out[key] == '?' else dt_out[key] for key in dt_out}
							if len(dt_result) == len(all_data):
								handle_write.write(json.dumps(dt_result))
							else:
								logger.error('DatasetID: {}  different number of lines processing Sample {} -> JSON {}'.format(
									dataset_id, len(dt_result), len(all_data)))
								string_file_content = None	## return error
								
					if (string_file_content != None and len(string_file_content) > 0):
						data['is_ok'] = True
						data['tree'] = string_file_content
						data['root'] = dataset.get_first_reference_name()
						data['url_sample'] = file_name_url_json			
						data['tree_nwk_id'] = mark_safe('<strong>Tree (.nwk):</strong> <a href="{}" download="{}"> {}</a>'.format(file_name_nwk, 
														os.path.basename(file_name_nwk), os.path.basename(file_name_nwk)))
						data['tree_tree_id'] = mark_safe('<strong>Tree (.tree):</strong> <a href="{}" download="{}"> {}</a>'.format(file_name_tree,
														os.path.basename(file_name_tree), os.path.basename(file_name_tree)))
			except Dataset.DoesNotExist:
				pass
		return JsonResponse(data)
