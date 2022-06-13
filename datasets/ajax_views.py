'''
Created on Dec 6, 2017
@author: mmp
'''

import os
import csv, json, logging

from utils.utils import Utils
from django.http import JsonResponse
from django.utils.safestring import mark_safe
from django.views.decorators.csrf import csrf_protect
from extend_user.models import Profile
from django.conf import settings
from datetime import datetime
from django.db import transaction
from datasets.models import Dataset, Consensus, DatasetConsensus

### Logger
logger_debug = logging.getLogger("fluWebVirus.debug")
logger_production = logging.getLogger("fluWebVirus.production")
	
	
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
			
			data = { 
				'is_ok' : True,
				'dataset_name' : dataset_name_str,
				'id' : dataset.pk,
				'date_created' : dataset.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE),
				'message' : "The Dataset '{}' was created".format(dataset_name_str)
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

