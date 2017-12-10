'''
Created on Dec 6, 2017

@author: mmp
'''

import os
from managing_files.models import Project, ProjectSample, DataSet, VacineStatus
from constants.constants import Constants, TypePath, FileExtensions
from utils.utils import Utils
from django.http import JsonResponse
from django.utils.safestring import mark_safe
from django.views.decorators.csrf import csrf_protect
from django.utils.translation import ugettext_lazy as _

######################################
###
###		AJAX methods for check box in session
###

@csrf_protect
def set_check_box_values(request):
	"""
	manage check boxes through ajax
	"""
	if request.is_ajax():
		data = { 'is_ok' : False }
		utils = Utils()
		if (Constants.CHECK_BOX_ALL in request.GET):
			request.session[Constants.CHECK_BOX_ALL] = utils.str2bool(request.GET.get(Constants.CHECK_BOX_ALL))
			## check all unique
			for key in request.session.keys():
				if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and utils.is_integer(key.split('_')[2])):
					request.session[key] = request.session[Constants.CHECK_BOX_ALL]
			data = {
				'is_ok': True
			}
		### one check box is pressed
		elif (Constants.CHECK_BOX in request.GET):
			request.session[Constants.CHECK_BOX_ALL] = False	### set All false
			key_name = "{}_{}".format(Constants.CHECK_BOX, request.GET.get(Constants.CHECK_BOX_VALUE))
			request.session[key_name] = utils.str2bool(request.GET.get(Constants.CHECK_BOX))
			total_checked = 0
			for key in request.session.keys():
				if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and utils.is_integer(key.split('_')[2])):
					if (request.session[key]): total_checked += 1
			data = {
				'is_ok': True,
				'total_checked' : total_checked,
			}
		## get the status of a check_box_single
		elif (Constants.GET_CHECK_BOX_SINGLE in request.GET):
			data['is_ok'] = True
			for key in request.session.keys():
				if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and utils.is_integer(key.split('_')[2])):
					data[key] = request.session[key]
		## change single status of a check_box_single
		elif (Constants.GET_CHANGE_CHECK_BOX_SINGLE in request.GET):
			data['is_ok'] = True
			key_name = "{}_{}".format(Constants.CHECK_BOX, request.GET.get(Constants.CHECK_BOX_VALUE))
			for key in request.session.keys():
				if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and utils.is_integer(key.split('_')[2])):
					if (request.session[key]): data[key] = False
					if (key == key_name): request.session[key] = utils.str2bool(request.GET.get(Constants.GET_CHANGE_CHECK_BOX_SINGLE))
					else: request.session[key] = False
		elif (Constants.COUNT_CHECK_BOX in request.GET):
			count = 0
			for key in request.session.keys():
				if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and utils.is_integer(key.split('_')[2])):
					if (request.session[key]): count += 1
				elif (key == Constants.CHECK_BOX_ALL and request.session[Constants.CHECK_BOX_ALL]): count += 1
			data = {
				'is_ok': True,
				Constants.COUNT_CHECK_BOX : count,
			}
		return JsonResponse(data)


###
###		END AJAX methods for check box in session
###
######################################

@csrf_protect
def show_phylo_canvas(request):
	"""
	manage check boxes through ajax
	"""
	if request.is_ajax():
		data = { 'is_ok' : False }
		utils = Utils()
		key_with_project_id = 'project_id'
		if (key_with_project_id in request.GET):
			project_id = int(request.GET.get(key_with_project_id))
			element_name = 'all_together'
			key_element_name = 'key_element_name'
			if (key_element_name in request.GET): element_name = int(request.GET.get(key_element_name))
			try:
				project = Project.objects.get(id=project_id)
				if (element_name == 'all_together'): file_name = project.get_global_file_by_project(TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_FASTTREE_tree)
				else: file_name = project.get_global_file_by_element(TypePath.MEDIA_ROOT, element_name, Project.PROJECT_FILE_NAME_FASTTREE_tree)
				if (os.path.exists(file_name)):
					string_file_content = utils.read_file_to_string(file_name).strip()
					if (string_file_content != None and len(string_file_content) > 0):
						data['is_ok'] = True
						data['tree'] = string_file_content
			except Project.DoesNotExist:
				pass
		return JsonResponse(data)

@csrf_protect
def get_image_coverage(request):
	"""
	get image coverage
	"""
	if request.is_ajax():
		data = { 'is_ok' : False }
		key_with_project_sample_id = 'project_sample_id'
		key_element = 'element'
		if (key_with_project_sample_id in request.GET and key_element in request.GET):
			try:
				project_sample = ProjectSample.objects.get(id=request.GET.get(key_with_project_sample_id))
				path_name = project_sample.get_global_file_by_element(TypePath.MEDIA_URL, ProjectSample.PREFIX_FILE_COVERAGE, request.GET.get(key_element), FileExtensions.FILE_PNG)
				data['is_ok'] = True
				data['image'] = mark_safe('<img src="{}" style="width: 100%;">'.format(path_name))
				data['text'] = _("Coverage for element '{}'".format(request.GET.get(key_element)))
			except ProjectSample.DoesNotExist as e:
				pass
		return JsonResponse(data)


@csrf_protect
def validate_project_reference_name(request):
	"""
	test if exist this project name
	"""
	if request.is_ajax():
		project_name = request.GET.get('project_name')
		data = {
			'is_taken': Project.objects.filter(name__iexact=project_name, owner__username=request.user.username).exists()
		}
		if (data['is_taken']): data['error_message'] = _('Exists a project with this name.')
		return JsonResponse(data)


@csrf_protect
def add_single_value_database(request):
	"""
	add a single value to a table in database
	possible tables to add: TagName, DataSet, VacineStatus
	"""
	if request.is_ajax():
		data = { 'is_ok' : False }
		key_type_data = 'type_data'
		key_value = 'value'
		value = request.GET[key_value].strip()
		if (key_value in request.GET and key_type_data in request.GET and len(value) > 0):
			## add to data_set
			if (request.GET[key_type_data] == 'id_data_set_add_modal'):
				try:
					DataSet.objects.get(name__iexact=value)
				except DataSet.DoesNotExist as e:
					dataset = DataSet()
					dataset.name = value
					dataset.owner = request.user
					dataset.save()
					
					### set the data for jscript
					data['is_ok'] = True
					data['text'] = dataset.name
					data['value'] = dataset.id
			## add to vaccine status
			elif (request.GET[key_type_data] == 'id_vaccine_add_modal'):
				try:
					VacineStatus.objects.get(name__iexact=value)
				except VacineStatus.DoesNotExist as e:
					vacine_status = VacineStatus()
					vacine_status.name = value
					vacine_status.owner = request.user
					vacine_status.save()
					
					### set the data for jscript
					data['is_ok'] = True
					data['text'] = vacine_status.name
					data['value'] = vacine_status.id
		return JsonResponse(data)


@csrf_protect
def remove_single_value_database(request):
	"""
	test if is prossible to remove
	"""
	if request.is_ajax():
		data = { 'is_ok' : False }
		key_type_data = 'type_data'
		key_value = 'value'
		key_is_to_test = 'is_to_test'
		
		if (key_is_to_test in request.GET and key_value in request.GET and key_type_data in request.GET and len(request.GET[key_value].strip()) > 0):
			utils = Utils()
			value = request.GET[key_value].strip()
			is_to_test = utils.str2bool(request.GET[key_is_to_test])
			
			## add to data_set
			if (request.GET[key_type_data] == 'id_data_set_remove_modal'):
				if (value == Constants.DATA_SET_GENERIC):
					data['is_ok'] = True;
					data['is_can_remove'] = False;
					data['message'] = "You can't remove 'Generic' data set";
				else:
					try:
						data_set = DataSet.objects.get(name__iexact=value)
						if (data_set.sample.count() > 0):
							if (is_to_test):
								data['is_ok'] = True;
								data['is_can_remove'] = False;
								data['message'] = "You can't remove '{}' name because as a relation in database.".format(value);
							else:
								data['is_ok'] = True;
								data['is_remove'] = False;
								data['message'] = "You can't remove '{}' name because as a relation in database.".format(value);
						else:
							if (is_to_test):
								data['is_ok'] = True;
								data['is_can_remove'] = True;
							else:
								data['is_ok'] = True;
								data['is_remove'] = True;
								data['value_to_remove'] = data_set.id;
								data['message'] = "'{}' was removed.".format(value);
								data_set.delete()
					except DataSet.DoesNotExist as e:
						pass

			## add to vaccine status
			elif (request.GET[key_type_data] == 'id_vaccine_remove_modal'):
				try:
					vacine_status = VacineStatus.objects.get(name__iexact=value)
					if (vacine_status.sample.count() > 0):
						if (is_to_test):
							data['is_ok'] = True;
							data['is_can_remove'] = False;
							data['message'] = "You can't remove '{}' name because has a relation in database.".format(value);
						else:
							data['is_ok'] = True;
							data['is_remove'] = False;
							data['message'] = "You can't remove '{}' name because has a relation in database.".format(value);
					else:
						if (is_to_test):
							data['is_ok'] = True;
							data['is_can_remove'] = True;
						else:
							data['is_ok'] = True;
							data['is_remove'] = False;
							data['message'] = "'{}' was removed.".format(value);
							data['value_to_remove'] = vacine_status.id;
							vacine_status.delete()
				except VacineStatus.DoesNotExist as e:
					data['is_ok'] = True;
					data['is_remove'] = False;
					data['message'] = "You can't remove '{}'.".format(value);
					pass
					
		return JsonResponse(data)


