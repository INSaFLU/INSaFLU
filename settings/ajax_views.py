'''
Created on Dec 6, 2017

@author: mmp
'''
from settings.models import Software
from settings.default_software import DefaultSoftware
from settings.default_software_project_sample import DefaultProjectSoftware
from managing_files.models import Project, ProjectSample, Sample
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_protect, csrf_exempt
from extend_user.models import Profile
from constants.meta_key_and_values import MetaKeyAndValue
from managing_files.manage_database import ManageDatabase
from settings.constants_settings import ConstantsSettings
from utils.process_SGE import ProcessSGE
from utils.result import MaskingConsensus, DecodeObjects


@csrf_protect
def set_default_parameters(request):
	"""
	remove a reference. It can only be removed if not belongs to any deleted project
	"""
	if request.is_ajax():
		data = { 'is_ok' : False }
		software_id_a = 'software_id'
		project_id_a = 'project_id'
		project_sample_id_a = 'project_sample_id'
		sample_id_a = 'sample_id'
		
		## some pre-requisites
		if (not request.user.is_active or not request.user.is_authenticated): return JsonResponse(data)
		try:
			profile = Profile.objects.get(user__pk=request.user.pk)
		except Profile.DoesNotExist:
			return JsonResponse(data)
		if (profile.only_view_project): return JsonResponse(data)
			
		if (software_id_a in request.GET):
			software_id = request.GET[software_id_a]
			b_change = False
			try:
				software = Software.objects.get(pk=software_id)
				if (project_id_a in request.GET):
					project_id = request.GET[project_id_a]
					project = Project.objects.get(pk=project_id)
					
					default_project_software = DefaultProjectSoftware()
					default_project_software.set_default_software(software, request.user, Software.TYPE_OF_USE_project, project,
											None, None)
					## set a new default
					data['default'] = default_project_software.get_parameters(software.name, request.user, Software.TYPE_OF_USE_project,
													project, None, None, software.technology.name)
					
					b_change = default_project_software.is_change_values_for_software(software.name, software.technology.name)
				elif (project_sample_id_a in request.GET):
					project_sample_id = request.GET[project_sample_id_a]
					project_sample = ProjectSample.objects.get(pk=project_sample_id)
					
					default_project_software = DefaultProjectSoftware()
					default_project_software.set_default_software(software, request.user, Software.TYPE_OF_USE_project_sample,
										None, project_sample, None)
					## set a new default
					data['default'] = default_project_software.get_parameters(software.name, request.user, Software.TYPE_OF_USE_project_sample,\
													None, project_sample, None, software.technology.name)
					
					### need to re-run this sample with snippy if the values change
					if (default_project_software.is_change_values_for_software(software.name, ConstantsSettings.TECHNOLOGY_illumina \
								if project_sample.is_sample_illumina() else ConstantsSettings.TECHNOLOGY_minion)):
						b_change = True
						
						### re-run data
						metaKeyAndValue = MetaKeyAndValue()
						manageDatabase = ManageDatabase()
						process_SGE = ProcessSGE()
						
						### change flag to nor finished
						project_sample.is_finished = False
						project_sample.save()
		
						### create a task to perform the analysis of snippy and freebayes
						try:
							(job_name_wait, job_name) = request.user.profile.get_name_sge_seq(Profile.SGE_GLOBAL)
							if (project_sample.is_sample_illumina()):
								taskID = process_SGE.set_second_stage_snippy(project_sample, request.user, job_name, job_name_wait)
							else:
								taskID = process_SGE.set_second_stage_medaka(project_sample, request.user, job_name, job_name_wait)
								
							### set project sample queue ID
							manageDatabase.set_project_sample_metakey(project_sample, request.user,\
											metaKeyAndValue.get_meta_key_queue_by_project_sample_id(project_sample.id),\
											MetaKeyAndValue.META_VALUE_Queue, taskID)
							
							### need to collect global files again
							taskID = process_SGE.set_collect_global_files(project, request.user)
							manageDatabase.set_project_metakey(project, request.user, metaKeyAndValue.get_meta_key(\
									MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id),
									MetaKeyAndValue.META_VALUE_Queue, taskID)
						except:
							pass
						
				elif (sample_id_a in request.GET):
					sample_id = request.GET[sample_id_a]
					sample = Sample.objects.get(pk=sample_id)
					
					default_project_software = DefaultProjectSoftware()
					default_project_software.set_default_software(software, request.user, Software.TYPE_OF_USE_sample,
										None, None, sample)
					## set a new default
					data['default'] = default_project_software.get_parameters(software.name, request.user, Software.TYPE_OF_USE_sample,\
													None, None, sample, software.technology.name)
					
					### need to re-run this sample with NanoFilt if the values change
					if (default_project_software.is_change_values_for_software(software.name, ConstantsSettings.TECHNOLOGY_illumina \
								if sample.is_type_fastq_gz_sequencing() else ConstantsSettings.TECHNOLOGY_minion)):
						b_change = True
						
						### re-run data
						manageDatabase = ManageDatabase()
						process_SGE = ProcessSGE()
						
						### change flag to nor finished
						sample.is_ready_for_projects = False
						sample.save()
		
						### create a task to perform the analysis of NanoFilt
						try:
							(job_name_wait, job_name) = request.user.profile.get_name_sge_seq(Profile.SGE_GLOBAL)
							if (sample.is_type_fastq_gz_sequencing()):
								taskID = process_SGE.set_run_trimmomatic_species(sample, request.user, job_name)
							else:
								taskID = process_SGE.set_run_clean_minion(sample, request.user, job_name)
								
							### set sample queue ID
							manageDatabase.set_sample_metakey(sample, sample.owner, MetaKeyAndValue.META_KEY_Queue_TaskID,
											MetaKeyAndValue.META_VALUE_Queue, taskID)
						except:
							pass
					
				else:
					default_software = DefaultSoftware()
					default_software.set_default_software(software)
					
					## set a new default
					data['default'] = default_software.get_parameters(software.name, request.user,
												software.technology.name)
					b_change = default_software.is_change_values_for_software(software)
					
				### message to show
				if b_change: data['message'] = "were set to default values."
				else: data['message'] = "already had the default values."
					
			except Software.DoesNotExist:
				return JsonResponse(data)
			except Project.DoesNotExist:
				return JsonResponse(data)
			data['is_ok'] = True
		return JsonResponse(data)

## @csrf_protect
@csrf_exempt
def mask_consensus(request):
	"""
	mask consensus
	"""
	if request.is_ajax():
		data = { 'is_ok' : False }
		
		## some pre-requisites
		if (not request.user.is_active or not request.user.is_authenticated): return JsonResponse(data)
		try:
			profile = Profile.objects.get(user__pk=request.user.pk)
		except Profile.DoesNotExist:
			return JsonResponse(data)
		if (profile.only_view_project): return JsonResponse(data)
		
		project_id_a = 'project_id'
		all_data_a = 'all_data'
		
		if (project_id_a in request.POST):
			
			project_id = request.POST[project_id_a]
			try:
				project = Project.objects.get(pk=project_id)
			except ProjectSample.DoesNotExist:
				return JsonResponse(data)

			manageDatabase = ManageDatabase()
			decode_masking_consensus = DecodeObjects()
			meta_value = manageDatabase.get_project_metakey_last(project, MetaKeyAndValue.META_KEY_Masking_consensus, MetaKeyAndValue.META_VALUE_Success)
			masking_consensus_original = None
			if not meta_value is None:
				masking_consensus_original = decode_masking_consensus.decode_result(meta_value.description)

			masking_consensus = decode_masking_consensus.decode_result(request.POST[all_data_a])
			masking_consensus.cleaning_mask_results()	## clean data
			
			if masking_consensus_original is None or masking_consensus_original != masking_consensus:
				manageDatabase.set_project_metakey(project, project.owner, MetaKeyAndValue.META_KEY_Masking_consensus,
							MetaKeyAndValue.META_VALUE_Success, masking_consensus.to_json())
				data['message'] = "The project '{}' is going to mask/unmask consensus. ".format(project.name)
				data['is_going_to_mask'] = True
			else: 
				data['message'] = "Masking regions are the same, nothing to do for project '{}'".format(project.name)
				data['is_going_to_mask'] = False
			
			data['new_title_i'] = masking_consensus.get_message_mask_to_show_in_web_site()
			data['new_class_i'] = "padding-button-table {} fa fa-superpowers padding-button-table tip".format(
					"warning_fa_icon" if masking_consensus.has_masking_data() else "")

			## check if they have projects			
			count = ProjectSample.objects.filter(project=project, is_deleted=False, is_error=False, is_finished=True).count()
			if count > 0:	### need to send a message to recalculate the global files
				metaKeyAndValue = MetaKeyAndValue()
				manageDatabase = ManageDatabase()
				try:
					process_SGE = ProcessSGE()
					taskID = process_SGE.set_collect_global_files(project, request.user)
					manageDatabase.set_project_metakey(project, request.user, metaKeyAndValue.get_meta_key(\
								MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id),
								MetaKeyAndValue.META_VALUE_Queue, taskID)
					data['is_ok'] = True
				except:
					data = { 'is_ok' : False }
			else:	### there's no projects to apply 
				data['is_ok'] = True
				data['message'] = "Mask/unmask consensus for the project '{}' are set.".format(project.name)
				
		return JsonResponse(data)


@csrf_protect
def get_mask_consensus_actual_values(request):
	"""
	return mask consensus of actual values
	"""
	data = { 'is_ok' : False }
	if request.is_ajax():
		
		## some pre-requisites
		if (not request.user.is_active or not request.user.is_authenticated): return JsonResponse(data)
		try:
			profile = Profile.objects.get(user__pk=request.user.pk)
		except Profile.DoesNotExist:
			return JsonResponse(data)
		if (profile.only_view_project): return JsonResponse(data)
		
		project_id_a = 'project_id'
		
		if (project_id_a in request.GET):
			
			project_id = request.GET[project_id_a]
			try:
				project = Project.objects.get(pk=project_id)
			except ProjectSample.DoesNotExist:
				return JsonResponse(data)

			manageDatabase = ManageDatabase()
			genetic_element = DecodeObjects()
			meta_value = manageDatabase.get_project_metakey_last(project, MetaKeyAndValue.META_KEY_Masking_consensus, MetaKeyAndValue.META_VALUE_Success)
			if meta_value is None:
				meta_value = manageDatabase.get_reference_metakey_last(project.reference, MetaKeyAndValue.META_KEY_Elements_And_CDS_Reference,
														MetaKeyAndValue.META_VALUE_Success)
				masking_consensus_original = genetic_element.decode_result(meta_value.description)
				for element in masking_consensus_original.get_sorted_elements():
					masking_consensus_original.dt_elements_mask[element] = MaskingConsensus()
			else:
				masking_consensus_original = genetic_element.decode_result(meta_value.description)
			
			### passing data
			data['all_data'] = masking_consensus_original.to_json()
			data['is_ok'] = True
		return JsonResponse(data)


@csrf_protect
def turn_on_off_software(request):
	"""
	Denies is_to_run in main software description.
	Don't do this if the software already run
	"""
	if request.is_ajax():
		data = { 'is_ok' : False }
		software_id_a = 'software_id'
		
		## some pre-requisites
		if (not request.user.is_active or not request.user.is_authenticated): return JsonResponse(data)
		try:
			profile = Profile.objects.get(user__pk=request.user.pk)
		except Profile.DoesNotExist:
			return JsonResponse(data)
		if (profile.only_view_project): return JsonResponse(data)
			
		if (software_id_a in request.GET):
			software_id = request.GET[software_id_a]
			try:
				software = Software.objects.get(pk=software_id)
				if software.can_be_on_off_in_pipeline:
					software.is_to_run = not software.is_to_run 
					software.save()
					
					## set a new default
					data['is_to_run'] = software.is_to_run
					data['message'] = "The '{}' in '{}' technology was turned '{}'.".format(
						software.name_extended, software.technology.name,
						"OFF" if software.is_to_run else "ON")
					
			except Software.DoesNotExist:
				return JsonResponse(data)
			except Project.DoesNotExist:
				return JsonResponse(data)
			data['is_ok'] = True
		return JsonResponse(data)

@csrf_protect
def get_software_name_to_turn_on_off(request):

	data = { 'is_ok' : False, 'message' : "You are not allow to do this operation." }
	
	if request.is_ajax():

		software_id_a = 'software_id'
		
		## some pre-requisites
		if (not request.user.is_active or not request.user.is_authenticated): return JsonResponse(data)
		try:
			profile = Profile.objects.get(user__pk=request.user.pk)
		except Profile.DoesNotExist:
			return JsonResponse(data)
		if (profile.only_view_project): return JsonResponse(data)
			
		if (software_id_a in request.GET):
			software_id = request.GET[software_id_a]
			try:
				software = Software.objects.get(pk=software_id)
				if software.can_be_on_off_in_pipeline:
					data['is_ok'] = True
					data['message'] = "Do you want to turn '{}' the '{}' in '{}' technology?.".format(
						"OFF" if software.is_to_run else "ON",
						software.name_extended, software.technology.name)
				else:
					data['message'] = "The '{}' in '{}' technology can not be ON/OFF.".format(software.name_extended, software.technology.name)
					
			except Software.DoesNotExist:
				return JsonResponse(data)
		return JsonResponse(data)
	
