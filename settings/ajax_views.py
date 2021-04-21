'''
Created on Dec 6, 2017

@author: mmp
'''


from settings.models import Software
from settings.default_software import DefaultSoftware
from settings.default_software_project_sample import DefaultProjectSoftware
from managing_files.models import Project, ProjectSample, Sample
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_protect
from extend_user.models import Profile
from constants.meta_key_and_values import MetaKeyAndValue
from managing_files.manage_database import ManageDatabase
from utils.process_SGE import ProcessSGE
from constants.software_names import SoftwareNames

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
			try:
				software = Software.objects.get(pk=software_id)
				technology_name = SoftwareNames.TECHNOLOGY_illumina if software.technology is None else\
								software.technology.name
				if (project_id_a in request.GET):
					project_id = request.GET[project_id_a]
					project = Project.objects.get(pk=project_id)
					
					default_project_software = DefaultProjectSoftware()
					default_project_software.set_default_software(software, request.user, Software.TYPE_OF_USE_project, project,
											None, None)
					## set a new default
					data['default'] = default_project_software.get_parameters(software.name, request.user, Software.TYPE_OF_USE_project,
													project, None, None, technology_name)
				elif (project_sample_id_a in request.GET):
					project_sample_id = request.GET[project_sample_id_a]
					project_sample = ProjectSample.objects.get(pk=project_sample_id)
					
					default_project_software = DefaultProjectSoftware()
					default_project_software.set_default_software(software, request.user, Software.TYPE_OF_USE_project_sample,
										None, project_sample, None)
					## set a new default
					data['default'] = default_project_software.get_parameters(software.name, request.user, Software.TYPE_OF_USE_project_sample,\
													None, project_sample, None, technology_name)
					
					### need to re-run this sample with snippy if the values change
					if (default_project_software.is_change_values_for_software(software.name, SoftwareNames.TECHNOLOGY_illumina \
								if project_sample.is_sample_illumina() else SoftwareNames.TECHNOLOGY_minion)):
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
													None, None, sample, technology_name)
					
					### need to re-run this sample with NanoFilt if the values change
					if (default_project_software.is_change_values_for_software(software.name, SoftwareNames.TECHNOLOGY_illumina \
								if sample.is_type_fastq_gz_sequencing() else SoftwareNames.TECHNOLOGY_minion)):
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
					default_software.set_default_software(software, request.user)
					
					
					## set a new default
					data['default'] = default_software.get_parameters(software.name, request.user,
												technology_name)
					
			except Software.DoesNotExist:
				return JsonResponse(data)
			except Project.DoesNotExist:
				return JsonResponse(data)
			data['is_ok'] = True
		return JsonResponse(data)

