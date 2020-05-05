'''
Created on Dec 6, 2017

@author: mmp
'''


from settings.models import Software
from settings.default_software import DefaultSoftware
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_protect
from extend_user.models import Profile

@csrf_protect
def set_default_parameters(request):
	"""
	remove a reference. It can only be removed if not belongs to any deleted project
	"""
	if request.is_ajax():
		data = { 'is_ok' : False }
		software_id_a = 'software_id'
		
		if (software_id_a in request.GET):
			
			## some pre-requisites
			if (not request.user.is_active or not request.user.is_authenticated): return JsonResponse(data)
			try:
				profile = Profile.objects.get(user__pk=request.user.pk)
			except Profile.DoesNotExist:
				return JsonResponse(data)
			if (profile.only_view_project): return JsonResponse(data)
			
			software_id = request.GET[software_id_a]
			try:
				software = Software.objects.get(pk=software_id)
				default_software = DefaultSoftware()
				default_software.set_default_software(software, request.user)
				## set a new default
				data['default'] = default_software.get_parameters(software.name, request.user)
			except Profile.DoesNotExist:
				return JsonResponse(data)
			data['is_ok'] = True
		return JsonResponse(data)

