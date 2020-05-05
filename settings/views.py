from django.views.generic import ListView, UpdateView
from braces.views import LoginRequiredMixin
from settings.models import Software, Parameter
from settings.default_software import DefaultSoftware
from settings.tables import SoftwaresTable
from settings.forms import SoftwareForm
from utils.utils import ShowInfoMainPage
from django.contrib import messages
from django.urls import reverse_lazy
from django.db import transaction


# Create your views here.
class SettingsView(LoginRequiredMixin, ListView):
	"""
	Home page
	"""
	model = Software
	template_name = 'settings/settings.html'
	context_object_name = 'software'
	
	def get_context_data(self, **kwargs):
		context = super(SettingsView, self).get_context_data(**kwargs)
		
		### test all defaults first, if exist in database
		default_software = DefaultSoftware()
		default_software.test_all_defaults(self.request.user) ## the user can have defaults yet
		
		query_set = Software.objects.filter(owner=self.request.user)
		table = SoftwaresTable(query_set)
		context['nav_settings'] = True
		context['table'] = table	
		context['show_paginatior'] = False
		context['show_info_main_page'] = ShowInfoMainPage()		## show main information about the institute	
		return context

	
class UpdateParametersView(LoginRequiredMixin, UpdateView):
	model = Software
	form_class = SoftwareForm
	success_url = reverse_lazy('settings')
	template_name = 'settings/software_update.html'

	## Other solution to get the reference
	## https://pypi.python.org/pypi?%3aaction=display&name=django-contrib-requestprovider&version=1.0.1
	def get_form_kwargs(self):
		"""
		Set the request to pass in the form
		"""
		kw = super(UpdateParametersView, self).get_form_kwargs()
		kw['request'] = self.request # the trick!
		return kw
	
	def get_context_data(self, **kwargs):
		context = super(UpdateParametersView, self).get_context_data(**kwargs)
		
		context['error_cant_see'] = self.request.user != context['software'].owner
		context['nav_settings'] = True
		context['nav_modal'] = True	## short the size of modal window
		context['show_info_main_page'] = ShowInfoMainPage()		## show main information about the institute	
		return context
	
	def form_valid(self, form):
		"""
		form update 
		"""

		## save it...
		with transaction.atomic():
			software = form.save(commit=False)
			paramers = Parameter.objects.filter(software=software)
			
			for parameter in paramers:
				if (not parameter.can_change): continue
				if (parameter.get_unique_id() in form.cleaned_data):
					parameter.parameter = "{}".format(form.cleaned_data[parameter.get_unique_id()])
					parameter.save()
			
		messages.success(self.request, "Software '" + software.name + "' parameters was updated successfully", fail_silently=True)
		return super(UpdateParametersView, self).form_valid(form)


	## static method, not need for now.
	form_valid_message = ""		## need to have this


