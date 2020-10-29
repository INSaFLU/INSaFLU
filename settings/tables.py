'''
Created on 03/05/2020

@author: mmp
'''
import django_tables2 as tables

from django.utils.safestring import mark_safe
from constants.constants import Constants
from django.urls import reverse
from managing_files.models import ProjectSample
from settings.models import Software
from settings.default_software import DefaultSoftware
from settings.default_software_project_sample import DefaultProjectSoftware

class SoftwaresTable(tables.Table):
	
#   Renders a normal value as an internal hyperlink to another page.
#   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
	software = tables.Column('Software', empty_values=())
	version = tables.Column('Version', empty_values=())
	parameters = tables.Column('Parameters', empty_values=())
	options = tables.Column('Options', empty_values=())
	constants = Constants()

	def __init__(self, query_set, project = None, project_sample = None, b_enable_options = True):
		tables.Table.__init__(self, query_set)
		self.project = project
		self.project_sample = project_sample
		self.b_enable_options = b_enable_options
		
		self.count_project_sample = 0
		### get number of samples inside of this project, if project exist
		if (not project is None):
			self.count_project_sample = ProjectSample.objects.filter(project=project, is_deleted=False).count()

	class Meta:
		model = Software()
		fields = ('software', 'version', 'parameters', 'options')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There are no Softwares to show..."

	def render_software(self, record):
		return record.name if record.name_extended is None else record.name_extended

	def render_parameters(self, **kwargs):
		"""
		render parameters for the software
		"""
		from crequest.middleware import CrequestMiddleware
		current_request = CrequestMiddleware.get_request()
		user = current_request.user
		
		record = kwargs.pop("record")
		if (self.project is None and self.project_sample is None):
			default_software = DefaultSoftware()
			return default_software.get_parameters(record.name, user)
		elif (self.project_sample is None):
			default_software_projects = DefaultProjectSoftware()
			return default_software_projects.get_parameters(record.name, user, Software.TYPE_OF_USE_project, self.project, None)
		elif (self.project is None):
			default_software_projects = DefaultProjectSoftware()
			return default_software_projects.get_parameters(record.name, user, Software.TYPE_OF_USE_project_sample, None, self.project_sample)
		return ""

	def render_options(self, record):
		
		### if project 
		## Edit
		if (not self.project is None):
			str_links = '<a href=' + reverse('software-project-update', args=[record.pk, self.project.pk]) + ' data-toggle="tooltip" title="Edit parameters" ' +\
				'{}'.format("" if self.b_enable_options else 'onclick=\'return false;\' disable') + '><span><i class="fa fa-2x fa-pencil padding-button-table ' +\
				'{}'.format("" if self.b_enable_options else 'disable_fa_icon') + '"></i></span></a>'
			str_links += '<a href="{}"'.format('#id_set_default_modal' if self.b_enable_options else '') +\
				' id="id_default_parameter" data-toggle="modal" data-toggle="tooltip" title="Set default parameters"' +\
				'{}'.format("" if self.b_enable_options else 'onclick=\'return false;\' disable') +\
				' ref_name="' + record.name + '" pk="' + str(record.pk) + '" pk_proj="' + str(self.project.pk) +\
				'" type_software="{}'.format('software' if record.is_software() else 'INSaFLU') +\
				'" proj_name="' + str(self.project.name) + '"><span ><i class="fa fa-2x fa-power-off padding-button-table ' +\
				'{}'.format("" if self.b_enable_options else 'disable_fa_icon') + '"></i></span></a>'
		elif (not self.project_sample is None):
			str_links = '<a href=' + reverse('software-project-sample-update', args=[record.pk, self.project_sample.pk]) + ' data-toggle="tooltip" title="Edit parameters" ' +\
				'{}'.format("" if self.b_enable_options else 'onclick=\'return false;\' disable') + '><span><i class="fa fa-2x fa-pencil padding-button-table ' +\
				'{}'.format("" if self.b_enable_options else 'disable_fa_icon') + '"></i></span></a>'
				
			str_links += '<a href="{}"'.format('#id_set_default_modal' if self.b_enable_options else '') +\
				' id="id_default_parameter" data-toggle="modal" data-toggle="tooltip" title="Set default parameters"' +\
				'{}'.format("" if self.b_enable_options else 'onclick=\'return false;\' disable') +\
				' ref_name="' + record.name + '" pk="' + str(record.pk) + '" pk_proj_sample="' + str(self.project_sample.pk) +\
				'" type_software="{}'.format('software' if record.is_software() else 'INSaFLU') +\
				'" proj_name="' + str(self.project_sample.project.name) + '"><span ><i class="fa fa-2x fa-power-off padding-button-table ' +\
				'{}'.format("" if self.b_enable_options else 'disable_fa_icon') + '"></i></span></a>'
		else:
			str_links = '<a href=' + reverse('software-update', args=[record.pk]) + ' data-toggle="tooltip" title="Edit parameters" ' +\
				'><span><i class="fa fa-2x fa-pencil padding-button-table ' +\
				'"></i></span></a>'
				## Remove
			str_links += '<a href="{}"'.format('#id_set_default_modal' if self.b_enable_options else '') +\
				' id="id_default_parameter" data-toggle="modal" data-toggle="tooltip" title="Set default parameters"' +\
				' ref_name="' + record.name +\
				'" type_software="{}'.format('software' if record.is_software() else 'INSaFLU') +\
				'" pk="' + str(record.pk) + '"><span ><i class="fa fa-2x fa-power-off padding-button-table ' +\
				'"></i></span></a>'
		return mark_safe(str_links)

class INSaFLUParametersTable(tables.Table):
	
#   Renders a normal value as an internal hyperlink to another page.
#   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
	software = tables.Column('Software', empty_values=())
	parameters = tables.Column('Parameters', empty_values=())
	options = tables.Column('Options', empty_values=())
	constants = Constants()

	def __init__(self, query_set, project = None, project_sample = None, b_enable_options = True):
		tables.Table.__init__(self, query_set)
		self.project = project
		self.project_sample = project_sample
		self.b_enable_options = b_enable_options
		
		self.count_project_sample = 0
		### get number of samples inside of this project, if project exist
		if (not project is None):
			self.count_project_sample = ProjectSample.objects.filter(project=project, is_deleted=False).count()

	class Meta:
		model = Software()
		fields = ('software', 'parameters', 'options')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There are no INSaFLU parameters to show..."

	def render_software(self, record):
		return record.name if record.name_extended is None else record.name_extended

	def render_parameters(self, **kwargs):
		"""
		render parameters for the software
		"""
		from crequest.middleware import CrequestMiddleware
		current_request = CrequestMiddleware.get_request()
		user = current_request.user
		
		record = kwargs.pop("record")
		if (self.project is None and self.project_sample is None):
			default_software = DefaultSoftware()
			return default_software.get_parameters(record.name, user)
		elif (self.project_sample is None):
			default_software_projects = DefaultProjectSoftware()
			return default_software_projects.get_parameters(record.name, user, Software.TYPE_OF_USE_project, self.project, None)
		elif (self.project is None):
			default_software_projects = DefaultProjectSoftware()
			return default_software_projects.get_parameters(record.name, user, Software.TYPE_OF_USE_project_sample, None, self.project_sample)
		return ""

	def render_options(self, record):
		
		### if project 
		## Edit
		if (not self.project is None):
			str_links = '<a href=' + reverse('software-project-update', args=[record.pk, self.project.pk]) + ' data-toggle="tooltip" title="Edit parameters" ' +\
				'{}'.format("" if self.b_enable_options else 'onclick=\'return false;\' disable') + '><span><i class="fa fa-2x fa-pencil padding-button-table ' +\
				'{}'.format("" if self.b_enable_options else 'disable_fa_icon') + '"></i></span></a>'
			str_links += '<a href="{}"'.format('#id_set_default_modal' if self.b_enable_options else '') +\
				' id="id_default_parameter" data-toggle="modal" data-toggle="tooltip" title="Set default parameters"' +\
				'{}'.format("" if self.b_enable_options else 'onclick=\'return false;\' disable') +\
				' ref_name="' + record.name + '" pk="' + str(record.pk) + '" pk_proj="' + str(self.project.pk) +\
				'" type_software="{}'.format('software' if record.is_software() else 'INSaFLU') +\
				'" proj_name="' + str(self.project.name) + '"><span ><i class="fa fa-2x fa-power-off padding-button-table ' +\
				'{}'.format("" if self.b_enable_options else 'disable_fa_icon') + '"></i></span></a>'
		elif (not self.project_sample is None):
			str_links = '<a href=' + reverse('software-project-sample-update', args=[record.pk, self.project_sample.pk]) + ' data-toggle="tooltip" title="Edit parameters" ' +\
				'{}'.format("" if self.b_enable_options else 'onclick=\'return false;\' disable') + '><span><i class="fa fa-2x fa-pencil padding-button-table ' +\
				'{}'.format("" if self.b_enable_options else 'disable_fa_icon') + '"></i></span></a>'
				
			str_links += '<a href="{}"'.format('#id_set_default_modal' if self.b_enable_options else '') +\
				' id="id_default_parameter" data-toggle="modal" data-toggle="tooltip" title="Set default parameters"' +\
				'{}'.format("" if self.b_enable_options else 'onclick=\'return false;\' disable') +\
				' ref_name="' + record.name + '" pk="' + str(record.pk) + '" pk_proj_sample="' + str(self.project_sample.pk) +\
				'" type_software="{}'.format('software' if record.is_software() else 'INSaFLU') +\
				'" proj_name="' + str(self.project_sample.project.name) + '"><span ><i class="fa fa-2x fa-power-off padding-button-table ' +\
				'{}'.format("" if self.b_enable_options else 'disable_fa_icon') + '"></i></span></a>'
		else:
			str_links = '<a href=' + reverse('software-update', args=[record.pk]) + ' data-toggle="tooltip" title="Edit parameters" ' +\
				'><span><i class="fa fa-2x fa-pencil padding-button-table ' +\
				'"></i></span></a>'
				## Remove
			str_links += '<a href="{}"'.format('#id_set_default_modal' if self.b_enable_options else '') +\
				' id="id_default_parameter" data-toggle="modal" data-toggle="tooltip" title="Set default parameters"' +\
				' ref_name="' + record.name +\
				'" type_software="{}'.format('software' if record.is_software() else 'INSaFLU') +\
				'" pk="' + str(record.pk) + '"><span ><i class="fa fa-2x fa-power-off padding-button-table ' +\
				'"></i></span></a>'
		print(str_links)
		return mark_safe(str_links)




