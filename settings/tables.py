'''
Created on 03/05/2020

@author: mmp
'''
import django_tables2 as tables
from settings.models import Software
from django.utils.safestring import mark_safe
from constants.constants import Constants
from django.urls import reverse
from settings.default_software import DefaultSoftware

class SoftwaresTable(tables.Table):
#   Renders a normal value as an internal hyperlink to another page.
#   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
	software = tables.Column('Software', empty_values=())
	version = tables.Column('Version', empty_values=())
	parameters = tables.Column('Parameters', empty_values=())
	options = tables.Column('Options', empty_values=())
	constants = Constants()
	
	class Meta:
		model = Software()
		fields = ('software', 'version', 'parameters', 'options')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There are no Softwares to show..."

	def render_software(self, record):
		return record.name
	
	def render_parameters(self, **kwargs):
		"""
		render parameters for the software
		"""
		from crequest.middleware import CrequestMiddleware
		current_request = CrequestMiddleware.get_request()
		user = current_request.user
		
		record = kwargs.pop("record")
		default_software = DefaultSoftware()
		return default_software.get_parameters(record.name, user)

	def render_options(self, record):
		## Edit
		str_links = '<a href=' + reverse('software-update', args=[record.pk]) + ' data-toggle="tooltip" title="Edit parameters"><span><i class="fa fa-2x fa-pencil padding-button-table"></i></span></a>'
		## Remove
		str_links += '<a href="#id_set_default_modal" id="id_default_parameter" data-toggle="modal" data-toggle="tooltip" title="Set default parameters"' +\
 					' ref_name="' + record.name + '" pk="' + str(record.pk) + '"><span ><i class="fa fa-2x fa-power-off padding-button-table">' +\
					'</i></span> </a>'
		return mark_safe(str_links)




