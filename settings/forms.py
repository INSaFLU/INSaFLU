'''
Created on 04/05/2020

@author: mmp
'''
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Div, Layout, ButtonHolder, Submit, Button, Fieldset
from django.urls import reverse
from django import forms
from settings.models import Software, Parameter
from managing_files.models import Project, ProjectSample
from django.utils.html import escape

## https://kuanyui.github.io/2015/04/13/django-crispy-inline-form-layout-with-bootstrap/
class SoftwareForm(forms.ModelForm):
	"""
	Reference form, name, isolate_name and others
	"""
	error_css_class = 'error'
	
	class Meta:
		model = Software
		fields = ()
		
	def __init__(self, *args, **kwargs):
		self.request = kwargs.pop('request')
		## try to get project or project sample, for settings update
		pk_project = kwargs.get('pk_project')
		project = None
		if (not pk_project is None): 
			kwargs.pop('pk_project')
			project = Project.objects.get(pk=pk_project)

		project_sample = None
		pk_project_sample = kwargs.get('pk_project_sample')
		if (not pk_project_sample is None):
			kwargs.pop('pk_project_sample')
			project_sample = ProjectSample.objects.get(pk=pk_project_sample)
			
		## end
		super(SoftwareForm, self).__init__(*args, **kwargs)

		### return the parameters that is possible to change
		paramers = Parameter.objects.filter(software=self.instance, project=project, project_sample=project_sample)
		dt_fields = {}
		vect_divs = []
		for parameter in paramers:
			if (not parameter.can_change): continue
			if (parameter.is_integer()):
				dt_fields[parameter.get_unique_id()] = forms.IntegerField(max_value=int(parameter.range_max),\
							min_value=int(parameter.range_min), required = True)
				help_text = parameter.description + " Range: {}.".format(parameter.range_available)
				if (not parameter.not_set_value is None):
					help_text += " If value equal to {} this parameter is excluded.".format(parameter.not_set_value)
				dt_fields[parameter.get_unique_id()].help_text = escape(help_text) 

			elif (parameter.is_float()):
				dt_fields[parameter.get_unique_id()] = forms.FloatField(max_value=float(parameter.range_max),\
							min_value=float(parameter.range_min), required = True)
				help_text = parameter.description + " Range: {}.".format(parameter.range_available)
				if (not parameter.not_set_value is None):
					help_text += " If value equal to {} this parameter is excluded.".format(parameter.not_set_value)
				dt_fields[parameter.get_unique_id()].help_text = escape(help_text) 
			else:
				dt_fields[parameter.get_unique_id()] = forms.CharField(required = True)
				help_text = parameter.description + " Range: {}.".format(parameter.range_available)
				if (not parameter.not_set_value is None):
					help_text += " If value equal to {} this parameter is excluded.".format(parameter.not_set_value)
				dt_fields[parameter.get_unique_id()].help_text = escape(help_text)

			dt_fields[parameter.get_unique_id()].label = parameter.name
			dt_fields[parameter.get_unique_id()].initial = parameter.parameter
			vect_divs.append(
					Div(parameter.get_unique_id(),
						css_class = "col-sm-4"))
		### set all fields
		self.fields.update(dt_fields)

		vect_rows_divs = []
		for _ in range(0, len(vect_divs), 3):
			if (_ + 2 < len(vect_divs)): vect_rows_divs.append(Div(vect_divs[_], vect_divs[_ + 1], vect_divs[_ + 2], css_class = 'row'))
			elif (_ + 1 < len(vect_divs)): vect_rows_divs.append(Div(vect_divs[_], vect_divs[_ + 1], css_class = 'row'))
			else: vect_rows_divs.append(Div(vect_divs[_], css_class = 'row'))
		
		self.helper = FormHelper()
		self.helper.form_method = 'POST'
		
		#### message in form
		form_message = 'Update {} parameters for -{}-'.format("software" if self.instance.is_software() else "INSaFLU",\
					self.instance.name)
		if (not project is None):
			form_message = "Update {} parameters for -{}- project:'{}'".format(\
				"software" if self.instance.is_software() else "INSaFLU",\
				self.instance.name, project.name)
		if (not project_sample is None):
			form_message = "Update {} parameters for -{}- project:'{}' for sample:'{}'.".format(\
				"software" if self.instance.is_software() else "INSaFLU",\
				self.instance.name, project_sample.project.name, project_sample.sample.name)

		if (len(vect_rows_divs) == 1):
			self.helper.layout = Layout(
				Fieldset(
					form_message,
					vect_rows_divs[0],
					css_class = 'article-content'
				),
				ButtonHolder(
					Submit('save', 'Save', css_class='btn-primary'),
					Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(
						self._get_reverse(project, project_sample)))
				)
			)
		elif (len(vect_rows_divs) == 2):
			self.helper.layout = Layout(
				Fieldset(
					form_message,
					vect_rows_divs[0],
					vect_rows_divs[1],
					css_class = 'article-content'
				),
				ButtonHolder(
					Submit('save', 'Save', css_class='btn-primary'),
					Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(
						self._get_reverse(project, project_sample)))
				)
			)
		elif (len(vect_rows_divs) == 3):
			self.helper.layout = Layout(
				Fieldset(
					form_message,
					vect_rows_divs[0],
					vect_rows_divs[1],
					vect_rows_divs[2],
					css_class = 'article-content'
				),
				ButtonHolder(
					Submit('save', 'Save', css_class='btn-primary'),
					Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(
						self._get_reverse(project, project_sample)))
				)
			)
		elif (len(vect_rows_divs) == 4):
			self.helper.layout = Layout(
				Fieldset(
					form_message,
					vect_rows_divs[0],
					vect_rows_divs[1],
					vect_rows_divs[2],
					vect_rows_divs[3],
					css_class = 'article-content'
				),
				ButtonHolder(
					Submit('save', 'Save', css_class='btn-primary'),
					Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(
						self._get_reverse(project, project_sample)))
				)
			)
		
	def _get_reverse(self, project, project_sample):
		"""
		"""
		if (not project is None):
			return reverse('project-settings', args=[project.pk])
		if (not project_sample is None):
			return reverse('sample-project-settings', args=[project_sample.pk])
		return reverse('settings')


