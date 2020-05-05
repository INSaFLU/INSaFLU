'''
Created on 04/05/2020

@author: mmp
'''
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Div, Layout, ButtonHolder, Submit, Button, Fieldset
from django.urls import reverse
from settings.models import Software, Parameter
from django import forms

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
		super(SoftwareForm, self).__init__(*args, **kwargs)

		paramers = Parameter.objects.filter(software=self.instance)
		dt_fields = {}
		vect_divs = []
		for parameter in paramers:
			if (not parameter.can_change): continue
			if (parameter.is_integer()):
				dt_fields[parameter.get_unique_id()] = forms.IntegerField(max_value=int(parameter.range_max),\
							min_value=int(parameter.range_min), required = True)
				dt_fields[parameter.get_unique_id()].help_text = parameter.description + " Range: {}.".format(parameter.range_available)
				if (not parameter.not_set_value is None):
					dt_fields[parameter.get_unique_id()].help_text += " If value equal to {} this parameter is excluded.".format(parameter.not_set_value)

			elif (parameter.is_float()):
				dt_fields[parameter.get_unique_id()] = forms.FloatField(max_value=int(parameter.range_max),\
							min_value=int(parameter.range_min), required = True)
				dt_fields[parameter.get_unique_id()].help_text = parameter.description + " Range: {}.".format(parameter.range_available)
				if (not parameter.not_set_value is None):
					dt_fields[parameter.get_unique_id()].help_text += " If value equal to {} this parameter is excluded.".format(parameter.not_set_value)

			else:
				dt_fields[parameter.get_unique_id()] = forms.CharField(required = True)
				dt_fields[parameter.get_unique_id()].help_text = parameter.description + " Range: {}.".format(parameter.range_available)
				if (not parameter.not_set_value is None):
					dt_fields[parameter.get_unique_id()].help_text += " If value equal to {} this parameter is excluded.".format(parameter.not_set_value)

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
		self.helper.layout = Layout(
			Fieldset(
				'Update software parameters for -{}-'.format(self.instance.name),
				## TODO need to improve
				vect_rows_divs[0],
				vect_rows_divs[1],
				vect_rows_divs[2],
				css_class = 'article-content'
			),
			ButtonHolder(
				Submit('save', 'Save', css_class='btn-primary'),
				Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(reverse('settings')))
			)
		)
		


