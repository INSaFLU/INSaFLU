
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, ButtonHolder, Submit
from crispy_forms.layout import Div, Layout, ButtonHolder, Submit, Button
from django import forms
from django.urls import reverse
from .models import Reference
from django.utils.translation import ugettext_lazy as _
from django.core.files.temp import NamedTemporaryFile
from constants.Constants import Constants
import os

## https://kuanyui.github.io/2015/04/13/django-crispy-inline-form-layout-with-bootstrap/
class ReferenceForm(forms.ModelForm):
	class Meta:
		model = Reference
		# specify what fields should be used in this form.
		fields = ('name', 'scientific_name', 'reference_fasta', 'reference_genbank')
		
	def __init__(self, *args, **kwargs):
		self.request = kwargs.pop('request')
		super(ReferenceForm, self).__init__(*args, **kwargs)

		## can exclude explicitly
		## exclude = ('md5',)
		field_text= [
			# (field_name, Field title label, Detailed field description, requiered)
			('name', 'Name', 'Regular name for this reference', True),
			('scientific_name', 'Scientific name', 'Scientific name for this reference', False),
			('reference_fasta', 'Reference (fasta)', 'Reference file in fasta format', True),
			('reference_genbank', 'Reference (genBank)', 'Reference file in genBank format, all locus must have the same name of fasta locus file', True),
		]
		for x in field_text:
			self.fields[x[0]].label = x[1]
			self.fields[x[0]].help_text = x[2]
			self.fields[x[0]].required = x[3]

## in case you have to undo to specific ID
##		cancel_url = reverse('references')
##		if self.instance: cancel_url = reverse('references', kwargs={'pk': self.instance.id})

		self.helper = FormHelper()
		self.helper.form_method = 'POST'
		self.helper.layout = Layout(
			Div(
				Div('name', css_class="col-sm-3"),
				Div('scientific_name', css_class="col-sm-7"),
				css_class = 'row'
			),
			Div('reference_fasta', css_class = 'show-for-sr'),
			Div('reference_genbank', css_class = 'show-for-sr'),
			ButtonHolder(
				Submit('save', 'Save', css_class='btn-primary'),
				Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(reverse('references')))
			)
		)

	def clean(self):
		"""
		Clean all together because it's necessary to compare the genbank and fasta files
		"""
		cleaned_data = super(ReferenceForm, self).clean()
		name = self.cleaned_data['name']
		try:
			Reference.objects.get(name=name, owner__username=self.request.user.username)
			self.add_error('name', _("This name '" + name +"' already exist in database, please choose other."))
		except Reference.DoesNotExist:
			pass
		
		### testing fasta
		reference_fasta = self.cleaned_data['reference_fasta']
		reference_genbank = self.cleaned_data['reference_genbank']
		if (reference_fasta.name == reference_genbank.name):
			self.add_error('reference_fasta', _("Error: both files has the same name. Please, two different files."))
			self.add_error('reference_genbank', _("Error: both files has the same name. Please, two different files."))
			return cleaned_data
		
		## testing fasta
		some_error_in_files = False
		reference_fasta_temp_file_name = NamedTemporaryFile(prefix='flu_fa_', delete=False)
		reference_fasta_temp_file_name.write(reference_fasta.file.read())
		reference_fasta_temp_file_name.flush()
		reference_fasta_temp_file_name.close()
		constants = Constants()
		try:
			number_locus = constants.is_fasta(reference_fasta_temp_file_name.name)
			self.request.session[Constants.NUMBER_LOCUS_FASTA_FILE] = number_locus
		except IOError as e:	## (e.errno, e.strerror)
			os.unlink(reference_fasta_temp_file_name.name)
			some_error_in_files = True
			self.add_error('reference_fasta', e.args[0])
		
		### testing genbank
		reference_genbank_temp_file_name = NamedTemporaryFile(prefix='flu_gb_', delete=False)
		reference_genbank = self.cleaned_data['reference_genbank']
		reference_genbank_temp_file_name.write(reference_genbank.read())
		reference_genbank_temp_file_name.flush()
		reference_genbank_temp_file_name.close()
		constants = Constants()
		try:
			constants.is_genbank(reference_genbank_temp_file_name.name)
		except IOError as e:
			some_error_in_files = True
			os.unlink(reference_genbank_temp_file_name.name)
			self.add_error('reference_genbank', e.args[0])
		
		## if some errors in the files, fasta or genBank, return
		if (some_error_in_files): return cleaned_data
		
		## test locus names and length of sequences
		try:
			constants.compare_locus_fasta_gb(reference_fasta_temp_file_name.name, reference_genbank_temp_file_name.name)
		except ValueError as e:
			self.add_error('reference_fasta', e.args[0])
			self.add_error('reference_genbank', e.args[0])
		
		## remove temp files
		os.unlink(reference_genbank_temp_file_name.name)
		os.unlink(reference_fasta_temp_file_name.name)
		return cleaned_data

		
# 	def clean_name(self):
# 		"""
# 		Clean only the name
# 		"""
# 		name = self.cleaned_data['name']
# 		try:
# 			Reference.objects.get(name=name, owner__username=self.request.user.username)
# 			raise forms.ValidationError(_("This name '" + name +"' already exist in database, please choose other."))
# 		except Reference.DoesNotExist:
# 			pass
# 		return name
#	
# 	def clean_reference_fasta(self):
# 		"""
# 		Test Fasta file
# 		"""
# 		reference_fasta = self.cleaned_data['reference_fasta']
# 		#reference_genbank = self.cleaned_data['reference_genbank']
# 		#if (reference_fasta.name == reference_genbank): raise forms.ValidationError("Error: both files has the same name. Please two different files.")
# 		
# 		reference_fasta_temp_file_name = NamedTemporaryFile(prefix='flu_fa_', delete=False)
# 		reference_fasta_temp_file_name.write(reference_fasta.file.read())
# 		reference_fasta_temp_file_name.flush()
# 		reference_fasta_temp_file_name.close()
# 		constants = Constants()
# 		try:
# 			number_locus = constants.is_fasta(reference_fasta_temp_file_name.name)
# 			self.request.session[Constants.NUMBER_LOCUS_FASTA_FILE] = number_locus
# 			os.unlink(reference_fasta_temp_file_name.name)
# 		except IOError as e:	## (e.errno, e.strerror)
# 			os.unlink(reference_fasta_temp_file_name.name)
# 			raise forms.ValidationError(e.args[0])
# 		return reference_fasta
# 	
# 	def clean_reference_genbank(self):
# 		"""
# 		Test genBank file
# 		"""
# 		reference_genbank_temp_file_name = NamedTemporaryFile(prefix='flu_gb_', delete=False)
# 		reference_genbank = self.cleaned_data['reference_genbank']
# 		reference_genbank_temp_file_name.write(reference_genbank.read())
# 		reference_genbank_temp_file_name.flush()
# 		reference_genbank_temp_file_name.close()
# 		constants = Constants()
# 		try:
# 			constants.is_genbank(reference_genbank_temp_file_name.name)
# 			os.unlink(reference_genbank_temp_file_name.name)
# 		except IOError as e:
# 			os.unlink(reference_genbank_temp_file_name.name)
# 			raise forms.ValidationError(e.args[0])
# 		return reference_genbank
	
