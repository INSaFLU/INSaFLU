
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Div, Layout, ButtonHolder, Submit, Button, HTML, Fieldset
from django.conf import settings 
from django import forms
from django.urls import reverse
from django.utils.translation import ugettext_lazy as _
from django.core.files.temp import NamedTemporaryFile
from django.contrib.gis.geos import Point
from django.forms.models import inlineformset_factory
from django.forms import widgets
from django.utils.safestring import mark_safe
from django.core.exceptions import ValidationError
from utils.utils import Utils
from utils.constants import Constants
from managing_files.models import Reference, Sample, DataSet, VacineStatus, Project
import os

## https://kuanyui.github.io/2015/04/13/django-crispy-inline-form-layout-with-bootstrap/
class ReferenceForm(forms.ModelForm):
	"""
	Reference form, name, isolate_name and others
	"""
	utils = Utils()
	error_css_class = 'error'
	
	class Meta:
		model = Reference
		# specify what fields should be used in this form.
		fields = ('name', 'isolate_name', 'reference_fasta', 'reference_genbank')
		
	def __init__(self, *args, **kwargs):
		self.request = kwargs.pop('request')
		super(ReferenceForm, self).__init__(*args, **kwargs)

		## can exclude explicitly
		## exclude = ('md5',)
		field_text= [
			# (field_name, Field title label, Detailed field description, requiered)
			('name', 'Name', 'Regular name for this reference', True),
			('isolate_name', 'Isolate name', 'Isolate name for this reference', False),
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
				Div('isolate_name', css_class="col-sm-7"),
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
		name = cleaned_data['name']
		try:
			Reference.objects.get(name=name, owner__username=self.request.user.username)
			self.add_error('name', _("This name '" + name +"' already exist in database, please choose other."))
		except Reference.DoesNotExist:
			pass
		
		### testing file names
		reference_fasta = cleaned_data['reference_fasta']
		reference_genbank = cleaned_data['reference_genbank']
		if (reference_fasta.name == reference_genbank.name):
			self.add_error('reference_fasta', _("Error: both files has the same name. Please, different files."))
			self.add_error('reference_genbank', _("Error: both files has the same name. Please, different files."))
			return cleaned_data
		
		## testing fasta
		some_error_in_files = False
		reference_fasta_temp_file_name = NamedTemporaryFile(prefix='flu_fa_', delete=False)
		reference_fasta_temp_file_name.write(reference_fasta.file.read())
		reference_fasta_temp_file_name.flush()
		reference_fasta_temp_file_name.close()
		try:
			number_locus = self.utils.is_fasta(reference_fasta_temp_file_name.name)
			self.request.session[Constants.NUMBER_LOCUS_FASTA_FILE] = number_locus
		except IOError as e:	## (e.errno, e.strerror)
			os.unlink(reference_fasta_temp_file_name.name)
			some_error_in_files = True
			self.add_error('reference_fasta', e.args[0])
		
		### testing genbank
		reference_genbank_temp_file_name = NamedTemporaryFile(prefix='flu_gb_', delete=False)
		reference_genbank = cleaned_data['reference_genbank']
		reference_genbank_temp_file_name.write(reference_genbank.read())
		reference_genbank_temp_file_name.flush()
		reference_genbank_temp_file_name.close()
		try:
			self.utils.is_genbank(reference_genbank_temp_file_name.name)
		except IOError as e:
			some_error_in_files = True
			os.unlink(reference_genbank_temp_file_name.name)
			self.add_error('reference_genbank', e.args[0])
		
		## if some errors in the files, fasta or genBank, return
		if (some_error_in_files): return cleaned_data
		
		## test locus names and length of sequences
		try:
			self.utils.compare_locus_fasta_gb(reference_fasta_temp_file_name.name, reference_genbank_temp_file_name.name)
		except ValueError as e:
			self.add_error('reference_fasta', e.args[0])
			self.add_error('reference_genbank', e.args[0])
		
		## remove temp files
		os.unlink(reference_genbank_temp_file_name.name)
		os.unlink(reference_fasta_temp_file_name.name)
		return cleaned_data


class DataSetForm(forms.ModelForm):

	class Meta:
		model = DataSet
		fields = ('name',)
		
	def __init__(self, *args, **kwargs):
		self.request = kwargs.pop('request')
		self.fields['data_set'].queryset = DataSet.objects.filter(owner__id = self.request.user.id)

		field_text= [
			# (field_name, Field title label, Detailed field description, requiered)
			('data_set', 'Dataset', 'Define a specific dataset, can be used in the future to filter samples', True),
		]
		for x in field_text:
			self.fields[x[0]].label = x[1]
			self.fields[x[0]].help_text = x[2]
			self.fields[x[0]].required = x[3]

		self.helper = FormHelper()
		self.helper.form_method = 'POST'
		self.helper.layout = Layout(
			Div(
				Div('data_set', css_class="col-sm-3"),
				css_class = 'row'
			),
			ButtonHolder(
				Submit('save', 'Save', css_class='btn-primary'),
				Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(reverse('sample-add')))
			)
		)


class RelatedFieldWidgetCanAdd(widgets.Select):

	def __init__(self, related_model, related_url=None, *args, **kw):
	
		super(RelatedFieldWidgetCanAdd, self).__init__(*args, **kw)
	
		if not related_url:
			rel_to = related_model
			info = (rel_to._meta.app_label, rel_to._meta.object_name.lower())
			related_url = 'admin:%s_%s_add' % info
	
		# Be careful that here "reverse" is not allowed
		self.related_url = related_url
	
	def render(self, name, value, *args, **kwargs):
		self.related_url = reverse(self.related_url)
		output = [super(RelatedFieldWidgetCanAdd, self).render(name, value, *args, **kwargs)]
		output.append('<a href="%s" class="add-another" id="add_id_%s" onclick="return showAddAnotherPopup(this);"> ' % (self.related_url, name))
		output.append('<img src={% static \'admin/img/icon_addlink.gif\' %} width="10" height="10" alt="Add Another"/></a>')
		return mark_safe(''.join(output))


class DateInput(forms.DateInput):
	input_type = 'date'


## https://stackoverflow.com/questions/4497684/django-class-based-views-with-inline-model-form-or-formset
## https://gist.github.com/neara/6209563
## https://stackoverflow.com/questions/11198638/django-inline-formset-error/11199031
class SampleForm(forms.ModelForm):
	"""
	Sample form, name, generic dataset and files
	"""
	utils = Utils()
	
	DATE_CHOICES=[('date_of_onset','Onset date'),
		 ('date_of_collection','Collection date'),
		 ('date_of_receipt_lab','Lab reception date'),]
	
	## set the date format
	# 38.7223° N, 9.1393° 
	lat = forms.FloatField(required=True)
	lng = forms.FloatField(required=True)
	like_dates = forms.ChoiceField(choices=DATE_CHOICES, widget=forms.RadioSelect())
	error_css_class = 'error'

##	geo_local = PointField(widget=CustomPointWidget(), required=False, srid=4326)
	
	class Meta:
		model = Sample
		# specify what fields should be used in this form.
		fields = ('name', 'date_of_collection', 'date_of_onset', 'date_of_receipt_lab', 'vaccine_status', 'data_set', 'path_name_1', 'path_name_2')
		date_of_onset = forms.DateField(input_formats = settings.DATE_INPUT_FORMATS)
		date_of_collection = forms.DateField(input_formats = settings.DATE_INPUT_FORMATS)
		date_of_receipt_lab = forms.DateField(input_formats = settings.DATE_INPUT_FORMATS)
		widgets = {
			'date_of_onset': DateInput(attrs={'class':'datepicker', 'placeHolder' : 'dd/mm/yyyy'}),
			'date_of_collection': DateInput(attrs={'class':'datepicker', 'placeHolder' : 'dd/mm/yyyy'}),
			'date_of_receipt_lab': DateInput(attrs={'class':'datepicker', 'placeHolder' : 'dd/mm/yyyy'}),
		}

	def __init__(self, *args, **kwargs):
		self.request = kwargs.pop('request')
		super(SampleForm, self).__init__(*args, **kwargs)
		geo_local = self.initial.get("geo_local", None)
		if isinstance(geo_local, Point):
			self.initial["lng"], self.initial["lat"] = geo_local.tuple
		
		## define a specific query set
		self.fields['data_set'].queryset = DataSet.objects.filter(owner__id = self.request.user.id)
		self.fields['data_set'].empty_label = None
		self.fields['vaccine_status'].queryset = VacineStatus.objects.filter(owner__id = self.request.user.id)
		
		## can exclude explicitly
		## exclude = ('md5',)
		field_text= [
			# (field_name, Field title label, Detailed field description, requiered)
			('name', 'Name', 'Unique identify for this sample', True),
			('date_of_onset', 'Date of onset', 'Date of onset', False),
			('date_of_collection', 'Date of collection', 'Date of collection', False),
			('date_of_receipt_lab', 'Date on Lab', 'Date receipt on the lab', False),
			('vaccine_status', 'Vaccine status', 'Discrimination of vaccination status', False),
			('data_set', 'Dataset', 'Specific dataset, can be used to organize samples', False),
		##	('geo_local', 'Global position', 'Geo position where the sample was collected', False),
			('lat', 'Latitude', 'Geo position where the sample was collected', False),
			('lng', 'Longitude', 'Geo position where the sample was collected', False),
			('path_name_1', 'Raw fastq.gz (R1)', 'Raw file R1 with fastq gzip file (< 30MB)', True),
			('path_name_2', 'Raw fastq.gz (R2)', 'Raw file R2 with fastq gzip file (< 30MB)', False),
			('like_dates', 'Choose a date', 'Choose a date to collect the Day, Week and Year', False),
		]
		for x in field_text:
			self.fields[x[0]].label = x[1]
			self.fields[x[0]].help_text = x[2]
			self.fields[x[0]].required = x[3]

		self.helper = FormHelper()
		self.helper.form_method = 'POST'
		self.helper.layout = Layout(
			Div(
				Div('name', css_class="col-sm-4"),
				Div( 
					HTML('<a href="{% url "sample-dataset" %}" id="data_set_add_modal" <span ><i class="fa fa-plus-square"></i></span></a>'),
					Div('data_set'), 
					css_class = "row col-sm-4"),
				Div( 
					HTML('<a href="{% url "sample-vaccine" %}" id="vaccine_add_modal" <span ><i class="fa fa-plus-square"></i></span></a>'),
					Div('vaccine_status'), 
					css_class = "row col-sm-4"),
				css_class = 'row'
			),
			Fieldset(
				'Dates',
				Div(
					Div('like_dates', css_class="col-sm-3"),
					Div('date_of_onset', css_class="col-sm-3"),
					Div('date_of_collection', css_class="col-sm-3"),
					Div('date_of_receipt_lab', css_class="col-sm-3"),
					css_class = 'row '),
				css_class = 'rounded-group-box'
			),
			Fieldset(
				'Global position',
				Div(
					Div('lat', css_class="col-sm-4"),
					Div('lng', css_class="col-sm-4"),
					css_class = 'row'
				),
			),
			Div('path_name_1', css_class = 'show-for-sr'),
			Div('path_name_2', css_class = 'show-for-sr'),
# 			Div(
# 				Div('geo_local', css_class="col-sm-8"),
# 				css_class = 'row'
# 			),
			ButtonHolder(
				Submit('save', 'Save', css_class='btn-primary'),
				Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(reverse('samples')))
			)
		)

	def clean(self):
		"""
		Clean all together because it's necessary to compare the genbank and fasta files
		"""
		cleaned_data = super(SampleForm, self).clean()
		name = self.cleaned_data['name']
		try:
			Sample.objects.get(name=name, owner__username=self.request.user.username)
			self.add_error('name', ValidationError(_("This name '" + name +"' already exist in database, please choose other."), code='invalid'))
		except Sample.DoesNotExist:
			pass
		
		## latitude in degrees is -90 and +90 for the southern and northern hemisphere respectively. Longitude is in the range -180 and +180
		lat = self.cleaned_data.get('lat')
		if (lat != None and (lat > 90 or lat < -90)): self.add_error('lat', _("Latitute must have values between +90 and -90."))
		lng = self.cleaned_data.get('lng')
		if (lng != None and (lng > 180 or lng < -180)): self.add_error('lng', _("Longitude must have values between +180 and -180."))
		
		### testing file names
		path_name_1 = self.cleaned_data.get('path_name_1')
		path_name_2 = self.cleaned_data.get('path_name_2')
		if (path_name_2 != None and path_name_1.name == path_name_2.name):
			self.add_error('path_name_1', _("Error: both files has the same name. Please, different files."))
			self.add_error('path_name_2', _("Error: both files has the same name. Please, different files."))
			return cleaned_data
		
		## testing fastq
		fastaq_temp_file_name = NamedTemporaryFile(prefix='flu_fq_', suffix='.fastq.gz', delete=False)
		fastaq_temp_file_name.write(path_name_1.file.read())
		fastaq_temp_file_name.flush()
		fastaq_temp_file_name.close()
		
		like_dates = self.cleaned_data.get('like_dates')
		date_of_onset = self.cleaned_data.get('date_of_onset')
		date_of_collection = self.cleaned_data.get('date_of_collection')
		date_of_receipt_lab = self.cleaned_data.get('date_of_receipt_lab')
		
		#####
		if (like_dates == None and date_of_onset != None and date_of_collection != None and date_of_receipt_lab != None):
			self.add_error('like_dates', _("Please, choose a data to collect the day, week and year."))
		elif (like_dates == 'date_of_onset' and date_of_onset == None):
			self.add_error('like_dates', _("Error, the Onset date is null."))
		elif (like_dates == 'date_of_collection' and date_of_collection == None):
			self.add_error('date_of_collection', _("Error, the Collection date is null."))
		elif (like_dates == 'date_of_receipt_lab' and date_of_receipt_lab == None):
			self.add_error('date_of_receipt_lab', _("Error, the Lab Receipt date is null."))
			
		try:
			self.utils.is_fastq_gz(fastaq_temp_file_name.name)
		except Exception as e:	## (e.errno, e.strerror)
			os.unlink(fastaq_temp_file_name.name)
			self.add_error('path_name_1', e.args[0])
			return cleaned_data
		
		## testing fastq
		if (path_name_2 != None):
			fastaq_temp_file_name_2 = NamedTemporaryFile(prefix='flu_fq_', suffix='.fastq.gz', delete=False)
			fastaq_temp_file_name_2.write(path_name_2.file.read())
			fastaq_temp_file_name_2.flush()
			fastaq_temp_file_name_2.close()
			
			try:
				self.utils.is_fastq_gz(fastaq_temp_file_name_2.name)
			except Exception as e:	## (e.errno, e.strerror)
				os.unlink(fastaq_temp_file_name_2.name)
				self.add_error('path_name_1', e.args[0])
				return cleaned_data

		## remove temp files
		os.unlink(fastaq_temp_file_name.name)
		if (path_name_2 != None): os.unlink(fastaq_temp_file_name_2.name)
		return cleaned_data
		

class ReferenceProjectForm(forms.ModelForm):
	"""
	Create a new project 
	"""
	error_css_class = 'error'
	
	class Meta:
		model = Project
		exclude = ()

	def __init__(self, *args, **kwargs):
		super(ReferenceProjectForm, self).__init__(*args, **kwargs)
		
		field_text= [
			('name', 'Name', 'Unique identify for this project', True),
		]
		for x in field_text:
			self.fields[x[0]].label = x[1]
			self.fields[x[0]].help_text = x[2]
			self.fields[x[0]].required = x[3]

	def clean(self):
		"""
		Clean all together because it's necessary to compare the genbank and fasta files
		"""
		cleaned_data = super(ReferenceProjectForm, self).clean()
		name = cleaned_data['name']
		
		return cleaned_data
	

ReferenceProjectFormSet = inlineformset_factory(Reference, Project, form=ReferenceProjectForm, extra=1)


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
# 		utils = Constants()
# 		try:
# 			number_locus = utils.is_fasta(reference_fasta_temp_file_name.name)
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
# 		utils = Constants()
# 		try:
# 			utils.is_genbank(reference_genbank_temp_file_name.name)
# 			os.unlink(reference_genbank_temp_file_name.name)
# 		except IOError as e:
# 			os.unlink(reference_genbank_temp_file_name.name)
# 			raise forms.ValidationError(e.args[0])
# 		return reference_genbank
	
