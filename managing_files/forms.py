
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
from utils.parse_in_files import ParseInFiles
from utils.software import Software
from constants.constants import Constants, TypeFile
from managing_files.models import Reference, Sample, DataSet, VaccineStatus, Project, UploadFiles
import os, re, logging, humanfriendly

## https://kuanyui.github.io/2015/04/13/django-crispy-inline-form-layout-with-bootstrap/
class ReferenceForm(forms.ModelForm):
	"""
	Reference form, name, isolate_name and others
	"""
	utils = Utils()
	software = Software()
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
			('reference_genbank', 'Reference (genBank)', 
					"""Reference file in genBank format.<br>
					Locus designations in the multi-GenBank file must have the same name as in the respective fasta file.<br>
					If you do not upload a Genbank file, INSaFLU will annotate the upload fasta file for you""", False),
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
				Submit('save', 'Save', css_class='btn-primary', onclick="this.disabled=true,this.form.submit();"),
				Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(reverse('references')))
			)
		)

	def clean(self):
		"""
		Clean all together because it's necessary to compare the genbank and fasta files
		"""
		cleaned_data = super(ReferenceForm, self).clean()
		name = self.cleaned_data.get('name', '').strip()
		if (len(name) == 0):
			self.add_error('name', _("Error: You must give a unique name."))
			return cleaned_data

		try:
			Reference.objects.get(name__iexact=name, owner=self.request.user, is_obsolete=False, is_deleted=False)
			self.add_error('name', _("This name '" + name +"' already exist in database, please choose other."))
		except Reference.DoesNotExist:
			pass
		
		## test reference_fasta
		if ('reference_fasta' not in cleaned_data):
			self.add_error('reference_fasta', _("Error: Must have a file."))
			return cleaned_data
		
		### testing file names
		reference_fasta = cleaned_data['reference_fasta']
		reference_genbank = cleaned_data['reference_genbank']
		if (reference_genbank != None and reference_fasta.name == reference_genbank.name):
			self.add_error('reference_fasta', _("Error: both files has the same name. Please, different files."))
			self.add_error('reference_genbank', _("Error: both files has the same name. Please, different files."))
			return cleaned_data
		
		## testing fasta
		some_error_in_files = False
		reference_fasta_temp_file_name = NamedTemporaryFile(prefix='flu_fa_', delete=False)
		reference_fasta_temp_file_name.write(reference_fasta.read())
		reference_fasta_temp_file_name.flush()
		reference_fasta_temp_file_name.close()
		self.software.dos_2_unix(reference_fasta_temp_file_name.name)
		try:
			number_locus = self.utils.is_fasta(reference_fasta_temp_file_name.name)
			self.request.session[Constants.NUMBER_LOCUS_FASTA_FILE] = number_locus
			
			## test the max numbers
			if (number_locus > Constants.MAX_SEQUENCES_FROM_FASTA):
				self.add_error('reference_fasta', _('Max allow number of sequences in fasta: {}'.format(Constants.MAX_SEQUENCES_FROM_FASTA)))
				some_error_in_files = True
			total_length_fasta = self.utils.get_total_length_fasta(reference_fasta_temp_file_name.name)
			if (not some_error_in_files and total_length_fasta > Constants.MAX_LENGTH_SEQUENCE_TOTAL_FROM_FASTA):
				some_error_in_files = True
				self.add_error('reference_fasta', _('The length sum of the sequences in fasta: {}'.format(Constants.MAX_LENGTH_SEQUENCE_TOTAL_FROM_FASTA)))
			
			n_seq_name_bigger_than = self.utils.get_number_seqs_names_bigger_than(reference_fasta_temp_file_name.name, Constants.MAX_LENGTH_SEQ_NAME)
			if (not some_error_in_files and n_seq_name_bigger_than > 0):
				some_error_in_files = True
				if (n_seq_name_bigger_than == 1):
					self.add_error('reference_fasta', _('There is one sequence name length bigger than {0}. The max. length name is {0}. Prokka constrainments.'.format(Constants.MAX_LENGTH_SEQ_NAME)))
				else:
					self.add_error('reference_fasta', _('There are {0} sequences with name length bigger than {1}. The max. length name is {1}. Prokka constrainments.'.format(n_seq_name_bigger_than, Constants.MAX_LENGTH_SEQ_NAME)))
					
				## if some errors in the files, fasta or genBank, return
				if (some_error_in_files): return cleaned_data
			
			if (not self.utils.test_sequences_same_length(reference_fasta_temp_file_name.name)):
				self.add_error('reference_fasta', _('There are sequences that have not the same length. This produce errors for samtools faidx.'))
				return cleaned_data
		
		except IOError as e:	## (e.errno, e.strerror)
			os.unlink(reference_fasta_temp_file_name.name)
			some_error_in_files = True
			self.add_error('reference_fasta', e.args[0])
		except: 
			os.unlink(reference_fasta_temp_file_name.name)
			some_error_in_files = True
			self.add_error('reference_fasta', "Not a valid 'fasta' file.")
			
		### test if it has degenerated bases
		if (os.path.exists(reference_fasta_temp_file_name.name)):
			try:
				self.utils.has_degenerated_bases(reference_fasta_temp_file_name.name)
			except Exception as e:
				os.unlink(reference_fasta_temp_file_name.name)
				some_error_in_files = True
				self.add_error('reference_fasta', e.args[0])
			
		### testing genbank
		reference_genbank_temp_file_name = NamedTemporaryFile(prefix='flu_gb_', delete=False)
		reference_genbank = cleaned_data['reference_genbank']
		if (reference_genbank != None):
			reference_genbank_temp_file_name.write(reference_genbank.read())
			reference_genbank_temp_file_name.flush()
			reference_genbank_temp_file_name.close()
			self.software.dos_2_unix(reference_genbank_temp_file_name.name)
			try:
				self.utils.is_genbank(reference_genbank_temp_file_name.name)
			except IOError as e:
				some_error_in_files = True
				os.unlink(reference_genbank_temp_file_name.name)
				self.add_error('reference_genbank', e.args[0])
			except: 
				os.unlink(reference_genbank_temp_file_name.name)
				some_error_in_files = True
				self.add_error('reference_genbank', "Not a valid 'genbank' file.")
			
		## if some errors in the files, fasta or genBank, return
		if (some_error_in_files): return cleaned_data
		
		## test locus names and length of sequences
		if (reference_genbank != None):
			try:
				self.utils.compare_locus_fasta_gb(reference_fasta_temp_file_name.name, reference_genbank_temp_file_name.name)
			except ValueError as e:
				self.add_error('reference_fasta', e.args[0])
				self.add_error('reference_genbank', e.args[0])
		
		## remove temp files
		os.unlink(reference_genbank_temp_file_name.name)
		os.unlink(reference_fasta_temp_file_name.name)
		return cleaned_data


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
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	class Meta:
		model = Sample
		# specify what fields should be used in this form.
		fields = ('name', 'vaccine_status', 'data_set', 'date_of_onset', 'date_of_collection', 'date_of_receipt_lab', 'path_name_1', 'path_name_2')

	def __init__(self, *args, **kwargs):
		self.request = kwargs.pop('request')
		super(SampleForm, self).__init__(*args, **kwargs)
		geo_local = self.initial.get("geo_local", None)
		if isinstance(geo_local, Point):
			self.initial["lng"], self.initial["lat"] = geo_local.tuple
		
		## define a specific query set
		self.fields['data_set'].queryset = DataSet.objects.filter(owner__id = self.request.user.id)
		self.fields['data_set'].empty_label = None		## to remove empty label in combo box
		self.fields['vaccine_status'].queryset = VaccineStatus.objects.filter(owner__id = self.request.user.id)
		
		##// <input placeholder="dd/mm/yyyy" name="date_of_onset" id="id_date_of_onset"/>
		self.fields['date_of_onset'].widget.attrs.update({
			'placeholder': 'dd/mm/yyyy',
		})
		self.fields['date_of_receipt_lab'].widget.attrs.update({
			'placeholder': 'dd/mm/yyyy',
		})
		self.fields['date_of_collection'].widget.attrs.update({
			'placeholder': 'dd/mm/yyyy',
		})
		
		#placeholder="dd/mm/yyyy"
		
		## can exclude explicitly
		## exclude = ('md5',)
		
		#### mesages
		if (settings.DOWN_SIZE_FASTQ_FILES):
			message_r1 = "Max raw file R1 with fastq gzip file (< {}). Files between {}-{} will be downsized randomly to ~{} before analysis.".format(
					humanfriendly.format_size(int(settings.MAX_FASTQ_FILE_WITH_DOWNSIZE)),
					humanfriendly.format_size(int(settings.MAX_FASTQ_FILE_UPLOAD)),
					humanfriendly.format_size(int(settings.MAX_FASTQ_FILE_WITH_DOWNSIZE)),
					humanfriendly.format_size(int(settings.MAX_FASTQ_FILE_UPLOAD))
					)
			message_r2 = "Max raw file R2 with fastq gzip file (< {}). Files between {}-{} will be downsized randomly to ~{} before analysis.".format(
					humanfriendly.format_size(int(settings.MAX_FASTQ_FILE_WITH_DOWNSIZE)),
					humanfriendly.format_size(int(settings.MAX_FASTQ_FILE_UPLOAD)),
					humanfriendly.format_size(int(settings.MAX_FASTQ_FILE_WITH_DOWNSIZE)),
					humanfriendly.format_size(int(settings.MAX_FASTQ_FILE_UPLOAD))
					)
		else:
			message_r1 = "Max raw file R1 with fastq gzip file (< {}).".format(
					humanfriendly.format_size(int(settings.MAX_FASTQ_FILE_UPLOAD))
					)
			message_r2 = "Max raw file R2 with fastq gzip file (< {}).".format(
					humanfriendly.format_size(int(settings.MAX_FASTQ_FILE_UPLOAD))
					)
		field_text= [
			# (field_name, Field title label, Detailed field description, requiered)
			('name', 'Name', 'Unique identifier for this sample. Only letters, numbers and underscores are allowed.', True),
			('date_of_onset', 'Onset date', 'Date of onset', False),
			('date_of_collection', 'Collection date', 'Date of collection', False),
			('date_of_receipt_lab', 'Lab reception date', 'Date receipt on the lab.', False),
			('vaccine_status', 'Vaccine status', 'Discrimination of vaccination status', False),
			('data_set', 'Dataset', 'Specific dataset (useful for grouping sets of samples)', False),
		##	('geo_local', 'Global position', 'Geo position where the sample was collected', False),
			('lat', 'Latitude', 'Geolocation where the sample was collected. Decimal Degrees (DD) format.', False),
			('lng', 'Longitude', 'Geolocation where the sample was collected. Decimal Degrees (DD) format.', False),
			('path_name_1', 'Raw fastq.gz (R1)', message_r1, True),
			('path_name_2', 'Raw fastq.gz (R2)', message_r2, False),
			('like_dates', 'Choose a date', 'Choose the option you want to be used to calculate the calendar week number.', False),
		]
		for x in field_text:
			self.fields[x[0]].label = x[1]
			self.fields[x[0]].help_text = x[2]
			self.fields[x[0]].required = x[3]

		### set formats of date
		self.fields['date_of_onset'].input_formats = settings.DATETIME_INPUT_FORMATS
		self.fields['date_of_collection'].input_formats = settings.DATETIME_INPUT_FORMATS
		self.fields['date_of_receipt_lab'].input_formats = settings.DATETIME_INPUT_FORMATS
		
		self.helper = FormHelper()
		self.helper.form_method = 'POST'
		self.helper.layout = Layout(
			Fieldset(
				'General data',
				Div(
					Div('name', css_class="col-sm-4"),
					Div('data_set', css_class="col-sm-4"), 
					Div('vaccine_status', css_class="col-sm-4"), 
					css_class = 'row'
				),
				css_class = 'article-content'
			),
			Fieldset(
				'Dates',
				Div(
					Div('like_dates', css_class="col-sm-3"),
					Div('date_of_onset', css_class="col-sm-3"),
					Div('date_of_collection', css_class="col-sm-3"),
					Div('date_of_receipt_lab', css_class="col-sm-3"),
				#	HTML('<div id="div_id_date_of_onset" class="form-group col-sm-3"> <label for="id_date_of_onset" class="form-control-label ">Date of onset</label>  <input placeholder="dd/mm/yyyy" name="date_of_onset" id="id_date_of_onset"/> <small id="hint_id_date_of_onset" class="text-muted">Date of onset</small> </div>'),
				#	HTML('<div id="div_id_date_of_collection" class="form-group col-sm-3"> <label for="id_date_of_collection" class="form-control-label ">Date of collection</label> <input placeholder="dd/mm/yyyy" name="date_of_collection" id="id_date_of_collection"/> <small id="hint_id_date_of_collection" class="text-muted">Date of collection</small> </div>'),
				#	HTML('<div id="div_id_date_on_lab" class="form-group col-sm-3"> <label for="id_date_on_lab" class="form-control-label ">Date on Lab</label> <input placeholder="dd/mm/yyyy" name="date_on_lab" id="id_date_on_lab"/> <small id="hint_id_date_on_lab" class="text-muted">Date on Lab</small> </div>'),
					css_class = 'row '),
				css_class = 'article-content'
			),
			Fieldset(
				'Global position',
				Div(
					Div('lat', css_class="col-sm-4"),
					Div('lng', css_class="col-sm-4"),
					css_class = 'row'
				),
				css_class = 'article-content'
			),
			Fieldset(
				'Fastq files',
				Div(
					Div('path_name_1', css_class = 'show-for-sr'),
					Div('path_name_2', css_class = 'show-for-sr')),
				css_class = 'article-content'
			),
# 			Div(
# 				Div('geo_local', css_class="col-sm-8"),
# 				css_class = 'row'
# 			),
			ButtonHolder(
				Submit('save', 'Save', css_class='btn-primary', onclick="this.disabled=true,this.form.submit();"),
				Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(reverse('samples'))),
			),
		)

	def clean(self):
		"""
		Clean all together because it's necessary to compare the genbank and fasta files
		"""
		cleaned_data = super(SampleForm, self).clean()
		name = self.cleaned_data.get('name', '').strip()
		if (len(name) == 0):
			self.add_error('name', _("Error: You must give a unique name."))
			return cleaned_data

		try:
			result_filer_name = re.sub('[^A-Za-z0-9_]+', '', name)
			if (len(result_filer_name) != len(name)):
				self.add_error('name', _("Error: Only letters, numbers and underscores are allowed."))
				return cleaned_data
				
			try:
				Sample.objects.get(name__iexact=name, owner__username=self.request.user.username, is_obsolete=False, is_deleted=False)
				self.add_error('name', ValidationError(_("This name '" + name +"' already exist in database, please choose other."), code='invalid'))
			except Sample.DoesNotExist:
				pass
			
			## latitude in degrees is -90 and +90 for the southern and northern hemisphere respectively. Longitude is in the range -180 and +180
			lat = self.cleaned_data.get('lat')
			if (lat != None and (lat > 90 or lat < -90)): self.add_error('lat', _("Latitute must have values between +90 and -90."))
			lng = self.cleaned_data.get('lng')
			if (lng != None and (lng > 180 or lng < -180)): self.add_error('lng', _("Longitude must have values between +180 and -180."))
			
			### test file name
			if ('path_name_1' not in self.cleaned_data):
				self.add_error('path_name_1', _("Error: must have a file."))
				return cleaned_data
			
			### testing file names
			path_name_1 = self.cleaned_data.get('path_name_1')
			path_name_2 = self.cleaned_data.get('path_name_2')
			
			## verbose log...
			self.logger_production.info('New Sample: {}  Path name1: {}'.format(name, path_name_1))
			self.logger_production.info('New Sample: {}  Path name2: {}'.format(name, path_name_2))
			self.logger_debug.info('New Sample: {}  Path name1: {}'.format(name, path_name_1))
			self.logger_debug.info('New Sample: {}  Path name2: {}'.format(name, path_name_2))

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
			
			## verbose log...
			self.logger_production.info('New Sample: {} Pass dates...'.format(name))
			self.logger_debug.info('New Sample: {} Pass dates...'.format(name))
			
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
				self.add_error('path_name_1', _(e.args[0]))
				return cleaned_data
			
			## verbose log...
			self.logger_production.info('New Sample: {} Pass is_fastq_gz 1...'.format(name))
			self.logger_debug.info('New Sample: {} Pass is_fastq_gz 1...'.format(name))
			
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
					self.add_error('path_name_2', _(e.args[0]))
					return cleaned_data
			
			## verbose log...
			self.logger_production.info('New Sample: {} end...'.format(name))
			self.logger_debug.info('New Sample: {} end...'.format(name))
			
			## remove temp files
			os.unlink(fastaq_temp_file_name.name)
			if (path_name_2 != None): os.unlink(fastaq_temp_file_name_2.name)
		except:
			self.logger_production.error("New Sample: {}  Path name2: {} Can't reach second file".format(name, path_name_1))
			self.logger_debug.error("New Sample: {}  Path name2: {} Can't reach second file".format(name, path_name_2))
			self.add_error('path_name_2', "Can't reach second file")
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
			('name', 'Name', 'Unique identifier for this project', True),
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
		return cleaned_data


## to put both inline
ReferenceProjectFormSet = inlineformset_factory(Reference, Project, form=ReferenceProjectForm, extra=1)

	
class AddSampleProjectForm(forms.ModelForm):
	"""
	Add samples to project 
	"""
	error_css_class = 'error'
	
	class Meta:
		model = Sample
		exclude = ()

	def __init__(self, *args, **kwargs):
		super(AddSampleProjectForm, self).__init__(*args, **kwargs)
		
	def clean(self):
		"""
		Clean all together because it's necessary to compare the genbank and fasta files
		"""
		cleaned_data = super(AddSampleProjectForm, self).clean()
		return cleaned_data
	


class SamplesUploadDescriptionForm(forms.ModelForm):
	"""
	Add samples to project without files 
	"""
	error_css_class = 'error'
	
	class Meta:
		model = UploadFiles
		exclude = ()
		fields = ('path_name', )

	def __init__(self, *args, **kwargs):
		self.request = kwargs.pop('request')
		super(SamplesUploadDescriptionForm, self).__init__(*args, **kwargs)
		
		field_text= [
			('path_name', 'File name', '"csv" or "tsv" file with several samples to upload fastq.gz files later', True),
		]
		for x in field_text:
			self.fields[x[0]].label = x[1]
			self.fields[x[0]].help_text = x[2]
			self.fields[x[0]].required = x[3]
		
		self.helper = FormHelper()
		self.helper.form_method = 'POST'
		self.helper.layout = Layout(
			HTML('<p> </p>'),
			HTML('<div class="alert alert-dark"> <a href="' + mark_safe(os.path.join(getattr(settings, "STATIC_URL", None), Constants.DIR_TEMPLATE_INPUT,\
					Constants.FILE_TEMPLATE_INPUT_csv)) + '" download> <span> <i class="fa fa-download"></i></span> Template file \'csv\'</a> </div>'),
			HTML('<p> </p>'),
			HTML('<div class="alert alert-dark"> <a href="' + mark_safe(os.path.join(getattr(settings, "STATIC_URL", None), Constants.DIR_TEMPLATE_INPUT,\
					Constants.FILE_TEMPLATE_INPUT_tsv)) + '" download> <span> <i class="fa fa-download"></i></span> Template file \'tsv\'</a> </div>'),
			HTML('<p> </p>'),
			HTML('<div class="alert alert-dark"> <a href="' + mark_safe(os.path.join(getattr(settings, "STATIC_URL", None), Constants.DIR_TEMPLATE_INPUT,\
					Constants.FILE_TEMPLATE_INPUT_data_csv)) + '" download> <span> <i class="fa fa-download"></i></span> Example Template file \'csv\'</a> </div>'),
			HTML('<p> </p>'),
			HTML('<div class="alert alert-dark"> <a href="' + mark_safe(os.path.join(getattr(settings, "STATIC_URL", None), Constants.DIR_TEMPLATE_INPUT,\
					Constants.FILE_TEMPLATE_INPUT_data_tsv)) + '" download> <span> <i class="fa fa-download"></i></span> Example Template file \'tsv\'</a> </div>'),
			HTML('<p> </p>'),
			Div('path_name', css_class="col-lm-3"),
			HTML('<p> </p>'),
			ButtonHolder(
				Submit('save', 'Upload', css_class='btn-primary', onclick="this.disabled=true,this.form.submit();"),
				Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(reverse('sample-add-file')))
			),
		)

	def clean(self):
		"""
		Clean all 
		"""
		cleaned_data = super(SamplesUploadDescriptionForm, self).clean()
		
		### get path name
		path_name = self.cleaned_data['path_name']
		
		## testing fastq
		temp_file_name = NamedTemporaryFile(prefix='flu_fq_', suffix='.csv', delete=False)
		temp_file_name.write(path_name.file.read())
		temp_file_name.flush()
		temp_file_name.close()
		
		parse_in_files = ParseInFiles()
		b_test_char_encoding = True
		parse_in_files.parse_sample_files(temp_file_name.name, self.request.user, b_test_char_encoding, ParseInFiles.STATE_READ_only_detect_errors)
		
		os.unlink(temp_file_name.name)
		if (parse_in_files.get_errors().has_errors()):
			self.add_error('path_name', _('There are errors in the file'))
			self.error_in_file = str(parse_in_files.get_errors())		## pass error_file kwargs for context
		else:
			self.number_files_to_process = parse_in_files.get_number_samples()
		return cleaned_data


class SamplesUploadDescriptionMetadataForm(forms.ModelForm):
	"""
	Upload csv file to update metadata 
	"""
	error_css_class = 'error'
	
	class Meta:
		model = UploadFiles
		exclude = ()
		fields = ('path_name', )

	def __init__(self, *args, **kwargs):
		self.request = kwargs.pop('request')
		super(SamplesUploadDescriptionMetadataForm, self).__init__(*args, **kwargs)
		
		field_text= [
			('path_name', 'File name', '"csv" or "tsv" file with metadata to update existing samples.', True),
		]
		for x in field_text:
			self.fields[x[0]].label = x[1]
			self.fields[x[0]].help_text = x[2]
			self.fields[x[0]].required = x[3]
		
		self.helper = FormHelper()
		self.helper.form_method = 'POST'
		self.helper.layout = Layout(
			HTML('<p> </p>'),
			HTML('<div class="alert alert-dark"> <a href="' + mark_safe(os.path.join(getattr(settings, "STATIC_URL", None), Constants.DIR_TEMPLATE_INPUT,\
					Constants.FILE_TEMPLATE_INPUT_METADATA_csv)) + '" download> <span> <i class="fa fa-download"></i></span> Template Metadata file \'csv\'</a> </div>'),
			HTML('<p> </p>'),
			HTML('<div class="alert alert-dark"> <a href="' + mark_safe(os.path.join(getattr(settings, "STATIC_URL", None), Constants.DIR_TEMPLATE_INPUT,\
					Constants.FILE_TEMPLATE_INPUT_METADATA_tsv)) + '" download> <span> <i class="fa fa-download"></i></span> Template Metadata file \'tsv\'</a> </div>'),
			HTML('<p> </p>'),
			HTML('<div class="alert alert-dark"> <a href="' + mark_safe(os.path.join(getattr(settings, "STATIC_URL", None), Constants.DIR_TEMPLATE_INPUT,\
					Constants.FILE_TEMPLATE_INPUT_METADATA_data_csv)) + '" download> <span> <i class="fa fa-download"></i></span> Example Template Metadata file \'csv\'</a> </div>'),
			HTML('<p> </p>'),
			HTML('<div class="alert alert-dark"> <a href="' + mark_safe(os.path.join(getattr(settings, "STATIC_URL", None), Constants.DIR_TEMPLATE_INPUT,\
					Constants.FILE_TEMPLATE_INPUT_METADATA_data_tsv)) + '" download> <span> <i class="fa fa-download"></i></span> Example Template Metadata file \'tsv\'</a> </div>'),
			HTML('<p> </p>'),
			Div('path_name', css_class="col-lm-3"),
			HTML('<p> </p>'),
			ButtonHolder(
				Submit('save', 'Upload', css_class='btn-primary', onclick="this.disabled=true,this.form.submit();"),
				Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(reverse('sample-update-metadata')))
			),
		)

	def clean(self):
		"""
		Clean all 
		"""
		cleaned_data = super(SamplesUploadDescriptionMetadataForm, self).clean()
		
		### get path name
		path_name = self.cleaned_data['path_name']
		
		## testing fastq
		temp_file_name = NamedTemporaryFile(prefix='flu_fq_', suffix='.csv', delete=False)
		temp_file_name.write(path_name.file.read())
		temp_file_name.flush()
		temp_file_name.close()
		
		parse_in_files = ParseInFiles()
		b_test_char_encoding = True
		parse_in_files.parse_sample_files(temp_file_name.name, self.request.user, b_test_char_encoding,\
									ParseInFiles.STATE_READ_metadata_only_detect_errors_and_chech_samples)
		
		os.unlink(temp_file_name.name)
		if (parse_in_files.get_errors().has_errors()):
			self.add_error('path_name', _('There are errors in the file'))
			self.error_in_file = str(parse_in_files.get_errors())		## pass error_file kwargs for context
		else:
			self.number_files_to_process = parse_in_files.get_number_samples()
		return cleaned_data


class SamplesUploadMultipleFastqForm(forms.ModelForm):
	"""
	Add multiple fastq files 
	"""
	
	class Meta:
		model = UploadFiles
		fields = ('path_name', )

	def __init__(self, *args, **kwargs):
		self.request = kwargs.pop('request')
		super(SamplesUploadMultipleFastqForm, self).__init__(*args, **kwargs)
		
	def clean(self):
		"""
		Clean all 
		"""
		cleaned_data = super(SamplesUploadMultipleFastqForm, self).clean()
		
		### get path name
		if ('path_name' not in self.cleaned_data):
			self.add_error('path_name', "There's no file to upload")
			return cleaned_data
		path_name = self.cleaned_data['path_name']
		
		### test the file name if exist
		number_files = UploadFiles.objects.filter(file_name__iexact=os.path.basename(path_name.name), owner=self.request.user,\
 										is_processed=False, is_deleted=False, type_file__name=TypeFile.TYPE_FILE_fastq_gz).count()

		if (number_files > 0):
			self.add_error('path_name', "There's one file with this name already")
			return cleaned_data
		
		utils = Utils()
		## testing fastq
		temp_file_name = NamedTemporaryFile(prefix='flu_fq_', suffix='.fastq.gz', delete=False)
		if (str(type(path_name.file)) == "<class '_io.BytesIO'>"):
			with open(temp_file_name.name, 'wb') as out: ## Open temporary file as bytes
				out.write(path_name.file.read())                ## Read bytes into file
		else: utils.link_file(path_name.file.name, temp_file_name.name)

		try:
			utils.is_fastq_gz(temp_file_name.name)
		except Exception as e:	## (e.errno, e.strerror)
			os.unlink(temp_file_name.name)
			self.add_error('path_name', e.args[0])
			return cleaned_data

		if (settings.DOWN_SIZE_FASTQ_FILES):
			if (path_name.size > settings.MAX_FASTQ_FILE_WITH_DOWNSIZE):
				os.unlink(temp_file_name.name)
				self.add_error('path_name', "Max file size is: {}".format( humanfriendly.format_size(int(settings.MAX_FASTQ_FILE_WITH_DOWNSIZE)) ))
				return cleaned_data
		elif (path_name.size > settings.MAX_FASTQ_FILE_UPLOAD):
			os.unlink(temp_file_name.name)
			self.add_error('path_name', "Max file size is: {}".format( humanfriendly.format_size(int(settings.MAX_FASTQ_FILE_UPLOAD)) ))
			return cleaned_data
		
		os.unlink(temp_file_name.name)
		return cleaned_data

