import django_tables2 as tables
from django_tables2.utils import A
from managing_files.models import Reference, Sample, Project, ProjectSample
from django.utils.safestring import mark_safe
from managing_files.manage_database import ManageDatabase
from constants.meta_key_and_values import MetaKeyAndValue
from utils.result import DecodeResultAverageAndNumberReads
from constants.constants import Constants, TypePath
from django.utils.translation import ugettext_lazy as _
from django.urls import reverse
from utils.result import DecodeCoverage
from django.conf import settings

class CheckBoxColumnWithName(tables.CheckBoxColumn):
	@property
	def header(self):
		default = {'type': 'checkbox'}
		return self.verbose_name


class ReferenceTable(tables.Table):
#   Renders a normal value as an internal hyperlink to another page.
#   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
	reference_fasta_name = tables.LinkColumn('reference_fasta_name', args=[tables.A('pk')], verbose_name='Fasta file')
	reference_genbank_name = tables.LinkColumn('reference_genbank_name', args=[tables.A('pk')], verbose_name='GenBank file')
	owner = tables.Column("Owner", orderable=True, empty_values=())
	number_of_locus = tables.Column(orderable=False)

	class Meta:
		model = Reference
		fields = ('name', 'isolate_name', 'reference_fasta_name', 'reference_genbank_name',
				  'creation_date', 'owner', 'number_of_locus')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There are no References to show..."

	def render_owner(self, **kwargs):
		record = kwargs.pop("record")
		return record.owner.username
	
	def render_reference_fasta_name(self, **kwargs):
		record = kwargs.pop("record")
		href = record.get_reference_fasta(TypePath.MEDIA_URL)		
		return mark_safe('<a href="' + href + '" download>' + record.reference_fasta_name + '</a>')

	def render_reference_genbank_name(self, **kwargs):
		record = kwargs.pop("record")
		href = record.get_reference_gbk(TypePath.MEDIA_URL)		
		return mark_safe('<a href="' + href + '" download>' + record.reference_genbank_name + '</a>')

	def render_creation_date(self, **kwargs):
		record = kwargs.pop("record")
		return record.creation_date.strftime(settings.DATE_FORMAT_FOR_TABLE)

class ReferenceProjectTable(tables.Table):
	"""
	create a new project with a reference
	"""
	select_ref = CheckBoxColumnWithName(verbose_name=('Select One'), accessor="pk", orderable=False)
#	number_of_locus = tables.Column(orderable=False)

	class Meta:
		model = Reference
		fields = ('select_ref', 'name', 'isolate_name', 'creation_date', 'owner', 'number_of_locus')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There are no References to add..."

	def render_creation_date(self, **kwargs):
		record = kwargs.pop("record")
		return record.creation_date.strftime(settings.DATE_FORMAT_FOR_TABLE)

	def render_select_ref(self, value, record):
		return mark_safe('<input name="select_ref" id="{}_{}" type="checkbox" value="{}"/>'.format(Constants.CHECK_BOX, record.id, value))

class SampleToProjectsTable(tables.Table):
	"""
	To add samples to projects
	"""
	type_subtype = tables.Column("Type-Subtype", orderable=True, empty_values=())
	select_ref = tables.CheckBoxColumn(accessor="pk", attrs = { "th__input": {"id": "checkBoxAll"}}, orderable=False)
	
	class Meta:
		model = Sample
		fields = ('select_ref', 'name', 'creation_date', 'type_subtype', 'week', 'data_set')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There are no Samples to add..."
	
	def render_select_ref(self, value, record):
		return mark_safe('<input name="select_ref" id="{}_{}" type="checkbox" value="{}"/>'.format(Constants.CHECK_BOX, record.id, value))

	def render_creation_date(self, **kwargs):
		record = kwargs.pop("record")
		return record.creation_date.strftime(settings.DATE_FORMAT_FOR_TABLE)

class SampleTable(tables.Table):
#   Renders a normal value as an internal hyperlink to another page.
#   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
	number_quality_sequences = tables.Column('#Quality Seq.', orderable=False, empty_values=())
#	extra_info = tables.LinkColumn('sample-description', args=[tables.A('pk')], orderable=False, verbose_name='Extra Information', empty_values=())
	extra_info = tables.LinkColumn('Extra Information', orderable=False, empty_values=())
	type_and_subtype = tables.Column('Type and SubType', empty_values=())
	fastq_files = tables.Column('#Fastq Files', empty_values=())
	data_set = tables.Column('Data set', empty_values=())
	
	class Meta:
		model = Sample
		fields = ('name', 'creation_date', 'fastq_files', 'type_and_subtype', 'data_set', 'number_quality_sequences', 'extra_info')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There are no Samples to show..."
	
	def render_fastq_files(self, record):
		"""
		number of fastqFiles, to show if there is 
		"""
		if (record.has_files):
			if (record.file_name_2 is None or len(record.file_name_2) == 0): return "1"
			return "2"
		return "0"
	
	def render_creation_date(self, **kwargs):
		record = kwargs.pop("record")
		return record.creation_date.strftime(settings.DATE_FORMAT_FOR_TABLE)

	def render_type_and_subtype(self, record):
		"""
		get type and sub type
		"""
		result = record.get_type_sub_type()
		return _('Not yet') if result == Constants.EMPTY_VALUE_TABLE else record.get_type_sub_type()
	
	def render_data_set(self, record):
		"""
		return name of dataset
		"""
		return record.name
		
		
	def render_number_quality_sequences(self, record):
		"""
		number of quality sequences and average
		"""
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(record, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None)
		if (list_meta.count() > 0 and list_meta[0].value == MetaKeyAndValue.META_VALUE_Success):
			decodeResultAverageAndNumberReads = DecodeResultAverageAndNumberReads()
			result_average = decodeResultAverageAndNumberReads.decode_result(list_meta[0].description)
			if (result_average.number_file_2 is None): return _('%s/%s' % (result_average.number_file_1, result_average.average_file_1 ))
			return _('%s/%s-%s/%s' % (result_average.number_file_1,\
					result_average.average_file_1, result_average.number_file_2,\
					result_average.average_file_2) )
		elif (list_meta.count() > 0 and list_meta[0].value.equals(MetaKeyAndValue.META_VALUE_Error)):
			return _("Error")
		return _('Not yet')

	def render_extra_info(self, record):
		"""
		icon with link to extra info
		"""
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(record, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
		if (list_meta.count() > 0 and list_meta[0].value == MetaKeyAndValue.META_VALUE_Success):
			return mark_safe('<a href=' + reverse('sample-description', args=[record.pk]) + '><span ><i class="fa fa-plus-square"></i></span> More Info</a>')
		elif (list_meta.count() > 0 and list_meta[0].value == MetaKeyAndValue.META_VALUE_Error): return _("Error")
		return _('Not yet')
	
	
class ProjectTable(tables.Table):
#   Renders a normal value as an internal hyperlink to another page.
#   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
	reference = tables.Column('Reference', empty_values=())
	samples = tables.Column('#Samples (P/W/E)', orderable=False, empty_values=())
	creation_date = tables.Column('Creation date', empty_values=())
	results = tables.LinkColumn('Results', orderable=False, empty_values=())
	
	class Meta:
		model = Project
		fields = ('name', 'reference', 'creation_date', 'samples', 'results')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There are no Projects to show..."
	
	def render_samples(self, record):
		"""
		return a reference name
		"""
		add_remove = ""
		if (ProjectSample.objects.filter(project__id=record.id, is_deleted=False).count() > 0):
			add_remove = ' <a href=' + reverse('remove-sample-project', args=[record.pk]) + '><span ><i class="fa fa-trash"></i></span> Remove</a>'
		return mark_safe("({}/{}/{}) ".format("0", "0", "0") + '<a href=' + reverse('add-sample-project', args=[record.pk]) +\
						 '><i class="fa fa-plus-square"></i> Add</a>' + add_remove)
	
	def render_creation_date(self, **kwargs):
		record = kwargs.pop("record")
		return record.creation_date.strftime(settings.DATE_FORMAT_FOR_TABLE)
	
	def render_results(self, record):
		"""
		icon with link to extra info
		"""
		## there's nothing to show
		count = ProjectSample.objects.filter(project__id=record.id, is_deleted=False, is_error=False, is_finished=True).count()
		if (count > 0): return mark_safe('<a href=' + reverse('show-sample-project-results', args=[record.pk]) + '><span ><i class="fa fa-info-circle"></i></span> More info</a>')
		count_not_finished = ProjectSample.objects.filter(project__id=record.id, is_deleted=False, is_error=False, is_finished=False).count()
		if (count_not_finished > 0): return _("{} processing".format(count_not_finished))
		return "-"
	

class ShowProjectSamplesResults(tables.Table):
	"""
	Has the results from the result samples processed against references
	"""
	sample_name = tables.Column('Sample name', empty_values=())
	coverage = tables.Column('Coverage >9', orderable=False, empty_values=())
	alerts = tables.Column('Alerts', empty_values=())
	results = tables.LinkColumn('Results', orderable=False, empty_values=())
	type_subtype = tables.LinkColumn('Type-Subtype', empty_values=())
	dataset = tables.LinkColumn('Dataset', empty_values=())
	results = tables.LinkColumn('Results', orderable=False, empty_values=())
		
	class Meta:
		model = ProjectSample
		fields = ('sample_name', 'type_subtype', 'dataset', 'coverage', 'alerts', 'results')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There are no samples processed to show..."
	
	def render_sample_name(self, record):
		"""
		return sample name
		"""
		return record.sample.name
	
	def render_coverage(self, record):
		"""
		return icons about coverage
		"""
		manageDatabase = ManageDatabase()
		meta_value = manageDatabase.get_project_sample_metakey(record, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
		decode_coverage = DecodeCoverage()
		coverage = decode_coverage.decode_result(meta_value.description)
		return_html = ""
		for key in coverage.get_sorted_elements_name():
			return_html += '<a href="#coverageModal" id="showImageCoverage" data-toggle="modal" project_sample_id="{}" sequence="{}"><img title="{}" class="tip" src="{}"></a>'.format(\
					record.id, key, coverage.get_message_to_show_in_web_site(key), coverage.get_icon(key))
		return mark_safe(return_html)
		
	def render_alerts(self, record):
		"""
		return number
		"""
		return "{}".format(record.alert_first_level + record.alert_second_level)
	
	def render_type_subtype(self, record):
		"""
		return number
		"""
		return record.sample.type_subtype
	
	def render_dataset(self, record):
		"""
		return number
		"""
		return record.sample.data_set.name
	
	def render_results(self, record):
		"""
		icon with link to extra info
		"""
		return mark_safe('<a href=' + reverse('show-sample-project-results', args=[record.sample.pk]) + '><span ><i class="fa fa-info-circle"></i> More info</a>')

