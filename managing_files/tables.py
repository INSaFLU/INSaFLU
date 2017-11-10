import django_tables2 as tables
from django_tables2.utils import A
from .models import Reference
from django.utils.safestring import mark_safe
from django.conf import settings
import os

class ReferenceTable(tables.Table):
#   Renders a normal value as an internal hyperlink to another page.
#   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
	reference_fasta_name = tables.LinkColumn('reference_fasta_name', args=[tables.A('pk')], verbose_name='Fasta file')
	reference_genbank_name = tables.LinkColumn('reference_genbank_name', args=[tables.A('pk')], verbose_name='GenBank file')
	
	class Meta:
		model = Reference
		fields = ('name', 'scentific_name', 'reference_fasta_name', 'reference_genbank_name',
				  'creation_date', 'is_obsolete', 'number_of_locus')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There are no References to show..."
		
	def render_reference_fasta_name(self, **kwargs):
		record = kwargs.pop("record")
		href = os.path.join(getattr(settings, "MEDIA_URL", None), record.reference_fasta.name)		
		return mark_safe('<a href="' + href + '">' + record.reference_fasta_name + '</a>')
	
	def render_reference_genbank_name(self, **kwargs):
		record = kwargs.pop("record")
		href = os.path.join(getattr(settings, "MEDIA_URL", None), record.reference_genbank.name)		
		return mark_safe('<a href="' + href + '">' + record.reference_genbank_name + '</a>')
	

class SampleTable(tables.Table):
#   Renders a normal value as an internal hyperlink to another page.
#   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
	number_quality_sequences = tables.Column('#Quality Seq.')
	extra_info = tables.Column('Extra Information')
	type_and_subtype = tables.Column('Type and SubType')
	
	class Meta:
		model = Reference
		fields = ('name', 'creation_date', 'fastq_files', 'type_and_subtype', 'number_quality_sequences', 'extra_info')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There are no Samples to show..."
	
	def render_fastq_files(self, **kwargs):
		"""
		number of fastqFiles, to show if there is 
		"""
		record = kwargs.pop("record")
		if (record.has_files()):
			if (record.file_name_2 is None or len(record.file_name_2) == 0): return "1"
			return "2"
		return "0"
	
	def render_type_and_subtype(self, **kwargs):
		"""
		get type and sub type
		"""
		record = kwargs.pop("record")
		return 'Not yet'
	
	def render_number_quality_sequences(self, **kwargs):
		"""
		number of quality sequences and average
		"""
		record = kwargs.pop("record")
		
		href = os.path.join(getattr(settings, "MEDIA_URL", None), record.reference_genbank.name)		
		return mark_safe('<a href="' + href + '">' + record.reference_genbank_name + '</a>')

	def render_extra_info(self, **kwargs):
		"""
		icon with link to extra info
		"""
		record = kwargs.pop("record")
		return 'Not yet'
	
	

	
