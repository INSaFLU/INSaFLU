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
		empty_text = "There aren't no References to show..."
		
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
	tag_names = tables.Column('Tags')
	quality = tables.Column('Quality Available')
	
	class Meta:
		model = Reference
		fields = ('name', 'sample_date', 'number_quality_sequences', 'tag_names', 'quality')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There aren't no Samples to show..."
		
	def render_number_quality_sequences(self, **kwargs):
		record = kwargs.pop("record")
		
		href = os.path.join(getattr(settings, "MEDIA_URL", None), record.reference_genbank.name)		
		return mark_safe('<a href="' + href + '">' + record.reference_genbank_name + '</a>')

	def render_quality(self, **kwargs):
		record = kwargs.pop("record")
		return 'Not yet'

	def render_tag_names(self, **kwargs):
		record = kwargs.pop("record")
		text_to_draw = record.tag_names.join(', ')
		return text_to_draw
