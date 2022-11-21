'''
Created on 22/05/2022

@author: mmp
'''
from constants.software_names import SoftwareNames
import django_tables2 as tables
from django.urls import reverse
from django.conf import settings
from django.db.models import F
from datasets.models import Dataset, DatasetConsensus
from django.utils.safestring import mark_safe
from constants.constants import Constants
from settings.constants_settings import ConstantsSettings
from managing_files.models import Reference, Project, UploadFiles
from settings.default_software_project_sample import DefaultProjectSoftware
from settings.models import Parameter


class DatasetTable(tables.Table):
#   Renders a normal value as an internal hyperlink to another page.
#   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
	### P -> sequences from projects
	### C -> sequences from consensus
	### R -> sequences from references

	# build cannot be orderable because it is not (yet?) an explicit field of the Dataset model
	build = tables.Column('NextStrain Build', orderable=False, empty_values=())
	sequences = tables.Column('#Sequences (P/C/R)', orderable=False, empty_values=())
	last_change_date = tables.Column('Last Change date', empty_values=())
	creation_date = tables.Column('Creation date', empty_values=())
	results = tables.LinkColumn('Options', orderable=False, empty_values=())
	
	class Meta:
		model = Dataset
		fields = ('name',  'build', 'last_change_date','creation_date', 'totla_alerts', 'sequences', 'results')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There are no Datasets to show..."
	
	def render_name(self, record):

		from crequest.middleware import CrequestMiddleware
		current_request = CrequestMiddleware.get_request()
		user = current_request.user

		## there's nothing to show
		count = record.number_of_sequences_from_projects + record.number_of_sequences_from_consensus + record.number_of_sequences_from_references
		if (count > 0):
			project_sample = '<a href=' + reverse('show-dataset-consensus', args=[record.pk]) + ' data-toggle="tooltip" title="Show sequences">' +\
				'{}</a>'.format(record.name)
		else:
			project_sample = record.name
		
		if (user.username == Constants.USER_ANONYMOUS): return mark_safe(project_sample);
		if (user.username == record.owner.username):
			return mark_safe('<a href="#id_remove_modal" id="id_remove_dataset_modal" data-toggle="modal" data-toggle="tooltip" title="Delete"' +\
					' ref_name="' + record.name + '" pk="' + str(record.pk) + '"><i class="fa fa-trash"></i></span> </a>' + project_sample)
		return project_sample


	def render_build(self, **kwargs):
		
		record = kwargs.pop("record")
		
		build = 'NA'
		parameters = Parameter.objects.filter(dataset__pk=record.pk)
		if(len(parameters) > 0):
			build = parameters[0].parameter
		
		for [id, desc] in SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_DESC:
			if id == build:
				build = desc
				break

		return build


	#def order_build(self, queryset, is_descending):
	#	queryset = queryset.annotate(
	#		amount=F("shirts") + F("pants")
	#	).order_by(("-" if is_descending else "") + "amount")
	#	return (queryset, True)


	def render_sequences(self, record):
		"""
		return a reference name
		"""
		tip_info = '<span ><i class="tip fa fa-info-circle" title="Consensus from projects (P): {}\nConsensus uploaded (C): {}\nReferences (R): {}"></i></span>'.format(\
			record.number_of_sequences_from_projects, record.number_of_sequences_from_consensus, record.number_of_sequences_from_references)
		add_sequence = '<div class="btn-group">' + \
			'<button type="button" class="btn btn-primary dropdown-toggle" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false"' +\
			'data-toggle="tooltip" title="Add Consensus/References/Projects" style="margin-left: 8px;">Add Sequences</button>' + \
			'<div class="dropdown-menu" x-placement="bottom-start" style="position: absolute; transform: translate3d(0px, 38px, 0px); top: 0px; left: 0px; will-change: transform;">' + \
			'<a rel="nofollow" class="dropdown-item" href=' + reverse('add-references-dataset', args=[record.pk]) + '> Add References</a>' + \
			'<a rel="nofollow" class="dropdown-item" href=' + reverse('add-projects-dataset', args=[record.pk]) + ' > Add Consensus from Projects</a>' + \
			'<a rel="nofollow" class="dropdown-item" href=' + reverse('add-consensus-dataset', args=[record.pk]) + ' > Add your own Consensus</a>' + \
			'</div> </div>'
			
		return mark_safe(tip_info + " ({}/{}/{})   ".format(record.number_of_sequences_from_projects, record.number_of_sequences_from_consensus, \
						record.number_of_sequences_from_references) + \
						add_sequence)
	
	def render_creation_date(self, **kwargs):
		record = kwargs.pop("record")
		return record.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)
	
	def render_last_change_date(self, **kwargs):
		record = kwargs.pop("record")
		if (record.last_change_date is None): return "Not set yet"
		return record.last_change_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)
	
	def render_results(self, record):
		"""
		icon with link to extra info
		"""

		from crequest.middleware import CrequestMiddleware
		current_request = CrequestMiddleware.get_request()
		user = current_request.user
		
		## there's nothing to show
		count = record.number_of_sequences_from_projects + record.number_of_sequences_from_consensus + record.number_of_sequences_from_references
		sz_project_sample = ""
		if (count > 0):
			sz_project_sample = '<a href=' + reverse('show-dataset-consensus', args=[record.pk]) + ' data-toggle="tooltip" title="See Results"> ' +\
				'<span ><i class="padding-button-table fa fa-info-circle padding-button-table"></i></span></a> '
			## Can also launch all the sequeneces in the NextTrain
			sz_project_sample += '<a href=' + reverse('dataset-settings', args=[record.pk]) + ' data-toggle="tooltip" title="Dataset settings">' +\
				'<span ><i class="fa fa-magic padding-button-table"></i></span></a>'
			# TODO Make this a query to ProcessSGE instead
			if(record.is_processed):
				if(user.username != Constants.USER_ANONYMOUS):
					# Do not change this data-toggle="modal"...
					sz_project_sample += '<a href="#id_rebuild_modal" id="id_rebuild_dataset_modal" data-toggle="modal" title="Rebuild Dataset Results"' +\
						' ref_name="' + record.name + '" pk="' + str(record.pk) + '"><span ><i class="fa fa-flask" style="color:#55aa55;"></i></span> </a>'
			else:
				sz_project_sample += '<a href="#id_rebuild_modal" id="id_rebuild_dataset_modal" data-toggle="modal" title="Waiting for Dataset Results (click to rebuild at your own risk)"' +\
					' ref_name="' + record.name + '" pk="' + str(record.pk) + '"><span ><i class="fa fa-spinner fa-spin" style="color:#ff0000;"></i></span> </a>'				
		else:	## can change settings
			sz_project_sample += '<a href=' + reverse('dataset-settings', args=[record.pk]) + ' data-toggle="tooltip" title="Dataset settings">' +\
				 '<span ><i class="padding-button-table fa fa-magic padding-button-table"></i></span></a>'
	
		return mark_safe(sz_project_sample)


class ReferenceTable(tables.Table):
#   Renders a normal value as an internal hyperlink to another page.
#   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
	select_ref = tables.CheckBoxColumn(accessor="pk", attrs = { "th__input": {"id": "checkBoxAll"}}, orderable=False)
	reference_fasta_name = tables.LinkColumn('reference_fasta_name', args=[tables.A('pk')], verbose_name='Fasta file')
	reference_genbank_name = tables.LinkColumn('reference_genbank_name', args=[tables.A('pk')], verbose_name='GenBank file')
	owner = tables.Column("Owner", orderable=True, empty_values=())
	constants = Constants()
	
	SHORT_NAME_LENGTH = 20
	
	class Meta:
		model = Reference
		fields = ('select_ref', 'name', 'reference_fasta_name', 'reference_genbank_name',
				  'creation_date', 'owner', 'number_of_locus')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There are no References to show..."

	def render_select_ref(self, value, record):
		return mark_safe('<input name="select_ref" id="{}_{}" type="checkbox" value="{}"/>'.format(Constants.CHECK_BOX, record.id, value))
	
	def render_owner(self, record):
		return record.owner.username
	
	def render_name(self, record):
		from crequest.middleware import CrequestMiddleware
		current_request = CrequestMiddleware.get_request()
		user = current_request.user
		if (user.username == Constants.USER_ANONYMOUS): return record.name
		return record.name
	
	def render_reference_fasta_name(self, **kwargs):
		record = kwargs.pop("record")
		return record.get_reference_fasta_web()

	def render_reference_genbank_name(self, **kwargs):
		record = kwargs.pop("record")
		return record.get_reference_gb_web()

	def render_creation_date(self, **kwargs):
		record = kwargs.pop("record")
		return record.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)

	def order_owner(self, queryset, is_descending):
		queryset = queryset.annotate(owner_name = F('owner__username')).order_by(('-' if is_descending else '') + 'owner_name')
		return (queryset, True)


class ConsensusTable(tables.Table):
#   Renders a normal value as an internal hyperlink to another page.
#   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
	select_ref = tables.CheckBoxColumn(accessor="pk", attrs = { "th__input": {"id": "checkBoxAll"}}, orderable=False)
	consensus_fasta_name = tables.LinkColumn('consensus_fasta_name', args=[tables.A('pk')], verbose_name='Fasta file')
	owner = tables.Column("Owner", orderable=True, empty_values=())
	constants = Constants()
	
	SHORT_NAME_LENGTH = 20
	
	class Meta:
		model = Reference
		fields = ('select_ref', 'name', 'consensus_fasta_name',
				  'creation_date', 'owner', 'number_of_locus')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There are no Consensus to show..."

	def render_select_ref(self, value, record):
		return mark_safe('<input name="select_ref" id="{}_{}" type="checkbox" value="{}"/>'.format(Constants.CHECK_BOX, record.id, value))
	
	def render_owner(self, record):
		return record.owner.username
	
	def render_name(self, record):
		from crequest.middleware import CrequestMiddleware
		current_request = CrequestMiddleware.get_request()
		user = current_request.user
		
		if (user.username == Constants.USER_ANONYMOUS): return record.owner.username
		if (user.username == record.owner.username and 
			DatasetConsensus.objects.filter(consensus__pk=record.pk, dataset__is_deleted=False, is_deleted=False).count() == 0):
			return mark_safe('<a href="#id_remove_modal" id="id_remove_consensus_modal" data-toggle="modal" data-toggle="tooltip" title="Delete"' +\
					' ref_name="' + record.name + '" pk="' + str(record.pk) + '"><i class="fa fa-trash"></i></span> </a>' + record.name)
		return record.name
	
	def render_consensus_fasta_name(self, **kwargs):
		record = kwargs.pop("record")
		return mark_safe(record.get_consensus_fasta_web())

	def render_creation_date(self, **kwargs):
		record = kwargs.pop("record")
		return record.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)

	def order_owner(self, queryset, is_descending):
		queryset = queryset.annotate(owner_name = F('owner__username')).order_by(('-' if is_descending else '') + 'owner_name')
		return (queryset, True)


class ProjectTable(tables.Table):
#   Renders a normal value as an internal hyperlink to another page.
#   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
	select_ref = tables.CheckBoxColumn(accessor="pk", attrs = { "th__input": {"id": "checkBoxAll"}}, orderable=False)
	reference = tables.Column('Reference', empty_values=())
	samples = tables.Column('#Consensus (S/A)', orderable=False, empty_values=())
	last_change_date = tables.Column('Last Change date', empty_values=())
	creation_date = tables.Column('Creation date', empty_values=())
	results = tables.LinkColumn('Options', orderable=False, empty_values=())
	
	class Meta:
		model = Project
		fields = ('select_ref', 'name', 'reference', 'last_change_date','creation_date', 'samples', 'results')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There are no Projects to show..."
	
	def render_name(self, record):
		return record.name
			
	def render_samples(self, record):
		"""
		return a reference name
		"""
		tip_info = '<span ><i class="tip fa fa-info-circle" title="Selected: {}\nAvailable: {}"></i></span>'.format(
			record.number_passed_sequences, record.number_passed_sequences)
		return mark_safe(tip_info + " ({}/{}) ".format(record.number_passed_sequences,
				record.number_passed_sequences))
	
	def render_creation_date(self, **kwargs):
		record = kwargs.pop("record")
		return record.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)
	
	def render_last_change_date(self, **kwargs):
		record = kwargs.pop("record")
		if (record.last_change_date is None): return "Not set yet"
		return record.last_change_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)
	
	def render_results(self, record):
		"""
		icon with link to extra info
		"""
		## there's nothing to show
		# count = ProjectSample.objects.filter(project__id=record.id, is_deleted=False, is_error=False, is_finished=True).count()
		
		sz_project_sample = "---"
		# if (count > 0):
		#	 sz_project_sample = '<a href=' + reverse('show-dataset-consensus', args=[record.pk]) + ' data-toggle="tooltip" title="See Results"> ' +\
		#		 '<span ><i class="padding-button-table fa fa-info-circle padding-button-table"></i></span></a> '
		#	 ## only can change settings when has projects finished
		#	 #sz_project_sample += '<a href=' + reverse('project-settings', args=[record.pk]) + ' data-toggle="tooltip" title="Software settings">' +\
		#	 #	'<span ><i class="fa fa-magic padding-button-table"></i></span></a>'
		#	 #sz_project_sample += '<a href=' + reverse('project-settings', args=[record.pk]) + ' data-toggle="tooltip" title="Software settings">' +\
		#	 #	'<span ><i class="padding-button-table fa fa-magic padding-button-table"></i></span></a>'
		# else:	## can change settings
		#	 #sz_project_sample += '<a href=' + reverse('project-settings', args=[record.pk]) + ' data-toggle="tooltip" title="Software settings">' +\
		#	 #	'<span ><i class="padding-button-table fa fa-magic padding-button-table"></i></span></a>'
		
		return mark_safe(sz_project_sample)

class DatasetConsensusTable(tables.Table):
	"""
	Has the results from the result samples processed against references
	"""
	
	name = tables.Column('Name', empty_values=())
	project_name = tables.Column('Project Name', orderable=False, empty_values=())		## when came from projects
	source = tables.Column('Source', orderable=False, empty_values=())
	type_and_subtype = tables.LinkColumn('Classification', orderable=False, empty_values=())
	alerts = tables.Column('Alerts', empty_values=())
	technology = tables.Column('Technology', orderable=False, empty_values=())
	consensus_file = tables.LinkColumn('Consensus File', orderable=False, empty_values=())
	results = tables.LinkColumn('Options', orderable=False, empty_values=())
	EMPTY = "---"
	
	class Meta:
		model = DatasetConsensus
		fields = ('name', 'project_name', 'source', 'type_and_subtype', 'technology',\
			'alerts', 'consensus_file', 'results')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There are no consensus to show..."
	
	def render_project_name(self, record):
		"""
		return project name
		"""
		return record.get_project_name() if len(record.get_project_name()) > 0 else self.EMPTY
	
	def render_name(self, record):
		"""
		return sample name
		"""
		from crequest.middleware import CrequestMiddleware
		current_request = CrequestMiddleware.get_request()
		
		user = current_request.user
		if (user.username == Constants.USER_ANONYMOUS): return record.name;
		if (user.username == record.dataset.owner.username):
			return mark_safe('<a href="#id_remove_modal" id="id_remove_consensus_modal" data-toggle="modal"' +\
					' ref_name="' + record.name + '" pk="' + str(record.pk) + '" +\
					" ref_project="' + record.name + '" data-toggle="tooltip" title="Remove consensus">' +\
					'<i class="fa fa-trash"></i></span></a> {}'.format(record.name))
					## '<a href=' + reverse('show-sample-project-single-detail', args=[record.pk]) + ' data-toggle="tooltip" title="Show more information">' +\
					## '{}</a>'.format(record.sample.name))
		return record.name

	def render_technology(self, record):
		""" shows if it is Illumina or Minion """
		if  not record.project_sample is None:
			return ConstantsSettings.TECHNOLOGY_illumina if record.project_sample.sample.is_type_fastq_gz_sequencing()\
				else ConstantsSettings.TECHNOLOGY_minion
		return self.EMPTY

	def render_alerts(self, record):
		"""
		return number
		"""
		return "{}".format(record.alert_first_level + record.alert_second_level)
	
	def render_source(self, record):
		"""
		return source of choise
		"""
		if  not record.project_sample is None: 
			return "Con. Project"
		elif not record.consensus is None:
			return "Private Con."
		elif not record.reference is None:
			return "Reference"
		return DatasetConsensusTable.EMPTY
	
	def render_type_and_subtype(self, record):
		"""
		return number
		"""
		if  not record.project_sample is None: 
			return record.project_sample.sample.type_subtype
		return self.EMPTY
	
		return record.sample.type_subtype
	
	def render_consensus_file(self, record):
		"""
		return link to consensus file
		"""
		default_software = DefaultProjectSoftware()
		if  not record.project_sample is None: 
			return record.project_sample.get_consensus_file_web(not default_software.include_consensus(record.project_sample))
		elif not record.consensus is None:
			return record.consensus.get_consensus_fasta_web()
		elif not record.reference is None:
			return record.reference.get_reference_fasta_web()
		return DatasetConsensusTable.EMPTY

	def render_results(self, record):
		"""
		icon with link to extra info
		"""
		# from crequest.middleware import CrequestMiddleware
		# current_request = CrequestMiddleware.get_request()
		# user = current_request.user
		# str_links = ""
		# if (user.username != Constants.USER_ANONYMOUS):
		#	 str_links = '<a href=' + reverse('sample-project-settings', args=[record.pk]) + ' data-toggle="tooltip" title="Software settings">' +\
		#		 '<span ><i class="padding-button-table fa fa-magic padding-button-table"></i></span></a>'
		#
		# b_space_add = False
		# if (record.is_mask_consensus_sequences):
		#	 str_links += '<a href="javascript:void(0);" data-toggle="tooltip" title="Sample run with user-selected parameters (see update 30 Oct 2020)">' +\
		#		 '<span ><i class="padding-button-table fa fa-calendar-check-o"></i></a>'
		# else:
		#	 b_space_add = True
		#	 str_links += '  '
		#
		# ### check if it was downsized
		# if (self.manage_database.is_sample_downsized(record.sample)):
		#	 str_links += '<a href="javascript:void(0);" data-toggle="tooltip" title="Sample was downsized">' +\
		#		 '<span ><i class="padding-button-table warning_fa_icon fa fa-level-down"></i></a>'
		# elif (not b_space_add):
		#	 str_links += '  '
		#
		# str_links += '<a href=' + reverse('show-sample-project-single-detail', args=[record.pk]) + ' data-toggle="tooltip" title="Show more information" class="padding-button-table">' +\
		#		 '<span ><i class="fa fa-info-circle"></i> More info</a>'
		# return mark_safe(str_links)
		return DatasetConsensusTable.EMPTY

	def order_sample_name(self, queryset, is_descending):
		queryset = queryset.annotate(sample_name = F('sample__name')).order_by(('-' if is_descending else '') + 'sample_name')
		return (queryset, True)

	def order_type_and_subtype(self, queryset, is_descending):
		queryset = queryset.annotate(type_subtype = F('sample__type_subtype')).order_by(('-' if is_descending else '') + 'type_subtype')
		return (queryset, True)

	def order_alerts(self, queryset, is_descending):
		queryset = queryset.annotate(alerts = F('alert_first_level') + F('alert_second_level')).order_by(('-' if is_descending else '') + 'alerts')
		return (queryset, True)
	
	def order_technology(self, queryset, is_descending):
		""" shows if it is Illumina or Minion """
		queryset = queryset.annotate(technology = F('sample__type_of_fastq')).order_by(('-' if is_descending else '') + 'technology')
		return (queryset, True)


class AddDatasetFromCvsFileTableMetadata(tables.Table):
	"""
	To add metadata to datasets
	"""
	is_processed = tables.Column('Is processed', empty_values=())

	class Meta:
		model = UploadFiles
		fields = ('file_name', 'creation_date', 'owner', 'is_valid', 'is_completed')
		attrs = {"class": "table-striped table-bordered"}
		empty_text = "There are no (csv) or (tsv) files with metadata to add..."
	
	def render_creation_date(self, value, record):
		return record.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)

	def render_is_processed(self, value, record):
		return 'True' if record.is_processed else 'False'
	
	def render_file_name(self, value, record):
		"""
		return file name
		"""
		return record.get_metadata_fasta_web()
	
	def order_is_processed(self, queryset, is_descending):
		queryset = queryset.annotate(is_processed = F('is_processed')).order_by(('-' if is_descending else '') + 'is_processed')
		return (queryset, True)
	
