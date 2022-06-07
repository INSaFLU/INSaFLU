'''
Created on 22/05/2022

@author: mmp
'''
import django_tables2 as tables
from django.urls import reverse
from django.conf import settings
from django.db.models import F
from datasets.models import Dataset, DatasetConsensus
from django.utils.safestring import mark_safe
from constants.constants import Constants, TypePath
from managing_files.models import Reference, Project, ProjectSample

class DatasetTable(tables.Table):
#   Renders a normal value as an internal hyperlink to another page.
#   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
    ### P -> sequences from projects
    ### C -> sequences from consensus
    ### R -> sequences from references
    sequences = tables.Column('#Sequences (P/C/R)', orderable=False, empty_values=())
    last_change_date = tables.Column('Last Change date', empty_values=())
    creation_date = tables.Column('Creation date', empty_values=())
    results = tables.LinkColumn('Options', orderable=False, empty_values=())
    
    class Meta:
        model = Dataset
        fields = ('name', 'last_change_date','creation_date', 'sequences', 'results')
        attrs = {"class": "table-striped table-bordered"}
        empty_text = "There are no Datasets to show..."
    
    def render_name(self, record):
        from crequest.middleware import CrequestMiddleware
        current_request = CrequestMiddleware.get_request()
        user = current_request.user
    
        ## there's nothing to show
        count = record.number_of_sequences_from_projects + record.number_of_sequences_from_consensus + record.number_of_sequences_from_references
        if (count > 0):
            project_sample = '<a href=' + reverse('show-sequences-dataset', args=[record.pk]) + ' data-toggle="tooltip" title="Show sequences">' +\
                '{}</a>'.format(record.name)
        else:
            project_sample = record.name
        
        if (user.username == Constants.USER_ANONYMOUS): return mark_safe(project_sample);
        if (user.username == record.owner.username):
            return mark_safe('<a href="#id_remove_modal" id="id_remove_dataset_modal" data-toggle="modal" data-toggle="tooltip" title="Delete"' +\
                    ' ref_name="' + record.name + '" pk="' + str(record.pk) + '"><i class="fa fa-trash"></i></span> </a>' + project_sample)
        return project_sample;
    
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
        ## there's nothing to show
        count = record.number_of_sequences_from_projects + record.number_of_sequences_from_consensus + record.number_of_sequences_from_references
        sz_project_sample = ""
        if (count > 0):
            sz_project_sample = '<a href=' + reverse('show-sequences-dataset', args=[record.pk]) + ' data-toggle="tooltip" title="See sequences"> ' +\
                '<span ><i class="padding-button-table fa fa-info-circle padding-button-table"></i></span></a> '
            ## Can also launch all the sequeneces in the NextTrain
            # sz_project_sample += '<a href=' + reverse('project-settings', args=[record.pk]) + ' data-toggle="tooltip" title="Software settings">' +\
            #    '<span ><i class="fa fa-magic padding-button-table"></i></span></a>'
            # sz_project_sample += '<a href=' + reverse('project-settings', args=[record.pk]) + ' data-toggle="tooltip" title="Software settings">' +\
            #    '<span ><i class="padding-button-table fa fa-magic padding-button-table"></i></span></a>'
        # else:    ## can change settings
        #     sz_project_sample += '<a href=' + reverse('project-settings', args=[record.pk]) + ' data-toggle="tooltip" title="Software settings">' +\
        #         '<span ><i class="padding-button-table fa fa-magic padding-button-table"></i></span></a>'
    
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
        href = record.get_reference_fasta(TypePath.MEDIA_URL)        
        return mark_safe('<a href="' + href + '" download>' + self.constants.short_name(record.reference_fasta_name, self.SHORT_NAME_LENGTH) + '</a>')

    def render_reference_genbank_name(self, **kwargs):
        record = kwargs.pop("record")
        href = record.get_reference_gbk(TypePath.MEDIA_URL)        
        return mark_safe('<a href="' + href + '" download>' + self.constants.short_name(record.reference_genbank_name, self.SHORT_NAME_LENGTH) + '</a>')

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
        href = record.get_consensus_fasta(TypePath.MEDIA_URL)        
        return mark_safe('<a href="' + href + '" download>' + self.constants.short_name(record.consensus_fasta_name, self.SHORT_NAME_LENGTH) + '</a>')

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
        # count_not_finished = ProjectSample.objects.filter(project__id=record.id, is_deleted=False, is_error=False, is_finished=False).count()
        #
        # sz_project_sample = ""
        # if (count > 0):
        #     sz_project_sample = '<a href=' + reverse('show-sample-project-results', args=[record.pk]) + ' data-toggle="tooltip" title="See Results"> ' +\
        #         '<span ><i class="padding-button-table fa fa-info-circle padding-button-table"></i></span></a> '
        #     ## only can change settings when has projects finished
        #     #sz_project_sample += '<a href=' + reverse('project-settings', args=[record.pk]) + ' data-toggle="tooltip" title="Software settings">' +\
        #     #    '<span ><i class="fa fa-magic padding-button-table"></i></span></a>'
        #     sz_project_sample += '<a href=' + reverse('project-settings', args=[record.pk]) + ' data-toggle="tooltip" title="Software settings">' +\
        #         '<span ><i class="padding-button-table fa fa-magic padding-button-table"></i></span></a>'
        # elif (count_not_finished > 0): 
        #     sz_project_sample = _("{} processing ".format(count_not_finished))
        # else:    ## can change settings
        #     sz_project_sample += '<a href=' + reverse('project-settings', args=[record.pk]) + ' data-toggle="tooltip" title="Software settings">' +\
        #         '<span ><i class="padding-button-table fa fa-magic padding-button-table"></i></span></a>'
        #
        # return mark_safe(sz_project_sample)
        return ""

