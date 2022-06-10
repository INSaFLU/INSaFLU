
import logging, datetime, ntpath, os
from django.urls import reverse_lazy
from django.views import generic
from braces.views import LoginRequiredMixin, FormValidMessageMixin
from django.conf import settings
from managing_files.models import Reference, Project, ProjectSample
from datasets.forms import ConsensusForm
from django.views.generic import ListView
from utils.utils import Utils
from datasets.models import Dataset, DatasetConsensus, Consensus
from datasets.tables import DatasetTable, ReferenceTable, ConsensusTable, ProjectTable, DatasetConsensusTable
from datasets.forms import AddReferencesDatasetForm, AddConsensusDatasetForm, AddProjectsDatasetForm
from django.http import HttpResponseRedirect
from extend_user.models import Profile
from django.contrib import messages
from django.db.models import Q
from django_tables2 import RequestConfig
from constants.constants import Constants, TypePath
from utils.utils import ShowInfoMainPage
from django.template.defaultfilters import pluralize
from operator import attrgetter
from itertools import chain
from utils.session_variables import clean_check_box_in_session, is_all_check_box_in_session
from utils.software import Software
from django.db import transaction



class DatasetsView(LoginRequiredMixin, ListView):
    """
    List all datasets
    """
    
    utils = Utils()
    model = Dataset
    template_name = 'datasets/datasets.html'
    context_object_name = 'datasets'
##    group_required = u'company-user' security related with GroupRequiredMixin
    
    def get_context_data(self, **kwargs):
        context = super(DatasetsView, self).get_context_data(**kwargs)
        tag_search = 'search_projects'
        query_set = Dataset.objects.filter(owner__id=self.request.user.id, is_deleted=False).order_by('-creation_date')
        if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)):
            query_set = query_set.filter(Q(name__icontains=self.request.GET.get(tag_search))).\
                            distinct()
                            
        table = DatasetTable(query_set)
        RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
        if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
        
        ### clean check box in the session
        clean_check_box_in_session(self.request) ## for both samples and references

        ### clean project name session
        if Constants.PROJECT_NAME_SESSION in self.request.session:
            del self.request.session[Constants.PROJECT_NAME_SESSION]

        context['table'] = table
        context['nav_dataset'] = True
        context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
        context['query_set_count'] = query_set.count()
        context['show_info_main_page'] = ShowInfoMainPage()        ## show main information about the institute
        return context



class AddDatasetsReferencesView(LoginRequiredMixin, FormValidMessageMixin, generic.CreateView):
    """
    Add references to Dataset (one or more)
    """
    utils = Utils()
    model = Dataset
    fields = ['name']
    template_name = 'datasets/datasets_references.html'
    success_url = reverse_lazy('datasets')
    context_object_name = 'datasets'

    if settings.DEBUG: logger = logging.getLogger("fluWebVirus.debug")
    else: logger = logging.getLogger("fluWebVirus.production")

    def get_context_data(self, **kwargs):
        context = super(AddDatasetsReferencesView, self).get_context_data(**kwargs)
        
        ### test if the user is the same of the page
        dataset = Dataset.objects.get(pk=self.kwargs['pk'])
        context['nav_dataset'] = True
        if (dataset.owner.id != self.request.user.id): 
            context['error_cant_see'] = "1"
            return context
        
        tag_search = 'search_references'
        query_set = Reference.objects.filter(owner__id=self.request.user.id, is_obsolete=False, is_deleted=False).order_by('-name')
        if (not self.request.GET.get(tag_search) is None and self.request.GET.get(tag_search)): 
            query_set = query_set.filter(Q(name__icontains=self.request.GET.get(tag_search)) |\
                                        Q(owner__username__icontains=self.request.GET.get(tag_search)) |\
                                        Q(reference_genbank_name__icontains=self.request.GET.get(tag_search)) |\
                                        Q(reference_fasta_name__icontains=self.request.GET.get(tag_search)) |\
                                        Q(isolate_name__icontains=self.request.GET.get(tag_search)))
        
        ### get the references from the system
        query_set_system = Reference.objects.filter(owner__username=Constants.DEFAULT_USER, is_obsolete=False, is_deleted=False).order_by('-name')
        if (not self.request.GET.get(tag_search) is None and self.request.GET.get(tag_search)): 
            query_set_system = query_set_system.filter(Q(name__icontains=self.request.GET.get(tag_search)) |\
                                    Q(owner__username__icontains=self.request.GET.get(tag_search)) |\
                                    Q(reference_genbank_name__icontains=self.request.GET.get(tag_search)) |\
                                    Q(reference_fasta_name__icontains=self.request.GET.get(tag_search)) |\
                                    Q(isolate_name__icontains=self.request.GET.get(tag_search)))
            
        ### get references all already set
        list_ids = Reference.objects.filter(dataset_reference__dataset=dataset, dataset_reference__is_deleted=False).values('pk')
         
        ### final dataset
        query_set_result = sorted(chain(query_set.exclude(id__in=list_ids), \
                query_set_system.exclude(id__in=list_ids)), key=attrgetter('creation_date'))
        
        ## there's no references
        if len(query_set_result) == 0:
            context['error_cant_see'] = "1"
            return context
        
        ## start table building
        table = ReferenceTable(query_set_result)
        RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
        if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
        
        ### set the check_box
        if (Constants.CHECK_BOX_ALL not in self.request.session):
            self.request.session[Constants.CHECK_BOX_ALL] = False
            is_all_check_box_in_session(["{}_{}".format(Constants.CHECK_BOX, key.id) for key in query_set_result], self.request)

        context[Constants.CHECK_BOX_ALL] = self.request.session[Constants.CHECK_BOX_ALL]
        ## need to clean all the others if are reject in filter
        dt_reference_id_add_temp = {}
        if (context[Constants.CHECK_BOX_ALL]):
            for reference in query_set_result: dt_reference_id_add_temp[reference.id] = 1    ## add the ids that are in the tables
            for key in self.request.session.keys():
                if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and self.utils.is_integer(key.split('_')[2])):
                    ### this is necessary because of the search. Can occur some checked box that are out of filter.
                    if (int(key.split('_')[2]) not in dt_reference_id_add_temp):
                        self.request.session[key] = False
                    else: self.request.session[key] = True
        ## END need to clean all the others if are reject in filter

        context['table'] = table
        context['show_paginatior'] = len(query_set_result) > Constants.PAGINATE_NUMBER
        context['query_set_count'] = len(query_set_result)
        context['dataset_name'] = dataset.name
        context['add_all_references_message'] = "Add {} reference{}".format(len(query_set_result), pluralize(len(query_set_result), ',s'))
        context['show_info_main_page'] = ShowInfoMainPage()        ## show main information about the institute
        
        ## Add references to DataSet
        if self.request.POST: 
            context['reference_dataset'] = AddReferencesDatasetForm(self.request.POST)
        else: 
            context['reference_dataset'] = AddReferencesDatasetForm()
        return context

    def form_valid(self, form):
        """
        Validate the form
        """

        ### test anonymous account
        try:
            profile = Profile.objects.get(user=self.request.user)
            if (profile.only_view_project):
                messages.warning(self.request, "'{}' account can not add references to a dataset.".format(self.request.user.username), fail_silently=True)
                return super(AddDatasetsReferencesView, self).form_invalid(form)
        except Profile.DoesNotExist:
            return super(AddDatasetsReferencesView, self).form_invalid(form)
        
        if form.is_valid():
            
            ### get project sample..
            context = self.get_context_data()
    
            ### get project
            try:
                dataset = Dataset.objects.get(pk=self.kwargs['pk'])
            except Dataset.DoesNotExist:
                return super(AddDatasetsReferencesView, self).form_invalid(form)
            
            vect_sample_id_add_temp = []
            for sample in context['table'].data:
                vect_sample_id_add_temp.append(sample.id)
                            
            vect_sample_id_add = []
            if ("submit_checked" in self.request.POST):
                for key in self.request.session.keys():
                    if (self.request.session[key] and key.startswith(Constants.CHECK_BOX) and\
                        len(key.split('_')) == 3 and self.utils.is_integer(key.split('_')[2])):
                        ### this is necessary because of the search. Can occur some checked box that are out of filter.
                        if (int(key.split('_')[2]) in vect_sample_id_add_temp):
                            vect_sample_id_add.append(int(key.split('_')[2]))
            elif ("submit_all" in self.request.POST):
                vect_sample_id_add = vect_sample_id_add_temp
            
            ### start adding...
            reference_add = 0
            for id_reference in vect_sample_id_add:
                
                try:
                    dataset_consensus = DatasetConsensus.objects.get(reference__pk=id_reference,
                                        dataset=dataset)
                    
                    if dataset_consensus.is_deleted or dataset_consensus.is_error:
                        dataset_consensus.is_deleted = False
                        dataset_consensus.is_error = False
                        dataset_consensus.save()
                        reference_add += 1
                    continue
                except DatasetConsensus.DoesNotExist:
                    
                    try:
                        reference = Reference.objects.get(pk=id_reference)
                        dataset_consensus = DatasetConsensus()
                        dataset_consensus.name = reference.name
                        dataset_consensus.dataset = dataset
                        dataset_consensus.reference = reference
                        dataset_consensus.save() 
                        reference_add += 1
                    except reference.DoesNotExist:
                        self.logger.error('Fail to get reference_id {} in DataSet {}'.format(id_reference, dataset.name))
                    
            ### necessary to calculate the global results again 
            if (reference_add > 0):
                dataset.last_change_date = datetime.datetime.now()
                dataset.number_of_sequences_from_references += reference_add
                dataset.totla_alerts = 1 if dataset.get_number_different_references() > 1 else 0
                dataset.save()
            
            if (reference_add == 0):
                messages.warning(self.request, "No references was added to the data set '{}'".format(dataset.name))
            else:
                if (reference_add > 1):
                    messages.success(self.request, "'{}' references were added to your data set '{}'.".format(\
                                reference_add, dataset.name), fail_silently=True)
                else:
                    messages.success(self.request, "One reference was added to your data set '{}'.".format(dataset.name), fail_silently=True)
            return HttpResponseRedirect(reverse_lazy('datasets'))
        else:
            return super(AddDatasetsReferencesView, self).form_invalid(form)

    form_valid_message = ""        ## need to have this, even empty


class AddDatasetsConsensusView(LoginRequiredMixin, FormValidMessageMixin, generic.CreateView):
    """
    Add consensus to Dataset, can add consensus to his/her account
    """
    utils = Utils()
    model = Dataset
    fields = ['name']
    template_name = 'datasets/datasets_consensus.html'
    success_url = reverse_lazy('datasets')
    context_object_name = 'datasets'
    
    if settings.DEBUG: logger = logging.getLogger("fluWebVirus.debug")
    else: logger = logging.getLogger("fluWebVirus.production")
    
    def get_context_data(self, **kwargs):
        context = super(AddDatasetsConsensusView, self).get_context_data(**kwargs)
        
        ### test if the user is the same of the page
        dataset = Dataset.objects.get(pk=self.kwargs['pk'])
        context['nav_dataset'] = True
        if (dataset.owner.id != self.request.user.id): 
            context['error_cant_see'] = "1"
            return context
        
        tag_search = 'search_consensus'
        query_set = Consensus.objects.filter(owner__id=self.request.user.id, is_deleted=False).order_by('-name')
        if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)):
            query_set = query_set.filter(Q(name__icontains=self.request.GET.get(tag_search))).\
                            distinct()
        
        ### get consensus all already set
        list_ids = Consensus.objects.filter(dataset_consensus__dataset=dataset, dataset_consensus__is_deleted=False).values('pk')
         
        ### final dataset
        query_set_result = sorted(query_set.exclude(id__in=list_ids), key=attrgetter('creation_date'))
        
        table = ConsensusTable(query_set_result)
        RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
        if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
        
        ### set the check_box
        if (Constants.CHECK_BOX_ALL not in self.request.session):
            self.request.session[Constants.CHECK_BOX_ALL] = False
            is_all_check_box_in_session(["{}_{}".format(Constants.CHECK_BOX, key.id) for key in query_set_result], self.request)

        context[Constants.CHECK_BOX_ALL] = self.request.session[Constants.CHECK_BOX_ALL]
        ## need to clean all the others if are reject in filter
        dt_reference_id_add_temp = {}
        if (context[Constants.CHECK_BOX_ALL]):
            for consensus in query_set_result: dt_reference_id_add_temp[consensus.id] = 1    ## add the ids that are in the tables
            for key in self.request.session.keys():
                if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and self.utils.is_integer(key.split('_')[2])):
                    ### this is necessary because of the search. Can occur some checked box that are out of filter.
                    if (int(key.split('_')[2]) not in dt_reference_id_add_temp):
                        self.request.session[key] = False
                    else: self.request.session[key] = True
        ## END need to clean all the others if are reject in filter

        context['pk'] = dataset.pk
        context['table'] = table
        context['show_paginatior'] = len(query_set_result) > Constants.PAGINATE_NUMBER
        context['query_set_count'] = len(query_set_result)
        context['dataset_name'] = dataset.name
        context['add_all_references_message'] = "Add {} consensus{}".format(len(query_set_result), pluralize(len(query_set_result), ',s'))
        context['show_info_main_page'] = ShowInfoMainPage()        ## show main information about the institute
        
        ## Add references to DataSet
        if self.request.POST: 
            context['consensus_dataset'] = AddConsensusDatasetForm(self.request.POST)
        else: 
            context['consensus_dataset'] = AddConsensusDatasetForm()
        return context
    
    def form_valid(self, form):
        """
        Validate the form
        """

        ### test anonymous account
        try:
            profile = Profile.objects.get(user=self.request.user)
            if (profile.only_view_project):
                messages.warning(self.request, "'{}' account can not add Consensus to a dataset.".format(self.request.user.username), fail_silently=True)
                return super(AddDatasetsConsensusView, self).form_invalid(form)
        except Profile.DoesNotExist:
            return super(AddDatasetsConsensusView, self).form_invalid(form)
        
        if form.is_valid():
            
            ### get project sample..
            context = self.get_context_data()
    
            ### get project
            try:
                dataset = Dataset.objects.get(pk=self.kwargs['pk'])
            except Dataset.DoesNotExist:
                return super(AddDatasetsConsensusView, self).form_invalid(form)
            
            vect_sample_id_add_temp = []
            for sample in context['table'].data:
                vect_sample_id_add_temp.append(sample.id)
                            
            vect_sample_id_add = []
            if ("submit_checked" in self.request.POST):
                for key in self.request.session.keys():
                    if (self.request.session[key] and key.startswith(Constants.CHECK_BOX) and\
                        len(key.split('_')) == 3 and self.utils.is_integer(key.split('_')[2])):
                        ### this is necessary because of the search. Can occur some checked box that are out of filter.
                        if (int(key.split('_')[2]) in vect_sample_id_add_temp):
                            vect_sample_id_add.append(int(key.split('_')[2]))
            elif ("submit_all" in self.request.POST):
                vect_sample_id_add = vect_sample_id_add_temp
            
            ### start adding...
            reference_add = 0
            for id_consensus in vect_sample_id_add:
                
                try:
                    dataset_consensus = DatasetConsensus.objects.get(consensus__pk=id_consensus,
                                        dataset=dataset)
                    
                    if dataset_consensus.is_deleted or dataset_consensus.is_error:
                        dataset_consensus.is_deleted = False
                        dataset_consensus.is_error = False
                        dataset_consensus.save()
                        reference_add += 1
                    continue
                except DatasetConsensus.DoesNotExist:
                    
                    try:
                        consensus = Consensus.objects.get(pk=id_consensus)
                        dataset_consensus = DatasetConsensus()
                        dataset_consensus.name = consensus.name
                        dataset_consensus.dataset = dataset
                        dataset_consensus.consensus = consensus
                        dataset_consensus.save() 
                        reference_add += 1
                    except consensus.DoesNotExist:
                        self.logger.error('Fail to get consensus_id {} in DataSet {}'.format(id_consensus, dataset.name))
                    
            ### necessary to calculate the global results again 
            if (reference_add > 0):
                dataset.last_change_date = datetime.datetime.now()
                dataset.number_of_sequences_from_consensus += reference_add
                dataset.save()
            
            if (reference_add == 0):
                messages.warning(self.request, "No consensus was added to the data set '{}'".format(dataset.name))
            else:
                if (reference_add > 1):
                    messages.success(self.request, "'{}' consensus were added to your data set '{}'.".format(\
                                reference_add, dataset.name), fail_silently=True)
                else:
                    messages.success(self.request, "One consensus was added to your data set '{}'.".format(dataset.name), fail_silently=True)
            return HttpResponseRedirect(reverse_lazy('datasets'))
        else:
            return super(AddDatasetsConsensusView, self).form_invalid(form)

    form_valid_message = ""        ## need to have this, even empty


class AddDatasetsProjectsView(LoginRequiredMixin, FormValidMessageMixin, generic.CreateView):
    """
    Add new projects/consensus to Dataset
    """
    utils = Utils()
    model = Dataset
    fields = ['name']
    template_name = 'datasets/datasets_projects.html'
    success_url = reverse_lazy('datasets')
    context_object_name = 'datasets'
    
    if settings.DEBUG: logger = logging.getLogger("fluWebVirus.debug")
    else: logger = logging.getLogger("fluWebVirus.production")
##    group_required = u'company-user' security related with GroupRequiredMixin
    
    def get_context_data(self, **kwargs):
        context = super(AddDatasetsProjectsView, self).get_context_data(**kwargs)
        
        ### test if the user is the same of the page
        dataset = Dataset.objects.get(pk=self.kwargs['pk'])
        context['nav_dataset'] = True
        if (dataset.owner.id != self.request.user.id): 
            context['error_cant_see'] = "1"
            return context
        
        tag_search = 'search_consensus'
        query_set = Project.objects.filter(owner__id=self.request.user.id, is_deleted=False,\
                                    number_passed_sequences__gt=0).order_by('-name')
        if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)):
            query_set = query_set.filter(Q(name__icontains=self.request.GET.get(tag_search))).\
                            distinct()
        
        ### get consensus all already set
        list_ids = Project.objects.filter(project_samples__dataset_project_samples__dataset=dataset,\
                    project_samples__dataset_project_samples__is_deleted=False).\
                    distinct().values('pk')
         
        ### final dataset
        query_set_result = sorted(query_set.exclude(id__in=list_ids), key=attrgetter('creation_date'))
        
        table = ProjectTable(query_set_result)
        RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
        if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
        
        ### set the check_box
        if (Constants.CHECK_BOX_ALL not in self.request.session):
            self.request.session[Constants.CHECK_BOX_ALL] = False
            is_all_check_box_in_session(["{}_{}".format(Constants.CHECK_BOX, key.id) for key in query_set_result], self.request)

        context[Constants.CHECK_BOX_ALL] = self.request.session[Constants.CHECK_BOX_ALL]
        ## need to clean all the others if are reject in filter
        dt_reference_id_add_temp = {}
        if (context[Constants.CHECK_BOX_ALL]):
            for consensus in query_set_result: dt_reference_id_add_temp[consensus.id] = 1    ## add the ids that are in the tables
            for key in self.request.session.keys():
                if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and self.utils.is_integer(key.split('_')[2])):
                    ### this is necessary because of the search. Can occur some checked box that are out of filter.
                    if (int(key.split('_')[2]) not in dt_reference_id_add_temp):
                        self.request.session[key] = False
                    else: self.request.session[key] = True
        ## END need to clean all the others if are reject in filter

        context['pk'] = dataset.pk
        context['table'] = table
        context['show_paginatior'] = len(query_set_result) > Constants.PAGINATE_NUMBER
        context['query_set_count'] = len(query_set_result)
        context['dataset_name'] = dataset.name
        context['add_all_references_message'] = "Add {} project{}".format(len(query_set_result), pluralize(len(query_set_result), ',s'))
        context['show_info_main_page'] = ShowInfoMainPage()        ## show main information about the institute
        
        ## Add references to DataSet
        if self.request.POST: 
            context['consensus_dataset'] = AddProjectsDatasetForm(self.request.POST)
        else: 
            context['consensus_dataset'] = AddProjectsDatasetForm()
        return context


    def form_valid(self, form):
        """
        Validate the form
        """
        ### test anonymous account
        try:
            profile = Profile.objects.get(user=self.request.user)
            if (profile.only_view_project):
                messages.warning(self.request, "'{}' account can not add Projects to a dataset.".format(self.request.user.username), fail_silently=True)
                return super(AddDatasetsProjectsView, self).form_invalid(form)
        except Profile.DoesNotExist:
            return super(AddDatasetsProjectsView, self).form_invalid(form)
        
        if form.is_valid():
            
            ### get project sample..
            context = self.get_context_data()
    
            ### get dataset
            try:
                dataset = Dataset.objects.get(pk=self.kwargs['pk'])
            except Dataset.DoesNotExist:
                return super(AddDatasetsProjectsView, self).form_invalid(form)
            
            vect_sample_id_add_temp = []
            for sample in context['table'].data:
                vect_sample_id_add_temp.append(sample.id)
                            
            vect_sample_id_add = []
            if ("submit_checked" in self.request.POST):
                for key in self.request.session.keys():
                    if (self.request.session[key] and key.startswith(Constants.CHECK_BOX) and\
                        len(key.split('_')) == 3 and self.utils.is_integer(key.split('_')[2])):
                        ### this is necessary because of the search. Can occur some checked box that are out of filter.
                        if (int(key.split('_')[2]) in vect_sample_id_add_temp):
                            vect_sample_id_add.append(int(key.split('_')[2]))
            elif ("submit_all" in self.request.POST):
                vect_sample_id_add = vect_sample_id_add_temp
            
            ### start adding...
            reference_add, consensus_add = 0, 0
            for id_project in vect_sample_id_add:
                
                try:
                    project = Project.objects.get(pk=id_project)
                except Project.DoesNotExist:
                    self.logger.error("Project with id '{}' does not exist".format(id_project))
                    continue
                
                ## run hover all project samples
                for project_sample in ProjectSample.objects.filter(project=project, is_deleted= False):
                    ## don't add project samples not added
                    if not project_sample.is_finished: continue
                    try:
                        dataset_consensus = DatasetConsensus.objects.get(project_sample=project_sample,
                                            dataset=dataset)
                        
                        if dataset_consensus.is_deleted or dataset_consensus.is_error:
                            dataset_consensus.is_deleted = False
                            dataset_consensus.is_error = False
                            dataset_consensus.save()
                            consensus_add += 1
                        ## now is finished and before is not finished yet
                        elif (project_sample.is_finished and not dataset_consensus.is_project_sample_finished):
                            dataset_consensus.is_project_sample_finished = True
                            consensus_add += 1
                        continue
                    except DatasetConsensus.DoesNotExist:
                        dataset_consensus = DatasetConsensus()
                        dataset_consensus.name = project_sample.sample.name
                        dataset_consensus.type_subtype = project_sample.sample.type_subtype
                        dataset_consensus.dataset = dataset
                        dataset_consensus.project_sample = project_sample
                        dataset_consensus.is_project_sample_finished = project_sample.is_finished
                        dataset_consensus.save() 
                        consensus_add += 1
                
                ### Add the reference of this project if not there yet
                try:
                    dataset_consensus = DatasetConsensus.objects.get(reference=project.reference,
                                        dataset=dataset)
                    if dataset_consensus.is_deleted or dataset_consensus.is_error:
                        dataset_consensus.is_deleted = False
                        dataset_consensus.is_error = False
                        dataset_consensus.save()
                        reference_add += 1
                except DatasetConsensus.DoesNotExist:
                    dataset_consensus = DatasetConsensus()
                    dataset_consensus.name = project.reference.name
                    dataset_consensus.dataset = dataset
                    dataset_consensus.reference = project.reference
                    dataset_consensus.save() 
                    reference_add += 1

            ### necessary to calculate the global results again 
            if (reference_add > 0):
                dataset.last_change_date = datetime.datetime.now()
                dataset.number_of_sequences_from_projects += consensus_add
                dataset.number_of_sequences_from_references += reference_add
                dataset.totla_alerts = 1 if dataset.get_number_different_references() > 1 else 0
                dataset.save()
            
            if (reference_add == 0):
                messages.warning(self.request, "No consensus was added to the data set '{}'".format(dataset.name))
            else:
                if (reference_add > 1):
                    messages.success(self.request, "'{}' projects with samples were added to your data set '{}'.".format(\
                                reference_add, dataset.name), fail_silently=True)
                else:
                    messages.success(self.request, "One consensus from a sample was added to your data set '{}'.".format(dataset.name), fail_silently=True)
            return HttpResponseRedirect(reverse_lazy('datasets'))
        else:
            return super(AddDatasetsProjectsView, self).form_invalid(form)

    form_valid_message = ""        ## need to have this, even empty

class ShowDatasetsSequencesView(LoginRequiredMixin, ListView):
    """
    Show sequences of Dataset
    """
    utils = Utils()
    model = Dataset
    template_name = 'datasets/datasets.html'
    context_object_name = 'datasets'
##    group_required = u'company-user' security related with GroupRequiredMixin
    
    def get_context_data(self, **kwargs):
        context = super(ShowDatasetsSequencesView, self).get_context_data(**kwargs)
        tag_search = 'search_projects'
        query_set = Dataset.objects.filter(owner__id=self.request.user.id, is_deleted=False).order_by('-creation_date')
        if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)):
            query_set = query_set.filter(Q(name__icontains=self.request.GET.get(tag_search))).\
                            distinct()
                            
        table = DatasetTable(query_set)
        RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
        if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
        
        ### clean check box in the session
        #clean_check_box_in_session(self.request) ## for both samples and references

        ### clean project name session
        if Constants.PROJECT_NAME_SESSION in self.request.session:
            del self.request.session[Constants.PROJECT_NAME_SESSION]

        context['table'] = table
        context['nav_dataset'] = True
        context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
        context['query_set_count'] = query_set.count()
        
        context['show_info_main_page'] = ShowInfoMainPage()        ## show main information about the institute
        return context

class UploadNewConsensusView(LoginRequiredMixin, FormValidMessageMixin, generic.FormView):
    
    form_class = ConsensusForm
    template_name = 'datasets/datasets_new_consensus.html'

    ## Other solution to get the Consensus
    ## https://pypi.python.org/pypi?%3aaction=display&name=django-contrib-requestprovider&version=1.0.1
    def get_form_kwargs(self):
        """
        Set the request to pass in the form
        """
        kw = super(UploadNewConsensusView, self).get_form_kwargs()
        kw['request'] = self.request # the trick!
        kw['pk'] = self.kwargs['pk']
        return kw
    
    def get_success_url(self):
        """
        get source_pk from consensus dataset, need to pass it in context
        """
        return reverse_lazy('add-consensus-dataset', kwargs={'pk': self.kwargs.get('pk')})


    def get_context_data(self, **kwargs):
        context = super(UploadNewConsensusView, self).get_context_data(**kwargs)
        context['nav_dataset'] = True
        context['nav_modal'] = True    ## short the size of modal window
        context['show_info_main_page'] = ShowInfoMainPage()        ## show main information about the institute
        return context
    
    @transaction.atomic
    def form_valid(self, form):
        software = Software()
        utils = Utils()
        
        ### test anonymous account
        try:
            profile = Profile.objects.get(user=self.request.user)
            if (profile.only_view_project):
                messages.warning(self.request, "'{}' account can not add references.".format(self.request.user.username), fail_silently=True)
                return super(UploadNewConsensusView, self).form_invalid(form)
        except Profile.DoesNotExist:
            pass
        
        if form.is_valid():
            name = form.cleaned_data['name']
            consensus_fasta = form.cleaned_data['consensus_fasta']
            
            consensus = form.save(commit=False)
            ## set other data
            consensus.display_name = consensus.name
            consensus.owner = self.request.user
            consensus.is_obsolete = False
            consensus.number_of_locus = self.request.session[Constants.NUMBER_LOCUS_FASTA_FILE]
            consensus.consensus_fasta_name = utils.clean_name(ntpath.basename(consensus_fasta.name))
            consensus.save()
    
            ## move the files to the right place
            sz_file_to = os.path.join(settings.MEDIA_ROOT, utils.get_path_to_consensus_file(self.request.user.id, consensus.id), consensus.consensus_fasta_name)
            software.dos_2_unix(os.path.join(settings.MEDIA_ROOT, consensus.consensus_fasta.name))
            ## test if bases all lower
            software.fasta_2_upper(os.path.join(settings.MEDIA_ROOT, consensus.consensus_fasta.name))
            utils.move_file(os.path.join(settings.MEDIA_ROOT, consensus.consensus_fasta.name), sz_file_to)
            consensus.hash_reference_fasta = utils.md5sum(sz_file_to)
            consensus.consensus_fasta.name = os.path.join(utils.get_path_to_consensus_file(self.request.user.id, consensus.id), consensus.consensus_fasta_name)
            consensus.save()
            
            ## create the index before commit in database, throw exception if something goes wrong
            # don't create this index, not necessary
            # do not test the same sequence length in from validation
            #software.create_fai_fasta(os.path.join(getattr(settings, "MEDIA_ROOT", None), consensus.consensus_fasta.name))
            
            messages.success(self.request, "Consensus '" + name + "' was created successfully", fail_silently=True)
            return super(UploadNewConsensusView, self).form_valid(form)
        else:
            return super(UploadNewConsensusView, self).form_invalid(form)

    form_valid_message = ""        ## need to have this


class ShowDatasetsConsensusView(LoginRequiredMixin, ListView):
    model = Project
    template_name = 'datasets/dataset_consensus_show.html'
    context_object_name = 'dataset_consensus'
    
    def get_context_data(self, **kwargs):
        context = super(ShowDatasetsConsensusView, self).get_context_data(**kwargs)
        dataset = Dataset.objects.get(pk=self.kwargs['pk'])
        
        ### can't see this project
        context['nav_dataset'] = True
        if (dataset.owner.id != self.request.user.id): 
            context['error_cant_see'] = "1"
            return context
        
        query_set = DatasetConsensus.objects.filter(dataset=dataset, is_deleted=False, is_error=False).order_by('creation_date')
        tag_search = 'search_add_project_sample'
        ### filter the search
        if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)): 
            query_set = query_set.filter(Q(name__icontains=self.request.GET.get(tag_search)) |\
                                        Q(type_subtype__icontains=self.request.GET.get(tag_search)))
        table = DatasetConsensusTable(query_set)
        RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)

        if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get('search_add_project_sample')        
        context['table'] = table
        context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
        context['query_set_count'] = query_set.count()
        context['dataset'] = dataset
        context['different_references'] = dataset.get_number_different_references()
        context['spinner_url'] = os.path.join("/" + Constants.DIR_STATIC, Constants.DIR_ICONS, Constants.AJAX_LOADING_GIF)
        context['show_info_main_page'] = ShowInfoMainPage()        ## show main information about the institute

        ## Files
            
        return context
