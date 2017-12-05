# Create your views here.

from django.views import generic
from braces.views import LoginRequiredMixin, FormValidMessageMixin
from django.core.urlresolvers import reverse_lazy
from django.views.generic import ListView, DetailView
from django_tables2 import RequestConfig
from managing_files.models import Reference, Sample, Project, ProjectSample
from managing_files.tables import ReferenceTable, SampleTable, ProjectTable, ReferenceProjectTable, SampleToProjectsTable, ShowProjectSamplesResults
from managing_files.forms import ReferenceForm, SampleForm, ReferenceProjectFormSet, AddSampleProjectForm
from managing_files.manage_database import ManageDatabase
from constants.constants import Constants, TypePath, FileExtensions
from constants.software_names import SoftwareNames
from constants.meta_key_and_values import MetaKeyAndValue
from utils.software import Software
from utils.collect_extra_data import CollectExtraData
from utils.utils import Utils
from utils.result import DecodeResult
import hashlib, ntpath, os, logging
from django.contrib import messages
from django.conf import settings
from django.contrib.gis.geos import Point
from django_modalview.generic.base import ModalTemplateView
from django_q.tasks import async
from django.utils.safestring import mark_safe
from django.utils.translation import ugettext_lazy as _
from django.http import JsonResponse
from django.db import transaction
from django.db.models import Q

# http://www.craigderington.me/generic-list-view-with-django-tables/
	
#class ReferencesView(LoginRequiredMixin, GroupRequiredMixin, ListView):
class ReferenceView(LoginRequiredMixin, ListView):
	model = Reference
	template_name = 'references/references.html'
	context_object_name = 'reference'
	ordering = ['id']
##	group_required = u'company-user' security related with GroupRequiredMixin
	
	def get_context_data(self, **kwargs):
		context = super(ReferenceView, self).get_context_data(**kwargs)
		query_set = Reference.objects.filter(owner__id=self.request.user.id).order_by('-name')
		table = ReferenceTable(query_set)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
		context['table'] = table
		context['nav_reference'] = True
		context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
		return context


class ReferenceAddView(LoginRequiredMixin, FormValidMessageMixin, generic.FormView):
	"""
	Create a new reference
	"""
	form_class = ReferenceForm
	success_url = reverse_lazy('references')
	template_name = 'references/reference_add.html'

	## Other solution to get the reference
	## https://pypi.python.org/pypi?%3aaction=display&name=django-contrib-requestprovider&version=1.0.1
	def get_form_kwargs(self):
		"""
		Set the request to pass in the form
		"""
		kw = super(ReferenceAddView, self).get_form_kwargs()
		kw['request'] = self.request # the trick!
		return kw
	
	
	def get_context_data(self, **kwargs):
		context = super(ReferenceAddView, self).get_context_data(**kwargs)
		context['nav_reference'] = True
		context['nav_modal'] = True	## short the size of modal window
		return context
	
	def form_valid(self, form):
		software = Software()
		utils = Utils()
		
		name = form.cleaned_data['name']
		scentific_name = form.cleaned_data['isolate_name']
		reference_fasta = form.cleaned_data['reference_fasta']
		reference_genbank = form.cleaned_data['reference_genbank']
		if (scentific_name is None): scentific_name = ""
		
		hash_value_fasta = hashlib.md5(form.files.get('reference_fasta').read()).hexdigest()
		hash_value_genbank = hashlib.md5(form.files.get('reference_genbank').read()).hexdigest()
		
		reference = form.save(commit=False)
		## set other data
		reference.owner = self.request.user
		reference.is_obsolete = False
		reference.number_of_locus = self.request.session[Constants.NUMBER_LOCUS_FASTA_FILE]
		reference.hash_reference_fasta = hash_value_fasta
		reference.reference_fasta_name = ntpath.basename(reference_fasta.name)
		reference.scentific_name = scentific_name
		reference.hash_reference_genbank = hash_value_genbank
		reference.reference_genbank_name = ntpath.basename(reference_genbank.name)
		reference.save()

		## move the files to the right place
		sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), utils.get_path_to_reference_file(self.request.user.id, reference.id), reference.reference_fasta_name)
		utils.move_file(os.path.join(getattr(settings, "MEDIA_ROOT", None), reference.reference_fasta.name), sz_file_to)
		reference.reference_fasta.name = os.path.join(utils.get_path_to_reference_file(self.request.user.id, reference.id), reference.reference_fasta_name)
		
		sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), utils.get_path_to_reference_file(self.request.user.id, reference.id), reference.reference_genbank_name)
		utils.move_file(os.path.join(getattr(settings, "MEDIA_ROOT", None), reference.reference_genbank.name), sz_file_to)
		reference.reference_genbank.name = os.path.join(utils.get_path_to_reference_file(self.request.user.id, reference.id), reference.reference_genbank_name)
		reference.save()
		
		## create the index before commit in database, throw exception if something goes wrong
		software.create_fai_fasta(reference.reference_fasta.name)
		
		messages.success(self.request, "Reference '" + name + "' was created successfully", fail_silently=True)
		return super(ReferenceAddView, self).form_valid(form)

	## static method, not need for now.
	form_valid_message = ""		## need to have this
	

class SamplesView(LoginRequiredMixin, ListView):
	model = Sample
	template_name = 'samples/samples.html'
	context_object_name = 'samples'
##	group_required = u'company-user' security related with GroupRequiredMixin
	
	def get_context_data(self, **kwargs):
		context = super(SamplesView, self).get_context_data(**kwargs)
		query_set = Sample.objects.filter(owner__id=self.request.user.id).order_by('-creation_date')
		table = SampleTable(query_set)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
		context['table'] = table
		context['nav_sample'] = True
		context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
		return context
	
class SamplesAddView(LoginRequiredMixin, FormValidMessageMixin, generic.FormView):
	"""
	Create a new reference
	"""
	form_class = SampleForm
	success_url = reverse_lazy('samples')
	template_name = 'samples/sample_add.html'


	def get_form_kwargs(self):
		"""
		Set the request to pass in the form
		"""
		kw = super(SamplesAddView, self).get_form_kwargs()
		kw['request'] = self.request 	# the trick!
		return kw


	def get_context_data(self, **kwargs):
		context = super(SamplesAddView, self).get_context_data(**kwargs)
		context['nav_sample'] = True
		context['nav_modal'] = True	## short the size of modal window
		return context


	def form_valid(self, form):
		"""
		Validate the form
		"""
		software = Software()
		utils = Utils()
		
		name = form.cleaned_data['name']
		lat = form.cleaned_data['lat']
		lng = form.cleaned_data['lng']
		like_dates = form.cleaned_data['like_dates']
			
		sample = form.save(commit=False)
		## set other data
		sample.owner = self.request.user
		sample.is_rejected = False
		sample.is_obsolete = False
		sample.file_name_1 = os.path.basename(sample.path_name_1.name)
		sample.is_valid_1 = True
		if (sample.exist_file_2()):
			sample.file_name_2 = os.path.basename(sample.path_name_2.name)
			sample.is_valid_2 = True 
		else: sample.is_valid_2 = False
		sample.has_files = True
		
		if (like_dates == 'date_of_onset'):
			sample.day = int(sample.date_of_onset.strftime("%d"))
			sample.week = int(sample.date_of_onset.strftime("%W")) + 1
			sample.year = int(sample.date_of_onset.strftime("%Y"))
			sample.month = int(sample.date_of_onset.strftime("%m"))
		elif (like_dates == 'date_of_collection'):
			sample.day = int(sample.date_of_collection.strftime("%d"))
			sample.week = int(sample.date_of_collection.strftime("%W")) + 1
			sample.year = int(sample.date_of_collection.strftime("%Y"))
			sample.month = int(sample.date_of_collection.strftime("%m"))
		elif (like_dates == 'date_of_receipt_lab'):
			sample.day = int(sample.date_of_receipt_lab.strftime("%d"))
			sample.week = int(sample.date_of_receipt_lab.strftime("%W")) + 1
			sample.year = int(sample.date_of_receipt_lab.strftime("%Y"))
			sample.month = int(sample.date_of_receipt_lab.strftime("%m"))
		
		sample.geo_local = Point(lat, lng)
		sample.save()

		## move the files to the right place
		sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), utils.get_path_to_fastq_file(self.request.user.id, sample.id), sample.file_name_1)
		utils.move_file(os.path.join(getattr(settings, "MEDIA_ROOT", None), sample.path_name_1.name), sz_file_to)
		sample.path_name_1.name = os.path.join(utils.get_path_to_fastq_file(self.request.user.id, sample.id), sample.file_name_1)
		
		if (sample.exist_file_2()):
			sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), utils.get_path_to_fastq_file(self.request.user.id, sample.id), sample.file_name_2)
			utils.move_file(os.path.join(getattr(settings, "MEDIA_ROOT", None), sample.path_name_2.name), sz_file_to)
			sample.path_name_2.name = os.path.join(utils.get_path_to_fastq_file(self.request.user.id, sample.id), sample.file_name_2)
		sample.save()

		### create a task to perform the analysis of fastq and trimmomatic
		taskID = async(software.run_fastq_and_trimmomatic_and_identify_species, sample, self.request.user)
		
		### 
		manageDatabase = ManageDatabase()
		manageDatabase.set_metakey(sample, self.request.user, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue, taskID)
		
		messages.success(self.request, "Sample '" + name + "' was created successfully", fail_silently=True)
		return super(SamplesAddView, self).form_valid(form)

	form_valid_message = ""		## need to have this, even empty

class AddValueModal(ModalTemplateView):
	'''
		 This modal inherit of ModalTemplateView, so it just display a text without logic.
	'''
	def __init__(self, *args, **kwargs):
		'''
			You have to call the init method of the parent, before to overide the values:
				- title: The title display in the modal-header
				- icon: The css class that define the modal's icon
				- description: The content of the modal.
				- close_button: A button object that has several attributes.(explain below)
		'''
		super(AddValueModal, self).__init__(*args, **kwargs)
		self.title = "Add value"
		self.description = "This is my description"
		self.icon = "icon-mymodal"


class SamplesDetailView(LoginRequiredMixin, DetailView):
	"""
	Sample detail view
	"""
	model = Sample
	template_name = "samples/sample_detail.html"

	def get_context_data(self, **kwargs):
		context = super(SamplesDetailView, self).get_context_data(**kwargs)
		sample = kwargs['object']
		if (sample.owner.id != self.request.user.id): context['error_cant_see'] = "1"
		context['nav_sample'] = True
		context['virus_identify'] = sample.get_type_sub_type()
		context['href_fastq_1'] = mark_safe('<a href="' + sample.get_fastq(TypePath.MEDIA_URL, True) + '">' + sample.file_name_1 + '</a>')
		context['href_fastq_2'] = mark_safe('<a href="' + sample.get_fastq(TypePath.MEDIA_URL, False) + '">' + sample.file_name_2 + '</a>')
		context['href_fastq_quality_1'] = mark_safe('<a target="_blank" href="' + sample.get_fastq_output(TypePath.MEDIA_URL, True) + '">' + sample.file_name_1 + '.html</a>')
		context['href_fastq_quality_2'] = mark_safe('<a target="_blank" href="' + sample.get_fastq_output(TypePath.MEDIA_URL, True) + '">' + sample.file_name_2 + '.html</a>')
		context['href_trimmonatic_1'] = mark_safe('<a href="' + sample.get_trimmomatic_file(TypePath.MEDIA_URL, True) + '">' + os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_URL, True)) + '</a>')
		context['href_trimmonatic_2'] = mark_safe('<a href="' + sample.get_trimmomatic_file(TypePath.MEDIA_URL, False) + '">' + os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_URL, False)) + '</a>')
		context['href_trimmonatic_quality_1'] = mark_safe('<a target="_blank" href="' + sample.get_fastq_trimmomatic(TypePath.MEDIA_URL, True) + '">' + os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_URL, True)) + '.html</a>')
		context['href_trimmonatic_quality_2'] = mark_safe('<a target="_blank" href="' + sample.get_fastq_trimmomatic(TypePath.MEDIA_URL, False) + '">' + os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_URL, False)) + '.html</a>')
	
		### software
		manageDatabase = ManageDatabase()
		meta_sample = manageDatabase.get_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software, MetaKeyAndValue.META_VALUE_Success)
		if (meta_sample == None):
			meta_sample = manageDatabase.get_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, MetaKeyAndValue.META_VALUE_Success)
			if (meta_sample != None):
				lst_data = meta_sample.description.split(',')
				context['fastq_software'] = lst_data[1].strip()
				context['trimmomatic_software'] = lst_data[2].strip()
		else:
			decodeResult = DecodeResult()
			result = decodeResult.decode_result(meta_sample.description)
			context['fastq_software'] = result.get_software(SoftwareNames.SOFTWARE_FASTQ_name)
			context['trimmomatic_software'] = result.get_software(SoftwareNames.SOFTWARE_TRIMMOMATIC_name)

		meta_sample = manageDatabase.get_metakey(sample, MetaKeyAndValue.META_KEY_Identify_Sample_Software, MetaKeyAndValue.META_VALUE_Success)
		if (meta_sample == None):
			meta_sample = manageDatabase.get_metakey(sample, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Success)
			if (meta_sample != None):
				lst_data = meta_sample.description.split(',')
				context['spades_software'] = lst_data[1].strip()
				context['abricate_software'] = lst_data[2].strip()
		else:
			decodeResult = DecodeResult()
			result = decodeResult.decode_result(meta_sample.description)
			context['spades_software'] = result.get_software(SoftwareNames.SOFTWARE_SPAdes_name)
			context['abricate_software'] = result.get_software(SoftwareNames.SOFTWARE_ABRICATE_name)
		
		return context


class ProjectsView(LoginRequiredMixin, ListView):
	model = Project
	template_name = 'project/projects.html'
	context_object_name = 'projects'
##	group_required = u'company-user' security related with GroupRequiredMixin
	
	def get_context_data(self, **kwargs):
		context = super(ProjectsView, self).get_context_data(**kwargs)
		query_set = Project.objects.filter(owner__id=self.request.user.id, is_deleted=False).order_by('-creation_date')
		table = ProjectTable(query_set)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
		### clean check box in the session
		clean_check_box_in_session(self.request)

		context['table'] = table
		context['nav_project'] = True
		context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
		return context


class ReferenceProjectList(ListView, LoginRequiredMixin):
	"""
	
	"""
	model = Reference
	context_object_name = 'reference'
	ordering = ['id']
##	group_required = u'company-user' security related with GroupRequiredMixin
	
	def get_context_data(self, **kwargs):
		context = super(ReferenceProjectList, self).get_context_data(**kwargs)
		query_set = Reference.objects.filter(owner__id=self.request.user.id, is_obsolete=False).order_by('-name')
		table = ReferenceTable(query_set)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
		
		context['table'] = table
		context['nav_reference'] = True
		context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
		return context

class ProjectCreateView(LoginRequiredMixin, FormValidMessageMixin, generic.CreateView):
	"""
	Create a new reference
	"""
	model = Project
	fields = ['name']
	success_url = reverse_lazy('projects')
	template_name = 'project/project_add.html'


	def get_form_kwargs(self):
		"""
		Set the request to pass in the form
		"""
		kw = super(ProjectCreateView, self).get_form_kwargs()
##		kw['request'] = self.request 	# get error
		return kw


	def get_context_data(self, **kwargs):
		context = super(ProjectCreateView, self).get_context_data(**kwargs)
		query_set = Reference.objects.filter(owner__id=self.request.user.id, is_obsolete=False).order_by('-name')
		table = ReferenceProjectTable(query_set)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER_SMALL}).configure(table)
		context['table'] = table
		context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER_SMALL
			
		context['nav_project'] = True
		context['nav_modal'] = True	## short the size of modal window
		if self.request.POST: 
			context['project_references'] = ReferenceProjectFormSet(self.request.POST)
		else: 
			context['project_references'] = ReferenceProjectFormSet()
		return context


	def form_valid(self, form):
		"""
		Validate the form
		"""
		context = self.get_context_data()
		project_references = context['project_references']
		name = form.cleaned_data['name']
		try:
			Project.objects.get(name=name, owner__username=self.request.user.username)
			return super(ProjectCreateView, self).form_invalid(form)
		except Project.DoesNotExist:
			pass
		if ('select_ref' in project_references.data ): select_ref = project_references.data['select_ref']
		else: return super(ProjectCreateView, self).form_invalid(form)
		if (name is None): name = ""
		
		try:
			reference = Reference.objects.get(pk=select_ref)
		except Sample.DoesNotExist:
			return super(ProjectCreateView, self).form_invalid(form)
		
		with transaction.atomic():
			project = Project()
			project.name = name
			project.reference = reference
			project.owner = self.request.user
			project.is_deleted = False
			project.save()

		messages.success(self.request, "Project '" + name + "' was created successfully", fail_silently=True)
		return super(ProjectCreateView, self).form_valid(form)
	form_valid_message = ""		## need to have this, even empty


class AddSamplesProjectsView(LoginRequiredMixin, FormValidMessageMixin, generic.CreateView):
	"""
	Create a new reference
	"""
	utils = Utils()
	model = Sample
	fields = ['name']
	success_url = reverse_lazy('projects')
	template_name = 'project_sample/project_sample_add.html'
	
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	def get_form_kwargs(self):
		"""
		Set the request to pass in the form
		"""
		kw = super(AddSamplesProjectsView, self).get_form_kwargs()
##		kw['request'] = self.request 	# get error
		return kw

	def get_context_data(self, **kwargs):
		context = super(AddSamplesProjectsView, self).get_context_data(**kwargs)
		
		### test if the user is the same of the page
		project = Project.objects.get(pk=self.kwargs['pk'])
		if (project.owner.id != self.request.user.id): context['error_cant_see'] = "1"

		tag_search = 'search_add_project_sample'
		## catch everything that is not in connection with project 
		query_set = Sample.objects.filter(Q(owner__id=self.request.user.id) & Q(is_obsolete=False) & Q(is_rejected=False) & Q(is_ready_for_projects=True) &\
									(Q(project_sample__isnull=True) | ~Q(project_sample__project__id=project.id) |\
									(Q(project_sample__project__id=project.id) & Q(project_sample__is_deleted=True))) ).distinct()
		table = SampleToProjectsTable(query_set)

		### set the check_box
		if (Constants.CHECK_BOX_ALL not in self.request.session):
			self.request.session[Constants.CHECK_BOX_ALL] = False
			is_all_check_box_in_session(["{}_{}".format(Constants.CHECK_BOX, key.id) for key in query_set], self.request)
		elif ("search_in_table" not in self.request.GET and not is_all_check_box_in_session(["{}_{}".format(Constants.CHECK_BOX, key.id) for key in query_set], self.request)):
			self.request.session[Constants.CHECK_BOX_ALL] = False

		context[Constants.CHECK_BOX_ALL] = self.request.session[Constants.CHECK_BOX_ALL]
		## need to clean all the others if are reject in filter
		dt_sample_id_add_temp = {}
		if (context[Constants.CHECK_BOX_ALL]):
			for sample in query_set: dt_sample_id_add_temp[sample.id] = 1	## add the ids that are in the tables
			for key in self.request.session.keys():
				if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and self.utils.is_integer(key.split('_')[2])):
					### this is necessary because of the search. Can occur some checked box that are out of filter.
					if (int(key.split('_')[2]) not in dt_sample_id_add_temp):
						self.request.session[key] = False
					else: self.request.session[key] = True
		## END need to clean all the others if are reject in filter
			
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER_SMALL}).configure(table)
		if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get('search_add_project_sample')
		context['table'] = table
		context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER_SMALL
		context['project_name'] = project.name
		context['nav_project'] = True
		context['nav_modal'] = True	## short the size of modal window
		if self.request.POST: 
			context['project_sample'] = AddSampleProjectForm(self.request.POST)
		else: 
			context['project_sample'] = AddSampleProjectForm()
		return context


	def form_valid(self, form):
		"""
		Validate the form
		"""
		metaKeyAndValue = MetaKeyAndValue()
		manageDatabase = ManageDatabase()
		software = Software()
		collect_extra_data = CollectExtraData()
		
		### get project sample..
		context = self.get_context_data()

		### get project 		
		project = Project.objects.get(pk=self.kwargs['pk'])
		
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
		project_sample_add = 0
		vect_task_id_submited = []
		for id_sample in vect_sample_id_add:
			try:
				sample = Sample.objects.get(pk=id_sample)
			except Sample.DoesNotExist:
				## log
				self.logger_production.error('Fail to get sample_id {} in ProjectSample'.format(key.split('_')[2]))
				self.logger_debug.error('Fail to get sample_id {} in ProjectSample'.format(key.split('_')[2]))
				continue
			
			## get project sample
			try:
				project_sample = ProjectSample.objects.get(project__id=project.id, sample__id=sample.id)
				
				### if exist can be deleted, active
				if (project_sample.is_deleted and not project_sample.is_error):
					project_sample.is_deleted = False
					project_sample.save()
					project_sample_add += 1
				
			except ProjectSample.DoesNotExist:
				project_sample = ProjectSample()
				project_sample.project = project
				project_sample.sample = sample
				project_sample.save()
				project_sample_add += 1
				
				### create a task to perform the analysis of fastq and trimmomatic
				taskID = async(software.process_second_stage_snippy_coverage_freebayes, project_sample, self.request.user)
				vect_task_id_submited.append(taskID)

				### set project sample queue ID
				manageDatabase.set_project_sample_metakey(project_sample, self.request.user,\
									metaKeyAndValue.get_meta_key_queue_by_project_sample_id(project_sample.id),\
									MetaKeyAndValue.META_VALUE_Queue, taskID)

		### necessary to calculate the global results again 
		if (project_sample_add > 0):
			taskID = async(collect_extra_data.collect_extra_data_for_project, project, self.request.user, vect_task_id_submited)
			manageDatabase.set_project_metakey(project, self.request.user, metaKeyAndValue.get_meta_key_by_project_id(\
					MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id), MetaKeyAndValue.META_VALUE_Queue, taskID)
		
		if (len(vect_task_id_submited) == 0):
			messages.warning(self.request, _("None sample was added to the project '{}'".format(project.name)))
		else:
			messages.success(self.request, _("'{}' {} added to your project '{}'".format(\
				len(vect_task_id_submited), "samples were" if len(vect_task_id_submited) > 1 else "sample is", project.name)), fail_silently=True)
			
		return super(AddSamplesProjectsView, self).form_valid(form)
	form_valid_message = ""		## need to have this, even empty


class ShowSampleProjectsView(LoginRequiredMixin, ListView):
	model = Project
	template_name = 'project_sample/project_sample_show.html'
	context_object_name = 'project_sample'
	
	def get_context_data(self, **kwargs):
		context = super(ShowSampleProjectsView, self).get_context_data(**kwargs)
		project = Project.objects.get(pk=self.kwargs['pk'])
		query_set = ProjectSample.objects.filter(project__id=project.id, is_finished=True, is_deleted=False).order_by('-creation_date')
		
		tag_search = 'search_add_project_sample'
		### filter the search
		if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)): 
			query_set = query_set.filter(Q(sample__name__icontains=self.request.GET.get(tag_search)) |\
										Q(sample__data_set__name__icontains=self.request.GET.get(tag_search)) |\
										Q(sample__type_subtype__icontains=self.request.GET.get(tag_search)))
		table = ShowProjectSamplesResults(query_set)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)

		### can't see this project
		if (project.owner.id != self.request.user.id): context['error_cant_see'] = "1"

		if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get('search_add_project_sample')		
		context['table'] = table
		context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
		context['project_name'] = project.name
		context['project_id'] = project.id
		context['spinner_url'] = os.path.join("/" + Constants.DIR_STATIC, Constants.DIR_ICONS, Constants.AJAX_LOADING_GIF)
		context['nav_project'] = True
		context['nav_modal'] = True	## short the size of modal window
		
		## tODO , set this in database 
		utils = Utils()
		dict_genes = utils.get_elements_and_genes(project.reference.get_reference_gbk(TypePath.MEDIA_ROOT))
		if (dict_genes != None): context['elements'] = sorted(dict_genes.keys())
		
		return context
	
######################################
###
###		AJAX methods for check box in session
###
def is_all_check_box_in_session(vect_check_to_test, request):
	"""
	test if all check boxes are in session 
	If not remove the ones that are in and create the new ones all False
	"""
	utils = Utils()
	dt_data = {}
	
	## get the dictonary
	for key in request.session.keys():
		if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and utils.is_integer(key.split('_')[2])):
			dt_data[key] = True
	
	b_different = False
	if (len(vect_check_to_test) != len(dt_data)): 
		b_different = True
	
	## test the vector
	if (not b_different):
		for key in vect_check_to_test:
			if (key not in dt_data):
				b_different = True
				break
		
	if (b_different):
		## remove all
		for key in dt_data:
			del request.session[key] 
	
		## create new
		for key in vect_check_to_test:
			request.session[key] = False
		return False
	return True
	
	
def clean_check_box_in_session(request):
	"""
	check all check boxes
	"""
	utils = Utils()
	request.session[Constants.CHECK_BOX_ALL] = False
	## clean all check unique
	vect_keys_to_remove = [Constants.CHECK_BOX_ALL]
	for key in request.session.keys():
		if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and utils.is_integer(key.split('_')[2])):
			vect_keys_to_remove.append(key)
	for key in vect_keys_to_remove:
		del request.session[key] 

def set_check_box_values(request):
	"""
	manage check boxes through ajax
	"""
	data = { 'is_ok' : False }
	utils = Utils()
	if (Constants.CHECK_BOX_ALL in request.GET):
		request.session[Constants.CHECK_BOX_ALL] = utils.str2bool(request.GET.get(Constants.CHECK_BOX_ALL))
		## check all unique
		for key in request.session.keys():
			if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and utils.is_integer(key.split('_')[2])):
				request.session[key] = request.session[Constants.CHECK_BOX_ALL]
		data = {
			'is_ok': True
		}
	### one check box is pressed
	elif (Constants.CHECK_BOX in request.GET):
		request.session[Constants.CHECK_BOX_ALL] = False	### set All false
		key_name = "{}_{}".format(Constants.CHECK_BOX, request.GET.get(Constants.CHECK_BOX_VALUE))
		request.session[key_name] = utils.str2bool(request.GET.get(Constants.CHECK_BOX))
		total_checked = 0
		for key in request.session.keys():
			if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and utils.is_integer(key.split('_')[2])):
				if (request.session[key]): total_checked += 1
		data = {
			'is_ok': True,
			'total_checked' : total_checked,
		}
	## get the status of a check_box_single
	elif (Constants.GET_CHECK_BOX_SINGLE in request.GET):
		data['is_ok'] = True
		for key in request.session.keys():
			if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and utils.is_integer(key.split('_')[2])):
				data[key] = request.session[key]
	elif (Constants.COUNT_CHECK_BOX in request.GET):
		count = 0
		for key in request.session.keys():
			if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and utils.is_integer(key.split('_')[2])):
				if (request.session[key]): count += 1
			elif (key == Constants.CHECK_BOX_ALL and request.session[Constants.CHECK_BOX_ALL]): count += 1
		data = {
			'is_ok': True,
			Constants.COUNT_CHECK_BOX : count,
		}
	return JsonResponse(data)
###
###		END AJAX methods for check box in session
###
######################################

	
def show_phylo_canvas(request):
	"""
	manage check boxes through ajax
	"""
	data = { 'is_ok' : False }
	utils = Utils()
	key_with_project_id = 'project_id'
	if (key_with_project_id in request.GET):
		project_id = int(request.GET.get(key_with_project_id))
		element_name = 'all_together'
		key_element_name = 'key_element_name'
		if (key_element_name in request.GET): element_name = int(request.GET.get(key_element_name))
		try:
			project = Project.objects.get(id=project_id)
			if (element_name == 'all_together'): file_name = project.get_global_file_by_project(TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_FASTTREE_tree)
			else: file_name = project.get_global_file_by_element(TypePath.MEDIA_ROOT, element_name, Project.PROJECT_FILE_NAME_FASTTREE_tree)
			if (os.path.exists(file_name)):
				string_file_content = utils.read_file_to_string(file_name).strip()
				if (string_file_content != None and len(string_file_content) > 0):
					data['is_ok'] = True
					data['tree'] = string_file_content
		except Project.DoesNotExist:
			pass
	return JsonResponse(data)


def get_image_coverage(request):
	"""
	get image coverage
	"""
	data = { 'is_ok' : False }
	key_with_project_sample_id = 'project_sample_id'
	key_element = 'element'
	if (key_with_project_sample_id in request.GET and key_element in request.GET):
		try:
			project_sample = ProjectSample.objects.get(id=request.GET.get(key_with_project_sample_id))
			path_name = project_sample.get_global_file_by_element(TypePath.MEDIA_URL, ProjectSample.PREFIX_FILE_COVERAGE, request.GET.get(key_element), FileExtensions.FILE_PNG)
			data['is_ok'] = True
			data['image'] = mark_safe('<img src="{}" style="width: 100%;">'.format(path_name))
			data['text'] = _("Coverage for element '{}'".format(request.GET.get(key_element)))
		except ProjectSample.DoesNotExist as e:
			pass
	return JsonResponse(data)


#####
#####
#####		VALIDATION METHODS WITH AJAX
#####
#####
def validate_project_reference_name(request):
	"""
	test if a name is the same
	"""
	project_name = request.GET.get('project_name')
	user_name = request.GET.get('user_name').strip()
	if (user_name.endswith(" - Logout")): user_name = user_name.replace(" - Logout", "").strip()
	data = {
		'is_taken': Project.objects.filter(name=project_name, owner__username=user_name).exists()
	}
	if data['is_taken']: data['error_message'] = _('Exists a project with this name.')
	return JsonResponse(data)




