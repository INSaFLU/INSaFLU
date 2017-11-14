# Create your views here.

from django.views import generic
from braces.views import LoginRequiredMixin, FormValidMessageMixin
from django.core.urlresolvers import reverse_lazy
from django.views.generic import ListView, DetailView
from django_tables2 import RequestConfig
from managing_files.models import Reference, Sample, Project
from managing_files.tables import ReferenceTable, SampleTable, ProjectTable, ReferenceProjectTable
from managing_files.forms import ReferenceForm, SampleForm, ReferenceProjectFormSet
from managing_files.manage_database import ManageDatabase
from utils.constants import Constants, TypePath
from utils.meta_key_and_values import MetaKeyAndValue
from utils.software import Software
from utils.utils import Utils
from utils.result import DecodeResult
import hashlib, ntpath, os
from django.contrib import messages
from django.conf import settings
from django.contrib.gis.geos import Point
from django_modalview.generic.base import ModalTemplateView
from django_q.tasks import async
from django.utils.safestring import mark_safe
from django.utils.translation import ugettext_lazy as _
from django.http import JsonResponse
from django.db import transaction

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
		software.createFaiToFastaFile(reference.reference_fasta.name)
		
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
		manageDatabase.set_metakey(sample, self.request.user, MetaKeyAndValue.META_KEY_Import_Sample_Import_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue, taskID)
		
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
			context['fastq_software'] = result.get_software(Software.SOFTWARE_FASTQ_name)
			context['trimmomatic_software'] = result.get_software(Software.SOFTWARE_TRIMMOMATIC_name)

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
			context['spades_software'] = result.get_software(Software.SOFTWARE_SPAdes_name)
			context['abricate_software'] = result.get_software(Software.SOFTWARE_ABRICATE_name)
		
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

		
#####
#####
#####		VALIDATION METHODS WITH AJAX
#####
#####
def validate_project_reference_name(request):
	project_name = request.GET.get('project_name')
	user_name = request.GET.get('user_name').strip()
	if (user_name.endswith(" - Logout")): user_name = user_name.replace(" - Logout", "").strip()
	data = {
		'is_taken': Project.objects.filter(name=project_name, owner__username=user_name).exists()
	}
	if data['is_taken']: data['error_message'] = _('Exists a project with this name.')
	return JsonResponse(data)



