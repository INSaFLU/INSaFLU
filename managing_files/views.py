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
from django_q.tasks import async
from django.utils.safestring import mark_safe
from django.utils.translation import ugettext_lazy as _
from django.core.files.temp import NamedTemporaryFile
from django.db import transaction
from django.db.models import Q

# http://www.craigderington.me/generic-list-view-with-django-tables/
	
#class ReferencesView(LoginRequiredMixin, GroupRequiredMixin, ListView):
class ReferenceView(LoginRequiredMixin, ListView):
	model = Reference
	template_name = 'references/references.html'
	context_object_name = 'reference'
	ordering = ['id']
	
	def get_context_data(self, **kwargs):
		context = super(ReferenceView, self).get_context_data(**kwargs)
		
		tag_search = 'search_references'
		query_set = Reference.objects.filter(owner__id=self.request.user.id).order_by('-name')
		if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)): 
			query_set = query_set.filter(Q(name__icontains=self.request.GET.get(tag_search)) |\
										Q(owner__username__icontains=self.request.GET.get(tag_search)) |\
										Q(reference_genbank_name__icontains=self.request.GET.get(tag_search)) |\
										Q(reference_fasta_name__icontains=self.request.GET.get(tag_search)) |\
										Q(isolate_name__icontains=self.request.GET.get(tag_search)))
		
		table = ReferenceTable(query_set)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
		if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
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
	
	@transaction.atomic
	def form_valid(self, form):
		software = Software()
		utils = Utils()
		
		name = form.cleaned_data['name']
		scentific_name = form.cleaned_data['isolate_name']
		reference_fasta = form.cleaned_data['reference_fasta']
		reference_genbank = form.cleaned_data['reference_genbank']
		if (scentific_name is None): scentific_name = ""
		
		## create a genbank file
		temp_genbank_dir = None
		if (reference_genbank == None):
			reference_fasta_temp_file_name = NamedTemporaryFile(prefix='flu_fa_', delete=False)
			reference_fasta.file.seek(0)
			reference_fasta_temp_file_name.write(reference_fasta.read())
			reference_fasta_temp_file_name.flush()
			reference_fasta_temp_file_name.close()
			
			software = Software()
			try:
				temp_genbank_dir = software.run_prokka(reference_fasta_temp_file_name.name, ntpath.basename(reference_fasta.name))
			except Exception as e:
				os.unlink(reference_fasta_temp_file_name.name)
				messages.error(self.request, "Error creating the genbank file", fail_silently=True)
				return super(ReferenceAddView, self).form_invalid(form)
			os.unlink(reference_fasta_temp_file_name.name)
		
			temp_genbank_file = os.path.join(temp_genbank_dir, utils.clean_extension(ntpath.basename(reference_fasta.name)) + FileExtensions.FILE_GBK)
			if (not os.path.exists(temp_genbank_file)):
				utils.remove_dir(temp_genbank_dir)
				messages.error(self.request, "Error creating the genbank file", fail_silently=True)
				return super(ReferenceAddView, self).form_invalid(form)

		hash_value_fasta = hashlib.md5(form.files.get('reference_fasta').read()).hexdigest()
		if (reference_genbank == None): hash_value_genbank = utils.md5sum(temp_genbank_file)
		else: hash_value_genbank = hashlib.md5(form.files.get('reference_genbank').read()).hexdigest()
		
		## remove genbank temp dir if exist 
		if (temp_genbank_dir != None): utils.remove_dir(temp_genbank_dir)
			
		reference = form.save(commit=False)
		## set other data
		reference.owner = self.request.user
		reference.is_obsolete = False
		reference.number_of_locus = self.request.session[Constants.NUMBER_LOCUS_FASTA_FILE]
		reference.hash_reference_fasta = hash_value_fasta
		reference.reference_fasta_name = ntpath.basename(reference_fasta.name)
		reference.scentific_name = scentific_name
		reference.hash_reference_genbank = hash_value_genbank
		if (reference_genbank == None): reference.reference_genbank_name = ntpath.basename(temp_genbank_file)
		else: reference.reference_genbank_name = ntpath.basename(reference_genbank.name)
		reference.save()

		## move the files to the right place
		sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), utils.get_path_to_reference_file(self.request.user.id, reference.id), reference.reference_fasta_name)
		utils.move_file(os.path.join(getattr(settings, "MEDIA_ROOT", None), reference.reference_fasta.name), sz_file_to)
		reference.reference_fasta.name = os.path.join(utils.get_path_to_reference_file(self.request.user.id, reference.id), reference.reference_fasta_name)
		
		sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), utils.get_path_to_reference_file(self.request.user.id, reference.id), reference.reference_genbank_name)
		if (reference_genbank == None): utils.move_file(temp_genbank_file, sz_file_to)
		else: utils.move_file(os.path.join(getattr(settings, "MEDIA_ROOT", None), reference.reference_genbank.name), sz_file_to)
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
		
	def get_context_data(self, **kwargs):
		context = super(SamplesView, self).get_context_data(**kwargs)
		tag_search = 'search_samples'
		query_set = Sample.objects.filter(owner__id=self.request.user.id).order_by('-creation_date')
		if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)): 
			query_set = query_set.filter(Q(name__icontains=self.request.GET.get(tag_search)) |\
										Q(type_subtype__icontains=self.request.GET.get(tag_search)) |\
										Q(data_set__name__icontains=self.request.GET.get(tag_search)))
			
		table = SampleTable(query_set)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
		if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
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
		tag_search = 'search_projects'
		query_set = Project.objects.filter(owner__id=self.request.user.id, is_deleted=False).order_by('-creation_date')
		if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)): 
			query_set = query_set.filter(Q(name__icontains=self.request.GET.get(tag_search)) |\
										Q(reference__name__icontains=self.request.GET.get(tag_search)))
		table = ProjectTable(query_set)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
		if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
		
		### clean check box in the session
		clean_check_box_in_session(self.request) ## for both samples and references

		context['table'] = table
		context['nav_project'] = True
		context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
		return context


class ProjectCreateView(LoginRequiredMixin, FormValidMessageMixin, generic.CreateView):
	"""
	Create a new reference
	"""
	utils = Utils()
	model = Project
	fields = ['name']
	success_url = reverse_lazy('projects')
	template_name = 'project/project_add.html'

	def get_form_kwargs(self):
		"""
		Set the request to pass in the form
		"""
		kw = super(ProjectCreateView, self).get_form_kwargs()
		return kw

	def get_context_data(self, **kwargs):
		context = super(ProjectCreateView, self).get_context_data(**kwargs)
		
		tag_search = 'search_references'
		query_set = Reference.objects.filter(owner__id=self.request.user.id, is_obsolete=False).order_by('-name')
		if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)): 
			query_set = query_set.filter(Q(name__icontains=self.request.GET.get(tag_search)) |\
							Q(owner__username__icontains=self.request.GET.get(tag_search)) |\
							Q(isolate_name__icontains=self.request.GET.get(tag_search)))
		table = ReferenceProjectTable(query_set)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
		
		### set the check_box, check_box_all is only to control if is the first time or not
		if (Constants.CHECK_BOX_ALL not in self.request.session):
			self.request.session[Constants.CHECK_BOX_ALL] = False
			is_all_check_box_in_session(["{}_{}".format(Constants.CHECK_BOX, key.id) for key in query_set], self.request)
		elif ("search_references" in self.request.GET):
			# clean check boxes in search
			dt_sample_id_add_temp = {}
			for reference in query_set: dt_sample_id_add_temp[reference.id] = 1	## add the ids that are in the tables
			for key in self.request.session.keys():
				if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and self.utils.is_integer(key.split('_')[2])):
					### this is necessary because of the search. Can occur some checked box that are out of filter.
					if (int(key.split('_')[2]) not in dt_sample_id_add_temp):
						self.request.session[key] = False
			
		## clean the 
		if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
		if (Constants.PROJECT_NAME in self.request.session): 
			context[Constants.PROJECT_NAME] = self.request.session[Constants.PROJECT_NAME]
			del self.request.session[Constants.PROJECT_NAME]
		if (Constants.ERROR_PROJECT_NAME in self.request.session): 
			context[Constants.ERROR_PROJECT_NAME] = self.request.session[Constants.ERROR_PROJECT_NAME]
			del self.request.session[Constants.ERROR_PROJECT_NAME]
		if (Constants.ERROR_REFERENCE in self.request.session): 
			context[Constants.ERROR_REFERENCE] = self.request.session[Constants.ERROR_REFERENCE]
			del self.request.session[Constants.ERROR_REFERENCE]
		
		context['table'] = table
		context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
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
		b_error = False
		try:
			Project.objects.get(name=name, owner__username=self.request.user.username)
			self.request.session[Constants.ERROR_PROJECT_NAME] = _("Exists a project with this name.")
			self.request.session[Constants.PROJECT_NAME] = name
			b_error = True
		except Project.DoesNotExist:
			pass
		
		select_ref = None
		if ('select_ref' in project_references.data ): select_ref = project_references.data['select_ref']
		else:
			self.request.session[Constants.ERROR_REFERENCE] = "You need to select a reference." 
			b_error = True
		if (name is None): name = ""
		
		if (select_ref != None):
			try:
				reference = Reference.objects.get(pk=select_ref)
			except Reference.DoesNotExist:
				self.request.session[Constants.ERROR_REFERENCE] = "You need to select a reference."
				self.request.session[Constants.PROJECT_NAME] = name
				return super(ProjectCreateView, self).form_invalid(form)
		
		### exists an error
		if (b_error): return super(ProjectCreateView, self).form_invalid(form)
		
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

		## catch everything that is not in connection with project 
		query_set = Sample.objects.filter(Q(owner__id=self.request.user.id) & Q(is_obsolete=False) & Q(is_rejected=False) & Q(is_ready_for_projects=True) &\
									(Q(project_sample__isnull=True) | ~Q(project_sample__project__id=project.id) |\
									(Q(project_sample__project__id=project.id) & Q(project_sample__is_deleted=True))) ).distinct()
		tag_search = 'search_add_project_sample'
		if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)): 
			query_set = query_set.filter(Q(name__icontains=self.request.GET.get(tag_search)) |\
							Q(type_subtype__icontains=self.request.GET.get(tag_search)) |\
							Q(data_set__name__icontains=self.request.GET.get(tag_search)) |\
							Q(week__icontains=self.request.GET.get(tag_search)))
			tag_search
		table = SampleToProjectsTable(query_set)

		### set the check_box
		if (Constants.CHECK_BOX_ALL not in self.request.session):
			self.request.session[Constants.CHECK_BOX_ALL] = False
			is_all_check_box_in_session(["{}_{}".format(Constants.CHECK_BOX, key.id) for key in query_set], self.request)

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
		if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
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
			manageDatabase.set_project_metakey(project, self.request.user, metaKeyAndValue.get_meta_key(\
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
		manage_database = ManageDatabase()
		vect_genes = manage_database.get_elements_and_genes_from_db(project, self.request.user)
		if (vect_genes != None): context['elements'] = vect_genes
		
		return context
	


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
	check all check boxes on samples/references to add samples
	"""
	utils = Utils()
	## clean all check unique
	if (Constants.CHECK_BOX_ALL in request.session): del request.session[Constants.CHECK_BOX_ALL]
	vect_keys_to_remove = []
	for key in request.session.keys():
		if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and utils.is_integer(key.split('_')[2])):
			vect_keys_to_remove.append(key)
	for key in vect_keys_to_remove:
		del request.session[key] 
