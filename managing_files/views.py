# Create your views here.

import hashlib, ntpath, os, logging, sys
from django.views import generic
from braces.views import LoginRequiredMixin, FormValidMessageMixin
from django.core.urlresolvers import reverse_lazy
from django.views.generic import ListView, DetailView
from django_tables2 import RequestConfig
from managing_files.models import Reference, Sample, Project, ProjectSample, UploadFiles, MetaKey
from managing_files.tables import ReferenceTable, SampleTable, ProjectTable, ReferenceProjectTable, SampleToProjectsTable
from managing_files.tables import ShowProjectSamplesResults, AddSamplesFromCvsFileTable, AddSamplesFromFastqFileTable
from managing_files.forms import ReferenceForm, SampleForm, ReferenceProjectFormSet, AddSampleProjectForm, SamplesUploadMultipleFastqForm
from managing_files.forms import SamplesUploadDescriptionForm
from managing_files.manage_database import ManageDatabase
from constants.constants import Constants, TypePath, FileExtensions, TypeFile, FileType
from constants.software_names import SoftwareNames
from constants.meta_key_and_values import MetaKeyAndValue
from utils.software import Software
from utils.collect_extra_data import CollectExtraData
from utils.utils import Utils
from utils.parse_in_files import UploadFilesByDjangoQ, ParseInFiles
from utils.result import DecodeObjects
from utils.process_SGE import ProcessSGE
from django.contrib import messages
from django.conf import settings
from django.contrib.gis.geos import Point
from django_q.tasks import async
from django.utils.safestring import mark_safe
from django.utils.translation import ugettext_lazy as _
from django.core.files.temp import NamedTemporaryFile
from django.db import transaction
from django.db.models import Q
from django.http import JsonResponse
from operator import attrgetter
from itertools import chain
from extend_user.models import Profile
from django.http import HttpResponseRedirect

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
		query_set = Reference.objects.filter(owner__id=self.request.user.id, is_obsolete=False, is_deleted=False).order_by('-name')
		if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)): 
			query_set = query_set.filter(Q(name__icontains=self.request.GET.get(tag_search)) |\
										Q(owner__username__icontains=self.request.GET.get(tag_search)) |\
										Q(reference_genbank_name__icontains=self.request.GET.get(tag_search)) |\
										Q(reference_fasta_name__icontains=self.request.GET.get(tag_search)) |\
										Q(isolate_name__icontains=self.request.GET.get(tag_search)))
		
		### get the references from the system
		query_set_system = Reference.objects.filter(owner__username=Constants.DEFAULT_USER, is_obsolete=False, is_deleted=False).order_by('-name')
		if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)): 
			query_set_system = query_set_system.filter(Q(name__icontains=self.request.GET.get(tag_search)) |\
									Q(owner__username__icontains=self.request.GET.get(tag_search)) |\
									Q(reference_genbank_name__icontains=self.request.GET.get(tag_search)) |\
									Q(reference_fasta_name__icontains=self.request.GET.get(tag_search)) |\
									Q(isolate_name__icontains=self.request.GET.get(tag_search)))
		query_set_result = sorted(chain(query_set, query_set_system), key=attrgetter('creation_date'))
		table = ReferenceTable(query_set_result)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
		if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
		context['table'] = table
		context['nav_reference'] = True
		context['show_paginatior'] = len(query_set_result) > Constants.PAGINATE_NUMBER
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
		context['user_mmp'] = (self.request.user.username == "mmp")	## short the size of modal window
		return context
	
	@transaction.atomic
	def form_valid(self, form):
		
		### test anonymous account
		try:
			profile = Profile.objects.get(user=self.request.user)
			if (profile.only_view_project):
				messages.warning(self.request, "'{}' account can not add references.".format(self.request.user.username), fail_silently=True)
				return super(ReferenceAddView, self).form_invalid(form)
		except Profile.DoesNotExist:
			pass
		
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
		
		reference = form.save(commit=False)
		## set other data
		reference.display_name = reference.name
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
		software.dos_2_unix(sz_file_to)
		reference.reference_fasta.name = os.path.join(utils.get_path_to_reference_file(self.request.user.id, reference.id), reference.reference_fasta_name)
		
		sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), utils.get_path_to_reference_file(self.request.user.id, reference.id), reference.reference_genbank_name)
		if (reference_genbank == None): utils.move_file(temp_genbank_file, sz_file_to)
		else: utils.move_file(os.path.join(getattr(settings, "MEDIA_ROOT", None), reference.reference_genbank.name), sz_file_to)
		software.dos_2_unix(sz_file_to)
		reference.reference_genbank.name = os.path.join(utils.get_path_to_reference_file(self.request.user.id, reference.id), reference.reference_genbank_name)
		reference.save()
		
		### create bed and index for genbank
		utils.from_genbank_to_bed(sz_file_to, reference.get_reference_bed(TypePath.MEDIA_ROOT))
		software.create_index_files_from_igv_tools(reference.get_reference_bed(TypePath.MEDIA_ROOT))
				
		### save in database the elements and coordinates
		utils.get_elements_from_db(reference, self.request.user)
		utils.get_elements_and_cds_from_db(reference, self.request.user)

		## create the index before commit in database, throw exception if something goes wrong
		software.create_fai_fasta(os.path.join(getattr(settings, "MEDIA_ROOT", None), reference.reference_fasta.name))
		
		## remove genbank temp dir if exist 
		if (temp_genbank_dir != None): utils.remove_dir(temp_genbank_dir)
		
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
		query_set = Sample.objects.filter(owner__id=self.request.user.id, is_deleted=False).order_by('-creation_date')
		if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)): 
			query_set = query_set.filter(Q(name__icontains=self.request.GET.get(tag_search)) |\
										Q(type_subtype__icontains=self.request.GET.get(tag_search)) |\
										Q(data_set__name__icontains=self.request.GET.get(tag_search)))
			
		table = SampleTable(query_set)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
		if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
		context['table'] = table
		context['nav_sample'] = True
		context['total_itens'] = query_set.count()
		context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
		return context


class SamplesAddView(LoginRequiredMixin, FormValidMessageMixin, generic.FormView):
	"""
	Create a new reference
	"""
	form_class = SampleForm
	success_url = reverse_lazy('samples')
	template_name = 'samples/sample_add.html'

	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")

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


	@transaction.atomic
	def form_valid(self, form):
		"""
		Validate the form
		"""

		### test anonymous account
		try:
			profile = Profile.objects.get(user=self.request.user)
			if (profile.only_view_project):
				messages.warning(self.request, "'{}' account can not add samples.".format(self.request.user.username), fail_silently=True)
				return super(SamplesAddView, self).form_invalid(form)
		except Profile.DoesNotExist:
			pass

		software = Software()
		utils = Utils()
		name = form.cleaned_data['name']
		lat = form.cleaned_data['lat']
		lng = form.cleaned_data['lng']
		like_dates = form.cleaned_data['like_dates']
			
		sample = form.save(commit=False)
		## set other data
		sample.owner = self.request.user
		sample.is_deleted = False
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
		
		### test geo spacing
		if (lat != None and lng != None): sample.geo_local = Point(lat, lng)
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
		try:
			if (settings.RUN_SGE):
				process_SGE = ProcessSGE()
				if (settings.RUN_SGE_INTO_DJANGOQ):
					taskID = async(process_SGE.set_run_trimmomatic_species, sample, self.request.user)
				else:
					taskID = process_SGE.set_run_trimmomatic_species(sample, self.request.user)
			else:
				taskID = async(software.run_fastq_and_trimmomatic_and_identify_species, sample, self.request.user)
		except Exception as e:
			self.logger_production.error('Fail to run: ProcessSGE - ' + str(e))
			self.logger_debug.error('Fail to run: ProcessSGE - ' + str(e))
			return super(SamplesAddView, self).form_invalid(form)
		### 
		manageDatabase = ManageDatabase()
		manageDatabase.set_sample_metakey(sample, self.request.user, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue, taskID)
		
		messages.success(self.request, "Sample '" + name + "' was created successfully", fail_silently=True)
		return super(SamplesAddView, self).form_valid(form)

	form_valid_message = ""		## need to have this, even empty


class SamplesAddDescriptionFileView(LoginRequiredMixin, FormValidMessageMixin, generic.CreateView):
	"""
	Create a new reference
	"""
	utils = Utils()
	success_url = reverse_lazy('samples')
	template_name = 'samples/sample_description_file.html'
	model = UploadFiles
	fields = ['file_name']
	
	def get_form_kwargs(self):
		"""
		Set the request to pass in the form
		"""
		kw = super(SamplesAddDescriptionFileView, self).get_form_kwargs()
		return kw
	
	def get_context_data(self, **kwargs):
		context = super(SamplesAddDescriptionFileView, self).get_context_data(**kwargs)
		tag_search = 'search_samples'
		query_set = UploadFiles.objects.filter(owner__id=self.request.user.id, is_deleted=False,\
				type_file__name=TypeFile.TYPE_FILE_sample_file).order_by('-creation_date')
		if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)): 
			query_set = query_set.filter(Q(file_name__icontains=self.request.GET.get(tag_search)) |\
										Q(owner__username__icontains=self.request.GET.get(tag_search)))
		table = AddSamplesFromCvsFileTable(query_set)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
		if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
		context['table'] = table
		context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
		context['nav_sample'] = True
		
		### test if exists files to process to match with (csv/tsv) file
		context['does_not_exists_fastq_files_to_process'] = UploadFiles.objects.filter(owner__id=self.request.user.id, is_deleted=False,\
				type_file__name=TypeFile.TYPE_FILE_sample_file).order_by('-creation_date').count() == 0
				
		### test if can add other csv file
		count_not_complete = UploadFiles.objects.filter(owner__id=self.request.user.id, is_deleted=False,\
				type_file__name=TypeFile.TYPE_FILE_sample_file, is_processed=False).count()
		if (count_not_complete > 0): context['can_add_other_file'] = "You can't add other file because there's a file not completed." 
		if (count_not_complete > 0): context['can_add_other_file'] = "You cannot add a new file because you must first upload NGS data regarding another file." 
		return context

	def form_valid(self, form):
		"""
		Validate the form
		"""
		return super(SamplesAddDescriptionFileView, self).form_valid(form)

	form_valid_message = ""		## need to have this, even empty


class SamplesUploadDescriptionFileView(LoginRequiredMixin, FormValidMessageMixin, generic.FormView):
	"""
	Create a new reference
	"""
	form_class = SamplesUploadDescriptionForm
	success_url = reverse_lazy('sample-add-file')
	template_name = 'samples/samples_upload_description_file.html'

	def get_form_kwargs(self):
		"""
		Set the request to pass in the form
		"""
		kw = super(SamplesUploadDescriptionFileView, self).get_form_kwargs()
		kw['request'] = self.request 	# get error
		return kw
	
	def get_context_data(self, **kwargs):
		context = super(SamplesUploadDescriptionFileView, self).get_context_data(**kwargs)
		if ('form' in kwargs and hasattr(kwargs['form'], 'error_in_file')):
			context['error_in_file'] = mark_safe(kwargs['form'].error_in_file.replace('\n', "<br>")) ## pass a list
		context['nav_sample'] = True
		context['nav_modal'] = True	## short the size of modal window
		return context

	def form_valid(self, form):
		
		### test anonymous account
		try:
			profile = Profile.objects.get(user=self.request.user)
			if (profile.only_view_project):
				messages.warning(self.request, "'{}' account can not add file with samples.".format(self.request.user.username), fail_silently=True)
				return super(SamplesUploadDescriptionFileView, self).form_invalid(form)
		except Profile.DoesNotExist:
			pass

		utils = Utils()
		path_name = form.cleaned_data['path_name']

		## create a genbank file
		if (path_name != None):
			upload_files = form.save(commit=False)
			upload_files.is_valid = True
			upload_files.is_processed = False
			upload_files.is_deleted = False
			upload_files.number_errors = 0
			upload_files.number_files_processed = 0
			upload_files.number_files_to_process = form.number_files_to_process
			
			try:
				type_file = MetaKey.objects.get(name=TypeFile.TYPE_FILE_sample_file)
			except MetaKey.DoesNotExist:
				type_file = MetaKey()
				type_file.name = TypeFile.TYPE_FILE_sample_file
				type_file.save()
			
			upload_files.type_file = type_file
			upload_files.file_name = ntpath.basename(path_name.name)
			upload_files.owner = self.request.user
			
			upload_files.description = ""
			upload_files.save()
		
			## move the files to the right place
			sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), utils.get_path_upload_file(self.request.user.id,\
													TypeFile.TYPE_FILE_sample_file), upload_files.file_name)
			sz_file_to = utils.get_unique_file(sz_file_to)		## get unique file name, user can upload files with same name...
			utils.move_file(os.path.join(getattr(settings, "MEDIA_ROOT", None), upload_files.path_name.name), sz_file_to)
			upload_files.path_name.name = os.path.join(utils.get_path_upload_file(self.request.user.id,\
									TypeFile.TYPE_FILE_sample_file),  ntpath.basename(sz_file_to))
			upload_files.save()
			
			### send message to upload samples in the system
			upload_files_by_djangoq = UploadFilesByDjangoQ()
			
			try:
				if (settings.RUN_SGE):
					process_SGE = ProcessSGE()
					if (settings.RUN_SGE_INTO_DJANGOQ):
						taskID = async(process_SGE.set_read_sample_file, upload_files, self.request.user)
					else:
						taskID =  process_SGE.set_read_sample_file(upload_files, self.request.user)
				else:
					b_testing = False
					taskID = async(upload_files_by_djangoq.read_sample_file, self.request.user, upload_files, b_testing)
			except:
				return super(SamplesUploadDescriptionFileView, self).form_invalid(form)
			
			messages.success(self.request, "File '" + upload_files.file_name + "' with samples was uploaded successfully", fail_silently=True)
			return super(SamplesUploadDescriptionFileView, self).form_valid(form)
		return super(SamplesUploadDescriptionFileView, self).form_invalid(form)

	## static method, not need for now.
	form_valid_message = ""		## need to have this

class SamplesAddFastQView(LoginRequiredMixin, FormValidMessageMixin, generic.FormView):
	"""
	Create a new reference
	"""
	form_class = SampleForm
	success_url = reverse_lazy('samples')
	template_name = 'samples/sample_fastq_file.html'


	def get_form_kwargs(self):
		"""
		Set the request to pass in the form
		"""
		kw = super(SamplesAddFastQView, self).get_form_kwargs()
		kw['request'] = self.request 	# the trick!
		return kw


	def get_context_data(self, **kwargs):
		context = super(SamplesAddFastQView, self).get_context_data(**kwargs)
		
		### test anonymous account
		disable_upload_files = False
		try:
			profile = Profile.objects.get(user=self.request.user)
			if (profile.only_view_project): disable_upload_files = True
		except Profile.DoesNotExist:
			pass
		
		### test to show only the processed or all
		if ('show-not-only-checked' in self.request.GET):
			b_show_all = self.request.GET.get('show-not-only-checked') != 'on'
		else: b_show_all = True
		
		### 
		tag_search = 'search_samples'
		query_set = UploadFiles.objects.filter(owner__id=self.request.user.id, is_deleted=False,\
				type_file__name=TypeFile.TYPE_FILE_fastq_gz).order_by('-creation_date')
		if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)):
			query_set = query_set.filter(Q(file_name__icontains=self.request.GET.get(tag_search)) |\
							Q(owner__username__icontains=self.request.GET.get(tag_search)))
		if (not b_show_all):
			query_set = query_set.filter(Q(is_processed=b_show_all))
			
		table = AddSamplesFromFastqFileTable(query_set)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
		if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
		context['table'] = table
		context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
		context['total_itens'] = query_set.count() 
		context['nav_sample'] = True
		context['disable_upload_files'] = disable_upload_files
		context['check_box_not_show_processed_files'] = not b_show_all
		return context


	def form_valid(self, form):
		"""
		Validate the form
		"""
		return super(SamplesAddFastQView, self).form_valid(form)

	form_valid_message = ""		## need to have this, even empty
	
class SamplesUploadFastQView(LoginRequiredMixin, FormValidMessageMixin, generic.FormView):
	"""
	Create a new reference
	"""
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	form_class = SamplesUploadMultipleFastqForm
	success_url = reverse_lazy('sample-add-fastq')
	template_name = 'samples/samples_upload_fastq_files.html'
	utils = Utils()

	def get_form_kwargs(self):
		"""
		Set the request to pass in the form
		"""
		kw = super(SamplesUploadFastQView, self).get_form_kwargs()
		kw['request'] = self.request 	# get error
		return kw
	
	def get_context_data(self, **kwargs):
		context = super(SamplesUploadFastQView, self).get_context_data(**kwargs)
		context['nav_sample'] = True
		context['nav_modal'] = True	## short the size of modal window
		return context

	def post(self, request):
		form = SamplesUploadMultipleFastqForm(request.POST, request.FILES, request=request)
		
		data = {}	## return data
		try:
			if form.is_valid():
				## doesn't work like that
				#upload_files = form.save()
							
				### get the temporary variable
				path_name = form.cleaned_data['path_name']
				
				if (path_name == None):
					data = {'is_valid': False, 'name': self.request.FILES['path_name'].name, 'message' : 'Internal server error, path not found.' }
					return JsonResponse(data)
				
				upload_files = UploadFiles()
				upload_files.file_name = path_name.name
				## move the files to the right place
				sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), self.utils.get_path_upload_file(self.request.user.id,\
														TypeFile.TYPE_FILE_fastq_gz), upload_files.file_name)
				sz_file_to = self.utils.get_unique_file(sz_file_to)		## get unique file name, user can upload files with same name...
				self.utils.copy_file(path_name.file.name, sz_file_to)
				upload_files.path_name.name = os.path.join(self.utils.get_path_upload_file(self.request.user.id,\
										TypeFile.TYPE_FILE_fastq_gz), ntpath.basename(sz_file_to))
				try:
					type_file = MetaKey.objects.get(name=TypeFile.TYPE_FILE_fastq_gz)
				except MetaKey.DoesNotExist:
					type_file = MetaKey()
					type_file.name = TypeFile.TYPE_FILE_fastq_gz
					type_file.save()
	
				upload_files.is_valid = True
				upload_files.is_processed = False			## True when all samples are set
				upload_files.owner = request.user
				upload_files.type_file = type_file
				upload_files.number_files_to_process = 1
				upload_files.number_files_processed = 0
				upload_files.description = ""
				upload_files.save()
	
				data = {'is_valid': True, 'name': upload_files.file_name, 'url': mark_safe(upload_files.get_path_to_file(TypePath.MEDIA_URL)) }
			else:
				data = {'is_valid': False, 'name': self.request.FILES['path_name'].name, 'message' : str(form.errors['path_name'][0]) }
		except:
			self.logger_debug.error(sys.exc_info())
			self.logger_production.error(sys.exc_info())
			data = {'is_valid': False, 'name': self.request.FILES['path_name'].name, 'message' : 'Internal server error, unknown error.' }
			return JsonResponse(data)
	
		## if is last file send a message to link files with sample csv file
		if ('is_valid' in data and data['is_valid']): 
			parse_in_files = ParseInFiles()
			
			try:
				if (settings.RUN_SGE):
					process_SGE = ProcessSGE()
					if (settings.RUN_SGE_INTO_DJANGOQ):
						taskID = async(process_SGE.set_link_files, self.request.user)
					else:
						taskID = process_SGE.set_link_files(self.request.user)
				else:
					b_testing = False
					taskID = async(parse_in_files.link_files, self.request.user, b_testing)
			except:
				data = {'is_valid': False, 'name': self.request.FILES['path_name'].name, 'message' : 'Fail to submit SGE job.' }
				return JsonResponse(data)
		return JsonResponse(data)
	
	form_valid_message = ""		## need to have this, even empty



class SamplesDetailView(LoginRequiredMixin, DetailView):
	"""
	Sample detail view
	"""
	utils = Utils()
	model = Sample
	template_name = "samples/sample_detail.html"

	def get_context_data(self, **kwargs):
		context = super(SamplesDetailView, self).get_context_data(**kwargs)
		sample = kwargs['object']
		if (sample.owner.id != self.request.user.id): context['error_cant_see'] = "1"
		context['nav_sample'] = True
		
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
		## the info can be called from two distinct sources
		## main info
		if (list_meta.count() > 0 and list_meta[0].value == MetaKeyAndValue.META_VALUE_Success):
			context['virus_identify'] = sample.get_type_sub_type()
			context['href_fastq_1'] = mark_safe('<a rel="nofollow" href="' + sample.get_fastq(TypePath.MEDIA_URL, True) + '" download="' + sample.file_name_1 + '">' + sample.file_name_1 + '</a>')
			context['href_fastq_quality_1'] = mark_safe('<a rel="nofollow" target="_blank" href="' + sample.get_fastq_output(TypePath.MEDIA_URL, True) + '">' + sample.file_name_1 + '.html</a>')
			context['href_trimmonatic_1'] = mark_safe('<a rel="nofollow" href="' + sample.get_trimmomatic_file(TypePath.MEDIA_URL, True) + '" download="'\
				+ os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_URL, True)) + '">' + os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_URL, True)) + '</a>')
			context['href_trimmonatic_quality_1'] = mark_safe('<a rel="nofollow" target="_blank" href="' + sample.get_fastq_trimmomatic(TypePath.MEDIA_URL, True) + '">' +\
				os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_URL, True)) + '.html</a>')
			
			trimmomatic_file_name = sample.get_trimmomatic_file(TypePath.MEDIA_URL, False)
			if (trimmomatic_file_name != None):
				context['href_trimmonatic_2'] = mark_safe('<a rel="nofollow" href="' + sample.get_trimmomatic_file(TypePath.MEDIA_URL, False) + '" download="'\
					+ os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_URL, False)) + '">' + os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_URL, False)) + '</a>')
				context['href_fastq_quality_2'] = mark_safe('<a rel="nofollow" target="_blank" href="' + sample.get_fastq_output(TypePath.MEDIA_URL, True) + '"">' + sample.file_name_2 + '.html</a>')
				context['href_fastq_2'] = mark_safe('<a rel="nofollow" href="' + sample.get_fastq(TypePath.MEDIA_URL, False) + '" download="' + sample.file_name_2 + '">' + sample.file_name_2 + '</a>')
				context['href_trimmonatic_quality_2'] = mark_safe('<a rel="nofollow" target="_blank" href="' + sample.get_fastq_trimmomatic(TypePath.MEDIA_URL, False) + '">' +\
					os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_URL, False)) + '.html</a>')
			else:	## there's no second file
				context['href_trimmonatic_2'] = _("Not available")
				context['href_fastq_quality_2'] = _("Not available")
				context['href_fastq_2'] = _("Not available")
				context['href_trimmonatic_quality_2'] = _("Not available")

			### software
			meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software, MetaKeyAndValue.META_VALUE_Success)
			if (meta_sample == None):
				meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, MetaKeyAndValue.META_VALUE_Success)
				if (meta_sample != None):
					lst_data = meta_sample.description.split(',')
					context['fastq_software'] = lst_data[1].strip()
					context['trimmomatic_software'] = lst_data[2].strip()
			else:
				decodeResult = DecodeObjects()
				result = decodeResult.decode_result(meta_sample.description)
				context['fastq_software'] = result.get_software(SoftwareNames.SOFTWARE_FASTQ_name)
				context['trimmomatic_software'] = result.get_software(SoftwareNames.SOFTWARE_TRIMMOMATIC_name)
	
			meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Identify_Sample_Software, MetaKeyAndValue.META_VALUE_Success)
			if (meta_sample == None):
				meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Success)
				if (meta_sample != None):
					lst_data = meta_sample.description.split(',')
					context['spades_software'] = lst_data[1].strip()
					context['abricate_software'] = lst_data[2].strip()
			else:
				decodeResult = DecodeObjects()
				result = decodeResult.decode_result(meta_sample.description)
				context['spades_software'] = result.get_software(SoftwareNames.SOFTWARE_SPAdes_name)
				if len(context['spades_software']) == 0: context['spades_software'] = result.get_software(SoftwareNames.SOFTWARE_ABYSS_name)
				context['abricate_software'] = result.get_software(SoftwareNames.SOFTWARE_ABRICATE_name)
				
			##### extra data sample, columns added by the user
			## [[header1, value1], [header2, value2], [header3, value3], ...]
			### if it's to big expand button is better
			tag_names = sample.get_tag_names()
			context['extra_data_sample_expand'] = (tag_names != None and tag_names.count()  > (Constants.START_EXPAND_SAMPLE_TAG_NAMES_ROWS))
			if (tag_names != None): context['extra_data_sample'] = self.utils.grouped(tag_names, 4)
			
			### get different types of alerts
			metaKeyAndValue = MetaKeyAndValue()
			alert_out = []
			for key in metaKeyAndValue.get_keys_show_alerts_in_sample_details_view():
				meta_data = manageDatabase.get_sample_metakey(sample, key, MetaKeyAndValue.META_VALUE_Success)
				if (meta_data != None): alert_out.append(meta_data.description)
			context['alerts'] = alert_out
			context['has_type_subtype'] = sample.identify_virus.all().count() > 0
			
		elif (sample.candidate_file_name_1 != None and len(sample.candidate_file_name_1) > 0):
			context['candidate_file_name_1'] = sample.candidate_file_name_1
			if (sample.candidate_file_name_2 != None and len(sample.candidate_file_name_2) > 0):
				context['candidate_file_name_2'] = sample.candidate_file_name_2
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

		### clean project name session
		if Constants.PROJECT_NAME_SESSION in self.request.session:
			del self.request.session[Constants.PROJECT_NAME_SESSION]

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
		
		### get project name
		project_name = self.request.session[Constants.PROJECT_NAME_SESSION] if Constants.PROJECT_NAME_SESSION in self.request.session else ''
		tag_search = 'search_references'
		query_set = Reference.objects.filter(owner__id=self.request.user.id, is_obsolete=False, is_deleted=False).order_by('-name')
		if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)): 
			query_set = query_set.filter(Q(name__icontains=self.request.GET.get(tag_search)) |\
							Q(owner__username__icontains=self.request.GET.get(tag_search)) |\
							Q(isolate_name__icontains=self.request.GET.get(tag_search)))
		
		### get the references from the system
		query_set_system = Reference.objects.filter(owner__username=Constants.DEFAULT_USER, is_obsolete=False, is_deleted=False).order_by('-name')
		if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)): 
			query_set_system = query_set_system.filter(Q(name__icontains=self.request.GET.get(tag_search)) |\
							Q(owner__username__icontains=self.request.GET.get(tag_search)) |\
							Q(isolate_name__icontains=self.request.GET.get(tag_search)))
		query_set_result = sorted(chain(query_set, query_set_system), key=attrgetter('creation_date'))
		table = ReferenceProjectTable(query_set_result)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
		
		### set the check_box, check_box_all is only to control if is the first time or not
		if (Constants.CHECK_BOX_ALL not in self.request.session):
			self.request.session[Constants.CHECK_BOX_ALL] = False
			is_all_check_box_in_session(["{}_{}".format(Constants.CHECK_BOX, key.id) for key in query_set_result], self.request)
		elif ("search_references" in self.request.GET):
			# clean check boxes in search
			dt_sample_id_add_temp = {}
			for reference in query_set_result: dt_sample_id_add_temp[reference.id] = 1	## add the ids that are in the tables
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
		
		context['project_name'] = project_name
		context['table'] = table
		context['show_paginatior'] = len(query_set_result) > Constants.PAGINATE_NUMBER
		context['nav_project'] = True
		if self.request.POST: 
			context['project_references'] = ReferenceProjectFormSet(self.request.POST)
		else: 
			context['project_references'] = ReferenceProjectFormSet()
		return context


	def form_valid(self, form):
		"""
		Validate the form
		"""
		### test anonymous account
		try:
			profile = Profile.objects.get(user=self.request.user)
			if (profile.only_view_project):
				messages.warning(self.request, "'{}' account can not create projects.".format(self.request.user.username), fail_silently=True)
				return super(ProjectCreateView, self).form_invalid(form)
		except Profile.DoesNotExist:
			pass

		
		context = self.get_context_data()
		project_references = context['project_references']
		name = form.cleaned_data['name']
		b_error = False
		try:
			Project.objects.get(name__iexact=name, is_deleted=False, owner__username=self.request.user.username)
			self.request.session[Constants.ERROR_PROJECT_NAME] = _("Exists a project with this name.")
			self.request.session[Constants.PROJECT_NAME] = name
			b_error = True
		except Project.DoesNotExist:
			pass
		
		select_ref = None
		if ('select_ref' in project_references.data ): select_ref = project_references.data['select_ref']
		else:
			### test if it's in the session 
			select_ref = get_unique_pk_from_session(self.request)
			if (select_ref == None):
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
			
			if (reference.is_deleted):
				self.request.session[Constants.ERROR_REFERENCE] = "The reference '{}' was removed.".format(reference.name)
				self.request.session[Constants.PROJECT_NAME] = name
				return super(ProjectCreateView, self).form_invalid(form)
		else:
			self.request.session[Constants.ERROR_REFERENCE] = "You need to select a reference."
			self.request.session[Constants.PROJECT_NAME] = name
			return super(ProjectCreateView, self).form_invalid(form)
			
		### exists an error
		if (b_error): return super(ProjectCreateView, self).form_invalid(form)
		
		with transaction.atomic():
			project = form.save()
			project.reference = reference
			project.owner = self.request.user
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
	
	def get_context_data(self, **kwargs):
		context = super(AddSamplesProjectsView, self).get_context_data(**kwargs)
		
		### test if the user is the same of the page
		project = Project.objects.get(pk=self.kwargs['pk'])
		if (project.owner.id != self.request.user.id): context['error_cant_see'] = "1"

		## catch everything that is not in connection with project 
		samples_out = ProjectSample.objects.filter(Q(project=project) & ~Q(is_deleted=True) & ~Q(is_error=True)).values('sample__pk')
		query_set = Sample.objects.filter(owner__id=self.request.user.id, is_obsolete=False, is_deleted=False,\
					is_ready_for_projects=True).exclude(pk__in=samples_out)
					
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
			
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
		if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
		context['table'] = table
		context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
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

		### test anonymous account
		try:
			profile = Profile.objects.get(user=self.request.user)
			if (profile.only_view_project):
				messages.warning(self.request, "'{}' account can not add samples to a project.".format(self.request.user.username), fail_silently=True)
				return super(AddSamplesProjectsView, self).form_invalid(form)
		except Profile.DoesNotExist:
			pass
		
		if form.is_valid():
			metaKeyAndValue = MetaKeyAndValue()
			manageDatabase = ManageDatabase()
			software = Software()
			collect_extra_data = CollectExtraData()
			process_SGE = ProcessSGE()
			
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
			for id_sample in vect_sample_id_add:
				try:
					sample = Sample.objects.get(pk=id_sample)
				except Sample.DoesNotExist:
					## log
					self.logger_production.error('Fail to get sample_id {} in ProjectSample'.format(key.split('_')[2]))
					self.logger_debug.error('Fail to get sample_id {} in ProjectSample'.format(key.split('_')[2]))
					continue
				
				## the sample can be deleted by other session
				if (sample.is_deleted): continue
				
				## get project sample
				try:
					project_sample = ProjectSample.objects.get(project__id=project.id, sample__id=sample.id)
					
					### if exist can be deleted, pass to active
					if (project_sample.is_deleted and not project_sample.is_error and not project_sample.is_deleted_in_file_system):
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
					try:
						if (settings.RUN_SGE):
							if (settings.RUN_SGE_INTO_DJANGOQ):
								taskID = async(process_SGE.set_second_stage_snippy, project_sample, self.request.user)
							else:
								taskID = process_SGE.set_second_stage_snippy(project_sample, self.request.user)
						else:
							taskID = async(software.process_second_stage_snippy_coverage_freebayes, project_sample, self.request.user)
							
						### set project sample queue ID
						manageDatabase.set_project_sample_metakey(project_sample, self.request.user,\
										metaKeyAndValue.get_meta_key_queue_by_project_sample_id(project_sample.id),\
										MetaKeyAndValue.META_VALUE_Queue, taskID)
					except:
						pass
					
	
			### necessary to calculate the global results again 
			if (project_sample_add > 0):
				try:
					if (settings.RUN_SGE):
						if (settings.RUN_SGE_INTO_DJANGOQ):
							taskID = async(process_SGE.set_collect_global_files, project, self.request.user)
						else:
							taskID = process_SGE.set_collect_global_files(project, self.request.user)
					else:
						taskID = async(collect_extra_data.collect_extra_data_for_project, project, self.request.user)
				
					manageDatabase.set_project_metakey(project, self.request.user, metaKeyAndValue.get_meta_key(\
							MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id), MetaKeyAndValue.META_VALUE_Queue, taskID)
				except:
						pass
			
			if (project_sample_add == 0):
				messages.warning(self.request, _("No sample was added to the project '{}'".format(project.name)))
			else:
				if (project_sample_add > 1):
					messages.success(self.request, _("'{}' samples were added to your project.".format(\
								project_sample_add)), fail_silently=True)
				else:
					messages.success(self.request, _("One sample was added to your project."), fail_silently=True)
			return HttpResponseRedirect(reverse_lazy('projects'))
		else:
			return super(AddSamplesProjectsView, self).form_invalid(form)


	form_valid_message = ""		## need to have this, even empty


class ShowSampleProjectsView(LoginRequiredMixin, ListView):
	model = Project
	template_name = 'project_sample/project_sample_show.html'
	context_object_name = 'project_sample'
	
	def get_context_data(self, **kwargs):
		context = super(ShowSampleProjectsView, self).get_context_data(**kwargs)
		project = Project.objects.get(pk=self.kwargs['pk'])
		query_set = ProjectSample.objects.filter(project__id=project.id, is_finished=True, is_deleted=False, is_error=False).order_by('-creation_date')
		
		tag_search = 'search_add_project_sample'
		### filter the search
		if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)): 
			query_set = query_set.filter(Q(sample__name__icontains=self.request.GET.get(tag_search)) |\
										Q(sample__data_set__name__icontains=self.request.GET.get(tag_search)) |\
										Q(mixed_infections__tag__name__icontains=self.request.GET.get(tag_search)) |\
										Q(sample__type_subtype__icontains=self.request.GET.get(tag_search)))
		table = ShowProjectSamplesResults(query_set)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)

		### can't see this project
		if (project.owner.id != self.request.user.id): context['error_cant_see'] = "1"

		if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get('search_add_project_sample')		
		context['table'] = table
		context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
		context['project_id'] = project.id
		context['spinner_url'] = os.path.join("/" + Constants.DIR_STATIC, Constants.DIR_ICONS, Constants.AJAX_LOADING_GIF)
		context['nav_project'] = True
		
		context['project_name'] = project.name
		context['reference_name'] = project.reference.name
		context['number_of_samples'] = ProjectSample.objects.filter(project=project, is_deleted=False, is_error=False, is_finished=True).count()
		context['number_of_alerts'] = ProjectSample.objects.filter(Q(project=project) & Q(is_deleted=False) &\
								Q(is_error=False) & Q(is_finished=True) & (Q(alert_first_level__gte=1) |\
								Q(alert_second_level__gte=1))).count()
		context['samples_in_process'] = ProjectSample.objects.filter(project=project, is_deleted=False, is_error=False, is_finished=False).count()
		context['samples_error'] = ProjectSample.objects.filter(project=project, is_deleted=False, is_error=True, is_finished=False).count()
		
		## Files
		context['coverage_file'] = project.get_global_file_by_project_web(Project.PROJECT_FILE_NAME_COVERAGE)
		context['main_variations_snippy_file'] = project.get_global_file_by_project_web(Project.PROJECT_FILE_NAME_TAB_VARIATIONS_SNIPPY)
		context['freebays_variations_50_file'] = project.get_global_file_by_project_web(Project.PROJECT_FILE_NAME_TAB_VARIATIONS_FREEBAYES)
		context['sample_file_result_csv'] = project.get_global_file_by_project_web(Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV)
		context['sample_file_result_tsv'] = project.get_global_file_by_project_web(Project.PROJECT_FILE_NAME_SAMPLE_RESULT_TSV)
		
		## tODO , set this in database 
		utils = Utils()
		vect_elements = utils.get_elements_from_db(project.reference, self.request.user)
		if (vect_elements != None and len(vect_elements) > 0): 
			context['elements'] = vect_elements
			### get a vect of genes name
			context['genes'] = utils.get_vect_cds_from_element_from_db(vect_elements[0], project.reference, self.request.user)
		return context

class ShowSampleProjectsDetailsView(LoginRequiredMixin, ListView):
	"""
	"""
	utils = Utils()
	model = ProjectSample
	template_name = 'project_sample/project_sample_single_detail.html'
	context_object_name = 'project_sample'
	
	def get_context_data(self, **kwargs):
		context = super(ShowSampleProjectsDetailsView, self).get_context_data(**kwargs)
		try:
			project_sample = ProjectSample.objects.get(pk=self.kwargs['pk'])
			
			## several data to show
			context['project_sample'] = project_sample
			context['num_alerts'] = project_sample.alert_first_level + project_sample.alert_second_level
			context['nav_project'] = True
			
			## collect alerts
			alert_out = []
			manageDatabase = ManageDatabase()
			metaKeyAndValue = MetaKeyAndValue()
			vect_elements = self.utils.get_elements_from_db(project_sample.project.reference, project_sample.project.owner)
			for element_temp in vect_elements:
				meta_key = metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_ALERT_COVERAGE_0, element_temp)
				meta_data = manageDatabase.get_project_sample_metakey(project_sample, meta_key, MetaKeyAndValue.META_VALUE_Success)
				if (meta_data != None):
					alert_out.append(meta_data.description)
				meta_key = metaKeyAndValue.get_meta_key(MetaKeyAndValue.META_KEY_ALERT_COVERAGE_9, element_temp)
				meta_data = manageDatabase.get_project_sample_metakey(project_sample, meta_key, MetaKeyAndValue.META_VALUE_Success)
				if (meta_data != None):
					alert_out.append(meta_data.description)
			
			### get different types of alerts
			for key in metaKeyAndValue.get_keys_show_alerts_in_sample_projects_details_view():
				meta_data = manageDatabase.get_project_sample_metakey(project_sample, key, MetaKeyAndValue.META_VALUE_Success)
				if (meta_data != None): alert_out.append(meta_data.description)
			context['alerts'] = alert_out
			
			##### extra data sample, columns added by the user
			## [[header1, value1], [header2, value2], [header3, value3], ...]
			### if it's to big expand button is better
			tag_names = project_sample.sample.get_tag_names()
			context['extra_data_sample_expand'] = (tag_names != None) ## (tag_names != None and tag_names.count()  > (Constants.START_EXPAND_SAMPLE_TAG_NAMES_ROWS))
			if (tag_names != None): context['extra_data_sample'] = self.utils.grouped(tag_names, 4)
			
			context['consensus_file'] = project_sample.get_consensus_file_web()
			context['snippy_variants_file'] = project_sample.get_file_web(FileType.FILE_TAB ,SoftwareNames.SOFTWARE_SNIPPY_name)
			context['freebayes_variants_file'] = project_sample.get_file_web(FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)
			
		except ProjectSample.DoesNotExist:
			context['error_cant_see'] = 1
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


def get_first_pk_from_session(request):
	"""
	return the first pk selected, if none, return none 
	"""
	utils = Utils()
	for key in request.session.keys():
		if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and utils.is_integer(key.split('_')[2])\
				and request.session[key]):
			return key.split('_')[2]
	return None

def get_unique_pk_from_session(request):
	"""
	return the unique pk selected. If exists more than one return none 
	"""
	utils = Utils()
	return_pk = None
	for key in request.session.keys():
		if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and utils.is_integer(key.split('_')[2])\
				and request.session[key]):
			if (return_pk == None): return_pk = key.split('_')[2]
			else: return None
	return return_pk



