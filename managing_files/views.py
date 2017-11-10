# Create your views here.

from django.views import generic
from braces.views import LoginRequiredMixin, FormValidMessageMixin
from django.core.urlresolvers import reverse_lazy
from django.views.generic import ListView
from django_tables2 import RequestConfig
from managing_files.models import Reference, Sample
from managing_files.tables import ReferenceTable, SampleTable
from managing_files.forms import ReferenceForm, SampleForm
from utils.constants import Constants
from utils.software import Software
from utils.utils import Utils
import hashlib, ntpath, os
from django.contrib import messages
from django.conf import settings
from django.contrib.gis.geos import Point
from django_modalview.generic.base import ModalTemplateView
from django_q.tasks import async

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
		context['show_paginatior'] = True
		context['nav_reference'] = query_set.count() > Constants.PAGINATE_NUMBER
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
		scentific_name = form.cleaned_data['scientific_name']
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
		
		messages.success(self.request, "Reference '" + name + "'was created successfully", fail_silently=True)
		return super(ReferenceAddView, self).form_valid(form)

	## static method, not need for now.
	#form_valid_message = "You have created a new Reference."
	form_valid_message = ""
	

class SamplesView(LoginRequiredMixin, ListView):
	model = Reference
	template_name = 'samples/samples.html'
	context_object_name = 'samples'
	ordering = ['id']
##	group_required = u'company-user' security related with GroupRequiredMixin
	
	def get_context_data(self, **kwargs):
		context = super(SamplesView, self).get_context_data(**kwargs)
		query_set = Sample.objects.filter(owner__id=self.request.user.id).order_by('id')
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
		constants = Constants()
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
			sample.week = int(sample.date_of_onset.strftime("%m"))
			sample.year = int(sample.date_of_onset.strftime("%Y"))
		elif (like_dates == 'date_of_collection'):
			sample.day = int(sample.date_of_collection.strftime("%d"))
			sample.week = int(sample.date_of_collection.strftime("%m"))
			sample.year = int(sample.date_of_collection.strftime("%Y"))
		elif (like_dates == 'date_of_receipt_lab'):
			sample.day = int(sample.date_of_receipt_lab.strftime("%d"))
			sample.week = int(sample.date_of_receipt_lab.strftime("%m"))
			sample.year = int(sample.date_of_receipt_lab.strftime("%Y"))
		
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

		### create a task to perform the anlysis of fastq and trimmomatic
		async(software.run_fastq_and_trimmomatic, sample, self.request.user)
		
		### queue the quality check and
		async(software.identify_type_and_sub_type, sample.get_trimmomatic_file_1(), sample.get_trimmomatic_file_2(), sample, self.request.user)
		
		messages.success(self.request, "Sample '" + name + "'was created successfully", fail_silently=True)
		return super(SamplesAddView, self).form_valid(form)



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
