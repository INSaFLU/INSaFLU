# Create your views here.

from django.views import generic
from braces.views import LoginRequiredMixin, FormValidMessageMixin
from django.core.urlresolvers import reverse_lazy
from django.views.generic import ListView
from django_tables2 import RequestConfig
from .models import Reference, Sample
from .tables import ReferenceTable, SampleTable
from .forms import ReferenceForm
from constants.Constants import Constants
from constants.Software import Software
import hashlib, ntpath, os
from django.contrib import messages
from django.conf import settings 

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
		table = ReferenceTable(Reference.objects.filter(owner__id=self.request.user.id).order_by('-name'))
		RequestConfig(self.request, paginate={'per_page': 15}).configure(table)
		context['table'] = table
		context['show_paginatior'] = True
		context['nav_reference'] = True
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
		constants = Constants()
		
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
		sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), constants.get_path_to_reference_file(self.request.user.id, reference.id), reference.reference_fasta_name)
		constants.move_file(os.path.join(getattr(settings, "MEDIA_ROOT", None), reference.reference_fasta.name), sz_file_to)
		reference.reference_fasta.name = os.path.join(constants.get_path_to_reference_file(self.request.user.id, reference.id), reference.reference_fasta_name)
		
		sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), constants.get_path_to_reference_file(self.request.user.id, reference.id), reference.reference_genbank_name)
		constants.move_file(os.path.join(getattr(settings, "MEDIA_ROOT", None), reference.reference_genbank.name), sz_file_to)
		reference.reference_genbank.name = os.path.join(constants.get_path_to_reference_file(self.request.user.id, reference.id), reference.reference_genbank_name)
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
		table = SampleTable(Sample.objects.filter(owner__id=self.request.user.id).order_by('id'))
		RequestConfig(self.request, paginate={'per_page': 15}).configure(table)
		context['table'] = table
		context['nav_sample'] = True
		context['show_paginatior'] = True
		return context

class SamplesAddView(LoginRequiredMixin, FormValidMessageMixin, generic.FormView):
	"""
	Create a new reference
	"""
	form_class = ReferenceForm
	success_url = reverse_lazy('samples')
	template_name = 'samples/sample_add.html'

	def get_context_data(self, **kwargs):
		context = super(SamplesAddView, self).get_context_data(**kwargs)
		context['nav_sample'] = True
		context['nav_modal'] = True	## short the size of modal window
		return context
