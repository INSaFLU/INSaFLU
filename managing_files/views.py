# Create your views here.

from django.views import generic
from django import forms
from braces.views import LoginRequiredMixin, FormValidMessageMixin
from django.core.urlresolvers import reverse_lazy
from django.views.generic import ListView
from django_tables2 import RequestConfig
from .models import Reference
from .tables import ReferenceTable
from .forms import ReferenceForm
from django.core.files.temp import NamedTemporaryFile
from constants.Constants import Constants
import hashlib
from django.contrib import messages

# http://www.craigderington.me/generic-list-view-with-django-tables/
    
#class ReferencesView(LoginRequiredMixin, GroupRequiredMixin, ListView):
class ReferenceView(LoginRequiredMixin, ListView):
    model = Reference
    template_name = 'references/references.html'
    context_object_name = 'reference'
    ordering = ['id']
##    group_required = u'company-user' security related with GroupRequiredMixin
    
    def get_context_data(self, **kwargs):
        context = super(ReferenceView, self).get_context_data(**kwargs)
        context['nav_reference'] = True
        table = ReferenceTable(Reference.objects.filter(owner__id=self.request.user.id).order_by('-name'))
        RequestConfig(self.request, paginate={'per_page': 15}).configure(table)
        context['table'] = table
        return context


class ReferenceAddView(LoginRequiredMixin, FormValidMessageMixin, generic.FormView):

    form_class = ReferenceForm
    success_url = reverse_lazy('reference')
    template_name = 'references/reference_add.html'

    ## Other solution
    ## https://pypi.python.org/pypi?%3aaction=display&name=django-contrib-requestprovider&version=1.0.1
    def get_form_kwargs(self):
        """
        Set the request to pass in the form
        """
        kw = super(ReferenceAddView, self).get_form_kwargs()
        kw['request'] = self.request # the trick!
        return kw
    
    def form_valid(self, form):
        name = form.cleaned_data['name']
        scentific_name = form.cleaned_data['scientific_name']
        reference_fasta = form.cleaned_data['reference_fasta']
        reference_genbank = form.cleaned_data['reference_genbank']
        
        hash_value_fasta = hashlib.md5(form.files.get('reference_fasta').read()).hexdigest()
        hash_value_genbank = hashlib.md5(form.files.get('reference_genbank').read()).hexdigest()
        
        reference = form.save(commit=False)
        reference.owner = self.request.user
        reference.is_obsolete = False
        reference.number_locus = self.request.session[Constants.NUMBER_LOCUS_FASTA_FILE]
        reference.hash_reference_fasta = hash_value_fasta
        reference.hash_reference_genbank = hash_value_genbank
        reference.save()
        messages.success(self.request, "Reference '" + name + "'was created successfully", fail_silently=True)
        return super(ReferenceAddView, self).form_valid(form)
#         user = authenticate(username=username, password=password)
#         if user is not None and user.is_active:
#             login(self.request, user)
#             return super(LoginView, self).form_valid(form)
#         else:
#             return self.form_invalid(form)

    ## static method
    form_valid_message = "You have a new Reference."
    


