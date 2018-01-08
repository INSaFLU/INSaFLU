from __future__ import absolute_import

# Create your views here.

from django.shortcuts import render
from django.views import generic
from braces.views import LoginRequiredMixin


class PhylogenyView(LoginRequiredMixin, generic.TemplateView):
	"""
	Home page
	"""
	template_name = 'phylogeny/phylogeny.html'
	
	def get_context_data(self, **kwargs):
		context = super(PhylogenyView, self).get_context_data(**kwargs)
		context['nav_phylogeny'] = True			
		return context
