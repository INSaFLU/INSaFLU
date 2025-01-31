from __future__ import absolute_import

from braces.views import LoginRequiredMixin
from django.shortcuts import render
from django.urls import reverse
from django.utils.functional import cached_property
from django.views import generic
from view_breadcrumbs import BaseBreadcrumbMixin

# Create your views here.


class PhylogenyView(BaseBreadcrumbMixin, LoginRequiredMixin, generic.TemplateView):
    """
    Home page
    """

    template_name = "phylogeny/phylogeny.html"

    add_home = True

    @cached_property
    def crumbs(self):
        return [("Phylogeny", reverse("phylogeny"))]

    def get_context_data(self, **kwargs):
        context = super(PhylogenyView, self).get_context_data(**kwargs)
        context["nav_phylogeny"] = True
        return context
