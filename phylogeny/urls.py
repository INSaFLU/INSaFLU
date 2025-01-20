'''
Created on Jan 7, 2018

@author: mmp
'''
from django.urls import re_path
from phylogeny import views

urlpatterns = [
	re_path(r'phylogeny', views.PhylogenyView.as_view(), name='phylogeny'),
]
