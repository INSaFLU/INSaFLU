'''
Created on Jan 7, 2018

@author: mmp
'''
from django.conf.urls import url
from phylogeny import views

urlpatterns = [
	url(r'phylogeny', views.PhylogenyView.as_view(), name='phylogeny'),
]
