'''
Created on Jan 7, 2018

@author: mmp
'''
from django.conf.urls import url
from settings import views
from settings import ajax_views

urlpatterns = [
	url(r'settings', views.SettingsView.as_view(), name='settings'),
	url(r'software-update/(?P<pk>\d+)', views.UpdateParametersView.as_view(), name='software-update'),
	
	url(r'^ajax/default_parameters', ajax_views.set_default_parameters, name='default_parameters'),			## remove a project
]
