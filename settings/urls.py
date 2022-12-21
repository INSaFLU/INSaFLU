"""
Created on Jan 7, 2018

@author: mmp
"""
from django.conf.urls import url

from settings import ajax_views, views

urlpatterns = [
    url(
        r"pathogenID-pipeline/(?P<level>\d+)$",
        views.PISettingsView.as_view(),
        name="pathogenID_pipeline",
    ),
    url(
        r"pathogenID-pipeline-reset$",
        ajax_views.reset_project_settings,
        name="reset_project_parameters",
    ),
    url(
        r"software_dataset_update/soft_(?P<pk>\d+)/dataset_(?P<pk_dataset>\d+)$",
        views.UpdateParametersDatasetView.as_view(),
        name="software-dataset-update",
    ),
    url(
        r"software_televir_project_update/soft_(?P<pk>\d+)/televir_project_(?P<pk_televir_project>\d+)$",
        views.UpdateParametersTelevirProjView.as_view(),
        name="software-televir-project-update",
    ),
    url(r"^$", views.index.as_view(), name="settings-index"),
    url(r"set_quality", views.QCSettingsView.as_view(), name="settings_qc"),
    url(r"403/$", views.Maintenance.as_view(), name="under_construction"),
    url(r"settings", views.SettingsView.as_view(), name="settings"),
    url(
        r"software_update/(?P<pk>\d+)$",
        views.UpdateParametersView.as_view(),
        name="software-update",
    ),
    url(
        r"software_project_update/soft_(?P<pk>\d+)/proj_(?P<pk_proj>\d+)$",
        views.UpdateParametersProjView.as_view(),
        name="software-project-update",
    ),
    url(
        r"software_project_sample_update/soft_(?P<pk>\d+)/proj_sample_(?P<pk_proj_sample>\d+)$",
        views.UpdateParametersProjSampleView.as_view(),
        name="software-project-sample-update",
    ),
    url(
        r"software_sample_update/soft_(?P<pk>\d+)/sample_(?P<pk_sample>\d+)$",
        views.UpdateParametersSampleView.as_view(),
        name="software-sample-update",
    ),
    url(
        r"^ajax/default_parameters$",
        ajax_views.set_default_parameters,
        name="default_parameters",
    ),
    url(
        r"^ajax/turn_on_off_software$",
        ajax_views.turn_on_off_software,
        name="turn_on_off_software",
    ),
    url(
        r"^ajax/mask_consensus$", ajax_views.mask_consensus, name="mask_consensus"
    ),  ## set positions to mask consensus
    url(
        r"^ajax/mask_consensus_actual_values$",
        ajax_views.get_mask_consensus_actual_values,
        name="mask_consensus_actual_values",
    ),  ## get positions to mask consensus
    url(
        r"^ajax/get_software_name_to_turn_on_off$",
        ajax_views.get_software_name_to_turn_on_off,
        name="get_software_name_to_turn_on_off",
    ),  ## message to toggle on/off software
]
