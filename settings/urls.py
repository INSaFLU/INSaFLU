"""
Created on Jan 7, 2018

@author: mmp
"""

from django.conf.urls import url

from settings import ajax_views, views

urlpatterns = [
    path(
        "pathogenID-pipeline/<int:level>",
        views.PISettingsGroupsView.as_view(),
        name="pathogenID_pipeline",
    ),
    path(
        "pathogenID-pipeline-sample/<int:sample>",
        views.PIMetagenSampleView.as_view(),
        name="pathogenID_sample_settings",
    ),
    path(
        "pathogenID-pipeline-reset",
        ajax_views.reset_project_settings,
        name="reset_project_parameters",
    ),
    path(
        "software_dataset_update/soft_<int:pk>/dataset_<int:pk_dataset>",
        views.UpdateParametersDatasetView.as_view(),
        name="software-dataset-update",
    ),
    path(
        "software_televir_project_update/soft_<int:pk>/televir_project_<int:pk_televir_project>",
        views.UpdateParametersTelevirProjView.as_view(),
        name="software-televir-project-update",
    ),
    path("", views.index.as_view(), name="settings-index"),
    re_path(r"set_quality", views.QCSettingsView.as_view(), name="settings_qc"),
    path("403/", views.Maintenance.as_view(), name="under_construction"),
    re_path(r"settings", views.SettingsView.as_view(), name="settings"),
    path(
        "software_update/<int:pk>",
        views.UpdateParametersView.as_view(),
        name="software-update",
    ),
    path(
        "software_project_update/soft_<int:pk>/proj_<int:pk_proj>",
        views.UpdateParametersProjView.as_view(),
        name="software-project-update",
    ),
    path(
        "software_project_sample_update/soft_<int:pk>/proj_sample_<int:pk_proj_sample>",
        views.UpdateParametersProjSampleView.as_view(),
        name="software-project-sample-update",
    ),
    path(
        "software_sample_update/soft_<int:pk>/sample_<int:pk_sample>",
        views.UpdateParametersSampleView.as_view(),
        name="software-sample-update",
    ),
    path(
        "ajax/default_parameters",
        ajax_views.set_default_parameters,
        name="default_parameters",
    ),
    url(
        r"^ajax/get_mdcg_project_software$",
        ajax_views.get_mdcg_project_software,
        name="get_mdcg_project_software",
    ),
    url(
        r"^ajax/turn_on_off_software$",
        ajax_views.turn_on_off_software,
        name="turn_on_off_software",
    ),
    path(
        "ajax/mask_consensus", ajax_views.mask_consensus, name="mask_consensus"
    ),  ## set positions to mask consensus
    path(
        "ajax/mask_consensus_actual_values",
        ajax_views.get_mask_consensus_actual_values,
        name="mask_consensus_actual_values",
    ),  ## get positions to mask consensus
    path(
        "ajax/get_software_name_to_turn_on_off",
        ajax_views.get_software_name_to_turn_on_off,
        name="get_software_name_to_turn_on_off",
    ),  ## message to toggle on/off software
]
