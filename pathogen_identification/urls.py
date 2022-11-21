from django.conf.urls import url
from managing_files import ajax_views, views

import pathogen_identification.ajax_views as PIajax_views
import pathogen_identification.views as PIviews

urlpatterns = [
    url(
        r"projects$",
        PIviews.PathId_ProjectsView.as_view(),
        name="PIprojects_main",
    ),
    url(
        r"Projects/(?P<pk>\d+)$",
        PIviews.MainPage.as_view(),
        name="PIproject_samples",
    ),
    url(
        r"project_add$",
        PIviews.PathID_ProjectCreateView.as_view(),
        name="PIproject-add",
    ),
    url(
        r"Project_samples/(?P<pk>\d+)/add_sample_project$",
        PIviews.AddSamples_PIProjectsView.as_view(),
        name="add-sample-PIproject",
    ),
    url(
        r"Project/(?P<pk>\d+)/show_project_settings$",
        views.ProjectsSettingsView.as_view(),
        name="PIproject-settings",
    ),
    url(
        r"Project_samples/(?P<pk>\d+)$",
        views.SamplesDetailView.as_view(),
        name="remove-sample-PIproject",
    ),
    url(
        r"Projects/(?P<project>[a-zA-Z0-9_]+)/all_reports$",
        PIviews.Project_reports,
        name="all_PIproject_reports",
    ),
    url(
        r"Projects/project_(?P<pk1>\d+)/sample_(?P<pk2>\d+)",
        PIviews.Sample_main.as_view(),
        name="sample_main",
    ),
    url(
        r"Summary/project_(?P<pk1>\d+)/sample_(?P<pk2>\d+)/run_(?P<pk3>\d+)",
        PIviews.Sample_detail.as_view(),
        name="sample_detail",
    ),
    url(
        "igv_display",
        PIviews.IGV_display,
        name="igv_browser",
    ),  ## get values for IGV
    url(
        "show_igv_<slug:sample_name>/<slug:run_name>/<slug:reference>",
        ajax_views.show_igv,
        name="show_igv",
    ),  ## get values for IGV
    url("download_file", PIviews.download_file, name="download_file"),  ##
    url("download_file_igv", PIviews.download_file_igv, name="download_file_igv"),
    url(
        "download_refmap_files", PIviews.download_file_ref, name="download_refmap_files"
    ),
    url(
        r"Scaffold/project_(?P<project>[a-zA-Z0-9_]+)/sample_(?P<sample>[a-zA-Z0-9_]+)/run_(?P<run>[a-zA-Z0-9_]+)/scaffold_(?P<reference>[a-zA-Z0-9_]+)",
        PIviews.Scaffold_Remap.as_view(),
        name="scaffold_remap",
    ),
    url(
        r"ajax/deploy_ProjectPI$",
        PIajax_views.deploy_ProjectPI,
        name="deploy_ProjectPI",
    ),
    url(
        r"ajax/deploy_televir_map$",
        PIajax_views.deploy_televir_map,
        name="deploy_televir_map",
    ),
]
