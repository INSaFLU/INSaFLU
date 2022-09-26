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
        "igv_display",
        PIviews.IGV_display,
        name="igv_browser",
    ),  ## get values for IGV
    url(
        r"projects/(?P<project>[a-zA-Z0-9]+)/all_reports$",
        PIviews.Project_reports,
        name="all_PIproject_reports",
    ),
    url(
        "show_igv_<slug:sample_name>/<slug:run_name>/<slug:reference>",
        ajax_views.show_igv,
        name="show_igv",
    ),  ## get values for IGV
    url("download_file", PIviews.download_file, name="download_file"),  ##
    url("download_file_igv", PIviews.download_file_igv, name="download_file_igv"),
    url(
        r"Projects/(?P<project_name>[a-zA-Z0-9]+)/(?P<sample_name>[a-zA-Z0-9]+)",
        PIviews.Sample_main,
        name="sample_main",
    ),
    url(
        "<slug:project>/sample_<slug:sample>/<slug:name>",
        PIviews.Sample_detail,
        name="sample_detail",
    ),
    url(
        "<slug:project>/sample_<slug:sample>/<slug:run>/<slug:reference>",
        PIviews.Scaffold_Remap,
        name="scaffold_remap",
    ),
    url(
        r"ajax/deploy_ProjectPI$",
        PIajax_views.deploy_ProjectPI,
        name="deploy_ProjectPI",
    ),
]
