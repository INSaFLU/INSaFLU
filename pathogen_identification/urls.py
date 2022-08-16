from django.conf.urls import url
from managing_files import ajax_views, views

import pathogen_identification.views as PIviews

urlpatterns = [
    url(
        r"projects$",
        PIviews.PathId_ProjectsView.as_view(),
        name="PIprojects_main",
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
        r"Project_samples/(?P<pk>\d+)/remove_sample_project$",
        views.SamplesDetailView.as_view(),
        name="remove-sample-PIproject",
    ),
    url(
        "igv_display",
        PIviews.IGV_display,
        name="igv_browser",
    ),  ## get values for IGV
    url(
        "project_<project>/all_reports",
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
    url("<slug:project_name>", PIviews.MainPage, name="PIproject_samples"),
    url(
        "<slug:project_name>/sample_<slug:sample_name>",
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
        r"^ajax/add_single_value_database$",
        ajax_views.add_single_value_database,
        name="add_single_value_database",
    ),  ## add a single value to a table in database
    url(
        r"^ajax/remove_single_value_database$",
        ajax_views.remove_single_value_database,
        name="remove_single_value_database",
    ),  ## add a single value to a table in database
    url(
        r"^ajax/show_igv$", ajax_views.show_igv, name="show_igv"
    ),  ## get values for IGV
    url(
        r"^ajax/remove_project$", ajax_views.remove_project, name="remove_project"
    ),  ## remove a project
    url(
        r"^ajax/unlock_sample_file$",
        ajax_views.unlock_sample_file,
        name="unlock_sample_file",
    ),  ## unlock sample list files
    url(
        r"^ajax/get_process_running$",
        ajax_views.get_process_running,
        name="get_process_running",
    ),  ## get process to run
    url(r"^ajax/submit_sge$", ajax_views.submit_sge, name="submit-sge"),
]
