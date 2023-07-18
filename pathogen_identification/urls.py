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
        r"check_project_name",
        PIajax_views.validate_project_name,
        name="check_project_name",
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
        r"Projects/(?P<pk1>\d+)/all_reports$",
        PIviews.Project_reports,
        name="all_PIproject_reports",
    ),
    url(
        r"Projects/project_(?P<pk1>\d+)/sample_(?P<pk2>\d+)",
        PIviews.Sample_main.as_view(),
        name="sample_main",
    ),
    url(
        r"Project_Samples/(?P<pk1>\d+)/(?P<pk2>\d+)/all_reports$",
        PIviews.Sample_reports,
        name="all_PIsample_reports",
    ),
    url(
        r"Project_Samples/(?P<pk1>\d+)/(?P<pk2>\d+)/sample_report$",
        PIviews.Sample_ReportCombined.as_view(),
        name="televir_sample_compound_report",
    ),
    url(
        r"Summary/project_(?P<pk1>\d+)/sample_(?P<pk2>\d+)/run_(?P<pk3>\d+)",
        PIviews.Sample_detail.as_view(),
        name="sample_detail",
    ),
    url(
        "igv_display",
        PIajax_views.IGV_display,
        name="igv_browser",
    ),  ## get values for IGV
    url("download_file", PIviews.download_file, name="download_file"),  ##
    url(
        "download_intermediate_reports",
        PIviews.download_intermediate_reports_zipfile,
        name="download_intermediate_reports",
    ),  ##
    url("download_file_igv", PIviews.download_file_igv, name="download_file_igv"),
    url(
        "download_refmap_files", PIviews.download_file_ref, name="download_refmap_files"
    ),
    url(
        r"Scaffold/project_(?P<pk1>\d+)/sample_(?P<pk2>\d+)/run_(?P<pk3>\d+)/scaffold_(?P<reference>[a-zA-Z0-9_]+)",
        PIviews.Scaffold_Remap.as_view(),
        name="scaffold_remap",
    ),
    url(
        r"ajax/deploy_ProjectPI$",
        PIajax_views.deploy_ProjectPI,
        name="deploy_ProjectPI",
    ),
    url(
        r"ajax/deploy_ProjectPI$",
        PIajax_views.deploy_ProjectPI_runs,
        name="deploy_ProjectPI_runs",
    ),
    url(
        r"^ajax/submit_televir_sample$",
        PIajax_views.submit_televir_project_sample,
        name="submit_televir_project_sample",
    ),  ## remove a televir project
    url(
        r"^ajax/submit_televir_sample$",
        PIajax_views.submit_televir_project_sample_runs,
        name="submit_televir_project_sample_runs",
    ),  ## remove a televir project
    url(
        r"ajax/sort_reports$",
        PIajax_views.sort_report_projects,
        name="sort_project_reports",
    ),
    url(
        r"^ajax/kill_televir_sample$",
        PIajax_views.kill_televir_project_tree_sample,
        name="kill_televir_project_sample",
    ),  ## remove a televir project
    url(
        r"ajax/deploy_televir_map$",
        PIajax_views.deploy_televir_map,
        name="deploy_televir_map",
    ),
    url(
        r"ajax/set_sample_control$",
        PIajax_views.set_sample_reports_control,
        name="set_control_televir_project_sample",
    ),
]
