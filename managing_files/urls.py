from django.urls import path, re_path

from managing_files import ajax_views, views

urlpatterns = [
    path("references/references", views.ReferenceView.as_view(), name="references"),
    path(
        "references/reference_add",
        views.ReferenceAddView.as_view(),
        name="reference-add",
    ),
    path("samples/samples", views.SamplesView.as_view(), name="samples"),
    path(
        "samples/<int:pk>/show_sample_settings",
        views.SampleSettingsView.as_view(),
        name="sample-settings",
    ),
    path("samples/sample_add", views.SamplesAddView.as_view(), name="sample-add"),
    path(
        "samples/sample_add_file",
        views.SamplesAddDescriptionFileView.as_view(),
        name="sample-add-file",
    ),  ## add xls file with several samples
    path(
        "samples/sample_add_single_csv_file",
        views.SamplesUploadDescriptionFileView.as_view(),
        name="sample-add-single-csv-file",
    ),  ## upload xls file with several samples
    path(
        "samples/<int:pk>/sample_description",
        views.SamplesDetailView.as_view(),
        name="sample-description",
    ),
    path(
        "samples/sample_add_fastq",
        views.SamplesAddFastQView.as_view(),
        name="sample-add-fastq",
    ),  ## add several fastq.gz
    path(
        "samples/sample_update_metadata",
        views.SamplesUpdateMetadata.as_view(),
        name="sample-update-metadata",
    ),  ## upload new matadata to replace the exist one
    path(
        "samples/sample_add_single_csv_file_metadata",
        views.SamplesUploadDescriptionFileViewMetadata.as_view(),
        name="sample-add-single-csv-file-metadata",
    ),  ## upload xls file with new metadata for several samples
    path(
        "samples/sample_upload_fastq",
        views.SamplesUploadFastQView.as_view(),
        name="sample-upload-fastq",
    ),  ## upload several fastq.gz
    re_path(r"project-index", views.ProjectIndex.as_view(), name="project-index"),
    re_path(
        r"references-index", views.ReferencesIndex.as_view(), name="references-index"
    ),
    path("project/projects", views.ProjectsView.as_view(), name="projects"),
    path("project/project_add", views.ProjectCreateView.as_view(), name="project-add"),
    path(
        "project/<int:pk>/show_project_settings",
        views.ProjectsSettingsView.as_view(),
        name="project-settings",
    ),
    path(
        r"project/(?P<pk>\d+)/show_project_settings_setup$",
        views.ProjectsSettingsSetupView.as_view(),
        name="project-settings-setup",
    ),
    path(
        r"project_samples/(?P<pk>\d+)/(?P<tf>\d+)/add_sample_project$",
        views.AddSamplesProjectsView.as_view(),
        name="add-sample-project",
    ),
    path(
        "project_samples/<int:pk>/remove_sample_project",
        views.SamplesDetailView.as_view(),
        name="remove-sample-project",
    ),
    path(
        "project_samples/<int:pk>/show_sample_project_results",
        views.ShowSampleProjectsView.as_view(),
        name="show-sample-project-results",
    ),
    path(
        "project_samples/<int:pk>/show_project_sample_settings",
        views.SampleProjectsSettingsView.as_view(),
        name="sample-project-settings",
    ),
    path(
        "project_samples/<int:pk>/show_sample_project_single_details",
        views.ShowSampleProjectsDetailsView.as_view(),
        name="show-sample-project-single-detail",
    ),
    ### ajax functions
    path(
        "ajax/validate_project_reference_name",
        ajax_views.validate_project_reference_name,
        name="validate-project-reference",
    ),
    path(
        "ajax/validate_reference_name",
        ajax_views.validate_reference_name,
        name="validate-reference-name",
    ),
    path(
        "ajax/set_check_box_values",
        ajax_views.set_check_box_values,
        name="set-check-box-values",
    ),
    path(
        "ajax/show_phylo_canvas",
        ajax_views.show_phylo_canvas,
        name="show-phylo-canvas",
    ),
    path(
        "ajax/show_variants_as_a_table",
        ajax_views.show_variants_as_a_table,
        name="show-variants-as-a-table",
    ),
    path(r"^ajax/show_aln2pheno$", ajax_views.show_aln2pheno, name="show-aln2pheno"),
    path(r"^ajax/show_flumut$", ajax_views.show_flumut, name="show-flumut"),
    path(
        r"^ajax/show_coverage_as_a_table$",
        ajax_views.show_coverage_as_a_table,
        name="show-coverage-as-a-table",
    ),
    path(
        "ajax/show_msa_nucleotide",
        ajax_views.show_msa_nucleotide,
        name="show-msa-nucleotide",
    ),
    path("ajax/show_msa_protein", ajax_views.show_msa_protein, name="show-msa-protein"),
    path(
        "ajax/show_count_variations",
        ajax_views.show_count_variations,
        name="show-count-variations",
    ),
    path(
        "ajax/get_cds_from_element",
        ajax_views.get_cds_from_element,
        name="get-cds-from-element",
    ),
    path(
        "ajax/get_image_coverage",
        ajax_views.get_image_coverage,
        name="get_image_coverage",
    ),
    re_path(
        r"^ajax/update_project_pangolin",
        ajax_views.update_project_pangolin,
        name="update_project_pangolin",
    ),
    path(
        r"^ajax/update_project_mutation_report",
        ajax_views.update_project_mutation_report,
        name="update_project_mutation_report",
    ),
    path(
        r"^ajax/add_single_value_database$",
        ajax_views.add_single_value_database,
        name="add_single_value_database",
    ),  ## add a single value to a table in database
    path(
        "ajax/remove_single_value_database",
        ajax_views.remove_single_value_database,
        name="remove_single_value_database",
    ),  ## add a single value to a table in database
    path("ajax/show_igv", ajax_views.show_igv, name="show_igv"),  ## get values for IGV
    ## remove
    path(
        "ajax/remove_reference", ajax_views.remove_reference, name="remove_reference"
    ),  ## remove a reference
    path(
        "ajax/remove_sample", ajax_views.remove_sample, name="remove_sample"
    ),  ## remove a sample
    path(
        r"^ajax/swap_technology$", ajax_views.swap_technology, name="swap_technology"
    ),  ## swap technology in sample
    path(
        r"^ajax/remove_project$", ajax_views.remove_project, name="remove_project"
    ),  ## remove a project
    path(
        "ajax/remove_televir_project",
        ajax_views.remove_televir_project,
        name="remove_televir_project",
    ),  ## remove a televir project
    path(
        "ajax/remove_project_sample",
        ajax_views.remove_project_sample,
        name="remove_project_sample",
    ),  ## remove a project sample
    path(
        "ajax/remove_televir_project_sample",
        ajax_views.remove_televir_project_sample,
        name="remove_televir_project_sample",
    ),  ## remove a project sample
    path(
        "ajax/remove_uploaded_file",
        ajax_views.remove_uploaded_file,
        name="remove_uploaded_file",
    ),  ## remove remove_uploaded_file
    path(
        "ajax/remove_uploaded_files",
        ajax_views.remove_uploaded_files,
        name="remove_uploaded_files",
    ),  ## remove remove_uploaded_files, several at once
    path(
        "ajax/remove_unattached_samples",
        ajax_views.remove_unattached_samples,
        name="remove_unattached_samples",
    ),  ## remove unattached samples
    path(
        "ajax/relink_uploaded_file",
        ajax_views.relink_uploaded_files,
        name="relink_uploaded_files",
    ),  ## remove remove_uploaded_file
    path(
        r"^ajax/unlock_sample_file$",
        ajax_views.unlock_sample_file,
        name="unlock_sample_file",
    ),  ## unlock sample list files
    path(
        "ajax/get_process_running",
        ajax_views.get_process_running,
        name="get_process_running",
    ),  ## get process to run
    path(
        "ajax/submit_sge", ajax_views.submit_sge, name="submit-sge"
    ),  ## only for tests because of sge/apache in centos
]
