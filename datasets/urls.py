from django.urls import path

from datasets import ajax_views, views

urlpatterns = [
    path("datasets", views.DatasetsView.as_view(), name="datasets"),
    path(
        "<int:pk>/add_references_dataset",
        views.AddDatasetsReferencesView.as_view(),
        name="add-references-dataset",
    ),
    path(
        "<int:pk>/add_consensus_dataset",
        views.AddDatasetsConsensusView.as_view(),
        name="add-consensus-dataset",
    ),
    path(
        "<int:pk>/upload_new_consensus",
        views.UploadNewConsensusView.as_view(),
        name="upload-new-consensus",
    ),
    path(
        "<int:pk>/add_projects_dataset",
        views.AddDatasetsProjectsView.as_view(),
        name="add-projects-dataset",
    ),
    path(
        "<int:pk>/show_dataset_consensus",
        views.ShowDatasetsConsensusView.as_view(),
        name="show-dataset-consensus",
    ),
    path(
        "<int:pk>/show_dataset_settings",
        views.DatasetsSettingsView.as_view(),
        name="dataset-settings",
    ),
    path(
        "<int:pk>/dataset_update_metadata",
        views.UpdateMetadataDataset.as_view(),
        name="dataset-update-metadata",
    ),  ## upload new matadata to replace the exist
    path(
        "<int:pk>/dataset_update_metadata_file",
        views.AddSingleMetadataDatasetFile.as_view(),
        name="dataset-add-single-file-metadata",
    ),  ## upload new matadata to replace the exist
    path(
        "ajax/remove_dataset", ajax_views.remove_dataset, name="remove_dataset"
    ),  ## remove a dataset
    path(
        "ajax/add_dataset_name", ajax_views.add_dataset_name, name="add_dataset_name"
    ),  ## add a dataset
    path(
        "ajax/test_dataset_name",
        ajax_views.test_dataset_name,
        name="test_dataset_name",
    ),  ## test if a dataset name exists
    path(
        "ajax/test_consensus_name",
        ajax_views.test_consensus_name,
        name="test_consensus_name",
    ),  ## test if a consensus name exists
    path(
        "ajax/add_consensus_name",
        ajax_views.add_consensus_name,
        name="add_consensus_name",
    ),  ## add a consensus name and file
    path(
        "ajax/remove_consensus", ajax_views.remove_consensus, name="remove_consensus"
    ),  ## remove a consensus name and file
    path(
        "ajax/remove_consensus_in_dataset",
        ajax_views.remove_consensus_in_dataset,
        name="remove_consensus_in_dataset",
    ),  ## remove a consensus dataset
    path(
        "ajax/validate_consensus_name",
        ajax_views.validate_consensus_name,
        name="validate-consensus-name",
    ),  ## Validate consensus name
    path(
        "ajax/dataset_rebuild", ajax_views.dataset_rebuild, name="dataset_rebuild"
    ),  ## rebuild results of a dataset
    path(
        "ajax/show_phylo_canvas",
        ajax_views.show_phylo_canvas,
        name="show-phylo-canvas-datasets",
    ),
    path(
        "ajax/show_msa_nucleotide",
        ajax_views.show_msa_nucleotide,
        name="show-msa-nucleotide-datasets",
    ),
]
