
from django.conf.urls import url
from datasets import views, ajax_views

urlpatterns = [
    # url(r'datasets/consensus$', views.ConsensusView.as_view(), name='consensus'),
    # url(r'datasets/consensus_add$', views.ConsensusAddView.as_view(), name='consensus-add'),
    url(r'datasets/datasets$', views.DatasetsView.as_view(), name='datasets'),
    url(r'datasets/(?P<pk>\d+)/add_references_dataset$', views.AddDatasetsReferencesView.as_view(), name='add-references-dataset'),
    url(r'datasets/(?P<pk>\d+)/add_consensus_dataset$', views.AddDatasetsConsensusView.as_view(), name='add-consensus-dataset'),
    url(r'datasets/(?P<pk>\d+)/upload_new_consensus$', views.UploadNewConsensusView.as_view(), name='upload-new-consensus'),
    url(r'datasets/(?P<pk>\d+)/add_projects_dataset$', views.AddDatasetsProjectsView.as_view(), name='add-projects-dataset'),
    url(r'datasets/(?P<pk>\d+)/show_sequences_dataset$', views.ShowDatasetsSequencesView.as_view(), name='show-sequences-dataset'),
    url(r'datasets/(?P<pk>\d+)/show_dataset_consensus$', views.ShowDatasetsConsensusView.as_view(), name='show-dataset-consensus'),
    # url(r'datasets/datasets_projects_add$', views.DatasetsProjectsAddView.as_view(), name='datasets-projects-add'),
    # url(r'datasets/datasets_consensus_add$', views.DatasetsConsensusAddView.as_view(), name='datasets-consensus-add'),
    
    url(r'^ajax/remove_dataset$', ajax_views.remove_dataset, name='remove_dataset'),            ## remove a dataset
    url(r'^ajax/add_dataset_name$', ajax_views.add_dataset_name, name='add_dataset_name'),      ## add a dataset 
    url(r'^ajax/test_dataset_name$', ajax_views.test_dataset_name, name='test_dataset_name'),   ## test if a dataset name exists
    url(r'^ajax/test_consensus_name$', ajax_views.test_consensus_name, name='test_consensus_name'),   ## test if a consensus name exists
    url(r'^ajax/add_consensus_name$', ajax_views.add_consensus_name, name='add_consensus_name'),      ## add a consensus name and file
    url(r'^ajax/remove_consensus$', ajax_views.remove_consensus, name='remove_consensus'),      ## add a consensus name and file
    url(r'^ajax/remove_consensus_in_dataset$', ajax_views.remove_consensus_in_dataset, name='remove_consensus_in_dataset'),      ## add a consensus name and file
    url(r'^ajax/validate_consensus_name$', ajax_views.validate_consensus_name, name='validate-consensus-name'), ## Validate consensus name 
    
    
]