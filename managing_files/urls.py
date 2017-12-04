
from django.conf.urls import url
from managing_files import views

urlpatterns = [
	url(r'references/references$', views.ReferenceView.as_view(), name='references'),
	url(r'references/reference_add$', views.ReferenceAddView.as_view(), name='reference-add'),
	url(r'samples/samples$', views.SamplesView.as_view(), name='samples'),
	url(r'samples/sample_add$', views.SamplesAddView.as_view(), name='sample-add'),
	url(r'samples/sample_file$', views.AddValueModal.as_view(), name='sample-file'),
	url(r'samples/sample_fastq$', views.SamplesAddView.as_view(), name='sample-fastq'),
	url(r'samples/(?P<pk>\d+)/sample_description$', views.SamplesDetailView.as_view(), name='sample-description'),
	url(r'project/projects$', views.ProjectsView.as_view(), name='projects'),
	url(r'project/project_add', views.ProjectCreateView.as_view(), name='project-add'),
	url(r'samples/(?P<pk>\d+)/add_sample_project$', views.AddSamplesProjectsView.as_view(), name='add-sample-project'),
	url(r'samples/(?P<pk>\d+)/remove_sample_project$', views.SamplesDetailView.as_view(), name='remove-sample-project'),
	url(r'samples/(?P<pk>\d+)/show_sample_project_result$', views.ShowSampleProjectsView.as_view(), name='show-sample-project-results'),

	###  dialog boxes
	url(r'samples/(?P<pk>\d+)/show-image-coverage$', views.AddValueModal.as_view(), name='show-image-coverage'),
	url(r'samples/sample_dataset$', views.AddValueModal.as_view(), name='sample-dataset'),
	url(r'samples/sample_vaccine$', views.SamplesAddView.as_view(), name='sample-vaccine'),
		
	### ajax functions
	url(r'^ajax/validate_project_reference_name$', views.validate_project_reference_name, name='validate-project-reference'),
	url(r'^ajax/set_check_box_values', views.set_check_box_values, name='set-check-box-values'),
	url(r'^ajax/show_phylo_canvas', views.show_phylo_canvas, name='show-phylo-canvas'),
] 
