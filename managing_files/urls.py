
from django.conf.urls import url
from managing_files import views, ajax_views

urlpatterns = [
	url(r'references/references$', views.ReferenceView.as_view(), name='references'),
	url(r'references/reference_add$', views.ReferenceAddView.as_view(), name='reference-add'),
	url(r'samples/samples$', views.SamplesView.as_view(), name='samples'),
	url(r'samples/sample_add$', views.SamplesAddView.as_view(), name='sample-add'),
	url(r'samples/sample_file$', views.SamplesAddView.as_view(), name='sample-file'),	## add xls file with several samples 
	url(r'samples/sample_fastq$', views.SamplesAddView.as_view(), name='sample-fastq'),	## add several fastq.gz
	url(r'samples/(?P<pk>\d+)/sample_description$', views.SamplesDetailView.as_view(), name='sample-description'),
	url(r'project/projects$', views.ProjectsView.as_view(), name='projects'),
	url(r'project/project_add', views.ProjectCreateView.as_view(), name='project-add'),
	url(r'samples/(?P<pk>\d+)/add_sample_project$', views.AddSamplesProjectsView.as_view(), name='add-sample-project'),
	url(r'samples/(?P<pk>\d+)/remove_sample_project$', views.SamplesDetailView.as_view(), name='remove-sample-project'),
	url(r'samples/(?P<pk>\d+)/show_sample_project_result$', views.ShowSampleProjectsView.as_view(), name='show-sample-project-results'),

	### ajax functions
	url(r'^ajax/validate_project_reference_name$', ajax_views.validate_project_reference_name, name='validate-project-reference'),
	url(r'^ajax/set_check_box_values', ajax_views.set_check_box_values, name='set-check-box-values'),
	url(r'^ajax/show_phylo_canvas', ajax_views.show_phylo_canvas, name='show-phylo-canvas'),
	url(r'^ajax/get_image_coverage', ajax_views.get_image_coverage, name='get_image_coverage'),
	url(r'^ajax/add_single_value_database', ajax_views.add_single_value_database, name='add_single_value_database'),	## add a single value to a table in database
	url(r'^ajax/remove_single_value_database', ajax_views.remove_single_value_database, name='remove_single_value_database'),	## add a single value to a table in database
] 
