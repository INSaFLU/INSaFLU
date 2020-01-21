
from django.conf.urls import url
from managing_files import views, ajax_views

urlpatterns = [
	url(r'references/references$', views.ReferenceView.as_view(), name='references'),
	url(r'references/reference_add$', views.ReferenceAddView.as_view(), name='reference-add'),
	url(r'samples/samples$', views.SamplesView.as_view(), name='samples'),
	url(r'samples/sample_add$', views.SamplesAddView.as_view(), name='sample-add'),
	url(r'samples/sample_add_file$', views.SamplesAddDescriptionFileView.as_view(), name='sample-add-file'),	## add xls file with several samples 
	url(r'samples/sample_add_single_csv_file$', views.SamplesUploadDescriptionFileView.as_view(), name='sample-add-single-csv-file'),	## upload xls file with several samples 
	url(r'samples/(?P<pk>\d+)/sample_description$', views.SamplesDetailView.as_view(), name='sample-description'),
	url(r'samples/sample_add_fastq$', views.SamplesAddFastQView.as_view(), name='sample-add-fastq'),			## add several fastq.gz
	url(r'samples/sample_update_metadata$', views.SamplesUpdateMetadata.as_view(), name='sample-update-metadata'),				## upload new matadata to replace the exist one
	url(r'samples/sample_add_single_csv_file_metadata$', views.SamplesUploadDescriptionFileViewMetadata.as_view(), name='sample-add-single-csv-file-metadata'),	## upload xls file with new metadata for several samples
	url(r'samples/sample_upload_fastq$', views.SamplesUploadFastQView.as_view(), name='sample-upload-fastq'),		## upload several fastq.gz
	
	url(r'project/projects$', views.ProjectsView.as_view(), name='projects'),
	url(r'project/project_add$', views.ProjectCreateView.as_view(), name='project-add'),
	url(r'project_samples/(?P<pk>\d+)/add_sample_project$', views.AddSamplesProjectsView.as_view(), name='add-sample-project'),
	url(r'project_samples/(?P<pk>\d+)/remove_sample_project$', views.SamplesDetailView.as_view(), name='remove-sample-project'),
	url(r'project_samples/(?P<pk>\d+)/show_sample_project_results$', views.ShowSampleProjectsView.as_view(), name='show-sample-project-results'),
	url(r'project_samples/(?P<pk>\d+)/show_sample_project_single_details$', views.ShowSampleProjectsDetailsView.as_view(), name='show-sample-project-single-detail'),

	### ajax functions
	url(r'^ajax/validate_project_reference_name$', ajax_views.validate_project_reference_name, name='validate-project-reference'),
	url(r'^ajax/set_check_box_values$', ajax_views.set_check_box_values, name='set-check-box-values'),
	url(r'^ajax/show_phylo_canvas$', ajax_views.show_phylo_canvas, name='show-phylo-canvas'),
	url(r'^ajax/show_msa_nucleotide$', ajax_views.show_msa_nucleotide, name='show-msa-nucleotide'),
	url(r'^ajax/show_msa_protein$', ajax_views.show_msa_protein, name='show-msa-protein'),
	url(r'^ajax/show_count_variations$', ajax_views.show_count_variations, name='show-count-variations'),
	url(r'^ajax/get_cds_from_element$', ajax_views.get_cds_from_element, name='get-cds-from-element'),
	url(r'^ajax/get_image_coverage$', ajax_views.get_image_coverage, name='get_image_coverage'),
	url(r'^ajax/add_single_value_database$', ajax_views.add_single_value_database, name='add_single_value_database'),	## add a single value to a table in database
	url(r'^ajax/remove_single_value_database$', ajax_views.remove_single_value_database, name='remove_single_value_database'),	## add a single value to a table in database
	url(r'^ajax/show_igv$', ajax_views.show_igv, name='show_igv'),	## get values for IGV
	
	## remove
	url(r'^ajax/remove_reference$', ajax_views.remove_reference, name='remove_reference'),	## remove a reference
	url(r'^ajax/remove_sample$', ajax_views.remove_sample, name='remove_sample'),			## remove a sample
	url(r'^ajax/remove_project$', ajax_views.remove_project, name='remove_project'),			## remove a project
	url(r'^ajax/remove_project_sample$', ajax_views.remove_project_sample, name='remove_project_sample'),	## remove a project sample
	url(r'^ajax/remove_uploaded_file$', ajax_views.remove_uploaded_file, name='remove_uploaded_file'),	## remove remove_uploaded_file
	url(r'^ajax/get_process_running$', ajax_views.get_process_running, name='get_process_running'),		## get process to run
	url(r'^ajax/submit_sge', ajax_views.submit_sge, name='submit-sge'),		## only for tests because of sge/apache in centos
] 
