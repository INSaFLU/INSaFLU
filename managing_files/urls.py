
from django.conf.urls import url
from managing_files.views import ReferenceView, ReferenceAddView
from managing_files.views import SamplesView, SamplesAddView, AddValueModal, SamplesDetailView

urlpatterns = [
	url(r'references/references$', ReferenceView.as_view(), name='references'),
	url(r'references/reference_add$', ReferenceAddView.as_view(), name='reference-add'),
	url(r'samples/samples$', SamplesView.as_view(), name='samples'),
	url(r'samples/sample_add$', SamplesAddView.as_view(), name='sample-add'),
	url(r'samples/sample_file$', AddValueModal.as_view(), name='sample-file'),
	url(r'samples/sample_fastq$', SamplesAddView.as_view(), name='sample-fastq'),
	url(r'samples/sample_dataset$', AddValueModal.as_view(), name='sample-dataset'),
	url(r'samples/sample_vaccine$', SamplesAddView.as_view(), name='sample-vaccine'),
	url(r'samples/(?P<pk>\d+)/sample_description$', SamplesDetailView.as_view(), name='sample-description'),
] 
