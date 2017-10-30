
from django.conf.urls import url
from .views import ReferenceView, ReferenceAddView
from .views import SamplesView, SamplesAddView

urlpatterns = [
	url(r'references/references$', ReferenceView.as_view(), name='references'),
	url(r'references/reference_add$', ReferenceAddView.as_view(), name='reference-add'),
	url(r'samples/samples$', SamplesView.as_view(), name='samples'),
	url(r'samples/sample_add$', SamplesAddView.as_view(), name='sample-add'),

] 
