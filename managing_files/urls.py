
from django.conf.urls import url

from .views import ReferenceView, ReferenceAddView

urlpatterns = [
    url(r'references/references$', ReferenceView.as_view(), name='references'),
    url(r'references/reference_add$', ReferenceAddView.as_view(), name='reference-add'),
]