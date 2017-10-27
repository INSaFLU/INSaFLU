import django_tables2 as tables
from django_tables2.utils import A
from .models import Reference


class ReferenceTable(tables.Table):
#   Renders a normal value as an internal hyperlink to another page.
#   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
    
    class Meta:
        model = Reference
        fields = ('name', 'scentific_name', 'file_name',
                  'creation_date', 'is_obsolete', 'number_of_locus')
        attrs = {"class": "table-striped table-bordered"}
        empty_text = "There are no References matching the search criteria..."