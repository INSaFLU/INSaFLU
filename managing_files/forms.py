
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, ButtonHolder, Submit
from crispy_forms.layout import Div, Layout, ButtonHolder, Submit
from django import forms
from .models import Reference
from django.core.exceptions import NON_FIELD_ERRORS
from django.utils.translation import ugettext_lazy as _
from django.core.files.temp import NamedTemporaryFile
from constants.Constants import Constants
import os

## https://kuanyui.github.io/2015/04/13/django-crispy-inline-form-layout-with-bootstrap/
class ReferenceForm(forms.ModelForm):
    error_css_class = 'has-error'
    
    class Meta:
        model = Reference
        # specify what fields should be used in this form.
        fields = ('name', 'scientific_name', 'reference_fasta', 'reference_genbank')
        
    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop('request')
        super(ReferenceForm, self).__init__(*args, **kwargs)
 
        ## can exclude explicitly
        ## exclude = ('md5',)
        field_text= [
            # (field_name, Field title label, Detailed field description, requiered)
            ('name', 'Name', 'Regular name for this reference', True),
            ('scientific_name', 'Scientific name', 'Scientific name for this reference', False),
            ('reference_fasta', 'Reference (fasta)', 'File reference in fasta format', True),
            ('reference_genbank', 'Reference (genBank)', 'File reference in genBank format, to have the genes annotated', True),
        ]
        for x in field_text:
            self.fields[x[0]].label = x[1]
            self.fields[x[0]].help_text = x[2]
            self.fields[x[0]].required = x[3]
              
        self.helper = FormHelper()
        self.helper.form_method = 'POST'
        self.helper.layout = Layout(
            Div(
                Div('name', css_class="col-sm-3"),
                Div('scientific_name', css_class="col-sm-7"),
                css_class = 'row'
            ),
            Div('reference_fasta', css_class = 'show-for-sr'),
            Div('reference_genbank', css_class = 'show-for-sr'),
            ButtonHolder(
                Submit('save', 'Save', css_class='btn-primary'),
                Submit('cancel', 'Cancel', css_class='btn-secondary')
            )
        )

    def clean_name(self):
        name = self.cleaned_data['name']
        try:
            Reference.objects.get(name=name, owner__username=self.request.user.username)
            raise forms.ValidationError("This name '" + name +"' already exist in database, please choose other.")
        except Reference.DoesNotExist:
            pass
        return name
    
    def clean_reference_fasta(self):
        """
        Test the Fasta file
        """
        reference_fasta = self.cleaned_data['reference_fasta']
        reference_genbank = self.cleaned_data['reference_genbank']
        if (reference_fasta.name == reference_genbank): raise forms.ValidationError("Error: both files has the same name. Please two different files.")
        
        reference_fasta_temp_file_name = NamedTemporaryFile(prefix='flu_', delete=False)
        reference_fasta_temp_file_name.write(reference_fasta.file.read())
        reference_fasta_temp_file_name.flush()
        reference_fasta_temp_file_name.close()
        constants = Constants()
        try:
            number_locus = constants.is_fasta(reference_fasta_temp_file_name.name)
            self.request.session[Constants.NUMBER_LOCUS_FASTA_FILE] = number_locus
            os.unlink(reference_fasta_temp_file_name.name)
        except IOError as e:    ## (e.errno, e.strerror)
            os.unlink(reference_fasta_temp_file_name.name)
            raise forms.ValidationError(e.args[0])
        return reference_fasta
    
    def clean_reference_genbank(self):
        """
        Test the genBank file
        """
        reference_genbank_temp_file_name = NamedTemporaryFile(prefix='flu_', delete=False)
        reference_genbank = self.cleaned_data['reference_genbank']
        reference_genbank_temp_file_name.write(reference_genbank.read())
        reference_genbank_temp_file_name.flush()
        reference_genbank_temp_file_name.close()
        constants = Constants()
        try:
            constants.is_genbank(reference_genbank_temp_file_name.name)
            os.unlink(reference_genbank_temp_file_name.name)
        except IOError as e:
            os.unlink(reference_genbank_temp_file_name.name)
            raise forms.ValidationError(e.args[0])
        return reference_genbank
    
