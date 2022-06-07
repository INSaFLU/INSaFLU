'''
Created on 31/05/2022

@author: mmp
'''
import os
from django import forms
from managing_files.models import Reference, Project
from datasets.models import Consensus
from utils.utils import Utils
from django.conf import settings
from django.template.defaultfilters import filesizeformat
from crispy_forms.helper import FormHelper
from django.core.files.temp import NamedTemporaryFile
from crispy_forms.layout import Div, Layout, ButtonHolder, Submit, Button
from django.urls import reverse
from constants.constants import Constants
from utils.software import Software

class AddReferencesDatasetForm(forms.ModelForm):
    """
    References model to dataset 
    """
    error_css_class = 'error'
    
    class Meta:
        model = Reference
        exclude = ()

    def __init__(self, *args, **kwargs):
        super(AddReferencesDatasetForm, self).__init__(*args, **kwargs)
        
    def clean(self):
        """
        """
        cleaned_data = super(AddReferencesDatasetForm, self).clean()
        return cleaned_data
    
    
class AddConsensusDatasetForm(forms.ModelForm):
    """
    Add consensus to Insaflu 
    """
    error_css_class = 'error'
    
    class Meta:
        model = Consensus
        exclude = ()

    def __init__(self, *args, **kwargs):
        super(AddConsensusDatasetForm, self).__init__(*args, **kwargs)
        
    def clean(self):
        """
        """
        cleaned_data = super(AddConsensusDatasetForm, self).clean()
        return cleaned_data

class AddProjectsDatasetForm(forms.ModelForm):
    """
    Add project to Insaflu 
    """
    error_css_class = 'error'
    
    class Meta:
        model = Project
        exclude = ()

    def __init__(self, *args, **kwargs):
        super(AddProjectsDatasetForm, self).__init__(*args, **kwargs)
        
    def clean(self):
        """
        """
        cleaned_data = super(AddProjectsDatasetForm, self).clean()
        return cleaned_data
    
## https://kuanyui.github.io/2015/04/13/django-crispy-inline-form-layout-with-bootstrap/
class ConsensusForm(forms.ModelForm):
    """
    Reference form, name, isolate_name and others
    """
    software = Software()
    utils = Utils()
    error_css_class = 'error'
    
    class Meta:
        model = Consensus
        # specify what fields should be used in this form.
        fields = ('name', 'consensus_fasta')
        
    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop('request')
        self.pk = kwargs.pop('pk')
        super(ConsensusForm, self).__init__(*args, **kwargs)

        ## can exclude explicitly
        ## exclude = ('md5',)
        field_text= [
            # (field_name, Field title label, Detailed field description, requiered)
            ('name', 'Name', 'Regular name for this reference', True),
            ('consensus_fasta', 'Consensus (FASTA)', "Consensus file in fasta format.<br>" +\
                "Max total sequence length: {}<br>".format(filesizeformat(settings.MAX_LENGTH_SEQUENCE_TOTAL_FROM_FASTA))  +\
                "Max FASTA file size: {}".format(filesizeformat(settings.MAX_REF_FASTA_FILE)), True),
        ]
        for x in field_text:
            self.fields[x[0]].label = x[1]
            self.fields[x[0]].help_text = x[2]
            self.fields[x[0]].required = x[3]

## in case you have to undo to specific ID
##        cancel_url = reverse('references')
##        if self.instance: cancel_url = reverse('references', kwargs={'pk': self.instance.id})

        self.helper = FormHelper()
        self.helper.form_method = 'POST'
        self.helper.attrs["data-validate-consensus-url"] = reverse('validate-consensus-name')
        self.helper.attrs["id"] = "id_form_consensus"
        self.helper.layout = Layout(
            Div(
                Div('name', css_class="col-sm-7"),
                css_class = 'row'
            ),
            Div('consensus_fasta', css_class = 'show-for-sr'),
            ButtonHolder(
                Submit('save', 'Save', css_class='btn-primary', onclick="this.disabled=true,this.form.submit();"),
                Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(
                    reverse('add-consensus-dataset', args=[self.pk])))
            )
        )

    def clean(self):
        """
        Clean all together because it's necessary to compare the genbank and fasta files
        """
        cleaned_data = super(ConsensusForm, self).clean()
        name = self.cleaned_data.get('name', '').strip()
        if (len(name) == 0):
            self.add_error('name', "Error: You must give a unique name.")
            return cleaned_data

        try:
            Consensus.objects.get(name__iexact=name, owner=self.request.user, is_obsolete=False, is_deleted=False)
            self.add_error('name', "This name '" + name +"' already exist in database, please choose other.")
        except Consensus.DoesNotExist:
            pass
        
        ## test reference_fasta
        if ('consensus_fasta' not in cleaned_data):
            self.add_error('consensus_fasta', "Error: Must have a file.")
            return cleaned_data
        
        ### testing file names
        reference_fasta = cleaned_data['consensus_fasta']
        
        ## testing fasta
        some_error_in_files = False
        reference_fasta_temp_file_name = NamedTemporaryFile(prefix='flu_fa_', delete=False)
        reference_fasta_temp_file_name.write(reference_fasta.read())
        reference_fasta_temp_file_name.flush()
        reference_fasta_temp_file_name.close()
        self.software.dos_2_unix(reference_fasta_temp_file_name.name)
        try:
            number_locus = self.utils.is_fasta(reference_fasta_temp_file_name.name)
            self.request.session[Constants.NUMBER_LOCUS_FASTA_FILE] = number_locus
            
            ## test the max numbers
            if (number_locus > Constants.MAX_SEQUENCES_FROM_FASTA):
                self.add_error('consensus_fasta', _('Max allow number of sequences in fasta: {}'.format(Constants.MAX_SEQUENCES_FROM_FASTA)))
                some_error_in_files = True
            total_length_fasta = self.utils.get_total_length_fasta(reference_fasta_temp_file_name.name)
            if (not some_error_in_files and total_length_fasta > settings.MAX_LENGTH_SEQUENCE_TOTAL_FROM_FASTA):
                some_error_in_files = True
                self.add_error('consensus_fasta', _('The max sum length of the sequences in fasta: {}'.format(settings.MAX_LENGTH_SEQUENCE_TOTAL_FROM_FASTA)))
            
            n_seq_name_bigger_than = self.utils.get_number_seqs_names_bigger_than(reference_fasta_temp_file_name.name, Constants.MAX_LENGTH_SEQ_NAME)
            if (not some_error_in_files and n_seq_name_bigger_than > 0):
                some_error_in_files = True
                if (n_seq_name_bigger_than == 1):
                    self.add_error('consensus_fasta', _('There is one sequence name length bigger than {0}. The max. length name is {0}.'.format(Constants.MAX_LENGTH_SEQ_NAME)))
                else:
                    self.add_error('consensus_fasta', _('There are {0} sequences with name length bigger than {1}. The max. length name is {1}.'.format(n_seq_name_bigger_than, Constants.MAX_LENGTH_SEQ_NAME)))
                    
                ## if some errors in the files, fasta or genBank, return
                if (some_error_in_files): return cleaned_data
            
        except IOError as e:    ## (e.errno, e.strerror)
            os.unlink(reference_fasta_temp_file_name.name)
            some_error_in_files = True
            self.add_error('consensus_fasta', e.args[0])
        except ValueError as e:    ## (e.errno, e.strerror)
            os.unlink(reference_fasta_temp_file_name.name)
            some_error_in_files = True
            self.add_error('consensus_fasta', e.args[0])
        except: 
            os.unlink(reference_fasta_temp_file_name.name)
            some_error_in_files = True
            self.add_error('consensus_fasta', "Not a valid 'fasta' file.")

        ### test if it has degenerated bases
        # if (os.path.exists(reference_fasta_temp_file_name.name)):
        #     try:
        #         self.utils.has_degenerated_bases(reference_fasta_temp_file_name.name)
        #     except Exception as e:
        #         os.unlink(reference_fasta_temp_file_name.name)
        #         some_error_in_files = True
        #         self.add_error('consensus_fasta', e.args[0])
            
        ## if some errors in the files, fasta or genBank, return
        if (some_error_in_files): return cleaned_data
        
        ## remove temp files
        os.unlink(reference_fasta_temp_file_name.name)
        return cleaned_data

