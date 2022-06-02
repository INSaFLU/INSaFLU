'''
Created on 31/05/2022

@author: mmp
'''
from django import forms
from managing_files.models import Reference
from datasets.models import Consensus

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