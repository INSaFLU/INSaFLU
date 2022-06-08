from django.db import models

import os

from django.conf import settings
from django.utils.safestring import mark_safe
from django.contrib.auth.models import User
from managing_files.models import ProjectSample, Reference
from fluwebvirus.formatChecker import ContentTypeRestrictedFileField
from constants.constants import Constants, TypePath, FileExtensions
# Create your models here.

def consensus_directory_path(instance, filename):
    # file will be uploaded to MEDIA_ROOT/<filename>
    return 'uploads/generic_data/user_{0}/consensus/{1}'.format(instance.owner.id, filename)

class Consensus(models.Model):
    constants = Constants()
    
    name = models.CharField(max_length=200, db_index=True, verbose_name='Consensus name')
    display_name = models.CharField(max_length=200, db_index=True, default='', verbose_name='Display name')
    creation_date = models.DateTimeField(auto_now_add=True, verbose_name='Uploaded Date')
    
    ## Size 100K
    consensus_fasta = ContentTypeRestrictedFileField(upload_to=consensus_directory_path, content_types=['application/octet-stream',
                                                                                                        'text/plain'],\
                                        max_upload_size=settings.MAX_REF_FASTA_FILE, blank=True, null=True, max_length=500)
    consensus_fasta_name = models.CharField(max_length=200, default='', verbose_name='Fasta file')
    hash_consensus_fasta = models.CharField(max_length=50, blank=True, null=True)

    owner = models.ForeignKey(User, related_name='consensus', blank=True, null=True, on_delete=models.CASCADE)
    is_obsolete = models.BooleanField(default=False, verbose_name='Obsolete')
    is_deleted = models.BooleanField(default=False, verbose_name='Deleted')
    is_deleted_in_file_system = models.BooleanField(default=False, verbose_name='Deleted in FileSystem')
    description = models.CharField(max_length=500, default='', blank=True, null=True, verbose_name='Description')
    number_of_locus = models.IntegerField(default=0, verbose_name='#Locus')

    ### if is deleted in file system
    is_deleted_in_file_system = models.BooleanField(default=False)            ## if this file was removed in file system
    date_deleted = models.DateTimeField(blank=True, null=True, verbose_name='Date attached') ## this date has the time of deleted by web page
    
    def __str__(self):
        return self.name

    def get_consensus_fasta(self, type_path):
        """
        get a path, type_path, from MEDIA_URL or MEDIA_ROOT
        """
        path_to_find = self.consensus_fasta.name
        if (type_path == TypePath.MEDIA_ROOT): 
            if not path_to_find.startswith('/'): path_to_find = os.path.join(getattr(settings, "MEDIA_ROOT", None), path_to_find)
        else:
            path_to_find = os.path.join(getattr(settings, "MEDIA_URL", None), path_to_find)
        return path_to_find
    
    def get_consensus_fasta_index(self, type_path):
        """
        get a fasta fai path, type_path, from MEDIA_URL or MEDIA_ROOT
        """
        return self.get_consensus_fasta(type_path) + FileExtensions.FILE_FAI

    def get_consensus_fasta_web(self):
        """
        return web link for consensus
        """
        out_file = self.get_consensus_fasta(TypePath.MEDIA_ROOT)
        if (os.path.exists(out_file)):
            return mark_safe('<a href="{}" download="{}"> {}</a>'.format(self.get_consensus_fasta(\
                        TypePath.MEDIA_URL), os.path.basename(self.get_consensus_fasta(TypePath.MEDIA_ROOT)),
                        self.constants.short_name(self.name, Constants.SHORT_NAME_LENGTH)))
        return _('File not available.')

 
    class Meta:
        verbose_name = 'Consensus'
        verbose_name_plural = 'Consensus'
        ordering = ['-creation_date', ]
        indexes = [
            models.Index(fields=['name'], name='name_consensus_idx'),
        ]

class Dataset(models.Model):
    """
    Dataset, join several consensus from different projects or individual consensus
    """
    ### has the path for main result
    PATH_MAIN_RESULT = 'main_result'
    
    name = models.CharField(max_length=200, db_index=True, blank=True, null=True, verbose_name='Dataset name')
    owner = models.ForeignKey(User, related_name='dataset', blank=True, null=True, on_delete=models.CASCADE)
    creation_date = models.DateTimeField('Uploaded date', auto_now_add=True)
    last_change_date = models.DateTimeField('Last change date', blank=True, null=True)
    is_deleted = models.BooleanField(default=False)
    number_of_sequences_from_projects = models.SmallIntegerField(default=0)      ### has the number of sequences from projects
    number_of_sequences_from_consensus = models.SmallIntegerField(default=0)      ### has the number of sequences from consensus
    number_of_sequences_from_references = models.SmallIntegerField(default=0)      ### has the number of sequences from references
    is_processed = models.BooleanField(default=False)               ### if some sequence was added and not processed yet
    
    ### if is deleted in file system
    is_deleted_in_file_system = models.BooleanField(default=False)            ## if this file was removed in file system
    date_deleted = models.DateTimeField(blank=True, null=True, verbose_name='Date attached')    ## this date has the time of deleted by web page
    
    ## has the number of alerts
    totla_alerts = models.IntegerField(default=0, verbose_name='Alerts')

    def __str__(self):
        return self.name
    
    class Meta:
        verbose_name = 'Dataset'
        verbose_name_plural = 'Datasets'
        ordering = ['-creation_date', ]

    def get_number_different_references(self):
        """ return the number of different references 
        :param user """
        dt_id_ref = {}
        for dataset_consensus in DatasetConsensus.objects.filter(dataset=self, is_deleted=False):
            if not dataset_consensus.consensus is None: continue
            if not dataset_consensus.reference is None and not dataset_consensus.reference.pk in dt_id_ref:
                dt_id_ref[dataset_consensus.reference.pk] = 1
            if not dataset_consensus.project_sample is None and not dataset_consensus.project_sample.project.reference.pk in dt_id_ref:
                dt_id_ref[dataset_consensus.project_sample.project.reference.pk] = 1
        return len(dt_id_ref)

class DatasetConsensus(models.Model):
    
    dataset = models.ForeignKey(Dataset, related_name='dataset_consensus', blank=True, null=True, on_delete=models.CASCADE)
    ## only can have one of this three
    project_sample = models.ForeignKey(ProjectSample, related_name='dataset_project_samples', blank=True, null=True, on_delete=models.CASCADE)
    reference = models.ForeignKey(Reference, related_name='dataset_reference', blank=True, null=True, on_delete=models.CASCADE)
    consensus = models.ForeignKey(Consensus, related_name='dataset_consensus', blank=True, null=True, on_delete=models.CASCADE)
    is_project_sample_finished = models.BooleanField(default=False)    ## True if all process in ProjectSample Over
    ##
    creation_date = models.DateTimeField('uploaded date', auto_now_add=True)
    is_finished = models.BooleanField(default=False)
    is_error = models.BooleanField(default=False)        ## if some problem occurs
    alert_first_level = models.IntegerField(default=0)    ## has the number of alerts for high errors
    alert_second_level = models.IntegerField(default=0)    ## has the number of alerts for low errors
    seq_name_all_consensus = models.CharField(blank=True, null=True, max_length=200)    ## name of the sequence when saved in AllConsensus.fasta file

    ### remove
    is_deleted = models.BooleanField(default=False)
    date_deleted = models.DateTimeField(blank=True, null=True, verbose_name='Date attached')    ## this date has the time of deleted by web page
    
    def __str__(self):
        return self.dataset.name
    
    class Meta:
        verbose_name = 'Consensus'
        verbose_name_plural = 'Consensus'
        ordering = ['-creation_date', ]




