from django.db import models

import os

from django.conf import settings
from django.utils.safestring import mark_safe
from django.contrib.auth.models import User
from fluwebvirus.formatChecker import ContentTypeRestrictedFileField
from constants.constants import Constants, TypePath, FileExtensions, FileType
# Create your models here.

def consensus_directory_path(instance, filename):
    # file will be uploaded to MEDIA_ROOT/<filename>
    return 'uploads/generic_data/user_{0}/consensus/{1}'.format(instance.owner.id, filename)

class Consensus(models.Model):
    name = models.CharField(max_length=200, db_index=True, verbose_name='Consensus name')
    display_name = models.CharField(max_length=200, db_index=True, default='', verbose_name='Display name')
    creation_date = models.DateTimeField(auto_now_add=True, verbose_name='Uploaded Date')
    
    ## Size 100K
    consensus_fasta = ContentTypeRestrictedFileField(upload_to=consensus_directory_path, content_types=['application/octet-stream'],\
                                        max_upload_size=settings.MAX_REF_FASTA_FILE, blank=True, null=True, max_length=500)
    consensus_fasta_name = models.CharField(max_length=200, default='', verbose_name='Fasta file')
    hash_consensus_fasta = models.CharField(max_length=50, blank=True, null=True)

    owner = models.ForeignKey(User, related_name='consensus', blank=True, null=True, on_delete=models.CASCADE)
    is_obsolete = models.BooleanField(default=False, verbose_name='Obsolete')
    is_deleted = models.BooleanField(default=False, verbose_name='Deleted')
    description = models.CharField(max_length=500, default='', blank=True, null=True, verbose_name='Description')

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
                        self.name))
        return _('File not available.')

 
    class Meta:
        verbose_name = 'Consensus'
        verbose_name_plural = 'Consensus'
        ordering = ['-creation_date', ]
        indexes = [
            models.Index(fields=['name'], name='name_consensus_idx'),
        ]
