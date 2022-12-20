from django.db import models

import os

from django.conf import settings
from django.utils.safestring import mark_safe
from django.contrib.auth.models import User
from managing_files.models import ProjectSample, Reference, user_directory_path
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
                                        max_upload_size=settings.MAX_CONSENSUS_FASTA_FILE, blank=True, null=True, max_length=500)
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
        return 'File not available.'

 
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
    
    DATASET_FILE_NAME_MAFFT = "Nextstrain_Alignment_All.fasta"
    DATASET_FILE_NAME_FASTTREE = "Nextstrain_Tree_All.nwk"
    DATASET_FILE_NAME_FASTTREE_tree = "Nextstrain_Tree_All.tree"
    DATASET_FILE_NAME_nex = "Nextstrain_Alignment_All.nex"

    DATASET_FILE_NAME_RESULT_TSV = "Dataset_list.tsv"     ### first column ID instead of 'sample name' to be compatible with Phandango e Microreact
    DATASET_FILE_NAME_RESULT_NEXTSTRAIN_TSV = "Nextstrain_metadata.tsv"     ### nextstrain list, mandatory fields are in db/nextstrain/ncov/references_metadata.tsv
    DATASET_FILE_NAME_RESULT_NEXTSTRAIN_CSV = "Nextstrain_metadata.csv"     ### nextstrain list, mandatory fields are in db/nextstrain/ncov/references_metadata.tsv
    DATASET_FILE_NAME_RESULT_CSV = "Dataset_list.csv"     ### first column ID instead of 'sample name' to be compatible with Phandango e Microreact
    DATASET_FILE_NAME_RESULT_json = "Dataset_list_simple.json"     ### first column ID instead of 'sample name' to be compatible with Phandango e Microreact, to download to 
    DATASET_FILE_NAME_RESULT_all_consensus = "AllConsensus.fasta"     ### all consensus sequences for a project sample
    
    ###NextStrain Expected, apear inside 'auspice' folder whenrun NextStrain
    REFERENCE_NAME = "Wuhan/Hu-1/2019"      ### need to change in the future
    RUN_out_path = 'auspice'
    DATASET_FILE_NAME_nextstrain_auspice_zip = "auspice_json.zip"
    #DATASET_FILE_NAME_nextstrain_default_build = "ncov_default-build.json"
    #DATASET_FILE_NAME_nextstrain_build_root = "ncov_default-build_root-sequence.json"
    #DATASET_FILE_NAME_nextstrain_build_tip = "ncov_default-build_tip-frequencies.json"
    #DATASET_FILE_NAME_nextstrain_log = "ncov_default-build.json"
    
    DATASET_FILE_NAME_nextstrain_error = "NextStrainError.txt"      ## has the error of the nextStrain
    
    ####
    #VECT_files_next_strain = [
    #    DATASET_FILE_NAME_nextstrain_default_build,
    #    DATASET_FILE_NAME_nextstrain_build_root,
    #   DATASET_FILE_NAME_nextstrain_build_tip
    #]
    
    ### files to zip
    VECT_files_to_zip = [ 
        DATASET_FILE_NAME_RESULT_TSV,
        DATASET_FILE_NAME_RESULT_CSV,
        DATASET_FILE_NAME_RESULT_NEXTSTRAIN_TSV,
        DATASET_FILE_NAME_RESULT_NEXTSTRAIN_CSV,
        DATASET_FILE_NAME_RESULT_all_consensus,
        DATASET_FILE_NAME_MAFFT,
        DATASET_FILE_NAME_FASTTREE,
        DATASET_FILE_NAME_FASTTREE_tree,
        DATASET_FILE_NAME_nex,
        DATASET_FILE_NAME_nextstrain_auspice_zip,
        DATASET_FILE_NAME_nextstrain_error 
    ]
    #] + VECT_files_next_strain + [DATASET_FILE_NAME_nextstrain_error]
        
    DATASET_FILE_NAME_all_files_zipped = "AllFiles.zip"                    ### Several files zipped
    
    ## put the type file here to clean if there isn't enough sequences to create the trees and alignments
    vect_clean_file = [DATASET_FILE_NAME_MAFFT, DATASET_FILE_NAME_FASTTREE,\
                    DATASET_FILE_NAME_FASTTREE_tree,\
                    DATASET_FILE_NAME_nex]
    
    name = models.CharField(max_length=200, db_index=True, blank=True, null=True, verbose_name='Dataset name')
    owner = models.ForeignKey(User, related_name='dataset', blank=True, null=True, on_delete=models.CASCADE)
    creation_date = models.DateTimeField('Uploaded date', auto_now_add=True)
    last_change_date = models.DateTimeField('Last change date', blank=True, null=True)
    is_deleted = models.BooleanField(default=False)
    number_of_sequences_from_projects = models.SmallIntegerField(default=0)     ### has the number of sequences from projects
    number_of_sequences_from_consensus = models.SmallIntegerField(default=0)    ### has the number of sequences from consensus
    number_of_sequences_from_references = models.SmallIntegerField(default=0)   ### has the number of sequences from references
    is_processed = models.BooleanField(default=False)                           ### if some sequence was added and not processed yet
    
    number_passed_sequences = models.SmallIntegerField(default=0)     ### has the number of consensus that pass in join consensus
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
        for dataset_consensus in self.dataset_consensus.all():
            if dataset_consensus.is_deleted or dataset_consensus.is_error: continue
            if not dataset_consensus.consensus is None: continue
            if not dataset_consensus.reference is None and not dataset_consensus.reference.pk in dt_id_ref:
                dt_id_ref[dataset_consensus.reference.pk] = 1
            if not dataset_consensus.project_sample is None and not dataset_consensus.project_sample.project.reference.pk in dt_id_ref:
                dt_id_ref[dataset_consensus.project_sample.project.reference.pk] = 1
        return len(dt_id_ref)

    def _clean_name(self, name_to_clean, dict_to_clean = { ' ' : '_', '(' : '', ')' : '', '$' : '', '#' : '', '&' : '', '/' : '', '\\' : '', '-' : '_' }):
        """
        clean a name based on dictionary, dict_to_clean = { ' ' : '_', '(' : '' , ')' : '' }
        
        """
        for key in dict_to_clean:
            name_to_clean = name_to_clean.replace(key, dict_to_clean[key])
        return name_to_clean
    
    def get_clean_dataset_name(self):
        return self._clean_name(self.name)
    
    def get_global_file_by_dataset(self, type_path, file_name):
        """
        type_path: constants.TypePath -> MEDIA_ROOT, MEDIA_URL
        file_name:    Project.DATASET_FILE_NAME_MAFFT, ....  
        """
        return os.path.join(self.__get_global_path__(type_path, None), file_name)

    def get_global_file_by_dataset_web(self, file_name):
        
        out_file = self.get_global_file_by_dataset(TypePath.MEDIA_ROOT, file_name)
        if (os.path.exists(out_file)):
            return mark_safe('<a href="{}" download="{}"> {}</a>'.format(self.get_global_file_by_dataset(\
                        TypePath.MEDIA_URL, file_name), file_name, file_name))
        return 'File not available yet.'
        
        
    def __get_user_result_global_directory_path__(self, element = None):
        # file will be uploaded to MEDIA_ROOT/<filename>
        if (not element is None and len(element) > 0): return Constants.DIR_PROCESSED_FILES_DATASETS + \
            '/user_{0}/dataset_{1}/{2}/{3}'.\
            format(self.owner.id, self.pk, self.PATH_MAIN_RESULT, element)
        return Constants.DIR_PROCESSED_FILES_DATASETS + '/user_{0}/dataset_{1}/{2}'.format(
            self.owner.id, self.pk, self.PATH_MAIN_RESULT)
    
    def __get_global_path__(self, type_path, element = None):
        """
        get a path, from MEDIA_URL or MEDIA_ROOT
        """
        path_to_find = self.__get_user_result_global_directory_path__(element)
        if (type_path == TypePath.MEDIA_ROOT): 
            if not path_to_find.startswith('/'): path_to_find = os.path.join(getattr(settings, "MEDIA_ROOT", None), path_to_find)
        ### URL
        else: path_to_find = os.path.join(getattr(settings, "MEDIA_URL", None), path_to_find)
        return path_to_find
    
    def get_first_reference_name(self):
        """
        return first reference name 
        """
        for dataset_consensus in self.dataset_consensus.all():
            if dataset_consensus.is_deleted or dataset_consensus.is_error: continue
            if not dataset_consensus.reference is None: return dataset_consensus.reference.name
        return ""

    def get_first_reference(self):
        """
        return first reference 
        """
        for dataset_consensus in self.dataset_consensus.all():
            if dataset_consensus.is_deleted or dataset_consensus.is_error: continue
            if not dataset_consensus.reference is None: return dataset_consensus.reference
        return None

class DatasetConsensus(models.Model):
    
    ## Name from sample, reference or consensus, to improve the search by name 
    name = models.CharField(max_length=200, db_index=True, blank=True, null=True, verbose_name='Name')
    ## from sample only
    type_subtype = models.CharField(max_length=50, blank=True, null=True)    ## has the type/subtype collected
    
    dataset = models.ForeignKey(Dataset, related_name='dataset_consensus', blank=True, null=True, on_delete=models.CASCADE)
    ## only can have one of this three
    project_sample = models.ForeignKey(ProjectSample, related_name='dataset_project_samples', blank=True, null=True, on_delete=models.CASCADE)
    reference = models.ForeignKey(Reference, related_name='dataset_reference', blank=True, null=True, on_delete=models.CASCADE)
    consensus = models.ForeignKey(Consensus, related_name='dataset_consensus', blank=True, null=True, on_delete=models.CASCADE)
    is_project_sample_finished = models.BooleanField(default=False)    ## True if all process in ProjectSample Over
    ##
    creation_date = models.DateTimeField('uploaded date', auto_now_add=True)
    is_error = models.BooleanField(default=False)           ## if some problem occurs
    alert_first_level = models.IntegerField(default=0)      ## has the number of alerts for high errors
    alert_second_level = models.IntegerField(default=0)     ## has the number of alerts for low errors
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
        
    def get_name(self):
        if not self.project_sample is None: return self.project_sample.sample.name
        elif not self.consensus is None: return self.consensus.name
        elif not self.reference is None: return self.reference.name
        return ""

    def is_ready_to_proccess(self):
        if not self.project_sample is None: return self.project_sample.get_is_ready_to_proccess()
        return True
    
    def get_consensus_file(self, type_data):
        if not self.project_sample is None: return self.project_sample.get_consensus_file(type_data)
        elif not self.consensus is None: return self.consensus.get_consensus_fasta(type_data)
        elif not self.reference is None: return self.reference.get_reference_fasta(type_data)
        return ""
       
    def get_project_name(self):
        if not self.project_sample is None: return self.project_sample.project.name
        return ""
        
class MetaKey(models.Model):
    """
    Has meta tags to put values, for example, quality in the files, or samples
    """
    name = models.CharField(max_length=200, db_index=True, blank=True, null=True)
    def __str__(self):
        return self.name
    
    class Meta:
        ordering = ['name', ]


class UploadFiles(models.Model):
    """
    this class has the files that the user can upload, has he want,
    then the system make the relations with the samples
    """
    is_valid = models.BooleanField(default=False, verbose_name="Is valid")          ## true if everything is OK with the file, without errors
    is_processed = models.BooleanField(default=False, verbose_name='Is processed')  ## True if the New metadata is set to consensus
                                                                            ## if fastq.gz -> True when the file is attributed
    is_deleted = models.BooleanField(default=False)                         ## if this file is removed
    number_errors = models.IntegerField(default=0)          ## if has errors don't do anything, need to remove and upload again.
    
    type_file = models.ForeignKey(MetaKey, related_name='upload_metadata_files', blank=True, null=True, on_delete=models.CASCADE)    ## has the type of file
                                                            ## constants.TYPE_FILE.TYPE_FILE_dataset_file_metadata
    file_name = models.CharField(max_length=300, blank=True, null=True)     ## in fastq file, must have the same name in samples_list file
    creation_date = models.DateTimeField(auto_now_add=True, verbose_name='Uploaded Date')
    
    ### if is deleted in file system
    is_deleted_in_file_system = models.BooleanField(default=False)            ## if this file was removed in file system
    date_deleted = models.DateTimeField(blank=True, null=True, verbose_name='Date attached') ## this date has the time of deleted by web page
    
    owner = models.ForeignKey(User, related_name='upload_metadata_files', on_delete=models.CASCADE)
    
    ### need to create a random name for this file
    path_name = ContentTypeRestrictedFileField(upload_to=user_directory_path, blank=True, null=True,\
                    content_types=['application/octet-stream', 'application/gzip', 'application/x-gzip', 'text/csv', 'text/txt', 'text/tsv',\
                                'application/vnd.ms-excel', 'text/tab-separated-values', 'text/plain'],\
                    max_upload_size=settings.MAX_FASTQ_FILE_WITH_DOWNSIZE if settings.DOWN_SIZE_FASTQ_FILES else settings.MAX_FASTQ_FILE_UPLOAD,\
                    max_length=500)

    dataset = models.ForeignKey(Dataset, related_name='dataset_upload_metadata_file', blank=True, null=True, on_delete=models.CASCADE)
    description = models.TextField(default="")                ## has a json result.ProcessResults instance with errors or successes
    
    ## constants
    constants = Constants()
    
    class Meta:
        ordering = ['-creation_date']

    def __str__(self):
        return self.file_name
    
    def get_path_to_file(self, type_path):
        """
        get a path, type_path, from MEDIA_URL or MEDIA_ROOT
        """
        path_to_find = self.path_name.name
        if (type_path == TypePath.MEDIA_ROOT): 
            if not path_to_find.startswith('/'): path_to_find = os.path.join(getattr(settings, "MEDIA_ROOT", None), path_to_find)
        else:
            path_to_find = os.path.join(getattr(settings, "MEDIA_URL", None), path_to_find)
        return path_to_find

    def get_metadata_fasta_web(self):
        """
        return web link for metadata
        """
        out_file = self.get_path_to_file(TypePath.MEDIA_ROOT)
        if (os.path.exists(out_file)):
            return mark_safe('<a href="{}" download="{}"> {}</a>'.format(self.get_path_to_file(\
                        TypePath.MEDIA_URL), os.path.basename(self.get_path_to_file(TypePath.MEDIA_ROOT)),
                        self.constants.short_name(self.file_name, Constants.SHORT_NAME_LENGTH)))
        return 'File not available.'

class MetaKeyDataset(models.Model):
    """
    Relation ManyToMany in 
    """
    meta_tag = models.ForeignKey(MetaKey, related_name='meta_key_dataset', on_delete=models.CASCADE)
    dataset = models.ForeignKey(Dataset, related_name='meta_key_dataset', on_delete=models.CASCADE)
    owner = models.ForeignKey(User, related_name='meta_key_dataset', on_delete=models.CASCADE)
    creation_date = models.DateTimeField('uploaded date', db_index=True, auto_now_add=True)
    value = models.CharField(default=Constants.META_KEY_VALUE_NOT_NEED, max_length=200)
    description = models.TextField(default="")
    
    class Meta:
        ordering = ['dataset__id', '-creation_date']
    
    def __str__(self):
        return self.meta_tag.name + " " + self.value + " " + self.description


