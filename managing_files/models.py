from django.db import models

# Create your models here.
from django.contrib.gis.db.models import PointField
from django.contrib.auth.models import User
from constants.Constants import Constants


def reference_directory_path(instance, filename):
    # file will be uploaded to MEDIA_ROOT/<filename>
    return 'uploads/reference/user_{0}/{1}'.format(instance.user.id, filename)


def user_directory_path(instance, filename):
    # file will be uploaded to MEDIA_ROOT/<filename>
    return 'uploads/raw_files/user_{0}/%Y/%m/%d/{1}'.format(instance.user.id, filename)


class Reference(models.Model):
    name = models.CharField(max_length=200, default='New reference')
    scientific_name = models.CharField(max_length=200, default='')
    creation_date = models.DateTimeField('uploaded date', auto_now_add=True)
    reference_fasta = models.FileField(upload_to=reference_directory_path, max_length=600, blank=True, null=True)
    hash_reference_fasta = models.CharField(max_length=50, blank=True, null=True)
    reference_genbank = models.FileField(upload_to=reference_directory_path, max_length=600, blank=True, null=True)
    hash_reference_genbank = models.CharField(max_length=50, blank=True, null=True)
    
    owner = models.ForeignKey(User, related_name='reference', blank=True, null=True, on_delete=models.CASCADE)
    is_obsolete = models.BooleanField(default=False)
    number_of_locus = models.IntegerField(default=0)
         
    def __str__(self):
        return self.file_name

    class Meta:
        verbose_name = 'Reference'
        verbose_name_plural = 'References'
        ordering = ['-creation_date', ]
        indexes = [
            models.Index(fields=['name'], name='name_idx'),
        ]


class Files(models.Model):
    """
    Files with raw data, only fastq.gz
    """
    is_valid_1 = models.BooleanField(default=True)
    file_name_1 = models.CharField(max_length=300, blank=True, null=True)
    path_name_1 = models.FileField(upload_to=user_directory_path, max_length=400)
    is_valid_2 = models.BooleanField(default=True)
    file_name_2 = models.CharField(max_length=300, blank=True, null=True)
    path_name_2 = models.FileField(upload_to=user_directory_path, max_length=400)
    owner = models.ForeignKey(User, related_name='files', blank=True, null=True, on_delete=models.CASCADE)
    
    def __str__(self):
        return self.file_name_1 + " " + self.file_name_2


class SampleTag(models.Model):
    """
    Has a simple tag to put for the samples of one year
    """
    name = models.CharField(max_length=200, blank=True, null=True)
    owner = models.ForeignKey(User, related_name='sample_tag', blank=True, null=True, on_delete=models.CASCADE)
    
    def __str__(self):
        return self.name
    
    class Meta:
        ordering = ['name', ]
    

class MetaKey(models.Model):
    """
    Has meta tags to put values, for example, quality in the files, or samples
    """
    name = models.CharField(max_length=200, blank=True, null=True)
    def __str__(self):
        return self.name
    class Meta:
        ordering = ['name', ]

    
class Sample(models.Model):
    """
    Sample, each sample has one or two files...
    """
    name = models.CharField(max_length=200, blank=True, null=True)
    date_sample = models.DateField('sample date', default='New sample')
    creation_date = models.DateTimeField('uploaded date', auto_now_add=True)
    is_rejected = models.BooleanField(default=False)
    owner = models.ForeignKey(User, related_name='sample', blank=True, null=True, on_delete=models.CASCADE)
    files = models.ForeignKey(Files, related_name='sample', on_delete=models.CASCADE)
    tag_name = models.ForeignKey(SampleTag, related_name='sample', blank=True, null=True, on_delete=models.CASCADE)
    geo_local = PointField(null=True, blank=True,);
        
    def __str__(self):
        return self.name

    class Meta:
        ordering = ['-creation_date', ]


class MetaKeySample(models.Model):
    """
    Relation ManyToMany in 
    """
    meta_tag = models.ForeignKey(MetaKey, related_name='meta_key_sample', on_delete=models.CASCADE)
    sample = models.ForeignKey(Sample, related_name='meta_key_sample', on_delete=models.CASCADE)
    owner = models.ForeignKey(User, related_name='meta_key_sample', on_delete=models.CASCADE)
    creation_date = models.DateTimeField('uploaded date', auto_now_add=True)
    value = models.CharField(default=Constants.META_KEY_VALUE_NOT_NEED, max_length=200)
    description = models.TextField(default="")
    
    class Meta:
        ordering = ['sample__id', 'meta_tag__id', 'creation_date']
    

class Project(models.Model):
    """
    Project for run a pipeline in the data,
    It's possible to run several times and add data to the project
    """
    name = models.CharField(max_length=200, blank=True, null=True)
    owner = models.ForeignKey(User, related_name='project', blank=True, null=True, on_delete=models.CASCADE)
    samples = models.ManyToManyField(Sample, related_name='project')
    creation_date = models.DateTimeField('uploaded date', auto_now_add=True)
    is_finished = models.BooleanField(default=False)
    
    def __str__(self):
        return self.name
    
    class Meta:
        ordering = ['-creation_date', ]


class Version(models.Model):
    """
    """
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name
    
    class Meta:
        ordering = ['name', ]
        

class Software(models.Model):
    """
    Has all the softwares
    It can be used for all users 
    """
    name = models.CharField(max_length=100, blank=True, null=True)
    path_to_run = models.CharField(max_length=300)
    version = models.ForeignKey(Version, related_name='software', on_delete=models.CASCADE)
    
    class Meta:
        ordering = ['name', 'version__name']





