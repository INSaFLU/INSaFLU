from django.db import models

# Create your models here.
from django.contrib.gis.db.models import PointField
from django.contrib.gis.db.models import GeoManager
from django.contrib.auth.models import User
from utils.constants import Constants
from fluwebvirus.formatChecker import ContentTypeRestrictedFileField
from manage_virus.models import IdentifyVirus

def reference_directory_path(instance, filename):
	# file will be uploaded to MEDIA_ROOT/<filename>
	return 'uploads/generic_data/user_{0}/{1}'.format(instance.owner.id, filename)


def user_directory_path(instance, filename):
	# file will be uploaded to MEDIA_ROOT/<filename>
	return 'uploads/generic_data/user_{0}/{1}'.format(instance.owner.id, filename)

class Reference(models.Model):
	name = models.CharField(max_length=200, default='New reference')
	scientific_name = models.CharField(max_length=200, default='')
	creation_date = models.DateTimeField(auto_now_add=True, verbose_name='Up.Date')
	## Size 100K
	reference_fasta = ContentTypeRestrictedFileField(upload_to=reference_directory_path, content_types=['application/octet-stream'], max_upload_size=100000, blank=True, null=True)
	reference_fasta_name = models.CharField(max_length=200, default='', verbose_name='Fasta file')
	hash_reference_fasta = models.CharField(max_length=50, blank=True, null=True)

	## Size 200K
	reference_genbank = ContentTypeRestrictedFileField(upload_to=reference_directory_path, content_types=['application/octet-stream'], max_upload_size=200000, blank=True, null=True)
	reference_genbank_name = models.CharField(max_length=200, default='', verbose_name='Genbank file')
	hash_reference_genbank = models.CharField(max_length=50, blank=True, null=True)

	owner = models.ForeignKey(User, related_name='reference', blank=True, null=True, on_delete=models.CASCADE)
	is_obsolete = models.BooleanField(default=False, verbose_name='Obsolete')
	number_of_locus = models.IntegerField(default=0, verbose_name='#Segments')

	def __str__(self):
		return self.file_name

	class Meta:
		verbose_name = 'Reference'
		verbose_name_plural = 'References'
		ordering = ['-creation_date', ]
		indexes = [
			models.Index(fields=['name'], name='name_idx'),
		]



class MetaKey(models.Model):
	"""
	Has meta tags to put values, for example, quality in the files, or samples
	"""
	name = models.CharField(max_length=200, blank=True, null=True)
	def __str__(self):
		return self.name
	class Meta:
		ordering = ['name', ]


class TagName(models.Model):
	"""
	Has the tags to ticked the samples
	"""
	name = models.CharField(max_length=100, blank=True, null=True)
	owner = models.ForeignKey(User, related_name='tag_name', blank=True, null=True, on_delete=models.CASCADE)
	def __str__(self):
		return self.name
	class Meta:
		ordering = ['name', ]


class DataSet(models.Model):
	"""
	Each sample needs a dataset 
	"""
	name = models.CharField(max_length=100, blank=True, null=True)
	owner = models.ForeignKey(User, related_name='data_set', blank=True, null=True, on_delete=models.CASCADE)
	def __str__(self):
		return self.name
	class Meta:
		ordering = ['name', ]
		
			
class Sample(models.Model):
	"""
	Sample, each sample has one or two files...
	"""
	name = models.CharField(max_length=200, blank=True, null=True)
	sample_date = models.DateField('sample date', blank=True, null=True)
	creation_date = models.DateTimeField('uploaded date', auto_now_add=True)
	is_rejected = models.BooleanField(default=False)
	is_obsolete = models.BooleanField(default=False)
	owner = models.ForeignKey(User, related_name='sample', blank=True, null=True, on_delete=models.CASCADE)
	tag_names = models.ManyToManyField(TagName)
	data_set = models.ForeignKey(DataSet, related_name='sample', blank=True, null=True, on_delete=models.CASCADE)
	geo_local = PointField(null=True, blank=True, srid=4326);  ## 4326 which means latitude and longitude
	objects = GeoManager()
	identify_virus = models.ManyToManyField(IdentifyVirus)

	### files
	is_valid_1 = models.BooleanField(default=True)
	file_name_1 = models.CharField(max_length=300, blank=True, null=True)
	path_name_1 = ContentTypeRestrictedFileField(upload_to=reference_directory_path, content_types=['application/octet-stream'], max_upload_size=30971520, blank=True, null=True)
	is_valid_2 = models.BooleanField(default=True)
	file_name_2 = ContentTypeRestrictedFileField(upload_to=reference_directory_path, content_types=['application/octet-stream'], max_upload_size=30971520, blank=True, null=True)
	path_name_2 = models.FileField(upload_to=user_directory_path, blank=True, null=True, max_length=400)

	def __str__(self):
		return self.name

	class Meta:
		ordering = ['-creation_date', ]

	def get_file_names(self):
		sz_return = "" if self.file_name_1 == None else self.file_name_1
		sz_return += "/" if len(sz_return) > 0 and self.file_name_2 is not None else ""
		sz_return += "" if self.file_name_2 == None else self.file_name_2
		return sz_return

	def exit_file_2(self):
		if (self.path_name_2 is None): return False
		return True

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


