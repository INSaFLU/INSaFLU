from django.db import models

# Create your models here.
from django.contrib.gis.db.models import PointField
from django.contrib.gis.db.models import GeoManager
from django.contrib.auth.models import User
from constants.constants import Constants, TypePath
from fluwebvirus.formatChecker import ContentTypeRestrictedFileField
from manage_virus.models import IdentifyVirus
from django.conf import settings
import os

def reference_directory_path(instance, filename):
	# file will be uploaded to MEDIA_ROOT/<filename>
	return 'uploads/generic_data/user_{0}/{1}'.format(instance.owner.id, filename)


def user_directory_path(instance, filename):
	# file will be uploaded to MEDIA_ROOT/<filename>
	return 'uploads/generic_data/user_{0}/{1}'.format(instance.owner.id, filename)

class SeasonReference(models.Model):
	"""
	Each sample needs a dataset 
	"""
	name = models.CharField(max_length=100, db_index=True, blank=True, null=True)
	owner = models.ForeignKey(User, related_name='season_reference', blank=True, null=True, on_delete=models.CASCADE)
	creation_date = models.DateTimeField(auto_now_add=True, verbose_name='Date creation')
	def __str__(self):
		return self.name
	class Meta:
		ordering = ['name', ]

	
class MetaKey(models.Model):
	"""
	Has meta tags to put values, for example, quality in the files, or samples
	"""
	name = models.CharField(max_length=200, db_index=True, blank=True, null=True)
	def __str__(self):
		return self.name
	
	class Meta:
		ordering = ['name', ]
		
class Reference(models.Model):
	name = models.CharField(max_length=200, db_index=True, verbose_name='Reference name')
	isolate_name = models.CharField(max_length=200, default='', verbose_name='Isolate Name')
	creation_date = models.DateTimeField(auto_now_add=True, verbose_name='Uploaded Date')
	
	## Size 100K
	reference_fasta = ContentTypeRestrictedFileField(upload_to=reference_directory_path, content_types=['application/octet-stream', 'application/gzip'],\
										max_upload_size=Constants.MAX_REF_FASTA_FILE, blank=True, null=True, max_length=500)
	reference_fasta_name = models.CharField(max_length=200, default='', verbose_name='Fasta file')
	hash_reference_fasta = models.CharField(max_length=50, blank=True, null=True)

	## Size 200K
	reference_genbank = ContentTypeRestrictedFileField(upload_to=reference_directory_path, content_types=['application/octet-stream', 'application/gzip'],\
									max_upload_size=Constants.MAX_REF_GENBANK_FILE, blank=True, null=True, max_length=500)
	reference_genbank_name = models.CharField(max_length=200, default='', verbose_name='Genbank file')
	hash_reference_genbank = models.CharField(max_length=50, blank=True, null=True)

	owner = models.ForeignKey(User, related_name='reference', blank=True, null=True, on_delete=models.CASCADE)
	is_obsolete = models.BooleanField(default=False, verbose_name='Obsolete')
	number_of_locus = models.IntegerField(default=0, verbose_name='#Sequences')
	
	season = models.ManyToManyField(SeasonReference)		## can have the season
	description = models.CharField(max_length=500, default='', blank=True, null=True, verbose_name='Description')

	def __str__(self):
		return self.name

	def get_reference_gbk(self, type_path):
		"""
		get a path, type_path, from MEDIA_URL or MEDIA_ROOT
		"""
		path_to_find = self.reference_genbank.name
		if (type_path == TypePath.MEDIA_ROOT): 
			if not path_to_find.startswith('/'): path_to_find = os.path.join(getattr(settings, "MEDIA_ROOT", None), path_to_find)
		else:
			path_to_find = os.path.join(getattr(settings, "MEDIA_URL", None), path_to_find)
		return path_to_find
	
	def get_reference_fasta(self, type_path):
		"""
		get a path, type_path, from MEDIA_URL or MEDIA_ROOT
		"""
		path_to_find = self.reference_fasta.name
		if (type_path == TypePath.MEDIA_ROOT): 
			if not path_to_find.startswith('/'): path_to_find = os.path.join(getattr(settings, "MEDIA_ROOT", None), path_to_find)
		else:
			path_to_find = os.path.join(getattr(settings, "MEDIA_URL", None), path_to_find)
		return path_to_find

	class Meta:
		verbose_name = 'Reference'
		verbose_name_plural = 'References'
		ordering = ['-creation_date', ]
		indexes = [
			models.Index(fields=['name'], name='name_idx'),
		]

class MetaKeyReference(models.Model):
	"""
	Relation ManyToMany in 
	"""
	meta_tag = models.ForeignKey(MetaKey, related_name='meta_key_reference', on_delete=models.CASCADE)
	reference = models.ForeignKey(Reference, related_name='meta_key_reference', on_delete=models.CASCADE)
	owner = models.ForeignKey(User, related_name='meta_key_reference', on_delete=models.CASCADE)
	creation_date = models.DateTimeField('uploaded date', db_index=True, auto_now_add=True)
	value = models.CharField(default=Constants.META_KEY_VALUE_NOT_NEED, max_length=200)
	description = models.TextField(default="")
	
	class Meta:
		ordering = ['reference__id', '-creation_date']
	
	def __str__(self):
		return self.meta_tag.name + " " + self.value + " " + self.description



class TagName(models.Model):
	"""
	Has the tags to ticked the samples
	"""
	name = models.CharField(max_length=100, db_index=True, blank=True, null=True)
	owner = models.ForeignKey(User, related_name='tag_name', blank=True, null=True, on_delete=models.CASCADE)
	creation_date = models.DateTimeField(auto_now_add=True, verbose_name='Uploaded Date')
	is_meta_data = models.BooleanField(default=False)	## if this tag belongs to meta data or not.
	def __str__(self):
		return self.name
	class Meta:
		ordering = ['name', ]


class DataSet(models.Model):
	"""
	Each sample needs a dataset 
	"""
	name = models.CharField(max_length=100, db_index=True, blank=True, null=True)
	owner = models.ForeignKey(User, related_name='data_set', blank=True, null=True, on_delete=models.CASCADE)
	creation_date = models.DateTimeField(auto_now_add=True, verbose_name='Date creation')
	def __str__(self):
		return self.name
	class Meta:
		ordering = ['creation_date', 'name', ]
		
class VaccineStatus(models.Model):
	"""
	Each sample needs a dataset 
	"""
	name = models.CharField(max_length=100, db_index=True, blank=True, null=True)
	owner = models.ForeignKey(User, related_name='vaccine_status', blank=True, null=True, on_delete=models.CASCADE)
	creation_date = models.DateTimeField(auto_now_add=True, verbose_name='Date creation')
	
	def __str__(self):
		return self.name
	
	class Meta:
		ordering = ['creation_date', 'name', ]

	
class Sample(models.Model):
	"""
	Sample, each sample has one or two files...
	"""
	OUT_FILE_ABRICATE = "abricate.txt"
	
	## to remove in future
	objects = models.Manager()	## need to check this
	
	name = models.CharField(max_length=200, db_index=True, blank=True, null=True)  ## This Id should match the prefix of the reads files (i.e. prefix_R1_001.fastq.gz /  
																	##	prefix_R2_001.fastq.gz),
	date_of_onset = models.DateField('date of onset', blank=True, null=True)
	date_of_collection = models.DateField('date of collection', blank=True, null=True)
	date_of_receipt_lab = models.DateField('date of receipt lab', blank=True, null=True)
	week = models.IntegerField(blank=True, null=True)
	day = models.IntegerField(blank=True, null=True)	## from "Date of onset” or “Date of collection” 
	month = models.IntegerField(blank=True, null=True)
	year = models.IntegerField(blank=True, null=True)
	vaccine_status = models.ForeignKey(VaccineStatus, related_name='sample', blank=True, null=True, on_delete=models.CASCADE)
	creation_date = models.DateTimeField(auto_now_add=True, verbose_name='Uploaded Date')
	is_deleted = models.BooleanField(default=False)
	is_obsolete = models.BooleanField(default=False)
	owner = models.ForeignKey(User, related_name='sample', blank=True, null=True, on_delete=models.CASCADE)
	data_set = models.ForeignKey(DataSet, related_name='sample', blank=True, null=True, on_delete=models.CASCADE)
	geo_local = PointField(null=True, blank=True, srid=4326);  ## 4326 which means latitude and longitude
	geo_manager = GeoManager()
	identify_virus = models.ManyToManyField(IdentifyVirus)
	type_subtype = models.CharField(max_length=50, blank=True, null=True)	## has the type/subtype collected

	## many to many relation	
	tag_names = models.ManyToManyField(TagName, through='TagNames')
	
	### files
	is_valid_1 = models.BooleanField(default=False)
	file_name_1 = models.CharField(max_length=200, blank=True, null=True)
	candidate_file_name_1 = models.CharField(max_length=200, blank=True, null=True)
	
	## 30M
	path_name_1 = ContentTypeRestrictedFileField(upload_to=user_directory_path, blank=True, null=True,\
					content_types=['application/octet-stream', 'application/gzip'], max_upload_size=Constants.MAX_FASTQ_FILE, max_length=500)
	is_valid_2 = models.BooleanField(default=False)
	file_name_2 = models.CharField(max_length=200, blank=True, null=True)
	candidate_file_name_2 = models.CharField(max_length=200, blank=True, null=True)
	path_name_2 = ContentTypeRestrictedFileField(upload_to=user_directory_path, blank=True, null=True,\
					content_types=['application/octet-stream', 'application/gzip'], max_upload_size=Constants.MAX_FASTQ_FILE, max_length=500)

	## has files, the user can upload the files after
	has_files = models.BooleanField(default=False)
	
	###	has the flag indicating that the sample can be processed by projects
	is_ready_for_projects = models.BooleanField(default=False)
	
	def __str__(self):
		return self.name

	class Meta:
		ordering = ['-creation_date', ]

	def get_file_names(self):
		sz_return = "" if self.file_name_1 == None else self.file_name_1
		sz_return += "/" if len(sz_return) > 0 and self.file_name_2 is not None else ""
		sz_return += "" if self.file_name_2 == None else self.file_name_2
		return sz_return

	def exist_file_2(self):
		if (self.path_name_2 is None or self.path_name_2.name is None): return False
		return True

	def get_abricate_output(self, type_path):
		"""
		type_path = [MEDIA_ROOT, MEDIA_URL]
		Return the file name of the abricate output base on fastq File input
		path it's a FileField instance, or a string
		"""
		return os.path.join(self.__get_path__(type_path, True), self.OUT_FILE_ABRICATE)
	
	def get_trimmomatic_file(self, type_path, b_first_file):
		"""
		get the trimmomatic files, it's going to be use for all processing
		type_path = [MEDIA_ROOT, MEDIA_URL]
		"""
		if (not b_first_file and not self.exist_file_2()): return None
		path_out = os.path.join(self.__get_path__(type_path, b_first_file), Constants.DIR_PROCESSED_PROCESSED)
		if (type_path == TypePath.MEDIA_ROOT): os.makedirs(path_out, mode=0o755, exist_ok=True)
		return os.path.join(path_out, self.name + ("_1P.fastq.gz" if b_first_file else "_2P.fastq.gz"))

	def get_fastq_trimmomatic(self, type_path, b_first_file):
		"""
		return fastq output first step
		"""
		if (not b_first_file and not self.exist_file_2()): return None
		
		path_out = os.path.join(self.__get_path__(type_path, b_first_file), Constants.DIR_PROCESSED_PROCESSED)
		if (type_path == TypePath.MEDIA_ROOT): os.makedirs(path_out, mode=0o755, exist_ok=True)
		return os.path.join(path_out, self.name + ("_1P_fastqc.html" if b_first_file else "_2P_fastqc.html"))

	def get_fastq(self, type_path, b_first_file):
		"""
		return fastq output first step
		"""
		if (not b_first_file and not self.exist_file_2()): return None
		return os.path.join(self.__get_path__(type_path, b_first_file), self.file_name_1 if b_first_file else self.file_name_2)

	def get_fastq_output(self, type_path, b_first_file):
		"""
		return fastq output first step
		"""
		if (not b_first_file and not self.exist_file_2()): return None
		return os.path.join(self.__get_path__(type_path, b_first_file), self.file_name_1.replace(".fastq.gz", "_fastqc.html") if b_first_file else self.file_name_2.replace(".fastq.gz", "_fastqc.html"))

	def __get_path__(self, type_path, b_first_file):
		"""
		get a path, from MEDIA_URL or MEDIA_ROOT
		"""
		if (b_first_file):
			path_to_find = os.path.dirname(self.path_name_1.name)
			if (type_path == TypePath.MEDIA_ROOT): 
				if not path_to_find.startswith('/'): path_to_find = os.path.join(getattr(settings, "MEDIA_ROOT", None), path_to_find)
			else:
				path_to_find = os.path.join(getattr(settings, "MEDIA_URL", None), path_to_find)
		else:
			path_to_find = os.path.dirname(self.path_name_2.name)
			if (type_path == TypePath.MEDIA_ROOT): 
				if not path_to_find.startswith('/'): path_to_find = os.path.join(getattr(settings, "MEDIA_ROOT", None), path_to_find)
			else:
				path_to_find = os.path.join(getattr(settings, "MEDIA_URL", None), path_to_find)
		return path_to_find

	def get_type_sub_type(self):
		vect_identify_virus = self.identify_virus.all()
		sz_return = ""
		if (vect_identify_virus.count() > 0):
			sz_return = self.__get_type__(vect_identify_virus, Constants.SEQ_VIRUS_TYPE)
			sz_type = self.__get_type__(vect_identify_virus, Constants.SEQ_VIRUS_SUB_TYPE)
			if (len(sz_type) > 0): sz_return += "" if len(sz_return) == 0 else "-" + sz_type
			sz_type = self.__get_type__(vect_identify_virus, Constants.SEQ_VIRUS_LINEAGE)
			if (len(sz_type) > 0): sz_return += "" if len(sz_return) == 0 else "-" + sz_type
			return sz_return
		else: return Constants.EMPTY_VALUE_TABLE
		
	def __get_type__(self, vect_identify_virus, type_to_test):
		vect_return = []
		for identify_virus in vect_identify_virus:
			if (identify_virus.seq_virus.kind_type.name == type_to_test):
				vect_return.append(identify_virus.seq_virus.name)
		return ''.join(sorted(vect_return))
				
	def get_is_ready_for_projects(self):
		"""
		need to be true to be ready for projects
		"""
		return self.is_ready_for_projects and not self.is_obsolete and not self.is_deleted

class TagNames(models.Model):
	value = models.CharField(max_length=150)
	tag_name = models.ForeignKey(TagName, on_delete=models.CASCADE)
	sample = models.ForeignKey(Sample, on_delete=models.CASCADE)

	def __str__(self):
		return self.value
	
class MetaKeySample(models.Model):
	"""
	Relation ManyToMany in 
	"""
	meta_tag = models.ForeignKey(MetaKey, related_name='meta_key_sample', on_delete=models.CASCADE)
	sample = models.ForeignKey(Sample, related_name='meta_key_sample', on_delete=models.CASCADE)
	owner = models.ForeignKey(User, related_name='meta_key_sample', on_delete=models.CASCADE)
	creation_date = models.DateTimeField('uploaded date', db_index=True, auto_now_add=True)
	value = models.CharField(default=Constants.META_KEY_VALUE_NOT_NEED, max_length=200)
	description = models.TextField(default="")
	
	class Meta:
		ordering = ['sample__id', '-creation_date']
	
	def __str__(self):
		return self.meta_tag.name + " " + self.value + " " + self.description
	
class Project(models.Model):
	"""
	Project for run a pipeline in the data,
	It's possible to run several times and add data to the project
	"""
	### has the path for main result
	PATH_MAIN_RESULT = 'main_result'
	
	PROJECT_FILE_NAME_MAFFT = "Alignment_whole_genome.fasta"
	PROJECT_FILE_NAME_FASTTREE = "Tree_ML_WG.nwk"
	PROJECT_FILE_NAME_FASTTREE_tree = "Tree_ML_WG.tree"
	PROJECT_FILE_NAME_nex = "Alignment_whole_genome.nex"
	
	## put the type file here to clean if there isn't enough sequences to create the trees and alignments
	vect_clean_file = [PROJECT_FILE_NAME_MAFFT, PROJECT_FILE_NAME_FASTTREE,\
					PROJECT_FILE_NAME_FASTTREE_tree, PROJECT_FILE_NAME_nex]

	## obsolete
	PROJECT_FILE_NAME_GRAPH_MINO_VAR_HTML = "graph_minor_var.html"
	PROJECT_FILE_NAME_GRAPH_MINO_VAR_PNG = "graph_minor_var.png"
	## end obsolete
	
	### this is only to join with other names
	PROJECT_FILE_NAME_FASTTREE_element = "Tree"
	PROJECT_FILE_NAME_MAFFT_element = "Alignment"
	
	name = models.CharField(max_length=200, db_index=True, blank=True, null=True)
	owner = models.ForeignKey(User, related_name='project', blank=True, null=True, on_delete=models.CASCADE)
	reference = models.ForeignKey(Reference, related_name='project', blank=True, null=True, on_delete=models.CASCADE)
	creation_date = models.DateTimeField('uploaded date', auto_now_add=True)
	is_deleted = models.BooleanField(default=False)
	
	def __str__(self):
		return self.name
	
	class Meta:
		ordering = ['-creation_date', ]

	def get_global_file_by_element(self, type_path, element, file_name):
		"""
		type_path: constants.TypePath -> MEDIA_ROOT, MEDIA_URL
		element: element name or None
		file_name: Project.PROJECT_FILE_NAME_MAFFT, ....   
		"""
		if (self.PROJECT_FILE_NAME_MAFFT == file_name):
			return os.path.join(self.__get_global_path__(type_path, element), "{}_{}.fasta".format(self.PROJECT_FILE_NAME_MAFFT_element, element))
		if (self.PROJECT_FILE_NAME_FASTTREE == file_name):
			return os.path.join(self.__get_global_path__(type_path, element), "{}_{}.nwk".format(self.PROJECT_FILE_NAME_FASTTREE_element, element))
		if (self.PROJECT_FILE_NAME_FASTTREE_tree == file_name):
			return os.path.join(self.__get_global_path__(type_path, element), "{}_{}.tree".format(self.PROJECT_FILE_NAME_FASTTREE_element, element))
		if (self.PROJECT_FILE_NAME_nex == file_name):
			return os.path.join(self.__get_global_path__(type_path, element), "{}_{}.nex".format(self.PROJECT_FILE_NAME_MAFFT_element, element))
		return None

	def get_global_file_by_element_and_cds(self, type_path, element, CDS, file_name):
		"""
		get file names for proteins
		type_path: constants.TypePath -> MEDIA_ROOT, MEDIA_URL
		element: element name or None
		file_name: Project.PROJECT_FILE_NAME_MAFFT, ....   
		"""
		if (self.PROJECT_FILE_NAME_MAFFT == file_name):		## protein alignement
			return os.path.join(self.__get_global_path__(type_path, element), "{}_{}_{}.faa".format(self.PROJECT_FILE_NAME_MAFFT_element, element, CDS))
		if (self.PROJECT_FILE_NAME_FASTTREE == file_name):
			return os.path.join(self.__get_global_path__(type_path, element), "{}_{}_{}.nwk".format(self.PROJECT_FILE_NAME_FASTTREE_element, element, CDS))
		if (self.PROJECT_FILE_NAME_FASTTREE_tree == file_name):
			return os.path.join(self.__get_global_path__(type_path, element), "{}_{}_{}.tree".format(self.PROJECT_FILE_NAME_FASTTREE_element, element, CDS))
		if (self.PROJECT_FILE_NAME_nex == file_name):
			return os.path.join(self.__get_global_path__(type_path, element), "{}_{}_{}.nex".format(self.PROJECT_FILE_NAME_MAFFT_element, element, CDS))
		return None

	def get_global_file_by_project(self, type_path, file_name):
		"""
		type_path: constants.TypePath -> MEDIA_ROOT, MEDIA_URL
		file_name:	Project.PROJECT_FILE_NAME_MAFFT, ....  
		"""
		return os.path.join(self.__get_global_path__(type_path, None), file_name)
	
	def __get_user_result_global_directory_path__(self, element):
		# file will be uploaded to MEDIA_ROOT/<filename>
		if (element != None and len(element) > 0): return 'projects/result/user_{0}/project_{1}/{2}/{3}'.\
				format(self.owner.id, self.pk, self.PATH_MAIN_RESULT, element)
		return 'projects/result/user_{0}/project_{1}/{2}'.format(self.owner.id, self.pk, self.PATH_MAIN_RESULT)
	
	def __get_global_path__(self, type_path, element):
		"""
		get a path, from MEDIA_URL or MEDIA_ROOT
		"""
		path_to_find = self.__get_user_result_global_directory_path__(element)
		if (type_path == TypePath.MEDIA_ROOT): 
			if not path_to_find.startswith('/'): path_to_find = os.path.join(getattr(settings, "MEDIA_ROOT", None), path_to_find)
		else: path_to_find = os.path.join(getattr(settings, "MEDIA_URL", None), path_to_find)
		return path_to_find
	
class MetaKeyProject(models.Model):
	"""
	Relation ManyToMany in 
	"""
	meta_tag = models.ForeignKey(MetaKey, related_name='meta_key_project', on_delete=models.CASCADE)
	project = models.ForeignKey(Project, related_name='meta_key_project', on_delete=models.CASCADE)
	owner = models.ForeignKey(User, related_name='meta_key_project', on_delete=models.CASCADE)
	creation_date = models.DateTimeField('uploaded date', db_index=True, auto_now_add=True)
	value = models.CharField(default=Constants.META_KEY_VALUE_NOT_NEED, max_length=200)
	description = models.TextField(default="")
	
	class Meta:
		ordering = ['project__id', '-creation_date']
	
	def __str__(self):
		return self.meta_tag.name + " " + self.value + " " + self.description

class MixedInfectionsTag(models.Model):
	"""
	Used to tag mixed infections
	"""
	name = models.CharField(max_length=50, db_index=True, blank=True, null=True)
	def __str__(self):
		return self.name
	
	class Meta:
		ordering = ['name', ]

class MixedInfections(models.Model):
	"""
	Used to identify mixed infections
	"""
	tag = models.ForeignKey(MixedInfectionsTag, related_name='mixed_infections', blank=True, null=True, on_delete=models.CASCADE)
	average_value = models.FloatField(default=0.0)
	description = models.TextField(default="")
	creation_date = models.DateTimeField('uploaded date', auto_now_add=True)
	has_master_vector = models.BooleanField(default=False)  ## if it has the master vector, has the vector to compare to all others
															## and is not used in projectSample,
															## It can change across time
															## to trace the change of tag is set a metaValue in ProjectSample
	class Meta:
		ordering = ['tag', ]

class CountVariations(models.Model):
	"""
	has the number of variations for a sample
	"""
	total = models.PositiveIntegerField(default=0)
	var_less_50 = models.PositiveIntegerField(default=0)
	var_bigger_50_90 = models.PositiveIntegerField(default=0)
	var_bigger_90 = models.PositiveIntegerField(default=0)
	
	def __str__(self):
		return 'Total: {} Less 50:{} 50<Value<90 :{}  Bigger >90:{}'.format(self.total, self.var_less_50, self.var_bigger_50_90, self.var_bigger_90)


class ProjectSample(models.Model):
	
	PATH_MAIN_RESULT = 'main_result'
	PREFIX_FILE_COVERAGE = 'coverage'
		
	project = models.ForeignKey(Project, related_name='project_sample', blank=True, null=True, on_delete=models.CASCADE)
	sample = models.ForeignKey(Sample, related_name='project_sample', blank=True, null=True, on_delete=models.CASCADE)
	mixed_infections = models.ForeignKey(MixedInfections, related_name='project_sample', blank=True, null=True, on_delete=models.CASCADE)
	count_variations = models.ForeignKey(CountVariations, related_name='project_sample', blank=True, null=True, on_delete=models.CASCADE)
	creation_date = models.DateTimeField('uploaded date', auto_now_add=True)
	is_finished = models.BooleanField(default=False)
	is_deleted = models.BooleanField(default=False)
	is_error = models.BooleanField(default=False)		## if some problem occurs
	alert_first_level = models.IntegerField(default=0)	## has the number of alerts for high errors
	alert_second_level = models.IntegerField(default=0)	## has the number of alerts for low errors
	
	class Meta:
		ordering = ['project__id', '-creation_date']
	
	def __str__(self):
		return self.project.name

	def get_global_file_by_element(self, type_path, prefix_file_name, sequence_name, extension):
		"""
		type_path: constants.TypePath -> MEDIA_ROOT, MEDIA_URL
		prefix_file_name: ProjectSample.PREFIX_FILE_COVERAGE
		sequence_name: sequence name  
		extension: FileExtensions.FILE_PNG
		"""
		return os.path.join(self.__get_path__(type_path, ProjectSample.PATH_MAIN_RESULT), "{}_{}{}".format(prefix_file_name, sequence_name, extension))
	
	def get_file_output(self, type_path, file_type, software):
		"""
		return file path output by software
		type_path: constants.TypePath -> MEDIA_ROOT, MEDIA_URL
		file_type: constants.FileType -> FILE_BAM, FILE_BAM_BAI, FILE_CONSENSUS_FA, ...
		software: Software.SOFTWARE_FREEBAYES_name, Software.SOFTWARE_SNIPPY_name
		"""
		constants = Constants()
		return os.path.join(self.__get_path__(type_path, software.lower() if software != None else None), constants.get_extensions_by_file_type(self.sample.name, file_type))

	def __get_user_result_directory_path__(self, software):
		# file will be uploaded to MEDIA_ROOT/<filename>
		if (software is None or not software): return 'projects/result/user_{0}/project_{1}/sample_{2}'.format(self.project.owner.id, self.project.id, self.sample.id)
		return 'projects/result/user_{0}/project_{1}/sample_{2}/{3}'.format(self.project.owner.id, self.project.id, self.sample.id, software)


	def __get_path__(self, type_path, software):
		"""
		get a path, from MEDIA_URL or MEDIA_ROOT
		"""
		path_to_find = self.__get_user_result_directory_path__(software)
		if (type_path == TypePath.MEDIA_ROOT): 
			if not path_to_find.startswith('/'): path_to_find = os.path.join(getattr(settings, "MEDIA_ROOT", None), path_to_find)
		else: path_to_find = os.path.join(getattr(settings, "MEDIA_URL", None), path_to_find)
		return path_to_find
	

	def get_is_ready_to_proccess(self):
		"""
		test if is ready to process
		"""
		return self.is_finished and not self.is_deleted and not self.is_error and self.sample.get_is_ready_for_projects()

	
	
class MetaKeyProjectSample(models.Model):
	"""
	Relation ManyToMany in 
	"""
	meta_tag = models.ForeignKey(MetaKey, related_name='meta_key_project_sample', on_delete=models.CASCADE)
	project_sample = models.ForeignKey(ProjectSample, related_name='meta_key_project_sample', on_delete=models.CASCADE)
	owner = models.ForeignKey(User, related_name='meta_key_project_sample', on_delete=models.CASCADE)
	creation_date = models.DateTimeField('uploaded date', db_index=True, auto_now_add=True)
	value = models.CharField(default=Constants.META_KEY_VALUE_NOT_NEED, max_length=200)
	description = models.TextField(default="")
	
	class Meta:
		ordering = ['project_sample__id', '-creation_date']
	
	def __str__(self):
		return self.meta_tag.name + " " + self.value + " " + self.description

class Statistics(models.Model):
	"""
	Has several percentils about the table CountVariations 
	"""
	tag = models.ForeignKey(TagName, related_name='statistics', on_delete=models.CASCADE)
	value = models.FloatField(default=0.0)
	
	def __str__(self):
		return 'Tag: {} Value: {}'.format(self.tag.name, self.value)

# class CosineDistance(models.Model):
# 	"""
# 	Has the values of cosine distance 
# 	"""
# 	pass

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


class UploadFiles(models.Model):
	"""
	this class has the files that the user can upload, has he want,
	then the system make the relations with the samples
	"""
	is_valid = models.BooleanField(default=False)			## true if everything is OK with the file, without errors
	is_processed = models.BooleanField(default=False)		## if samples file -> True when all files is attributed
															## if fastq.gz -> True when the file is attributed
	is_deleted = models.BooleanField(default=False)			## if this file is removed
	number_errors = models.IntegerField(default=0)			## if has errors don't do anything, need to remove and upload again.
	number_files_processed = models.IntegerField(default=0)	## samples_list, has the number of files already processed
															## in fastq files, if this file is associated to a sample or not
	number_files_to_process = models.IntegerField(default=0)	## samples_list, has the number of files to process. At the end this number must be equal
	
	type_file = models.ForeignKey(MetaKey, related_name='upload_files', blank=True, null=True, on_delete=models.CASCADE)	## has the type of file
															## constants.TYPE_FILE.TYPE_FILE_fastq_gz
															## constants.TYPE_FILE.TYPE_FILE_sample_file
	file_name = models.CharField(max_length=300, blank=True, null=True)	## in fastq file, must have the same name in samples_list file
	creation_date = models.DateTimeField('upload_files', auto_now_add=True)
	owner = models.ForeignKey(User, related_name='upload_files', on_delete=models.CASCADE)
	
	### need to create a random name for this file
	path_name = ContentTypeRestrictedFileField(upload_to=user_directory_path, blank=True, null=True,\
					content_types=['application/octet-stream', 'application/gzip', 'text/csv', 'text/tsv'], max_upload_size=Constants.MAX_FASTQ_FILE, max_length=500)
#	path_name = models.FileField(upload_to=user_directory_path, blank=True, null=True)

	samples = models.ManyToManyField(Sample) 	## if fastq file has the sample where it belongs
												## if samples_file has all the relations with samples. Must be all created, files attributed, or deleted
												##   to add other samples file
	upload_file = models.ForeignKey('self', blank=True, null=True) 			## in fastq file has the sample list where it belongs
	description = models.TextField(default="")				## has a json result.ProcessResults instance with errors or successes
	
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


