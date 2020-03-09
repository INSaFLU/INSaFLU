from django.db import models

# Create your models here.
from django.contrib.gis.db.models import PointField
from django.contrib.gis.db.models import GeoManager    ##  change to django  2.x
#from django.db.models import Manager as GeoManager
from django.contrib.auth.models import User
from constants.constants import Constants, TypePath, FileExtensions, FileType
from fluwebvirus.formatChecker import ContentTypeRestrictedFileField
from manage_virus.models import IdentifyVirus
from django.conf import settings
from django.utils.translation import ugettext_lazy as _
from django.utils.safestring import mark_safe
from operator import itemgetter
import os
from constants.software_names import SoftwareNames
from manage_virus.constants_virus import ConstantsVirus
from constants.constants_mixed_infection import ConstantsMixedInfection

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
	display_name = models.CharField(max_length=200, db_index=True, default='', verbose_name='Display name')
	isolate_name = models.CharField(max_length=200, default='', verbose_name='Isolate Name')
	creation_date = models.DateTimeField(auto_now_add=True, verbose_name='Uploaded Date')
	
	## Size 100K
	reference_fasta = ContentTypeRestrictedFileField(upload_to=reference_directory_path, content_types=['application/octet-stream'],\
										max_upload_size=Constants.MAX_REF_FASTA_FILE, blank=True, null=True, max_length=500)
	reference_fasta_name = models.CharField(max_length=200, default='', verbose_name='Fasta file')
	hash_reference_fasta = models.CharField(max_length=50, blank=True, null=True)

	## Size 200K
	## application/x-gameboy-rom because of 'gb' extension file of gbk
	reference_genbank = ContentTypeRestrictedFileField(upload_to=reference_directory_path, content_types=['application/octet-stream', 'application/x-gameboy-rom'],\
									max_upload_size=Constants.MAX_REF_GENBANK_FILE, blank=True, null=True, max_length=500)
	reference_genbank_name = models.CharField(max_length=200, default='', verbose_name='Genbank file')
	hash_reference_genbank = models.CharField(max_length=50, blank=True, null=True)

	owner = models.ForeignKey(User, related_name='reference', blank=True, null=True, on_delete=models.CASCADE)
	is_obsolete = models.BooleanField(default=False, verbose_name='Obsolete')
	is_deleted = models.BooleanField(default=False, verbose_name='Deleted')
	number_of_locus = models.IntegerField(default=0, verbose_name='#Locus')
	
	season = models.ManyToManyField(SeasonReference)		## can have the season
	description = models.CharField(max_length=500, default='', blank=True, null=True, verbose_name='Description')

	### if is deleted in file system
	is_deleted_in_file_system = models.BooleanField(default=False)			## if this file was removed in file system
	date_deleted = models.DateTimeField(blank=True, null=True, verbose_name='Date attached') ## this date has the time of deleted by web page
	
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
	
	def get_reference_fasta_index(self, type_path):
		"""
		get a fasta fai path, type_path, from MEDIA_URL or MEDIA_ROOT
		"""
		return self.get_reference_fasta(type_path) + FileExtensions.FILE_FAI

	def get_reference_bed(self, type_path):
		"""
		get a fasta bed path, type_path, from MEDIA_URL or MEDIA_ROOT
		"""
		file_name = self.get_reference_gbk(type_path)
		return file_name[:file_name.rfind('.')] + FileExtensions.FILE_BED
	
	def get_reference_bed_index(self, type_path):
		"""
		get a fasta bed.idx path, type_path, from MEDIA_URL or MEDIA_ROOT
		"""
		return self.get_reference_bed(type_path) + FileExtensions.FILE_IDX


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
	last_change_date = models.DateTimeField('uploaded date', blank=True, null=True)
	has_master_vector = models.BooleanField(default=False)  ## if it has the master vector, has the vector to compare to all others
															## and is not used in projectSample,
															## It can change across time
															## to trace the change of tag is set a metaValue in ProjectSample
	class Meta:
		ordering = ['tag', ]


class Sample(models.Model):
	"""
	Sample, each sample has one or two files...
	"""
	OUT_FILE_ABRICATE = "abricate.txt"
	DRAFT_CONTIGS = "_draft_contigs.fasta"
	SEGMENTS_CONTIGS = "_segments_to_contigs.tsv"
	
	### consensus file
	OUT_CONSENSUS_FILE = "consensus.fasta"
	
	## to remove in future
	objects = models.Manager()	## need to check this
	
	name = models.CharField(max_length=200, db_index=True, blank=True, null=True, verbose_name='Sample Name')  ## This Id should match the prefix of the reads files (i.e. prefix_R1_001.fastq.gz /  
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
	number_alerts = models.IntegerField(verbose_name='Alerts', default=0, blank=True, null=True)	## has the number of alerts
	mixed_infections_tag = models.ForeignKey(MixedInfectionsTag, verbose_name='Putative Mixed Infection', related_name='sample', blank=True, null=True, on_delete=models.CASCADE)
																			### has the tag of Yes/No mixed infection

	## many to many relation	
	tag_names = models.ManyToManyField(TagName, through='TagNames')

	### files
	is_valid_1 = models.BooleanField(default=False)
	file_name_1 = models.CharField(max_length=200, blank=True, null=True)
	candidate_file_name_1 = models.CharField(max_length=200, blank=True, null=True)
	
	## 30M
	path_name_1 = ContentTypeRestrictedFileField(upload_to=user_directory_path, blank=True, null=True,\
					content_types=['application/octet-stream', 'application/gzip', 'application/x-gzip'],\
					max_upload_size=settings.MAX_FASTQ_FILE_WITH_DOWNSIZE if settings.DOWN_SIZE_FASTQ_FILES else settings.MAX_FASTQ_FILE_UPLOAD,\
					max_length=500)
	is_valid_2 = models.BooleanField(default=False)
	file_name_2 = models.CharField(max_length=200, blank=True, null=True)
	candidate_file_name_2 = models.CharField(max_length=200, blank=True, null=True)
	path_name_2 = ContentTypeRestrictedFileField(upload_to=user_directory_path, blank=True, null=True,\
					content_types=['application/octet-stream', 'application/gzip', 'application/x-gzip'],\
					max_upload_size=settings.MAX_FASTQ_FILE_WITH_DOWNSIZE if settings.DOWN_SIZE_FASTQ_FILES else settings.MAX_FASTQ_FILE_UPLOAD,\
					max_length=500)

	## has files, the user can upload the files after
	has_files = models.BooleanField(default=False)
	
	###	has the flag indicating that the sample can be processed by projects
	is_ready_for_projects = models.BooleanField(default=False)
	
	### if is deleted in file system
	is_deleted_in_file_system = models.BooleanField(default=False)			## if this file was removed in file system
	date_deleted = models.DateTimeField(blank=True, null=True, verbose_name='Date attached')	## this date has the time of deleted by web page
	
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
		if (self.path_name_2 is None or self.path_name_2.name is None or len(self.path_name_2.name) == 0): return False
		return True

	def get_abricate_output(self, type_path):
		"""
		type_path = [MEDIA_ROOT, MEDIA_URL]
		Return the file name of the abricate output base on fastq File input
		path it's a FileField instance, or a string
		"""
		return os.path.join(self.__get_path__(type_path, True), self.OUT_FILE_ABRICATE)
	
	def get_draft_contigs_output(self, type_path):
		"""
		type_path = [MEDIA_ROOT, MEDIA_URL]
		Return the file name of the abricate output base on fastq File input
		path it's a FileField instance, or a string
		"""
		return os.path.join(self.__get_path__(type_path, True), "{}{}".format(self.name.replace(' ', '_'), self.DRAFT_CONTIGS))
	
	def get_draft_contigs_abricate_output(self, type_path):
		"""
		type_path = [MEDIA_ROOT, MEDIA_URL]
		Return the file name of the abricate output base on fastq File input
		path it's a FileField instance, or a string
		"""
		return os.path.join(self.__get_path__(type_path, True), "{}{}".format(self.name.replace(' ', '_'), self.SEGMENTS_CONTIGS))
	
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
		return fastq output first step, from MEDIA_URL or MEDIA_ROOT
		"""
		if (not b_first_file and not self.exist_file_2()): return None
		return os.path.join(self.__get_path__(type_path, b_first_file), self.file_name_1 if b_first_file else self.file_name_2)

	def is_original_fastq_removed(self):
		"""
		Test if original fastq were removed already
		"""
		return not os.path.exists(self.get_fastq(TypePath.MEDIA_ROOT, True))
	
	def get_fastq_output(self, type_path, b_first_file):
		"""
		return fastq output second step
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
		"""
		return a string with Type/Subtpye
		"""
		vect_identify_virus_temp = self.identify_virus.all()
		
		### clean some possible repeated
		vect_identify_virus = []
		for identify_virus in vect_identify_virus_temp:
			if (not identify_virus in vect_identify_virus):
				vect_identify_virus.append(identify_virus)

		if (len(vect_identify_virus) > 0):

			### Corona
			sz_return_c = self.__get_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_GENUS, [ConstantsVirus.TYPE_BetaCoV])
			sz_subtype = self.__get_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_HUMAN, [])
			if (len(sz_subtype) > 0): sz_return_c += sz_subtype if len(sz_return_c) == 0 else "-" + sz_subtype
						
			### Type A
			## the second flag is for corona
			sz_return_a = self.__get_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_TYPE, [ConstantsVirus.TYPE_A, \
												ConstantsVirus.TYPE_BetaCoV])
			sz_subtype = self.__get_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SUB_TYPE, [ConstantsVirus.TYPE_A])
			if (len(sz_subtype) > 0): sz_return_a += sz_subtype if len(sz_return_a) == 0 else "-" + sz_subtype
			
			### Type B
			sz_return_b = self.__get_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_TYPE, [ConstantsVirus.TYPE_B])
			sz_subtype = self.__get_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_LINEAGE, [ConstantsVirus.TYPE_B])
			if (len(sz_subtype) > 0): sz_return_b += sz_subtype if len(sz_return_b) == 0 else "-" + sz_subtype
			
			sz_return = sz_return_c
			if (len(sz_return_a) > 0): sz_return = sz_return_a if len(sz_return) == 0 else "{}; {}".format(sz_return, sz_return_a)
			if (len(sz_return_b) > 0): sz_return = sz_return_b if len(sz_return) == 0 else "{}; {}".format(sz_return, sz_return_b)
			return sz_return
		else: return Constants.EMPTY_VALUE_TYPE_SUBTYPE
		
	def __get_type__(self, vect_identify_virus, type_to_test, vect_virus_name):
		vect_return = []
		for identify_virus in vect_identify_virus:
			if (identify_virus.seq_virus.kind_type.name == type_to_test and (identify_virus.seq_virus.name in vect_virus_name \
							or len(vect_virus_name) == 0)):
				vect_return.append(identify_virus.seq_virus.name)
			elif (type_to_test == ConstantsVirus.SEQ_VIRUS_SUB_TYPE and identify_virus.seq_virus.kind_type.name == type_to_test and ConstantsVirus.TYPE_A in vect_virus_name):
				vect_return.append(identify_virus.seq_virus.name)
			elif (type_to_test == ConstantsVirus.SEQ_VIRUS_LINEAGE and identify_virus.seq_virus.kind_type.name == type_to_test and ConstantsVirus.TYPE_B in vect_virus_name):
				vect_return.append(identify_virus.seq_virus.name)
		if (type_to_test == ConstantsVirus.SEQ_VIRUS_SUB_TYPE and len(vect_return) > 2): return '|'.join(sorted(vect_return))
		if (type_to_test == ConstantsVirus.SEQ_VIRUS_HUMAN and len(vect_return) > 1): return '|'.join(sorted(vect_return))
		return ''.join(sorted(vect_return))
	
	def __exists_type(self, vect_identify_virus, virus_name, type_virus = ConstantsVirus.SEQ_VIRUS_TYPE):
		"""
		test if exist some specific type
		"""
		for identify_virus in vect_identify_virus:
			if (identify_virus.seq_virus.kind_type.name == type_virus and \
				( identify_virus.seq_virus.name == virus_name or len(virus_name) == 0)): return True
		return False
				
	def __get_number_type__(self, vect_identify_virus, type_to_test):
		"""
		get a number for a specific type
		"""
		n_return = 0
		for identify_virus in vect_identify_virus:
			if (identify_virus.seq_virus.kind_type.name == type_to_test): n_return += 1
		return n_return
	
	def __get_number_type_and_start_sub_type(self, vect_identify_virus, type_to_test, starts_sub_type):
		"""
		get a number for a specific type and name starts with a specific sub_type
		if (starts_sub_type == None) get all of a specific 
		"""
		n_return = 0
		for identify_virus in vect_identify_virus:
			if (identify_virus.seq_virus.kind_type.name == type_to_test and \
					(starts_sub_type == None or identify_virus.seq_virus.name.startswith(starts_sub_type))): 
				n_return += 1
		return n_return

	def get_mixed_infection(self):
		"""
		mixed infection based on the table static/mixed_infections/mixed_infections.xls
		return tuble (tag_mixed_infection, alert, message)
		tag_mixed_infection: ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO or ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES
		alert: positive number, or zero
		message: None, empty or the string with the error message
		"""
		vect_identify_virus_temp = self.identify_virus.all()
		if (vect_identify_virus_temp.count() == 0): return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, 1,\
					"Warning: no typing data was obtained (possible reason: low number of influenza reads).")
	
		vect_identify_virus = []
		for identify_virus in vect_identify_virus_temp:
			if (not identify_virus in vect_identify_virus):
				vect_identify_virus.append(identify_virus)
		
		## Only A not B
		if (self.__exists_type(vect_identify_virus, ConstantsVirus.TYPE_A) and not self.__exists_type(vect_identify_virus, ConstantsVirus.TYPE_B)):
			#  A; #any subtype; > 0 lineage
			if (self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_LINEAGE) > 0):
				return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, 1,\
					"Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.")
				
			## A; = 2 subtype
			if (self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SUB_TYPE) == 2): return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, 0, None)
			## A; > 2 subtype
			if (self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SUB_TYPE) > 2): return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, 1,\
					"Warning: more than two subtypes were detected for this sample, suggesting that may represent a 'mixed infection'.")
			## A < 2 subtype
			if (self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SUB_TYPE) < 2): return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, 1,\
					"Warning: an incomplete subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).")
		
		## Only B not A
		if (not self.__exists_type(vect_identify_virus, ConstantsVirus.TYPE_A) and self.__exists_type(vect_identify_virus, ConstantsVirus.TYPE_B)):
			#  B; #any subtype; > 0 lineage
			if (self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SUB_TYPE) > 0):
				return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, 1,\
					"Warning: more than one type/subtypes were detected for this sample, suggesting that may represent a 'mixed infection'.")
			
			## B; == 1 lineage
			if (self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_LINEAGE) == 1): return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, 0, None)
			## B; > 1 lineage
			if (self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_LINEAGE) > 1): return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, 1,\
						"Warning: more than one lineage were detected for this sample, suggesting that may represent a 'mixed infection'.")
			## B; < 1 lineage
			if (self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_LINEAGE) < 1): return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, 1,\
						"Warning: an incomplete lineage has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).")
			
		## A and B
		if (self.__exists_type(vect_identify_virus, ConstantsVirus.TYPE_A) and self.__exists_type(vect_identify_virus, ConstantsVirus.TYPE_B)):
			
			if (self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_LINEAGE) == 0 and\
					self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SUB_TYPE) == 0):
				return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, 1,\
						"Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.")
			return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, 1,\
						"Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.")
		
		## BEGIN corona
		if (self.__exists_type(vect_identify_virus, "", ConstantsVirus.SEQ_VIRUS_GENUS) or self.__exists_type(vect_identify_virus, "", ConstantsVirus.SEQ_VIRUS_HUMAN)):
			
			if (self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_GENUS) == 1):
				if (self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_HUMAN) == 1):
					return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, 0, None)
				if (self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_HUMAN) > 1):
					return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, 1,\
						"Warning: more than one human BetaCoV virus are likely present in this sample, suggesting that may represent a 'mixed infection'")
				if (self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_HUMAN) == 0):
					return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, 1,\
						"Warning: an incomplete human BetaCoV identification has been obtained (possible reasons: low number of  human BetaCoV reads, mixed infection, etc)")
			else:
				if (self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_HUMAN) == 1):
					return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, 1,\
						"Warning: an incomplete human BetaCoV identification has been obtained (possible reasons: low number of  human BetaCoV reads, mixed infection, etc)")
				else:
					return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, 1,\
						"Warning: more than one human BetaCoV virus are likely present in this sample, suggesting that may represent a 'mixed infection'")
		### END corona
					
		## not A and not B
		if (not self.__exists_type(vect_identify_virus, ConstantsVirus.TYPE_A) and not self.__exists_type(vect_identify_virus, ConstantsVirus.TYPE_B)):
			if ((self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_LINEAGE) == 1 and\
					self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SUB_TYPE) == 0) or\
					(self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SUB_TYPE) == 1 and\
					self.__get_number_type__(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_LINEAGE) == 0)):
				return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, 1,\
						"Warning: an incomplete type/subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).")
			else:
				count_type_N = self.__get_number_type_and_start_sub_type(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SUB_TYPE, ConstantsVirus.SUB_TYPE_STARTS_N)
				count_type_H = self.__get_number_type_and_start_sub_type(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SUB_TYPE, ConstantsVirus.SUB_TYPE_STARTS_H)
				count_type_other = self.__get_number_type_and_start_sub_type(vect_identify_virus, ConstantsVirus.SEQ_VIRUS_LINEAGE, None)
				
				if (count_type_N == 1 and count_type_H == 1 and count_type_other == 0):
					return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, 1,\
						"Warning: an incomplete type/subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).")
				elif (count_type_N > 1 or count_type_H > 1 or count_type_other > 0):
					return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, 1,\
						"Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.")
			return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES, 1,\
						"Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.")
		
		## default	
		return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, 0, None)
	
	def get_is_ready_for_projects(self):
		"""
		need to be true to be ready for projects
		"""
		return self.is_ready_for_projects and not self.is_obsolete and not self.is_deleted

	def get_vect_tag_names_and_value(self):
		"""
		return [[header1, value1], [header2, value2], [header3, value3], ...]
		"""
		vect_out = []
		query_set = TagNames.objects.filter(sample=self)
		for tag_names in query_set:
			vect_out.append([tag_names.tag_name.name, tag_names.value])
		return sorted(vect_out, key=itemgetter(1))

	def get_tag_names(self):
		"""
		get the tag names grouped by a number
		"""
		query_set = TagNames.objects.filter(sample=self)
		if (query_set.count() == 0): return None
		return query_set


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
	
	PROJECT_FILE_NAME_MAFFT = "Alignment_nt_All.fasta"
	PROJECT_FILE_NAME_FASTA = "All_nt.fasta"
	PROJECT_FILE_NAME_FASTTREE = "Tree_ML_All.nwk"
	PROJECT_FILE_NAME_FASTTREE_tree = "Tree_ML_All.tree"
	PROJECT_FILE_NAME_nex = "Alignment_nt_All.nex"
	PROJECT_FILE_NAME_COVERAGE = "coverage.tsv"
	PROJECT_FILE_NAME_TOTAL_VARIATIONS = "proportions_iSNVs_graph.tsv"
	PROJECT_FILE_NAME_TAB_VARIATIONS_SNIPPY = "validated_variants.tsv" 
	PROJECT_FILE_NAME_TAB_VARIATIONS_FREEBAYES = "validated_minor_iSNVs.tsv" 	## remove del and ins and everything bigger than >50
	PROJECT_FILE_NAME_SAMPLE_RESULT_TSV = "Sample_list.tsv" 	### first column ID instead of 'sample name' to be compatible with Phandango e Microreact
	PROJECT_FILE_NAME_SAMPLE_RESULT_CSV = "Sample_list.csv" 	### first column ID instead of 'sample name' to be compatible with Phandango e Microreact
	PROJECT_FILE_NAME_SAMPLE_RESULT_json = "Sample_list.json" 	### first column ID instead of 'sample name' to be compatible with Phandango e Microreact, to download to 
	
	## put the type file here to clean if there isn't enough sequences to create the trees and alignments
	vect_clean_file = [PROJECT_FILE_NAME_MAFFT, PROJECT_FILE_NAME_FASTTREE,\
					PROJECT_FILE_NAME_FASTTREE_tree,\
					PROJECT_FILE_NAME_nex, PROJECT_FILE_NAME_FASTA]
	
	vect_exclude_clean_file_from_proteins = [PROJECT_FILE_NAME_FASTA]

	## obsolete
	PROJECT_FILE_NAME_GRAPH_MINO_VAR_HTML = "graph_minor_var.html"
	PROJECT_FILE_NAME_GRAPH_MINO_VAR_PNG = "graph_minor_var.png"
	## end obsolete
	
	### this is only to join with other names
	PROJECT_FILE_NAME_FASTTREE_element = "Tree_ML"
	PROJECT_FILE_NAME_MAFFT_element_nt = "Alignment_nt"
	PROJECT_FILE_NAME_MAFFT_element_aa = "Alignment_aa"
	PROJECT_FILE_NAME_FASTA_element = "Sequences_nt"
	
	name = models.CharField(max_length=200, db_index=True, blank=True, null=True, verbose_name='Project name')
	owner = models.ForeignKey(User, related_name='project', blank=True, null=True, on_delete=models.CASCADE)
	reference = models.ForeignKey(Reference, related_name='project', blank=True, null=True, on_delete=models.CASCADE)
	creation_date = models.DateTimeField('uploaded date', auto_now_add=True)
	is_deleted = models.BooleanField(default=False)
	
	### if is deleted in file system
	is_deleted_in_file_system = models.BooleanField(default=False)			## if this file was removed in file system
	date_deleted = models.DateTimeField(blank=True, null=True, verbose_name='Date attached')	## this date has the time of deleted by web page
	
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
			return os.path.join(self.__get_global_path__(type_path, element), "{}_{}.fasta".format(self.PROJECT_FILE_NAME_MAFFT_element_nt, element))
		if (self.PROJECT_FILE_NAME_FASTA == file_name):
			return os.path.join(self.__get_global_path__(type_path, element), "{}_{}.fasta".format(self.PROJECT_FILE_NAME_FASTA_element, element))
		if (self.PROJECT_FILE_NAME_FASTTREE == file_name):
			return os.path.join(self.__get_global_path__(type_path, element), "{}_{}.nwk".format(self.PROJECT_FILE_NAME_FASTTREE_element, element))
		if (self.PROJECT_FILE_NAME_FASTTREE_tree == file_name):
			return os.path.join(self.__get_global_path__(type_path, element), "{}_{}.tree".format(self.PROJECT_FILE_NAME_FASTTREE_element, element))
		if (self.PROJECT_FILE_NAME_nex == file_name):
			return os.path.join(self.__get_global_path__(type_path, element), "{}_{}.nex".format(self.PROJECT_FILE_NAME_MAFFT_element_nt, element))
		return None

	def get_global_file_by_element_and_cds(self, type_path, element, CDS, file_name):
		"""
		get file names for proteins
		type_path: constants.TypePath -> MEDIA_ROOT, MEDIA_URL
		element: element name or None
		file_name: Project.PROJECT_FILE_NAME_MAFFT, ....   
		"""
		CDS = self._clean_name(CDS)
		if (self.PROJECT_FILE_NAME_MAFFT == file_name):		## protein alignement
			return os.path.join(self.__get_global_path__(type_path, element), "{}_{}_{}.fasta".format(self.PROJECT_FILE_NAME_MAFFT_element_aa, element, CDS))
		if (self.PROJECT_FILE_NAME_FASTTREE == file_name):
			return os.path.join(self.__get_global_path__(type_path, element), "{}_{}_{}.nwk".format(self.PROJECT_FILE_NAME_FASTTREE_element, element, CDS))
		if (self.PROJECT_FILE_NAME_FASTTREE_tree == file_name):
			return os.path.join(self.__get_global_path__(type_path, element), "{}_{}_{}.tree".format(self.PROJECT_FILE_NAME_FASTTREE_element, element, CDS))
		if (self.PROJECT_FILE_NAME_nex == file_name):
			return os.path.join(self.__get_global_path__(type_path, element), "{}_{}_{}.nex".format(self.PROJECT_FILE_NAME_MAFFT_element_aa, element, CDS))
		return None

	def _clean_name(self, name_to_clean, dict_to_clean = { ' ' : '_', '(' : '', ')' : '', '$' : '', '#' : '', '&' : '', '/' : '', '\\' : '', '-' : '_' }):
		"""
		clean a name based on dictionary, dict_to_clean = { ' ' : '_', '(' : '' , ')' : '' }
		
		"""
		for key in dict_to_clean:
			name_to_clean = name_to_clean.replace(key, dict_to_clean[key])
		return name_to_clean
	
	def get_global_file_by_project(self, type_path, file_name):
		"""
		type_path: constants.TypePath -> MEDIA_ROOT, MEDIA_URL
		file_name:	Project.PROJECT_FILE_NAME_MAFFT, ....  
		"""
		return os.path.join(self.__get_global_path__(type_path, None), file_name)

	def get_global_file_by_project_web(self, file_name):
		
		out_file = self.get_global_file_by_project(TypePath.MEDIA_ROOT, file_name)
		if (os.path.exists(out_file)):
			return mark_safe('<a href="{}" download> {}</a>'.format(self.get_global_file_by_project(\
						TypePath.MEDIA_URL, file_name), file_name))
		return _('File not available yet.')
		
		
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
	
	constants = Constants()
	
	PATH_MAIN_RESULT = 'main_result'
	PREFIX_FILE_COVERAGE = 'coverage'
	FILE_CONSENSUS_FILE = "Consensus_"
	FILE_SNIPPY_TAB = "validated_variants_sample_"
	FILE_FREEBAYES_TAB = "validated_minor_iSNVs_sample_"
	
	project = models.ForeignKey(Project, related_name='project_samples', blank=True, null=True, on_delete=models.CASCADE)
	sample = models.ForeignKey(Sample, related_name='project_samples', blank=True, null=True, on_delete=models.CASCADE)
	mixed_infections = models.ForeignKey(MixedInfections, related_name='project_samples', blank=True, null=True, on_delete=models.CASCADE)
	count_variations = models.ForeignKey(CountVariations, related_name='project_samples', blank=True, null=True, on_delete=models.CASCADE)
	creation_date = models.DateTimeField('uploaded date', auto_now_add=True)
	is_finished = models.BooleanField(default=False)
	is_deleted = models.BooleanField(default=False)
	is_error = models.BooleanField(default=False)		## if some problem occurs
	alert_first_level = models.IntegerField(default=0)	## has the number of alerts for high errors
	alert_second_level = models.IntegerField(default=0)	## has the number of alerts for low errors
	
	### if is deleted in file system
	is_deleted_in_file_system = models.BooleanField(default=False)			## if this file was removed in file system
	date_deleted = models.DateTimeField(blank=True, null=True, verbose_name='Date attached') ## this date has the time of deleted by web page

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
		software: SoftwareNames.SOFTWARE_FREEBAYES_name, SoftwareNames.SOFTWARE_SNIPPY_name
		"""
		constants = Constants()
		return os.path.join(self.__get_path__(type_path, software.lower() if not software is None else None), constants.get_extensions_by_file_type(self.sample.name, file_type))

	def get_file_output_human(self, type_path, file_type, software):
		"""
		return file path output by software, but with human name
		type_path: constants.TypePath -> MEDIA_ROOT, MEDIA_URL
		file_type: constants.FileType -> FILE_BAM, FILE_BAM_BAI, FILE_CONSENSUS_FA, ...
		software: SoftwareNames.SOFTWARE_FREEBAYES_name, SoftwareNames.SOFTWARE_SNIPPY_name
		"""
		constants = Constants()
		
		return os.path.join(self.__get_path__(type_path, software.lower() if software != None else None),\
				constants.get_extensions_by_file_type(self.__get_human_name_file__(software, file_type), file_type))

	def __get_human_name_file__(self, software, file_type):
		"""
		get human file name
		"""
		if (software == SoftwareNames.SOFTWARE_SNIPPY_name):
			if (file_type == FileType.FILE_TAB): return "{}{}".format(ProjectSample.FILE_SNIPPY_TAB, self.sample.name)
		if (software == SoftwareNames.SOFTWARE_FREEBAYES_name):
			if (file_type == FileType.FILE_TAB): return "{}{}".format(ProjectSample.FILE_FREEBAYES_TAB, self.sample.name)
		return self.sample.name


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

	def get_consensus_file(self, type_path):
		"""
		get clean consensus file name
		"""
		return os.path.join(self.__get_path__(type_path, ProjectSample.PATH_MAIN_RESULT), "{}{}{}".format(\
					ProjectSample.FILE_CONSENSUS_FILE, self.sample.name, FileExtensions.FILE_FASTA))
	
	def get_consensus_file_web(self):
		"""
		get consensus file web
		"""
		out_file = self.get_consensus_file(TypePath.MEDIA_ROOT)
		if (os.path.exists(out_file)):
			return mark_safe('<a href="{}" download> {}</a>'.format(self.get_consensus_file(\
						TypePath.MEDIA_URL), self.constants.short_name(os.path.basename(out_file), 15)))
		return _('Not available.')

	def get_file_web(self, file_type, software):
		"""
		get file web from different softwares
		type_path: constants.TypePath -> MEDIA_ROOT, MEDIA_URL
		file_type: constants.FileType -> FILE_BAM, FILE_BAM_BAI, FILE_CONSENSUS_FA, ...
		software: SoftwareNames.SOFTWARE_FREEBAYES_name, SoftwareNames.SOFTWARE_SNIPPY_name
		"""
		out_file = self.get_file_output_human(TypePath.MEDIA_ROOT, file_type, software)
		if (os.path.exists(out_file)):
			return mark_safe('<a href="{}" download> {}</a>'.format(\
				self.get_file_output_human(TypePath.MEDIA_URL, file_type, software),\
				self.constants.short_name(os.path.basename(out_file), 20)))
		return _('Not available.')
		
		
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
	is_valid = models.BooleanField(default=False)							## true if everything is OK with the file, without errors
	is_processed = models.BooleanField(default=False)						## if samples file -> True when all files is attributed
																			## if fastq.gz -> True when the file is attributed
	is_deleted = models.BooleanField(default=False)							## if this file is removed
	number_errors = models.IntegerField(default=0)			## if has errors don't do anything, need to remove and upload again.
	number_files_processed = models.IntegerField(default=0)	## samples_list, has the number of files already processed
															## in fastq files, if this file is associated to a sample or not
	number_files_to_process = models.IntegerField(default=0)	## samples_list, has the number of files to process. At the end this number must be equal to number_files_processed
	
	type_file = models.ForeignKey(MetaKey, related_name='upload_files', blank=True, null=True, on_delete=models.CASCADE)	## has the type of file
															## constants.TYPE_FILE.TYPE_FILE_fastq_gz
															## constants.TYPE_FILE.TYPE_FILE_sample_file
															## constants.TYPE_FILE.TYPE_FILE_sample_file_metadata
	file_name = models.CharField(max_length=300, blank=True, null=True)	## in fastq file, must have the same name in samples_list file
	creation_date = models.DateTimeField(auto_now_add=True, verbose_name='Uploaded Date')
	attached_date = models.DateTimeField(blank=True, null=True, verbose_name='Date attached')		## only used in fastq.gz files
	
	### if is deleted in file system
	is_deleted_in_file_system = models.BooleanField(default=False)			## if this file was removed in file system
	date_deleted = models.DateTimeField(blank=True, null=True, verbose_name='Date attached') ## this date has the time of deleted by web page
	
	owner = models.ForeignKey(User, related_name='upload_files', on_delete=models.CASCADE)
	
	### need to create a random name for this file
	path_name = ContentTypeRestrictedFileField(upload_to=user_directory_path, blank=True, null=True,\
					content_types=['application/octet-stream', 'application/gzip', 'application/x-gzip', 'text/csv', 'text/txt', 'text/tsv',\
								'application/vnd.ms-excel', 'text/tab-separated-values', 'text/plain'],\
					max_upload_size=settings.MAX_FASTQ_FILE_WITH_DOWNSIZE if settings.DOWN_SIZE_FASTQ_FILES else settings.MAX_FASTQ_FILE_UPLOAD,\
					max_length=500)

	samples = models.ManyToManyField(Sample) 	## if fastq file has the sample where it belongs
												## if samples_file has all the relations with samples. Must be all created, files attributed, or deleted
												##   to add other samples file
	upload_file = models.ForeignKey('self', blank=True, null=True, on_delete=models.CASCADE) 			## in fastq file has the sample list where it belongs
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


class ProcessControler(models.Model):
	"""
	this class has the info about the process
	"""
	PREFIX_SAMPLE = "sample_"
	PREFIX_PROJECT_SAMPLE = "project_sample_"
	PREFIX_PROJECT = "project_"
	PREFIX_UPLOAD_FILES = "upload_files_"
	PREFIX_LINK_FILES_USER = "link_files_user_"
	
	### flags of status
	FLAG_FINISHED = 'flag_finished'
	FLAG_RUNNING = 'flag_running'
	FLAG_ERROR = 'flag_error'
	
	owner = models.ForeignKey(User, related_name='process_controles', on_delete=models.CASCADE)
	creation_date = models.DateTimeField(auto_now_add=True, verbose_name='Submit job')
	close_date = models.DateTimeField(blank=True, null=True, verbose_name='Close date')
	is_finished = models.BooleanField(default=False)
	is_running = models.BooleanField(default=False)		## when is running
	is_error = models.BooleanField(default=False)
	name = models.CharField(max_length=50, db_index=True, blank=True, null=True)  ## has the name of the process
								## could be: sample_<pk>, project_sample_<pk>, project_<pk>
	name_sge_id = models.CharField(max_length=20, db_index=True, blank=True, null=True) ## has the SGE id about this process
	
	class Meta:
		ordering = ['creation_date']

	### get names for samples, project and project samples
	def get_name_sample(self, sample):
		return "{}{}".format(ProcessControler.PREFIX_SAMPLE, sample.pk)
	def get_name_project_sample(self, project_sample):
		return "{}{}".format(ProcessControler.PREFIX_PROJECT_SAMPLE, project_sample.pk)
	def get_name_project(self, project):
		return "{}{}".format(ProcessControler.PREFIX_PROJECT, project.pk)
	def get_name_upload_files(self, upload_files):
		return "{}{}".format(ProcessControler.PREFIX_UPLOAD_FILES, upload_files.pk)
	def get_name_link_files_user(self, user):
		return "{}{}".format(ProcessControler.PREFIX_LINK_FILES_USER, user.pk)

	def __str__(self):
		return "PK:{} name:{}  is_finished:{}  is_running:{}  is_error:{}".format(self.pk, self.name, self.is_finished, self.is_running, self.is_error)


