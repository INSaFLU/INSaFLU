from django.db import models
# Create your models here.


class Tags(models.Model):
	name = models.CharField(max_length=50, blank=True, null=True)
	def __str__(self):
		return self.name
	
	class Meta:
		ordering = ['name', ]

class UploadFile(models.Model):
	name = models.CharField(max_length=100, blank=True, null=True)
	path = models.CharField(max_length=500, blank=True, null=True)
	version = models.IntegerField(blank=True, null=True)
	creation_date = models.DateTimeField('uploaded date', auto_now_add=True)
	abricate_name = models.CharField(max_length=100, blank=True, null=True)
	def __str__(self):
		return self.name
	
	class Meta:
		ordering = ['version', ]

class SeqVirus(models.Model):
	"""
	Has the virus identification
	The file must have this name <name>_<V1>.fasta, to get the version of the file
	"""
	accession = models.CharField(max_length=50, blank=True, null=True)
	name = models.CharField(max_length=100, blank=True, null=True)
	kind_type = models.ForeignKey(Tags, related_name='seq_virus', blank=True, null=True, on_delete=models.CASCADE)
	file = models.ForeignKey(UploadFile, related_name='seq_virus', blank=True, null=True, on_delete=models.CASCADE)
	description = models.CharField(max_length=500, blank=True, null=True)
	def __str__(self):
		return self.name
	
	class Meta:
		ordering = ['accession', ]

class IdentifyVirus(models.Model):
	"""
	identify the virus
	"""
	seq_virus = models.ForeignKey(SeqVirus, related_name='identify_virus', blank=True, null=True, on_delete=models.CASCADE)
	rank = models.IntegerField(default=0)	## has the rank, first the type, then the subType
	coverage = models.CharField(max_length=10, blank=True, null=True)
	identity = models.CharField(max_length=10, blank=True, null=True)
	
	def __str__(self):
		return "name:{}     type:{}    rank:{}".format(self.seq_virus.name, self.seq_virus.kind_type, self.rank)
	
	def __eq__(self, other):
		return self.seq_virus.name == other.seq_virus.name and self.seq_virus.kind_type == other.seq_virus.kind_type



