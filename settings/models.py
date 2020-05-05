from django.db import models
from django.contrib.auth.models import User
# Create your models here.


class Software(models.Model):
	"""
	Each user has it software parameters 
	"""
	name = models.CharField(max_length=100, db_index=True, blank=True, null=True)
	version = models.CharField(max_length=100, db_index=True, blank=True, null=True)
	owner = models.ForeignKey(User, related_name='software_settings', blank=True, null=True, on_delete=models.PROTECT)

	def __str__(self):
		return self.name

	class Meta:
		ordering = ['name', ]
		

class Parameter(models.Model):
	"""
	Can have more than one, for example
	
	SOFTWARE_TRIMMOMATIC_PARAMETERS = "SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33"
	"""
	
	PARAMETER_int = 0 
	PARAMETER_float = 1 
	PARAMETER_char = 2
	PARAMETER_null = 3		### its like
	
	name = models.CharField(max_length=50, db_index=True, blank=True, null=True)
	parameter = models.CharField(max_length=50, db_index=True, blank=True, null=True)
	type_data = models.SmallIntegerField()
	software = models.ForeignKey(Software, related_name='parameter', on_delete=models.PROTECT)
	union_char = models.CharField(max_length=10, default="")	### ":" in the case LEADING:3
	can_change = models.BooleanField(default=True)	### TOPHRED33 can change, and it is not to show in change dialog
	sequence_out = models.SmallIntegerField()
	range_available = models.CharField(max_length=50, default="")		## description of the range
	
	range_max = models.CharField(max_length=50, default="")				## only used in int and float fields
	range_min = models.CharField(max_length=50, default="")				## only used in int and float fields
	description = models.CharField(max_length=500, default="")			## description of this size
	not_set_value = models.CharField(max_length=50, db_index=True, blank=True, null=True)  ## dont set value if not this value
	
	def __str__(self):
		return self.name
	
	class Meta:
		ordering = ['sequence_out', ]
		
	def is_integer(self):
		return self.type_data == Parameter.PARAMETER_int
	def is_float(self):
		return self.type_data == Parameter.PARAMETER_float
	def is_char(self):
		return self.type_data == Parameter.PARAMETER_char
	
	def get_unique_id(self):
		"""
		get unique ID for forms
		"""
		return "{}_{}".format(self.name, self.sequence_out)




