from django.db import models
from django.contrib.auth.models import User
# Create your models here.


## in INSaFLU database
## select * from log_login_loginhistory order by creation_date desc limit 50;
class LoginHistory(models.Model):
	"""
	Each sample needs a dataset 
	"""
	LOGIN_IN = 'in'
	LOGIN_OUT = 'out'
	
	owner = models.ForeignKey(User, related_name='login_history', blank=True, null=True, on_delete=models.CASCADE)
	
	ip = models.CharField(max_length=50, blank=True, null=True)
	operation = models.CharField(max_length=10)
	creation_date = models.DateTimeField(auto_now_add=True, verbose_name='Uploaded Date')
	description = models.TextField(default="")
	
	def __str__(self):
		return "{} {} {}".format(self.owner.username, self.ip, self.creation_date)
	class Meta:
		ordering = ['-creation_date', ]
	
