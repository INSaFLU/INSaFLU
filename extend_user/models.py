
from django.db import models
from django.contrib.auth.models import User
from django.db.models.signals import post_save
from django.dispatch import receiver

class Profile(models.Model):
	"""
	has the name of the institution of the account
	"""
	user = models.OneToOneField(User, on_delete=models.CASCADE)
	institution = models.TextField(max_length=100, blank=True)
	email_confirmed = models.BooleanField(default=False)
	
	## user only has a possibility to view a project 
	only_view_project = models.BooleanField(default=False)
	
	## some limits by user
	max_references = models.IntegerField(default=30)
	max_samples = models.IntegerField(default=500)
	max_file_size_fastq = models.IntegerField(default=50000000)
	max_length_reference_fasta = models.IntegerField(default=20000)
	max_sequence_reference = models.IntegerField(default=20)

@receiver(post_save, sender=User)
def create_user_profile(sender, instance, created, **kwargs):
	if created:
		Profile.objects.create(user=instance)

@receiver(post_save, sender=User)
def save_user_profile(sender, instance, **kwargs):
	if (instance != None): instance.profile.save()
	
# @receiver(post_save, sender=User)
# def update_user_profile(sender, instance, created, **kwargs):
# 	if created and Profile.objects.filter(user=instance).count() == 0:
# 		Profile.objects.create(user=instance)
# 	instance.profile.save()