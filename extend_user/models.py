
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