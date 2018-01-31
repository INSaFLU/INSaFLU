'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from django.contrib.auth.models import User

class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Show all users."

	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	
	# A command must define handle()
	def handle(self, *args, **options):
		
		list_users = User.objects.all().order_by('date_joined');
		self.stdout.write("Pk".ljust(10) + "Name".ljust(40) + "#projects".ljust(13) + "#samples".ljust(11) +  "Last access".ljust(20) + "Last IP".ljust(18) + "SuperUser".ljust(15) + "Active".ljust(10) + "Email Active".ljust(20))
		self.stdout.write("-" * (10 + 40 + 13 + 11 + 20 + 18 + 15 + 10 + 20))
		for user in list_users:
			
			username = user.username
			last_login = "" if user.last_login == None else self.format_date(user.last_login)
			pk = str(user.pk)
			number_projects = str(user.project.all().count())
			number_samples = str(user.sample.all().count())
			is_active = str(user.is_active)
			is_superuser = str(user.is_superuser)
			is_email_active = str(user.profile.email_confirmed)
			lst_login_set = user.login_history.all().order_by('-creation_date').first()
			lst_login_ip = lst_login_set.ip if lst_login_set != None else ""
			self.stdout.write(pk.ljust(10) + username.ljust(40) + number_projects.ljust(13) + number_samples.ljust(11) +\
						last_login.ljust(20) + lst_login_ip.ljust(18) + is_superuser.ljust(15) + is_active.ljust(10) + is_email_active.ljust(20))
		
		self.stdout.write("\n" + "=" * (10 + 40 + 40 + 40) + "\n\n")
		self.stdout.write("Pk".ljust(10) + "Name".ljust(40) + "Email".ljust(40) + "Institution".ljust(40))
		self.stdout.write("-" * (10 + 40 + 40 + 40))
		for user in list_users:
			username = user.username
			email = user.email
			if (email == None): email = ""
			pk = str(user.pk)
			self.stdout.write(pk.ljust(10) + username.ljust(40) + email.ljust(40) + user.profile.institution.ljust(40))
			
	def format_date(self, date_and_time):
		return date_and_time.strftime("%d/%m/%Y %H:%M")