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
	help = "Confirm a account by user name. If you want activate last account pass the argument 'last_account_created'"

	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	def add_arguments(self, parser):
		parser.add_argument('--username', type=str, help='User name to activate email...')

	# A command must define handle()
	def handle(self, *args, **options):
		
		username = options['username']
		
		### activate last account created
		if username == 'last_account_created':
			user = User.objects.all().order_by('-id')[0]
			username = user.username
		else:
			try:
				user = User.objects.get(username__iexact=username.strip())
			except User.DoesNotExist as e:
				self.stdout.write("Error: User login '{}' does not exist.".format(username))
				return
		
		is_email_active = user.profile.email_confirmed
		if (is_email_active): self.stdout.write("User name email '{}' is already confirmed.".format(username))
		else:
			profile = user.profile
			profile.email_confirmed = True
			profile.save()
			self.stdout.write("User name email '{}' is confirmed.".format(username))

