'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from django.contrib.auth.models import User
from utils.collect_extra_data import CollectExtraData
import logging

class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Create a list of samples by user."
	
	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")

	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	## https://docs.djangoproject.com/en/dev/howto/custom-management-commands/
	def add_arguments(self, parser):
		parser.add_argument('--user_id', type=int, help='User id to create a sample list')
	
	# A command must define handle()
	def handle(self, *args, **options):
		
		user_id = options['user_id']
		self.stdout.write("Start creating sample list for user_id: " + str(user_id))
		self.logger_production.info("Start creating sample list for user_id: " + str(user_id))
		self.logger_debug.info("Start creating sample list for user_id: " + str(user_id))
		try:
			user = User.objects.get(pk=user_id)
			self.stdout.write("User name: " + user.username)
			collect_extra_data = CollectExtraData()
			collect_extra_data.collect_sample_list(user)
			self.stdout.write("End")
		except User.DoesNotExist as e:
			self.stdout.write("Error: User id '{}' does not exist.".format(user_id))
			self.logger_production.info("Error: User id '{}' does not exist.".format(user_id))
			self.logger_debug.info("Error: User id '{}' does not exist.".format(user_id))
		