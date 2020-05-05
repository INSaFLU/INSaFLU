'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from utils.parse_in_files import ParseInFiles
from django.contrib.auth.models import User
import logging

class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Run link files from fastq to samples."
	
	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	## https://docs.djangoproject.com/en/dev/howto/custom-management-commands/
	def add_arguments(self, parser):
		parser.add_argument('--user_id', type=int, help='User id to attached fastq to samples')
	
	# A command must define handle()
	def handle(self, *args, **options):
		
		parse_in_files = ParseInFiles()
		user_id = options['user_id']
		self.stdout.write("Starting for user_id: " + str(user_id))
		self.logger_production.info("Starting for user_id: " + str(user_id))
		self.logger_debug.info("Starting for user_id: " + str(user_id))
		try:
			user = User.objects.get(pk=user_id)
			parse_in_files.link_files(user)
			self.stdout.write("End")
		except User.DoesNotExist as e:
			self.stdout.write("Error: User id '{}' does not exist.".format(user_id))
		