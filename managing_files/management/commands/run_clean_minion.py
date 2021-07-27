'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from managing_files.models import Sample
from utils.software_minion import SoftwareMinion
from django.contrib.auth.models import User
import logging

class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Run clean minion and identify species."
	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	## https://docs.djangoproject.com/en/dev/howto/custom-management-commands/
	def add_arguments(self, parser):
		parser.add_argument('--sample_id', type=int, help='Sample id to process')
		parser.add_argument('--user_id', nargs='?', type=int, help='User id to join to this process')
	
	# A command must define handle()
	def handle(self, *args, **options):
		
		software_minion = SoftwareMinion()
		sample_id = options['sample_id']
		user_id = options['user_id']
		self.stdout.write("Starting for sample_id: " + str(sample_id))
		self.logger_production.info("Starting for sample_id: " + str(sample_id))
		self.logger_debug.info("Starting for sample_id: " + str(sample_id))
		try:
			sample = Sample.objects.get(pk=sample_id)
			if (user_id == None): user = sample.owner
			else: user = User.objects.get(pk=user_id)
			## identify species
			b_make_identify_species = True
			b_return = software_minion.run_clean_minion(sample, user, b_make_identify_species)
			self.stdout.write("Resulting: " + str(b_return))
			self.stdout.write("End")
		except Sample.DoesNotExist as e:
			self.stdout.write("Error: Sample id '{}' does not exist.".format(sample_id))
		except User.DoesNotExist as e:
			self.stdout.write("Error: User id '{}' does not exist.".format(user_id))
		