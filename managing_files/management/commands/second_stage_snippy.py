'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from managing_files.models import ProjectSample
from utils.software import Software
from django.contrib.auth.models import User
import logging

class Command(BaseCommand):
	'''
	classdocs
	Ex: python3 manage.py second_stage_snippy --project_sample_id 19 --user_id 1
	
	'''
	help = "Run second stage snippy and freebayes."

	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	## https://docs.djangoproject.com/en/dev/howto/custom-management-commands/
	def add_arguments(self, parser):
		parser.add_argument('--project_sample_id', type=int, help='Project sample id to process')
		parser.add_argument('--user_id', nargs='?', type=int, help='User id to join to this process')
	
	# A command must define handle()
	def handle(self, *args, **options):
		
		software = Software()
		project_sample_id = options['project_sample_id']
		user_id = options['user_id']
		self.stdout.write("Starting for project_sample_id: " + str(project_sample_id))
		self.logger_production.info("Starting for project_sample_id: " + str(project_sample_id))
		self.logger_debug.info("Starting for project_sample_id: " + str(project_sample_id))
		try:
			project_sample = ProjectSample.objects.get(pk=project_sample_id)
			if (user_id == None): user = project_sample.project.owner
			else: user = User.objects.get(pk=user_id)
			software.process_second_stage_snippy_coverage_freebayes(project_sample, user)
			self.stdout.write("End")
		except ProjectSample.DoesNotExist as e:
			self.stdout.write("Error: ProjectSample id '{}' does not exist.".format(project_sample_id))
		except User.DoesNotExist as e:
			self.stdout.write("Error: User id '{}' does not exist.".format(user_id))
		