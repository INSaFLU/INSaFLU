'''
Created on Jan 5, 2018

@author: mmp
'''
import logging

from django.contrib.auth.models import User
from django.core.management import BaseCommand

from managing_files.models import Project
from settings.default_software_project_sample import DefaultProjectSoftware
from utils.coverage import DrawAllCoverage


class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Create coverage images to all samples in a specific project."

	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	## https://docs.djangoproject.com/en/dev/howto/custom-management-commands/
	def add_arguments(self, parser):
		parser.add_argument('--project_id', type=int, help='Project id to process')
	
	# A command must define handle()
	def handle(self, *args, **options):
		
		draw_all_coverage = DrawAllCoverage()
		project_id = options['project_id']
		self.stdout.write("Starting for project_id: " + str(project_id))
		self.logger_production.info("Starting for project_id: " + str(project_id))
		self.logger_debug.info("Starting for project_id: " + str(project_id))
		try:
			project = Project.objects.get(pk=project_id)
			default_project_software = DefaultProjectSoftware()
			for project_sample in project.project_samples.all():
				software_mdcg_name= default_project_software.get_software_project_sample_mdcg_illumina(
					project_sample,
					)
				if (not project_sample.get_is_ready_to_proccess()): continue
				self.stdout.write("Processing sample: {}".format(project_sample.sample.name))
				draw_all_coverage.draw_all_coverages(project_sample, software_name = software_mdcg_name)
			self.stdout.write("End")
		except Project.DoesNotExist as e:
			self.stdout.write("Project id '{}' does not exist.".format(project_id))
		self.stdout.write("Finished")			self.stdout.write("End")
		except Project.DoesNotExist as e:
			self.stdout.write("Project id '{}' does not exist.".format(project_id))
		self.stdout.write("Finished")