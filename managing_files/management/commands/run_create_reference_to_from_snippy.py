'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from managing_files.models import Project
from utils.software import Software
import logging

class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Create new reference for snippy."

	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	
	# A command must define handle()
	def handle(self, *args, **options):
		
		software = Software()
		for project in Project.objects.all():
			if (project.is_deleted): continue
			self.stdout.write("Processing project: {}".format(project.name))
			for project_sample in project.project_samples.all():
				if (not project_sample.get_is_ready_to_proccess()): continue
				self.stdout.write("Processing sample: {}".format(project_sample.sample.name))
				software.creat_new_reference_to_snippy(project_sample)
		self.stdout.write("End")
		self.stdout.write("Finished")