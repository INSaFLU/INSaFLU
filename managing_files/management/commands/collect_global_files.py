'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from managing_files.models import Project
from constants.constants import TypePath
from utils.utils import Utils
from utils.software import Software
from utils.collect_extra_data import CollectExtraData
import os

class Command(BaseCommand):
	'''
	classdocs
	'''
	## args = '<poll_id poll_id ...>'
	help = "Create global files by project sample."

	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	# A command must define handle()
	def handle(self, *args, **options):
		self.stdout.write("Stating")
		
		collect_extra_data = CollectExtraData()
		project_id = 1
		try:
			project = Project.objects.get(pk=project_id)
			collect_extra_data.collect_extra_data_for_project(project, project.owner, None)
			self.stdout.write("End")
		except Project.DoesNotExist as e:
			self.stdout.write("Project id '{}' does not exist".format(project_id))
		