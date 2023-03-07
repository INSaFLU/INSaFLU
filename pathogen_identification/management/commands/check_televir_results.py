'''
Created on January 30, 2023

@author: daniel.sobral
'''
import os
import logging
import pandas as pd
from django.core.management import BaseCommand
from django.contrib.auth.models import User
from django.db import transaction
from managing_files.models import Sample
from pathogen_identification.models import Projects, RunMain, ParameterSet, FinalReport
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager
from settings.constants_settings import ConstantsSettings

class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Checks if a given televir project finished."

	# logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")

	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)

	def add_arguments(self, parser):
		parser.add_argument('--project_name', nargs='?', type=str,
							required=True, help='Project Name')				
		parser.add_argument('--user_login', nargs='?', type=str, required=True,
							help='User login of the project and sample owner')	
		#parser.add_argument('--report_file', nargs='?', type=str, required=True,
		#					help='Path for the report file to be saved (if project finished)')								

	# A command must define handle()
	def handle(self, *args, **options):

		project_name = options['project_name']
		account = options['user_login']
		#report_file = options['report_file']

		try:

			user = User.objects.get(username=account)

			project = Projects.objects.get(
                    name__iexact=project_name,
                    is_deleted=False,
                    owner__username=user.username,
                )
			
			runs = RunMain.objects.filter(project = project)
			if(len(runs)>0):
				finished = True
				for run in runs:
					if(run.parameter_set.status != ParameterSet.STATUS_FINISHED): 
						finished = False
				if(finished):

					# Experiments to save file etc...
					#runids = map(lambda x: x.id, runs)
					#all_reports = FinalReport.objects.filter(
            		#	run__project__pk=int(pk1), sample__pk=int(pk2)
        			#).order_by("-coverage")
					#all_reports = FinalReport.objects.filter(
            		#	run__project__pk__in=runids
        			#).order_by("-coverage")
					#report_pd = pd.DataFrame(all_reports.values())
					#report_pd.to_csv(report_file, sep="\t")
					
					self.stdout.write("Project finished.")
					
				else:
					self.stdout.write("Error: Project not finished yet.")

			else:
				self.stdout.write("Error: no runs found for this project.")

		except User.DoesNotExist as e:
			self.stdout.write("Error: User '{}' does not exist.".format(account))                
		except Exception as e:
			self.stdout.write("Error: {}.".format(e))
           

