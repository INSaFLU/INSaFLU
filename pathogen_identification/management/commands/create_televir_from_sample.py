'''
Created on January 30, 2023

@author: daniel.sobral
'''
import os
import logging
from django.core.management import BaseCommand
from django.contrib.auth.models import User
from django.db import transaction
from managing_files.models import Sample
from pathogen_identification.models import Projects, PIProject_Sample
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager
from settings.constants_settings import ConstantsSettings
from utils.process_SGE import ProcessSGE

class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Create a TELEVIR project and add a sample to it. Returns project id."

	# logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")

	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)

	def add_arguments(self, parser):
		parser.add_argument('--name', nargs='?', type=str,
							required=True, help='Sample Name')
		parser.add_argument('--user_login', nargs='?', type=str, required=True,
							help='User login of the sample owner')
		parser.add_argument('--technology', nargs='?', type=str, required=False,
							help='Technology: Illumina or ONT (ONT by default)')							

	# A command must define handle()
	def handle(self, *args, **options):

		sample_name = options['name']
		account = options['user_login']
		technology = ConstantsSettings.TECHNOLOGY_minion
		if(options['technology']):
			# Todo check if technology is one of the ones in ConstantsSettings
			technology = options['technology']

		try:
			
			user = User.objects.get(username=account)

			project_count = Projects.objects.filter(
				name__iexact=sample_name,
                is_deleted=False,
                owner__username=user.username,
            ).count()
			
			if(project_count > 0):
				self.stdout.write("Error: Project '{}' already exists.".format(sample_name))
			else:

				sample = Sample.objects.get(name=sample_name, owner=user)

				with transaction.atomic():
					project = Projects()
					project.name = sample_name
					project.owner = user
					project.owner_id = user.id
					project.technology = technology
					project.save()
					project_sample_input = sample.file_name_1
					if sample.is_valid_2:
						project_sample_input += ";" + sample.file_name_2
					project_sample = PIProject_Sample()
					project_sample.project = project
					project_sample.sample = sample
					project_sample.name = sample.name
					project_sample.input = project_sample_input
					project_sample.technology = sample.type_of_fastq
					project_sample.report = "report"
					project_sample.save()
					self.stdout.write("Project created with id {}.".format(project.id))   

				process_SGE = ProcessSGE()

				utils = Utils_Manager()
				runs_to_deploy = utils.check_runs_to_deploy(user, project)

				if runs_to_deploy:
					taskID = process_SGE.set_submit_televir_job(
                    	user=user,
                    	project_pk=project.pk,
                	)
					self.stdout.write("Project submitted as task {}.".format(project.id, taskID))   

		except User.DoesNotExist as e:
			self.stdout.write("Error: User '{}' does not exist.".format(account))                
		except Exception as e:
			self.stdout.write("Error: {}.".format(e))
           

