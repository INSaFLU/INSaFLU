'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from datasets.models import Dataset
from utils.collect_dataset_data import CollectExtraDatasetData
from django.contrib.auth.models import User
from django.conf import settings
import logging, datetime

class Command(BaseCommand):
	'''
	classdocs
	Ex: python3 manage.py collect_global_dataset_files --dataset_id 11 --user_id 1
	
	'''
	help = "Create global files by dataset sample."

	## logging
	if settings.DEBUG: logger = logging.getLogger("fluWebVirus.debug")
	else: logger = logging.getLogger("fluWebVirus.production")
	
	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	## https://docs.djangodataset.com/en/dev/howto/custom-management-commands/
	def add_arguments(self, parser):
		parser.add_argument('--dataset_id', type=int, help='Dataset id to process')
		parser.add_argument('--user_id', nargs='?', type=int, help='User id to join to this process')
	
	# A command must define handle()
	def handle(self, *args, **options):
		
		collect_extra_data = CollectExtraDatasetData()
		dataset_id = options['dataset_id']
		user_id = options['user_id']
		self.stdout.write("Starting for dataset_id: " + str(dataset_id))
		self.logger.info("Starting for dataset_id: " + str(dataset_id))
		try:
			dataset = Dataset.objects.get(pk=dataset_id)
			if (user_id == None): user = dataset.owner
			else: user = User.objects.get(pk=user_id)
			
			### change date
			dataset.last_change_date = datetime.datetime.now()
			dataset.save()
			
			### collect extra data for datasets
			collect_extra_data.collect_extra_data_for_dataset(dataset, user)
			self.stdout.write("End")
		except Dataset.DoesNotExist as e:
			self.stdout.write("Dataset id '{}' does not exist.".format(dataset_id))
		except User.DoesNotExist as e:
			self.stdout.write("Error: User id '{}' does not exist.".format(user_id))
