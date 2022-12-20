'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from utils.software_pangolin import SoftwarePangolin
import logging

class Command(BaseCommand):
	'''
	classdocs
	Ex: python3 manage.py update_pangolin
	
	'''
	help = "Update pangolin learn."

	## logging
	logger = logging.getLogger("fluWebVirus.update_pangolin")
	
	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	# A command must define handle()
	def handle(self, *args, **options):
		
		software_pangolin = SoftwarePangolin()
		self.stdout.write("Starting update pangolin")
		self.logger.info("Starting update pangolin")

		try:
			software_pangolin.run_pangolin_update()
			self.logger.info("End update pangolin")
			self.stdout.write("Success update pangolin")
		except Exception as e:
			message = "Fail to update pangolin: " + str(e)
			self.logger.info(message)
			self.stdout.write(message)
