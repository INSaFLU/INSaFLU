'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from utils.software_pangolin import SoftwarePangolin

class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Check software versions."

	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	
	# A command must define handle()
	def handle(self, *args, **options):

		self.stdout.write("Check pangolin...")
		software_pangolin = SoftwarePangolin()
		software_pangolin.run_pangolin_update()
		
		self.stdout.write("Done")
