'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from django_q.cluster import Cluster

class Command(BaseCommand):
	'''
	classdocs
	'''
	## args = '<poll_id poll_id ...>'
	help = "Stop django-q if its running."

	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	# A command must define handle()
	def handle(self, *args, **options):
		self.stdout.write("Hello world!")
