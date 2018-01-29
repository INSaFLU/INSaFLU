'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from managing_files.models import Sample
from utils.software import Software
from django.contrib.auth.models import User

class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Run trimmomatic and identify species."

	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	## https://docs.djangoproject.com/en/dev/howto/custom-management-commands/
	def add_arguments(self, parser):
		parser.add_argument('--sample_id', type=int, help='Sample id to process')
		parser.add_argument('--user_id', nargs='?', type=int, help='User id to join to this process')
	
	# A command must define handle()
	def handle(self, *args, **options):
		
		software = Software()
		sample_id = options['sample_id']
		user_id = options['user_id']
		self.stdout.write("Starting for sample_id: " + str(sample_id))
		try:
			sample = Sample.objects.get(pk=sample_id)
			if (user_id == None): user = sample.owner
			else: user = User.objects.get(pk=user_id)
			software.run_fastq_and_trimmomatic_and_identify_species(sample, user)
			self.stdout.write("End")
		except Sample.DoesNotExist as e:
			self.stdout.write("Error: Sample id '{}' does not exist.".format(sample_id))
		except User.DoesNotExist as e:
			self.stdout.write("Error: User id '{}' does not exist.".format(user_id))
		