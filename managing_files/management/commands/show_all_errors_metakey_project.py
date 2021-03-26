'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from managing_files.models import MetaKeyProject

class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Show all users."

	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	def add_arguments(self, parser):
		parser.add_argument('--limit', nargs='?', type=int, help='Number of errors to print')
	
	# A command must define handle()
	def handle(self, *args, **options):
		limit = options['limit']
		if (limit is None): limit = 10
		
		metakey_project = MetaKeyProject.objects.filter(value="Error").order_by('id').reverse();
		self.stdout.write("Pk".ljust(10) + "Value".ljust(20) + "ProjectId".ljust(15) + \
						"Date".ljust(20) + "MetaTag".ljust(20) + "Description".ljust(20))
		for _, meta_key in enumerate(metakey_project):
			self.stdout.write(str(meta_key.pk).ljust(10) + str(meta_key.value).ljust(20) + str(meta_key.project_id).ljust(15) +\
						self.format_date(meta_key.creation_date).ljust(20) + meta_key.meta_tag.name.ljust(20) +\
						str(meta_key.description).ljust(20))
			self.stdout.write("# {} ".format(_) + "#" * 40)
			if (_ > limit): break

	def format_date(self, date_and_time):
		return date_and_time.strftime("%d/%m/%Y %H:%M")