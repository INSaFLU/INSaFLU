'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from managing_files.models import Reference
from constants.constants import TypePath
from utils.utils import Utils
from utils.software import Software
import os, logging

class Command(BaseCommand):
	'''
	classdocs
	'''
	## args = '<poll_id poll_id ...>'
	help = "Create a bed file for all genbank."

	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	# A command must define handle()
	def handle(self, *args, **options):
		self.stdout.write("Stating")
		
		utils = Utils()
		software = Software()
		count = 0
		for reference in Reference.objects.all():
			count += 1
			bed_path = reference.get_reference_bed(TypePath.MEDIA_ROOT)
			self.stdout.write("{} {}: ".format(count, reference.reference_genbank_name))
			self.logger_production.info("{} {}: ".format(count, reference.reference_genbank_name))
			self.logger_debug.info("{} {}: ".format(count, reference.reference_genbank_name))
			if os.path.exists(bed_path):
				self.stdout.write("already exist\n")
				continue
		
			### create bed and index for genbank
			utils.from_genbank_to_bed(reference.get_reference_gbk(TypePath.MEDIA_ROOT), reference.get_reference_bed(TypePath.MEDIA_ROOT))
			software.create_index_files_from_igv_tools(reference.get_reference_bed(TypePath.MEDIA_ROOT))
			self.stdout.write("created\n")
			
		self.stdout.write("End")