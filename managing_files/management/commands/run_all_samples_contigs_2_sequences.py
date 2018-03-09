'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from managing_files.models import Sample
from utils.software import Software, Contigs2Sequences
from constants.constants import TypePath
from utils.utils import Utils
import os

class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Run all samples to collect contigs and copy to sequences."

	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	
	# A command must define handle()
	def handle(self, *args, **options):
		
		utils = Utils()
		software = Software()
		for sample in Sample.objects.all():
			if (sample.is_deleted): continue
			if (not sample.is_ready_for_projects): continue
			if (os.path.exists(sample.get_draft_contigs_output(TypePath.MEDIA_ROOT))): continue
			
			if (not os.path.exists(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True))): continue
			
			out_dir = utils.get_temp_dir()
			cmd = software.run_spades(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True),\
					sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False), out_dir)
			file_out = os.path.join(out_dir, "contigs.fasta")

			if (os.path.exists(file_out)):
				b_run_tests = False
				contigs_2_sequences = Contigs2Sequences(b_run_tests)
				(out_file_clean, clean_abricate_file) = contigs_2_sequences.identify_contigs(file_out)
				## copy the contigs from spades
				utils.copy_file(out_file_clean, sample.get_draft_contigs_output(TypePath.MEDIA_ROOT))
				utils.copy_file(clean_abricate_file, sample.get_draft_contigs_abricate_output(TypePath.MEDIA_ROOT))
				
				if (os.path.exists(out_file_clean)): os.unlink(out_file_clean)
				if (os.path.exists(clean_abricate_file)): os.unlink(clean_abricate_file)
			utils.remove_dir(out_dir)
		