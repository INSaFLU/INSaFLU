'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from django.db import transaction
from utils.lock_atomic_transaction import LockedAtomicTransaction
from managing_files.models import MetaKey

class Command(BaseCommand):
	'''
	classdocs
	'''
	help = "Test lock tables with LockedAtomicTransaction. nothing else"
	## logging
	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	def add_arguments(self, parser):
		parser.add_argument('--key', help='Sample id to process')
		
	def get_meta_key(self, meta_key_name_):
		"""
		
		"""
		with LockedAtomicTransaction(MetaKey):
			try:
				metaKey = MetaKey.objects.get(name=meta_key_name_)
			except MetaKey.DoesNotExist:
				metaKey = MetaKey()
				metaKey.name = meta_key_name_
				metaKey.save()
			return metaKey

	def get_meta_key_without(self, meta_key_name_):
		"""
		
		"""
		try:
			metaKey = MetaKey.objects.get(name=meta_key_name_)
		except MetaKey.DoesNotExist:
			metaKey = MetaKey()
			metaKey.name = meta_key_name_
			metaKey.save()
		return metaKey
	
	@transaction.atomic
	def transaction_atomic(self):
		
		meatkey = self.get_meta_key(self.meta_key_name)
#		meatkey = self.get_meta_key_without(self.meta_key_name)
		
		# Ask for `input` so execution will pause and wait for input.
		input('Execution is paused and you can now inspect the database.\n'
            'Press return/enter key to continue:')
		
		
	# A command must define handle()
	def handle(self, *args, **options):
		
		self.meta_key_name = options['key']
		
		#### set default fields
		self.stdout.write("Key to process: {}".format(self.meta_key_name))
		
		print("Last wait...")
		self.transaction_atomic()
		self.stdout.write("count: {}".format(MetaKey.objects.filter(name=self.meta_key_name).count()))
		
		self.stdout.write("End")
		