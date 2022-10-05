'''
Created on Jan 5, 2018

@author: mmp
'''
from django.core.management import BaseCommand
from managing_files.models import Reference, Sample, Project, ProjectSample, UploadFiles
from constants.constants import TypePath, TypeFile
from constants.software_names import SoftwareNames
from django.contrib.auth.models import User
from utils.utils import Utils
from utils.software import Software
import datetime, dateutil.relativedelta
import argparse
import logging, sys

### Define date for NO-DATE parameter 
NO_DATE = datetime.datetime(1980, 1, 1)

def str2bool(v):
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')

def str2date(v):
	if len(v) == 0: return NO_DATE
	try:
		return datetime.datetime.strptime(v, '%Y-%m-%d')
	except ValueError:
		pass
	
	try:
		return datetime.datetime.strptime(v, '%Y/%m/%d')
	except ValueError:
		pass
	
	try:
		return datetime.datetime.strptime(v, '%y-%m-%d')
	except ValueError:
		pass
	
	try:
		return datetime.datetime.strptime(v, '%y/%m/%d')
	except ValueError:
		raise argparse.ArgumentTypeError("\nNot a valid date '{}'. Must be the format 'YYYY-MM-DD'".format(v))

class Command(BaseCommand):
	'''
	Remove files that are marked to be removed.
	
	This command can run at least once a week, trigger by the crontab
	This is to prevent deleting files in runtime increasing the risk of attack
	'''
	help = "Run remove files when they are marked to be removed."
	
	## logging
	logger_remove_files = logging.getLogger("fluWebVirus.remove_files")
	utils = Utils()
	
	def __init__(self, *args, **kwargs):
		super(Command, self).__init__(*args, **kwargs)
	
	def add_arguments(self, parser):
		parser.add_argument('--only_identify_files', type=str2bool, nargs='?',
			const=True, default=True, 
			help='Only to identify files to remove. Does not remove them. The output it is in log file.')
		parser.add_argument('--remove_files_after_days', type=int, nargs='?',
			const=10, default=10, 
			help='Remove files after X days passed from been removed in web site by the users.\n' +
			"If the user set the parameter 'end_date' this parameter will be deactivated.")
		parser.add_argument('--login_to_pass', type=str, nargs='?',
			const="", default="",
			help='Login accounts to pass the remove files. Can be more than one separated by comma. Ex: mmp,demo')
		parser.add_argument('--login_to_process', type=str, nargs='?',
			const="", default="",
			help='Login accounts to remove files. Can be more than one separated by comma. Ex: mmp,demo')
		parser.add_argument('--delete_original_fastq', type=str2bool, nargs='?',
			const=False, default=False,
			help='Delete fastq original files if True. It is possible to remove these files because the system only works Trimmomatic output files.')
		parser.add_argument('--delete_trimmomatic_nanofilt_fastq', type=str2bool, nargs='?',
			const=False, default=False,
			help='Delete Trimmomatic/Nanofilt result files if True. If they remove it is not possible to not include the sample in a project anymore.')
		parser.add_argument('--from_date', type=str2date, nargs='?',
			const="", default="",
			help='From this date till more recent. It is greater and equal. Format (YYYY-MM-DD)')
		parser.add_argument('--end_date', type=str2date, nargs='?',
			const="", default="",
			help='From this date till older. It is lower and equal. Format (YYYY-MM-DD).')


	# A command must define handle()
	def handle(self, *args, **options):
		
		### remove files from the disk if were removed more than X days
		self.REMOVE_FILES_AFTER_DAYS = options['remove_files_after_days']
		
		### if true does not remove files, only save it in the log file
		only_identify_files = options['only_identify_files']
		delete_original_fastq = options['delete_original_fastq']
		delete_trimmomatic_nanofilt_fastq = options['delete_trimmomatic_nanofilt_fastq']
		self.out_message("\nMode - {}\nDelete Original fastq files: {}\n".format(
			"Only identify files" if only_identify_files else "Identify files and remove them.",
			str(delete_original_fastq)) +
			"Delete Trimmomatic/Nanofilt fastq files: {}".format(
			str(delete_trimmomatic_nanofilt_fastq)), False)
		### dates
		if options['from_date'] != NO_DATE and options['end_date'] != NO_DATE and options['end_date'] <= options['from_date']:
			self.out_message("\n'From Date' must be bigger to 'End Date'", True)
			sys.exit(1)

		if (options['from_date'] != NO_DATE): self.out_message("From date: {}".format(options['from_date']), False)
		if (options['end_date'] != NO_DATE): self.out_message("End date: {}".format(options['end_date']), False)
		if (options['from_date'] == NO_DATE and options['end_date'] == NO_DATE):
			self.out_message("Remove files after {} days passed from been removed in web site".format(self.REMOVE_FILES_AFTER_DAYS) + \
					" by the users.", False)
		
		## logins to pass
		logins_to_pass = options['login_to_pass']
		lst_accounts_to_pass = []
		if (len(logins_to_pass) > 0): lst_accounts_to_pass = logins_to_pass.strip().split(',')
		
		## logins to process
		logins_to_process = options['login_to_process']
		lst_accounts_to_process = []
		if (len(logins_to_process) > 0): lst_accounts_to_process = logins_to_process.strip().split(',')
		
		if (len(lst_accounts_to_process) > 0 and len(lst_accounts_to_pass) > 0):
			self.out_message("\nYou can not use '--login_to_process' and '--login_to_pass' simultaneously.", True)
			sys.exit(1)
			
		### info about the logins
		if len(lst_accounts_to_pass) > 0: self.out_message("Login to pass: {}".format(",".join(lst_accounts_to_pass)), False)
		if len(logins_to_process) > 0: self.out_message("Login to process: {}".format(",".join(lst_accounts_to_process)), False)
		if options['end_date'] != NO_DATE and options['from_date'] != NO_DATE:
			self.out_message("Warning: the parameter 'remove_files_after_days = {}' will be deactivated.".format(self.REMOVE_FILES_AFTER_DAYS), False)
					
		### message to user, to confirm all the information
		read_key = input("\nDo you want to proceed [Yy/Nn]? ")
		if (not read_key in ['Y', 'y']): sys.exit(0)
		
		### set star and end dates
		self.has_end_date = False 
		if (options['end_date'] != NO_DATE):
			self.has_end_date = True		###if true the parameter 'remove_files_after_days' will be disable
			self.end_date = options['end_date']
		self.has_from_date = False 
		if (options['from_date'] != NO_DATE):
			self.has_from_date = True
			self.from_date = options['from_date']
		
		## test the existence of accounts
		self.test_existence_accounts(lst_accounts_to_process, "to process")
		
		## test the existence of accounts
		self.test_existence_accounts(lst_accounts_to_pass, "to pass")
				
		### delete REFERENCES
		total_files_downloaded = self.delete_references(only_identify_files, lst_accounts_to_pass,
												lst_accounts_to_process)

		### delete SAMPLES
		total_files_downloaded += self.delete_all_samples(only_identify_files, lst_accounts_to_pass,
												lst_accounts_to_process)
		
		### delete original fastq SAMPLES
		if (delete_original_fastq):
			total_files_downloaded += self.delete_original_fastq_samples(only_identify_files, lst_accounts_to_pass,
												lst_accounts_to_process)

		if (delete_trimmomatic_nanofilt_fastq):
			total_files_downloaded += self.delete_trimmomatic_nanofilt_fastq_samples(only_identify_files,
						lst_accounts_to_pass, lst_accounts_to_process)
			
		### remove fastq.gz files not attached and already deleted
		total_files_downloaded += self.delete_upload_deleted_files_by_batch_and_not_attached(only_identify_files,
												lst_accounts_to_pass, lst_accounts_to_process)
		
		### remove project_samples "sample" is the sample and project sample was deleted
		total_files_downloaded += self.delete_project_samples(only_identify_files, lst_accounts_to_pass,
												lst_accounts_to_process)

		if (only_identify_files):
			if (total_files_downloaded == 0): self.out_message("No file were identified", False)
			else: self.out_message("Files identified: {}".format(total_files_downloaded), False)
		else:
			if (total_files_downloaded == 0): self.out_message("No files were removed", False)
			else: self.out_message("Files removed: {}".format(total_files_downloaded), False)
			
	def out_message(self, message, b_error):
		if (b_error): self.logger_remove_files.error(message.strip())
		else: self.logger_remove_files.info(message.strip())
		self.stdout.write(message)

	def test_existence_accounts(self, lst_accounts_to_pass, sz_message = ""):
		if (len(lst_accounts_to_pass) > 0):
			self.out_message("\nTesting existence of the accounts {}- {}".format(sz_message,
				",".join(lst_accounts_to_pass)), False)
			lst_not_exist = []
			for account in lst_accounts_to_pass:
				try:
					User.objects.all().get(username=account)
				except User.DoesNotExist:
					lst_not_exist.append(account)
			
			### testing accounts not exist
			if (len(lst_not_exist) > 0):
				for account in lst_not_exist:
					self.out_message("\nAccount does not exist: {}".format(account), True)
				sys.exit(1)
			self.out_message("Accounts OK", False)
				
	def is_out_date_to_delete(self, date_to_evaluate):
		## both dates are defined
		if (not date_to_evaluate is None and self.has_end_date and self.has_from_date):
			if (self.from_date > date_to_evaluate or self.end_date < date_to_evaluate): return True
			return False
		## only from_date not null
		if (not date_to_evaluate is None and self.has_from_date and \
			self.from_date > date_to_evaluate): return True
		## only end_date not null
		if (not date_to_evaluate is None and self.has_end_date and \
			self.end_date < date_to_evaluate): return True
		
		if (date_to_evaluate != None): removed_days = int(divmod((datetime.datetime.now() - date_to_evaluate).total_seconds(), 86400)[0])
		else: removed_days = 100000	## big number, older versions doesn't have this table field
		if (removed_days < self.REMOVE_FILES_AFTER_DAYS): return True
		return False

	def delete_references(self, only_identify_files, lst_accounts_to_pass, lst_accounts_to_process):
		"""
		param: only_identify_files
		param: lst_accounts_to_pass
		param: lst_accounts_to_process
		delete file references deleted by the users 
		"""
		### REFERENCES - removing references files
		self.out_message("\n### REFERENCES\n", False)
		references = Reference.objects.all().filter(is_deleted=True, is_deleted_in_file_system=False)
		count = 0
		for reference in references:
			### test the owners
			if (reference.owner.username in lst_accounts_to_pass): continue
			if (len(lst_accounts_to_process) > 0 and \
				not reference.owner.username in lst_accounts_to_process): continue
				
			### start deleting
			if (Project.objects.all().filter(reference=reference, is_deleted=False).count() == 0):
				files_removed = []
				files_to_remove = []
				files_to_remove.append(reference.get_reference_fasta(TypePath.MEDIA_ROOT))
				
				## test the days/date removed
				if (self.is_out_date_to_delete(reference.date_deleted)):
					#self.out_message("Not remove physically: {}; fail the date/remove days test {}.".format(files_to_remove[0], str(reference.date_deleted)), False)
					continue
				
				try:
					### try to remove extra files like gbk and index
					files_to_remove.append(reference.get_reference_gbk(TypePath.MEDIA_ROOT))
					files_to_remove.append(reference.get_reference_fasta_index(TypePath.MEDIA_ROOT))
					files_to_remove.append(reference.get_reference_bed(TypePath.MEDIA_ROOT))
					files_to_remove.append(reference.get_reference_bed_index(TypePath.MEDIA_ROOT))
						
					if (only_identify_files):
						files_removed = files_to_remove.copy()
					else:
						for path_to_remove in files_to_remove:
							if (self.utils.remove_file(path_to_remove)):
								files_removed.append(path_to_remove)
						
						if (len(files_removed) > 0):
							### save the flag in database
							reference.is_deleted_in_file_system = True
							reference.save()
				except Exception as e:
					self.out_message("Fail to remove: {}".format(str(e)), False)
					continue
				for file_path in files_removed:
					if (only_identify_files): self.out_message("Identified file: " + file_path, False)
					else: self.out_message("Remove file: " + file_path, False)
					count += 1
		
		if (only_identify_files):
			if (count == 0): self.out_message("None references were identified", False)
			else: self.out_message("Files identified: {}".format(count), False)
		else:
			if (count == 0): self.out_message("None references were removed", False)
			else: self.out_message("Files removed: {}".format(count), False)
		self.out_message("### END  REFERENCES\n", False)
		return count

	def delete_all_samples(self, only_identify_files, lst_accounts_to_pass, lst_accounts_to_process):
		"""
		delete all samples deleted by the users in the website
		param: only_identify_files
		param: lst_accounts_to_pass
		param: lst_accounts_to_process
		"""
		self.out_message("\n### Samples Deleted\n", False)
		samples_deleted = 0
		try:
			samples = Sample.objects.all().filter(is_deleted=True, is_deleted_in_file_system=False, is_valid_1=True, is_ready_for_projects=True)
			samples_deleted = self._delete_all_samples(samples, only_identify_files, lst_accounts_to_pass,
											lst_accounts_to_process)
		except Sample.DoesNotExist:
			pass
		
		### can not be ready for a project but it is in pipeline to run. It's better to wait one month 
		date_remove = datetime.datetime.now() - dateutil.relativedelta.relativedelta(months=1)
		samples = Sample.objects.all().filter(is_deleted=True, is_deleted_in_file_system=False, is_valid_1=True, 
								is_ready_for_projects=False, creation_date__lte=date_remove)
		samples_deleted += self._delete_all_samples(samples, only_identify_files, lst_accounts_to_pass,
											lst_accounts_to_process)
		
		if (samples_deleted == 0):
			if (only_identify_files): self.out_message("None samples were removed.", False)
			else: self.out_message("None samples were removed", False)
		self.out_message("### END  Samples Deleted\n", False)
		return samples_deleted

	def _delete_all_samples(self, samples, only_identify_files, lst_accounts_to_pass, lst_accounts_to_process):
		"""
		Delete all samples by array
		param: only_identify_files
		param: lst_accounts_to_pass
		param: lst_accounts_to_process
		"""
		count = 0
		for sample in samples:
			### test users
			if (sample.owner.username in lst_accounts_to_pass): continue
			if (len(lst_accounts_to_process) > 0 and \
				not sample.owner.username in lst_accounts_to_process): continue
				
			### start identify
			### don't remove samples that are in active projects
			if (ProjectSample.objects.all().filter(sample=sample, is_deleted=False).count() == 0):
				files_removed = []
				files_to_remove = []
				
				## test the days/date removed
				if (self.is_out_date_to_delete(sample.date_deleted)):
					#self.out_message("Not remove physically: {}; fail the date/remove days test {}.".format(sample.name, str(sample.date_deleted)), False)
					continue
				
				try:
					original_file_not_removed = False 
					if (sample.is_type_fastq_gz_sequencing()):	### illumina
						files_to_remove.append(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True))
						files_to_remove.append(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False))
						if (not sample.is_original_fastq_removed()):
							files_to_remove.append(sample.get_fastq(TypePath.MEDIA_ROOT, True))
							files_to_remove.append(sample.get_fastq(TypePath.MEDIA_ROOT, False))
							files_to_remove.append(sample.get_fastqc_output(TypePath.MEDIA_ROOT, True))
							files_to_remove.append(sample.get_fastqc_output(TypePath.MEDIA_ROOT, False))
							original_file_not_removed = True
						files_to_remove.append(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True))
						files_to_remove.append(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False))
						files_to_remove.append(sample.get_abricate_output(TypePath.MEDIA_ROOT))
					else:	### minion
						if (not sample.is_original_fastq_removed()):
							files_to_remove.append(sample.get_fastq(TypePath.MEDIA_ROOT, True))
							files_to_remove.append(sample.get_rabbitQC_output(TypePath.MEDIA_ROOT))
							original_file_not_removed = True
						files_to_remove.append(sample.get_nanofilt_file(TypePath.MEDIA_ROOT))
						files_to_remove.append(sample.get_rabbitQC_nanofilt(TypePath.MEDIA_ROOT))
						
					if (only_identify_files):
						files_removed = files_to_remove.copy()
					else:
						for path_to_remove in files_to_remove:
							if (self.utils.remove_file(path_to_remove)):
								files_removed.append(path_to_remove)
						
						if (len(files_removed) > 0):
							### save the flag in database
							sample.is_deleted_in_file_system = True
							sample.is_deleted_processed_fastq = True
							sample.date_deleted_processed_fastq = datetime.datetime.now()
							sample.save()
							
							### need to set the flag deleted in upload files by batch if the file was uploaded by batch
							if (original_file_not_removed):
								upload_files = UploadFiles.objects.all().filter(samples__id=sample.pk, 
											type_file__name=TypeFile.TYPE_FILE_fastq_gz, is_deleted_in_file_system=False)
								for upload_file in upload_files:
									upload_file.is_deleted = True
									if (upload_file.date_deleted == None):
										upload_file.date_deleted = datetime.datetime.now()
									## don't set this flag were, the next method will remove the files in "upload_file"
									#upload_file.is_deleted_in_file_system = True
									upload_file.save()
							
				except Exception as e:
					self.out_message("Fail to remove: {}".format(str(e)), False)
					continue
				for file_path in files_removed:
					if (file_path == None): continue
					if (only_identify_files): self.out_message("Identified file: " + file_path, False)
					else: self.out_message("Remove file: " + file_path, False)
					count += 1
		
		if (only_identify_files):
			if (count == 0): self.out_message("None samples were identified", False)
			else: self.out_message("Files identified: {}".format(count), False)
		else:
			if (count == 0): self.out_message("None samples were removed", False)
			else: self.out_message("Files removed: {}".format(count), False)
		return count


	def delete_original_fastq_samples(self, only_identify_files, lst_accounts_to_pass,
									lst_accounts_to_process):
		"""
		delete all samples deleted by the users 
		"""
		self.out_message("\n### Original fastq files\n", False)
		samples = Sample.objects.all().filter(is_deleted=False, is_deleted_in_file_system=False, is_valid_1=True, is_ready_for_projects=True)
		samples_deleted = self._delete_original_fastq(samples, only_identify_files, lst_accounts_to_pass,
									lst_accounts_to_process)
		
		if (samples_deleted == 0):
			if (only_identify_files): self.out_message("None original fastq.gz were identified", False)
			else: self.out_message("None original fastq.gz were removed", False)
		self.out_message("### END  Original fastq files\n", False)
		return samples_deleted

	def _delete_original_fastq(self, samples, only_identify_files, lst_accounts_to_pass,
							lst_accounts_to_process):
		"""
		Delete all samples by array
		"""
		count = 0
		for sample in samples:
			files_removed = []
			files_to_remove = []
			
			if (sample.is_original_fastq_removed()): continue
			if (sample.owner.username in lst_accounts_to_pass): continue
			if (len(lst_accounts_to_process) > 0 and \
				not sample.owner.username in lst_accounts_to_process): continue
				
			## can be removed already
			try:
				## Do not remove Quality Control outputs
				files_to_remove.append(sample.get_fastq(TypePath.MEDIA_ROOT, True))
				#if (sample.is_type_fastq_gz_sequencing()): files_to_remove.append(sample.get_fastqc_output(TypePath.MEDIA_ROOT, True))		## illumina default reads
				#else: files_to_remove.append(sample.get_rabbitQC_output(TypePath.MEDIA_ROOT))

				if (sample.is_valid_2): 
					files_to_remove.append(sample.get_fastq(TypePath.MEDIA_ROOT, False))
					#files_to_remove.append(sample.get_fastqc_output(TypePath.MEDIA_ROOT, False))

				## test the days/date removed
				if (self.is_out_date_to_delete(sample.creation_date)):
					#self.out_message("Not remove physically: {}; fail the date/remove days test {}.".format(files_to_remove[0], str()), False)
					continue
				
				if (only_identify_files):
					files_removed = files_to_remove.copy()
				else:
					for path_to_remove in files_to_remove:
						if (self.utils.remove_file(path_to_remove)):
							files_removed.append(path_to_remove)
					
					if (len(files_removed) > 0):
						### save the flag in database
						sample.is_deleted_in_file_system = True
						sample.save()
						
						### need to set the flag deleted in upload files by batch if the file was uploaded by batch
						upload_files = UploadFiles.objects.all().filter(samples__id=sample.pk, 
									type_file__name=TypeFile.TYPE_FILE_fastq_gz, is_deleted_in_file_system=False)
						for upload_file in upload_files:
							upload_file.is_deleted = True
							if (upload_file.date_deleted == None):
								upload_file.date_deleted = datetime.datetime.now()
							## don't set this flag were, the next method will remove the files in "upload_file"
							#upload_file.is_deleted_in_file_system = True
							upload_file.save()
							
			except Exception as e:
				self.out_message("Fail to remove: {}".format(str(e)), False)
				continue
			for file_path in files_removed:
				if (file_path == None): continue
				if (only_identify_files): self.out_message("Identified file: " + file_path, False)
				else: self.out_message("Remove file: " + file_path, False)
				count += 1
		
		if (only_identify_files):
			if (count == 0): self.out_message("None samples were identified", False)
			else: self.out_message("Files identified: {}".format(count), False)
		else:
			if (count == 0): self.out_message("None samples were removed", False)
			else: self.out_message("Files removed: {}".format(count), False)
		return count

	def delete_trimmomatic_nanofilt_fastq_samples(self, only_identify_files, lst_accounts_to_pass,
									lst_accounts_to_process):
		"""
		Delete all files processed from fastq Trimmomatic/Rabbit
		"""
		self.out_message("\n### trimmomatic-nanofilt-fastq files\n", False)
		samples = Sample.objects.all().filter(is_deleted=False, is_deleted_processed_fastq=False, is_valid_1=True, is_ready_for_projects=True)
		samples_deleted = self._delete_processed_fastq(samples, only_identify_files, lst_accounts_to_pass,
									lst_accounts_to_process)
		
		if (samples_deleted == 0):
			if (only_identify_files): self.out_message("None trimmomatic-nanofilt-fastq were identified", False)
			else: self.out_message("None trimmomatic-nanofilt-fastq were removed", False)
		self.out_message("### END trimmomatic-nanofilt-fastq files\n", False)
		return samples_deleted
		
	def _delete_processed_fastq(self, samples, only_identify_files, lst_accounts_to_pass,
							lst_accounts_to_process):
		"""
		Delete all files processed from fastq Trimmomatic/Rabbit
		"""
		count = 0
		for sample in samples:
			files_removed = []
			files_to_remove = []
			
			if (sample.is_processed_fastq_deleted()): continue
			if (sample.owner.username in lst_accounts_to_pass): continue
			if (len(lst_accounts_to_process) > 0 and \
				not sample.owner.username in lst_accounts_to_process): continue
				
			## test the days/date removed
			if (self.is_out_date_to_delete(sample.creation_date)):
				#self.out_message("Not remove physically: {}; fail the date/remove days test {}.".format(sample.name, str(sample.creation_date)), False)
				continue
			try:
				## Do not remove Quality Control outputs
				if (sample.is_type_fastq_gz_sequencing()):	### illumina
					files_to_remove.append(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True))
					files_to_remove.append(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False))
					#files_to_remove.append(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True))
					#files_to_remove.append(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False))
					# files_to_remove.append(sample.get_abricate_output(TypePath.MEDIA_ROOT))
				else:	### minion
					files_to_remove.append(sample.get_nanofilt_file(TypePath.MEDIA_ROOT))
					#files_to_remove.append(sample.get_rabbitQC_nanofilt(TypePath.MEDIA_ROOT))
					
				if (only_identify_files):
					files_removed = files_to_remove.copy()
				else:
					for path_to_remove in files_to_remove:
						if (self.utils.remove_file(path_to_remove)):
							files_removed.append(path_to_remove)
					
					if (len(files_removed) > 0):
						### save the flag in database
						sample.is_deleted_processed_fastq = True
						sample.date_deleted_processed_fastq = datetime.datetime.now()
						sample.save()
						
			except Exception as e:
				self.out_message("Fail to remove: {}".format(str(e)), False)
				continue
			### remove files
			for file_path in files_removed:
				if (file_path == None): continue
				if (only_identify_files): self.out_message("Identified file: " + file_path, False)
				else: self.out_message("Remove file: " + file_path, False)
				count += 1
		
		if (only_identify_files):
			if (count == 0): self.out_message("None samples were identified", False)
			else: self.out_message("Files identified: {}".format(count), False)
		else:
			if (count == 0): self.out_message("None samples were removed", False)
			else: self.out_message("Files removed: {}".format(count), False)
		return count


	def delete_upload_deleted_files_by_batch_and_not_attached(self, only_identify_files, lst_accounts_to_pass,
															lst_accounts_to_process):
		"""
		delete deleted files that are uploaded by batch
		"""
		self.out_message("\n### Uploaded fastq files by batch\n", False)
		count = 0
		
		upload_files = UploadFiles.objects.all().filter(is_deleted=True, is_deleted_in_file_system=False,
											type_file__name=TypeFile.TYPE_FILE_fastq_gz)
		for upload_file in upload_files:
			files_removed = []
			files_to_remove = []
			
			### test the owner
			if (upload_file.owner.username in lst_accounts_to_pass): continue
			if (len(lst_accounts_to_process) > 0 and \
				not upload_file.owner.username in lst_accounts_to_process): continue
				
			## can be removed already
			try:
				files_to_remove.append(upload_file.get_path_to_file(TypePath.MEDIA_ROOT))
				
				## test the days removed
				if (self.is_out_date_to_delete(upload_file.date_deleted)):
					#self.out_message("Not remove physically: {}; fail the date/remove days test {}.".format(files_to_remove[0], str(upload_file.date_deleted)), False)
					continue

				if (only_identify_files):
					files_removed = files_to_remove.copy()
				else:
					for path_to_remove in files_to_remove:
						if (self.utils.remove_file(path_to_remove)):
							files_removed.append(path_to_remove)
					
					if (len(files_removed) > 0):
						### save the flag in database
						upload_file.is_deleted_in_file_system = True
						upload_file.save()
							
			except Exception as e:
				self.out_message("Fail to remove: {}".format(str(e)), False)
				continue
			for file_path in files_removed:
				if (file_path == None): continue
				if (only_identify_files): self.out_message("Identified file: " + file_path, False)
				else: self.out_message("Remove file: " + file_path, False)
				count += 1
		
		if (only_identify_files):
			if (count == 0): self.out_message("None Uploaded fastq files by batch were identified", False)
			else: self.out_message("Files identified: {}".format(count), False)
		else:
			if (count == 0): self.out_message("None Uploaded fastq files by batch were removed", False)
			else: self.out_message("Files removed: {}".format(count), False)
		self.out_message("### END Uploaded fastq files by batch\n", False)
		return count
	

	def delete_project_samples(self, only_identify_files, lst_accounts_to_pass, lst_accounts_to_process):
		"""
		delete project samples if samples was deleted
		"""
		self.out_message("\n### Project samples\n", False)
		count = 0
		
		software = Software()
		software_names = SoftwareNames()
		project_samples = ProjectSample.objects.all().filter(is_deleted=True, is_deleted_in_file_system=False, sample__is_deleted=True)
		for project_sample in project_samples:
			files_removed = []
			files_to_remove = []

			## test the owner
			if (project_sample.project.owner.username in lst_accounts_to_pass): continue
			if (len(lst_accounts_to_process) > 0 and \
				not project_sample.project.owner.username in lst_accounts_to_process): continue
				
			## can be removed already
			try:
				### files from snippy
				for type_file in software.get_vect_type_files_to_copy(software_names.get_snippy_name()):
					files_to_remove.append(project_sample.get_file_output(TypePath.MEDIA_ROOT, type_file, software_names.get_snippy_name()))
				for type_file in software.get_vect_type_files_to_copy(software_names.get_freebayes_name()):
					files_to_remove.append(project_sample.get_file_output(TypePath.MEDIA_ROOT, type_file, software_names.get_freebayes_name()))
				for type_file in software.get_vect_type_files_to_copy(software_names.get_medaka_name()):
					files_to_remove.append(project_sample.get_file_output(TypePath.MEDIA_ROOT, type_file, software_names.get_medaka_name()))

				## test the days/date removed
				if (self.is_out_date_to_delete(project_sample.date_deleted)):
					#self.out_message("Not remove physically: {}; fail the date/remove days test {}.".format(files_to_remove[0], str(project_sample.date_deleted)), False)
					continue
				
				if (only_identify_files):
					files_removed = files_to_remove.copy()
				else:
					for path_to_remove in files_to_remove:
						if (self.utils.remove_file(path_to_remove)):
							files_removed.append(path_to_remove)
					
					if (len(files_removed) > 0):
						### save the flag in database
						project_sample.is_deleted_in_file_system = True
						project_sample.save()
							
			except Exception as e:
				self.out_message("Fail to remove: {}".format(str(e)), False)
				continue
			for file_path in files_removed:
				if (file_path == None): continue
				if (only_identify_files): self.out_message("Identified file: " + file_path, False)
				else: self.out_message("Remove file: " + file_path, False)
				count += 1
		
		if (only_identify_files):
			if (count == 0): self.out_message("None Uploaded fastq files by batch were identified", False)
			else: self.out_message("Files identified: {}".format(count), False)
		else:
			if (count == 0): self.out_message("None Uploaded fastq files by batch were removed", False)
			else: self.out_message("Files removed: {}".format(count), False)
		self.out_message("### END Project samples\n", False)
		return count




