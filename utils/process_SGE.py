#!/usr/bin/env python

import os, logging
from utils.utils import Utils
from constants.constants import Constants, FileExtensions
from django.conf import settings
from managing_files.models import ProcessControler
from extend_user.models import Profile
from datetime import datetime

# http://www.socher.org/index.php/Main/HowToInstallSunGridEngineOnUbuntu
# https://peteris.rocks/blog/sun-grid-engine-installation-on-ubuntu-server/
# http://biohpc.blogspot.pt/2016/10/sge-installation-of-son-of-grid-engine.html ## centos 7
# http://star.mit.edu/cluster/docs/0.93.3/guides/sge.html		### explain who to use SGE

# /usr/share/gridengine/scripts/init_cluster

#  => SGE_ROOT: /var/lib/gridengine
# => SGE_CELL: default
# => Spool directory: /var/spool/gridengine/spooldb
# => Initial manager user: sgeadmin

## logs
# <qmaster_spool_dir>/messages
# <qmaster_spool_dir>/schedd/messages
# <execd_spool_dir>/<hostname>/messages
# <sge_root>/<sge_cell>/common/accounting
# <sge_root>/<sge_cell>/common/statistics

## default configuration
# /etc/default/gridengine
class ProcessSGE(object):
	
	utils = Utils()
	
	FILE_NAME_SCRIPT_SGE = "launch_job_insa.sh"
	SGE_JOB_ID_PROCESSING = 1
	SGE_JOB_ID_QUEUE = 2
	SGE_JOB_ID_FINISH = 3
	DEFAULT_QUEUE_NAME = "all.q"
	
	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	def __init__(self):
		pass

	###########################################
	###
	###		IMPORTANT
	###			Put qsub in /usr/bin/qsub
	###
	def submitte_job(self, file_name):
		"""
		job submission
		raise exception if something wrong
		"""
		temp_file = self.utils.get_temp_file('qsub_out', FileExtensions.FILE_TXT)
		cmd = 'export SGE_ROOT={}; qsub {} > {}'.format(settings.SGE_ROOT, file_name, temp_file)
		exist_status = os.system(cmd)
		if (exist_status != 0):
			if (os.path.exists(temp_file)): os.unlink(temp_file)
			self.logger_production.error('Fail to run: ' + cmd + " - exit code: " + str(exist_status))
			self.logger_debug.error('Fail to run: ' + cmd + " - exit code: " + str(exist_status))
			raise Exception("Fail to submit qsub")
		## read output
		vect_out = self.utils.read_text_file(temp_file)
		if (os.path.exists(temp_file)): os.unlink(temp_file)
		b_found = False
		for line in vect_out:
			if (line.find("has been submitted") != -1):
				lst_line = line.split(' ')
				if (len(lst_line) > 4 and self.utils.is_integer(lst_line[2])): return int(lst_line[2])
				return None		## don't rise exception... 
		if (not b_found): raise Exception("\n".join(vect_out))


	def set_script_run_sge(self, out_dir, queue_name, vect_cmd, job_name, b_remove_out_dir = False, job_name_wait = "", nPriority = 0):
		"""
		create the script to run SGE
		"""
		if (len(vect_cmd) == 0): return None
		
		file_name_out = os.path.join(out_dir, ProcessSGE.FILE_NAME_SCRIPT_SGE)
		with open(file_name_out, 'w') as handleSGE:
			handleSGE.write("#!/bin/bash\n")
			handleSGE.write("#$ -V\n")	# Specifies  that  all  environment  variables active
										# within the qsub utility be exported to the context of the job.
			handleSGE.write("#$ -S /bin/bash\n") 	# interpreting shell
			if len(job_name_wait) > 0: handleSGE.write("#$ -hold_jid {}\n".format(job_name_wait))	# need to wait until all this jobs names finished
			handleSGE.write("#$ -j y\n")	# merge the standard error with standard output
			handleSGE.write("#$ -N {}\n".format(job_name))	# job name
			handleSGE.write("#$ -cwd\n")	# execute the job for the current work directory
			handleSGE.write("#$ -q {}\n".format(queue_name))	# queue name
			handleSGE.write("#$ -o {}\n".format(out_dir))		# out path file
			if (nPriority > 0): handleSGE.write("#$ -p {}\n".format(nPriority))	# execute the job for the current work directory
			for cline in vect_cmd: handleSGE.write("\n" + cline)
			if (b_remove_out_dir):
				handleSGE.write("\nif [ $? -eq 0 ]\nthen\n  rm -r {}\nfi\n".format(out_dir))
		return file_name_out
	
	def __get_sge_process__(self):
		"""
		#Job status - one of
	
		### test if all jobs submitted to the SGE are finish
		## return 0, if is end
		## return -1, error
		## other value, keeping running
		## also returns a vector with jobId already finish
		
		#	* d(eletion)
		#	* E(rror)
		#	* h(old)
		#	* r(unning)
		#	* R(estarted)
		#	* s(uspended),
		#	* S(uspended)
		#	* t(ransfering)
		#	* T(hreshold)
		#	* w(aiting)
		"""
		tagsSGERunning = ('r', 't')
		tagsSGEWaiting = ('qw', 'w')
		# test with qstat
		file_result = self.utils.get_temp_file('sge_stat', '.txt')
		cline = 'qstat > %s' % (file_result)
		os.system(cline)
			
		## read the FILE
		with open(file_result) as handle_result:
			vectRunning =[]
			vectWait =[]
			for line in handle_result:
				# pass header and other things
				if (line.find("job-ID") != -1 or len(line) < 3 or line.find("---") == 0): continue
				if (len(line.split()) > 0):
					## jobid is running
					if (line.split()[4] in tagsSGERunning): vectRunning.append(line.split()[0])
					elif (line.split()[4] in tagsSGEWaiting): vectWait.append(line.split()[0])
		
		## remove file
		if (os.path.exists(file_result)): os.unlink(file_result)
		return (vectRunning, vectWait)

	def get_status_process(self, n_SGE_id):
		(vectRunning, vectWait) = self.__get_sge_process__()
		if (str(n_SGE_id) in vectRunning): return self.SGE_JOB_ID_PROCESSING
		if (str(n_SGE_id) in vectWait): return self.SGE_JOB_ID_QUEUE
		return self.SGE_JOB_ID_FINISH
	
	def is_finished(self, n_SGE_id):
		"""
		is it finished
		"""
		return self.get_status_process(n_SGE_id) == self.SGE_JOB_ID_FINISH

	##### set collect global files
	def set_collect_global_files(self, project, user):
		"""
		job_name = "job_name_<user_id>_<seq_id>"
		only run this task after all second_stage_snippy
		"""
		process_controler = ProcessControler()
		vect_command = ['python3 {} collect_global_files --project_id {} --user_id {}'.format(\
				os.path.join(settings.BASE_DIR, 'manage.py'), project.pk, user.pk)]
		self.logger_production.info('Processing: ' + ";".join(vect_command))
		self.logger_debug.info('Processing: ' + ";".join(vect_command))
		out_dir = self.utils.get_temp_dir()
		
		queue_name = user.profile.queue_name_sge
		if (queue_name == None): queue_name = Constants.QUEUE_SGE_NAME_GLOBAL
		
		(job_name_wait, job_name) = user.profile.get_name_sge_seq(Profile.SGE_GLOBAL)
		path_file = self.set_script_run_sge(out_dir, queue_name, vect_command, job_name, True, job_name_wait)
		try:
			sge_id = self.submitte_job(path_file)
			if (sge_id != None): self.set_process_controlers(user, process_controler.get_name_project(project), sge_id)
		except:
			raise Exception('Fail to submit the job.')
		return sge_id
	
	##### set collect global files
	def set_collect_global_files_for_update_metadata(self, project, user):
		"""
		job_name = "job_name_<user_id>_<seq_id>"
		only run this task after all second_stage_snippy
		"""
		process_controler = ProcessControler()
		vect_command = ['python3 {} collect_global_files_for_update_metadata --project_id {} --user_id {}'.format(\
				os.path.join(settings.BASE_DIR, 'manage.py'), project.pk, user.pk)]
		self.logger_production.info('Processing: ' + ";".join(vect_command))
		self.logger_debug.info('Processing: ' + ";".join(vect_command))
		out_dir = self.utils.get_temp_dir()
		
		queue_name = user.profile.queue_name_sge
		if (queue_name == None): queue_name = Constants.QUEUE_SGE_NAME_GLOBAL
		(job_name_wait, job_name) = user.profile.get_name_sge_seq(Profile.SGE_GLOBAL)
		path_file = self.set_script_run_sge(out_dir, queue_name, vect_command, job_name, True, job_name_wait)
		try:
			sge_id = self.submitte_job(path_file)
			if (sge_id != None): self.set_process_controlers(user, process_controler.get_name_project(project), sge_id)
		except:
			raise Exception('Fail to submit the job.')
		return sge_id
	
	
	def set_second_stage_snippy(self, project_sample, user, job_name, job_name_wait):
		"""
		can make several in parallel but only after last collect_global_files
		"""
		process_controler = ProcessControler()
		vect_command = ['python3 {} second_stage_snippy --project_sample_id {} --user_id {}'.format(\
				os.path.join(settings.BASE_DIR, 'manage.py'), project_sample.pk, user.pk)]
		self.logger_production.info('Processing: ' + ";".join(vect_command))
		self.logger_debug.info('Processing: ' + ";".join(vect_command))
		out_dir = self.utils.get_temp_dir()
		queue_name = user.profile.queue_name_sge
		if (queue_name == None): queue_name = Constants.QUEUE_SGE_NAME_GLOBAL
		path_file = self.set_script_run_sge(out_dir, queue_name, vect_command, job_name, True, job_name_wait)
		try:
			sge_id = self.submitte_job(path_file)
			if (sge_id != None): self.set_process_controlers(user, process_controler.get_name_project_sample(project_sample), sge_id)
		except:
			raise Exception('Fail to submit the job.')
		return sge_id


	def set_run_trimmomatic_species(self, sample, user, job_name = "job_name_to_run"):
		"""
		Run trimmomatic and identify species 
		Can run free without wait for anything
		"""
		process_controler = ProcessControler()
		vect_command = ['python3 {} run_trimmomatic_species --sample_id {} --user_id {}'.format(\
				os.path.join(settings.BASE_DIR, 'manage.py'), sample.pk, user.pk)]
		self.logger_production.info('Processing: ' + ";".join(vect_command))
		self.logger_debug.info('Processing: ' + ";".join(vect_command))
		out_dir = self.utils.get_temp_dir()
		path_file = self.set_script_run_sge(out_dir, Constants.QUEUE_SGE_NAME_GLOBAL, vect_command, job_name, True)
		try:
			sge_id = self.submitte_job(path_file)
			if (sge_id != None): self.set_process_controlers(user, process_controler.get_name_sample(sample), sge_id)
		except:
			raise Exception('Fail to submit the job.')
		return sge_id
		
	def set_link_files(self, user):
		"""
		only can run one at a time
		"""
		process_controler = ProcessControler()
		vect_command = ['python3 {} link_files --user_id {}'.format(\
				os.path.join(settings.BASE_DIR, 'manage.py'), user.pk)]
		self.logger_production.info('Processing: ' + ";".join(vect_command))
		self.logger_debug.info('Processing: ' + ";".join(vect_command))
		out_dir = self.utils.get_temp_dir()
		
		(job_name_wait, job_name) = user.profile.get_name_sge_seq(Profile.SGE_LINK)
		path_file = self.set_script_run_sge(out_dir, Constants.QUEUE_SGE_NAME_FAST, vect_command, job_name, True, job_name_wait)
		try:
			sge_id = self.submitte_job(path_file)
			if (sge_id != None): self.set_process_controlers(user, process_controler.get_name_link_files_user(user), sge_id)
		except:
			raise Exception('Fail to submit the job.')
		return sge_id

	##### set read sample file 	### ultra fast queue
	def set_read_sample_file(self, upload_files, user):
		
		process_controler = ProcessControler()
		vect_command = ['python3 {} read_sample_file --upload_files_id {} --user_id {}'.format(\
				os.path.join(settings.BASE_DIR, 'manage.py'), upload_files.pk, user.pk)]
		self.logger_production.info('Processing: ' + ";".join(vect_command))
		self.logger_debug.info('Processing: ' + ";".join(vect_command))
		out_dir = self.utils.get_temp_dir()
		(job_name_wait, job_name) = user.profile.get_name_sge_seq(Profile.SGE_LINK)
		path_file = self.set_script_run_sge(out_dir, Constants.QUEUE_SGE_NAME_FAST, vect_command, job_name, True, job_name_wait)
		try:
			sge_id = self.submitte_job(path_file)
			if (sge_id != None): self.set_process_controlers(user, process_controler.get_name_upload_files(upload_files), sge_id)
		except:
			raise Exception('Fail to submit the job.')
		return sge_id

	def set_read_sample_file_with_metadata(self, upload_files, user, job_name = "job_name"):
		"""
		update metadata, normal queue, wait for all of other data
		"""
		process_controler = ProcessControler()
		vect_command = ['python3 {} update_metadata_sample_file --upload_files_id {} --user_id {}'.format(\
				os.path.join(settings.BASE_DIR, 'manage.py'), upload_files.pk, user.pk)]
		self.logger_production.info('Processing: ' + ";".join(vect_command))
		self.logger_debug.info('Processing: ' + ";".join(vect_command))
		out_dir = self.utils.get_temp_dir()
		queue_name = user.profile.queue_name_sge
		if (queue_name == None): queue_name = Constants.QUEUE_SGE_NAME_GLOBAL
		(job_name_wait, job_name) = user.profile.get_name_sge_seq(Profile.SGE_LINK)
		path_file = self.set_script_run_sge(out_dir, queue_name, vect_command, job_name, True, job_name_wait)
		try:
			sge_id = self.submitte_job(path_file)
			if (sge_id != None): self.set_process_controlers(user, process_controler.get_name_upload_files(upload_files), sge_id)
		except:
			raise Exception('Fail to submit the job.')
		return sge_id
	
	### only for tests
	def submit_dummy_sge(self, job_name = "job_name"):
		"""
		only for tests
		"""
		vect_command = ['echo "#####start" >> /tmp/sge.out', 'echo $HOSTNAME >> /tmp/sge.out', 'echo "start waiting" >> /tmp/sge.out',\
					'echo $PATH >> /tmp/sge.out', 'echo $SGE_ROOT >> /tmp/sge.out', 
					'date >> /tmp/sge.out', 'sleep 2s', 'date >> /tmp/sge.out',\
					'echo "end" >> /tmp/sge.out']
		out_dir = self.utils.get_temp_dir()
		path_file = self.set_script_run_sge(out_dir, Constants.QUEUE_SGE_NAME_FAST, vect_command, job_name, False)
		os.system("/bin/sh {}".format(path_file))
		try:
			sge_id = self.submitte_job(path_file)
		except:
			raise Exception('Fail to submit the job.')
		return sge_id
	
	def set_process_controlers(self, user, name_of_process, name_sge_id):
		"""
		Add a record in ProcessControlers
		"""
		process_controler = ProcessControler()
		process_controler.owner = user
		process_controler.name = name_of_process
		process_controler.name_sge_id = name_sge_id
		process_controler.save()

	
	def set_process_controler(self, user, name_of_process, flags):
		"""
		name_of_process:
			process_controler.get_name_upload_files(upload_files),
			process_controler.get_name_link_files_user(user),
			process_controler.get_name_sample(sample),
			process_controler.get_name_project(project), sge_id)
			process_controler.get_name_project_sample(project_sample)
			
		flags: ProcessControler.FLAG_FINISHED, ProcessControler.FLAG_RUNNING, ProcessControler.FLAG_ERROR
		"""
		
		if (flags == ProcessControler.FLAG_FINISHED):
			data_set = ProcessControler.objects.filter(owner__id=user.pk, name=name_of_process, is_running=True, is_finished=False, is_error=False)
		elif (flags == ProcessControler.FLAG_ERROR):
			data_set = ProcessControler.objects.filter(owner__id=user.pk, name=name_of_process, is_finished=False, is_error=False)
		else:
			data_set = ProcessControler.objects.filter(owner__id=user.pk, name=name_of_process, is_running=False, is_finished=False, is_error=False)
		
		if (data_set.count() > 0):
			process_controler = ProcessControler.objects.get(pk=data_set[0].pk)
			if (flags == ProcessControler.FLAG_FINISHED):
				process_controler.is_finished = True
				process_controler.is_running = False
				process_controler.close_date = datetime.now()
			elif (flags == ProcessControler.FLAG_ERROR):
				process_controler.is_finished = True
				process_controler.is_error = True
				process_controler.is_running = False
				process_controler.close_date = datetime.now()
			elif (flags == ProcessControler.FLAG_RUNNING):
				process_controler.is_running = True
			process_controler.save()

