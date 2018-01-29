#!/usr/bin/env python

import os, subprocess, logging
from utils.utils import Utils
from constants.constants import Constants
from django.conf import settings

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

	def submitte_job(self, file_name):
		"""
		job submission
		raise exception if something wrong
		"""
		sz_temp = os.getcwd()
		os.chdir(os.path.dirname(file_name))	## change dir
		cmd = 'qsub ' + file_name
		proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
		(out, err) = proc.communicate()
		os.chdir(sz_temp)		## change dir
		if (err != None):
			self.logger_production.error('Fail to run: ' + cmd)
			self.logger_debug.error('Fail to run: ' + cmd)
			raise Exception("Fail to submit qsub")
		out_str = out.decode("utf-8")
		b_found = False
		for line in out_str.split('\n'):
			print(line)
			if (line.find("has been submitted") != -1):
				lst_line = line.split(' ')
				if (len(lst_line) > 4 and self.utils.is_integer(lst_line[2])): return int(lst_line[2])
				return None		## don't rise exception... 
		if (not b_found): raise Exception(out_str)


	def set_script_run_sge(self, out_dir, queue_name, vect_cmd, b_remove_out_dir = False, nPriority = 0):
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
			handleSGE.write("#$ -j y\n")  # merge the standard error with standard output
			handleSGE.write("#$ -cwd\n")	# execute the job for the current work directory
			handleSGE.write("#$ -q {}\n".format(queue_name))	# queue name
			if (nPriority > 0): handleSGE.write("#$ -p %d\n" % (nPriority))	# execute the job for the current work directory
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

	##### set different process
	def set_collect_global_files(self, project, user):
		
		vect_command = ['python3 {} collect_global_files --project_id {} --user_id {}'.format(\
				os.path.join(settings.BASE_DIR, 'manage.py'), project.pk, user.pk)]
		out_dir = self.utils.get_temp_dir()
		queue_name = user.profile.queue_name_sge
		if (queue_name == None): queue_name = Constants.QUEUE_SGE_NAME_GLOBAL
		path_file = self.set_script_run_sge(out_dir, queue_name, vect_command, True)
		try:
			sge_id = self.submitte_job(path_file)
		except:
			raise Exception('Fail to submit the job.')
		return sge_id
	
	
		##### set different process
	def set_second_stage_snippy(self, project_sample, user):
		
		vect_command = ['python3 {} second_stage_snippy --project_sample_id {} --user_id {}'.format(\
				os.path.join(settings.BASE_DIR, 'manage.py'), project_sample.pk, user.pk)]
		out_dir = self.utils.get_temp_dir()
		queue_name = user.profile.queue_name_sge
		if (queue_name == None): queue_name = Constants.QUEUE_SGE_NAME_GLOBAL
		path_file = self.set_script_run_sge(out_dir, queue_name, vect_command, True)
		try:
			sge_id = self.submitte_job(path_file)
		except:
			raise Exception('Fail to submit the job.')
		return sge_id


	##### set different process
	def set_run_trimmomatic_species(self, sample, user):
		
		vect_command = ['python3 {} run_trimmomatic_species --sample_id {} --user_id {}'.format(\
				os.path.join(settings.BASE_DIR, 'manage.py'), sample.pk, user.pk)]
		out_dir = self.utils.get_temp_dir()
		path_file = self.set_script_run_sge(out_dir, Constants.QUEUE_SGE_NAME_GLOBAL, vect_command, True)
		try:
			sge_id = self.submitte_job(path_file)
		except:
			raise Exception('Fail to submit the job.')
		return sge_id
		

