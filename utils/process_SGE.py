#!/usr/bin/env python

import os, subprocess, logging
from utils.utils import Utils

# http://www.socher.org/index.php/Main/HowToInstallSunGridEngineOnUbuntu
# https://peteris.rocks/blog/sun-grid-engine-installation-on-ubuntu-server/
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
	
	szScriptSGEName = "launch_lvb.sh"
	utils = Utils()
	
	SGE_JOB_ID_PROCESSING = 1
	SGE_JOB_ID_QUEUE = 2
	SGE_JOB_ID_FINISH = 3
	DEFAULT_QUEUE_NAME = "all.q"
	
	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	def __init__(self):
		pass

	def submit_sge(self, sz_file, file_model, sz_path_directory):
	
		sz_out_dir = "%s/%s" % (sz_path_directory, file_model.key_id)
		if (not os.path.exists(sz_out_dir)):
			sz_cmd = "mkdir -p " + sz_out_dir
			if (os.system(sz_cmd)): return (True, "Error: fail to create project")
		
		sz_cmd = "cp %s %s/infile" % (sz_file, sz_out_dir)
		if (os.system(sz_cmd)): return (True, "Error: fail to copy fasta file")

		### create bash
		sz_cmd = "cd %s; %s -f fasta -s 509739986 -p %d > %s" % (sz_out_dir, self.constants.PATH_LVB_SOFTWARE, self.threadsHead, self.constants.FILE_NAME_OUT_REPORT_LVB)
		self.setScriptRunSGE(sz_out_dir, sz_cmd)
		n_job_id = int(self.submitteJob(sz_out_dir))
		if (n_job_id == 0):
			file_model.job_sge_id = 0 
			return (True, "Error: fail to submit OGE job")
		file_model.is_queue = True
		file_model.is_processing = False
		file_model.job_sge_id = n_job_id
		file_model.path_job = file_model.key_id
		return (False, "")


	def submitte_job(self, file_name):
		"""
		job submission
		raise exception if something wrong
		"""
		sz_temp = os.getcwd()
		os.chdir(os.path.dirname(file_name))	## change dir
		cmd = 'qsub %s/%s' % (file_name)
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
			if (line.find("has been submitted") != -1): 
				b_found = True
				break
		if (not b_found): raise Exception(out_str)


	def set_script_run_sge(self, file_name_out, queue_name, cline, nPriority = 0):
		"""
		create the script to run SGE
		"""
		with open(file_name_out, 'w') as handleSGE:
			handleSGE.write("#!/bin/bash\n")
			handleSGE.write("#$ -V\n")	# Specifies  that  all  environment  variables active
										# within the qsub utility be exported to the context of the job.
			handleSGE.write("#$ -S /bin/bash\n") 	# interpreting shell
			handleSGE.write("#$ -j y\n")  # merge the standard error with standard output
			handleSGE.write("#$ -cwd\n")	# execute the job for the current work directory
			handleSGE.write("#$ -q {}\n".format(queue_name))	# queue name
			if (nPriority > 0): handleSGE.write("#$ -p %d\n" % (nPriority))	# execute the job for the current work directory
			handleSGE.write("\n" + cline)
	
	
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
		szFileResultSGE = self.constants.get_file_name(self.constants.PATH_TEMP, "sge_stat")
		cline = 'qstat > %s' % (szFileResultSGE)
		os.system(cline)
			
		## read the FILE
		handleResult = open(szFileResultSGE)
		vectRunning =[]
		vectWait =[]
		for line in handleResult:
			# pass header and other things
			if (line.find("job-ID") != -1 or len(line) < 3 or line.find("---") == 0): continue
			if (len(line.split()) > 0):
				## jobid is running
				if (line.split()[4] in tagsSGERunning): vectRunning.append(line.split()[0])
				elif (line.split()[4] in tagsSGEWaiting): vectWait.append(line.split()[0])
		handleResult.close()
		
		## remove file
		cline = 'rm %s' % (szFileResultSGE)
		os.system(cline)
		return (vectRunning, vectWait)

	def get_status_process(self, n_SGE_id):
		(vectRunning, vectWait) = self.__get_sge_process__()
		if (str(n_SGE_id) in vectRunning): return self.SGE_JOB_ID_PROCESSING
		if (str(n_SGE_id) in vectWait): return self.SGE_JOB_ID_QUEUE
		return self.SGE_JOB_ID_FINISH
	

		