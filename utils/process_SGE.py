#!/usr/bin/env python

import logging
import os
import subprocess
import time
from datetime import datetime
from typing import List, Optional

from django.conf import settings
from django.contrib.auth.models import User
from django.db import transaction

from constants.constants import Constants, FileExtensions, TypePath
from extend_user.models import Profile
from managing_files.models import ProcessControler
from pathogen_identification.constants_settings import \
    ConstantsSettings as PICS
from utils.utils import Utils


class ProcessSched(object):
    utils: Utils = Utils()

    FILE_NAME_SCRIPT_SLURM = "launch_job_insa.sh"
    JOB_ID_PROCESSING = 1
    JOB_ID_QUEUE = 2
    JOB_ID_FINISH = 3
    DEFAULT_QUEUE_NAME = "all.q"

    ## logging
    logger_debug = logging.getLogger("fluWebVirus.debug")
    logger_production = logging.getLogger("fluWebVirus.production")

    def __init__(self):
        pass

    ###########################################
    ###
    ###        BEGIN main methods
    ###

    def submit_job(self, file_name):
        """
        job submission
        raise exception if something wrong
        """
        temp_file = self.utils.get_temp_file("qsub_out", FileExtensions.FILE_TXT)
        cmd = ""
        # Create a bash script to deploy the SLURM job
        # bash_script_path = "/data/tmp/submit_slurm_job.sh"
        bash_script_path = self.utils.get_temp_file(
            "submit_slurm_job", FileExtensions.FILE_BASH_SCRIPT
        )
        with open(bash_script_path, "w") as bash_script:
            bash_script.write("#!/bin/bash\n")
            bash_script.write("cd /data/tmp\n")
            bash_script.write("sbatch {} > {}\n".format(file_name, temp_file))
            bash_script.write("cd /insaflu_web/INSaFLU\n")

        # Make the script executable
        os.chmod(bash_script_path, 0o755)

        # Execute the bash script using subprocess
        cmd = bash_script_path
        process = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        exist_status = process.returncode

        # Check if error occurred
        if exist_status != 0:
            print("Error: ", stderr.decode())
            print("stdout: ", stdout.decode())
            print("Error: ", exist_status)
            print("cmd: ", cmd)

            ## remove file
            if os.path.exists(temp_file):
                os.unlink(temp_file)
            self.logger_production.error(
                "Fail to run: " + cmd + " - exit code: " + str(exist_status)
            )
            self.logger_debug.error(
                "Fail to run: " + cmd + " - exit code: " + str(exist_status)
            )
            raise Exception("Fail to submit qsub")

        ## read output

        b_found = False
        tries = 3
        time_wait_past_first = 3

        def read_file_once(temp_file) -> Optional[int]:
            vect_out = self.utils.read_text_file(temp_file)

            for line in vect_out:
                if line.find("Submitted batch job") != -1:
                    lst_line = line.split(" ")
                    if len(lst_line) > 2 and self.utils.is_integer(lst_line[3]):
                        return int(lst_line[3])
                    return None  ## don't rise exception...
        
        submitted = read_file_once(temp_file)
        if submitted is None:
            while tries > 0:
                time.sleep(time_wait_past_first)
                submitted = read_file_once(temp_file)
                if submitted is not None:
                    b_found = True
                    break
                tries -= 1
        else:
            b_found = True
        #if os.path.exists(temp_file):
        #    os.unlink(temp_file)

        if not b_found:
            raise Exception("Fail to submit job")
        
        return submitted

    def collect_jobname_jobid_slurm(self, job_name):
        """
        collect the job id from slurm squeue using the job name.
        if job_name is not found return None
        """
        cmd = "squeue -o '%j %i' | grep {}".format(job_name)

        process = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        exist_status = process.returncode

        if exist_status != 0 and stderr.decode().strip() != "": # counts as error if nothing is returned
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            return None

        output_lines = stdout.decode().strip().split("\n")
        #
        # JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
        for line in output_lines:
            parts = line.split()
            if len(parts) >= 2:
                job_name_found = parts[0]
                job_id = parts[1]
                if job_name_found == job_name:
                    return job_id

        return None

    def jobwait_ids(self, job_names_wait: list):
        ids = [
            self.collect_jobname_jobid_slurm(job_name) for job_name in job_names_wait
        ]
        return [job_id for job_id in ids if job_id is not None]

    def set_script_run_slurm(
        self,
        out_dir,
        queue_name,
        vect_cmd,
        job_name,
        b_remove_out_dir=False,
        job_name_wait=[],
        alternative_temp_dir=None,
        cpus = 1,
        memory = "4G"
    ):
        """
        create the script to run the job in slurm
        """
        b_remove_out_dir = False
        if len(vect_cmd) == 0:
            return None
        # b_remove_out_dir = False
        file_name_out = os.path.join(out_dir, ProcessSched.FILE_NAME_SCRIPT_SLURM)

        with open(file_name_out, "w") as handleSLURM:
            handleSLURM.write("#!/bin/bash\n")
            handleSLURM.write(
                "#SBATCH --export=ALL\n"
            )  # Specifies  that  all  environment  variables active
            ##
            handleSLURM.write("#SBATCH --cpus-per-task={}\n".format(cpus))
            handleSLURM.write("#SBATCH --mem={}\n".format(memory))
            # within the qsub utility be exported to the context of the job.
            # handleSLURM.write("#$ -S /bin/bash\n")  # interpreting shell
            ## hold_jid <comma separated list of job-ids, can also be a job id pattern such as 2722*> :
            ## will start the current job/job -array only after completion of all jobs in the comma separated list
            if isinstance(job_name_wait, str):
                job_name_wait = [job_name_wait]
            if len(job_name_wait) > 0:
                job_wait_ids = self.jobwait_ids(job_name_wait)
                if len(job_wait_ids) > 0:
                    handleSLURM.write(
                        "#$ --dependency=afterok:{}\n".format(",".join(job_wait_ids))
                    )  # need to wait until all this jobs names finished

            # handleSLURM.write(
            #    "#$ -j y\n"
            # )  # merge the standard error with standard output
            handleSLURM.write("#SBATCH -J {}\n".format(job_name))  # job name
            # handleSLURM.write(
            #    "#SBATCH --partition={}\n".format(queue_name)
            # )  # queue name
            # handleSLURM.write("#$ --output={}\n".format(out_dir))  # out path file
            # handleSLURM.write("#SBATCH --error={}\n".format(out_dir))
            handleSLURM.write("#SBATCH --ntasks=1\n")
            handleSLURM.write("#SBATCH --output={}/slurm_%j.out\n".format(out_dir))
            handleSLURM.write("#SBATCH --begin=now\n")
            handleSLURM.write("\n")

            for cline in vect_cmd:
                handleSLURM.write("\n" + cline)
            if b_remove_out_dir and not settings.RUN_TEST_IN_COMMAND_LINE:
                handleSLURM.write(
                    "\nif [ $? -eq 0 ]\nthen\n  rm -r {}\nfi\n".format(out_dir)
                )
                if alternative_temp_dir is not None:
                    handleSLURM.write(
                        "\nif [ $? -eq 0 ]\nthen\n  rm -r {}\nfi\n".format(
                            alternative_temp_dir
                        )
                    )

        return file_name_out

    def __get_slurm_process__(self):
        """
        #Job status - one of

        ### test if all jobs submitted are finished, if not return the job id that are still running or in waiting queue
        ## return 0, if is end
        ## return -1, error
        ## other value, keeping running
        ## also returns a vector with jobId already finish

        #    * d(eletion)
        #    * E(rror)
        #    * h(old)
        #    * r(unning)
        #    * R(estarted)
        #    * s(uspended),
        #    * S(uspended)
        #    * t(ransfering)
        #    * T(hreshold)
        #    * w(aiting)
        """
        tagsRunning = ("r", "t")
        tagsWaiting = ("hqw", "qw", "w")
        # test with squeue
        file_result = self.utils.get_temp_file("slurm_stat", ".txt")
        cline = "squeue > {}".format(file_result)
        os.system(cline)

        ## read the FILE
        with open(file_result) as handle_result:
            vectRunning = []
            vectWait = []
            for line in handle_result:
                # pass header and other things
                if line.find("job-ID") != -1 or len(line) < 3 or line.find("---") == 0:
                    continue
                if len(line.split()) > 0:
                    ## jobid is running
                    if line.split()[4] in tagsRunning:
                        vectRunning.append(line.split()[0])
                    elif line.split()[4] in tagsWaiting:
                        vectWait.append(line.split()[0])

        ## remove file
        if os.path.exists(file_result):
            os.unlink(file_result)
        return (vectRunning, vectWait)

    def get_status_process(self, n_job_id):
        (vectRunning, vectWait) = self.__get_slurm_process__()
        if str(n_job_id) in vectRunning:
            return self.JOB_ID_PROCESSING
        if str(n_job_id) in vectWait:
            return self.JOB_ID_QUEUE
        return self.JOB_ID_FINISH

    def is_finished(self, n_job_id):
        """
        is it finished
        """
        return self.get_status_process(n_job_id) == self.JOB_ID_FINISH

    def exists_taks_running(self):
        """
        test if there any tasks running...
        """
        file_result = self.utils.get_temp_file("slurm_stat", ".txt")
        cline = "squeue > {}".format(file_result)
        os.system(cline)
        ## read the FILE
        with open(file_result) as handle_result:
            for line in handle_result:
                if line.find("job-ID") != -1 or len(line) < 3 or line.find("---") == 0:
                    continue
                if len(line.strip()) > 0:
                    if os.path.exists(file_result):
                        os.unlink(file_result)
                    return True
        if os.path.exists(file_result):
            os.unlink(file_result)
        return False

    def _get_prefix_in_wait_queue(self, prefix_id):
        """
        check if a predefine ID is in the queue
        if in the waiting queue need to have this:
        scheduling info:            job dropped because of job dependencies

        """
        file_result = self.utils.get_temp_file("slurm_stat", ".txt")
        cline = "squeue -j {}* > {}".format(
            prefix_id, file_result
        )
        os.system(cline)
        ## read the FILE
        vect_job_ids = []
        job_candicate = None
        with open(file_result) as handle_result:
            for line in handle_result:
                if line.find("=========") == 0:
                    job_candicate = None
                elif line.find("job_number:") == 0:
                    job_candicate = line.split(":")[1].strip()
                elif (
                    line.find("scheduling info:") == 0 and not job_candicate is None
                ):  ## it is wainting in the queue
                    vect_job_ids.append(job_candicate)
        if os.path.exists(file_result):
            os.unlink(file_result)
        if len(vect_job_ids) > 0:
            return ",".join(vect_job_ids)
        return None

    def wait_until_finished(self, vect_slurm_ids):
        """
        wait till all end
        if len(vect_job_ids) == 0 wait till all are finished, doesn't matter the ID
        """
        ## expand vect_job_idsbecaus some lines can have morethan one ID
        vect_slurm_to_search = [c for b in vect_slurm_ids for c in str(b).split(",")]
        if len(vect_slurm_to_search) == 0:
            while self.exists_taks_running():
                print("=" * 50 + "\n  waiting for slurm\n" + str(datetime.now()))
                time.sleep(5)  ## wais 5 seconds
        else:
            while len(vect_slurm_to_search) > 0:
                print("=" * 50)
                print(
                    "   wait for these ids: {}".format(
                        ";".join([str(_) for _ in vect_slurm_to_search])
                    )
                )
                vect_remove = []
                for slurm_id in vect_slurm_to_search:
                    if self.is_finished(slurm_id):
                        vect_remove.append(slurm_id)

                ### remove 
                for slurm_id in vect_remove:
                    vect_slurm_to_search.remove(slurm_id)

                print("=" * 50)
                if len(vect_slurm_to_search) > 0:
                    time.sleep(5)  ## wais 5 seconds
        ## set the one still running
        return [int(_) for _ in vect_slurm_to_search]

    #### END MAIN files
    #############
    #############

    ##### set collect global files
    def set_collect_global_files(self, project, user):
        """
        job_name = "job_name_<user_id>_<seq_id>"
        only run this task after all second_stage_snippy
        """
        process_controler = ProcessControler()
        vect_command = [
            "python3 {} collect_global_files --project_id {} --user_id {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"), project.pk, user.pk
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()

        queue_name = user.profile.queue_name_slurm
        if queue_name == None:
            queue_name = Constants.QUEUE_NAME_GLOBAL

        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_projects, Constants.PROCESS_GLOBAL
        )
        path_file = self.set_script_run_slurm(
            out_dir, queue_name, vect_command, job_name, True, [job_name_wait],
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_project(project), job_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    ##### set collect global files
    def set_collect_global_files_for_update_metadata(self, project, user):
        """
        job_name = "job_name_<user_id>_<seq_id>"
        only run this task after all second_stage_snippy
        """
        process_controler = ProcessControler()
        vect_command = [
            "python3 {} collect_global_files_for_update_metadata --project_id {} --user_id {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"), project.pk, user.pk
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()

        queue_name = user.profile.queue_name_slurm
        if queue_name == None:
            queue_name = Constants.QUEUE_NAME_GLOBAL
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_dont_care, Constants.PROCESS_GLOBAL
        )
        path_file = self.set_script_run_slurm(
            out_dir, queue_name, vect_command, job_name, True, [job_name_wait]
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_project(project), job_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    ##### set collect global files
    def set_collect_update_pangolin_lineage(self, project, user):
        """
        job_name = "job_name_<user_id>_<seq_id>"
        only run this task after all second_stage_snippy
        """
        process_controler = ProcessControler()
        vect_command = [
            "python3 {} collect_update_pangolin_lineage --project_id {} --user_id {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"), project.pk, user.pk
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()

        queue_name = user.profile.queue_name_slurm
        if queue_name == None:
            queue_name = Constants.QUEUE_NAME_GLOBAL
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_dont_care, Constants.PROCESS_GLOBAL
        )
        path_file = self.set_script_run_slurm(
            out_dir, queue_name, vect_command, job_name, True, [job_name_wait]
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_project(project), job_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def set_collect_update_mutation_report(self, project, user):
        """
        job_name = "job_name_<user_id>_<seq_id>"
        only run this task after all second_stage_snippy
        """
        process_controler = ProcessControler()
        vect_command = [
            "python3 {} collect_update_mutation_report --project_id {} --user_id {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"), project.pk, user.pk
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()

        queue_name = user.profile.queue_name_slurm
        if queue_name == None:
            queue_name = Constants.QUEUE_NAME_GLOBAL
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_dont_care, Constants.PROCESS_GLOBAL
        )
        path_file = self.set_script_run_slurm(
            out_dir, queue_name, vect_command, job_name, True, [job_name_wait]
        )
        try:
            job_id = self.submit_job(path_file)

            if job_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_project(project), job_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def set_second_stage_snippy(
        self, project_sample, user, job_name, vect_job_name_wait
    ):
        """
        can make several in parallel but only after last collect_global_files
        """
        process_controler = ProcessControler()
        vect_command = [
            "python3 {} second_stage_snippy --project_sample_id {} --user_id {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                project_sample.pk,
                user.pk,
                (
                    "--settings fluwebvirus.settings_test"
                    if settings.RUN_TEST_IN_COMMAND_LINE
                    else ""
                ),
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()
        queue_name = user.profile.queue_name_slurm
        if queue_name == None:
            queue_name = Constants.QUEUE_NAME_GLOBAL
        path_file = self.set_script_run_slurm(
            out_dir, queue_name, vect_command, job_name, True, vect_job_name_wait
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user,
                    process_controler.get_name_project_sample(project_sample),
                    job_id,
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def set_second_stage_medaka(
        self, project_sample, user, job_name, vect_job_name_wait
    ):
        """
        can make several in parallel but only after last collect_global_files
        """
        process_controler = ProcessControler()
        vect_command = [
            "python3 {} second_stage_medaka --project_sample_id {} --user_id {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                project_sample.pk,
                user.pk,
                (
                    "--settings fluwebvirus.settings_test"
                    if settings.RUN_TEST_IN_COMMAND_LINE
                    else ""
                ),
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()
        queue_name = user.profile.queue_name_slurm
        if queue_name == None:
            queue_name = Constants.QUEUE_NAME_GLOBAL
        path_file = self.set_script_run_slurm(
            out_dir, queue_name, vect_command, job_name, True, vect_job_name_wait
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user,
                    process_controler.get_name_project_sample(project_sample),
                    job_id,
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def _remove_files_create_by_nanofilt_and_stat(self, sample):
        """
        remove all files create by run_nanofilt_and_stat
        """
        self.utils.remove_file(sample.get_rabbitQC_output(TypePath.MEDIA_ROOT))
        self.utils.remove_file(sample.get_nanofilt_file(TypePath.MEDIA_ROOT))
        self.utils.remove_file(sample.get_rabbitQC_nanofilt(TypePath.MEDIA_ROOT))

    def _remove_files_create_by_fastq_and_trimmomatic(self, sample):
        """
        remove all files create by run_nanofilt_and_stat
        """
        self.utils.remove_file(sample.get_fastqc_output(TypePath.MEDIA_ROOT, True))
        self.utils.remove_file(sample.get_fastqc_output(TypePath.MEDIA_ROOT, False))
        self.utils.remove_file(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True))
        self.utils.remove_file(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False))
        self.utils.remove_file(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True))
        self.utils.remove_file(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False))

    def _remove_files_create_by_identify_type_and_sub_type(self, sample):
        """
        remove all files create by run_nanofilt_and_stat
        """
        self.utils.remove_file(sample.get_abricate_output(TypePath.MEDIA_ROOT))
        self.utils.remove_file(sample.get_draft_contigs_output(TypePath.MEDIA_ROOT))
        self.utils.remove_file(
            sample.get_draft_contigs_abricate_output(TypePath.MEDIA_ROOT)
        )

    def set_run_trimmomatic_species(self, sample, user, job_name="job_name_to_run"):
        """
        Run trimmomatic and identify species
        Can run free without wait for anything
        """
        process_controler = ProcessControler()
        vect_command = [
            "python3 {} run_trimmomatic_species --sample_id {} --user_id {} {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                sample.pk,
                user.pk,
                (
                    "--settings fluwebvirus.settings_test"
                    if settings.RUN_TEST_IN_COMMAND_LINE
                    else ""
                ),
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()
        try:
            path_file = self.set_script_run_slurm(
                out_dir, Constants.QUEUE_NAME_GLOBAL, vect_command, job_name, True
            )
            job_id = self.submit_job(path_file)
        except Exception as e:
            print("Error: ", e)
            raise Exception("Fail to submit the job.")

        try:
            self._remove_files_create_by_identify_type_and_sub_type(sample)
            self._remove_files_create_by_fastq_and_trimmomatic(sample)
            if job_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_sample(sample), job_id
                )

            ### change flag to not finished
            sample.is_ready_for_projects = False
            sample.is_sample_in_the_queue = True
            sample.save()
        except:
            print("Error: ", e)
            raise Exception("Fail to submit the job.")
        return job_id

    def set_run_clean_minion(self, sample, user, job_name="job_name_to_run"):
        """
        Run trimmomatic and identify species
        Can run free without wait for anything
        """
        process_controler = ProcessControler()
        vect_command = [
            "python3 {} run_clean_minion --sample_id {} --user_id {} {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                sample.pk,
                user.pk,
                (
                    "--settings fluwebvirus.settings_test"
                    if settings.RUN_TEST_IN_COMMAND_LINE
                    else ""
                ),
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()
        path_file = self.set_script_run_slurm(
            out_dir, Constants.QUEUE_NAME_GLOBAL, vect_command, job_name, True
        )
        try:
            self._remove_files_create_by_identify_type_and_sub_type(sample)
            self._remove_files_create_by_nanofilt_and_stat(sample)
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_sample(sample), job_id
                )

            ### change flag to not finished
            sample.is_ready_for_projects = False
            sample.is_sample_in_the_queue = True
            sample.save()
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def set_link_files(self, user, b_test=False):
        """
        only can run one at a time
        """
        process_controler = ProcessControler()
        vect_command = [
            "python3 {} link_files --user_id {} {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                user.pk,
                "--settings fluwebvirus.settings_test" if b_test else "",
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()

        prefix_to_find = user.profile.get_prefix_name(
            Constants.PROCESS_link_files, Constants.PROCESS_LINK
        )
        job_id = self._get_prefix_in_wait_queue(prefix_to_find)
        if job_id is None:  ## if prefix does not exist in queue need to submit new one
            (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
                Constants.PROCESS_link_files, Constants.PROCESS_LINK
            )
            path_file = self.set_script_run_slurm(
                out_dir,
                Constants.QUEUE_NAME_FAST,
                vect_command,
                job_name,
                True,
                [job_name_wait],
            )
            try:
                job_id = self.submit_job(path_file)
                if job_id != None:
                    self.set_process_controlers(
                        user, process_controler.get_name_link_files_user(user), job_id
                    )
            except:
                raise Exception("Fail to submit the job.")
        return job_id

    def set_create_sample_list_by_user(self, user, vect_job_wait, b_test=False):
        """
        only can run one at a time
        :param user
        :param vect_job_wait wait for the job of the sample
        """
        process_controler = ProcessControler()
        vect_command = [
            "python3 {} create_sample_list_by_user --user_id {} {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                user.pk,
                "--settings fluwebvirus.settings_test" if b_test else "",
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()

        ## same has trimmomatic and minion
        prefix_to_find = user.profile.get_prefix_name(
            Constants.PROCESS_collect_all_samples, Constants.PROCESS_REGULAR
        )
        job_id = self._get_prefix_in_wait_queue(prefix_to_find)
        if job_id is None:  ## if prefix does not exist in queue need to submit new one
            (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
                Constants.PROCESS_collect_all_samples, Constants.PROCESS_REGULAR
            )
            path_file = self.set_script_run_slurm(
                out_dir,
                Constants.QUEUE_NAME_FAST,
                vect_command,
                job_name,
                True,
                vect_job_wait,
            )
            try:
                job_id = self.submit_job(path_file)
                if job_id != None:
                    self.set_process_controlers(
                        user,
                        process_controler.get_name_collect_all_samples_user(user),
                        job_id,
                    )
            except:
                raise Exception("Fail to submit the job.")
        return job_id

    def set_create_project_list_by_user(self, user, b_test=False):
        """
        only can run one at a time
        Usually it is called by CollectAllSamples
        """
        process_controler = ProcessControler()
        vect_command = [
            "python3 {} create_projects_list_by_user --user_id {} {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                user.pk,
                "--settings fluwebvirus.settings_test" if b_test else "",
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()

        prefix_to_find = user.profile.get_prefix_name(
            Constants.PROCESS_collect_all_projects, Constants.PROCESS_REGULAR
        )
        job_id = self._get_prefix_in_wait_queue(prefix_to_find)
        if job_id is None:  ## if prefix does not exist in queue need to submit new one
            (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
                Constants.PROCESS_collect_all_projects, Constants.PROCESS_REGULAR
            )
            path_file = self.set_script_run_slurm(
                out_dir,
                Constants.QUEUE_NAME_FAST,
                vect_command,
                job_name,
                True,
                [job_name_wait],
            )
            try:
                job_id = self.submit_job(path_file)
                if job_id != None:
                    self.set_process_controlers(
                        user,
                        process_controler.get_name_collect_all_projects_user(user),
                        job_id,
                    )
            except:
                raise Exception("Fail to submit the job.")
        return job_id

    def set_read_sample_file(self, upload_files, user, b_test=False):
        """
        set read sample file
        param: b_test == True is going to use other settings has also a test database
        """
        process_controler = ProcessControler()
        vect_command = [
            "python3 {} read_sample_file --upload_files_id {} --user_id {} {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                upload_files.pk,
                user.pk,
                (
                    "--settings fluwebvirus.settings_test"
                    if settings.RUN_TEST_IN_COMMAND_LINE or b_test
                    else ""
                ),
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_dont_care, Constants.PROCESS_LINK
        )
        path_file = self.set_script_run_slurm(
            out_dir,
            Constants.QUEUE_NAME_FAST,
            vect_command,
            job_name,
            True,
            [job_name_wait],
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_upload_files(upload_files), job_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def set_read_sample_file_with_metadata(self, upload_files, user):
        """
        update metadata, normal queue, wait for all of other data
        """
        process_controler = ProcessControler()
        vect_command = [
            "python3 {} update_metadata_sample_file --upload_files_id {} --user_id {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"), upload_files.pk, user.pk
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()
        queue_name = user.profile.queue_name_slurm
        if queue_name == None:
            queue_name = Constants.QUEUE_NAME_GLOBAL
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_dont_care, Constants.PROCESS_LINK
        )
        path_file = self.set_script_run_slurm(
            out_dir, queue_name, vect_command, job_name, True, [job_name_wait]
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_upload_files(upload_files), job_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id


    def set_submit_update_televir_project(self, project_id: int, user: User):
        """
        submit job to add references to sample
        """
        process_controler = ProcessControler()

        vect_command = [
            "python3 {} update_project_references --project_id {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                project_id,
            )
        ]

        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        queue_name = user.profile.queue_name_slurm
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_dont_care, Constants.PROCESS_LINK
        )
        outdir_job = self.utils.get_temp_dir()
        path_file = self.set_script_run_slurm(
            outdir_job,
            queue_name,
            vect_command,
            job_name,
            True,
            [job_name_wait],
        )

        try:
            job_id = self.submit_job(path_file)

            if job_id != None:
                pc_name = process_controler.get_name_update_televir_project(project_id)
                self.set_process_controlers(
                    user,
                    pc_name,
                    job_id,
                )
                self.set_process_controlers(
                    user,
                    pc_name,
                    job_id,
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def set_submit_televir_explify_merge(
        self,
        user: User,
        project_pk: int,
        rpip_filepath: str,
        upip_filepath: str,
        televir_report_filepath: str,
        out_dir: str,
    ):
        """
        submit the job to televir
        """
        user_pk = user.pk
        process_controler = ProcessControler()

        vect_command = [
            "python3 {} submit_televir_explify_merge --project_id {} --televir {} --rpip {} --upip {} -o {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                project_pk,
                televir_report_filepath,
                rpip_filepath,
                upip_filepath,
                out_dir,
            )
        ]

        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        queue_name = user.profile.queue_name_slurm
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_televir, Constants.PROCESS_LINK
        )
        outdir_job = self.utils.get_temp_dir()
        path_file = self.set_script_run_slurm(
            outdir_job,
            queue_name,
            vect_command,
            job_name,
            True,
            [job_name_wait],
            alternative_temp_dir=out_dir,
            cpus= Constants.get_process_cpu(Constants.PROCESS_LINK),
            memory= Constants.get_process_mem_string(Constants.PROCESS_LINK),
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                pc_name = process_controler.get_name_televir_project_merge_explify(
                    project_pk,
                )
                self.set_process_controlers(
                    user,
                    pc_name,
                    job_id,
                )

        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def set_submit_televir_explify_merge_external(
        self,
        user: User,
        rpip_filepath: str,
        upip_filepath: str,
        televir_report_filepath: str,
        out_dir: str,
    ):
        """
        submit the job to televir
        """
        user_pk = user.pk
        process_controler = ProcessControler()

        vect_command = [
            "python3 {} submit_televir_explify_merge_ext --user_id {} --televir {} --rpip {} --upip {} -o {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                user_pk,
                televir_report_filepath,
                rpip_filepath,
                upip_filepath,
                out_dir,
            )
        ]

        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        queue_name = user.profile.queue_name_slurm
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_televir, Constants.PROCESS_LINK
        )
        outdir_job = self.utils.get_temp_dir()
        path_file = self.set_script_run_slurm(
            outdir_job,
            queue_name,
            vect_command,
            job_name,
            True,
            [job_name_wait],
            alternative_temp_dir=out_dir,
        )
        try:
            job_id = self.submit_job(path_file)

            if job_id != None:
                pc_name = (
                    process_controler.get_name_televir_project_merge_explify_external(
                        user.pk,
                    )
                )
                self.set_process_controlers(
                    user,
                    pc_name,
                    job_id,
                )

        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def set_submit_televir_sample(self, user, project_pk: int, sample_pk: int, job_name, vect_job_name_wait):
        """
        submit the job to televir
        """
        user_pk = user.pk
        process_controler = ProcessControler()
        outdir_job = self.utils.get_temp_dir()

        vect_command = [
            "python3 {} submit_televir_job_tree_sample --user_id {} --project_id {} --sample_id {} -o {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                user_pk,
                project_pk,
                sample_pk,
                outdir_job,
            )
        ]

        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        queue_name = user.profile.queue_name_slurm

        path_file = self.set_script_run_slurm(
            outdir_job,
            queue_name,
            vect_command,
            job_name,
            True,
            vect_job_name_wait,
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user,
                    process_controler.get_name_televir_project_sample(
                        project_pk, sample_pk
                    ),
                    job_id,
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def set_submit_televir_sample_metagenomics(
        self,
        user,
        sample_pk: int,
        leaf_pk: int,
        combined_analysis: bool = False,
        mapping_request: bool = False,
        map_run_pk: int = None,
    ):
        """
        submit the job to televir
        """
        user_pk = user.pk
        process_controler = ProcessControler()
        outdir_job = self.utils.get_temp_dir()

        vect_command = [
            "python3 {} submit_televir_sample_metagenomics_run --user_id {} --sample_id {} --leaf_id {} {} {} {}-o {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                user_pk,
                sample_pk,
                leaf_pk,
                "--combined_analysis" if combined_analysis else "",
                "--mapping_request" if mapping_request else "",
                "--mapping_run_id {} ".format(map_run_pk) if map_run_pk else "",
                outdir_job,
            )
        ]

        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        queue_name = user.profile.queue_name_slurm
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_televir, Constants.PROCESS_REGULAR
        )
        path_file = self.set_script_run_slurm(
            outdir_job,
            queue_name,
            vect_command,
            job_name,
            True,
            [job_name_wait],
            cpus= Constants.get_process_cpu(Constants.PROCESS_televir),
            memory= Constants.get_process_mem_string(Constants.PROCESS_televir),
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user,
                    process_controler.get_name_televir_project_sample_metagenomics_run(
                        sample_pk,
                        leaf_pk,
                    ),
                    job_id,
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def set_submit_televir_sample_panel_map(
        self,
        user,
        sample_pk: int,
        leaf_pk: int,
        combined_analysis: bool = False,
        mapping_request: bool = False,
        panel_pk: int = None,
    ):
        """
        submit the job to televir
        """
        user_pk = user.pk
        process_controler = ProcessControler()
        outdir_job = self.utils.get_temp_dir()

        vect_command = [
            "python3 {} submit_televir_sample_panel_run --user_id {} --sample_id {} --leaf_id {} {} {} {} -o {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                user_pk,
                sample_pk,
                leaf_pk,
                "--combined_analysis" if combined_analysis else "",
                "--mapping_request" if mapping_request else "",
                "--panel_id {} ".format(panel_pk),
                outdir_job,
            )
        ]

        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        queue_name = user.profile.queue_name_slurm
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_televir, Constants.PROCESS_mapping
        )
        path_file = self.set_script_run_slurm(
            outdir_job,
            queue_name,
            vect_command,
            job_name,
            True,
            [job_name_wait],
            cpus= Constants.get_process_cpu(Constants.PROCESS_televir),
            memory= Constants.get_process_mem_string(Constants.PROCESS_televir),
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user,
                    process_controler.get_name_televir_project_sample_panel_map(
                        sample_pk,
                        leaf_pk,
                    ),
                    job_id,
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def set_submit_upload_reference_televir(
        self,
        user,
        file_id: int,
        fasta: str,
        metadata: str,
    ):
        """
        submit the job to televir
        """
        user_pk = user.pk
        process_controler = ProcessControler()

        vect_command = [
            "python3 {} submit_televir_upload_reference_file --user_id {} --file_id {} --fasta {} --metadata {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                user_pk,
                file_id,
                fasta,
                metadata,
            )
        ]

        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        queue_name = user.profile.queue_name_slurm
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_televir, Constants.PROCESS_LINK
        )
        outdir_job = self.utils.get_temp_dir()
        path_file = self.set_script_run_slurm(
            outdir_job,
            queue_name,
            vect_command,
            job_name,
            True,
            [job_name_wait],
            cpus= Constants.get_process_cpu(Constants.PROCESS_LINK),
            memory= Constants.get_process_mem_string(Constants.PROCESS_LINK),
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user,
                    process_controler.get_name_televir_file_upload(
                        file_id=file_id,
                    ),
                    job_id,
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def set_submit_televir_job(self, user, project_pk):
        """
        submit the job to televir
        """
        user_pk = user.pk
        process_controler = ProcessControler()
        outdir_job = self.utils.get_temp_dir()

        vect_command = [
            "python3 {} submit_televir_job --user_id {} --project_id {} -o {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                user_pk,
                project_pk,
                outdir_job,
            )
        ]

        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        queue_name = user.profile.queue_name_slurm
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_televir, Constants.PROCESS_LINK
        )
        path_file = self.set_script_run_slurm(
            outdir_job,
            queue_name,
            vect_command,
            job_name,
            True,
            [job_name_wait],
            cpus= Constants.get_process_cpu(Constants.PROCESS_televir),
            memory= Constants.get_process_mem_string(Constants.PROCESS_televir),
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_televir_project(project_pk), job_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def set_submit_televir_sort_pisample_reports(self, user, pisample_pk):
        """
        submit the job to televir
        """
        user_pk = user.pk
        process_controler = ProcessControler()
        out_dir = self.utils.get_temp_dir()

        vect_command = [
            "python3 {} submit_televir_job_sort_sample_reports --user_id {} --pisample_id {} -o {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                user_pk,
                pisample_pk,
                out_dir,
            )
        ]

        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        queue_name = user.profile.queue_name_slurm
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_dont_care, Constants.PROCESS_LINK
        )
        outdir_job = self.utils.get_temp_dir()
        path_file = self.set_script_run_slurm(
            outdir_job,
            queue_name,
            vect_command,
            job_name,
            True,
            [job_name_wait],
            alternative_temp_dir=out_dir,
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user,
                    process_controler.get_name_televir_project_sample_sort(pisample_pk),
                    job_id,
                )
        except:
            import traceback
            traceback.print_exc()
            raise Exception("Fail to submit the job.")
        return job_id

    def set_submit_file_televir_teleflu_create(self, user, ref_id):
        """
        submit the job to televir
        """
        user_pk = user.pk
        process_controler = ProcessControler()
        out_dir = self.utils.get_temp_dir()

        vect_command = [
            "python3 {} submit_televir_job_file_teleflu_ref_create --user_id {} --ref_id {} -o {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                user_pk,
                ref_id,
                out_dir,
            )
        ]

        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        queue_name = user.profile.queue_name_slurm
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_dont_care, Constants.PROCESS_LINK
        )
        outdir_job = self.utils.get_temp_dir()
        path_file = self.set_script_run_slurm(
            outdir_job,
            queue_name,
            vect_command,
            job_name,
            True,
            [job_name_wait],
            alternative_temp_dir=out_dir,
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user,
                    process_controler.get_name_raw_televir_teleflu_ref_create(ref_id),
                    job_id,
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def set_submit_televir_teleflu_project_create(self, user, project_pk):
        """
        submit the job to televir
        """
        user_pk = user.pk
        process_controler = ProcessControler()
        out_dir = self.utils.get_temp_dir()

        vect_command = [
            "python3 {} submit_televir_job_teleflu_project_create --user_id {} --project_id {} -o {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                user_pk,
                project_pk,
                out_dir,
            )
        ]

        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        queue_name = user.profile.queue_name_slurm
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_dont_care, Constants.PROCESS_LINK
        )
        outdir_job = self.utils.get_temp_dir()
        path_file = self.set_script_run_slurm(
            outdir_job,
            queue_name,
            vect_command,
            job_name,
            True,
            [job_name_wait],
            alternative_temp_dir=out_dir,
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user,
                    process_controler.get_name_televir_teleflu_project_create(
                        project_pk
                    ),
                    job_id,
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def set_submit_televir_map(self, user, reference_pk, project_pk):
        """
        submit the job to televir
        """
        user_pk = user.pk
        process_controler = ProcessControler()
        outdir_job = self.utils.get_temp_dir()

        vect_command = [
            "python3 {} submit_televir_map_specific --ref_id {} --project_id {} -o {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                reference_pk,
                project_pk,
                outdir_job,
            )
        ]

        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        queue_name = user.profile.queue_name_slurm
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_mapping, Constants.PROCESS_LINK
        )
        path_file = self.set_script_run_slurm(
            outdir_job,
            queue_name,
            vect_command,
            job_name,
            True,
            [job_name_wait],
            cpus= Constants.get_process_cpu(Constants.PROCESS_mapping),
            memory= Constants.get_process_mem_string(Constants.PROCESS_mapping),
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_televir_map(reference_pk), job_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def set_submit_teleflu_map(self, user, leaf_pk, project_pk):
        """
        submit the job to televir
        """
        user_pk = user.pk
        process_controler = ProcessControler()
        outdir_job = self.utils.get_temp_dir()

        vect_command = [
            "python3 {} submit_televir_job_teleflu_stacked_igv --leaf_id {} --project_id {} -o {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                leaf_pk,
                project_pk,
                outdir_job,
            )
        ]

        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        queue_name = user.profile.queue_name_slurm
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_televir, Constants.PROCESS_LINK
        )
        path_file = self.set_script_run_slurm(
            outdir_job,
            queue_name,
            vect_command,
            job_name,
            True,
            [job_name_wait],
            cpus= Constants.get_process_cpu(Constants.PROCESS_LINK),
            memory= Constants.get_process_mem_string(Constants.PROCESS_LINK),
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user,
                    process_controler.get_name_televir_teleflu_igv_stack(leaf_pk),
                    job_id,
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    ### only for tests
    def submit_dummy_job(self, job_name="job_name"):
        """
        only for tests
        """
        vect_command = [
            'echo "#####start" >> /tmp/job.out',
            "echo $HOSTNAME >> /tmp/job.out",
            'echo "start waiting" >> /tmp/job.out',
            "echo $PATH >> /tmp/job.out",
            "date >> /tmp/job.out",
            "sleep 2s",
            "date >> /tmp/job.out",
            'echo "end" >> /tmp/job.out',
        ]
        out_dir = self.utils.get_temp_dir()
        path_file = self.set_script_run_slurm(
            out_dir, Constants.QUEUE_NAME_FAST, vect_command, job_name, True
        )
        os.system("/bin/sh {}".format(path_file))
        try:
            job_id = self.submit_job(path_file)
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    def set_process_controlers(self, user, name_of_process, name_job_id):
        """
        Add a record in ProcessControlers
        """

        process_controler = ProcessControler()
        process_controler.owner = user
        process_controler.name = name_of_process
        process_controler.name_job_id = name_job_id
        process_controler.save()

    def kill_process(self, process_id: str):
        """
        kill the process
        """

        bash_command = (
            # SLURM: scancel process_id
            "scancel {}".format(
                process_id
            )
        )

        exit_status = os.system(bash_command)

        if exit_status != 0:
            print("Fail to kill the process")

        return exit_status

    @transaction.atomic
    def kill_televir_process_controler_runs(
        self, user_pk: int, project_pk: int, sample_pk: int, leaf_pk: int
    ):
        """
        Kill the process in process controler.
        """
        process_controler = ProcessControler()

        names_processes = [
            process_controler.get_name_televir_run(project_pk, sample_pk, leaf_pk),
            process_controler.get_name_televir_project_sample(
                project_pk=project_pk, sample_pk=sample_pk
            ),
            process_controler.get_name_televir_project_sample_metagenomics_run(
                sample_pk,
                leaf_pk,
            ),
            process_controler.get_name_televir_project_sample_panel_map(
                sample_pk=sample_pk, leaf_pk=leaf_pk
            ),
        ]

        processes = ProcessControler.objects.filter(
            owner__id=user_pk,
            name__in=names_processes,
            is_error=False,
            is_finished=False,
        )

        self.kill_processes(processes)

    @transaction.atomic
    def kill_televir_process_controler_samples(
        self, user_pk: int, project_pk: int, sample_pk: int, leaf_pk: int
    ):
        """
        Kill the process in process controler.
        """
        process_controler = ProcessControler()

        names_processes = [
            process_controler.get_name_televir_run(project_pk, sample_pk, leaf_pk),
            process_controler.get_name_televir_project_sample(
                project_pk=project_pk, sample_pk=sample_pk
            ),
            process_controler.get_name_televir_project_sample_metagenomics_run(
                sample_pk,
                leaf_pk,
            ),
            process_controler.get_name_televir_project_sample_panel_map(
                sample_pk=sample_pk, leaf_pk=leaf_pk
            ),
        ]

        processes = ProcessControler.objects.filter(
            owner__id=user_pk,
            name__in=names_processes,
            is_error=False,
            is_finished=False,
        )

        self.kill_processes(processes)

    @transaction.atomic
    def kill_project_samples(self, user_pk: int, project, project_sample_list):
        """
        Kill the processes for the given Project in process controler.
        """
        process_controler = ProcessControler()

        names_processes = []
        names_processes.append(process_controler.get_name_project(project))
        for project_sample in project_sample_list:
            names_processes.append(
                process_controler.get_name_project_sample(project_sample)
            )

        processes = ProcessControler.objects.filter(
            owner__id=user_pk,
            name__in=names_processes,
            is_error=False,
            is_finished=False,
        )

        self.kill_processes(processes)

    @transaction.atomic
    def kill_dataset(self, user_pk: int, dataset):
        """
        Kill the processes for the given Dataset in process controler.
        """
        process_controler = ProcessControler()

        names_processes = []
        names_processes.append(process_controler.get_name_dataset(dataset))

        processes = ProcessControler.objects.filter(
            owner__id=user_pk,
            name__in=names_processes,
            is_error=False,
            is_finished=False,
        )

        self.kill_processes(processes)

    def kill_processes(self, processes: List[ProcessControler]):
        """ """

        for process in processes:
            if process.name_job_id:
                self.kill_process(process.name_job_id)

            process.is_running = False
            process.is_finished = False
            process.is_error = True
            process.save()

    def set_specific_controler_flag(self, user, name_of_process, job_id, flags):
        try:
            process_controler = ProcessControler.objects.get(
                owner__id=user.pk,
                name=name_of_process,
                name_job_id=job_id,
            )

        except ProcessControler.DoesNotExist:
            process_controler = ProcessControler(
                owner=user,
                name=name_of_process,
                name_job_id=job_id,
            )

        if flags == ProcessControler.FLAG_FINISHED:
            process_controler.is_finished = True
            process_controler.is_running = False
            process_controler.is_error = False
        elif flags == ProcessControler.FLAG_RUNNING:
            process_controler.is_finished = False
            process_controler.is_running = True
            process_controler.is_error = False
        elif flags == ProcessControler.FLAG_ERROR:
            process_controler.is_finished = False
            process_controler.is_running = False
            process_controler.is_error = True

        process_controler.save()

    def set_process_controler(self, user, name_of_process, flags):
        """
        name_of_process:
                process_controler.get_name_upload_files(upload_files),
                process_controler.get_name_link_files_user(user),
                process_controler.get_name_sample(sample),
                process_controler.get_name_project(project), job_id)
                process_controler.get_name_project_sample(project_sample)

        flags: ProcessControler.FLAG_FINISHED, ProcessControler.FLAG_RUNNING, ProcessControler.FLAG_ERROR
        """

        if flags == ProcessControler.FLAG_FINISHED:
            data_set = ProcessControler.objects.filter(
                owner__id=user.pk,
                name=name_of_process,
                is_running=True,
                is_finished=False,
                is_error=False,
            )
        elif flags == ProcessControler.FLAG_ERROR:
            data_set = ProcessControler.objects.filter(
                owner__id=user.pk,
                name=name_of_process,
                is_finished=False,
                is_error=False,
            )

        else:
            data_set = ProcessControler.objects.filter(
                owner__id=user.pk,
                name=name_of_process,
                is_running=False,
                is_finished=False,
                is_error=False,
            )

        if data_set.count() > 0:
            for process_controler in data_set:
                # process_controler = ProcessControler.objects.get(pk=data_set[0].pk)
                if flags == ProcessControler.FLAG_FINISHED:
                    process_controler.is_finished = True
                    process_controler.is_running = False
                    process_controler.close_date = datetime.now()
                elif flags == ProcessControler.FLAG_ERROR:
                    process_controler.is_finished = True
                    process_controler.is_error = True
                    process_controler.is_running = False
                    process_controler.close_date = datetime.now()
                elif flags == ProcessControler.FLAG_RUNNING:
                    process_controler.is_running = True
                process_controler.save()

    ##### set collect global files
    def set_collect_dataset_global_files(self, dataset, user):
        """
        job_name = "job_name_<user_id>_<seq_id>"
        """

        # TODO: If dataset is running, do not run again... fail with error...

        process_controler = ProcessControler()
        vect_command = [
            "python3 {} collect_global_dataset_files --dataset_id {} --user_id {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"), dataset.pk, user.pk
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()

        queue_name = user.profile.queue_name_slurm
        if queue_name == None:
            queue_name = Constants.QUEUE_NAME_GLOBAL

        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_datasets, Constants.PROCESS_GLOBAL
        )
        path_file = self.set_script_run_slurm(
            out_dir, queue_name, vect_command, job_name, True, [job_name_wait]
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_dataset(dataset), job_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id

    ##### set collect global files
    def set_collect_dataset_global_files_for_update_metadata(self, dataset, user):
        """
        job_name = "job_name_<user_id>_<seq_id>"
        only run this task after all second_stage_snippy
        """

        dataset.is_processed = False
        dataset.save()

        process_controler = ProcessControler()
        vect_command = [
            "python3 {} collect_global_dataset_files_for_update_metadata --dataset_id {} --user_id {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"), dataset.pk, user.pk
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()

        queue_name = user.profile.queue_name_slurm
        if queue_name == None:
            queue_name = Constants.QUEUE_NAME_GLOBAL
        (job_name_wait, job_name) = user.profile.get_name_slurm_seq(
            Constants.PROCESS_datasets, Constants.PROCESS_GLOBAL
        )
        path_file = self.set_script_run_slurm(
            out_dir, queue_name, vect_command, job_name, True, [job_name_wait]
        )
        try:
            job_id = self.submit_job(path_file)
            if job_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_dataset(dataset), job_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return job_id
