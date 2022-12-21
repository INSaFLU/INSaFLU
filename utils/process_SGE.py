#!/usr/bin/env python

import logging
import os
import time
from datetime import datetime

from constants.constants import Constants, FileExtensions, TypePath
from django.conf import settings
from extend_user.models import Profile
from managing_files.models import ProcessControler

from utils.utils import Utils

# http://www.socher.org/index.php/Main/HowToInstallSunGridEngineOnUbuntu
# https://peteris.rocks/blog/sun-grid-engine-installation-on-ubuntu-server/
# http://biohpc.blogspot.pt/2016/10/sge-installation-of-son-of-grid-engine.html ## centos 7
# http://star.mit.edu/cluster/docs/0.93.3/guides/sge.html        ### explain who to use SGE

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
    ###        BEGIN main methods
    ###
    ###        IMPORTANT
    ###            Put qsub in /usr/bin/qsub
    ###

    def submitte_job(self, file_name):
        """
        job submission
        raise exception if something wrong
        """
        temp_file = self.utils.get_temp_file("qsub_out", FileExtensions.FILE_TXT)

        cmd = "export SGE_ROOT={}; qsub {} > {}".format(
            settings.SGE_ROOT, file_name, temp_file
        )
        exist_status = os.system(cmd)
        print("cmd: {}".format(cmd))

        if exist_status != 0:
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
        vect_out = self.utils.read_text_file(temp_file)
        if os.path.exists(temp_file):
            os.unlink(temp_file)
        b_found = False

        for line in vect_out:
            if line.find("has been submitted") != -1:
                lst_line = line.split(" ")
                if len(lst_line) > 4 and self.utils.is_integer(lst_line[2]):
                    return int(lst_line[2])
                return None  ## don't rise exception...
        if not b_found:
            raise Exception("\n".join(vect_out))

    def set_script_run_sge(
        self,
        out_dir,
        queue_name,
        vect_cmd,
        job_name,
        b_remove_out_dir=False,
        job_name_wait=[],
        alternative_temp_dir=None,
    ):
        """
        create the script to run SGE
        """
        if len(vect_cmd) == 0:
            return None

        file_name_out = os.path.join(out_dir, ProcessSGE.FILE_NAME_SCRIPT_SGE)
        with open(file_name_out, "w") as handleSGE:
            handleSGE.write("#!/bin/bash\n")
            handleSGE.write(
                "#$ -V\n"
            )  # Specifies  that  all  environment  variables active
            # within the qsub utility be exported to the context of the job.
            handleSGE.write("#$ -S /bin/bash\n")  # interpreting shell
            ## hold_jid <comma separated list of job-ids, can also be a job id pattern such as 2722*> :
            ## will start the current job/job -array only after completion of all jobs in the comma separated list
            if isinstance(job_name_wait, str):
                job_name_wait = [job_name_wait]
            if len(job_name_wait) > 0:
                handleSGE.write(
                    "#$ -hold_jid {}\n".format(",".join(job_name_wait))
                )  # need to wait until all this jobs names finished
            handleSGE.write(
                "#$ -j y\n"
            )  # merge the standard error with standard output
            handleSGE.write("#$ -N {}\n".format(job_name))  # job name
            handleSGE.write(
                "#$ -cwd\n"
            )  # execute the job for the current work directory
            handleSGE.write("#$ -q {}\n".format(queue_name))  # queue name
            handleSGE.write("#$ -o {}\n".format(out_dir))  # out path file
            for cline in vect_cmd:
                handleSGE.write("\n" + cline)
            if b_remove_out_dir and not settings.RUN_TEST_IN_COMMAND_LINE:
                handleSGE.write(
                    "\nif [ $? -eq 0 ]\nthen\n  rm -r {}\nfi\n".format(out_dir)
                )
                if alternative_temp_dir is not None:
                    handleSGE.write(
                        "\nif [ $? -eq 0 ]\nthen\n  rm -r {}\nfi\n".format(
                            alternative_temp_dir
                        )
                    )

        return file_name_out

    def __get_sge_process__(self):
        """
        #Job status - one of

        ### test if all jobs submitted to the SGE are finish
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
        tagsSGERunning = ("r", "t")
        tagsSGEWaiting = ("hqw", "qw", "w")
        # test with qstat
        file_result = self.utils.get_temp_file("sge_stat", ".txt")
        cline = "export SGE_ROOT={}; qstat > {}".format(settings.SGE_ROOT, file_result)
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
                    if line.split()[4] in tagsSGERunning:
                        vectRunning.append(line.split()[0])
                    elif line.split()[4] in tagsSGEWaiting:
                        vectWait.append(line.split()[0])

        ## remove file
        if os.path.exists(file_result):
            os.unlink(file_result)
        return (vectRunning, vectWait)

    def get_status_process(self, n_SGE_id):
        (vectRunning, vectWait) = self.__get_sge_process__()
        if str(n_SGE_id) in vectRunning:
            return self.SGE_JOB_ID_PROCESSING
        if str(n_SGE_id) in vectWait:
            return self.SGE_JOB_ID_QUEUE
        return self.SGE_JOB_ID_FINISH

    def is_finished(self, n_SGE_id):
        """
        is it finished
        """
        return self.get_status_process(n_SGE_id) == self.SGE_JOB_ID_FINISH

    def exists_taks_running(self):
        """
        test if there any tasks running...
        """
        file_result = self.utils.get_temp_file("sge_stat", ".txt")
        cline = "export SGE_ROOT={}; qstat > {}".format(settings.SGE_ROOT, file_result)
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
        file_result = self.utils.get_temp_file("sge_stat", ".txt")
        cline = "export SGE_ROOT={}; qstat -j {}* > {}".format(
            settings.SGE_ROOT, prefix_id, file_result
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

    def wait_until_finished(self, vect_sge_ids):
        """
        wait till all end
        if len(vect_sge_ids) == 0 wait till all are finished, doesn't matter the ID
        """
        ## expand vect_sge_idsbecaus some lines can have morethan one ID
        vect_sge_to_search = [c for b in vect_sge_ids for c in str(b).split(",")]
        if len(vect_sge_to_search) == 0:
            while self.exists_taks_running():
                print("=" * 50 + "\n  waiting for sge\n" + str(datetime.now()))
                time.sleep(5)  ## wais 5 seconds
        else:
            while len(vect_sge_to_search) > 0:
                print("=" * 50)
                print(
                    "   wait for these ids: {}".format(
                        ";".join([str(_) for _ in vect_sge_to_search])
                    )
                )
                vect_remove = []
                for sge_id in vect_sge_to_search:
                    if self.is_finished(sge_id):
                        vect_remove.append(sge_id)

                ### remove sge
                for sge_id in vect_remove:
                    vect_sge_to_search.remove(sge_id)

                print("=" * 50)
                if len(vect_sge_to_search) > 0:
                    time.sleep(5)  ## wais 5 seconds
        ## set the one still running
        return [int(_) for _ in vect_sge_to_search]

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

        queue_name = user.profile.queue_name_sge
        if queue_name == None:
            queue_name = Constants.QUEUE_SGE_NAME_GLOBAL

        (job_name_wait, job_name) = user.profile.get_name_sge_seq(
            Profile.SGE_PROCESS_projects, Profile.SGE_GLOBAL
        )
        path_file = self.set_script_run_sge(
            out_dir, queue_name, vect_command, job_name, True, [job_name_wait]
        )
        try:
            sge_id = self.submitte_job(path_file)
            if sge_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_project(project), sge_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return sge_id

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

        queue_name = user.profile.queue_name_sge
        if queue_name == None:
            queue_name = Constants.QUEUE_SGE_NAME_GLOBAL
        (job_name_wait, job_name) = user.profile.get_name_sge_seq(
            Profile.SGE_PROCESS_dont_care, Profile.SGE_GLOBAL
        )
        path_file = self.set_script_run_sge(
            out_dir, queue_name, vect_command, job_name, True, [job_name_wait]
        )
        try:
            sge_id = self.submitte_job(path_file)
            if sge_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_project(project), sge_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return sge_id

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

        queue_name = user.profile.queue_name_sge
        if queue_name == None:
            queue_name = Constants.QUEUE_SGE_NAME_GLOBAL
        (job_name_wait, job_name) = user.profile.get_name_sge_seq(
            Profile.SGE_PROCESS_dont_care, Profile.SGE_GLOBAL
        )
        path_file = self.set_script_run_sge(
            out_dir, queue_name, vect_command, job_name, True, [job_name_wait]
        )
        try:
            sge_id = self.submitte_job(path_file)
            if sge_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_project(project), sge_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return sge_id

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
                "--settings fluwebvirus.settings_test"
                if settings.RUN_TEST_IN_COMMAND_LINE
                else "",
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()
        queue_name = user.profile.queue_name_sge
        if queue_name == None:
            queue_name = Constants.QUEUE_SGE_NAME_GLOBAL
        path_file = self.set_script_run_sge(
            out_dir, queue_name, vect_command, job_name, True, vect_job_name_wait
        )
        try:
            sge_id = self.submitte_job(path_file)
            if sge_id != None:
                self.set_process_controlers(
                    user,
                    process_controler.get_name_project_sample(project_sample),
                    sge_id,
                )
        except:
            raise Exception("Fail to submit the job.")
        return sge_id

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
                "--settings fluwebvirus.settings_test"
                if settings.RUN_TEST_IN_COMMAND_LINE
                else "",
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()
        queue_name = user.profile.queue_name_sge
        if queue_name == None:
            queue_name = Constants.QUEUE_SGE_NAME_GLOBAL
        path_file = self.set_script_run_sge(
            out_dir, queue_name, vect_command, job_name, True, vect_job_name_wait
        )
        try:
            sge_id = self.submitte_job(path_file)
            if sge_id != None:
                self.set_process_controlers(
                    user,
                    process_controler.get_name_project_sample(project_sample),
                    sge_id,
                )
        except:
            raise Exception("Fail to submit the job.")
        return sge_id

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
                "--settings fluwebvirus.settings_test"
                if settings.RUN_TEST_IN_COMMAND_LINE
                else "",
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()
        path_file = self.set_script_run_sge(
            out_dir, Constants.QUEUE_SGE_NAME_GLOBAL, vect_command, job_name, True
        )
        try:
            self._remove_files_create_by_identify_type_and_sub_type(sample)
            self._remove_files_create_by_fastq_and_trimmomatic(sample)
            sge_id = self.submitte_job(path_file)
            if sge_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_sample(sample), sge_id
                )

            ### change flag to not finished
            sample.is_ready_for_projects = False
            sample.is_sample_in_the_queue = True
            sample.save()
        except:
            raise Exception("Fail to submit the job.")
        return sge_id

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
                "--settings fluwebvirus.settings_test"
                if settings.RUN_TEST_IN_COMMAND_LINE
                else "",
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()
        path_file = self.set_script_run_sge(
            out_dir, Constants.QUEUE_SGE_NAME_GLOBAL, vect_command, job_name, True
        )
        try:
            self._remove_files_create_by_identify_type_and_sub_type(sample)
            self._remove_files_create_by_nanofilt_and_stat(sample)
            sge_id = self.submitte_job(path_file)
            if sge_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_sample(sample), sge_id
                )

            ### change flag to not finished
            sample.is_ready_for_projects = False
            sample.is_sample_in_the_queue = True
            sample.save()
        except:
            raise Exception("Fail to submit the job.")
        return sge_id

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
            Profile.SGE_PROCESS_link_files, Profile.SGE_LINK
        )
        sge_id = self._get_prefix_in_wait_queue(prefix_to_find)
        if sge_id is None:  ## if prefix does not exist in queue need to submit new one
            (job_name_wait, job_name) = user.profile.get_name_sge_seq(
                Profile.SGE_PROCESS_link_files, Profile.SGE_LINK
            )
            path_file = self.set_script_run_sge(
                out_dir,
                Constants.QUEUE_SGE_NAME_FAST,
                vect_command,
                job_name,
                True,
                [job_name_wait],
            )
            try:
                sge_id = self.submitte_job(path_file)
                if sge_id != None:
                    self.set_process_controlers(
                        user, process_controler.get_name_link_files_user(user), sge_id
                    )
            except:
                raise Exception("Fail to submit the job.")
        return sge_id

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
            Profile.SGE_PROCESS_collect_all_samples, Profile.SGE_REGULAR
        )
        sge_id = self._get_prefix_in_wait_queue(prefix_to_find)
        if sge_id is None:  ## if prefix does not exist in queue need to submit new one
            (job_name_wait, job_name) = user.profile.get_name_sge_seq(
                Profile.SGE_PROCESS_collect_all_samples, Profile.SGE_REGULAR
            )
            path_file = self.set_script_run_sge(
                out_dir,
                Constants.QUEUE_SGE_NAME_FAST,
                vect_command,
                job_name,
                True,
                vect_job_wait,
            )
            try:
                sge_id = self.submitte_job(path_file)
                if sge_id != None:
                    self.set_process_controlers(
                        user,
                        process_controler.get_name_collect_all_samples_user(user),
                        sge_id,
                    )
            except:
                raise Exception("Fail to submit the job.")
        return sge_id

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
            Profile.SGE_PROCESS_collect_all_projects, Profile.SGE_REGULAR
        )
        sge_id = self._get_prefix_in_wait_queue(prefix_to_find)
        if sge_id is None:  ## if prefix does not exist in queue need to submit new one
            (job_name_wait, job_name) = user.profile.get_name_sge_seq(
                Profile.SGE_PROCESS_collect_all_projects, Profile.SGE_REGULAR
            )
            path_file = self.set_script_run_sge(
                out_dir,
                Constants.QUEUE_SGE_NAME_FAST,
                vect_command,
                job_name,
                True,
                [job_name_wait],
            )
            try:
                sge_id = self.submitte_job(path_file)
                if sge_id != None:
                    self.set_process_controlers(
                        user,
                        process_controler.get_name_collect_all_projects_user(user),
                        sge_id,
                    )
            except:
                raise Exception("Fail to submit the job.")
        return sge_id

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
                "--settings fluwebvirus.settings_test"
                if settings.RUN_TEST_IN_COMMAND_LINE or b_test
                else "",
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()
        (job_name_wait, job_name) = user.profile.get_name_sge_seq(
            Profile.SGE_PROCESS_dont_care, Profile.SGE_LINK
        )
        path_file = self.set_script_run_sge(
            out_dir,
            Constants.QUEUE_SGE_NAME_FAST,
            vect_command,
            job_name,
            True,
            [job_name_wait],
        )
        try:
            sge_id = self.submitte_job(path_file)
            if sge_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_upload_files(upload_files), sge_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return sge_id

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
        queue_name = user.profile.queue_name_sge
        if queue_name == None:
            queue_name = Constants.QUEUE_SGE_NAME_GLOBAL
        (job_name_wait, job_name) = user.profile.get_name_sge_seq(
            Profile.SGE_PROCESS_dont_care, Profile.SGE_LINK
        )
        path_file = self.set_script_run_sge(
            out_dir, queue_name, vect_command, job_name, True, [job_name_wait]
        )
        try:
            sge_id = self.submitte_job(path_file)
            if sge_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_upload_files(upload_files), sge_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return sge_id

    def set_submit_televir_job(self, user, project_pk):
        """
        submit the job to televir
        """
        user_pk = user.pk
        process_controler = ProcessControler()
        out_dir = self.utils.get_temp_dir()

        vect_command = [
            "python3 {} submit_televir_job --user_id {} --project_id {} -o {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                user_pk,
                project_pk,
                out_dir,
            )
        ]

        print(vect_command)

        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        queue_name = user.profile.queue_name_sge
        (job_name_wait, job_name) = user.profile.get_name_sge_seq(
            Profile.SGE_PROCESS_dont_care, Profile.SGE_LINK
        )
        outdir_sge = self.utils.get_temp_dir()
        path_file = self.set_script_run_sge(
            outdir_sge,
            queue_name,
            vect_command,
            job_name,
            True,
            [job_name_wait],
            alternative_temp_dir=out_dir,
        )
        try:
            sge_id = self.submitte_job(path_file)
            print("project submitted, sge_id: " + str(sge_id))
            if sge_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_televir_project(project_pk), sge_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return sge_id

    def set_submit_televir_map(self, user, reference_pk):
        """
        submit the job to televir
        """
        user_pk = user.pk
        process_controler = ProcessControler()
        out_dir = self.utils.get_temp_dir()

        vect_command = [
            "python3 {} submit_televir_map --ref_id {} -o {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"),
                reference_pk,
                out_dir,
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        queue_name = user.profile.queue_name_sge
        (job_name_wait, job_name) = user.profile.get_name_sge_seq(
            Profile.SGE_PROCESS_dont_care, Profile.SGE_LINK
        )
        outdir_sge = self.utils.get_temp_dir()
        path_file = self.set_script_run_sge(
            outdir_sge,
            queue_name,
            vect_command,
            job_name,
            True,
            [job_name_wait],
            alternative_temp_dir=out_dir,
        )
        try:
            sge_id = self.submitte_job(path_file)
            print("project submitted, sge_id: " + str(sge_id))
            if sge_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_televir_map(reference_pk), sge_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return sge_id

    ### only for tests
    def submit_dummy_sge(self, job_name="job_name"):
        """
        only for tests
        """
        vect_command = [
            'echo "#####start" >> /tmp/sge.out',
            "echo $HOSTNAME >> /tmp/sge.out",
            'echo "start waiting" >> /tmp/sge.out',
            "echo $PATH >> /tmp/sge.out",
            "echo $SGE_ROOT >> /tmp/sge.out",
            "date >> /tmp/sge.out",
            "sleep 2s",
            "date >> /tmp/sge.out",
            'echo "end" >> /tmp/sge.out',
        ]
        out_dir = self.utils.get_temp_dir()
        path_file = self.set_script_run_sge(
            out_dir, Constants.QUEUE_SGE_NAME_FAST, vect_command, job_name, False
        )
        os.system("/bin/sh {}".format(path_file))
        try:
            sge_id = self.submitte_job(path_file)
        except:
            raise Exception("Fail to submit the job.")
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
            process_controler = ProcessControler.objects.get(pk=data_set[0].pk)
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

        queue_name = user.profile.queue_name_sge
        if queue_name == None:
            queue_name = Constants.QUEUE_SGE_NAME_GLOBAL

        (job_name_wait, job_name) = user.profile.get_name_sge_seq(
            Profile.SGE_PROCESS_datasets, Profile.SGE_GLOBAL
        )
        path_file = self.set_script_run_sge(
            out_dir, queue_name, vect_command, job_name, True, [job_name_wait]
        )
        try:
            sge_id = self.submitte_job(path_file)
            if sge_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_dataset(dataset), sge_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return sge_id

    ##### set collect global files
    def set_collect_dataset_global_files_for_update_metadata(self, dataset, user):
        """
        job_name = "job_name_<user_id>_<seq_id>"
        only run this task after all second_stage_snippy
        """
        process_controler = ProcessControler()
        vect_command = [
            "python3 {} collect_global_dataset_files_for_update_metadata --dataset_id {} --user_id {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"), dataset.pk, user.pk
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()

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

        queue_name = user.profile.queue_name_sge
        if queue_name == None:
            queue_name = Constants.QUEUE_SGE_NAME_GLOBAL

        (job_name_wait, job_name) = user.profile.get_name_sge_seq(
            Profile.SGE_PROCESS_datasets, Profile.SGE_GLOBAL
        )
        path_file = self.set_script_run_sge(
            out_dir, queue_name, vect_command, job_name, True, [job_name_wait]
        )
        try:
            sge_id = self.submitte_job(path_file)
            if sge_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_dataset(dataset), sge_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return sge_id

    ##### set collect global files
    def set_collect_dataset_global_files_for_update_metadata(self, dataset, user):
        """
        job_name = "job_name_<user_id>_<seq_id>"
        only run this task after all second_stage_snippy
        """
        process_controler = ProcessControler()
        vect_command = [
            "python3 {} collect_global_dataset_files_for_update_metadata --dataset_id {} --user_id {}".format(
                os.path.join(settings.BASE_DIR, "manage.py"), dataset.pk, user.pk
            )
        ]
        self.logger_production.info("Processing: " + ";".join(vect_command))
        self.logger_debug.info("Processing: " + ";".join(vect_command))
        out_dir = self.utils.get_temp_dir()

        queue_name = user.profile.queue_name_sge
        if queue_name == None:
            queue_name = Constants.QUEUE_SGE_NAME_GLOBAL
        (job_name_wait, job_name) = user.profile.get_name_sge_seq(
            Profile.SGE_PROCESS_datasets, Profile.SGE_GLOBAL
        )
        path_file = self.set_script_run_sge(
            out_dir, queue_name, vect_command, job_name, True, [job_name_wait]
        )
        try:
            sge_id = self.submitte_job(path_file)
            if sge_id != None:
                self.set_process_controlers(
                    user, process_controler.get_name_dataset(dataset), sge_id
                )
        except:
            raise Exception("Fail to submit the job.")
        return sge_id
