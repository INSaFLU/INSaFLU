import os
import time
from datetime import date

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from django.db.models import Q

from pathogen_identification.models import ParameterSet, PIProject_Sample


def check_sample_available(sample: PIProject_Sample):
    parameter_set = ParameterSet.objects.filter(sample=sample)
    if parameter_set.exists():
        for ps in parameter_set:
            if ps.status in [
                ParameterSet.STATUS_FINISHED,
                ParameterSet.STATUS_RUNNING,
                ParameterSet.STATUS_QUEUED,
            ]:
                return False
    return True


def check_sample_future(sample: PIProject_Sample):
    parameter_set = ParameterSet.objects.filter(sample=sample).exclude(
        status__in=[ParameterSet.STATUS_FINISHED]
    )
    if parameter_set.exists():
        False

    return True


def count_samples_available(project_id: int):
    samples = PIProject_Sample.objects.filter(project__pk=project_id)
    count = 0
    for sample in samples:
        if check_sample_available(sample):
            count += 1

    return count


from abc import ABC, abstractmethod


class InsafluCommand(ABC):
    command = "submit_televir_job_tree_sample"

    def command_wrapper(self, args: dict):
        command = f"{self.command} {' '.join([f'{k} {v}' for k, v in args.items()])}"
        return command

    @abstractmethod
    def command_args(
        self, project_id: int, sample: PIProject_Sample, output_dir: str
    ) -> dict:
        pass

    def command_args_wrapper(
        self, project_id: int, sample: PIProject_Sample, output_dir: str
    ):
        args = self.command_args(project_id, sample, output_dir)
        return self.command_wrapper(args)


class TelevirTreeSample(InsafluCommand):
    command = "submit_televir_job_tree_sample"

    def command_args(
        self, project_id: int, sample: PIProject_Sample, output_dir: str
    ) -> dict:
        args = {
            "--project_id": project_id,
            "--sample_id": sample.pk,
            "-o": output_dir,
            "--user_id": sample.project.owner.pk,
        }

        return args


import logging


class DeploymentManager:
    python_bin = "/usr/bin/python3"

    def __init__(
        self, project_id: int, output_dir: str, log_dir: str, max_threads: int
    ):
        self.project_id = project_id
        self.output_dir = output_dir
        self.log_dir = log_dir
        self.max_threads = max_threads
        self.pid_deployed = []

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

        # create a file handler
        handler = logging.FileHandler(f"{self.log_dir}/deployment_manager.log")
        handler.setLevel(logging.INFO)
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)

    def update_pid_deployed(self):
        """Update pid_deployed list with samples that are running or queued, remove samples that are finished"""
        new_pid_list = []
        for pid in self.pid_deployed:
            if ParameterSet.objects.filter(
                sample__pk=pid,
                status__in=[
                    ParameterSet.STATUS_RUNNING,
                    ParameterSet.STATUS_QUEUED,
                ],
            ).exists():
                new_pid_list.append(pid)

        self.pid_deployed = new_pid_list

    def check_sample_deployed(self, sample: PIProject_Sample):
        if sample.pk in self.pid_deployed:
            return True

        parameter_set = ParameterSet.objects.filter(sample=sample)
        if parameter_set.exists():
            for ps in parameter_set:
                if ps.status in [
                    ParameterSet.STATUS_RUNNING,
                    ParameterSet.STATUS_QUEUED,
                ]:
                    return True

        return False

    def count_samples_deployed(self, project_id):
        samples = PIProject_Sample.objects.filter(project__pk=project_id)
        count = 0
        for sample in samples:
            if self.check_sample_deployed(sample):
                count += 1

        return count

    @staticmethod
    def nohup_wrapper(command: str, output_dir: str, id_job: int):
        nohup = f"nohup {command} > {output_dir}/nohup.{id_job}.out 2> {output_dir}/nohup.{id_job}.err &"

        return nohup

    def get_samples_available(self, project_id: int):
        samples = PIProject_Sample.objects.filter(project__pk=project_id).exclude(
            pk__in=self.pid_deployed
        )
        samples_available = []
        for sample in samples:
            if check_sample_available(sample):
                samples_available.append(sample)

        return samples_available

    def count_samples_future(self, project_id: int):
        samples = PIProject_Sample.objects.filter(project__pk=project_id).exclude(
            pk__in=self.pid_deployed
        )
        count = 0
        for sample in samples:
            if check_sample_available(sample) is True:
                count += 1

        return count

    def deploy_sample_in_background(self, sample: PIProject_Sample):
        project_id = self.project_id
        command_base = TelevirTreeSample().command_args_wrapper(
            project_id, sample, self.output_dir
        )

        command = f"{self.python_bin} manage.py {command_base} "

        nohup = self.nohup_wrapper(command, self.log_dir, sample.pk)

        sys_out = os.system(nohup)
        self.pid_deployed.append(sample.pk)

        return sys_out

    def find_sample_to_deploy(self):
        samples_deployed = self.count_samples_deployed(self.project_id)
        if samples_deployed >= self.max_threads:
            return None

        samples_available = self.get_samples_available(self.project_id)
        self.logger.info(f"Samples available: {len(samples_available)}")

        for sample in samples_available:
            if check_sample_available(sample):
                return sample

        return None

    def samples_remain(self):
        samples_remaining_n = self.count_samples_future(self.project_id)
        self.logger.info(f"Samples remaining: {samples_remaining_n}")
        if samples_remaining_n > 0:
            return True

        return False


class Command(BaseCommand):
    help = "nohup televir command to run in background"

    def add_arguments(self, parser):
        parser.add_argument(
            "--project_id",
            type=int,
            help="televir project to be merged (pk)",
        )

        parser.add_argument(
            "--out_dir",
            type=str,
            help="tmp directory",
        )

        parser.add_argument(
            "--log_dir",
            type=str,
            help="log directory",
        )

        parser.add_argument(
            "--max_threads",
            type=int,
            help="maximum number of threads",
        )

        parser.add_argument(
            "--wait_time",
            type=int,
            help="wait time between checks, default 5 minutes",
            default=5 * 60,
        )

    def handle(self, *args, **options):
        stop = False

        wait_time = options["wait_time"]

        # break if log_dir or out_dir does not exist
        if not os.path.exists(options["log_dir"]):
            print(f"log_dir {options['log_dir']} does not exist")
            stop = True

        if not os.path.exists(options["out_dir"]):
            print(f"out_dir {options['out_dir']} does not exist")
            stop = True

        manager = DeploymentManager(
            options["project_id"],
            options["out_dir"],
            options["log_dir"],
            options["max_threads"],
        )

        time_elased = 0

        while not stop:
            ##
            manager.update_pid_deployed()
            sample_to_deploy = manager.find_sample_to_deploy()

            while sample_to_deploy is not None:
                sys_out = manager.deploy_sample_in_background(sample_to_deploy)

                if sys_out == 0:
                    print(f"Sample {sample_to_deploy.pk} deployed")
                else:
                    print(f"Error, Sample {sample_to_deploy.pk} not deployed")
                    stop = True

                sample_to_deploy = manager.find_sample_to_deploy()

            else:
                if manager.samples_remain() is False:
                    print("No sample to deploy")
                    stop = True

            ## wait for 60 seconds

            time.sleep(wait_time)
            time_elased += wait_time
