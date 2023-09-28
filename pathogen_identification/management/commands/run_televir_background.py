import os
import time
from datetime import date

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

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


def get_samples_available(project_id: int):
    samples = PIProject_Sample.objects.filter(project__pk=project_id)
    samples_available = []
    for sample in samples:
        if check_sample_available(sample):
            samples_available.append(sample)

    return samples_available


def check_sample_deployed(sample: PIProject_Sample):
    parameter_set = ParameterSet.objects.filter(sample=sample)
    if parameter_set.exists():
        for ps in parameter_set:
            if ps.status in [
                ParameterSet.STATUS_RUNNING,
                ParameterSet.STATUS_QUEUED,
            ]:
                return True

    return False


def count_samples_available(project_id: int):
    samples = PIProject_Sample.objects.filter(project__pk=project_id)
    count = 0
    for sample in samples:
        if check_sample_available(sample):
            count += 1

    return count


def count_samples_deployed(project_id: int):
    samples = PIProject_Sample.objects.filter(project__pk=project_id)
    count = 0
    for sample in samples:
        if check_sample_deployed(sample):
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


class DeploymentManager:
    python_bin = "/usr/bin/python3"

    def __init__(
        self, project_id: int, output_dir: str, log_dir: str, max_threads: int
    ):
        self.project_id = project_id
        self.output_dir = output_dir
        self.log_dir = log_dir
        self.max_threads = max_threads

    @staticmethod
    def nohup_wrapper(command: str, output_dir: str):
        nohup = f"nohup {command} > {output_dir}/nohup.out 2> {output_dir}/nohup.err &"

        return nohup

    def deploy_sample_in_background(self, sample: PIProject_Sample):
        project_id = self.project_id
        command_base = TelevirTreeSample().command_args_wrapper(
            project_id, sample, self.output_dir
        )

        command = f"{self.python_bin} manage.py {command_base} "

        nohup = self.nohup_wrapper(command, self.log_dir)
        print(nohup)

        sys_out = os.system(nohup)

        return sys_out

    def find_sample_to_deploy(self):
        samples_deployed = count_samples_deployed(self.project_id)
        if samples_deployed >= self.max_threads:
            return None
        samples_available = get_samples_available(self.project_id)

        for sample in samples_available:
            if check_sample_available(sample):
                return sample

        return None


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

    def handle(self, *args, **options):
        stop = False
        wait_time = 60

        # break if log_dir or out_dir does not exist
        if not os.path.exists(options["log_dir"]):
            print(f"log_dir {options['log_dir']} does not exist")
            stop = True

        if not os.path.exists(options["out_dir"]):
            print(f"out_dir {options['out_dir']} does not exist")
            stop = True

        while not stop:
            ##
            manager = DeploymentManager(
                options["project_id"],
                options["out_dir"],
                options["log_dir"],
                options["max_threads"],
            )
            sample_to_deploy = manager.find_sample_to_deploy()

            if sample_to_deploy is not None:
                sys_out = manager.deploy_sample_in_background(sample_to_deploy)
                if sys_out == 0:
                    print(f"Sample {sample_to_deploy.pk} deployed")
                else:
                    print(f"Sample {sample_to_deploy.pk} not deployed")
                    break
            else:
                print("No sample to deploy")
                stop = True

            ## wait for 60 seconds

            time.sleep(wait_time)
