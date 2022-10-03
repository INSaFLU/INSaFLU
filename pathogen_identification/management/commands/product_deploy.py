import os
from datetime import date

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from pathogen_identification.deployment_main import Run_Main_from_Leaf
from pathogen_identification.models import (
    PIProject_Sample,
    Projects,
    SoftwareTree,
    SoftwareTreeNode,
)


class Command(BaseCommand):
    help = "deploy run"

    def add_arguments(self, parser):
        parser.add_argument(
            "--user_pk",
            type=int,
            help="user deploying the run",
        )
        parser.add_argument(
            "--sample_pk",
            type=int,
            help="sample to be run",
        )

        parser.add_argument(
            "--project_pk",
            type=int,
            help="project to be run",
        )

        parser.add_argument(
            "--pipeline_pk",
            type=int,
            help="pipeline to be run",
        )

    def handle(self, *args, **options):
        ###
        #
        user = User.objects.get(pk=options["user_pk"])
        sample = PIProject_Sample.objects.get(pk=options["sample_pk"])
        project = Projects.objects.get(pk=options["project_pk"])
        pipeline_leaf = SoftwareTreeNode.objects.get(pk=options["pipeline_pk"])
        pipeline_tree = SoftwareTree.objects.get(pk=pipeline_leaf.software_tree.pk)

        run = Run_Main_from_Leaf(
            user=user,
            sample=sample,
            project=project,
            pipeline_leaf=pipeline_leaf,
            pipeline_tree=pipeline_tree,
        )
        run.Submit()
