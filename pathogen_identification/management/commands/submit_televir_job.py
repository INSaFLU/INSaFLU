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
            "--user_id",
            type=int,
            help="user deploying the run (pk)",
        )
        parser.add_argument(
            "--sample_id",
            type=int,
            help="sample to be run (pk)",
        )

        parser.add_argument(
            "--project_id",
            type=int,
            help="project to be run (pk)",
        )

        parser.add_argument(
            "--pipeline_id",
            type=int,
            help="pipeline to be run (pk )",
        )

        parser.add_argument(
            "-o",
            "--outdir",
            type=str,
            help="output directory",
        )

    def handle(self, *args, **options):
        ###
        #
        user = User.objects.get(pk=options["user_id"])
        sample = PIProject_Sample.objects.get(pk=options["sample_id"])
        project = Projects.objects.get(pk=options["project_id"])
        pipeline_leaf = SoftwareTreeNode.objects.get(pk=options["pipeline_id"])
        pipeline_tree = SoftwareTree.objects.get(pk=pipeline_leaf.software_tree.pk)

        run = Run_Main_from_Leaf(
            user=user,
            input_data=sample,
            project=project,
            pipeline_leaf=pipeline_leaf,
            pipeline_tree=pipeline_tree,
            odir=options["outdir"],
        )
        run.Submit()
