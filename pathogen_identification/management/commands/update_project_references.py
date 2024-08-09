import traceback
from datetime import date

from django.conf import settings
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from managing_files.models import ProcessControler
from pathogen_identification.models import PIProject_Sample, Projects
from pathogen_identification.utilities.tree_deployment import TreeProgressGraph
from pathogen_identification.utilities.utilities_views import (
    RawReferenceUtils,
    calculate_reports_overlaps,
)
from utils.process_SGE import ProcessSGE


class Command(BaseCommand):
    help = "deploy run"

    def add_arguments(self, parser):

        parser.add_argument(
            "--project_id",
            type=int,
            help="project to be run (pk)",
        )

    def handle(self, *args, **options):

        project_id = options["project_id"]
        project = Projects.objects.get(pk=project_id)

        all_project_samples = PIProject_Sample.objects.filter(
            is_deleted=False, project__is_deleted=False, project=project
        )
        process_SGE = ProcessSGE()
        process_controler = ProcessControler()

        process_SGE.set_process_controler(
            project.owner,
            process_controler.get_name_update_televir_project(
                project_id=project.pk,
            ),
            ProcessControler.FLAG_RUNNING,
        )

        print(f"Updating references for {all_project_samples.count()} samples")

        try:
            for sample in all_project_samples:
                reference_utils = RawReferenceUtils(sample)

                sample_runs = reference_utils.filter_runs().count()
                if sample_runs == 0:
                    print(f"Sample {sample.name} has no runs")
                    continue

                calculate_reports_overlaps(sample, force=True)
                _ = reference_utils.retrieve_compound_references()

                graph_progress = TreeProgressGraph(sample)
                graph_progress.generate_graph()

            process_SGE.set_process_controler(
                project.owner,
                process_controler.get_name_update_televir_project(
                    project_id=project.pk,
                ),
                ProcessControler.FLAG_FINISHED,
            )

            project.updated_version = settings.APP_VERSION_NUMBER
            project.save()

        except Exception as e:
            print(f"Error: {e}")
            print(traceback.format_exc())
            process_SGE.set_process_controler(
                project.owner,
                process_controler.get_name_update_televir_project(
                    project_id=project.pk,
                ),
                ProcessControler.FLAG_ERROR,
            )
            return
