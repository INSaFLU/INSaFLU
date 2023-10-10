import pandas as pd
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from pathogen_identification.models import (
    FinalReport,
    ParameterSet,
    PIProject_Sample,
    Projects,
    RawReference,
    RunMain,
)
from pathogen_identification.utilities.utilities_pipeline import (
    PipelineTree,
    SoftwareTreeUtils,
    Utils_Manager,
)


def read_parameters(row):
    try:
        params = pd.read_csv(row.params_file_path, sep="\t")
        params["run"] = row.run_name
        params["run_id"] = row.run_id
        params["sample"] = row.sample_name_x
        params["sample_id"] = row.sample_id
        params["project"] = row.project_name
        params["project_id"] = row.project_id
        return params

    except FileNotFoundError:
        return None


def collect_parameters(queryset_df: pd.DataFrame) -> pd.DataFrame:
    all_parameters = []
    unique_runs_df = queryset_df.drop_duplicates(subset=["run_id"])
    for index, row in unique_runs_df.iterrows():
        params = read_parameters(row)
        if params is not None:
            all_parameters.append(params)
    all_parameters = pd.concat(all_parameters).rename(columns={"sample": "sample_name"})
    return all_parameters


def collect_parameters_project(project: Projects):
    utils = Utils_Manager()
    technology = project.technology
    user = project.owner
    pipeline_utils = Utils_Manager()

    samples = PIProject_Sample.objects.filter(project=project)
    project_params = []
    tree_utils = SoftwareTreeUtils(user=project.owner, project=project)

    for sample in samples:
        tree_utils.set_sample(sample=sample)

        parameterset = ParameterSet.objects.filter(
            project=project, sample=sample, status=ParameterSet.STATUS_FINISHED
        )
        mainruns = RunMain.objects.filter(parameter_set__in=parameterset)

        for run in mainruns:
            ps = run.parameter_set
            node_leaves = ps.get_leaf_descendants()
            leaf_node = node_leaves[0]
            leaf_node_index = leaf_node.index
            params = pipeline_utils.get_leaf_parameters(leaf_node)
            params["leaves"] = leaf_node_index

            if params is None:
                print("no params found for leaf: ")

            params["run"] = run.name
            params["run_id"] = run.pk
            params["sample"] = ps.sample.name
            params["sample_id"] = ps.sample.pk
            params["sample_name"] = ps.sample.name
            params["project"] = project.name
            params["project_id"] = project.pk

            project_params.append(params)

    if len(project_params):
        project_params = pd.concat(project_params)

    return project_params


def projects_params_summary(project_query):
    all_params = []

    for project in project_query:
        project_params = collect_parameters_project(project)
        if len(project_params):
            all_params.append(project_params)

    if len(all_params):
        all_params = pd.concat(all_params)
    return all_params


def projects_reports_summary(project_query):
    all_reports = []

    for project in project_query:
        samples = PIProject_Sample.objects.filter(project=project)
        for sample in samples:
            ps = ParameterSet.objects.filter(
                sample=sample, project=project, status=ParameterSet.STATUS_FINISHED
            )
            mainruns = RunMain.objects.filter(parameter_set__in=ps)
            for mainrun in mainruns:
                final_report = FinalReport.objects.filter(run=mainrun)

                if len(final_report) == 0:
                    continue
                final_report_dict = pd.DataFrame(final_report.values())
                final_report_dict["run"] = mainrun.name
                final_report_dict["run_id"] = mainrun.pk
                final_report_dict["sample"] = sample.name
                final_report_dict["sample_name"] = sample.name
                final_report_dict["sample_id"] = sample.pk
                final_report_dict["project"] = project.name
                final_report_dict["project_id"] = project.pk
                final_report_dict["runtime"] = mainrun.runtime
                all_reports.append(final_report_dict)

    if len(all_reports):
        all_reports = pd.concat(all_reports)

    return all_reports


def project_references(project_query):
    all_params = []

    for project in project_query:
        project_params = collect_references_project(project)
        if len(project_params):
            all_params.append(project_params)

    if len(all_params):
        all_params = pd.concat(all_params)
    return all_params


def collect_references_project(project: Projects):
    samples = PIProject_Sample.objects.filter(project=project)
    project_params = []

    for sample in samples:
        parameterset = ParameterSet.objects.filter(
            project=project, sample=sample, status=ParameterSet.STATUS_FINISHED
        )
        mainruns = RunMain.objects.filter(parameter_set__in=parameterset)

        for run in mainruns:
            reference_mapped = RawReference.objects.filter(run=run)
            reference_mapped = pd.DataFrame(reference_mapped.values())

            ps = run.parameter_set

            if reference_mapped is None:
                print("no references found for run: ", run.pk)

            reference_mapped["run"] = run.name
            reference_mapped["run_id"] = run.pk
            reference_mapped["sample"] = ps.sample.name
            reference_mapped["sample_id"] = ps.sample.pk
            reference_mapped["sample_name"] = ps.sample.name
            reference_mapped["project"] = project.name
            reference_mapped["project_id"] = project.pk

            project_params.append(reference_mapped)

    if len(project_params):
        project_params = pd.concat(project_params)

    return project_params


class Command(BaseCommand):
    help = "deploy run"

    def add_arguments(self, parser):
        parser.add_argument(
            "-u",
            "--user",
            type=str,
            help="user to get info from",
        )
        parser.add_argument(
            "-p",
            "--project",
            default="all",
            type=str,
            help="project to target",
        )

        parser.add_argument(
            "--tech",
            default="nanopore",
            type=str,
            help="project to target",
        )

    def handle(self, *args, **options):
        ###
        #

        user = User.objects.get(username=options.get("user"))
        project_get = options.get("project")

        if options.get("project") == "all":
            projects_query = Projects.objects.filter(owner=user)

        else:
            projects_query = Projects.objects.filter(
                name=options.get("project"), owner=user
            )

        # local_tree= self.generate_hardcoded_pipeline_tree(technology=technology)
        all_parameters = projects_params_summary(projects_query)
        reports_df = projects_reports_summary(projects_query)
        references_df = project_references(projects_query)

        ####
        reports_df.to_csv(
            f"{project_get}_u{user.pk}_reports.tsv", index=False, sep="\t"
        )

        ####
        all_parameters.to_csv(
            f"{project_get}_u{user.pk}_parameters.tsv", index=False, sep="\t"
        )
        #

        references_df.to_csv(
            f"{project_get}_u{user.pk}_references.tsv", index=False, sep="\t"
        )
