import pandas as pd

from pathogen_identification.models import (
    FinalReport,
    PIProject_Sample,
    Projects,
    RunMain,
)
from pathogen_identification.utilities.televir_parameters import TelevirParameters
from pathogen_identification.utilities.utilities_views import (  # ############################################
    EmptyRemapMain,
    RawReferenceUtils,
    ReportSorter,
    RunMainWrapper,
    final_report_best_cov_by_accid,
    recover_assembly_contigs,
)

project_pk = 123

project_main = Projects.objects.get(pk=123)
project_samples = PIProject_Sample.objects.filter(project=project_main)

all_data = []

for sample in project_samples:
    print(sample)
    final_report = FinalReport.objects.filter(
        sample=sample, run__project=project_main
    ).order_by("-coverage")
    unique_reports = final_report_best_cov_by_accid(final_report)
    #
    report_layout_params = TelevirParameters.get_report_layout_params(
        project_pk=project_main.pk
    )
    report_sorter = ReportSorter(sample, unique_reports, report_layout_params)
    sort_tree_plot_path = None
    if report_sorter.overlap_manager is not None:
        sort_tree_plot_path = report_sorter.overlap_manager.tree_plot_path_render
    sorted_reports = report_sorter.get_reports_compound()
    data = []
    print("reports: ", len(sorted_reports))
    for report_group in sorted_reports:
        group_name = report_group.name
        for report in report_group.group_list:
            data.append(
                {
                    "sample": report_sorter.sample.name,
                    "name": group_name,
                    "accid": report.accid,
                    "description": report.description,
                    "taxid": report.taxid,
                    "workflows": report.found_in,
                    "coverage": report.coverage,
                    "depth": report.depth,
                    "depthC": report.depthR,
                    "gaps": report.ngaps,
                    "mapped reads": report.mapped_reads,
                    "private_reads": report.private_reads,
                    "start prop": report.ref_proportion,
                    "mapped_proportion": report.mapped_proportion,
                    "windows_covered": report.windows_covered,
                    "error_rate": report.error_rate,
                    "class. success": report.classification_success,
                    "mapping success": report.mapping_success,
                }
            )
    df = pd.DataFrame(data)
    all_data.append(df)

df = pd.concat(all_data, ignore_index=True)
