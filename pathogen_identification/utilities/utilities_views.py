import os
from typing import List

import pandas as pd

from pathogen_identification.models import FinalReport, ReferenceMap_Main
from pathogen_identification.utilities.utilities_general import (
    infer_run_media_dir,
    simplify_name,
)

# import Django BaseManager
# from django.db.models import BaseManager


def set_control_reports(project_pk: int):
    """
    set control reports
    """

    try:
        project = Projects.objects.get(pk=project_pk)

        control_samples = PIProject_Sample.objects.filter(
            project=project, is_control=True
        )

        control_reports = FinalReport.objects.filter(sample__in=control_samples)

        control_report_taxids = control_reports.values_list("taxid", flat=True)
        control_report_taxids_set = set(control_report_taxids)
        print(control_report_taxids_set)

        other_reports = FinalReport.objects.filter(sample__project=project).exclude(
            sample__in=control_samples
        )

        for sample_report in other_reports:
            if sample_report.taxid in control_report_taxids_set:
                sample_report.control_flag = FinalReport.CONTROL_FLAG_PRESENT
            else:
                sample_report.control_flag = FinalReport.CONTROL_FLAG_NONE

            sample_report.save()

        for report in control_reports:
            report.control_flag = FinalReport.CONTROL_FLAG_NONE
            report.save()

    except Exception as e:
        print(e)
        pass
