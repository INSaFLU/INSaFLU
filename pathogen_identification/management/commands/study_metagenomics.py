from abc import ABC, abstractmethod
from typing import Any, List, Optional

from pathogen_identification.models import FinalReport, PIProject_Sample, RawReference


class SampleWrapper:

    sample: PIProject_Sample

    def __init__(self, sample):
        self.sample = sample

    @property
    def sample_class(self):
        if self.sample.sample.name is None:
            return "NONE"
        if "UPIP" in self.sample.sample.name:
            return "UPIP"
        elif "RPIP" in self.sample.sample.name:
            return "RPIP"
        else:
            return "OTHER"


class SamplesBenchCollection:
    name: str
    samples: List[SampleWrapper]

    def __init__(self, name):
        self.name = name
        self.samples = []

    def sample_add(self, sample):
        self.samples.append(SampleWrapper(sample))

    def sample_register(self, sample: PIProject_Sample):

        if sample.sample.name is None:
            return

        self.sample_add(sample)


class SampleCurator:

    collection: SamplesBenchCollection

    def __init__(self, project_id: int, sample_name: str):
        self.project_id = project_id
        self.samples = self.get_samples_by_project(sample_name)
        pass

    def get_samples_by_project(self, pattern=None) -> SamplesBenchCollection:
        project_samples = PIProject_Sample.objects.filter(project_id=self.project_id)
        collection = SamplesBenchCollection("project")

        if pattern:
            project_samples = project_samples.filter(sample__name__icontains=pattern)

        for sample in project_samples:
            collection.sample_register(sample)

        return collection


def match_name_score(name: str, reference) -> float:
    name_list = name.split(" ")

    score = 0
    for ix, string in enumerate(name_list):
        string_until_now = " ".join(name_list[: ix + 1])
        if reference.description is None:
            continue
        if string_until_now in reference.description:
            score += 1

    score = score / len(name_list)
    return score


class Hit:
    name: str
    taxid: Optional[int]
    raw_reference_id_list: List[int]

    def __init__(
        self,
        name,
        taxid: Optional[int] = None,
        raw_reference_id_list: List[int] = [],
        reports_list: List[int] = [],
    ):
        self.taxid = taxid
        self.raw_reference_id_list = raw_reference_id_list
        self.reports_list = reports_list
        self.name = name

        self.name_similarity = 0
        if len(self.raw_reference_id_list) > 0:
            self.name_similarity = match_name_score(
                self.name, RawReference.objects.get(pk=self.raw_reference_id_list[0])
            )

        if len(self.reports_list) > 0:
            self.name_similarity = match_name_score(
                self.name, FinalReport.objects.get(pk=self.reports_list[0])
            )

    @property
    def raw_reference_samples(self) -> List[SampleWrapper]:
        refs = (
            RawReference.objects.filter(id__in=self.raw_reference_id_list)
            .distinct("run__sample__id")
            .exclude(run__sample=None)
        )
        return [SampleWrapper(ref.run.sample) for ref in refs]

    @property
    def reported_samples(self) -> List[SampleWrapper]:
        refs = (
            FinalReport.objects.filter(id__in=self.reports_list)
            .distinct("run__sample__id")
            .exclude(run__sample=None)
        )
        return [SampleWrapper(ref.run.sample) for ref in refs]

    @property
    def reported_samples_classes(self) -> List[str]:
        return [sample.sample_class for sample in self.reported_samples]

    @property
    def intermediate_samples_classes(self) -> List[str]:
        return [sample.sample_class for sample in self.raw_reference_samples]


class HitFactory:

    def __init__(self, project_id: int, collection: SamplesBenchCollection):
        self.project_id = project_id
        self.sample = None
        self.collection = collection

    def hit_by_name(self, name: str) -> Hit:

        reference_hits = (
            RawReference.objects.filter(run__sample__in=self.collection.samples)
            .distinct("run__id")
            .exclude(run=None)
        )
        reference_hits_list = []

        if reference_hits.exists():
            reference_hits = reference_hits.all()
            reference_hits = list(reference_hits)
            reference_hits.sort(key=lambda x: match_name_score(name, x), reverse=True)
            if match_name_score(name, reference_hits[0]) == 1:
                reference_hits_list = [
                    x.pk for x in reference_hits if match_name_score(name, x) == 1
                ]
            else:
                reference_hits_list = [
                    x.pk for x in reference_hits if match_name_score(name, x) > 0
                ]

        report_hits = (
            FinalReport.objects.filter(run__sample__in=self.collection.samples)
            .distinct("run__id")
            .exclude(run=None)
        )

        report_hits_list = []
        if report_hits.exists():
            report_hits = report_hits.all()
            report_hits = list(report_hits)
            report_hits.sort(key=lambda x: match_name_score(name, x), reverse=True)
            if match_name_score(name, report_hits[0]) == 1:
                report_hits_list = [
                    x.pk for x in report_hits if match_name_score(name, x) == 1
                ]
            else:
                report_hits_list = [
                    x.pk for x in report_hits if match_name_score(name, x) > 0
                ]

        return Hit(
            name,
            raw_reference_id_list=reference_hits_list,
            reports_list=report_hits_list,
        )


import pandas as pd


def df_report_analysis(analysis_df_filename, project_id: int):

    # Python
    sheet_name = "Results"
    df = pd.read_excel(
        analysis_df_filename,
        sheet_name=sheet_name,
        header=1,
        index_col=0,
        engine="openpyxl",
    )

    bact_results = df[df["Class"] == "bacterial"]

    new_table = []

    for ix, row in bact_results.iterrows():
        print(f"############## {ix} ##############")
        sample_name = row["Sample_ID"]
        curator = SampleCurator(project_id, sample_name)
        hit_factory = HitFactory(project_id, curator.collection)

        expected_hit = hit_factory.hit_by_name(sample_name)

        row["name_similarity"] = expected_hit.name_similarity
        row["reported_samples"] = expected_hit.reported_samples_classes
        row["intermediate_samples"] = expected_hit.intermediate_samples_classes

        new_table.append(row)

    new_df = pd.DataFrame(new_table)

    return new_df


############################################################################################################
import os
from typing import List

import pandas as pd
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from managing_files.models import ProcessControler
from pathogen_identification.constants_settings import ConstantsSettings as CS
from pathogen_identification.models import FinalReport, PIProject_Sample, Projects
from pathogen_identification.templatetags.report_colors import flag_false_positive
from pathogen_identification.utilities.explify_merge import (
    get_illumina_found,
    merge_panels,
    process_televir,
    read_panel,
)
from pathogen_identification.utilities.utilities_general import get_project_dir
from utils.process_SGE import ProcessSGE


class Command(BaseCommand):
    help = "deploy run"

    def add_arguments(self, parser):
        parser.add_argument(
            "--project_id",
            type=int,
            help="televir project to be merged (pk)",
        )

        parser.add_argument(
            "--report",
            type=str,
            help="televir panel report",
        )

        parser.add_argument(
            "--output",
            type=str,
            help="output file",
        )

    def handle(self, *args, **options):
        ###
        # SETUP
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        project = Projects.objects.filter(
            is_deleted=False, pk=options["project_id"]
        ).first()

        report = options["report"]

        df = df_report_analysis(report, project.pk)

        df.to_excel(options["output"], index=False)
