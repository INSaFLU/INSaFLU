from typing import Any, List, Optional

import pandas as pd
from django.core.management.base import BaseCommand
from django.db.models import Q

from pathogen_identification.models import (
    FinalReport,
    PIProject_Sample,
    Projects,
    RawReference,
)


def process_class(r2, maxt=6):
    """
    Process classification results.
    """
    r2 = r2.drop_duplicates(subset=["taxid"], keep="first")
    r2 = r2.reset_index(drop=True)
    # r2 = r2.sort_values("counts", ascending=False)
    taxids_tokeep = []
    nr2 = []
    if "length" in r2.columns:
        r2["length"] = r2["length"].astype(int)
        r2c = r2.copy().sort_values("length", ascending=False)
        for i in range(r2.shape[0]):
            if i < maxt:
                taxids_tokeep.append(r2c.taxid[i])
                nr2.append(r2c.loc[i])
                if r2.taxid.tolist()[i] not in taxids_tokeep:
                    taxids_tokeep.remove(r2.taxid[i])
                    nr2.append(r2.loc[i])
            else:
                break
    else:
        r2 = r2.head(maxt)
    if len(nr2):
        r2 = pd.concat(nr2, axis=1).T
    return r2


def get_source(row):
    if row.counts_x > 0 and row.counts_y > 0:
        return 3
    elif row.counts_x > 0:
        return 1
    elif row.counts_y > 0:
        return 2


def descriptor_sources(fd, r1_raw, r2_raw):
    if len(r1_raw) == 0:
        fd["source"] = 2
    elif len(r2_raw) == 0:
        fd["source"] = 1
    else:
        fd["source"] = fd.apply(get_source, axis=1)
    return fd


def descriptor_counts(fd, r1_raw, r2_raw):
    if len(r1_raw) == 0 or len(r2_raw) == 0:
        if "counts" not in fd.columns:
            fd["counts"] = fd["counts_x"]
        if "counts_y" in fd.columns:
            fd["counts"] = fd["counts_y"] + fd["counts_x"]
        return fd
    else:
        fd["counts"] = fd.apply(get_counts, axis=1)
    return fd


def get_counts(row):
    if row.counts_x > 0 and row.counts_y > 0:
        return f"{int(row.counts_x)} / {int(row.counts_y)}"
    elif row.counts_x > 0:
        return str(row.counts_x)
    elif row.counts_y > 0:
        return str(row.counts_y)


def descriptor_description_remove(df: pd.DataFrame):
    if "description_x" in df.columns:
        df.drop("description_x", axis=1, inplace=True)
    if "description_y" in df.columns:
        df.drop("description_y", axis=1, inplace=True)
    if "description" in df.columns:
        df.drop("description", axis=1, inplace=True)
    return df


def merge_classes(r1: pd.DataFrame, r2: pd.DataFrame, maxt=6, exclude="phage"):
    """
    merge tables of taxids to columns.
    """
    ###
    if "description" in r1.columns:
        r1 = r1[~r1.description.str.contains(exclude)]
    if "description" in r2.columns:
        r2 = r2[~r2.description.str.contains(exclude)]
    ###
    r1 = r1[["taxid", "counts"]]
    r1 = r1.sort_values("counts", ascending=False)
    r2 = r2.sort_values("counts", ascending=False)
    r1_raw = r1.copy()
    r2_raw = r2.copy()
    ###
    if len(r2) > 0 and len(r1) > 0:
        full_descriptor = pd.merge(r1, r2, on="taxid", how="outer")
    elif len(r2) > 0:
        full_descriptor = r2.copy()
    else:
        full_descriptor = r1.copy()
    full_descriptor = full_descriptor.fillna(0)
    ###
    if len(r2) and len(r1):
        r1.taxid = r1.taxid.astype(str)
        r2.taxid = r2.taxid.astype(str)
        shared = pd.merge(r1, r2, on=["taxid"], how="inner").sort_values(
            "counts_x", ascending=False
        )
        maxt = maxt - shared.shape[0]
        if maxt <= 0:
            r1 = shared
        else:
            r2 = (
                pd.merge(r2, shared, indicator=True, how="outer")
                .query('_merge=="left_only"')
                .drop("_merge", axis=1)
            )
            r2 = process_class(r2, maxt=maxt)
            r1p = r1[~r1.taxid.isin(shared.taxid)]
            r1 = (
                pd.concat([shared, r2, r1p.head(maxt - r2.shape[0])], axis=0)
                .drop_duplicates(subset=["taxid"], keep="first")
                .reset_index(drop=True)
            )

    elif len(r2) == 0:
        r1 = r1.head(maxt)
    elif len(r1) == 0:
        r1 = r2.head(maxt)
    ###

    full_descriptor["taxid"] = full_descriptor["taxid"].astype(int)
    full_descriptor = descriptor_description_remove(full_descriptor)
    full_descriptor = descriptor_sources(full_descriptor, r1_raw, r2_raw)
    full_descriptor = descriptor_counts(full_descriptor, r1_raw, r2_raw)

    r1["taxid"] = r1.taxid.astype(int)
    merged_final = full_descriptor[full_descriptor.taxid.isin(r1.taxid.to_list())]
    # get taxid index in r1
    merged_final["taxid_index"] = merged_final.apply(
        lambda x: r1[r1.taxid == x.taxid].index[0], axis=1
    )
    merged_final = merged_final.sort_values("taxid_index")

    merged_final["source"] = merged_final.source.apply(
        lambda x: ["none", "reads", "contigs", "reads/contigs"][x]
    )
    merged_final = merged_final[["taxid", "counts", "source"]]

    return merged_final, full_descriptor


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

    @property
    def samples_televir(self):
        return [sample.sample for sample in self.samples]


from django.db.models import Q


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
            project_samples = project_samples.filter(
                Q(sample__name__icontains=pattern)
                | Q(sample__name__icontains=pattern.replace("-", "_"))
            )

        for sample in project_samples:
            collection.sample_register(sample)

        return collection

    def set_collection(self, name_pattern: str):

        print("##### Setting collection")
        print("Pattern", name_pattern)
        self.collection = self.get_samples_by_project(name_pattern)
        print(f"Collection has {len(self.collection.samples)} samples")
        print(
            "sample classes",
            [sample.sample_class for sample in self.collection.samples],
        )


def name_translator(name: str) -> str:

    if "EBV" in name:
        return "human herpesvirus 4"

    if "polyomavirus" in name:
        return "polyomavirus"

    if "Metapneumovírus" in name:
        return "human metapneumovirus"

    if "SARS-CoV-2" in name:
        return "Severe acute respiratory syndrome coronavirus"

    if "Varicella" in name:
        return "human herpesvirus 3"

    if "Cytomegalovirus" in name:
        return "human herpesvirus 5"

    if "Rinovírus" in name or "Rinovíru" in name:
        return "rhinovirus"

    if "Epstein-Barr virus (EBV)" in name:
        return "human herpesvirus 4"

    if "Herpes simplex virus 1" in name:
        return "human herpesvirus 1"

    if "RSV" in name:
        return "respiratory syncytial virus"

    return name


def match_name_score(name: str, reference: RawReference) -> float:
    name = name_translator(name)
    name_list = name.lower().split(" ")
    if reference.description is None:
        return 0

    description_lower = reference.description.lower()

    scores = []

    score = 0
    for ix, string in enumerate(name_list):
        string_until_now = " ".join(name_list[: ix + 1])

        if string_until_now in description_lower:
            score += 1

    score = score / len(name_list)
    scores.append(score)

    if "betaherpes" in name:
        name = name.replace("betaherpes", "herpes")

    if "BK polyomavirus" in name:
        name = "Betapolyomavirus hominis"

    score = 0
    for ix, string in enumerate(name_list):
        string_until_now = " ".join(name_list[: ix + 1])

        if string_until_now in description_lower:
            score += 1

    score = score / len(name_list)
    scores.append(score)

    return max(scores)


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
    def raw_references(self) -> List[RawReference]:
        return RawReference.objects.filter(id__in=self.raw_reference_id_list)

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
            RawReference.objects.filter(run__sample__in=self.collection.samples_televir)
            .exclude(run=None)
            .exclude(accid="-")
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

        report_hits = FinalReport.objects.filter(
            run__sample__in=self.collection.samples_televir
        ).exclude(run=None)

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


def get_hit_best_reference(hit: Hit) -> Optional[RawReference]:
    if len(hit.raw_reference_id_list) == 0:
        return None

    hit_references = RawReference.objects.filter(
        id__in=hit.raw_reference_id_list
    ).exclude(accid="-")
    hit_references = [(x, determine_raw_ref_index(x)) for x in hit_references]
    hit_references.sort(key=lambda x: x[1])

    return hit_references[0][0]


def get_hit_best_classifier_reference(
    hit: Hit, classifier: str, panel: str
) -> Optional[RawReference]:
    if len(hit.raw_reference_id_list) == 0:
        return None

    hit_references = RawReference.objects.filter(
        id__in=hit.raw_reference_id_list,
        run__read_classification=classifier,
        run__sample__name__icontains=panel,
    ).exclude(accid="-")

    if hit_references.exists():
        hit_references = [(x, determine_raw_ref_index(x)) for x in hit_references]
        hit_references.sort(key=lambda x: x[1])

        return hit_references[0][0]
    else:
        return None


def ref_index_by_pk(ref: RawReference):
    run = ref.run

    other_references = RawReference.objects.filter(run=run).order_by("id")
    other_references = list(other_references)
    for ix, other_ref in enumerate(other_references):
        if other_ref.pk == ref.pk:
            return ix

    return -1


def ref_index_by_counts(ref: RawReference):
    run = ref.run
    other_references = RawReference.objects.filter(run=run).exclude(accid="-")
    reference_counts = {
        ref.pk: ref.read_counts for ix, ref in enumerate(other_references)
    }
    reference_counts = {
        pk: read_counts
        for pk, read_counts in reference_counts.items()
        if read_counts is not None
    }
    reference_counts = {x: float(y) for x, y in reference_counts.items()}
    sorted_pks = sorted(reference_counts, key=reference_counts.get, reverse=True)

    for ix, pk in enumerate(sorted_pks):
        if pk == ref.pk:
            return ix
    return -1


def determine_raw_ref_index(ref: RawReference):

    # return ref_index_by_pk(ref)
    return ref_index_by_counts(ref)


def get_frames(ref: RawReference):
    contigs_df = []
    reads_df = []
    other_references = RawReference.objects.filter(run=ref.run).exclude(accid="-")
    for ref in other_references:
        if ref.contig_counts != "0":
            contigs_df.append([int(ref.taxid), int(float(ref.contig_counts))])

        if ref.read_counts != "0":
            reads_df.append([int(ref.taxid), int(float(ref.read_counts))])
    contigs_df = pd.DataFrame(contigs_df, columns=["taxid", "counts"])
    reads_df = pd.DataFrame(reads_df, columns=["taxid", "counts"])
    return contigs_df, reads_df


def determine_ref_televir_sort_index(ref: RawReference):
    """
    Determine the index of the reference in the televir sort."""
    contigs_df, reads_df = get_frames(ref)
    merged_final, full_descriptor = merge_classes(reads_df, contigs_df, maxt=3000)
    print(ref.accid)
    print(merged_final.head())
    print(ref.taxid in contigs_df.taxid)
    print(ref.taxid in reads_df.taxid)
    print(int(ref.taxid), merged_final.taxid)
    print(ref.taxid)
    ref_accid_index = merged_final[merged_final["taxid"] == int(ref.taxid)].index[0]
    return ref_accid_index


def determine_raw_ref_counts(ref: RawReference):
    return ref.read_counts


def determine_reads_stats(ref: RawReference, samples_list: List[SampleWrapper] = []):

    run = ref.run
    sample = run.sample

    accid_use = ref.accid
    accid_use = accid_use.split(".")[0]

    reports = FinalReport.objects.filter(
        taxid__icontains=ref.taxid,
        sample__in=[sample.sample for sample in samples_list],
    ).order_by("-coverage")
    mapped_reads = 0

    print(f"Found {reports.count()} reports for {ref.accid} in {sample.name}")
    print([f.mapped_reads for f in reports])

    if reports.exists():
        mapped_reads = reports.first().mapped_reads
        if mapped_reads is None:
            mapped_reads = 0
    else:
        mapped_reads = 0

    return mapped_reads


def determine_reads_map_qual(ref: RawReference, samples_list: List[SampleWrapper] = []):

    run = ref.run
    sample = run.sample

    accid_use = ref.accid
    accid_use = accid_use.split(".")[0]

    reports = FinalReport.objects.filter(
        taxid__icontains=ref.taxid,
        sample__in=[sample.sample for sample in samples_list],
    ).order_by("-coverage")
    error_rate = 0

    if reports.exists():
        try:
            error_rate = reports.first().error_rate
        except:
            error_rate = None

        if error_rate is None:
            error_rate = "-"
    else:
        error_rate = "-"

    return error_rate


def df_report_analysis(analysis_df_filename, project_id: int):

    # Python

    df = pd.read_csv(analysis_df_filename, sep="\t")

    bact_results = df[df["Class"] == "bacterial"]

    new_table = []

    for ix, row in df.iterrows():
        sample_name = str(row["Sample_ID"])
        hitname = str(row["Reporting Name"])
        org_class = str(row["Class"])
        curator = SampleCurator(project_id, sample_name)
        curator.set_collection(sample_name)
        hit_factory = HitFactory(project_id, curator.collection)

        expected_hit = hit_factory.hit_by_name(hitname)
        best_mapping = get_hit_best_reference(expected_hit)
        mapped = 0
        run_pk = -1
        sample = None

        best_mapping_rank = -1

        if best_mapping is not None:

            best_mapping_rank = determine_raw_ref_index(best_mapping)
            mapped = best_mapping.status
            run_pk = best_mapping.run.pk
            sample = best_mapping.run.sample

        new_row = {
            "sample": sample_name,
            "Class": org_class,
            "hitname": hitname,
            "name_similarity": expected_hit.name_similarity,
            "reported_samples": "/".join(expected_hit.reported_samples_classes),
            "intermediate_samples": "/".join(expected_hit.intermediate_samples_classes),
            "position": best_mapping_rank + 1,
            "mapped": mapped,
            "run_pk": run_pk,
            "sample_pk": sample.pk if sample is not None else -1,
        }

        new_row["mapped_reads UPIP / RPIP"] = "0 / 0"
        new_row["map quality UPIP / RPIP"] = "- / -"

        if best_mapping is not None:
            class_mapped_reads = {
                "UPIP": 0,
                "RPIP": 0,
            }

            class_map_quality = {
                "UPIP": "-",
                "RPIP": "-",
            }

            for sample in expected_hit.raw_reference_samples:
                sample_class = sample.sample_class
                class_mapped_reads[sample_class] += determine_reads_stats(
                    best_mapping, [sample]
                )

                class_map_quality[sample_class] = determine_reads_map_qual(
                    best_mapping, [sample]
                )

            new_row["mapped_reads UPIP / RPIP"] = (
                f"{class_mapped_reads['UPIP']} / {class_mapped_reads['RPIP']}"
            )

            new_row["map quality UPIP / RPIP"] = (
                f"{class_map_quality['UPIP']} / {class_map_quality['RPIP']}"
            )

        for classifier in ["kraken2", "centrifuge"]:

            for panel in ["UPIP", "RPIP"]:

                best_mapping = get_hit_best_classifier_reference(
                    expected_hit, classifier, panel
                )

                if best_mapping is not None:
                    best_mapping_rank = determine_raw_ref_index(best_mapping)
                    best_televir_sort_rank = determine_ref_televir_sort_index(
                        best_mapping
                    )
                    new_row[f"{classifier}_{panel}_position"] = best_mapping_rank + 1
                    new_row[f"{classifier}_{panel}_televir_sort_position"] = (
                        best_televir_sort_rank + 1
                    )
                else:
                    new_row[f"{classifier}_{panel}_position"] = -1
                    new_row[f"{classifier}_{panel}_televir_sort_position"] = -1

        new_table.append(new_row)

    new_df = pd.DataFrame(new_table)

    new_df["sample_pk"] = new_df.sample_pk.astype(int)
    new_df["run_pk"] = new_df.run_pk.astype(int)

    return new_df


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

        project = Projects.objects.filter(
            is_deleted=False, pk=options["project_id"]
        ).first()

        report = options["report"]

        df = df_report_analysis(report, project.pk)

        df.to_csv(options["output"], index=False, sep="\t")
