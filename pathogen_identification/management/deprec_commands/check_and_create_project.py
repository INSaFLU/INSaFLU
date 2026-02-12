import re
from typing import List, Optional

import pandas as pd
from django.contrib.auth.models import User
from django.db import transaction

from managing_files.models import Sample
from pathogen_identification.utilities.utilities_general import simplify_name_lower


def check_pattern_in_name(name, pattern):
    return bool(re.search(pattern, name))


def check_pattern_in_column(df, column, pattern) -> bool:
    return df[column].apply(lambda x: check_pattern_in_name(pattern, x)).any()


## Study samples
study_samples = "/insaflu_web/metagenomics_reads/metagenomics_study_samples.lst"
study_samples = pd.read_csv(study_samples, sep="\t", header=None).rename(
    columns={0: "name"}
)
study_samples["name_simple"] = study_samples["name"].apply(simplify_name_lower)

## Sample all
sample_all = Sample.objects.all().values()
sample_all = pd.DataFrame(sample_all)
sample_all = sample_all.dropna(subset=["name"])
sample_all = sample_all[["name", "id", "creation_date"]]
sample_all["name_simple"] = sample_all["name"].apply(simplify_name_lower)
sample_all.sort_values("creation_date", ascending=False, inplace=True)
sample_all["study"] = sample_all["name_simple"].apply(
    lambda x: check_pattern_in_column(study_samples, "name_simple", x)
)


## Check if study samples are in the database by checking if db name in study samples
def samples_with_pattern_in_column(df, column, pattern):
    return df[column].apply(lambda x: check_pattern_in_name(x, pattern))


def samples_with_pattern_in_column_ids(df, column, pattern) -> List[int]:
    ids = df[samples_with_pattern_in_column(df, column, pattern)]["id"]
    if ids.empty:
        return []
    return list(ids.astype(int))


study_samples["in_db"] = study_samples["name_simple"].apply(
    lambda x: samples_with_pattern_in_column_ids(sample_all, "name_simple", x)
)

study_samples_in_db = study_samples[study_samples["in_db"].apply(lambda x: len(x) > 0)]

print(study_samples_in_db.shape)
print(study_samples_in_db.head())

print(
    "Found {} samples out of {} in the database".format(
        study_samples_in_db.shape[0], study_samples_in_db.shape[0]
    )
)

### Check if there are multiple ids for the same sample
multiple_samples = study_samples_in_db[study_samples_in_db["in_db"].apply(len) > 1]
print(multiple_samples.shape)
print(multiple_samples.head())

print("Found {} samples with multiple ids".format(multiple_samples.shape[0]))


#### get sample names
def sample_name(id: int):
    return Sample.objects.get(id=id).name


def sample_names(ids: List[int]) -> List[str]:
    return [sample_name(x) for x in ids]


study_samples_in_db["names"] = study_samples_in_db["in_db"].apply(sample_names)

### select sample to keepd


def sample_pass(id: int):
    try:
        sample = Sample.objects.get(id=id)
        return sample.has_files
    except:
        return False


def select_sample(sample_list: List[int]) -> Optional[int]:

    samples_pass = [x for x in sample_list if sample_pass(x)]
    if not samples_pass:
        return None
    return samples_pass[0]


study_samples_in_db["id_select"] = study_samples_in_db["in_db"].apply(select_sample)
study_samples_in_db.head()

print(study_samples_in_db.shape)
print(study_samples_in_db.head())

print(
    "Found {} samples with data out of {} in the database".format(
        study_samples_in_db["id_select"].notnull().sum(), study_samples_in_db.shape[0]
    )
)

### Create project
from pathogen_identification.models import PIProject_Sample, Projects
from settings.constants_settings import ConstantsSettings

samples = [Sample.objects.get(pk=x) for x in study_samples_in_db["id_select"].dropna()]
samples_technology = [x.type_of_fastq for x in samples]
if len(list(set(samples_technology))) > 1:
    raise ValueError("Different technologies in samples")

technology = samples_technology[0]
if technology == Sample.TYPE_OF_FASTQ_illumina:
    technology = ConstantsSettings.TECHNOLOGY_illumina
else:
    technology = ConstantsSettings.TECHNOLOGY_minion

user_pk = 5
user = User.objects.get(pk=user_pk)
project_name = "metagenomics_study"
project_description = "Metagenomics study samples"

try:
    project = Projects.objects.get(name=project_name)
    print("Project already exists")
except Projects.DoesNotExist:
    project = Projects.objects.create(
        name=project_name,
        description=project_description,
        owner=user,
        technology=technology,
    )

for sample in samples:
    try:
        project_sample = PIProject_Sample.objects.get(project=project, sample=sample)
        print("Project sample already exists")
    except PIProject_Sample.DoesNotExist:
        with transaction.atomic():
            project_sample = PIProject_Sample()
            project_sample.project = project
            project_sample.sample = sample
            project_sample.name = sample.name
            project_sample_input = sample.file_name_1
            if sample.is_valid_2:
                project_sample_input += ";" + sample.file_name_2
            project_sample.input = project_sample_input
            sample_technology = ConstantsSettings.TECHNOLOGY_minion
            if sample.type_of_fastq == Sample.TYPE_OF_FASTQ_illumina:
                sample_technology = ConstantsSettings.TECHNOLOGY_illumina
            if project.technology != sample_technology:
                print(
                    "Project has different technology {} from sample technology {}...".format(
                        project.technology, sample_technology
                    )
                )
            else:
                project_sample.technology = sample.type_of_fastq
                project_sample.report = "report"
                project_sample.save()
