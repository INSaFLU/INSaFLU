import os
import shutil
import sys
from time import timezone

import pandas as pd
from django.contrib.auth.models import User

from pathogen_identification.constants_settings import ConstantsSettings
from pathogen_identification.install_registry import Deployment_Params
from pathogen_identification.models import ParameterSet, Processed, Submitted
from pathogen_identification.modules.run_main import RunMain_class
from pathogen_identification.utilities.update_DBs import (
    Update_QC_report,
    Update_Sample,
    Update_Sample_Runs,
)


class Run_Main_from_Leaf:
    user: User
    file_r1: str
    file_r2: str
    file_sample: str
    technology: str
    project_name: str
    description: str
    date_created: str
    date_modified: str
    pk: int

    container: PathogenIdentification_deployment

    def __init__(self, user, input_data, project, pipeline_leaf, pipeline_tree):
        self.user = user
        self.sample = input_data
        self.project = project
        self.pipeline_leaf = pipeline_leaf
        self.pipeline_tree = pipeline_tree

        self.file_r1 = input_data.sample.path_name_1
        self.file_r2 = input_data.sample.path_name_2

        self.technology = input_data.sample.technology
        # self.user = user.username
        self.project_name = project.name
        self.date_created = project.creation_date
        self.date_modified = project.last_change_date
        # self.pk = input_data.pk

        self.static_dir = os.path.join(
            ConstantsSettings.static_directory_deployment,
            self.user.username,
            self.project_name,
            simplify_name(os.path.basename(self.file_r1)),
        )

        self.static_dir = os.path.join(
            ConstantsSettings.static_directory, self.static_dir
        )

        self.container = PathogenIdentification_deployment(
            project=self.project_name,
            prefix=f"{simplify_name(input_data.name)}_{input_data.sample.pk}",
            username=self.user.username,
            static_dir=self.static_dir,
        )

        self.pk = self.register_parameter_set()

    def register_parameter_set(self):
        try:
            new_run = ParameterSet.objects.get(
                leaf=self.pipeline_leaf,
                sample=self.sample,
            )

            try:
                submitted = Submitted.objects.get(
                    parameter_set=new_run,
                )

                print("Sample Run already submitted. Please wait for results")
                sys.exit(1)

            except Submitted.DoesNotExist:
                print("Sample Run not submitted. Submitting now.")

            try:
                processed = Processed.objects.get(
                    parameter_set=new_run,
                )
                print("Sample Run already processed. Please wait for results")
                sys.exit(1)
            except Processed.DoesNotExist:
                print("Sample Run not processed. Processing now.")

            return new_run.pk

        except ParameterSet.DoesNotExist:
            new_run = ParameterSet.objects.create(
                leaf=self.pipeline_leaf,
                sample=self.sample,
            )

            return new_run.pk

    def configure(self):
        self.container.configure(
            self.file_r1,
            self.pipeline_leaf,
            r2_path=self.file_r2,
        )

    def Deploy(self):
        self.container.run_main_prep()
        self.container.run_engine.Run()
        self.container.run_engine.Summarize()
        self.container.run_engine.generate_output_data_classes()

    def Update_dbs(self):

        Update_Sample(
            self.container.run_engine.sample,
        )

        Update_QC_report(self.container.run_engine.sample)
        Update_Sample_Runs(self.container.run_engine)

    def register_submission(self):
        try:
            submitted = Submitted.objects.get(pk=self.pk)
        except Submitted.DoesNotExist:
            new_run = ParameterSet.objects.get(pk=self.pk)
            submitted = Submitted(
                parameter_set=new_run,
                date_submitted=timezone.now(),
            )
            submitted.save()

    def register_completion(self):

        try:
            processed = Processed.objects.get(pk=self.pk)
        except Processed.DoesNotExist:
            new_run = ParameterSet.objects.get(pk=self.pk)
            processed = Processed(
                parameter_set=new_run,
                date_processed=timezone.now(),
            )
            processed.save()

    def Submit(self):

        self.register_submission()
        self.configure()
        self.Deploy()
        self.Update_dbs()
        self.container.close()
        self.register_completion()


class PathogenIdentification_deployment:

    project: str
    prefix: str
    rdir: str
    threads: int = 3
    run_engine: RunMain_class
    params = dict
    run_params_db = pd.DataFrame()
    pk: int = 0
    username: str

    def __init__(
        self,
        project: str = "test",
        prefix: str = "main",
        username: str = "admin",
        technology: str = "ONT",
        pk: int = 0,
        static_dir: str = "static",
    ) -> None:

        self.username = username
        self.project = project
        self.prefix = prefix
        self.rdir = os.path.join(ConstantsSettings.project_directory, username, project)
        self.dir = os.path.join(self.rdir, prefix)

        os.makedirs(self.dir, exist_ok=True)

        self.prefix = prefix
        self.pk = pk
        self.static_dir = self.static_dir
        self.technology = technology
        self.install_registry = Deployment_Params

    def configure(self, r1_path: str, pipeline_leaf, r2_path: str = "") -> None:
        self.configure_params(pipeline_leaf)
        self.generate_config_file()
        self.prep_test_env()

        r1_name = os.path.basename(r1_path)
        new_r1_path = os.path.join(self.dir, "reads") + "/" + r1_name
        shutil.copy(r1_path, new_r1_path)

        if r2_path:
            r2_name = os.path.basename(r1_path)
            new_r2_path = os.path.join(self.dir, "reads") + "/" + r2_name
            shutil.copy(r2_path, new_r2_path)

        self.config["sample_name"] = simplify_name(r1_name)
        self.config["r1"] = new_r1_path
        self.config["r2"] = new_r2_path
        self.config["type"] = ["SE", "PE"][
            int(os.path.isfile(self.config["r2"] is not None))
        ]
        self.config["technology"] = self.technology

    def get_constants(self):

        if self.technology == "Illumina/IonTorrent":
            self.constants = ConstantsSettings.CONSTANTS_ILLUMINA
        if self.technology == "ONT":
            self.constants = ConstantsSettings.CONSTANTS_ONT

    def configure_params(self, pipeline_leaf):

        utils = Utils_Manager(owner=self.user)

        all_paths = utils.get_all_technology_pipelines(self.technology)

        leaf_index = self.pipeline_leaf.index

        if leaf_index not in all_paths.keys():
            raise ValueError("Pipeline not found")

        else:
            self.run_params_db = all_paths[leaf_index]

    def generate_config_file(self):

        self.config = {
            "project": self.project,
            "source": self.install_registry.SOURCE,
            "directories": {
                "root": self.rdir,
            },
            "static_dir": self.static_dir,
            "threads": self.threads,
            "prefix": self.prefix,
            "project_name": self.project,
            "metadata": {
                x: os.path.join(self.install_registry.METADATA["ROOT"], g)
                for x, g in self.install_registry.METADATA.items()
            },
            "technology": self.technology,
            "bin": self.install_registry.BINARIES,
            "actions": {},
        }

        for dr, g in ConstantsSettings.DIRS.items():
            self.config["directories"][dr] = os.path.join(self.dir, g)

        for dr, g in ConstantsSettings.ACTIONS.items():
            self.config["actions"][dr] = True
            self.config["actions"]["DEPLETE"] = False
            self.config["actions"]["CLEAN"] = False

        self.config.update(self.params_dict.CONSTANTS)

    def prep_test_env(self):
        """
        from main directory bearing scripts, params.py and main.sh, create metagenome run directory

        :return:
        """
        #

        for dir in self.config["directories"].values():
            os.makedirs(dir, exist_ok=True)

    def close(self):
        if os.path.exists(self.dir):
            shutil.rmtree(self.dir)

    def run_main_prep(self):

        os.makedirs(self.static_dir, exist_ok=True)

        self.run_engine = RunMain_class(self.config, self.run_params_db, self.username)
