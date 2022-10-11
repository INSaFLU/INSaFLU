import datetime
import os
import shutil
import traceback

import pandas as pd
from django.contrib.auth.models import User
from managing_files.models import ProcessControler
from utils.process_SGE import ProcessSGE

from pathogen_identification.constants_settings import ConstantsSettings
from pathogen_identification.install_registry import Deployment_Params
from pathogen_identification.models import (
    ParameterSet,
    PIProject_Sample,
    Processed,
    Projects,
    SoftwareTree,
    SoftwareTreeNode,
    Submitted,
)
from pathogen_identification.modules.run_main import RunMain_class
from pathogen_identification.utilities.update_DBs import (
    Update_QC_report,
    Update_Sample_Runs,
)
from pathogen_identification.utilities.utilities_general import simplify_name
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager


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
        pipeline_index: int,
        sample,  # sample name
        project: str = "test",
        prefix: str = "main",
        username: str = "admin",
        technology: str = "ONT",
        pk: int = 0,
        deployment_root_dir: str = "project",
        dir_branch: str = "deployment",
    ) -> None:

        self.pipeline_index = pipeline_index

        self.username = username
        self.project = project
        self.sample = sample
        self.prefix = prefix
        self.deployment_root_dir = deployment_root_dir
        self.dir_branch = dir_branch
        self.dir = os.path.join(self.deployment_root_dir, dir_branch)

        self.prefix = prefix
        self.pk = pk
        self.technology = technology
        self.install_registry = Deployment_Params

    def configure(self, r1_path: str, r2_path: str = "") -> None:
        self.get_constants()
        self.configure_params()
        self.generate_config_file()
        self.prep_test_env()

        r1_name = os.path.basename(r1_path)
        new_r1_path = os.path.join(self.dir, "reads") + "/" + r1_name
        shutil.copy(r1_path, new_r1_path)

        if r2_path:
            r2_name = os.path.basename(r1_path)
            new_r2_path = os.path.join(self.dir, "reads") + "/" + r2_name
            shutil.copy(r2_path, new_r2_path)
        else:
            new_r2_path = ""

        self.config["sample_name"] = self.sample
        self.config["r1"] = new_r1_path
        self.config["r2"] = new_r2_path
        self.config["type"] = ["SE", "PE"][int(os.path.isfile(self.config["r2"]))]

    def get_constants(self):

        if self.technology == "Illumina/IonTorrent":
            self.constants = ConstantsSettings.CONSTANTS_ILLUMINA
        if self.technology == "ONT":
            self.constants = ConstantsSettings.CONSTANTS_ONT

    def configure_params(self):
        user = User.objects.get(username=self.username)

        utils = Utils_Manager(user)

        all_paths = utils.get_all_technology_pipelines(self.technology)

        leaf_index = self.pipeline_index

        if leaf_index not in all_paths.keys():
            raise ValueError("Pipeline not found")

        else:
            self.run_params_db = all_paths[leaf_index]

    def generate_config_file(self):

        self.config = {
            "project": self.project,
            "source": self.install_registry.SOURCE,
            "technology": self.technology,
            "deployment_root_dir": self.deployment_root_dir,
            "sub_directory": self.dir_branch,
            "directories": {},
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
            self.config["actions"][dr] = g

        self.config.update(self.constants)

        print(self.config["directories"])

    def prep_test_env(self):
        """
        from main directory bearing scripts, params.py and main.sh, create metagenome run directory

        :return:
        """
        os.makedirs(self.dir, exist_ok=True)
        os.makedirs(
            os.path.join(ConstantsSettings.media_directory, self.dir_branch),
            exist_ok=True,
        )
        os.makedirs(
            os.path.join(ConstantsSettings.static_directory, self.dir_branch),
            exist_ok=True,
        )

        for directory in self.config["directories"].values():
            os.makedirs(directory, exist_ok=True)

    def close(self):
        if os.path.exists(self.dir):
            shutil.rmtree(self.dir)

    def run_main_prep(self):

        self.run_engine = RunMain_class(self.config, self.run_params_db, self.username)


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
    date_submitted = datetime.datetime

    container: PathogenIdentification_deployment

    def __init__(
        self,
        user: User,
        input_data: PIProject_Sample,
        project: Projects,
        pipeline_leaf: SoftwareTreeNode,
        pipeline_tree: SoftwareTree,
        odir: str,
    ):
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        self.user = user
        self.sample = input_data
        self.project = project
        self.pipeline_leaf = pipeline_leaf
        self.pipeline_tree = pipeline_tree
        prefix = f"{simplify_name(input_data.name)}_{user.pk}_{project.pk}_{pipeline_leaf.pk}"
        self.date_submitted = datetime.datetime.now()

        self.file_r1 = input_data.sample.path_name_1.path
        if input_data.sample.path_name_2:
            self.file_r2 = input_data.sample.path_name_2.path
        else:
            self.file_r2 = ""

        self.technology = input_data.sample.get_type_technology()
        # self.user = user.username
        self.project_name = project.name
        self.date_created = project.creation_date
        self.date_modified = project.last_change_date
        # self.pk = input_data.pk

        self.deployment_directory_structure = os.path.join(
            ConstantsSettings.televir_subdirectory,
            self.user.username,
            self.project_name,
            simplify_name(os.path.basename(input_data.name)),
            prefix,
        )

        self.unique_id = prefix

        self.deployment_directory = os.path.join(
            odir, self.deployment_directory_structure
        )

        self.container = PathogenIdentification_deployment(
            pipeline_index=pipeline_leaf.index,
            sample=input_data.name,
            project=self.project_name,
            prefix=prefix,
            username=self.user.username,
            deployment_root_dir=odir,
            dir_branch=self.deployment_directory_structure,
        )

        self.parameter_set = self.register_parameter_set()

        self.pk = self.parameter_set.pk

    def check_submission(self, parameter_set):
        try:
            submitted = Submitted.objects.get(
                parameter_set=parameter_set,
            )

            print("Sample Run already submitted. Please wait for results")
            True

        except Submitted.DoesNotExist:
            print("Sample Run not submitted. Submitting now.")
            return False

    def check_processed(self, parameter_set):
        try:
            processed = Processed.objects.get(
                parameter_set=parameter_set,
            )

            print("Sample Run already processed. Please wait for results")
            return True

        except Processed.DoesNotExist:
            print("Sample Run not processed. Processing now.")
            False

    def register_parameter_set(self):

        try:
            new_run = ParameterSet.objects.get(
                leaf=self.pipeline_leaf,
                sample=self.sample,
                project=self.project,
            )

            return new_run

        except ParameterSet.DoesNotExist:
            new_run = ParameterSet.objects.create(
                leaf=self.pipeline_leaf,
                sample=self.sample,
                project=self.project,
            )

            return new_run

    def configure(self):
        self.container.configure(
            self.file_r1,
            r2_path=self.file_r2,
        )

    def Deploy(self):

        try:
            self.container.run_main_prep()
            self.container.run_engine.Run()
            self.container.run_engine.export_sequences()
            self.container.run_engine.Summarize()
            self.container.run_engine.generate_output_data_classes()
            return True
        except Exception as e:
            print(traceback.format_exc())
            print(e)
            return False

    def Update_dbs(self):

        Update_QC_report(self.container.run_engine.sample, self.parameter_set)
        Update_Sample_Runs(self.container.run_engine, self.parameter_set)

    def register_submission(self):
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        process_SGE.set_process_controler(
            self.user,
            process_controler.get_name_televir_project(self.unique_id),
            ProcessControler.FLAG_RUNNING,
        )

        new_run = ParameterSet.objects.get(pk=self.pk)
        new_run.register_subprocess()

        try:
            submitted = Submitted.objects.get(
                parameter_set=new_run,
            )
        except Submitted.DoesNotExist:

            submitted = Submitted(
                parameter_set=new_run,
                date_submitted=self.date_submitted,
            )
            submitted.save()

    def delete_submission_error(self):
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        process_SGE.set_process_controler(
            self.user,
            process_controler.get_name_televir_project(self.unique_id),
            ProcessControler.FLAG_ERROR,
        )

        new_run = ParameterSet.objects.get(pk=self.pk)
        new_run.register_error()

        try:
            submitted = Submitted.objects.get(
                parameter_set=new_run,
            )
            submitted.delete()
        except Submitted.DoesNotExist:
            pass

    def delete_submission(self):
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        process_SGE.set_process_controler(
            self.user,
            process_controler.get_name_televir_project(self.unique_id),
            ProcessControler.FLAG_FINISHED,
        )

        new_run = ParameterSet.objects.get(pk=self.pk)
        new_run.register_finished()

        try:
            submitted = Submitted.objects.get(
                parameter_set=new_run,
            )
            submitted.delete()
        except Submitted.DoesNotExist:
            pass

    def register_completion(self):

        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        process_SGE.set_process_controler(
            self.user,
            process_controler.get_name_televir_project(self.unique_id),
            ProcessControler.FLAG_FINISHED,
        )

        new_run = ParameterSet.objects.get(pk=self.pk)
        date_processed = datetime.datetime.now()

        try:
            processed = Processed.objects.get(
                parameter_set=new_run,
            )
        except Processed.DoesNotExist:

            processed = Processed(
                parameter_set=new_run,
                date_processed=date_processed,
            )
            processed.save()

    def update_project_change_date(self):
        self.project.last_change_date = datetime.datetime.now()
        self.project.save()

    def Submit(self):

        if not self.check_submission(self.parameter_set) and not self.check_processed(
            self.parameter_set
        ):
            self.register_submission()
            self.configure()
            run_success = self.Deploy()
            if run_success:
                self.Update_dbs()
                self.container.close()
                self.register_completion()
                self.delete_submission()
                self.update_project_change_date()

            else:
                print("Error in run")
                self.delete_submission_error()
