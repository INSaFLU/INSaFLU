import datetime
import os
import shutil
import traceback
from typing import List, Optional

import pandas as pd
from django.contrib.auth.models import User

from constants.constants import Televir_Metadata_Constants as Televir_Metadata
from constants.constants import TypePath
from managing_files.models import ProcessControler
from pathogen_identification.constants_settings import ConstantsSettings
from pathogen_identification.models import (
    FinalReport,
    ParameterSet,
    PIProject_Sample,
    Projects,
    RunMain,
    SoftwareTree,
    SoftwareTreeNode,
)
from pathogen_identification.modules.object_classes import Remap_Target
from pathogen_identification.modules.run_main import RunMainTree_class
from pathogen_identification.utilities.televir_parameters import TelevirParameters
from pathogen_identification.utilities.update_DBs import (
    Update_Assembly,
    Update_Classification,
    Update_Metagenomics,
    Update_Remap,
    Update_RunMain_Initial,
    Update_RunMain_Secondary,
    get_run_parents,
)
from pathogen_identification.utilities.utilities_general import simplify_name_lower
from pathogen_identification.utilities.utilities_pipeline import (
    SoftwareTreeUtils,
    Utils_Manager,
)
from pathogen_identification.utilities.utilities_views import ReportSorter
from settings.constants_settings import ConstantsSettings as SettingsConstants
from utils.process_SGE import ProcessSGE


class PathogenIdentificationDeploymentCore:

    run_engine: RunMainTree_class

    def __init__(
        self,
        sample: PIProject_Sample,
        deployment_root_dir: str = "/tmp/insaflu/insaflu_something",
        dir_branch: str = "deployment",
        prefix: str = "main",
        threads: int = 3,
    ):
        self.sample = sample
        self.prefix = prefix

        self.deployment_root_dir = deployment_root_dir
        self.dir_branch = dir_branch
        self.dir = os.path.join(self.deployment_root_dir, dir_branch)

        self.install_registry = Televir_Metadata()

        self.threads = threads
        self.prepped = False
        self.run_params_db = pd.DataFrame()
        self.config = dict()

        self.file_r1 = sample.sample.get_fastq_available(TypePath.MEDIA_ROOT, True)
        if sample.sample.exist_file_2():
            self.file_r2 = sample.sample.get_fastq_available(TypePath.MEDIA_ROOT, False)
        else:
            self.file_r2 = ""

    @property
    def username(self):
        return self.sample.project.owner.username

    @property
    def user(self):
        return self.sample.project.owner

    @property
    def technology(self):
        return self.sample.project.technology

    @property
    def project_name(self):
        return self.sample.project.name

    @property
    def project_pk(self):
        return self.sample.project.pk

    @property
    def project(self):
        return self.sample.project

    def input_read_project_path(self, filepath) -> str:
        """copy input reads to project directory and return new path"""

        if not os.path.isfile(filepath):
            return ""

        rname = os.path.basename(filepath)
        new_rpath = os.path.join(self.dir, "reads") + "/" + rname
        shutil.copy(filepath, new_rpath)
        return new_rpath

    def delete_run_media(self):
        """delete project media directory"""

        if self.prepped:
            try:
                if os.path.isdir(self.run_engine.media_dir):
                    shutil.rmtree(self.run_engine.media_dir, ignore_errors=True)

            except AttributeError as e:
                pass

    def delete_run_static(self):
        """delete project static directory"""

        if self.prepped:

            try:
                if os.path.isdir(self.run_engine.static_dir):
                    shutil.rmtree(self.run_engine.static_dir, ignore_errors=True)

            except AttributeError as e:
                pass

    def get_constants(self):
        """set constants for technology"""
        if self.technology == SettingsConstants.TECHNOLOGY_illumina:
            self.constants = ConstantsSettings.CONSTANTS_ILLUMINA
        if self.technology == SettingsConstants.TECHNOLOGY_minion:
            self.constants = ConstantsSettings.CONSTANTS_ONT

    def generate_config_file(self):
        self.config = {
            "project": self.project_name,
            "source": self.install_registry.SOURCE,
            "technology": self.technology,
            "deployment_root_dir": self.deployment_root_dir,
            "sub_directory": self.dir_branch,
            "directories": {},
            "threads": self.threads,
            "prefix": self.prefix,
            "project_name": self.project_name,
            "metadata": self.install_registry.metadata_full_path,
            "bin": self.install_registry.BINARIES,
            "actions": {},
        }

        for dr, g in ConstantsSettings.DIRS.items():
            self.config["directories"][dr] = os.path.join(self.dir, g)

        for dr, g in ConstantsSettings.ACTIONS.items():
            self.config["actions"][dr] = g

        self.config.update(self.constants)

    def update_config_prefix(self):
        self.config["prefix"] = self.prefix

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

    def import_params(self, run_params_db: pd.DataFrame):

        self.run_params_db = run_params_db

    def configure_constants(self) -> bool:
        """generate config dictionary for run_main, and copy input reads to project directory."""

        try:
            self.get_constants()

            self.generate_config_file()
            self.prep_test_env()

            new_r1_path = self.input_read_project_path(self.file_r1)
            new_r2_path = self.input_read_project_path(self.file_r2)

            self.config["sample_registered"] = self.sample
            self.config["sample_name"] = self.sample.name
            self.config["r1"] = new_r1_path
            self.config["r2"] = new_r2_path

            self.config["type"] = [
                ConstantsSettings.SINGLE_END,
                ConstantsSettings.PAIR_END,
            ][int(os.path.isfile(self.config["r2"]))]

            return True

        except Exception as e:
            print(e)
            return False

    def run_main_prep_check_first(self):
        """
        check if run_main has been prepped and not root, if so, prep it
        by equiping it with the run_engine"""
        if self.prepped or self.run_params_db.empty:
            return

        self.run_engine = RunMainTree_class(
            self.config, self.run_params_db, self.project_pk
        )
        self.run_engine.Prep_deploy()
        self.run_engine.generate_output_data_classes()
        self.prepped = True

    def run_main(self):
        self.run_engine.Run_QC()
        self.run_engine.Run_PreProcess()
        self.run_engine.Sanitize_reads()
        self.run_engine.Run_Assembly()
        self.run_engine.Run_Contig_classification()
        self.run_engine.Run_Read_classification()
        self.run_engine.Run_Remapping()

    def update_engine(self):

        if "module" in self.run_params_db.columns:
            self.run_engine.Update(self.config, self.run_params_db)
            self.run_engine.Prep_deploy(remap_prep=False)

    def update_merged_targets(self, merged_targets: List[Remap_Target]):
        self.run_engine.update_merged_targets(merged_targets)

    def delete_run_record(self, parameter_set: ParameterSet):
        """delete project record in database"""

        if self.prepped:
            _, runmain, _ = get_run_parents(self.run_engine, parameter_set)

            if runmain is not None:
                runmain.delete()

    def delete_run(self):
        """delete project record in database"""

        self.delete_run_media()
        self.delete_run_static()


class PathogenIdentification_SingleDeployment(PathogenIdentificationDeploymentCore):
    project_name: str
    prefix: str
    threads: int
    run_engine: RunMainTree_class
    params = dict
    run_params_db = pd.DataFrame()
    pk: int
    username: str
    prepped: bool

    def __init__(
        self,
        pipeline_index: int,
        sample: PIProject_Sample,  # sample name
        prefix: str = "main",
        pk: int = 0,
        deployment_root_dir: str = "/tmp/insaflu/insaflu_something",
        dir_branch: str = "deployment",
        threads: int = 3,
    ) -> None:

        super().__init__(sample, deployment_root_dir, dir_branch, prefix, threads)

        self.pipeline_index = pipeline_index
        self.pk = pk
        self.parameter_set = ParameterSet.objects.get(pk=pk)
        self.tree_makup = self.parameter_set.leaf.software_tree.global_index

    def configure_params(self):
        """get pipeline parameters from database"""

        software_tree_utils = SoftwareTreeUtils(self.project.owner, self.project)

        all_paths = software_tree_utils.get_all_technology_pipelines(self.tree_makup)

        self.run_params_db = all_paths.get(self.pipeline_index, None)

        if self.run_params_db is None:
            print("Pipeline index not found")
            return False

        return True

    def configure_w_branch(self) -> bool:
        """generate config dictionary for run_main, and copy input reads to project directory."""

        branch_exists = self.configure_params()

        if not branch_exists:
            return False

        success = self.configure_constants()

        return success

    def run_main_prep_dump_tables(self):
        """prepare run_main object from config dictionary"""
        self.run_main_prep_check_first()

        utils = Utils_Manager()
        utils.dump_tables(self.run_engine.log_dir)

        self.prepped = True


class Run_Main_from_Leaf:
    user: User

    file_sample: str
    technology: str
    project_name: str
    description: str
    date_created: str
    date_modified: str
    pk: int
    date_submitted = datetime.datetime
    combined_analysis: bool
    mapping_request: bool
    container: PathogenIdentification_SingleDeployment

    def __init__(
        self,
        user: User,
        input_data: PIProject_Sample,
        project: Projects,
        pipeline_leaf: SoftwareTreeNode,
        pipeline_tree: SoftwareTree,
        odir: str,
        threads: int = 3,
        combined_analysis: bool = False,
        mapping_request: bool = False,
        run_pk: Optional[int] = None,
    ):
        self.user = user
        self.sample = input_data
        self.project = project
        self.combined_analysis = combined_analysis
        self.mapping_request = mapping_request
        self.run_pk = run_pk
        self.pipeline_leaf = pipeline_leaf
        self.pipeline_tree = pipeline_tree
        ########################################
        prefix = f"{simplify_name_lower(input_data.name)}_run{pipeline_leaf.index}"
        self.date_submitted = datetime.datetime.now()

        self.technology = input_data.sample.get_type_technology()
        self.project_name = project.name
        self.date_created = project.creation_date
        self.date_modified = project.last_change_date

        self.deployment_directory_structure = os.path.join(
            ConstantsSettings.televir_subdirectory,
            f"{self.user.pk}",
            f"{project.pk}",
            f"{input_data.sample.pk}",
            prefix,
        )

        self.unique_id = prefix

        self.deployment_directory = os.path.join(
            odir, self.deployment_directory_structure
        )

        self.parameter_set = self.register_parameter_set()
        self.pk = self.parameter_set.pk

        self.container = PathogenIdentification_SingleDeployment(
            pipeline_index=pipeline_leaf.index,
            sample=input_data,
            prefix=prefix,
            deployment_root_dir=odir,
            dir_branch=self.deployment_directory_structure,
            pk=self.pk,
            threads=threads,
        )

        self.is_available = self.check_availability()

    def get_status(self):
        return self.parameter_set.status

    def check_finished(self):
        return self.parameter_set.status == ParameterSet.STATUS_FINISHED

    def check_availability(self):
        return self.parameter_set.status not in [
            ParameterSet.STATUS_RUNNING,
            ParameterSet.STATUS_FINISHED,
        ]

    def set_to_queued(self):
        if self.is_available:
            self.parameter_set.status = ParameterSet.STATUS_QUEUED
            self.parameter_set.save()

    def check_submission(self):
        if self.parameter_set.status in [
            ParameterSet.STATUS_RUNNING,
        ]:
            return True

        else:
            return False

    def check_processed(self):
        if self.parameter_set.status in [
            ParameterSet.STATUS_FINISHED,
        ]:
            return True

        else:
            return False

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
        configured = self.container.configure_w_branch()

        if not configured:
            return False

        self.container.run_main_prep_dump_tables()

        try:

            if self.run_pk is not None:
                self.container.run_engine.run_pk = self.run_pk

            if (
                self.container.run_engine.run_type
                == RunMainTree_class.RUN_TYPE_SCREENING
            ):
                self.container.run_engine.remap_params.manual_references_include = True

            if self.mapping_request:
                self.container.run_engine.run_type = (
                    RunMainTree_class.RUN_TYPE_MAPPING_REQUEST
                )

                self.container.run_engine.metadata_tool.get_mapping_references(
                    self.run_pk,
                    max_accids=self.container.run_engine.remap_params.max_accids,
                )

            if self.combined_analysis:
                self.container.run_engine.run_type = (
                    RunMainTree_class.RUN_TYPE_COMBINED_MAPPING
                )

                if (
                    self.container.run_engine.remap_params.manual_references_include
                    is True
                ):
                    self.container.run_engine.metadata_tool.get_manual_references(
                        self.sample,
                        max_accids=self.container.run_engine.remap_params.max_accids,
                    )

        except Exception as e:
            print(e)
            print(traceback.format_exc())
            return False

        return True

    def set_run_process_running(self):
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        process_SGE.set_process_controler(
            self.user,
            process_controler.get_name_televir_run(
                self.project.pk,
                self.sample.pk,
                self.pipeline_leaf.pk,
            ),
            ProcessControler.FLAG_RUNNING,
        )

    def set_run_process_error(self):
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        process_SGE.set_process_controler(
            self.user,
            process_controler.get_name_televir_run(
                self.project.pk,
                self.sample.pk,
                self.pipeline_leaf.pk,
            ),
            ProcessControler.FLAG_ERROR,
        )

    def set_run_process_finished(self):
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        process_SGE.set_process_controler(
            self.user,
            process_controler.get_name_televir_run(
                self.project.pk,
                self.sample.pk,
                self.pipeline_leaf.pk,
            ),
            ProcessControler.FLAG_FINISHED,
        )

    def Deploy_Parts(self):
        try:

            self.container.run_engine.Prep_deploy()
            self.container.run_engine.Run_QC()
            db_updated = Update_RunMain_Initial(
                self.container.run_engine, self.parameter_set
            )
            self.register_running()
            if not db_updated:
                return False
        except Exception as e:
            print(traceback.format_exc())
            print(e)
            return False

        try:
            self.container.run_engine.Run_PreProcess()
            db_updated = Update_RunMain_Secondary(
                self.container.run_engine, self.parameter_set
            )
            if not db_updated:
                return False
        except Exception as e:
            print(traceback.format_exc())
            print(e)
            return False

        try:
            self.container.run_engine.Run_Assembly()
            self.container.run_engine.export_assembly()
            db_updated = Update_Assembly(self.container.run_engine, self.parameter_set)
            if not db_updated:
                return False
        except Exception as e:
            print(traceback.format_exc())
            print(e)
            return False

        try:
            self.container.run_engine.Run_Classification()
            self.container.run_engine.plan_remap_prep_safe()
            self.container.run_engine.export_intermediate_reports()
            self.container.run_engine.generate_output_data_classes()
            db_updated = Update_Classification(
                self.container.run_engine, self.parameter_set
            )
            if not db_updated:
                return False
        except Exception as e:
            print(traceback.format_exc())
            print(e)
            return False

        if self.combined_analysis:
            self.container.run_engine.plan_combined_remapping()

        try:
            self.container.run_engine.Run_Metagenomics_Classification()

            db_updated = Update_Metagenomics(
                self.container.run_engine, self.parameter_set
            )
            if not db_updated:
                return False

            if (
                self.container.run_engine.read_metagenomics_classification_performed
                is True
            ):
                return True

        except Exception as e:
            print(traceback.format_exc())
            print(e)
            return False

        try:
            #
            self.container.run_engine.Run_Remapping()
            self.container.run_engine.export_sequences()
            self.container.run_engine.export_intermediate_reports()
            self.container.run_engine.Summarize()
            self.container.run_engine.generate_output_data_classes()
            self.container.run_engine.export_logdir()

            db_updated = Update_Remap(self.container.run_engine, self.parameter_set)

            if not db_updated:
                return False

        except Exception as e:
            print(traceback.format_exc())
            print(e)
            return False

        return True

    def register_submission(self):
        self.set_run_process_running()
        self.parameter_set.register_subprocess()
        new_run = ParameterSet.objects.get(pk=self.pk)
        new_run.register_subprocess()
        print("registered_submission")

        if self.run_pk:
            run = RunMain.objects.get(pk=self.run_pk)
            run.status = RunMain.STATUS_PREP
            run.save()

    def register_running(self):
        if self.run_pk:
            run = RunMain.objects.get(pk=self.run_pk)
            run.status = RunMain.STATUS_RUNNING
            run.save()

    def register_error(self):
        self.set_run_process_error()
        print("REGISTERING ERROR")
        print("RUN PS PK", self.pk)

        new_run = ParameterSet.objects.get(pk=self.pk)
        new_run.register_error()

        try:
            if self.run_pk:
                run = RunMain.objects.get(pk=self.run_pk)
                run.status = RunMain.STATUS_ERROR
                run.save()
            else:
                run = RunMain.objects.get(parameter_set=new_run)
                run.delete()
        except RunMain.DoesNotExist:
            pass

        self.container.delete_run()

    def run_reference_overlap_analysis(self):
        final_reports = FinalReport.objects.filter(sample=self.sample).order_by(
            "-coverage"
        )
        #
        report_layout_params = TelevirParameters.get_report_layout_params(
            project_pk=self.parameter_set.project.pk
        )

        if final_reports.exists() is False:
            return

        report_sorter = ReportSorter(self.sample, final_reports, report_layout_params)

        try:
            report_sorter.sort_reports_save()
        except Exception as e:
            print(e)
            print(traceback.format_exc())

            print("Error in report sorter")
            return

    def register_completion(self):
        self.set_run_process_finished()
        new_run = ParameterSet.objects.get(pk=self.pk)
        new_run.register_finished()

        if self.run_pk:
            run = RunMain.objects.get(pk=self.run_pk)
            run.status = RunMain.STATUS_FINISHED
            run.save()

    def update_project_change_date(self):
        self.project.last_change_date = datetime.datetime.now()
        self.project.save()

    def Submit(self):
        if not self.check_submission() and not self.check_processed():
            self.register_submission()

            try:
                configured = self.configure()
            except Exception as e:
                print(e)
                self.register_error()
                return

            if configured:
                run_success = self.Deploy_Parts()
            else:
                print("Error in configuration")
                self.register_error()
                return

            if run_success:
                self.register_completion()
                self.update_project_change_date()

            else:
                print("Error in run")
                self.register_error()
                return
