import os
from random import sample
from typing import Dict, Tuple

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase, tag

from constants.constants import Constants
from constants.constants import Televir_Metadata_Constants as Deployment_Params
from constants.constantsTestsCase import ConstantsTestsCase
from constants.software_names import SoftwareNames
from fluwebvirus.settings import STATIC_ROOT
from pathogen_identification.constants_settings import \
    ConstantsSettings as PI_CS
from pathogen_identification.deployment_main import Run_Main_from_Leaf
from pathogen_identification.models import (ParameterSet, PIProject_Sample,
                                            Projects, SoftwareTree,
                                            SoftwareTreeNode)
from pathogen_identification.modules.object_classes import (
    Operation_Temp_Files, Read_class, RunCMD, Temp_File)
from pathogen_identification.utilities.utilities_pipeline import (
    Pipeline_Makeup, PipelineTree, SoftwareTreeUtils, Utils_Manager)
from settings.constants_settings import ConstantsSettings as CS
from settings.default_software import DefaultSoftware
from settings.models import Parameter, Sample, Software
from utils.software import Software as SoftwareUtils
from utils.utils import Utils

# Create your tests here.


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


def get_bindir_from_binaries(binaries, key, value: str = ""):
    if value == "":
        try:
            return os.path.join(binaries["ROOT"], binaries[key]["default"], "bin")
        except KeyError:
            return ""
    else:
        try:
            return os.path.join(binaries["ROOT"], binaries[key][value], "bin")
        except KeyError:
            return ""


def update_software_params_global_project(project, user):
    """
    update software global to project
    """
    ### get all global software
    query_set = Software.objects.filter(
        owner=user,
        type_of_use__in=[
            Software.TYPE_OF_USE_televir_global,
            Software.TYPE_OF_USE_televir_settings,
        ],
        type_of_software__in=[
            Software.TYPE_SOFTWARE,
            Software.TYPE_INSAFLU_PARAMETER,
        ],
        is_obsolete=False,
        technology__name=project.technology,
    )
    project = Projects.objects.get(pk=project.pk)
    for software in query_set:
        software_parameters = Parameter.objects.filter(
            software=software,
        )

        software.pk = None
        software.type_of_use = Software.TYPE_OF_USE_televir_project

        try:
            Software.objects.get(
                name=software.name,
                type_of_use=Software.TYPE_OF_USE_televir_project,
                parameter__televir_project=project,
                pipeline_step=software.pipeline_step,
            )

        except Software.MultipleObjectsReturned:
            pass

        except Software.DoesNotExist:
            software.save()
            for parameter in software_parameters:
                parameter.pk = None
                parameter.software = software
                parameter.televir_project = project
                parameter.save()


def duplicate_software_params_global_project(user, project):
    """
    duplicate software global to project
    """
    ### get all global software
    query_set = Software.objects.filter(
        owner=user,
        type_of_use__in=[
            Software.TYPE_OF_USE_televir_global,
            Software.TYPE_OF_USE_televir_settings,
        ],
        type_of_software__in=[
            Software.TYPE_SOFTWARE,
            Software.TYPE_INSAFLU_PARAMETER,
        ],
        is_obsolete=False,
        technology__name=project.technology,
        parameter__televir_project=None,
    )

    for software in query_set:
        software_parameters = Parameter.objects.filter(
            software=software,
        )

        software.pk = None
        if software.type_of_use == Software.TYPE_OF_USE_televir_global:
            software.type_of_use = Software.TYPE_OF_USE_televir_project
        else:
            software.type_of_use = Software.TYPE_OF_USE_televir_project_settings

        try:
            Software.objects.get(
                name=software.name,
                type_of_use=software.type_of_use,
                parameter__televir_project=project,
                pipeline_step=software.pipeline_step,
            )

        except Software.MultipleObjectsReturned:
            pass

        except Software.DoesNotExist:
            software.save()

            for parameter in software_parameters:
                parameter.pk = None
                parameter.software = software
                parameter.televir_project = project
                parameter.save()


def check_project_params_exist(project):
    """
    check if project parameters exist
    """

    query_set = Parameter.objects.filter(televir_project=project.pk)
    if query_set.count() == 0:
        return False
    return True


def create_project_params(user, televir_project):
    if not check_project_params_exist(televir_project):
        duplicate_software_params_global_project(user, televir_project)
    else:
        update_software_params_global_project(user, televir_project)


def reset_project_params(televir_project):
    if not check_project_params_exist(televir_project):
        return

    else:
        Parameter.objects.filter(televir_project=televir_project).delete()


def set_project_makeup(project: Projects, makeup: list):
    project_software = Software.objects.filter(
        type_of_use__in=[
            Software.TYPE_OF_USE_televir_project,
            Software.TYPE_OF_USE_televir_project_settings,
        ],
        parameter__televir_project=project,
        owner=project.owner,
    )

    for software in project_software:
        if software.pipeline_step.name not in makeup:
            software.is_to_run = False
            software.save()


def reset_project_makeup(project: Projects):
    project_software = Software.objects.filter(
        type_of_use=Software.TYPE_OF_USE_televir_project,
        parameter__televir_project=project,
        owner=project.owner,
    )

    for software in project_software:
        software.is_to_run = True
        software.save()


def get_global_tree(local_tree: PipelineTree, project: Projects):
    user = project.owner
    software_tree_manager = SoftwareTreeUtils(user, project)
    # local_tree = software_tree_manager.generate_project_tree()
    tree_makeup = local_tree.makeup
    pipeline_tree = software_tree_manager.generate_software_tree_extend(local_tree)
    pipeline_tree_index = software_tree_manager.get_software_tree_index(tree_makeup)
    pipeline_tree_query = SoftwareTree.objects.get(pk=pipeline_tree_index)

    return pipeline_tree, pipeline_tree_query


def determine_available_paths(project: Projects, user: User):
    utils_manager = Utils_Manager()
    softwaretree_manager = SoftwareTreeUtils(user, project)

    local_tree = softwaretree_manager.generate_project_tree()
    local_paths = local_tree.get_all_graph_paths_explicit()

    pipeline_tree, pipeline_tree_query = get_global_tree(local_tree, project)

    matched_paths = {
        leaf: utils_manager.utility_manager.match_path_to_tree_safe(path, pipeline_tree)
        for leaf, path in local_paths.items()
    }
    available_paths = {
        leaf: path for leaf, path in matched_paths.items() if path is not None
    }

    return pipeline_tree_query, available_paths


class Televir_Software_Test(TestCase):
    software = SoftwareUtils()
    utils = Utils()
    constants = Constants()
    software_names = SoftwareNames()
    constants_tests_case = ConstantsTestsCase()
    pipeline_makeup = Pipeline_Makeup()

    def setUp(self):
        self.baseDirectory = os.path.join(
            getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS
        )

        try:
            user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
        except User.DoesNotExist:
            user = User()
            user.username = ConstantsTestsCase.TEST_USER_NAME
            user.is_active = False
            user.password = ConstantsTestsCase.TEST_USER_NAME
            user.save()

        self.test_user = user

    def test_software_created(self):
        default_software = DefaultSoftware()
        default_software.test_all_defaults_pathogen_identification(self.test_user)

        blastn = Software.objects.filter(name="Blastn", owner=self.test_user)
        self.assertEqual(blastn.count(), 2)

        centrifuge = Software.objects.filter(name="Centrifuge", owner=self.test_user)

        self.assertEqual(centrifuge.count(), 4)

        kraken2 = Software.objects.filter(name="Kraken2", owner=self.test_user)
        self.assertEqual(kraken2.count(), 1)

        krakenuniq = Software.objects.filter(name="Krakenuniq", owner=self.test_user)
        self.assertEqual(krakenuniq.count(), 2)

        fastviromeexplorer = Software.objects.filter(
            name="FastViromeExplorer", owner=self.test_user
        )
        self.assertEqual(fastviromeexplorer.count(), 1)

        kaiju = Software.objects.filter(name="Kaiju", owner=self.test_user)
        self.assertEqual(kaiju.count(), 1)

        raven = Software.objects.filter(name="Raven", owner=self.test_user)
        self.assertEqual(raven.count(), 1)

        minimap2 = Software.objects.filter(
            name_extended="Minimap2", owner=self.test_user
        )
        self.assertEqual(minimap2.count(), 2)

        snippy = Software.objects.filter(name="Snippy_PI", owner=self.test_user)
        self.assertEqual(snippy.count(), 1)

        bamutil = Software.objects.filter(name="BamUtil", owner=self.test_user)
        self.assertEqual(bamutil.count(), 2)

        prinseqplusplus = Software.objects.filter(
            name="Prinseq++", owner=self.test_user
        )
        self.assertEqual(prinseqplusplus.count(), 2)

    def test_pipeline_makeup(self):
        """THIS TEST IS NOT COMPLETE"""
        default_software = DefaultSoftware()
        default_software.test_all_defaults_pathogen_identification(self.test_user)
        self.assertTrue(len(self.pipeline_makeup.MAKEUP), 16)

        self.assertEqual(
            self.pipeline_makeup.match_makeup_name_from_list(
                ["Contig classification", "Remapping"]
            ),
            None,
        )

        raven = Software.objects.get(
            name="Raven", type_of_use=5
        )  # type_of_use=5 is for pathogen identification

        pipeline_excluding_raven = (
            self.pipeline_makeup.get_software_pipeline_list_excluding(raven)
        )

        self.assertEqual(
            {
                CS.PIPELINE_NAME_viral_enrichment,
                CS.PIPELINE_NAME_contig_classification,
                CS.PIPELINE_NAME_read_classification,
                CS.PIPELINE_NAME_remapping,
                CS.PIPELINE_NAME_reporting,
            },
            set(pipeline_excluding_raven),
        )


def televir_test_project(user: User, project_ont_name: str = "project_televir"):
    ######### ONT ##########

    try:
        project_ont = Projects.objects.get(name=project_ont_name)
    except Projects.DoesNotExist:
        project_ont = Projects()
        project_ont.name = project_ont_name
        project_ont.owner = user
        project_ont.technology = CS.TECHNOLOGY_minion

        project_ont.save()

    return project_ont


def televir_test_sample(project_ont, sample_ont: Sample):
    try:
        ont_project_sample = PIProject_Sample.objects.get(
            project__id=project_ont.pk, sample__id=sample_ont.pk
        )
    except PIProject_Sample.DoesNotExist:
        ont_project_sample = PIProject_Sample()
        ont_project_sample.project = project_ont
        ont_project_sample.sample = sample_ont
        ont_project_sample.name = sample_ont.name
        ont_project_sample.save()

    return ont_project_sample


def test_user():
    try:
        user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
    except User.DoesNotExist:
        user = User()
        user.username = ConstantsTestsCase.TEST_USER_NAME
        user.is_active = False
        user.password = ConstantsTestsCase.TEST_USER_NAME
        user.save()

    return user


def test_fastq_file(
    baseDirectory, user: User, sample_name: str = "televir_sample_minion_1"
):
    utils = Utils()

    file_name = os.path.join(
        STATIC_ROOT,
        "tests",
        ConstantsTestsCase.DIR_FASTQ,
        ConstantsTestsCase.FASTQ_MINION_1,
    )

    utils.copy_file(file_name, os.path.join(baseDirectory, ConstantsTestsCase.FASTQ1_1))

    try:
        sample_ont = Sample.objects.get(name=sample_name)
    except Sample.DoesNotExist:
        sample_ont = Sample()
        sample_ont.name = sample_name
        sample_ont.is_valid_1 = True
        sample_ont.file_name_1 = ConstantsTestsCase.FASTQ1_1
        sample_ont.path_name_1.name = os.path.join(
            baseDirectory, ConstantsTestsCase.FASTQ1_1
        )
        sample_ont.is_valid_2 = False
        sample_ont.type_of_fastq = Sample.TYPE_OF_FASTQ_minion
        sample_ont.owner = user
        sample_ont.save()

    return sample_ont


class Televir_Objects_TestCase(TestCase):
    install_registry = Deployment_Params

    def setUp(self):
        self.baseDirectory = os.path.join(
            getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS
        )
        self.temp_directory = os.path.join(self.baseDirectory, "temp_objects_tests")
        os.makedirs(self.temp_directory, exist_ok=True)

        self.test_user = test_user()
        self.sample_ont = test_fastq_file(self.temp_directory, self.test_user)

    def tearDown(self) -> None:
        # os.system("rm -rf " + self.temp_directory)
        pass

    def test_temporary_operations(self):
        tempop = Operation_Temp_Files(self.baseDirectory)

        with tempop as tpf:
            open(tpf.script, "w").close()
            open(tpf.log, "w").close()
            open(tpf.flag, "w").close()

        self.assertFalse(os.path.exists(tpf.script))
        self.assertFalse(os.path.exists(tpf.log))
        self.assertFalse(os.path.exists(tpf.flag))

    def test_temp_file(self):
        tf = Temp_File(self.baseDirectory)

        with tf as tpf:
            self.assertTrue(os.path.exists(tpf))

        self.assertFalse(os.path.exists(tpf))

    def test_runCMD_strings(self):
        runcmd = RunCMD(
            self.baseDirectory,
            logdir=self.baseDirectory,
            prefix="test_runCMD",
        )

        cmd = "hello"
        self.assertEqual(runcmd.bash_cmd_string(cmd), cmd)
        self.assertEqual(
            runcmd.bash_software_cmd_string(cmd), f"{self.baseDirectory}/{cmd}"
        )
        self.assertEqual(
            runcmd.python_cmd_string(cmd), f"python {self.baseDirectory}/{cmd}"
        )
        java_bin = os.path.join(
            self.install_registry.BINARIES["ROOT"],
            self.install_registry.BINARIES["software"]["java"],
            "bin",
            "java",
        )
        self.assertEqual(
            runcmd.java_cmd_string(cmd), f"{java_bin} -cp {self.baseDirectory}/ {cmd}"
        )

    def test_runCMD_deployment(self):
        runcmd = RunCMD(
            self.baseDirectory,
            logdir=self.baseDirectory,
            prefix="test_runCMD",
        )
        tempf = Temp_File(self.baseDirectory)
        tempPython = Temp_File(self.baseDirectory, suffix=".py")
        with tempf as tmp:
            with tempPython as tpf:
                with open(tpf, "w") as f:
                    f.write("print('hello world')")

                cmd = f"{os.path.basename(tpf)} > {tmp}"

                runcmd.run_python(cmd)
                self.assertTrue(os.path.exists(tmp))
                with open(tmp, "r") as f:
                    self.assertEqual(f.read().strip(), "hello world")

        tempJava = Temp_File(self.baseDirectory, suffix=".java")

        with tempf as tmp:
            with tempJava as tpf:
                with open(tpf, "w") as f:
                    ###
                    f.write(
                        'public class Hello { public static void main(String[] args) { System.out.println("Hello World"); } }'
                    )
                cmd = f"{tpf} > {tmp}"
                runcmd.run_java(cmd)
                self.assertTrue(os.path.exists(tmp))
                with open(tmp, "r") as f:
                    self.assertEqual(f.read().strip(), "Hello World")

        with tempf as tmp:
            cmd = f"echo hello world > {tmp}"
            runcmd.run_bash(cmd)
            self.assertTrue(os.path.exists(tmp))
            with open(tmp, "r") as f:
                self.assertEqual(f.read().strip(), "hello world")

        cmd = "echo hello world"
        cmd_return = runcmd.run_bash_return(cmd)
        self.assertEqual(cmd_return.strip(), "hello world")

    def test_read_class(self):
        r1_path = self.sample_ont.path_name_1.name

        clean_dir = os.path.join(self.temp_directory, "clean")
        enriched_dir = os.path.join(self.temp_directory, "enriched")
        depleted_dir = os.path.join(self.temp_directory, "depleted")
        os.makedirs(clean_dir, exist_ok=True)
        os.makedirs(enriched_dir, exist_ok=True)
        os.makedirs(depleted_dir, exist_ok=True)

        preproc_library = os.path.join(
            self.install_registry.BINARIES["ROOT"],
            self.install_registry.BINARIES[CS.PIPELINE_NAME_read_quality_analysis][
                "default"
            ],
            "bin",
        )

        r1 = Read_class(
            r1_path,
            clean_dir,
            enriched_dir,
            depleted_dir,
            bin=preproc_library,
        )
        r2 = Read_class(
            os.path.join(self.temp_directory, "None"),
            clean_dir,
            enriched_dir,
            depleted_dir,
            bin=preproc_library,
        )

        self.assertEqual(r1.exists, True)
        self.assertEqual(r2.exists, False)
        self.assertEqual(r1.current_status, "raw")
        read_names = r1.get_read_names_fastq()
        self.assertEqual(len(read_names), r1.read_number_raw)

        read_names_subset = sample(read_names, 100)
        temp_read_file = Temp_File(self.temp_directory, suffix=".fq.gz")

        with temp_read_file as tpf:
            r1.read_filter_move(read_names_subset, tpf)

            temp_read = Read_class(
                tpf,
                clean_dir,
                enriched_dir,
                depleted_dir,
                bin=preproc_library,
            )

            self.assertEqual(temp_read.read_number_raw, len(read_names_subset))
            subsample = sample(read_names_subset, 10)
            temp_reads_file = Temp_File(self.temp_directory, suffix=".lst")
            with temp_reads_file as spf:
                with open(spf, "w") as f:
                    f.write("\n".join(subsample))
                temp_read.read_filter_inplace(spf)

            self.assertEqual(temp_read.get_current_fastq_read_number(), len(subsample))

            #
        r1.enrich(read_names_subset)  ## keep only the reads in the list
        self.assertEqual(r1.read_number_enriched, len(read_names_subset))
        self.assertEqual(r1.current_status, "enriched")
        r1.deplete(read_names_subset)  ## remove the reads in the list
        self.assertEqual(r1.current_status, "depleted")
        self.assertEqual(r1.read_number_depleted, 0)


class Televir_Project_Test(TestCase):
    software = SoftwareUtils()
    utils = Utils()
    constants = Constants()
    software_names = SoftwareNames()
    constants_tests_case = ConstantsTestsCase()
    pipeline_makeup = Pipeline_Makeup()

    def setUp(self):
        self.baseDirectory = os.path.join(
            getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS
        )

        self.test_user = test_user()
        self.sample_ont = test_fastq_file(self.baseDirectory, self.test_user)
        self.assertTrue(self.sample_ont.path_name_1.name)

        self.project_ont = televir_test_project(self.test_user)
        self.ont_project_sample = televir_test_sample(self.project_ont, self.sample_ont)

        ######
        default_software = DefaultSoftware()
        default_software.test_all_defaults_pathogen_identification(self.test_user)
        duplicate_software_params_global_project(self.test_user, self.project_ont)

    def test_project_trees_exist_ont(self):
        software_tree_utils = SoftwareTreeUtils(self.test_user, self.project_ont)

        for _, makeup in self.pipeline_makeup.MAKEUP.items():
            reset_project_makeup(self.project_ont)
            set_project_makeup(self.project_ont, makeup)

            local_tree = software_tree_utils.generate_project_tree()

            self.assertEqual(local_tree.__class__.__name__, "PipelineTree")

        reset_project_makeup(self.project_ont)

    @tag("slow")
    def test_all_project_paths_matched_ont(self):
        utils_manager = Utils_Manager()
        software_tree_utils = SoftwareTreeUtils(self.test_user, self.project_ont)
        for _, makeup in self.pipeline_makeup.MAKEUP.items():
            reset_project_makeup(self.project_ont)
            set_project_makeup(self.project_ont, makeup)

            local_tree = software_tree_utils.generate_project_tree()

            combined_table = (
                utils_manager.parameter_util.generate_merged_table_safe(
                    self.test_user, self.project_ont.technology
                )
            )

            if not utils_manager.check_pipeline_possible(
                combined_table, local_tree.makeup
            ):
                continue

            local_paths = local_tree.get_all_graph_paths_explicit()

            tree_makeup = local_tree.makeup

            pipeline_tree = software_tree_utils.generate_software_tree_extend(
                local_tree
            )
            pipeline_tree_index = software_tree_utils.get_software_tree_index(
                tree_makeup
            )

            matched_paths = {
                leaf: utils_manager.utility_manager.match_path_to_tree_safe(
                    path, pipeline_tree
                )
                for leaf, path in local_paths.items()
            }
            available_paths = {
                leaf: path for leaf, path in matched_paths.items() if path is not None
            }

            self.assertEqual(len(local_paths), len(available_paths))

            for _, path in available_paths.items():
                node = SoftwareTreeNode.objects.filter(
                    software_tree__pk=pipeline_tree_index, index=path
                ).exists()

                self.assertTrue(node)

    def contained_test_workflow_software_exist(self, run: Run_Main_from_Leaf):
        run.container.run_engine.software_check_map()

        for module in run.container.run_engine.modules:
            self.assertTrue(
                run.container.run_engine.module_software_check_map[module],
            )

    @tag("slow")
    def test_run_submit(self):
        pipeline_makeup_manager = Pipeline_Makeup()
        software_tree_utils = SoftwareTreeUtils(self.test_user, self.project_ont)

        for makeup, makeup_explicit in self.pipeline_makeup.MAKEUP.items():
            reset_project_makeup(self.project_ont)
            set_project_makeup(self.project_ont, makeup_explicit)

            pipeline_tree, available_paths = determine_available_paths(
                self.project_ont, self.test_user
            )

            pipeline_tree_index = software_tree_utils.get_software_tree_index(makeup)

            available_path_nodes = {
                leaf: SoftwareTreeNode.objects.get(
                    software_tree__pk=pipeline_tree_index, index=path
                )
                for leaf, path in available_paths.items()
            }

            for _, matched_path_node in available_path_nodes.items():
                run = Run_Main_from_Leaf(
                    user=self.test_user,
                    input_data=self.ont_project_sample,
                    project=self.project_ont,
                    pipeline_leaf=matched_path_node,
                    pipeline_tree=pipeline_tree,
                    odir=self.baseDirectory,
                    threads=3,
                )
                run.get_in_line()

                self.assertEqual(run.get_status(), ParameterSet.STATUS_QUEUED)
                self.assertEqual(run.check_finished(), False)
                self.assertEqual(run.check_availability(), True)
                self.assertEqual(run.is_available, True)
                self.assertEqual(run.container.prepped, False)
                self.assertEqual(run.check_processed(), False)
                self.assertEqual(run.check_submission(), False)
                self.assertEqual(
                    run.container.dir
                    == os.path.join(
                        self.baseDirectory, run.deployment_directory_structure
                    ),
                    True,
                )

                self.assertEqual(
                    makeup_explicit,
                    pipeline_makeup_manager.get_makeup(run.container.tree_makup),
                )

                ##### Test submission
                run.register_submission()

                self.assertEqual(run.check_submission(), True)
                configured = run.configure()

                self.assertEqual(configured, True)
                self.assertEqual(
                    run.container.config["sample_name"], run.container.sample.name
                )
                self.assertEqual(
                    os.path.join(
                        os.path.join(run.container.dir, "reads"),
                        os.path.basename(run.file_r1),
                    ),
                    run.container.config["r1"],
                )
                self.assertEqual(run.container.config["type"], PI_CS.SINGLE_END)
                self.assertEqual(os.path.exists(run.container.dir), True)
                self.assertEqual(
                    os.path.exists(
                        os.path.join(PI_CS.media_directory, run.container.dir_branch)
                    ),
                    True,
                )
                self.assertEqual(
                    os.path.exists(
                        os.path.join(PI_CS.static_directory, run.container.dir_branch)
                    ),
                    True,
                )

                for sbdir in run.container.config["directories"].values():
                    self.assertEqual(os.path.exists(sbdir), True)

                ##### test prepping
                run.container.run_main_prep()
                self.assertEqual(run.container.prepped, True)
                self.assertEqual(
                    os.path.exists(run.container.run_engine.media_dir_classification),
                    True,
                )
                self.assertEqual(
                    os.path.exists(run.container.run_engine.media_dir_igv), True
                )

                self.contained_test_workflow_software_exist(run)

                ##### test closing
                run.container.close()
                self.assertFalse(
                    os.path.exists(run.container.dir),
                )
