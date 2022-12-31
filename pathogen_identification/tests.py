from django.test import TestCase
from django.contrib.auth.models import User
from pathogen_identification.models import Projects
from pathogen_identification.models import PIProject_Sample
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager
from pathogen_identification.constants_settings import Pipeline_Makeup
from settings.default_software import DefaultSoftware
import os
from django.conf import settings
from settings.models import Software, Parameter, Sample
from settings.constants_settings import ConstantsSettings as CS
from pathogen_identification.models import SoftwareTree, ParameterSet, SoftwareTreeNode
from constants.constantsTestsCase import ConstantsTestsCase
from utils.software import Software as SoftwareUtils
from utils.utils import Utils
from constants.software_names import SoftwareNames
from constants.constants import Constants
from typing import Tuple, Dict
from pathogen_identification.deployment_main import Run_Main_from_Leaf
from pathogen_identification.modules.object_classes import (
    Operation_Temp_Files,
    Temp_File,
    RunCMD,
    Read_class,
    Sample_runClass,
    Software_detail,
    Bedgraph,
)
from pathogen_identification.install_registry import Deployment_Params

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
        type_of_use=Software.TYPE_OF_USE_televir_global,
        type_of_software=Software.TYPE_SOFTWARE,
        is_obsolete=False,
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
        type_of_use=Software.TYPE_OF_USE_televir_global,
        type_of_software=Software.TYPE_SOFTWARE,
        is_obsolete=False,
    )

    for software in query_set:

        software_parameters = Parameter.objects.filter(
            software=software,
        )

        software.pk = None
        software.type_of_use = Software.TYPE_OF_USE_televir_project

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
        type_of_use=Software.TYPE_OF_USE_televir_project,
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


def get_global_tree(technology, tree_makeup):
    utils_manager = Utils_Manager()
    pipeline_tree = utils_manager.generate_software_tree(technology, tree_makeup)
    pipeline_tree_index = utils_manager.get_software_tree_index(technology, tree_makeup)
    pipeline_tree_query = SoftwareTree.objects.get(pk=pipeline_tree_index)

    return pipeline_tree, pipeline_tree_query


def determine_available_paths(project: Projects, user: User):
    utils_manager = Utils_Manager()

    local_tree = utils_manager.generate_project_tree(project.technology, project, user)
    local_paths = local_tree.get_all_graph_paths_explicit()

    tree_makeup = local_tree.makeup
    pipeline_tree, pipeline_tree_query = get_global_tree(
        project.technology, tree_makeup
    )

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
        self.assertEqual(centrifuge.count(), 3)

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

        snippy = Software.objects.filter(name_extended="Snippy", owner=self.test_user)
        self.assertEqual(snippy.count(), 1)

        bowtie2 = Software.objects.filter(name_extended="Bowtie2", owner=self.test_user)
        self.assertEqual(bowtie2.count(), 1)

    def test_pipeline_makeup(self):
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
                "Remapping",
                "Contig classification",
                "Read classification",
                "Viral enrichment",
            },
            set(pipeline_excluding_raven),
        )

    def test_generate_default_trees(self):
        default_software = DefaultSoftware()
        default_software.test_all_defaults_pathogen_identification(self.test_user)
        utils_manager = Utils_Manager()
        utils_manager.generate_default_trees()
        all_trees = SoftwareTree.objects.all()

        self.assertEqual(all_trees.count(), len(self.pipeline_makeup.MAKEUP) * 2)

        ont_pipeline_tree = utils_manager.generate_software_tree(
            CS.TECHNOLOGY_minion, 0
        )

        self.assertEqual(ont_pipeline_tree.__class__.__name__, "PipelineTree")


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
        baseDirectory,
        ConstantsTestsCase.DIR_FASTQ,
        ConstantsTestsCase.FASTQ_MINION_1,
    )

    temp_dir = utils.get_temp_dir()
    utils.copy_file(file_name, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))

    try:
        sample_ont = Sample.objects.get(name=sample_name)
    except Sample.DoesNotExist:
        sample_ont = Sample()
        sample_ont.name = sample_name
        sample_ont.is_valid_1 = True
        sample_ont.file_name_1 = ConstantsTestsCase.FASTQ1_1
        sample_ont.path_name_1.name = os.path.join(
            temp_dir, ConstantsTestsCase.FASTQ1_1
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

        self.test_user = test_user()

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
        bindir = get_bindir_from_binaries(
            self.install_registry.BINARIES, CS.PIPELINE_NAME_read_quality_analysis
        )
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
                ## compile
                cmd = f"javac {os.path.basename(tpf)}"
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
        utils_manager = Utils_Manager()
        utils_manager.generate_default_trees()
        duplicate_software_params_global_project(self.test_user, self.project_ont)

    def test_project_trees_exist_ont(self):
        utils_manager = Utils_Manager()

        for tree, makeup in self.pipeline_makeup.MAKEUP.items():
            reset_project_makeup(self.project_ont)
            set_project_makeup(self.project_ont, makeup)

            local_tree = utils_manager.generate_project_tree(
                self.project_ont.technology, self.project_ont, self.test_user
            )

            self.assertEqual(local_tree.__class__.__name__, "PipelineTree")

        reset_project_makeup(self.project_ont)

    def test_all_project_paths_matched_ont(self):
        utils_manager = Utils_Manager()

        for tree, makeup in self.pipeline_makeup.MAKEUP.items():
            reset_project_makeup(self.project_ont)
            set_project_makeup(self.project_ont, makeup)

            local_tree = utils_manager.generate_project_tree(
                self.project_ont.technology, self.project_ont, self.test_user
            )
            local_paths = local_tree.get_all_graph_paths_explicit()

            tree_makeup = local_tree.makeup

            pipeline_tree = utils_manager.generate_software_tree(
                self.project_ont.technology, tree_makeup
            )
            pipeline_tree_index = utils_manager.get_software_tree_index(
                self.project_ont.technology, tree_makeup
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

            for leaf, path in available_paths.items():
                node = SoftwareTreeNode.objects.filter(
                    software_tree__pk=pipeline_tree_index, index=path
                ).exists()

                self.assertTrue(node)

    def test_run_submit(self):
        utils_manager = Utils_Manager()

        for tree, makeup in self.pipeline_makeup.MAKEUP.items():
            reset_project_makeup(self.project_ont)
            set_project_makeup(self.project_ont, makeup)

            pipeline_tree, available_paths = determine_available_paths(
                self.project_ont, self.test_user
            )
            pipeline_tree_index = utils_manager.get_software_tree_index(
                self.project_ont.technology, tree
            )

            available_path_nodes = {
                leaf: SoftwareTreeNode.objects.get(
                    software_tree__pk=pipeline_tree_index, index=path
                )
                for leaf, path in available_paths.items()
            }

            for leaf, matched_path_node in available_path_nodes.items():

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

                self.assertEqual(run.parameter_set.status, ParameterSet.STATUS_QUEUED)
