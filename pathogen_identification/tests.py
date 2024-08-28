import logging
import os
import random
import shutil
from typing import Dict, List, Optional, Tuple

import numpy as np
# Create your tests here.
import pandas as pd
from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase, tag

from constants.constants import Constants
from constants.constants import Televir_Metadata_Constants as Deployment_Params
from constants.constants import TypePath
from constants.constantsTestsCase import ConstantsTestsCase
from constants.software_names import SoftwareNames
from fluwebvirus.settings import STATIC_ROOT
from pathogen_identification.constants_settings import \
    ConstantsSettings as PI_CS
from pathogen_identification.deployment_main import Run_Main_from_Leaf
from pathogen_identification.models import (FinalReport, ParameterSet,
                                            PIProject_Sample, Projects,
                                            RawReference, ReferenceSource,
                                            ReferenceSourceFile,
                                            ReferenceSourceFileMap,
                                            ReferenceTaxid, RunMain,
                                            SoftwareTree, SoftwareTreeNode)
from pathogen_identification.modules.metadata_handler import RunMetadataHandler
from pathogen_identification.modules.object_classes import (
    Operation_Temp_Files, Read_class, RunCMD, Temp_File)
from pathogen_identification.utilities.overlap_manager import (
    MappingResultsParser, ReadOverlapManager, clade_private_proportions,
    pairwise_shared_count, pairwise_shared_reads,
    pairwise_shared_reads_distance, square_and_fill_diagonal,
    very_similar_groups_from_dataframe)
from pathogen_identification.utilities.reference_utils import extract_file
from pathogen_identification.utilities.televir_parameters import \
    TelevirParameters
from pathogen_identification.utilities.tree_deployment import Tree_Progress
from pathogen_identification.utilities.utilities_general import merge_classes
from pathogen_identification.utilities.utilities_pipeline import (
    Pipeline_Makeup, PipelineTree, SoftwareTreeUtils, Utils_Manager)
from pathogen_identification.utilities.utilities_views import (
    FinalReportCompound, RawReferenceUtils, ReportSorter,
    SampleReferenceManager)
from settings.constants_settings import ConstantsSettings as CS
from settings.default_software import DefaultSoftware
from settings.models import Parameter, Sample, Software
from utils.software import Software as SoftwareUtils
from utils.utils import Utils


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


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


def televir_test_project_ont(user: User, project_ont_name: str = "project_televir"):
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


def televir_test_project_illu(
    user: User, project_illu_name: str = "project_televir_illu"
):
    ######### ILLUMINA ##########

    try:
        project_illu = Projects.objects.get(
            name=project_illu_name, technology=CS.TECHNOLOGY_illumina
        )
    except Projects.DoesNotExist:
        project_illu = Projects()
        project_illu.name = project_illu_name
        project_illu.owner = user
        project_illu.technology = CS.TECHNOLOGY_illumina

        project_illu.save()

    return project_illu


def televir_test_sample(project, sample: Sample):
    try:
        project_sample = PIProject_Sample.objects.get(
            project__id=project.pk, sample__id=sample.pk
        )
    except PIProject_Sample.DoesNotExist:
        project_sample = PIProject_Sample()
        project_sample.project = project
        project_sample.sample = sample
        project_sample.name = sample.name
        project_sample.save()

    return project_sample


def test_fastq_ont_file(
    baseDirectory, user: User, sample_name: str = "televir_sample_minion_1"
):
    utils = Utils()

    file_name = os.path.join(
        STATIC_ROOT,
        ConstantsTestsCase.MANAGING_TESTS,
        ConstantsTestsCase.DIR_FASTQ,
        ConstantsTestsCase.FASTQ_MINION_1,
    )

    utils.copy_file(
        file_name, os.path.join(baseDirectory, ConstantsTestsCase.FASTQ_MINION_1)
    )

    try:
        sample_ont = Sample.objects.get(
            name=sample_name, type_of_fastq=Sample.TYPE_OF_FASTQ_minion
        )
    except Sample.DoesNotExist:
        sample_ont = Sample()
        sample_ont.name = sample_name
        sample_ont.is_valid_1 = True
        sample_ont.file_name_1 = ConstantsTestsCase.FASTQ_MINION_1
        sample_ont.path_name_1.name = os.path.join(
            baseDirectory, ConstantsTestsCase.FASTQ_MINION_1
        )
        sample_ont.is_valid_2 = False
        sample_ont.type_of_fastq = Sample.TYPE_OF_FASTQ_minion
        sample_ont.owner = user
        sample_ont.save()

    return sample_ont


def test_fastq_illu1_file(
    baseDirectory, user: User, sample_name: str = "televir_sample_illumina_1"
) -> Sample:
    utils = Utils()

    file_name = os.path.join(
        STATIC_ROOT,
        ConstantsTestsCase.MANAGING_TESTS,
        ConstantsTestsCase.DIR_FASTQ,
        ConstantsTestsCase.FASTQ1_1,
    )

    utils.copy_file(file_name, os.path.join(baseDirectory, ConstantsTestsCase.FASTQ1_1))

    file_name = os.path.join(
        STATIC_ROOT,
        ConstantsTestsCase.MANAGING_TESTS,
        ConstantsTestsCase.DIR_FASTQ,
        ConstantsTestsCase.FASTQ1_2,
    )

    utils.copy_file(file_name, os.path.join(baseDirectory, ConstantsTestsCase.FASTQ1_2))

    try:
        sample_illu1 = Sample.objects.get(
            name=sample_name, type_of_fastq=Sample.TYPE_OF_FASTQ_illumina
        )
    except Sample.DoesNotExist:
        sample_illu1 = Sample()
        sample_illu1.name = sample_name
        sample_illu1.is_valid_1 = True
        sample_illu1.file_name_1 = ConstantsTestsCase.FASTQ1_1
        sample_illu1.path_name_1.name = os.path.join(
            baseDirectory, ConstantsTestsCase.FASTQ1_1
        )

        sample_illu1.is_valid_2 = True
        sample_illu1.file_name_2 = ConstantsTestsCase.FASTQ1_2
        sample_illu1.path_name_2.name = os.path.join(
            baseDirectory, ConstantsTestsCase.FASTQ1_2
        )
        sample_illu1.type_of_fastq = Sample.TYPE_OF_FASTQ_illumina
        sample_illu1.owner = user
        sample_illu1.save()

    return sample_illu1


def generate_influenza_reference(
    output_file: str, accid_in_file: str, accid_to_be: Optional[str] = None
):
    deployment_params = Deployment_Params()
    output_file_gz = output_file + ".gz"

    source = os.path.join(
        STATIC_ROOT,
        ConstantsTestsCase.MANAGING_TESTS,
        ConstantsTestsCase.MANAGING_DIR,
        ConstantsTestsCase.MANAGING_FILES_FASTA,
    )

    samtools_bin = os.path.join(
        Deployment_Params.BINARIES["ROOT"],
        Deployment_Params.BINARIES["software"]["java"],
        "bin",
        "samtools",
    )

    if os.path.exists(output_file + ".gz"):
        os.remove(output_file + ".gz")

    cmd = f"{samtools_bin} faidx {source} {accid_in_file} > {output_file}"

    os.system(cmd)

    ### replace accid in file
    if accid_to_be is not None:
        cmd = f"sed -i 's/{accid_in_file}/{accid_to_be}/g' {output_file}"
        os.system(cmd)

    ## compress and index
    cmd = f"{deployment_params.get_software_binary('bgzip')} {output_file}"
    os.system(cmd)

    cmd = f"{deployment_params.get_software_binary('samtools')} faidx {output_file_gz}"
    os.system(cmd)

    return output_file_gz


def generate_NA_televir_reference(output_dir: str):

    taxid_str = "1674969"
    accid = "KT022278.1"
    accid_in_file = "NA"
    reference_file = "NA_H3N2.fasta"
    description = "Influenza A virus (A/chicken/Guangxi/125C8/2012(H3N2)) segment 6 neuraminidase (NA) gene, complete cds"
    output_file = os.path.join(output_dir, reference_file)

    test_reference_path = generate_influenza_reference(
        output_file, accid_in_file, accid_to_be=accid
    )
    test_reference = os.path.basename(test_reference_path)

    ####
    taxid = ReferenceTaxid()
    taxid.taxid = taxid_str
    taxid.save()
    ref_source = ReferenceSource()
    ref_source.accid = accid
    ref_source.taxid = taxid
    ref_source.description = description
    ref_source.save()

    ref_source_file = ReferenceSourceFile()
    ref_source_file.file = test_reference
    ref_source_file.save()

    ref_source_file = ReferenceSourceFile.objects.get(file=test_reference)

    ref_source_file_map = ReferenceSourceFileMap()
    ref_source_file_map.reference_source = ref_source
    ref_source_file_map.reference_source_file = ref_source_file
    ref_source_file_map.save()


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
        parameter__televir_project=None,
        parameter__televir_project_sample=None,
    )
    project = Projects.objects.get(pk=project.pk)
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
                parameter__televir_project_sample=None,
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


def duplicate_software_params_global_project(user, project: Projects):
    """
    duplicate software global to project
    """
    ### get all global software
    query_set = Software.objects.filter(
        owner=user,
        type_of_use__in=Software.TELEVIR_GLOBAL_TYPES,
        type_of_software__in=[
            Software.TYPE_SOFTWARE,
            Software.TYPE_INSAFLU_PARAMETER,
        ],
        is_obsolete=False,
        technology__name=project.technology,
        parameter__televir_project=None,
        parameter__televir_project_sample=None,
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
                parameter__televir_project_sample=None,
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


def user_project_turn_off_pipeline_steps(user, project, pipeline_steps):
    """
    turn off pipeline steps
    """

    software = Software.objects.filter(
        owner=user,
        parameter__televir_project=project,
        pipeline_step__name__in=pipeline_steps,
    )

    for software in software:
        software.is_to_run = False
        software.save()


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


def generate_compressed_tree(user, project, sample, makeup):
    software_tree_utils = SoftwareTreeUtils(user, project, sample)
    utils_manager = Utils_Manager()

    runs_to_deploy: Dict[PIProject_Sample, List[SoftwareTreeNode]] = (
        software_tree_utils.check_runs_to_deploy_sample(sample)
    )

    local_tree = software_tree_utils.generate_project_tree()
    local_paths = local_tree.get_all_graph_paths_explicit()

    pipeline_tree = software_tree_utils.generate_software_tree_extend(local_tree)
    pipeline_tree_index = local_tree.software_tree_pk

    matched_paths = {
        leaf: utils_manager.utility_manager.match_path_to_tree_safe(path, pipeline_tree)
        for leaf, path in local_paths.items()
    }

    assert set(list(matched_paths.keys())) == set(
        [x.index for x in runs_to_deploy[sample]]
    )

    available_path_nodes = {
        leaf: SoftwareTreeNode.objects.get(
            software_tree__pk=pipeline_tree_index, index=path
        )
        for leaf, path in matched_paths.items()
    }

    available_path_nodes = {
        leaf: utils_manager.parameter_util.check_ParameterSet_available_to_run(
            sample=sample,
            leaf=matched_path_node,
            project=project,
        )
        for leaf, matched_path_node in available_path_nodes.items()
    }

    matched_paths = {
        k: v for k, v in matched_paths.items() if available_path_nodes[k] == True
    }

    assert set(list(matched_paths.keys())) == set(
        [x.index for x in runs_to_deploy[sample]]
    )

    module_tree = utils_manager.module_tree(pipeline_tree, list(matched_paths.values()))

    return module_tree


#########################################################################
#########################################################################
###### Tests Begin Here


class OverlapManagerTests(TestCase):
    def setUp(self):
        # Setup test data
        self.baseDirectory = os.path.join(
            STATIC_ROOT, ConstantsTestsCase.MANAGING_TESTS
        )

        self.test_user = test_user()
        self.project_ont = televir_test_project_ont(self.test_user)

        self.temp_directory = os.path.join(self.baseDirectory, PI_CS.test_subdirectory)
        os.makedirs(self.temp_directory, exist_ok=True)

        self.test_matrix = pd.DataFrame(
            [[1, 0, 1, 1], [0, 1, 1, 0], [1, 1, 0, 1]],
            index=["id1", "id2", "id3"],
            columns=["feature1", "feature2", "feature3", "feature4"],
        )

        ## create files and metadata, iterate index
        metadata = []
        for idx in self.test_matrix.index:
            with open(os.path.join(self.temp_directory, f"{idx}.fasta"), "w") as f:
                for feature in self.test_matrix.columns:
                    if self.test_matrix.loc[idx, feature] == 1:
                        f.write(f">{feature}\n")

                metadata.append(
                    {
                        "filename": f"{idx}.fasta",
                        "file": os.path.join(self.temp_directory, f"{idx}.fasta"),
                        "accid": idx,
                        "description": f"desc{idx}",
                        "bam": f"{idx}.bam",
                    }
                )

        self.metadata_df = pd.DataFrame(metadata)

    def test_pairwise_shared_count(self):
        result = pairwise_shared_count(self.test_matrix)
        expected_result = pd.DataFrame(
            [[3, 1, 2], [1, 2, 1], [2, 1, 3]],
            index=["id1", "id2", "id3"],
            columns=["id1", "id2", "id3"],
        )

        pd.testing.assert_frame_equal(result, expected_result)

    def test_square_and_fill_diagonal(self):
        result = square_and_fill_diagonal(self.test_matrix)
        # Calculate expected result
        shared_counts = pairwise_shared_count(self.test_matrix)
        row_sums = self.test_matrix.sum(axis=1)
        expected_result = shared_counts.div(row_sums, axis=0)
        np.fill_diagonal(expected_result.values, 0)
        expected_result = expected_result.fillna(0)

        pd.testing.assert_frame_equal(result, expected_result)

    def test_very_similar_groups_from_dataframe(self):
        result = very_similar_groups_from_dataframe(self.test_matrix, threshold=0.6)

        expected_result = [
            ("id1", "id3"),
            ("id2",),
        ]
        self.assertEqual(sorted(result), sorted(expected_result))

    def test_clade_private_proportions(self):
        private_reads, total_reads, proportion_private = clade_private_proportions(
            self.test_matrix, ["id1", "id2"]
        )
        expected_private_reads = 1  # Replace with actual expected value
        expected_total_reads = 4  # Replace with actual expected value
        expected_proportion_private = 0.25  # Replace with actual expected value
        self.assertEqual(private_reads, expected_private_reads)
        self.assertEqual(total_reads, expected_total_reads)
        self.assertAlmostEqual(proportion_private, expected_proportion_private)

    def test_pairwise_shared_reads_distance(self):
        result = pairwise_shared_reads_distance(self.test_matrix)
        one_thirds = 1 / 3
        expected_result = pd.DataFrame(
            [[0, 0.5, one_thirds], [0.5, 0, 0.5], [one_thirds, 0.5, 0]],
            index=["id1", "id2", "id3"],
            columns=["id1", "id2", "id3"],
        )
        pd.testing.assert_frame_equal(result, expected_result)

    def test_pairwise_shared_reads(self):
        result = pairwise_shared_reads(self.test_matrix)

        one_thirds = 1 / 3
        two_thirds = 2 / 3
        expected_result = pd.DataFrame(
            [[0, one_thirds, two_thirds], [0.5, 0, 0.5], [two_thirds, one_thirds, 0]],
            index=["id1", "id2", "id3"],
            columns=["id1", "id2", "id3"],
        )
        pd.testing.assert_frame_equal(result, expected_result)

    def dont_test_parse_for_data(self):

        mapping_parser = MappingResultsParser(
            self.metadata_df, self.temp_directory, "test_pid"
        )

        mapping_parser.parse_for_data()

        self.assertTrue(mapping_parser.parsed)

        pd.testing.assert_frame_equal(
            mapping_parser.read_profile_matrix,
            mapping_parser.read_profile_matrix_filtered,
        )

        self.assertTrue(mapping_parser.total_read_counts.sum() == 8)

    def test_televir_parameters(self):

        report_layout_params = TelevirParameters.get_report_layout_params(
            project_pk=self.project_ont.pk
        )

        self.assertTrue(report_layout_params.shared_proportion_threshold == 0.2)
        self.assertTrue(report_layout_params.read_overlap_threshold == 0.5)

        self.assertTrue(
            report_layout_params.flag_str == PI_CS.FLAG_BUILD_DEFAULT.build_name
        )

    def test_overlap_manager(self):

        report_layout_params = TelevirParameters.get_report_layout_params(
            project_pk=self.project_ont.pk
        )

        reference_clade = ReportSorter.generate_reference_clade(report_layout_params)

        overlap_manager = ReadOverlapManager(
            self.metadata_df, reference_clade, self.temp_directory, "test_pid"
        )

        self.assertTrue(overlap_manager.all_accs_analyzed())

        distance_matrix = overlap_manager.generate_distance_matrix(
            force=overlap_manager.force_tree_rebuild
        )

        tree = overlap_manager.tree_from_distance_matrix(distance_matrix)

        self.assertTrue(tree is not None)
        self.assertFalse(tree.rooted)

        self.assertTrue(len(tree.root.clades) == 2)
        self.assertTrue(len(overlap_manager.all_clade_leaves_filtered) == 5)
        self.assertTrue(len(overlap_manager.get_leaf_clades()) == 3)


class MergeClassificationsTest(TestCase):
    def setUp(self):
        # Setup data frames similar to the pytest fixtures
        self.deployment_dir = os.path.join(
            STATIC_ROOT, ConstantsTestsCase.MANAGING_TESTS
        )
        self.test_user = test_user()
        self.data_frame_1 = pd.DataFrame(
            {
                "taxid": [1, 2, 3],
                "counts": [10, 20, 30],
                "description": ["btests_acteria", "virus", "phage"],
            }
        )

        self.data_frame_2 = pd.DataFrame(
            {
                "taxid": [3, 4, 5],
                "counts": [30, 40, 50],
                "description": ["phage", "fungi", "protozoa"],
            }
        )

    def test_non_overlapping_merge(self):
        merged, _ = merge_classes(
            self.data_frame_1.head(2), self.data_frame_2.tail(2), exclude="phage"
        )
        self.assertEqual(len(merged), 4, "Should merge non-overlapping data correctly.")

    def test_partial_overlap_merge(self):
        merged, _ = merge_classes(self.data_frame_1, self.data_frame_2, exclude="phage")
        self.assertEqual(
            len(merged),
            4,
            "Should correctly merge and handle partially overlapping data.",
        )

    def test_full_overlap_merge(self):
        merged, _ = merge_classes(self.data_frame_1, self.data_frame_1, exclude="phage")
        self.assertEqual(
            len(merged),
            2,
            "Should handle fully overlapping data correctly, excluding 'phage'.",
        )

    def test_exclusion_parameter(self):
        _, full_descriptor = merge_classes(self.data_frame_1, self.data_frame_2)
        self.assertNotIn(
            3,
            full_descriptor["taxid"].values,
            "Should exclude specified categories.",
        )

    def test_empty_dataframe_merge(self):
        empty_df = pd.DataFrame(columns=["taxid", "counts", "description"])
        merged, _ = merge_classes(self.data_frame_1, empty_df)
        self.assertEqual(
            len(merged),
            2,
            "Should handle empty data frames correctly, excluding 'phage'.",
        )

    def test_maxt_parameter(self):
        merged, _ = merge_classes(
            self.data_frame_1, self.data_frame_2, maxt=1, exclude="phage"
        )
        self.assertLessEqual(len(merged), 1, "Should respect the maxt parameter.")


class MetadataManagementTests(TestCase):

    def setUp(self):
        self.baseDirectory = os.path.join(
            STATIC_ROOT, ConstantsTestsCase.MANAGING_TESTS
        )
        self.nthreads = 3

        self.install_registry = Deployment_Params()

        self.test_user = test_user()
        self.temp_directory = os.path.join(self.baseDirectory, PI_CS.test_subdirectory)
        self.sample_illu = test_fastq_illu1_file(self.baseDirectory, self.test_user)
        self.sample_ont = test_fastq_ont_file(self.baseDirectory, self.test_user)

        self.project_ont = televir_test_project_ont(self.test_user, "project_ont")
        self.project_illu = televir_test_project_illu(self.test_user, "project_illu")
        self.project_sample_ont = televir_test_sample(self.project_ont, self.sample_ont)
        self.project_sample_illu = televir_test_sample(
            self.project_illu, self.sample_illu
        )
        self.test_mapping_reference_accid = "KT022278.1"
        generate_NA_televir_reference(self.temp_directory)

        refs = [
            {
                "accid": "ref1",
                "taxid": 1,
                "description": "virus1",
                "file": "tests_file1",
            },
            {
                "accid": "ref2",
                "taxid": 2,
                "description": "virus2",
                "file": "tests_file1",
            },
        ]

        self.reads_report = pd.DataFrame(
            {
                "taxid": [
                    "1674969",
                ],
                "counts": [
                    4000,
                ],
                "description": [
                    "Influenza A virus (A/chicken/Guangxi/125C8/2012(H3N2)) segment 6 neuraminidase (NA) gene, complete cds",
                ],
            }
        )

        self.contig_report = pd.DataFrame(
            [
                {
                    "taxid": "1674969",
                    "counts": 30,
                }
            ]
        )
        ### write files
        for ref in refs:
            with open(os.path.join(self.temp_directory, ref["file"]), "w") as f:
                f.write(f"{ref['accid']}\n")
                f.write(
                    "".join(np.random.choice(["A", "C", "G", "T"], 1000, replace=True))
                )

        files = list(set([ref["file"] for ref in refs]))
        for file in files:

            ref_source_file = ReferenceSourceFile()
            ref_source_file.file = file
            ref_source_file.save()

        ### generate references
        for ref in refs:
            taxid = ReferenceTaxid()
            taxid.taxid = ref["taxid"]
            taxid.save()
            ref_source = ReferenceSource()
            ref_source.accid = ref["accid"]
            ref_source.taxid = taxid
            ref_source.description = ref["description"]
            ref_source.save()

            ref_source_file = ReferenceSourceFile.objects.get(file=ref["file"])

            ref_source_file_map = ReferenceSourceFileMap()
            ref_source_file_map.reference_source = ref_source
            ref_source_file_map.reference_source_file = ref_source_file
            ref_source_file_map.save()

        self.dataframe_1 = pd.DataFrame(
            {
                "taxid": [1, 2, 3],
                "counts": [10, 20, 30],
            }
        )

        self.dataframe_2 = pd.DataFrame(
            {
                "taxid": [3, 4, 5],
                "counts": [30, 40, 50],
            }
        )

        ###########################
        ######
        default_software = DefaultSoftware()
        default_software.test_all_defaults_pathogen_identification(self.test_user)
        duplicate_software_params_global_project(self.test_user, self.project_ont)

    def test_add_reference(self):

        sample_ref_manager = SampleReferenceManager(self.project_sample_ont)

        self.assertTrue(RunMain.objects.filter(project=self.project_ont).exists())

        ref = ReferenceSourceFileMap.objects.get(reference_source__accid="KT022278.1")
        sample_ref_manager.add_reference(ref)

        self.assertTrue(
            RawReference.objects.filter(run=sample_ref_manager.storage_run).exists()
        )

        self.assertTrue(
            SoftwareTree.objects.filter(project=self.project_ont, model=-1).exists()
        )

        software_utils = SoftwareTreeUtils(
            self.test_user, self.project_ont, sample=self.project_sample_ont
        )
        runs_to_deploy = software_utils.check_runs_to_submit_metagenomics_sample(
            self.project_sample_ont
        )
        reference_manager = SampleReferenceManager(self.project_sample_ont)

        for leaf in runs_to_deploy[self.project_sample_ont]:

            for run_type in [
                RunMain.RUN_TYPE_COMBINED_MAPPING,
                RunMain.RUN_TYPE_MAP_REQUEST,
                RunMain.RUN_TYPE_SCREENING,
            ]:
                test_run = reference_manager.create_mapping_run(leaf, run_type)

                self.assertTrue(test_run.status == RunMain.STATUS_PREP)
                self.assertTrue(test_run.run_type == run_type)

    def dont_test_extract_file(self):
        tmp_fasta = extract_file("ref1")
        self.assertFalse(tmp_fasta is None)

        self.assertTrue(os.path.exists(tmp_fasta))

    def test_metadataHandler(self):

        config = {
            "metadata": self.install_registry.metadata_full_path,
            "bin": self.install_registry.BINARIES,
        }

        metadata_tool = RunMetadataHandler(
            self.test_user.username,
            config,
            sift_query="phage",
            prefix="test_run",
            rundir=self.temp_directory,
        )

        reads_report = pd.DataFrame(
            {
                "taxid": [
                    "1674969",
                    "1",
                ],
                "counts": [
                    4000,
                    1000,
                ],
                "description": [
                    "Influenza A virus (A/chicken/Guangxi/125C8/2012(H3N2)) segment 6 neuraminidase (NA) gene, complete cds",
                    "virus1",
                ],
            }
        )

        contig_report = pd.DataFrame(
            [
                {
                    "taxid": "1674969",
                    "counts": 30,
                }
            ]
        )

        metadata_tool.match_and_select_targets(
            reads_report, contig_report, max_remap=3, taxid_limit=3
        )

        self.assertTrue(metadata_tool.merged_targets.shape[0] == 2)
        self.assertTrue(
            metadata_tool.merged_targets.iloc[0]["source"] == "reads/contigs"
        )
        self.assertTrue(metadata_tool.merged_targets.iloc[1]["source"] == "reads")

        metadata_tool.merged_targets = pd.DataFrame(columns=["taxid", "counts"])
        metadata_tool.match_and_select_targets(
            pd.DataFrame(columns=["taxid", "counts"]),
            contig_report,
            max_remap=3,
            taxid_limit=3,
        )

        self.assertTrue(metadata_tool.merged_targets.iloc[0]["source"] == "contigs")
        self.assertTrue(metadata_tool.merged_targets.shape[0] == 1)
        ####
        metadata_tool.merged_targets = pd.DataFrame(columns=["taxid", "counts"])
        metadata_tool.match_and_select_targets(
            pd.DataFrame(columns=["taxid", "countsuy"]),
            contig_report,
            max_remap=3,
            taxid_limit=3,
        )

        self.assertTrue(metadata_tool.merged_targets.iloc[0]["source"] == "contigs")
        self.assertTrue(metadata_tool.merged_targets.shape[0] == 1)

        ###
        ####
        metadata_tool.merged_targets = pd.DataFrame(columns=["taxid", "counts"])
        metadata_tool.remap_targets = []
        metadata_tool.match_and_select_targets(
            self.reads_report,
            self.contig_report,
            max_remap=3,
            taxid_limit=3,
        )

        self.assertTrue(len(metadata_tool.remap_targets) == 1)

    def run_setup(self, sample_to_test: PIProject_Sample, project_to_test: Projects):

        reference_manager = SampleReferenceManager(sample_to_test)

        software_utils = SoftwareTreeUtils(
            self.test_user, project_to_test, sample=sample_to_test
        )

        runs_to_deploy, _ = software_utils.check_runs_to_submit_mapping_only(
            sample_to_test
        )

        first_run_node = runs_to_deploy[sample_to_test][0]

        references_to_add = ReferenceSourceFileMap.objects.get(
            reference_source__accid=self.test_mapping_reference_accid
        )

        reference_manager.add_reference(references_to_add)

        mapping_run = reference_manager.create_mapping_run(
            first_run_node, RunMain.RUN_TYPE_MAP_REQUEST
        )

        ref_storage = RawReference.objects.filter(
            run__sample=sample_to_test, run__run_type=RunMain.RUN_TYPE_STORAGE
        ).first()

        ref_storage.pk = None
        ref_storage.run = mapping_run
        ref_storage.save()

        run_engine = Run_Main_from_Leaf(
            self.test_user,
            sample_to_test,
            project_to_test,
            pipeline_leaf=first_run_node,
            pipeline_tree=first_run_node.software_tree,
            odir=self.temp_directory,
            combined_analysis=False,
            mapping_request=True,
            run_pk=mapping_run.pk,
        )

        ########
        ########
        run_engine.configure()
        ##
        self.assertTrue(
            run_engine.container.run_engine.metadata_tool.remap_targets[0].accid
            == ref_storage.accid
        )

        ## replace fasta file path with temp file
        for remap_target in run_engine.container.run_engine.metadata_tool.remap_targets:
            remap_target.file = os.path.join(
                self.temp_directory,
                os.path.basename(remap_target.file),
            )

        self.assertTrue(
            os.path.exists(
                run_engine.container.run_engine.metadata_tool.remap_targets[0].file
            )
        )

        return run_engine

    def run_mapping_to_reference(
        self, sample_to_test: PIProject_Sample, project_to_test: Projects
    ):

        run_engine = self.run_setup(sample_to_test, project_to_test)
        #### run mapping
        run_engine.container.run_engine.Prep_deploy()
        run_engine.container.run_engine.prep_REMAPPING()
        run_engine.container.run_engine.Run_Remapping()

        return run_engine

    @tag("slow")
    def test_mapping_to_reference(self):
        """
        Run for mappings against reference for both illumina and ont reads
        """

        run_engine = self.run_mapping_to_reference(
            self.project_sample_ont, self.project_ont
        )

        self.assertTrue(run_engine.container.run_engine.remap_manager.report.empty)

        #########################################################
        run_engine = self.run_mapping_to_reference(
            self.project_sample_illu, self.project_illu
        )

        self.assertTrue(
            run_engine.container.run_engine.remap_manager.report.shape[0] == 1
        )

    def test_deploy_mapping(self):
        """
        Run for mappings against reference for both illumina and ont reads"""

        #####
        user_project_turn_off_pipeline_steps(
            self.test_user, self.project_illu, [CS.PIPELINE_NAME_viral_enrichment]
        )
        run_engine = self.run_setup(self.project_sample_illu, self.project_illu)
        deployed_and_updated = run_engine.Deploy_Parts()

        self.assertTrue(deployed_and_updated)

        reset_project_makeup(self.project_illu)

        mapping_run = RunMain.objects.get(pk=run_engine.run_pk)

        final_report = FinalReport.objects.filter(
            sample=self.project_sample_illu, run=mapping_run
        ).order_by("-coverage")

        report_layout_params = TelevirParameters.get_report_layout_params(
            run_pk=run_engine.run_pk
        )
        reference_utils = RawReferenceUtils(project=self.project_illu)

        report_sorter = ReportSorter(
            self.project_sample_illu, final_report, report_layout_params
        )

        sorted_reports = report_sorter.get_reports()

        excluded_reports_exist = report_sorter.check_excluded_exist()
        # empty_reports = report_sorter.get_reports_empty()

        first_report_group = sorted_reports[0]
        self.assertEquals(first_report_group.max_coverage, 100)
        self.assertEquals(first_report_group.total_counts, "total counts 0")
        self.assertEquals(first_report_group.shared_proportion, 0)
        self.assertEquals(first_report_group.max_private_reads, 0)
        self.assertFalse(first_report_group.has_multiple)
        self.assertEquals(first_report_group.toggle, "on")

        self.assertFalse(excluded_reports_exist)

        sorted_reports = report_sorter.get_reports_compound()

        #############################
        query_string = "Influenza A virus"
        references = reference_utils.retrieve_compound_references(query_string)
        print("###################")
        print(references)

        first_report_group = sorted_reports[0]

        first_compound: FinalReportCompound = first_report_group.group_list[0]
        self.assertEquals(first_compound.found_in, "M")
        self.assertTrue(first_compound.run_main == mapping_run)
        self.assertTrue(first_compound.data_exists)
        self.assertEquals(first_compound.control_flag, FinalReport.CONTROL_FLAG_NONE)
        self.assertEquals(first_compound.private_reads, 0)


class Televir_Software_Test(TestCase):
    software = SoftwareUtils()
    utils = Utils()
    constants = Constants()
    software_names = SoftwareNames()
    constants_tests_case = ConstantsTestsCase()
    pipeline_makeup = Pipeline_Makeup()

    def setUp(self):
        self.baseDirectory = os.path.join(
            STATIC_ROOT, ConstantsTestsCase.MANAGING_TESTS
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
        self.assertEqual(kraken2.count(), 5)

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
                CS.PIPELINE_NAME_metagenomics_settings,
                CS.PIPELINE_NAME_request_mapping,
                CS.PIPELINE_NAME_viral_enrichment,
                CS.PIPELINE_NAME_contig_classification,
                CS.PIPELINE_NAME_read_classification,
                CS.PIPELINE_NAME_remapping,
            },
            set(pipeline_excluding_raven),
        )


class Televir_Objects_TestCase(TestCase):
    install_registry = Deployment_Params

    def setUp(self):
        self.baseDirectory = os.path.join(
            STATIC_ROOT, ConstantsTestsCase.MANAGING_TESTS
        )
        self.temp_directory = os.path.join(self.baseDirectory, PI_CS.test_subdirectory)
        os.makedirs(self.temp_directory, exist_ok=True)

        self.test_user = test_user()
        self.sample_ont = test_fastq_ont_file(self.temp_directory, self.test_user)

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

        read_names_subset = random.sample(read_names, 100)
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
            subsample = random.sample(read_names_subset, 10)
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
            STATIC_ROOT, ConstantsTestsCase.MANAGING_TESTS
        )

        self.test_user = test_user()
        self.project_ont = televir_test_project_ont(self.test_user)
        self.sample_ont = test_fastq_ont_file(self.baseDirectory, self.test_user)
        self.assertTrue(self.sample_ont.path_name_1.name)

        self.ont_project_sample = televir_test_sample(self.project_ont, self.sample_ont)

        ######
        default_software = DefaultSoftware()
        default_software.test_all_defaults_pathogen_identification(self.test_user)
        duplicate_software_params_global_project(self.test_user, self.project_ont)

    @tag("slow")
    def test_project_trees_exist_ont(self):
        software_tree_utils = SoftwareTreeUtils(self.test_user, self.project_ont)

        for _, makeup in self.pipeline_makeup.MAKEUP.items():
            reset_project_makeup(self.project_ont)
            set_project_makeup(self.project_ont, makeup)
            if (
                self.pipeline_makeup.check_makeuplist_has_classification(makeup)
                is False
            ):
                continue
            local_tree = software_tree_utils.generate_project_tree()

            self.assertEqual(local_tree.__class__.__name__, "PipelineTree")

        reset_project_makeup(self.project_ont)

    @tag("slow")
    def test_all_project_paths_matched_ont(self):
        utils_manager = Utils_Manager()
        software_tree_utils = SoftwareTreeUtils(self.test_user, self.project_ont)
        for _, makeup in self.pipeline_makeup.MAKEUP.items():
            if (
                self.pipeline_makeup.check_makeuplist_has_classification(makeup)
                is False
            ):
                continue

            reset_project_makeup(self.project_ont)
            set_project_makeup(self.project_ont, makeup)

            local_tree = software_tree_utils.generate_project_tree()

            combined_table = utils_manager.parameter_util.generate_merged_table_safe(
                self.test_user, self.project_ont.technology
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

    @tag("slow")
    def test_compress_tree(self):
        utils_manager = Utils_Manager()
        software_tree_utils = SoftwareTreeUtils(
            self.test_user, self.project_ont, self.ont_project_sample
        )
        for _, makeup in self.pipeline_makeup.MAKEUP.items():
            if (
                self.pipeline_makeup.check_makeuplist_has_classification(makeup)
                is False
            ):
                continue

            reset_project_makeup(self.project_ont)
            set_project_makeup(self.project_ont, makeup)

            runs_to_deploy: Dict[PIProject_Sample, List[SoftwareTreeNode]] = (
                software_tree_utils.check_runs_to_deploy_sample(self.ont_project_sample)
            )

            local_tree = software_tree_utils.generate_project_tree()
            local_paths = local_tree.get_all_graph_paths_explicit()

            pipeline_tree = software_tree_utils.generate_software_tree_extend(
                local_tree
            )
            pipeline_tree_index = local_tree.software_tree_pk

            matched_paths = {
                leaf: utils_manager.utility_manager.match_path_to_tree_safe(
                    path, pipeline_tree
                )
                for leaf, path in local_paths.items()
            }

            assert set(list(matched_paths.keys())) == set(
                [x.index for x in runs_to_deploy[self.ont_project_sample]]
            )

            available_path_nodes = {
                leaf: SoftwareTreeNode.objects.get(
                    software_tree__pk=pipeline_tree_index, index=path
                )
                for leaf, path in matched_paths.items()
            }

            available_path_nodes = {
                leaf: utils_manager.parameter_util.check_ParameterSet_available_to_run(
                    sample=self.ont_project_sample,
                    leaf=matched_path_node,
                    project=self.project_ont,
                )
                for leaf, matched_path_node in available_path_nodes.items()
            }

            matched_paths = {
                k: v
                for k, v in matched_paths.items()
                if available_path_nodes[k] == True
            }

            assert set(list(matched_paths.keys())) == set(
                [x.index for x in runs_to_deploy[self.ont_project_sample]]
            )

            module_tree = utils_manager.module_tree(
                pipeline_tree, list(matched_paths.values())
            )

            for leaf in module_tree.leaves:
                self.assertTrue(leaf in module_tree.compress_dag_dict)
                self.assertTrue(module_tree.compress_dag_dict[leaf] == [])

            modules_passed = []

            for node_index in module_tree.compress_dag_dict.keys():

                node_info = module_tree.node_index.node[node_index]

                if node_index == 0:
                    self.assertTrue(node_info[0] == "root")
                    continue

                if node_index in module_tree.leaves:
                    continue

                self.assertTrue(node_info[0] in makeup)
                modules_passed.append(node_info[0])

            self.assertEqual(set(modules_passed), set(makeup))

    @tag("slow")
    def test_progress_tree_tree(self):
        utils_manager = Utils_Manager()
        software_tree_utils = SoftwareTreeUtils(
            self.test_user, self.project_ont, self.ont_project_sample
        )
        for _, makeup in self.pipeline_makeup.MAKEUP.items():
            if (
                self.pipeline_makeup.check_makeuplist_has_classification(makeup)
                is False
            ):
                continue

            reset_project_makeup(self.project_ont)
            set_project_makeup(self.project_ont, makeup)

            module_tree = generate_compressed_tree(
                self.test_user, self.project_ont, self.ont_project_sample, makeup
            )

            deployment_tree = Tree_Progress(
                module_tree, self.ont_project_sample, self.project_ont
            )

            self.assertTrue(deployment_tree.get_current_module() == "root")

            deployment_tree.deploy_nodes()

            first_node = deployment_tree.current_nodes[0]

            self.assertFalse(
                deployment_tree.classification_monitor.ready_to_merge(first_node)
            )

            self.assertFalse(
                deployment_tree.classification_monitor.classification_performed(
                    first_node
                )
            )

            for leaf_index in first_node.leaves:

                leaf_node = deployment_tree.spawn_node_child_prepped(
                    first_node, leaf_index
                )

                registraction_success = deployment_tree.register_node_safe(leaf_node)
                self.assertTrue(registraction_success)

                self.assertTrue(
                    leaf_node.parameter_set.status == ParameterSet.STATUS_RUNNING
                )

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

            if (
                self.pipeline_makeup.check_makeuplist_has_classification(
                    makeup_explicit
                )
                is False
            ):
                continue

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

                run.set_to_queued()

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
                        os.path.basename(run.container.file_r1),
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
                run.container.run_main_prep_dump_tables()
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
