from django.test import TestCase
from django.contrib.auth.models import User
from pathogen_identification.models import Projects
from pathogen_identification.models import PIProject_Sample
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager
from pathogen_identification.constants_settings import Pipeline_Makeup
from settings.default_software import DefaultSoftware
import os, filecmp
from django.conf import settings
from settings.models import Software, Parameter, Sample
from settings.constants_settings import ConstantsSettings as CS
from pathogen_identification.models import SoftwareTree
from constants.constantsTestsCase import ConstantsTestsCase
from utils.software import Software as SoftwareUtils
from utils.utils import Utils
from constants.software_names import SoftwareNames
from constants.constants import Constants, TypePath, FileType, FileExtensions

# Create your tests here.


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


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

    def test_generate_trees(self):
        default_software = DefaultSoftware()
        default_software.test_all_defaults_pathogen_identification(self.test_user)
        utils_manager = Utils_Manager()
        utils_manager.generate_default_trees()

        all_trees = SoftwareTree.objects.all()
        self.assertEqual(all_trees.count(), len(self.pipeline_makeup.MAKEUP) * 2)


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
        file_name = os.path.join(
            self.baseDirectory,
            ConstantsTestsCase.DIR_FASTQ,
            ConstantsTestsCase.FASTQ_MINION_1,
        )
        self.assertTrue(file_name)

        try:
            user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
        except User.DoesNotExist:
            user = User()
            user.username = ConstantsTestsCase.TEST_USER_NAME
            user.is_active = False
            user.password = ConstantsTestsCase.TEST_USER_NAME
            user.save()

        self.test_user = user

        temp_dir = self.utils.get_temp_dir()
        self.utils.copy_file(
            file_name, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
        )

        sample_name = "televir_sample_minion_1"
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

        self.sample = sample_ont

        ######### ONT ##########
        project_ont_name = "project_televir"
        try:
            project_ont = Projects.objects.get(name=project_ont_name)
        except Projects.DoesNotExist:
            project_ont = Projects()
            project_ont.name = project_ont_name
            project_ont.owner = user
            project_ont.technology = CS.TECHNOLOGY_minion

            project_ont.save()

        self.project_ont = project_ont

        try:
            ont_project_sample = PIProject_Sample.objects.get(
                project__id=project_ont.pk, sample__id=sample_ont.pk
            )
        except PIProject_Sample.DoesNotExist:
            ont_project_sample = PIProject_Sample()
            ont_project_sample.project = project_ont
            ont_project_sample.sample = sample_ont
            ont_project_sample.save()

        ######## ILLUMINA #######
        project_illumina_name = "project_illumina"
        try:
            project_illumina = Projects.objects.get(name=project_illumina_name)
        except Projects.DoesNotExist:
            project_illumina = Projects()
            project_illumina.name = project_illumina_name
            project_illumina.owner = user
            project_illumina.technology = CS.TECHNOLOGY_illumina

            project_illumina.save()

        self.project_illumina = project_illumina

        ######
        default_software = DefaultSoftware()
        default_software.test_all_defaults_pathogen_identification(self.test_user)
        utils_manager = Utils_Manager()
        utils_manager.generate_default_trees()
