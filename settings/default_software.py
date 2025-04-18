"""
Created on 03/05/2020

@author: mmp
"""

from curses.ascii import SO

from django.contrib.auth.models import User

from constants.software_names import SoftwareNames
from pathogen_identification.constants_settings import ConstantsSettings as PICS
from pathogen_identification.utilities.utilities_pipeline import (
    Utility_Pipeline_Manager,
    Utils_Manager,
)
from settings.constants_settings import ConstantsSettings
from settings.default_parameters import DefaultParameters
from settings.models import Parameter, Software, SoftwareDefaultTest
from utils.lock_atomic_transaction import LockedAtomicTransaction


class DefaultSoftware(object):
    """
    classdocs
    """

    software_names = SoftwareNames()

    def __init__(self, test: int = 1):
        """change values"""
        self.default_parameters = DefaultParameters()
        self.televir_utiltity = Utility_Pipeline_Manager()
        # self.televir_utiltity.get_software_db_dict()

        self.change_values_software = {}  ### the key is the name of the software

    def test_televir_pipelines_available_once(self, user) -> bool:
        """
        test if exist, if not persist in database
        """
        try:
            software_test = SoftwareDefaultTest.objects.get(
                user=user, televir_pipelines_available=True
            )

            return True

        except SoftwareDefaultTest.DoesNotExist:
            utils = Utils_Manager()
            pipelines_available = utils.test_televir_pipelines_available(user)

            if pipelines_available:
                if not SoftwareDefaultTest.objects.filter(user=user).exists():
                    SoftwareDefaultTest.objects.create(
                        user=user, televir_pipelines_available=True
                    )
                else:
                    SoftwareDefaultTest.objects.filter(user=user).update(
                        televir_pipelines_available=True
                    )

            return pipelines_available

    def test_televir_software_available(self):
        """test if televir software is available"""
        user_system = User.objects.get(username="system")

        # self.test_all_defaults_pathogen_identification_once(user_system)

        televir_available = self.test_televir_pipelines_available_once(user_system)

        # self.remove_all_parameters(user_system)

        return televir_available

    def remove_all_parameters(self, user):
        """remove all parameters"""
        user_software = Software.objects.filter(owner=user)
        user_parameter = Parameter.objects.filter(software__in=user_software)
        with LockedAtomicTransaction(Parameter):
            user_parameter.delete()

        with LockedAtomicTransaction(Software):
            user_software.delete()

    def remove_all_software(self, user):
        """remove all software"""
        user_software = Software.objects.filter(owner=user)
        with LockedAtomicTransaction(Software):
            user_software.delete()

    def remove_all_televir_global_software(self, user):
        """remove all software"""
        user_software = Software.objects.filter(
            owner=user, type_of_use=Software.TYPE_OF_USE_televir_global
        )
        with LockedAtomicTransaction(Software):
            user_software.delete()

    def remove_all_televir_global_Parameters(self, user):
        """remove all software"""
        user_parameter = Parameter.objects.filter(
            software__type_of_use=Software.TYPE_OF_USE_televir_global,
            software__owner=user,
        )
        with LockedAtomicTransaction(Parameter):
            user_parameter.delete()

    def remove_all_televir_project_Parameters(self, user):
        user_project_parameter = Parameter.objects.filter(
            software__type_of_use=Software.TYPE_OF_USE_televir_project,
            software__owner=user,
        )
        with LockedAtomicTransaction(Parameter):
            user_project_parameter.delete()

    def remove_all_televir_project_software(self, user):
        """remove all software"""
        user_software = Software.objects.filter(
            owner=user, type_of_use=Software.TYPE_OF_USE_televir_project
        )
        with LockedAtomicTransaction(Software):
            user_software.delete()

    def remove_all_televir_global(self, user):
        """remove all software"""
        self.remove_all_televir_global_Parameters(user)
        self.remove_all_televir_global_software(user)

    def remove_all_televir_project(self, user):
        """remove all software"""
        self.remove_all_televir_project_Parameters(user)
        self.remove_all_televir_project_software(user)

    def remove_all_televir_software(self, user):
        """remove all software"""
        self.remove_all_televir_global(user)
        self.remove_all_televir_project(user)

    def reset_user_software(self, users):
        """
        for every user given, if has been active, remove all televir software and add the default televir software.
        Including the default parameters.
        Including Project software and Project parameters -> everything is reset.
        """

        for user in users:
            global_software = Software.objects.filter(
                type_of_use=Software.TYPE_OF_USE_televir_global
            )
            if global_software.exists():
                self.remove_all_televir_software(user)
                self.test_all_defaults_pathogen_identification(user)

    def test_all_defaults_once(self, user: User):

        try:
            SoftwareDefaultTest.objects.get(user=user, is_tested_all_defaults=True)
        except SoftwareDefaultTest.DoesNotExist:
            self.test_all_defaults(user)
            if not SoftwareDefaultTest.objects.filter(user=user).exists():
                SoftwareDefaultTest.objects.create(
                    user=user, is_tested_all_defaults=True
                )
            else:
                SoftwareDefaultTest.objects.filter(user=user).update(
                    is_tested_all_defaults=True
                )

    def test_all_defaults_pathogen_identification_once(self, user: User):

        try:
            SoftwareDefaultTest.objects.get(user=user, is_tested_televir_defaults=True)

        except SoftwareDefaultTest.DoesNotExist:
            self.test_all_defaults_pathogen_identification(user)
            if not SoftwareDefaultTest.objects.filter(user=user).exists():
                SoftwareDefaultTest.objects.create(
                    user=user, is_tested_televir_defaults=True
                )
            else:
                SoftwareDefaultTest.objects.filter(user=user).update(
                    is_tested_televir_defaults=True
                )

    def test_all_defaults(self, user: User):
        ### test all defaults
        self.test_default_db(
            SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
            self.default_parameters.get_trimmomatic_default(
                user, Software.TYPE_OF_USE_qc, ConstantsSettings.TECHNOLOGY_illumina
            ),
            user,
        )
        self.test_default_db(
            SoftwareNames.SOFTWARE_SNIPPY_name,
            self.default_parameters.get_snippy_default(
                user, Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_illumina
            ),
            user,
            SoftwareNames.SOFTWARE_SNIPPY_name_extended,
        )
        self.test_default_db(
            SoftwareNames.SOFTWARE_IVAR_name,
            self.default_parameters.get_ivar_default(
                user, Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_illumina
            ),
            user,
            SoftwareNames.SOFTWARE_IVAR_name_extended,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_IRMA_name,
            self.default_parameters.get_irma_default(
                user,
                Software.TYPE_OF_USE_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
            SoftwareNames.SOFTWARE_IRMA_name_extended,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_FREEBAYES_name,
            self.default_parameters.get_freebayes_default(
                user, Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_illumina
            ),
            user,
        )
        self.test_default_db(
            SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
            self.default_parameters.get_mask_consensus_threshold_default(
                user, Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_illumina
            ),
            user,
        )
        self.test_default_db(
            SoftwareNames.SOFTWARE_ABRICATE_name,
            self.default_parameters.get_abricate_default(
                user, Software.TYPE_OF_USE_qc, ConstantsSettings.TECHNOLOGY_illumina
            ),
            user,
        )
        #         self.test_default_db(SoftwareNames.SOFTWARE_CLEAN_HUMAN_READS_name,
        #                 self.default_parameters.get_clean_human_reads_default(user,
        #                 Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_illumina), user)

        ## ONT software
        self.test_default_db(
            SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
            self.default_parameters.get_mask_consensus_threshold_default(
                user, Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_minion
            ),
            user,
        )
        self.test_default_db(
            SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,
            self.default_parameters.get_limit_coverage_ONT_threshold_default(
                user, Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_minion
            ),
            user,
        )
        self.test_default_db(
            SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name,
            self.default_parameters.get_vcf_freq_ONT_threshold_default(
                user, Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_minion
            ),
            user,
        )
        self.test_default_db(
            SoftwareNames.SOFTWARE_Medaka_name_consensus,
            self.default_parameters.get_medaka_model_default(
                user, Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_minion
            ),
            user,
        )

        #         self.test_default_db(SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT,
        #                 self.default_parameters.get_samtools_depth_default_ONT(user,
        #                 Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_minion), user)

        self.test_default_db(
            SoftwareNames.SOFTWARE_NanoFilt_name,
            self.default_parameters.get_nanofilt_default(
                user, Software.TYPE_OF_USE_qc, ConstantsSettings.TECHNOLOGY_minion
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_ABRICATE_name,
            self.default_parameters.get_abricate_default(
                user, Software.TYPE_OF_USE_qc, ConstantsSettings.TECHNOLOGY_minion
            ),
            user,
        )

        #         self.test_default_db(SoftwareNames.SOFTWARE_CLEAN_HUMAN_READS_name,
        #                 self.default_parameters.get_clean_human_reads_default(user,
        #                 Software.TYPE_OF_USE_global, ConstantsSettings.TECHNOLOGY_minion), user)

        #############
        ############# PATHOGEN IDENTIFICATION SOFTWARE
        #############

        self.test_all_defaults_pathogen_identification(user)

    def test_all_defaults_pathogen_identification(self, user):
        self.test_default_db(
            SoftwareNames.SOFTWARE_REMAP_PARAMS_name,
            self.default_parameters.get_remap_defaults(
                user,
                Software.TYPE_OF_USE_televir_settings,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_REMAP_PARAMS_name,
            self.default_parameters.get_remap_defaults(
                user,
                Software.TYPE_OF_USE_televir_settings,
                ConstantsSettings.TECHNOLOGY_minion,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_BAMUTIL_name,
            self.default_parameters.get_bamutil_defaults(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_BAMUTIL_name,
            self.default_parameters.get_bamutil_defaults(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_MSAMTOOLS_name,
            self.default_parameters.get_msamtools_defaults(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_MSAMTOOLS_name,
            self.default_parameters.get_msamtools_defaults(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_map_filtering,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_MSAMTOOLS_name,
            self.default_parameters.get_msamtools_defaults(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_MSAMTOOLS_name,
            self.default_parameters.get_msamtools_defaults(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_map_filtering,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_DUSTMASKER_name,
            self.default_parameters.get_dustmasker_defaults(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_DUSTMASKER_name,
            self.default_parameters.get_dustmasker_defaults(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_remap_filtering,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_DUSTMASKER_name,
            self.default_parameters.get_dustmasker_defaults(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_DUSTMASKER_name,
            self.default_parameters.get_dustmasker_defaults(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_remap_filtering,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_PRINSEQ_name,
            self.default_parameters.get_prinseq_defaults(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_PRINSEQ_name,
            self.default_parameters.get_prinseq_defaults(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_televir_report_layout_name,
            self.default_parameters.get_televir_report_defaults(
                user,
                Software.TYPE_OF_USE_televir_settings,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_televir_report_layout_name,
            self.default_parameters.get_televir_report_defaults(
                user,
                Software.TYPE_OF_USE_televir_settings,
                ConstantsSettings.TECHNOLOGY_minion,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_CENTRIFUGE_name,
            self.default_parameters.get_centrifuge_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
        )
        self.test_default_db(
            SoftwareNames.SOFTWARE_CENTRIFUGE_name,
            self.default_parameters.get_centrifuge_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_read_classification,
                is_to_run=False,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_CENTRIFUGE_name,
            self.default_parameters.get_centrifuge_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_read_classification,
                is_to_run=False,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_CENTRIFUGE_name,
            self.default_parameters.get_centrifuge_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_BWA_name,
            self.default_parameters.get_bwa_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_KRAKEN2_name,
            self.default_parameters.get_kraken2_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_KAIJU_name,
            self.default_parameters.get_kaiju_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_read_classification,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_DIAMOND_name,
            self.default_parameters.get_diamond_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_DIAMOND_name,
            self.default_parameters.get_diamond_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_KRAKENUNIQ_name,
            self.default_parameters.get_krakenuniq_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_KRAKENUNIQ_name,
            self.default_parameters.get_krakenuniq_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_BLAST_name,
            self.default_parameters.get_blast_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_BLAST_name,
            self.default_parameters.get_blast_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_FASTVIROMEEXPLORER_name,
            self.default_parameters.get_fastviromeexplorer_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            ),
            user,
        )

        ########################### hslib19 in centos 7 is too old for this
        ########################### uncomment in ubuntu
        # self.test_default_db(
        #    SoftwareNames.SOFTWARE_DESAMBA_name,
        #    self.default_parameters.get_desamba_default(
        #        user,
        #        Software.TYPE_OF_USE_televir_global,
        #        ConstantsSettings.TECHNOLOGY_minion,
        #    ),
        #    user,
        # )

        self.test_default_db(
            SoftwareNames.SOFTWARE_SPAdes_name,
            self.default_parameters.get_spades_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_RAVEN_name,
            self.default_parameters.get_raven_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_SNIPPY_PI_name,
            self.default_parameters.get_snippy_pi_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_SNIPPY_PI_name,
            self.default_parameters.get_snippy_pi_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_request_mapping,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ONT_name,
            self.default_parameters.get_minimap2_remap_ONT_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_MINIMAP2_DEPLETE_ONT_name,
            self.default_parameters.get_minimap2_depletion_ONT_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_BOWTIE2_DEPLETE_name,
            self.default_parameters.get_bowtie2_deplete_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_KRAKEN2_name,
            self.default_parameters.get_kraken2_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_contig_classification,
                is_to_run=False,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_KRAKEN2_name,
            self.default_parameters.get_kraken2_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_read_classification,
                is_to_run=False,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_KRAKEN2_name,
            self.default_parameters.get_kraken2_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_contig_classification,
                is_to_run=False,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_KRAKEN2_name,
            self.default_parameters.get_kraken2_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_viral_enrichment,
                is_to_run=False,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_BOWTIE2_REMAP_name,
            self.default_parameters.get_bowtie2_remap_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
        )

        if PICS.TEST_SOFTWARE:
            self.test_defaults_test_televir(user)
        if PICS.METAGENOMICS:
            self.test_defaults_metagenomics(user)

    def test_defaults_metagenomics(self, user):
        """
        test if exist, if not persist in database, for metagenomics"""

        self.test_default_db(
            SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ONT_name,
            self.default_parameters.get_minimap2_remap_ONT_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_metagenomics_screening,
                is_to_run=True,
                job="screening",
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ONT_name,
            self.default_parameters.get_minimap2_remap_ONT_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_request_mapping,
                is_to_run=True,
                job="request_mapping",
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ILLU_name,
            self.default_parameters.get_minimap2_remap_illumina_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_metagenomics_screening,
                is_to_run=True,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_BOWTIE2_REMAP_name,
            self.default_parameters.get_bowtie2_remap_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_request_mapping,
                job="request_mapping",
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_METAGENOMICS_SETTINGS_name,
            self.default_parameters.get_metagenomics_settings_defaults(
                user,
                Software.TYPE_OF_USE_televir_settings,
                ConstantsSettings.TECHNOLOGY_illumina,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_METAGENOMICS_SETTINGS_name,
            self.default_parameters.get_metagenomics_settings_defaults(
                user,
                Software.TYPE_OF_USE_televir_settings,
                ConstantsSettings.TECHNOLOGY_minion,
            ),
            user,
        )

    def test_defaults_test_televir(self, user):
        """
        test if exist, if not persist in database, for televir"""

        self.test_default_db(
            SoftwareNames.SOFTWARE_BWA_name,
            self.default_parameters.get_bwa_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_request_mapping,
            ),
            user,
        )

        self.test_default_db(
            SoftwareNames.SOFTWARE_BWA_name,
            self.default_parameters.get_bwa_default(
                user,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_remapping,
            ),
            user,
        )

        #        self.test_default_db(
        #            SoftwareNames.SOFTWARE_MINIMAP2_MAP_ASSEMBLY_name,
        #            self.default_parameters.get_minimap2_map_assembly_default(
        #                user,
        #                Software.TYPE_OF_USE_televir_global,
        #                ConstantsSettings.TECHNOLOGY_minion,
        #            ),
        #            user,
        #        )
        #
        #        self.test_default_db(
        #            SoftwareNames.SOFTWARE_MINIMAP2_MAP_ASSEMBLY_name,
        #            self.default_parameters.get_minimap2_map_assembly_default(
        #                user,
        #                Software.TYPE_OF_USE_televir_global,
        #                ConstantsSettings.TECHNOLOGY_illumina,
        #            ),
        #            user,
        #        )

        # self.test_default_db(
        #    SoftwareNames.SOFTWARE_MINIMAP2_DEPLETE_ILLU_name,
        #    self.default_parameters.get_minimap2_depletion_illumina_default(
        #        user,
        #        Software.TYPE_OF_USE_televir_global,
        #        ConstantsSettings.TECHNOLOGY_illumina,
        #    ),
        #    user,
        # )

    def assess_db_dependency_met(self, vect_parameters, software_name):
        """for pipeline steps where sequence dbs are required, check that they exist."""
        if (
            vect_parameters[0].software.type_of_use
            == Software.TYPE_OF_USE_televir_global
        ):
            if (
                vect_parameters[0].software.pipeline_step.name
                in self.televir_utiltity.steps_db_dependant
            ):
                if not self.televir_utiltity.check_software_db_available(
                    software_name=software_name,
                ):
                    return False

        return True

    def test_default_db(self, software_name, vect_parameters, user, name_extended=None):
        """
        test if exist, if not persist in database
        """
        type_of_use = Software.TYPE_OF_USE_global

        if not self.assess_db_dependency_met(vect_parameters, software_name):
            return

        ## lock because more than one process can duplicate software names

        try:
            type_of_use = vect_parameters[0].software.type_of_use
        except:
            pass

        if name_extended is None:
            self.test_default_persist_general(
                software_name, vect_parameters, user, type_of_use
            )
        else:
            self.test_default_persist_specific(
                software_name, vect_parameters, user, name_extended, type_of_use
            )

    def test_default_persist_general(
        self, software_name, vect_parameters, user, type_of_use
    ):
        try:

            software_queried = Software.objects.get(
                name=software_name,
                owner=user,
                type_of_use=vect_parameters[0].software.type_of_use,
                technology__name=vect_parameters[0].software.technology.name,
                version_parameters=self.default_parameters.get_software_parameters_version(
                    software_name
                ),
                pipeline_step__name=vect_parameters[0].software.pipeline_step,
            )

        except Software.MultipleObjectsReturned:
            ## keep the first one, delete the rest
            software_query = Software.objects.filter(
                name=software_name,
                owner=user,
                type_of_use=vect_parameters[0].software.type_of_use,
                technology__name=vect_parameters[0].software.technology.name,
                version_parameters=self.default_parameters.get_software_parameters_version(
                    software_name
                ),
                pipeline_step__name=vect_parameters[0].software.pipeline_step,
                parameter__televir_project=None,
                parameter__televir_project_sample=None,
            ).order_by("id")

            if software_query.count() > 1:
                software = software_query.exclude(pk=software_query.last().pk)

                parameters = Parameter.objects.filter(software__in=software)
                with LockedAtomicTransaction(Parameter):
                    parameters.delete()
                with LockedAtomicTransaction(Software):
                    software.delete()

        except Software.DoesNotExist:  ### if not exist save it
            self.default_parameters.persist_parameters(vect_parameters, type_of_use)

    def test_default_persist_specific(
        self, software_name, vect_parameters, user, name_extended, type_of_use
    ):

        try:

            software_queried = Software.objects.get(
                name=software_name,
                owner=user,
                type_of_use=vect_parameters[0].software.type_of_use,
                technology__name=vect_parameters[0].software.technology.name,
                version_parameters=self.default_parameters.get_software_parameters_version(
                    software_name
                ),
                pipeline_step__name=vect_parameters[0].software.pipeline_step,
                name_extended=name_extended,
            )

        except Software.MultipleObjectsReturned:
            ## keep the first one, delete the rest
            software_query = Software.objects.filter(
                name=software_name,
                owner=user,
                type_of_use=vect_parameters[0].software.type_of_use,
                technology__name=vect_parameters[0].software.technology.name,
                version_parameters=self.default_parameters.get_software_parameters_version(
                    software_name
                ),
                pipeline_step__name=vect_parameters[0].software.pipeline_step,
                name_extended=name_extended,
                parameter__televir_project=None,
                parameter__televir_project_sample=None,
            ).order_by("id")

            if software_query.count() > 1:
                software = software_query.exclude(pk=software_query.last().pk)

                parameters = Parameter.objects.filter(software__in=software)
                with LockedAtomicTransaction(Parameter):
                    parameters.delete()
                with LockedAtomicTransaction(Software):
                    software.delete()

        except Software.DoesNotExist:
            self.default_parameters.persist_parameters(vect_parameters, type_of_use)

    def get_trimmomatic_parameters(self, user):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
            user,
            Software.TYPE_OF_USE_qc,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_illumina,
        )
        return "" if result is None else result

    def get_snippy_parameters(self, user, is_to_run=False):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_SNIPPY_name,
            user,
            Software.TYPE_OF_USE_global,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_illumina,
            is_to_run=is_to_run,
            software_name_extended=SoftwareNames.SOFTWARE_SNIPPY_name_extended,
        )
        return "" if result is None else result

    def get_ivar_parameters(self, user):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_IVAR_name,
            user,
            Software.TYPE_OF_USE_global,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_illumina,
            software_name_extended=SoftwareNames.SOFTWARE_IVAR_name_extended,
        )
        return "" if result is None else result

    def get_irma_parameters(self, user):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_IRMA_name,
            user,
            Software.TYPE_OF_USE_global,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_illumina,
            software_name_extended=SoftwareNames.SOFTWARE_IRMA_name_extended,
        )

        return "" if result is None else result

    def get_freebayes_parameters(self, user):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_FREEBAYES_name,
            user,
            Software.TYPE_OF_USE_global,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_illumina,
        )
        return "" if result is None else result

    def get_nanofilt_parameters(self, user):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_NanoFilt_name,
            user,
            Software.TYPE_OF_USE_qc,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        return "" if result is None else result

    def get_mask_consensus_threshold_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
            user,
            Software.TYPE_OF_USE_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_clean_human_reads_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_CLEAN_HUMAN_READS_name,
            user,
            Software.TYPE_OF_USE_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_limit_coverage_ONT_parameters(self, user):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,
            user,
            Software.TYPE_OF_USE_global,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        return "" if result is None else result

    def get_vcf_freq_ONT_parameters(self, user):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name,
            user,
            Software.TYPE_OF_USE_global,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        return "" if result is None else result

    def get_medaka_parameters_consensus(self, user):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_Medaka_name_consensus,
            user,
            Software.TYPE_OF_USE_global,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        return "" if result is None else result

    def get_samtools_parameters_depth_ONT(self, user):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT,
            user,
            Software.TYPE_OF_USE_global,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        return "" if result is None else result

    def get_abricate_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_ABRICATE_name,
            user,
            Software.TYPE_OF_USE_qc,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    ###
    ### PATHOGEN DETECTION PARAMETERS

    def get_remap_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_REMAP_PARAMS_name,
            user,
            Software.TYPE_OF_USE_televir_settings,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_prinseq_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_PRINSEQ_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_bamutil_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_BAMUTIL_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_dustmasker_parameters(self, user, technology_name, pipeline_step=None):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_DUSTMASKER_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
            pipeline_step=pipeline_step,
        )
        return "" if result is None else result

    def get_metagenomics_settings_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_METAGENOMICS_SETTINGS_name,
            user,
            Software.TYPE_OF_USE_televir_settings,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_msamtools_parameters(self, user, technology_name, pipeline_step=None):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_MSAMTOOLS_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
            pipeline_step=pipeline_step,
        )
        return "" if result is None else result

    def get_televir_report_layout_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_televir_report_layout_name,
            user,
            Software.TYPE_OF_USE_televir_settings,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_kaiju_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_KAIJU_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_centrifuge_parameters(self, user, technology_name, pipeline_step=None):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_CENTRIFUGE_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
            pipeline_step=pipeline_step,
        )
        return "" if result is None else result

    def get_bwa_parameters(self, user, technology_name, pipeline_step=None):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_BWA_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
            pipeline_step=pipeline_step,
        )
        return "" if result is None else result

    def get_kraken2_parameters(self, user, technology_name, pipeline_step=None):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_KRAKEN2_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name=technology_name,
            pipeline_step=pipeline_step,
        )
        return "" if result is None else result

    def get_diamond_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_DIAMOND_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_blast_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_BLAST_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_krakenuniq_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_KRAKENUNIQ_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_desamba_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_DESAMBA_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_raven_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_RAVEN_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_spades_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_SPAdes_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_fastviromeexplorer_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_FASTVIROMEEXPLORER_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_snippy_pi_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_SNIPPY_PI_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_minimap2_remap_ont_parameters(self, user):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ONT_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        return "" if result is None else result

    def get_bowtie2_deplete_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_BOWTIE2_DEPLETE_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_bowtie2_remap_parameters(self, user, technology_name, pipeline_step=None):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_BOWTIE2_REMAP_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
            pipeline_step=pipeline_step,
        )
        return "" if result is None else result

    def get_minimap2_deplete_ont_parameters(self, user):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_MINIMAP2_DEPLETE_ONT_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        return "" if result is None else result

    def get_minimap2_remap_illumina_parameters(self, user):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ILLU_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_illumina,
        )
        return "" if result is None else result

    def get_minimap2_map_assembly_parameters(self, user):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_MINIMAP2_MAP_ASSEMBLY_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        return "" if result is None else result

    def get_minimap2_deplete_illumina(self, user):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_MINIMAP2_DEPLETE_ILLU_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_illumina,
        )
        return "" if result is None else result

    def get_nextstrain_parameters(self, user):
        result = self.default_parameters.get_parameters_parsed(
            SoftwareNames.SOFTWARE_NEXTSTRAIN_name,
            user,
            Software.TYPE_OF_USE_global,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_generic,
        )
        return "" if result is None else result

    ####
    ####
    def set_default_software(self, software):
        """Set a default"""
        vect_parameters = self.default_parameters.get_vect_parameters(software)
        if vect_parameters is None:
            return

        key_value = "{}_{}".format(software.name, software.technology.name)
        parameters = Parameter.objects.filter(software=software)
        self.change_values_software[key_value] = False
        for parameter in parameters:
            if parameter.can_change:
                for parameter_to_set_default in vect_parameters:
                    if parameter_to_set_default.sequence_out == parameter.sequence_out:
                        ###   if change software name
                        if parameter.parameter != parameter_to_set_default.parameter:
                            self.change_values_software[key_value] = True
                            parameter.parameter = parameter_to_set_default.parameter
                            parameter.save()
                        break

    def is_change_values_for_software(self, software, user_name):
        """Return if the software has a value changed"""
        key_value = "{}_{}_{}".format(
            software.name, software.technology.name, user_name
        )
        return self.change_values_software.get(key_value, False)

    def get_parameters(
        self,
        software_name,
        user,
        technology_name=ConstantsSettings.TECHNOLOGY_illumina,
        pipeline_step=None,
    ):
        """
        Return the parameters for a software

        :param software_name: software name
        :param user: user
        :param technology_name: technology name
        :return: parameters

        """
        if software_name == SoftwareNames.SOFTWARE_TRIMMOMATIC_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
                self.default_parameters.get_trimmomatic_default(
                    user, Software.TYPE_OF_USE_qc, ConstantsSettings.TECHNOLOGY_illumina
                ),
                user,
            )
            return self.get_trimmomatic_parameters(user)

        if software_name == SoftwareNames.SOFTWARE_IVAR_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_IVAR_name,
                self.default_parameters.get_ivar_default(
                    user,
                    Software.TYPE_OF_USE_global,
                    ConstantsSettings.TECHNOLOGY_illumina,
                    pipeline_step=ConstantsSettings.PIPELINE_NAME_variant_detection,
                ),
                user,
                SoftwareNames.SOFTWARE_IVAR_name_extended,
            )
            return self.get_ivar_parameters(user)

        if software_name == SoftwareNames.SOFTWARE_IRMA_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_IRMA_name,
                self.default_parameters.get_irma_default(
                    user,
                    Software.TYPE_OF_USE_global,
                    ConstantsSettings.TECHNOLOGY_illumina,
                ),
                user,
                SoftwareNames.SOFTWARE_IRMA_name_extended,
            )
            return self.get_irma_parameters(user)

        if software_name == SoftwareNames.SOFTWARE_SNIPPY_name:

            self.test_default_db(
                SoftwareNames.SOFTWARE_SNIPPY_name,
                self.default_parameters.get_snippy_default(
                    user,
                    Software.TYPE_OF_USE_global,
                    ConstantsSettings.TECHNOLOGY_illumina,
                    pipeline_step=ConstantsSettings.PIPELINE_NAME_variant_detection,
                ),
                user,
                SoftwareNames.SOFTWARE_SNIPPY_name_extended,
            )

            return self.get_snippy_parameters(user)

        if software_name == SoftwareNames.SOFTWARE_FREEBAYES_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_FREEBAYES_name,
                self.default_parameters.get_freebayes_default(
                    user,
                    Software.TYPE_OF_USE_global,
                    ConstantsSettings.TECHNOLOGY_illumina,
                ),
                user,
            )

            return self.get_freebayes_parameters(user)
        if software_name == SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name:
            self.test_default_db(
                SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
                self.default_parameters.get_mask_consensus_threshold_default(
                    user, Software.TYPE_OF_USE_global, technology_name
                ),
                user,
            )
            return self.get_mask_consensus_threshold_parameters(user, technology_name)
        if software_name == SoftwareNames.SOFTWARE_CLEAN_HUMAN_READS_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_CLEAN_HUMAN_READS_name,
                self.default_parameters.get_clean_human_reads_default(
                    user, Software.TYPE_OF_USE_global, technology_name
                ),
                user,
            )
            return self.get_clean_human_reads_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_NanoFilt_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_NanoFilt_name,
                self.default_parameters.get_nanofilt_default(
                    user,
                    Software.TYPE_OF_USE_qc,
                    ConstantsSettings.TECHNOLOGY_minion,
                ),
                user,
            )
            return self.get_nanofilt_parameters(user)
        if software_name == SoftwareNames.SOFTWARE_Medaka_name_consensus:
            self.test_default_db(
                SoftwareNames.SOFTWARE_Medaka_name_consensus,
                self.default_parameters.get_medaka_model_default(
                    user,
                    Software.TYPE_OF_USE_global,
                    ConstantsSettings.TECHNOLOGY_minion,
                ),
                user,
            )
            return self.get_medaka_parameters_consensus(user)
        if software_name == SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT:
            self.test_default_db(
                software_name,
                self.default_parameters.get_samtools_depth_default_ONT(
                    user,
                    Software.TYPE_OF_USE_global,
                    ConstantsSettings.TECHNOLOGY_minion,
                ),
                user,
            )
            return self.get_samtools_parameters_depth_ONT(user)
        if software_name == SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name:
            self.test_default_db(
                SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,
                self.default_parameters.get_limit_coverage_ONT_threshold_default(
                    user,
                    Software.TYPE_OF_USE_global,
                    ConstantsSettings.TECHNOLOGY_minion,
                ),
                user,
            )
            return self.get_limit_coverage_ONT_parameters(user)
        if software_name == SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name:
            self.test_default_db(
                SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name,
                self.default_parameters.get_vcf_freq_ONT_threshold_default(
                    user,
                    Software.TYPE_OF_USE_global,
                    ConstantsSettings.TECHNOLOGY_minion,
                ),
                user,
            )
            return self.get_vcf_freq_ONT_parameters(user)
        if software_name == SoftwareNames.SOFTWARE_ABRICATE_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_ABRICATE_name,
                self.default_parameters.get_abricate_default(
                    user, Software.TYPE_OF_USE_qc, technology_name
                ),
                user,
            )
            return self.get_abricate_parameters(user, technology_name)
        if software_name == SoftwareNames.SOFTWARE_NEXTSTRAIN_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_NEXTSTRAIN_name,
                self.default_parameters.get_nextstrain_default(user),
                user,
            )
            return self.get_nextstrain_parameters(user)

        ################################################
        ############### TELEVIR SOFTWARE ###############

        if software_name == SoftwareNames.SOFTWARE_REMAP_PARAMS_name:

            return self.get_remap_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_PRINSEQ_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_PRINSEQ_name,
                self.default_parameters.get_prinseq_defaults(
                    user,
                    Software.TYPE_OF_USE_televir_global,
                    ConstantsSettings.TECHNOLOGY_illumina,
                ),
                user,
            )
            return self.get_prinseq_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_BAMUTIL_name:

            return self.get_bamutil_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_DUSTMASKER_name:

            return self.get_dustmasker_parameters(
                user, technology_name, pipeline_step=pipeline_step
            )

        if software_name == SoftwareNames.SOFTWARE_METAGENOMICS_SETTINGS_name:

            return self.get_metagenomics_settings_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_MSAMTOOLS_name:

            return self.get_msamtools_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_televir_report_layout_name:

            return self.get_televir_report_layout_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_CENTRIFUGE_name:

            return self.get_centrifuge_parameters(
                user, technology_name, pipeline_step=pipeline_step
            )

        if software_name == SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ONT_name:

            return self.get_minimap2_remap_ont_parameters(user)

        if software_name == SoftwareNames.SOFTWARE_BOWTIE2_DEPLETE_name:

            return self.get_bowtie2_deplete_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_BOWTIE2_REMAP_name:

            return self.get_bowtie2_remap_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_MINIMAP2_DEPLETE_ONT_name:

            return self.get_minimap2_deplete_ont_parameters(user)

        if software_name == SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ILLU_name:

            return self.get_minimap2_remap_illumina_parameters(user)

        if software_name == SoftwareNames.SOFTWARE_MINIMAP2_MAP_ASSEMBLY_name:

            return self.get_minimap2_map_assembly_parameters(user)

        if software_name == SoftwareNames.SOFTWARE_MINIMAP2_DEPLETE_ILLU_name:

            return self.get_minimap2_deplete_illumina(user)

        if software_name == SoftwareNames.SOFTWARE_KRAKEN2_name:

            return self.get_kraken2_parameters(
                user, technology_name, pipeline_step=pipeline_step
            )

        if software_name == SoftwareNames.SOFTWARE_BWA_name:

            return self.get_bwa_parameters(
                user, technology_name, pipeline_step=pipeline_step
            )

        if software_name == SoftwareNames.SOFTWARE_DIAMOND_name:

            return self.get_diamond_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_KRAKENUNIQ_name:

            return self.get_krakenuniq_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_KAIJU_name:

            return self.get_kaiju_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_BLAST_name:

            return self.get_blast_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_DESAMBA_name:

            return self.get_desamba_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_RAVEN_name:

            return self.get_raven_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_FASTVIROMEEXPLORER_name:

            return self.get_fastviromeexplorer_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_SPAdes_name:

            return self.get_spades_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_SNIPPY_PI_name:

            return self.get_snippy_pi_parameters(user, technology_name)

        return ""

    def get_all_software(self):
        """
        get all softwares available by this class
        """
        vect_software = []
        vect_software.append(self.software_names.get_trimmomatic_name())
        vect_software.append(self.software_names.get_snippy_name())
        vect_software.append(self.software_names.get_ivar_name())
        vect_software.append(self.software_names.get_irma_name())
        vect_software.append(self.software_names.get_freebayes_name())
        vect_software.append(self.software_names.get_NanoFilt_name())
        vect_software.append(
            self.software_names.get_insaflu_parameter_mask_consensus_name()
        )
        vect_software.append(self.software_names.get_medaka_name_extended_consensus())
        vect_software.append(self.software_names.get_samtools_name_depth_ONT())
        vect_software.append(
            self.software_names.get_insaflu_parameter_limit_coverage_name()
        )
        vect_software.append(self.software_names.get_insaflu_parameter_freq_vcf_name())
        vect_software.append(self.software_names.get_abricate_name())
        return vect_software

    def is_software_to_run(self, software_name, user, technology_name):
        """Test if it is necessary to run this software, By default return True"""
        return self.default_parameters.is_software_to_run(
            software_name,
            user,
            Software.TYPE_OF_USE_global,
            None,
            None,
            None,
            technology_name,
        )

    def set_software_to_run(self, software_name, user, technology_name, is_to_run):
        """Set True/False on software
        :output True if the is_to_run is changed"""
        return self.default_parameters.set_software_to_run(
            software_name,
            user,
            Software.TYPE_OF_USE_global,
            None,
            None,
            None,
            technology_name,
            is_to_run,
        )
