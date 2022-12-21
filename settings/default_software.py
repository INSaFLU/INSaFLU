"""
Created on 03/05/2020

@author: mmp
"""
from curses.ascii import SO

from constants.software_names import SoftwareNames
from pathogen_identification.utilities.utilities_pipeline import (
    Utility_Pipeline_Manager,
)
from utils.lock_atomic_transaction import LockedAtomicTransaction

from settings.constants_settings import ConstantsSettings
from settings.default_parameters import DefaultParameters
from settings.models import Parameter, Software


class DefaultSoftware(object):
    """
    classdocs
    """

    software_names = SoftwareNames()

    def __init__(self):
        """change values"""
        self.default_parameters = DefaultParameters()
        self.televir_utiltity = Utility_Pipeline_Manager()
        self.change_values_software = {}  ### the key is the name of the software

    def remove_all_parameters(self, user):
        """remove all parameters"""
        user_software = Software.objects.filter(owner=user)
        user_parameter = Parameter.objects.filter(software__in=user_software)
        with LockedAtomicTransaction(Parameter):
            user_parameter.delete()

        with LockedAtomicTransaction(Software):
            user_software.delete()

    def test_all_defaults(self, user):

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

    def test_default_db(self, software_name, vect_parameters, user):
        """
        test if exist, if not persist in database
        """
        type_of_use = Software.TYPE_OF_USE_global

        if not self.assess_db_dependency_met(vect_parameters, software_name):
            return

        ## lock because more than one process can duplicate software names
        with LockedAtomicTransaction(Software), LockedAtomicTransaction(Parameter):

            try:
                Software.objects.get(
                    name=software_name,
                    owner=user,
                    type_of_use=vect_parameters[0].software.type_of_use,
                    technology__name=vect_parameters[0].software.technology.name,
                    version_parameters=self.default_parameters.get_software_parameters_version(
                        software_name
                    ),
                    pipeline_step__name=vect_parameters[0].software.pipeline_step,
                )
            except Software.DoesNotExist:  ### if not exist save it
                self.default_parameters.persist_parameters(vect_parameters, type_of_use)

    def get_trimmomatic_parameters(self, user):
        result = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
            user,
            Software.TYPE_OF_USE_qc,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_illumina,
        )
        return "" if result is None else result

    def get_snippy_parameters(self, user):
        result = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_SNIPPY_name,
            user,
            Software.TYPE_OF_USE_global,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_illumina,
        )
        return "" if result is None else result

    def get_freebayes_parameters(self, user):
        result = self.default_parameters.get_parameters(
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
        result = self.default_parameters.get_parameters(
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
        result = self.default_parameters.get_parameters(
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
        result = self.default_parameters.get_parameters(
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
        result = self.default_parameters.get_parameters(
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
        result = self.default_parameters.get_parameters(
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
        result = self.default_parameters.get_parameters(
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
        result = self.default_parameters.get_parameters(
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
        result = self.default_parameters.get_parameters(
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

    def get_kaiju_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_KAIJU_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_centrifuge_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_CENTRIFUGE_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_bwa_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_BWA_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_kraken2_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_KRAKEN2_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_diamond_parameters(self, user, technology_name):
        result = self.default_parameters.get_parameters(
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
        result = self.default_parameters.get_parameters(
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
        result = self.default_parameters.get_parameters(
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
        result = self.default_parameters.get_parameters(
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
        result = self.default_parameters.get_parameters(
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
        result = self.default_parameters.get_parameters(
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
        result = self.default_parameters.get_parameters(
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
        result = self.default_parameters.get_parameters(
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
        result = self.default_parameters.get_parameters(
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
        result = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_BOWTIE2_DEPLETE_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            technology_name,
        )
        return "" if result is None else result

    def get_minimap2_deplete_ont_parameters(self, user):
        result = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_MINIMAP2_DEPLETE_ONT_name,
            user,
            Software.TYPE_OF_USE_televir_global,
            None,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        return "" if result is None else result

    def get_nextstrain_parameters(self, user):
        result = self.default_parameters.get_parameters(
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
        self, software_name, user, technology_name=ConstantsSettings.TECHNOLOGY_illumina
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
            print("get_abricate_parameters")
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

        ##########################################
        ############### TELEVIR SOFTWARE #########

        if software_name == SoftwareNames.SOFTWARE_CENTRIFUGE_name:
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
            return self.get_centrifuge_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ONT_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ONT_name,
                self.default_parameters.get_minimap2_remap_ONT_default(
                    user,
                    Software.TYPE_OF_USE_televir_global,
                    ConstantsSettings.TECHNOLOGY_minion,
                    pipeline_step=ConstantsSettings.PIPELINE_NAME_remapping,
                ),
                user,
            )

            return self.get_minimap2_remap_ont_parameters(user)

        if software_name == SoftwareNames.SOFTWARE_BOWTIE2_DEPLETE_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_BOWTIE2_DEPLETE_name,
                self.default_parameters.get_bowtie2_deplete_default(
                    user,
                    Software.TYPE_OF_USE_televir_global,
                    ConstantsSettings.TECHNOLOGY_illumina,
                ),
                user,
            )
            return self.get_bowtie2_deplete_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_MINIMAP2_DEPLETE_ONT_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_MINIMAP2_DEPLETE_ONT_name,
                self.default_parameters.get_minimap2_depletion_ONT_default(
                    user,
                    Software.TYPE_OF_USE_televir_global,
                    ConstantsSettings.TECHNOLOGY_minion,
                    pipeline_step=ConstantsSettings.PIPELINE_NAME_host_depletion,
                ),
                user,
            )

            return self.get_minimap2_deplete_ont_parameters(user)

        if software_name == SoftwareNames.SOFTWARE_KRAKEN2_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_KRAKEN2_name,
                self.default_parameters.get_kraken2_default(
                    user, Software.TYPE_OF_USE_televir_global, technology_name
                ),
                user,
            )
            return self.get_kraken2_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_BWA_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_BWA_name,
                self.default_parameters.get_bwa_default(
                    user,
                    Software.TYPE_OF_USE_televir_global,
                    ConstantsSettings.TECHNOLOGY_illumina,
                ),
                user,
            )
            return self.get_bwa_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_DIAMOND_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_DIAMOND_name,
                self.default_parameters.get_diamond_default(
                    user, Software.TYPE_OF_USE_televir_global, technology_name
                ),
                user,
            )
            return self.get_diamond_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_KRAKENUNIQ_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_KRAKENUNIQ_name,
                self.default_parameters.get_krakenuniq_default(
                    user, Software.TYPE_OF_USE_televir_global, technology_name
                ),
                user,
            )

            return self.get_krakenuniq_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_KAIJU_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_KAIJU_name,
                self.default_parameters.get_kaiju_default(
                    user,
                    Software.TYPE_OF_USE_televir_global,
                    technology_name,
                    pipeline_step=ConstantsSettings.PIPELINE_NAME_read_classification,
                ),
                user,
            )
            return self.get_kaiju_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_BLAST_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_BLAST_name,
                self.default_parameters.get_blast_default(
                    user, Software.TYPE_OF_USE_televir_global, technology_name
                ),
                user,
            )
            return self.get_blast_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_DESAMBA_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_DESAMBA_name,
                self.default_parameters.get_desamba_default(
                    user, Software.TYPE_OF_USE_televir_global, technology_name
                ),
                user,
            )
            return self.get_desamba_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_RAVEN_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_RAVEN_name,
                self.default_parameters.get_raven_default(
                    user, Software.TYPE_OF_USE_televir_global, technology_name
                ),
                user,
            )
            return self.get_raven_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_FASTVIROMEEXPLORER_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_FASTVIROMEEXPLORER_name,
                self.default_parameters.get_fastviromeexplorer_default(
                    user, Software.TYPE_OF_USE_televir_global, technology_name
                ),
                user,
            )
            return self.get_fastviromeexplorer_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_SPAdes_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_SPAdes_name,
                self.default_parameters.get_spades_default(
                    user, Software.TYPE_OF_USE_televir_global, technology_name
                ),
                user,
            )
            return self.get_spades_parameters(user, technology_name)

        if software_name == SoftwareNames.SOFTWARE_SNIPPY_PI_name:
            self.test_default_db(
                SoftwareNames.SOFTWARE_SNIPPY_PI_name,
                self.default_parameters.get_snippy_pi_default(
                    user, Software.TYPE_OF_USE_televir_global, technology_name
                ),
                user,
            )
            return self.get_snippy_pi_parameters(user, technology_name)

        return ""

    def get_all_software(self):
        """
        get all softwares available by this class
        """
        vect_software = []
        vect_software.append(self.software_names.get_trimmomatic_name())
        vect_software.append(self.software_names.get_snippy_name())
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
