"""
Created on 03/05/2020

@author: mmp
"""
import logging

from constants.software_names import SoftwareNames
from managing_files.models import Project, ProjectSample, Sample
from utils.lock_atomic_transaction import LockedAtomicTransaction

from settings.constants_settings import ConstantsSettings
from settings.default_parameters import DefaultParameters
from settings.default_software import DefaultSoftware
from settings.models import Parameter, Software


class DefaultProjectSoftware(object):
    """
    classdocs
    """

    software_names = SoftwareNames()

    def __init__(self):
        """change values"""
        self.default_parameters = DefaultParameters()
        self.change_values_software = {}  ### the key is the name of the software

    def test_all_defaults(self, user, project, project_sample, sample, dataset=None):
        """
        test all defaults for all software available
        """

        ## only for project and for all technology
        if not project is None:
            self.test_default_db(
                SoftwareNames.SOFTWARE_SNIPPY_name,
                user,
                Software.TYPE_OF_USE_project,
                project,
                None,
                None,
                ConstantsSettings.TECHNOLOGY_illumina,
            )
            self.test_default_db(
                SoftwareNames.SOFTWARE_FREEBAYES_name,
                user,
                Software.TYPE_OF_USE_project,
                project,
                None,
                None,
                ConstantsSettings.TECHNOLOGY_illumina,
            )
            self.test_default_db(
                SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
                user,
                Software.TYPE_OF_USE_project,
                project,
                None,
                None,
                ConstantsSettings.TECHNOLOGY_illumina,
            )
            self.test_default_db(
                SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
                user,
                Software.TYPE_OF_USE_project,
                project,
                None,
                None,
                ConstantsSettings.TECHNOLOGY_minion,
            )
            self.test_default_db(
                SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,
                user,
                Software.TYPE_OF_USE_project,
                project,
                None,
                None,
                ConstantsSettings.TECHNOLOGY_minion,
            )
            self.test_default_db(
                SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name,
                user,
                Software.TYPE_OF_USE_project,
                project,
                None,
                None,
                ConstantsSettings.TECHNOLOGY_minion,
            )
            self.test_default_db(
                SoftwareNames.SOFTWARE_Medaka_name_consensus,
                user,
                Software.TYPE_OF_USE_project,
                project,
                None,
                None,
                ConstantsSettings.TECHNOLOGY_minion,
            )

            ### SOFTWARE_MASK_CONSENSUS_BY_SITE_name can be both
            self.test_default_db(
                SoftwareNames.SOFTWARE_MASK_CONSENSUS_BY_SITE_name,
                user,
                Software.TYPE_OF_USE_project,
                project,
                None,
                None,
                ConstantsSettings.TECHNOLOGY_generic,
            )
            self.test_default_db(
                SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT,
                user,
                Software.TYPE_OF_USE_project,
                project,
                None,
                None,
                ConstantsSettings.TECHNOLOGY_minion,
            )
            ## TELEVIR SOFTWARE

            ############# PATHOGEN IDENTIFICATION SOFTWARE
            #############
            # self.test_default_db(
            #    SoftwareNames.SOFTWARE_CENTRIFUGE_name,
            #    user,
            #    Software.TYPE_OF_USE_televir_project,
            #    project,
            #    None,
            #    None,
            #    ConstantsSettings.TECHNOLOGY_illumina,
            # )
        elif not dataset is None:
            self.test_default_db(
                SoftwareNames.SOFTWARE_NEXTSTRAIN_name,
                user,
                Software.TYPE_OF_USE_dataset,
                project=None,
                project_sample=None,
                sample=None,
                technology_name=ConstantsSettings.TECHNOLOGY_generic,
                dataset=dataset,
            )

            ## only for project sample and by technology
        elif not project_sample is None:

            ### both technologies
            self.test_default_db(
                SoftwareNames.SOFTWARE_MASK_CONSENSUS_BY_SITE_name,
                user,
                Software.TYPE_OF_USE_project_sample,
                None,
                project_sample,
                None,
                ConstantsSettings.TECHNOLOGY_generic,
            )
            self.test_default_db(
                SoftwareNames.SOFTWARE_GENERATE_CONSENSUS_name,
                user,
                Software.TYPE_OF_USE_project_sample,
                None,
                project_sample,
                None,
                ConstantsSettings.TECHNOLOGY_generic,
            )

            ### illumina
            if project_sample.sample.is_type_fastq_gz_sequencing():
                self.test_default_db(
                    SoftwareNames.SOFTWARE_SNIPPY_name,
                    user,
                    Software.TYPE_OF_USE_project_sample,
                    None,
                    project_sample,
                    None,
                    ConstantsSettings.TECHNOLOGY_illumina,
                )
                self.test_default_db(
                    SoftwareNames.SOFTWARE_FREEBAYES_name,
                    user,
                    Software.TYPE_OF_USE_project_sample,
                    None,
                    project_sample,
                    None,
                    ConstantsSettings.TECHNOLOGY_illumina,
                )
                self.test_default_db(
                    SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
                    user,
                    Software.TYPE_OF_USE_project_sample,
                    None,
                    project_sample,
                    None,
                    ConstantsSettings.TECHNOLOGY_illumina,
                )

            else:  ## ONT
                self.test_default_db(
                    SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
                    user,
                    Software.TYPE_OF_USE_project_sample,
                    None,
                    project_sample,
                    None,
                    ConstantsSettings.TECHNOLOGY_minion,
                )
                self.test_default_db(
                    SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,
                    user,
                    Software.TYPE_OF_USE_project_sample,
                    None,
                    project_sample,
                    None,
                    ConstantsSettings.TECHNOLOGY_minion,
                )
                self.test_default_db(
                    SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name,
                    user,
                    Software.TYPE_OF_USE_project_sample,
                    None,
                    project_sample,
                    None,
                    ConstantsSettings.TECHNOLOGY_minion,
                )
                self.test_default_db(
                    SoftwareNames.SOFTWARE_Medaka_name_consensus,
                    user,
                    Software.TYPE_OF_USE_project_sample,
                    None,
                    project_sample,
                    None,
                    ConstantsSettings.TECHNOLOGY_minion,
                )
                self.test_default_db(
                    SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT,
                    user,
                    Software.TYPE_OF_USE_project_sample,
                    None,
                    project_sample,
                    None,
                    ConstantsSettings.TECHNOLOGY_minion,
                )

        ### only for sample and ONT technology
        elif not sample is None:
            if sample.is_type_fastq_gz_sequencing(Sample.TYPE_OF_FASTQ_minion):
                self.test_default_db(
                    SoftwareNames.SOFTWARE_NanoFilt_name,
                    user,
                    Software.TYPE_OF_USE_sample,
                    None,
                    None,
                    sample,
                    ConstantsSettings.TECHNOLOGY_minion,
                )
                self.test_default_db(
                    SoftwareNames.SOFTWARE_ABRICATE_name,
                    user,
                    Software.TYPE_OF_USE_sample,
                    None,
                    None,
                    sample,
                    ConstantsSettings.TECHNOLOGY_minion,
                )
            if sample.is_type_fastq_gz_sequencing(Sample.TYPE_OF_FASTQ_illumina):
                self.test_default_db(
                    SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
                    user,
                    Software.TYPE_OF_USE_sample,
                    None,
                    None,
                    sample,
                    ConstantsSettings.TECHNOLOGY_illumina,
                )
                self.test_default_db(
                    SoftwareNames.SOFTWARE_ABRICATE_name,
                    user,
                    Software.TYPE_OF_USE_sample,
                    None,
                    None,
                    sample,
                    ConstantsSettings.TECHNOLOGY_illumina,
                )

    def test_default_db(
        self,
        software_name,
        user,
        type_of_use,
        project,
        project_sample,
        sample,
        technology_name,
        dataset=None,
        televir_project=None,
    ):
        """
        test if exist, if not persist in database
        """
        ## lock because more than one process can duplicate software names
        with LockedAtomicTransaction(Software), LockedAtomicTransaction(Parameter):
            list_software = Software.objects.filter(
                name=software_name,
                owner=user,
                type_of_use=type_of_use,
                parameter__project=project,
                parameter__project_sample=project_sample,
                parameter__sample=sample,
                parameter__dataset=dataset,
                parameter__televir_project=televir_project,
                version_parameters=self.default_parameters.get_software_parameters_version(
                    software_name
                ),
                technology__name=technology_name,
            ).distinct("name")

            # logger = logging.getLogger("fluWebVirus.debug")
            # logger.debug("Test default db: {} ({})".format(list_software, len(list_software)))

            ### if not exist need to save
            if len(list_software) == 0:
                vect_parameters = self._get_default_parameters(
                    software_name,
                    user,
                    type_of_use,
                    project,
                    project_sample,
                    sample,
                    technology_name,
                    dataset,
                )
                if len(vect_parameters) > 0:  ### persist
                    self.default_parameters.persist_parameters(
                        vect_parameters, type_of_use
                    )

    def _get_default_parameters(
        self,
        software_name,
        user,
        type_of_use,
        project,
        project_sample,
        sample,
        technology_name,
        dataset=None,
    ):
        if software_name == SoftwareNames.SOFTWARE_SNIPPY_name:
            vect_parameters = self.default_parameters.get_snippy_default(
                user, type_of_use, technology_name, project, project_sample
            )  ### base values
            if not project is None:
                vect_parameters = self._get_default_project(
                    user,
                    SoftwareNames.SOFTWARE_SNIPPY_name,
                    None,
                    vect_parameters,
                    technology_name,
                )  ### base values
            if not project_sample is None:
                vect_parameters = self._get_default_project(
                    user,
                    SoftwareNames.SOFTWARE_SNIPPY_name,
                    project_sample.project,
                    vect_parameters,
                    technology_name,
                )  ### base values
            return vect_parameters
        elif software_name == SoftwareNames.SOFTWARE_FREEBAYES_name:
            vect_parameters = self.default_parameters.get_freebayes_default(
                user, type_of_use, technology_name, project, project_sample
            )
            if not project is None:
                vect_parameters = self._get_default_project(
                    user,
                    SoftwareNames.SOFTWARE_FREEBAYES_name,
                    None,
                    vect_parameters,
                    technology_name,
                )  ### base values
            if not project_sample is None:
                vect_parameters = self._get_default_project(
                    user,
                    SoftwareNames.SOFTWARE_FREEBAYES_name,
                    project_sample.project,
                    vect_parameters,
                    technology_name,
                )  ### base values
            return vect_parameters
        elif software_name == SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name:
            vect_parameters = (
                self.default_parameters.get_mask_consensus_threshold_default(
                    user, type_of_use, technology_name, project, project_sample
                )
            )
            if not project is None:
                vect_parameters = self._get_default_project(
                    user, software_name, None, vect_parameters, technology_name
                )  ### base values
            if not project_sample is None:
                vect_parameters = self._get_default_project(
                    user,
                    software_name,
                    project_sample.project,
                    vect_parameters,
                    technology_name,
                )  ### base values
            return vect_parameters
        elif software_name == SoftwareNames.SOFTWARE_MASK_CONSENSUS_BY_SITE_name:
            vect_parameters = (
                self.default_parameters.get_mask_consensus_by_site_default(
                    user, type_of_use, technology_name, project, project_sample
                )
            )
            if not project is None:
                vect_parameters = self._get_default_project(
                    user, software_name, None, vect_parameters, technology_name
                )  ### base values
            if not project_sample is None:
                vect_parameters = self._get_default_project(
                    user,
                    software_name,
                    project_sample.project,
                    vect_parameters,
                    technology_name,
                )  ### base values
            return vect_parameters
        elif software_name == SoftwareNames.SOFTWARE_GENERATE_CONSENSUS_name:
            vect_parameters = self.default_parameters.get_generate_consensus_default(
                user, project_sample
            )
            return vect_parameters
        elif software_name == SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name:
            vect_parameters = (
                self.default_parameters.get_limit_coverage_ONT_threshold_default(
                    user, type_of_use, technology_name, project, project_sample
                )
            )
            if not project is None:
                vect_parameters = self._get_default_project(
                    user, software_name, None, vect_parameters, technology_name
                )  ### base values
            if not project_sample is None:
                vect_parameters = self._get_default_project(
                    user,
                    software_name,
                    project_sample.project,
                    vect_parameters,
                    technology_name,
                )  ### base values
            return vect_parameters
        elif software_name == SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name:
            vect_parameters = (
                self.default_parameters.get_vcf_freq_ONT_threshold_default(
                    user, type_of_use, technology_name, project, project_sample
                )
            )
            if not project is None:
                vect_parameters = self._get_default_project(
                    user, software_name, None, vect_parameters, technology_name
                )  ### base values
            if not project_sample is None:
                vect_parameters = self._get_default_project(
                    user,
                    software_name,
                    project_sample.project,
                    vect_parameters,
                    technology_name,
                )  ### base values
            return vect_parameters
        elif software_name == SoftwareNames.SOFTWARE_Medaka_name_consensus:
            vect_parameters = self.default_parameters.get_medaka_model_default(
                user, type_of_use, technology_name, project, project_sample
            )
            if not project is None:
                vect_parameters = self._get_default_project(
                    user, software_name, None, vect_parameters, technology_name
                )  ### base values
            if not project_sample is None:
                vect_parameters = self._get_default_project(
                    user,
                    software_name,
                    project_sample.project,
                    vect_parameters,
                    technology_name,
                )  ### base values
            return vect_parameters
        elif software_name == SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT:
            vect_parameters = self.default_parameters.get_samtools_depth_default_ONT(
                user, type_of_use, technology_name, project, project_sample
            )
            if not project is None:
                vect_parameters = self._get_default_project(
                    user, software_name, None, vect_parameters, technology_name
                )  ### base values
            if not project_sample is None:
                vect_parameters = self._get_default_project(
                    user,
                    software_name,
                    project_sample.project,
                    vect_parameters,
                    technology_name,
                )  ### base values
            return vect_parameters
        elif software_name == SoftwareNames.SOFTWARE_NanoFilt_name:
            vect_parameters = self.default_parameters.get_nanofilt_default(
                user, type_of_use, technology_name, sample
            )
            if not sample is None:
                vect_parameters = self._get_default_project(
                    user, software_name, None, vect_parameters, technology_name
                )  ### base values
            return vect_parameters
        elif software_name == SoftwareNames.SOFTWARE_TRIMMOMATIC_name:
            vect_parameters = self.default_parameters.get_trimmomatic_default(
                user, type_of_use, technology_name, sample
            )
            if not sample is None:
                vect_parameters = self._get_default_project(
                    user, software_name, None, vect_parameters, technology_name
                )  ### base values
            return vect_parameters
        elif software_name == SoftwareNames.SOFTWARE_ABRICATE_name:
            vect_parameters = self.default_parameters.get_abricate_default(
                user, type_of_use, technology_name, sample
            )
            if not sample is None:
                vect_parameters = self._get_default_project(
                    user, software_name, None, vect_parameters, technology_name
                )  ### base values
            return vect_parameters
        elif software_name == SoftwareNames.SOFTWARE_NEXTSTRAIN_name:
            vect_parameters = self.default_parameters.get_nextstrain_default(
                user=user, dataset=dataset
            )
            return vect_parameters
        return []

    #####################################################
    #####
    #####        snippy
    #####

    def get_snippy_parameters(self, user, type_of_use, project, project_sample):
        """
        get snippy parameters
        """
        return self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_SNIPPY_name,
            user,
            type_of_use,
            project,
            project_sample,
            None,
            ConstantsSettings.TECHNOLOGY_illumina,
        )

    def get_snippy_parameters_all_possibilities(self, user, project_sample):
        """
        get snippy parameters for project_sample, project and default
        """

        ### Test project_sample first
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_SNIPPY_name,
            user,
            Software.TYPE_OF_USE_project_sample,
            None,
            project_sample,
            None,
            ConstantsSettings.TECHNOLOGY_illumina,
        )
        if not parameters is None:
            return parameters

        ### Test project
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_SNIPPY_name,
            user,
            Software.TYPE_OF_USE_project,
            project_sample.project,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_illumina,
        )
        if not parameters is None:
            return parameters

        ### can be a default one
        default_software = DefaultSoftware()
        parameters = default_software.get_snippy_parameters(user)
        if len(parameters) > 0:
            return parameters

        software_names = SoftwareNames()
        return software_names.get_snippy_parameters()

    def get_snippy_parameters_for_project(self, user, project):
        """
        get snippy parameters only for project or default
        """

        ### Test project
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_SNIPPY_name,
            user,
            Software.TYPE_OF_USE_project,
            project,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_illumina,
        )
        if not parameters is None:
            return parameters

        ### can be a default one
        default_software = DefaultSoftware()
        parameters = default_software.get_snippy_parameters(user)
        if len(parameters) > 0:
            return parameters
        return None

    def get_snippy_single_parameter_default(self, parameter_name):
        """
        :param parameter_name -> Only these two possibilities available SNIPPY_COVERAGE_NAME; SNIPPY_MAPQUAL_NAME
        :return value of parameter
        """
        vect_parameters = self.default_parameters.get_snippy_default(None, None, None)
        for parameters in vect_parameters:
            if parameters.name == parameter_name:
                return parameters.parameter
        return None

    def is_snippy_single_parameter_default(self, project_sample, parameter_name):
        """
        test if a specific parameter is default SNIPPY_COVERAGE_NAME; SNIPPY_MAPQUAL_NAME
        """

        value_default_parameter = self.get_snippy_single_parameter_default(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_snippy_single_parameter(
            project_sample, parameter_name
        )
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_snippy_single_parameter(self, project_sample, parameter_name):
        """
        get snippy single parameters
        :param parameter_name -> Only these two possibilities available SNIPPY_COVERAGE_NAME; SNIPPY_MAPQUAL_NAME
        """

        parameters_string = self.get_snippy_parameters_all_possibilities(
            project_sample.project.owner, project_sample
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1].split()[0]
        return None

    def is_snippy_single_parameter_default_for_project(self, project, parameter_name):
        """
        test if a specific parameter is default SNIPPY_COVERAGE_NAME; SNIPPY_MAPQUAL_NAME
        """

        value_default_parameter = self.get_snippy_single_parameter_default(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_snippy_single_parameter_for_project(
            project, parameter_name
        )
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_snippy_single_parameter_for_project(self, project, parameter_name):
        """
        get snippy single parameters
        :param parameter_name -> Only these two possibilities available SNIPPY_COVERAGE_NAME; SNIPPY_MAPQUAL_NAME
        """

        parameters_string = self.get_snippy_parameters_for_project(
            project.owner, project
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1].split()[0]
        return None

    #####
    #####        END snippy
    #####
    #####################################################

    #####################################################
    #####
    #####        BEGIN -- Freebayes
    #####

    def get_freebayes_parameters(self, user, type_of_use, project, project_sample):
        """
        get freebayes parameters
        Add extra -V to the end
        """
        return self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_FREEBAYES_name,
            user,
            type_of_use,
            project,
            project_sample,
            None,
            ConstantsSettings.TECHNOLOGY_illumina,
        )

    def is_to_run_freebayes(self, user, project_sample):
        """
        test if freebayes is to run
        """

        ### Test project_sample first
        return self.default_parameters.is_software_to_run(
            SoftwareNames.SOFTWARE_FREEBAYES_name,
            user,
            Software.TYPE_OF_USE_project_sample,
            None,
            project_sample,
            None,
            ConstantsSettings.TECHNOLOGY_illumina,
        )

    def set_freebayes_to_run(self, user, project_sample, is_to_run):
        """
        test if abricate is to run
        """
        return self.default_parameters.set_software_to_run(
            SoftwareNames.SOFTWARE_NanoFilt_name,
            user,
            Software.TYPE_OF_USE_project_sample,
            None,
            project_sample,
            None,
            ConstantsSettings.TECHNOLOGY_illumina,
            is_to_run,
        )

    #####
    #####        END -- Freebayes
    #####
    #####################################################

    #####################################################
    #####
    #####        nanofilt
    #####

    def get_nanofilt_parameters(self, user, type_of_use, sample):
        """
        get nanofilt parameters
        """
        return self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_NanoFilt_name,
            user,
            type_of_use,
            None,
            None,
            sample,
            ConstantsSettings.TECHNOLOGY_minion,
        )

    def get_nanofilt_parameters_all_possibilities(self, user, sample):
        """
        get nanofilt parameters for project_sample, project and default
        """

        ### Test project_sample first
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_NanoFilt_name,
            user,
            Software.TYPE_OF_USE_sample,
            None,
            None,
            sample,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        if not parameters is None:
            return parameters

        ### can be a default one
        default_software = DefaultSoftware()
        parameters = default_software.get_nanofilt_parameters(user)
        if len(parameters) > 0:
            return parameters

        software_names = SoftwareNames()
        return software_names.get_nanofilt_parameters()

    def is_to_run_nanofilt(self, user, sample):
        """
        test if nanofilt is to run
        """

        ### Test project_sample first
        return self.default_parameters.is_software_to_run(
            SoftwareNames.SOFTWARE_NanoFilt_name,
            user,
            Software.TYPE_OF_USE_sample,
            None,
            None,
            sample,
            ConstantsSettings.TECHNOLOGY_minion,
        )

    def set_nanofilt_to_run(self, user, sample, is_to_run):
        """
        test if abricate is to run
        """
        return self.default_parameters.set_software_to_run(
            SoftwareNames.SOFTWARE_NanoFilt_name,
            user,
            Software.TYPE_OF_USE_sample,
            None,
            None,
            sample,
            ConstantsSettings.TECHNOLOGY_minion,
            is_to_run,
        )

    def get_nanofilt_parameters_for_sample(self, user, sample):
        """
        get nanofilt parameters only for project or default
        """

        ### Test project
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_NanoFilt_name,
            user,
            Software.TYPE_OF_USE_sample,
            None,
            None,
            sample,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        if not parameters is None:
            return parameters

        ### can be a default one
        default_software = DefaultSoftware()
        parameters = default_software.get_nanofilt_parameters(user)
        if len(parameters) > 0:
            return parameters
        return None

    def get_nanofilt_single_parameter_default(self, parameter_name):
        """
        :param parameter_name -> NANOfilt_quality_read
        :return value of parameter
        """
        vect_parameters = self.default_parameters.get_nanofilt_default(None, None, None)
        for parameters in vect_parameters:
            if parameters.name == parameter_name:
                return parameters.parameter
        return None

    def is_nanofilt_single_parameter_default(self, sample, parameter_name):
        """
        test if a specific parameter is default NANOfilt_quality_read
        """

        value_default_parameter = self.get_nanofilt_single_parameter_default(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_nanofilt_single_parameter(sample, parameter_name)
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_nanofilt_single_parameter(self, sample, parameter_name):
        """
        get nanofilt single parameters
        :param parameter_name -> Only these two possibilities available NANOfilt_quality_read
        """

        parameters_string = self.get_nanofilt_parameters_all_possibilities(
            sample.owner, sample
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1].split()[0]
        return None

    def is_nanofilt_single_parameter_default_for_project(self, sample, parameter_name):
        """
        test if a specific parameter is default NANOfilt_quality_read
        """

        value_default_parameter = self.get_nanofilt_single_parameter_default(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_nanofilt_single_parameter_for_sample(
            sample, parameter_name
        )
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_nanofilt_single_parameter_for_sample(self, sample, parameter_name):
        """
        get nanofilt single parameters
        :param parameter_name -> Only these two possibilities available NANOfilt_quality_read
        """

        parameters_string = self.get_nanofilt_parameters_for_sample(
            sample.owner, sample
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1].split()[0]
        return None

    #####
    #####        END nanofilt
    #####
    #####################################################

    #####################################################
    #####
    #####        Trimmomatic
    #####

    def get_trimmomatic_parameters(self, user, type_of_use, sample):
        """
        get trimmomatic parameters
        """
        return self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
            user,
            type_of_use,
            None,
            None,
            sample,
            ConstantsSettings.TECHNOLOGY_illumina,
        )

    def get_trimmomatic_parameters_all_possibilities(self, user, sample):
        """
        get trimmomatic parameters for project_sample, project and default
        """

        ### Test project_sample first
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
            user,
            Software.TYPE_OF_USE_sample,
            None,
            None,
            sample,
            ConstantsSettings.TECHNOLOGY_illumina,
        )
        if not parameters is None:
            return parameters

        ### can be a default one
        default_software = DefaultSoftware()
        parameters = default_software.get_trimmomatic_parameters(user)
        if len(parameters) > 0:
            return parameters

        software_names = SoftwareNames()
        return software_names.get_trimmomatic_parameters()

    def is_to_run_trimmomatic(self, user, sample):
        """
        test if trimmomatic is to run
        """

        ### Test project_sample first
        return self.default_parameters.is_software_to_run(
            SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
            user,
            Software.TYPE_OF_USE_sample,
            None,
            None,
            sample,
            ConstantsSettings.TECHNOLOGY_illumina,
        )

    def set_trimmomatic_to_run(self, user, sample, is_to_run):
        """
        test if trimmomatic is to run
        """
        return self.default_parameters.set_software_to_run(
            SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
            user,
            Software.TYPE_OF_USE_sample,
            None,
            None,
            sample,
            ConstantsSettings.TECHNOLOGY_illumina,
            is_to_run,
        )

    def get_trimmomatic_parameters_for_sample(self, user, sample):
        """
        get trimmomatic parameters only for project or default
        """

        ### Test project
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_TRIMMOMATIC_name,
            user,
            Software.TYPE_OF_USE_sample,
            None,
            None,
            sample,
            ConstantsSettings.TECHNOLOGY_illumina,
        )
        if not parameters is None:
            return parameters

        ### can be a default one
        default_software = DefaultSoftware()
        parameters = default_software.get_trimmomatic_parameters(user)
        if len(parameters) > 0:
            return parameters
        return None

    def get_trimmomatic_single_parameter_default(self, parameter_name):
        """
        :param parameter_name -> NANOfilt_quality_read
        :return value of parameter
        """
        vect_parameters = self.default_parameters.get_trimmomatic_default(
            None, None, None
        )
        for parameters in vect_parameters:
            if parameters.name == parameter_name:
                return parameters.parameter
        return None

    def is_trimmomatic_single_parameter_default(self, sample, parameter_name):
        """
        test if a specific parameter is default NANOfilt_quality_read
        """

        value_default_parameter = self.get_trimmomatic_single_parameter_default(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_trimmomatic_single_parameter(
            sample, parameter_name
        )
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_trimmomatic_single_parameter(self, sample, parameter_name):
        """
        get trimmomatic single parameters
        :param parameter_name -> Only these two possibilities available NANOfilt_quality_read
        """

        parameters_string = self.get_trimmomatic_parameters_all_possibilities(
            sample.owner, sample
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1].split()[0]
        return None

    def is_trimmomatic_single_parameter_default_for_project(
        self, sample, parameter_name
    ):
        """
        test if a specific parameter is default NANOfilt_quality_read
        """

        value_default_parameter = self.get_trimmomatic_single_parameter_default(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_trimmomatic_single_parameter_for_sample(
            sample, parameter_name
        )
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_trimmomatic_single_parameter_for_sample(self, sample, parameter_name):
        """
        get trimmomatic single parameters
        :param parameter_name -> Only these two possibilities available NANOfilt_quality_read
        """

        parameters_string = self.get_trimmomatic_parameters_for_sample(
            sample.owner, sample
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1].split()[0]
        return None

    #####
    #####        END trimmomatic
    #####
    #####################################################

    #####################################################
    #####
    #####        Abricate
    #####

    def get_abricate_parameters(self, user, type_of_use, sample, technology_name):
        """
        get abricate parameters
        """
        return self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_ABRICATE_name,
            user,
            type_of_use,
            None,
            None,
            sample,
            technology_name,
        )

    def is_to_run_abricate(self, user, sample, technology_name):
        """
        test if abricate is to run
        """

        ### Test project_sample first
        return self.default_parameters.is_software_to_run(
            SoftwareNames.SOFTWARE_ABRICATE_name,
            user,
            Software.TYPE_OF_USE_sample,
            None,
            None,
            sample,
            technology_name,
        )

    def set_abricate_to_run(self, user, sample, technology_name, is_to_run):
        """
        test if abricate is to run
        """

        ### Test project_sample first
        return self.default_parameters.set_software_to_run(
            SoftwareNames.SOFTWARE_ABRICATE_name,
            user,
            Software.TYPE_OF_USE_sample,
            None,
            None,
            sample,
            technology_name,
            is_to_run,
        )

    def get_abricate_parameters_all_possibilities(self, user, sample, technology_name):
        """
        get abricate parameters for sample
        """

        ### Test project_sample first
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_ABRICATE_name,
            user,
            Software.TYPE_OF_USE_sample,
            None,
            None,
            sample,
            technology_name,
        )
        if not parameters is None:
            return parameters

        ### can be a default one
        default_software = DefaultSoftware()
        parameters = default_software.get_abricate_parameters(user, technology_name)
        if len(parameters) > 0:
            return parameters

        software_names = SoftwareNames()
        return software_names.get_trimmomatic_parameters()

    def get_abricatec_parameters_for_sample(self, user, sample, technology_name):
        """
        get abricate parameters only for project or default
        """

        ### Test project
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_ABRICATE_name,
            user,
            Software.TYPE_OF_USE_sample,
            None,
            None,
            sample,
            technology_name,
        )
        if not parameters is None:
            return parameters

        ### can be a default one
        default_software = DefaultSoftware()
        parameters = default_software.get_abricate_parameters(user, technology_name)
        if len(parameters) > 0:
            return parameters
        return None

    def get_abricate_single_parameter_default(self, parameter_name):
        """
        :return value of parameter
        """
        vect_parameters = self.default_parameters.get_abricate_default(None, None, None)
        for parameters in vect_parameters:
            if parameters.name == parameter_name:
                return parameters.parameter
        return None

    def is_abricate_single_parameter_default(
        self, sample, parameter_name, technology_name
    ):
        """ """

        value_default_parameter = self.get_abricate_single_parameter_default(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_abricate_single_parameter(
            sample, parameter_name, technology_name
        )
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_abricate_single_parameter(self, sample, parameter_name, technology_name):
        """
        get abricate single parameters
        :param parameter_name -> Only these two possibilities available NANOfilt_quality_read
        """

        parameters_string = self.get_abricate_parameters_all_possibilities(
            sample.owner, sample, technology_name
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1].split()[0]
        return None

    def is_abricate_single_parameter_default_for_project(
        self, sample, parameter_name, technology_name
    ):
        """
        test if a specific parameter is default NANOfilt_quality_read
        """

        value_default_parameter = self.get_abricate_single_parameter_default(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_abricate_single_parameter_for_sample(
            sample, parameter_name, technology_name
        )
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_abricate_single_parameter_for_sample(
        self, sample, parameter_name, technology_name
    ):
        """
        get trimmomatic single parameters
        :param parameter_name -> Only these two possibilities available NANOfilt_quality_read
        """

        parameters_string = self.get_abricatec_parameters_for_sample(
            sample.owner, sample, technology_name
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1].split()[0]
        return None

    #####
    #####        END abricate
    #####
    #####################################################

    #####################################################
    #####
    #####        Medaka
    #####

    def get_medaka_parameters(self, user, type_of_use, project, project_sample):
        """
        get medaka parameters
        """
        return self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_Medaka_name_consensus,
            user,
            type_of_use,
            project,
            project_sample,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )

    def get_medaka_parameters_all_possibilities(self, user, project_sample):
        """
        get medaka parameters for project_sample, project and default
        get the model to run, this is the only possibility
        """

        ### Test project_sample first
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_Medaka_name_consensus,
            user,
            Software.TYPE_OF_USE_project_sample,
            None,
            project_sample,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        if not parameters is None:
            return parameters

        ### Test project
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_Medaka_name_consensus,
            user,
            Software.TYPE_OF_USE_project,
            project_sample.project,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        if not parameters is None:
            return parameters

        ### can be a default one
        default_software = DefaultSoftware()
        parameters = default_software.get_medaka_parameters_consensus(user)
        if len(parameters) > 0:
            return parameters

        software_names = SoftwareNames()
        return software_names.get_medaka_parameters_consensus()

    def get_medaka_parameters_for_project(self, user, project):
        """
        get medaka parameters only for project or default
        """

        ### Test project
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_Medaka_name_consensus,
            user,
            Software.TYPE_OF_USE_project,
            project,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        if not parameters is None:
            return parameters

        ### can be a default one
        default_software = DefaultSoftware()
        parameters = default_software.get_medaka_parameters_consensus(user)
        if len(parameters) > 0:
            return parameters
        return None

    def get_medaka_single_parameter_default(self, parameter_name):
        """
        :param parameter_name -> Only these two possibilities available MEDAKA_model
        :return value of parameter
        """
        vect_parameters = self.default_parameters.get_medaka_model_default(
            None, None, None, None
        )
        for parameters in vect_parameters:
            if parameters.name == parameter_name:
                return parameters.parameter
        return None

    def is_medaka_single_parameter_default(self, project_sample, parameter_name):
        """
        test if a specific parameter is default MEDAKA_model
        """

        value_default_parameter = self.get_medaka_single_parameter_default(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_medaka_single_parameter(
            project_sample, parameter_name
        )
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_medaka_single_parameter(self, project_sample, parameter_name):
        """
        get medaka single parameters
        :param parameter_name -> Only these two possibilities available MEDAKA_model
        """

        parameters_string = self.get_medaka_parameters_all_possibilities(
            project_sample.project.owner, project_sample
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1].split()[0]
        return None

    def is_medaka_single_parameter_default_for_project(self, project, parameter_name):
        """
        test if a specific parameter is default MEDAKA_model
        """

        value_default_parameter = self.get_medaka_single_parameter_default(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_medaka_single_parameter_for_project(
            project, parameter_name
        )
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_medaka_single_parameter_for_project(self, project, parameter_name):
        """
        get medaka single parameters
        :param parameter_name -> Only these two possibilities available MEDAKA_model
        """

        parameters_string = self.get_medaka_parameters_for_project(
            project.owner, project
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1].split()[0]
        return None

    #####
    #####        END medaka
    #####
    #####################################################

    #####################################################
    #####
    #####        Samtools ONT
    #####

    def get_samtools_parameters_ONT(self, user, type_of_use, project, project_sample):
        """
        get samtools parameters
        """
        return self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT,
            user,
            type_of_use,
            project,
            project_sample,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )

    def get_samtools_parameters_all_possibilities_ONT(self, user, project_sample):
        """
        get samtools parameters ONT for project_sample, project and default
        get the model to run, this is the only possibility
        """

        ### Test project_sample first
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT,
            user,
            Software.TYPE_OF_USE_project_sample,
            None,
            project_sample,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        if not parameters is None:
            return parameters

        ### Test project
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT,
            user,
            Software.TYPE_OF_USE_project,
            project_sample.project,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        if not parameters is None:
            return parameters

        ### can be a default one
        default_software = DefaultSoftware()
        parameters = default_software.get_samtools_parameters_depth_ONT(user)
        if len(parameters) > 0:
            return parameters

        software_names = SoftwareNames()
        return software_names.get_samtools_parameters_depth_ONT()

    def get_samtools_parameters_for_project_ONT(self, user, project):
        """
        get samtools parameters only for project or default
        """

        ### Test project
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT,
            user,
            Software.TYPE_OF_USE_project,
            project,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        if not parameters is None:
            return parameters

        ### can be a default one
        default_software = DefaultSoftware()
        parameters = default_software.get_samtools_parameters_depth_ONT(user)
        if len(parameters) > 0:
            return parameters
        return None

    def get_samtools_single_parameter_default_ONT(self, parameter_name):
        """
        :param parameter_name -> Only these two possibilities available Samtools ONT
        :return value of parameter
        """
        vect_parameters = self.default_parameters.get_samtools_depth_default_ONT(
            None, None, None, None
        )
        for parameters in vect_parameters:
            if parameters.name == parameter_name:
                return parameters.parameter
        return None

    def is_samtools_single_parameter_default_ONT(self, project_sample, parameter_name):
        """
        test if a specific parameter is default Samtools ONT
        """

        value_default_parameter = self.get_samtools_single_parameter_default_ONT(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_samtools_single_parameter_ONT(
            project_sample, parameter_name
        )
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_samtools_single_parameter_ONT(self, project_sample, parameter_name):
        """
        get samtools single parameters
        :param parameter_name -> Only these two possibilities available Samtools ONT
        """

        parameters_string = self.get_samtools_parameters_all_possibilities_ONT(
            project_sample.project.owner, project_sample
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1].split()[0]
        return None

    def is_samtools_single_parameter_default_for_project_ONT(
        self, project, parameter_name
    ):
        """
        test if a specific parameter is default Samtools ONT
        """

        value_default_parameter = self.get_samtools_single_parameter_default_ONT(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_samtools_single_parameter_for_project_ONT(
            project, parameter_name
        )
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_samtools_single_parameter_for_project_ONT(self, project, parameter_name):
        """
        get samtools single parameters
        :param parameter_name -> Only these two possibilities available MEDAKA_model
        """

        parameters_string = self.get_samtools_parameters_for_project_ONT(
            project.owner, project
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1].split()[0]
        return None

    #####
    #####        END samtools ONT
    #####
    #####################################################

    #####################################################
    #####
    #####        Mask consensus
    #####

    def get_mask_consensus_parameters(
        self, user, type_of_use, project, project_sample, technology_name
    ):
        """
        get mask_consensus parameters
        """
        return self.default_parameters.get_parameters(
            SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
            user,
            type_of_use,
            project,
            project_sample,
            None,
            technology_name,
        )

    def get_mask_consensus_parameters_all_possibilities(
        self, user, project_sample, technology_name
    ):
        """
        get mask_consensus parameters for project and default
        """

        ### Test project_sample first
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
            user,
            Software.TYPE_OF_USE_project_sample,
            None,
            project_sample,
            None,
            technology_name,
        )
        if not parameters is None:
            return parameters

        ### Test project
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
            user,
            Software.TYPE_OF_USE_project,
            project_sample.project,
            None,
            None,
            technology_name,
        )
        if not parameters is None:
            return parameters

        ### can be a default one
        default_software = DefaultSoftware()
        parameters = default_software.get_mask_consensus_threshold_parameters(
            user, technology_name
        )
        if len(parameters) > 0:
            return parameters

        software_names = SoftwareNames()
        return software_names.get_insaflu_parameter_mask_consensus_parameters()

    def get_mask_consensus_parameters_for_project(self, user, project, technology_name):
        """
        get mask_consensus parameters only for project or default
        """

        ### Test project
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name,
            user,
            Software.TYPE_OF_USE_project,
            project,
            None,
            None,
            technology_name,
        )
        if not parameters is None:
            return parameters

        ### can be a default one
        default_software = DefaultSoftware()
        parameters = default_software.get_mask_consensus_threshold_parameters(
            user, technology_name
        )
        if len(parameters) > 0:
            return parameters
        return None

    def get_mask_consensus_single_parameter_default(self, parameter_name):
        """
        :param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
        :return value of parameter
        """
        vect_parameters = self.default_parameters.get_mask_consensus_threshold_default(
            None, None, None
        )
        for parameters in vect_parameters:
            if parameters.name == parameter_name:
                return parameters.parameter
        return None

    def is_mask_consensus_single_parameter_default(
        self, project_sample, parameter_name, technology_name
    ):
        """
        one possibility available MASK_CONSENSUS_threshold
        """

        value_default_parameter = self.get_mask_consensus_single_parameter_default(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_mask_consensus_single_parameter(
            project_sample, parameter_name, technology_name
        )
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_mask_consensus_single_parameter(
        self, project_sample, parameter_name, technology_name
    ):
        """
        get mask_consensus single parameters
        :param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
        """

        parameters_string = self.get_mask_consensus_parameters_all_possibilities(
            project_sample.project.owner, project_sample, technology_name
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1][1:]
        return None

    def is_mask_consensus_single_parameter_default_for_project(
        self, project, parameter_name, technology_name
    ):
        """
        :param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
        """

        value_default_parameter = self.get_mask_consensus_single_parameter_default(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_mask_consensus_single_parameter_for_project(
            project, parameter_name, technology_name
        )
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_mask_consensus_single_parameter_for_project(
        self, project, parameter_name, technology_name
    ):
        """
        get mask_consensus single parameters
        :param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
        """

        parameters_string = self.get_mask_consensus_parameters_for_project(
            project.owner, project, technology_name
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1][1:]
        return None

    #####
    #####        END Mask consensus
    #####
    #####################################################

    #####################################################
    #####
    #####        Generate consensus
    #####

    def include_consensus(self, project_sample):
        """
        get mask_consensus parameters for project and default
        """

        ### Test project_sample first
        vect_parameters = self.default_parameters.get_list_parameters(
            SoftwareNames.SOFTWARE_GENERATE_CONSENSUS_name,
            project_sample.sample.owner,
            Software.TYPE_OF_USE_project_sample,
            None,
            project_sample,
            None,
            ConstantsSettings.TECHNOLOGY_generic,
        )
        if not vect_parameters is None and len(vect_parameters) > 0:
            return vect_parameters[0].is_to_run
        return True

    #####
    #####        END Generate consensus
    #####
    #####################################################

    #####################################################
    #####
    #####        Coverage limit ONT
    #####

    def get_limit_coverage_ONT_parameters(
        self, user, type_of_use, project, project_sample
    ):
        """
        get limit_coverage_ONT parameters
        """
        return self.default_parameters.get_parameters(
            SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,
            user,
            type_of_use,
            project,
            project_sample,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )

    def get_limit_coverage_ONT_parameters_all_possibilities(self, user, project_sample):
        """
        get limit_coverage parameters for project and default
        """

        ### Test project_sample first
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,
            user,
            Software.TYPE_OF_USE_project_sample,
            None,
            project_sample,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        if not parameters is None:
            return parameters

        ### Test project
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,
            user,
            Software.TYPE_OF_USE_project,
            project_sample.project,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        if not parameters is None:
            return parameters

        ### can be a default one
        default_software = DefaultSoftware()
        parameters = default_software.get_limit_coverage_ONT_parameters(user)
        if len(parameters) > 0:
            return parameters

        software_names = SoftwareNames()
        return software_names.get_insaflu_parameter_limit_coverage_parameters()

    def get_limit_coverage_ONT_parameters_for_project(self, user, project):
        """
        get limit_coverage_ONT parameters only for project or default
        """

        ### Test project
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,
            user,
            Software.TYPE_OF_USE_project,
            project,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        if not parameters is None:
            return parameters

        ### can be a default one
        default_software = DefaultSoftware()
        parameters = default_software.get_limit_coverage_ONT_parameters(user)
        if len(parameters) > 0:
            return parameters
        return None

    def get_limit_coverage_ONT_single_parameter_default(self, parameter_name):
        """
        :param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
        :return value of parameter
        """
        vect_parameters = (
            self.default_parameters.get_limit_coverage_ONT_threshold_default(
                None, None, None, None
            )
        )
        for parameters in vect_parameters:
            if parameters.name == parameter_name:
                return parameters.parameter
        return None

    def is_limit_coverage_ONT_single_parameter_default(
        self, project_sample, parameter_name
    ):
        """
        one possibility available MASK_CONSENSUS_threshold
        """

        value_default_parameter = self.get_limit_coverage_ONT_single_parameter_default(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_limit_coverage_ONT_single_parameter(
            project_sample, parameter_name
        )
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_limit_coverage_ONT_single_parameter(self, project_sample, parameter_name):
        """
        get limit_coverage_ONT single parameters
        :param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
        """

        parameters_string = self.get_limit_coverage_ONT_parameters_all_possibilities(
            project_sample.project.owner, project_sample
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1][1:]
        return None

    def is_limit_coverage_ONT_single_parameter_default_for_project(
        self, project, parameter_name
    ):
        """
        :param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
        """

        value_default_parameter = self.get_limit_coverage_ONT_single_parameter_default(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_limit_coverage_ONT_single_parameter_for_project(
            project, parameter_name
        )
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_limit_coverage_ONT_single_parameter_for_project(
        self, project, parameter_name
    ):
        """
        get limit_coverage_ONT single parameters
        :param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
        """

        parameters_string = self.get_limit_coverage_ONT_parameters_for_project(
            project.owner, project
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1][1:]
        return None

    #####
    #####        END Coverage limit ONT
    #####
    #####################################################

    #####################################################
    #####
    #####        FREQ VCF ONT
    #####

    def get_freq_vcf_ONT_parameters(self, user, type_of_use, project, project_sample):
        """
        get freq_vcf_ONT parameters
        """
        return self.default_parameters.get_parameters(
            SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name,
            user,
            type_of_use,
            project,
            project_sample,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )

    def get_freq_vcf_ONT_parameters_all_possibilities(self, user, project_sample):
        """
        get freq_vcf parameters for project and default
        """

        ### Test project_sample first
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name,
            user,
            Software.TYPE_OF_USE_project_sample,
            None,
            project_sample,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        if not parameters is None:
            return parameters

        ### Test project
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name,
            user,
            Software.TYPE_OF_USE_project,
            project_sample.project,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        if not parameters is None:
            return parameters

        ### can be a default one
        default_software = DefaultSoftware()
        parameters = default_software.get_freq_vcf_ONT_parameters(user)
        if len(parameters) > 0:
            return parameters

        software_names = SoftwareNames()
        return software_names.get_insaflu_parameter_freq_vcf_parameters()

    def get_freq_vcf_ONT_parameters_for_project(self, user, project):
        """
        get freq_vcf_ONT parameters only for project or default
        """

        ### Test project
        parameters = self.default_parameters.get_parameters(
            SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name,
            user,
            Software.TYPE_OF_USE_project,
            project,
            None,
            None,
            ConstantsSettings.TECHNOLOGY_minion,
        )
        if not parameters is None:
            return parameters

        ### can be a default one_get_freq_vcf
        default_software = DefaultSoftware()
        parameters = default_software.get_freq_vcf_ONT_parameters(user)
        if len(parameters) > 0:
            return parameters
        return None

    def get_freq_vcf_ONT_single_parameter_default(self, parameter_name):
        """
        :param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
        :return value of parameter
        """
        vect_parameters = self.default_parameters.get_vcf_freq_ONT_threshold_default(
            None, None, None, None
        )
        for parameters in vect_parameters:
            if parameters.name == parameter_name:
                return parameters.parameter
        return None

    def is_freq_vcf_ONT_single_parameter_default(self, project_sample, parameter_name):
        """
        one possibility available MASK_CONSENSUS_threshold
        """

        value_default_parameter = self.get_freq_vcf_ONT_single_parameter_default(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_freq_vcf_ONT_single_parameter(
            project_sample, parameter_name
        )
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_freq_vcf_ONT_single_parameter(self, project_sample, parameter_name):
        """
        get freq_vcf_ONT single parameters
        :param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
        """

        parameters_string = self.get_freq_vcf_ONT_parameters_all_possibilities(
            project_sample.project.owner, project_sample
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1][1:]
        return None

    def is_freq_vcf_ONT_single_parameter_default_for_project(
        self, project, parameter_name
    ):
        """
        :param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
        """

        value_default_parameter = self.get_freq_vcf_ONT_single_parameter_default(
            parameter_name
        )
        if value_default_parameter is None:
            return False

        parameter_defined = self.get_freq_vcf_ONT_single_parameter_for_project(
            project, parameter_name
        )
        if (
            not parameter_defined is None
            and parameter_defined == value_default_parameter
        ):
            return True
        return False

    def get_freq_vcf_ONT_single_parameter_for_project(self, project, parameter_name):
        """
        get freq_vcf_ONT single parameters
        :param parameter_name -> Only one possibility available MASK_CONSENSUS_threshold
        """

        parameters_string = self.get_freq_vcf_ONT_parameters_for_project(
            project.owner, project
        )
        if parameters_string is None:
            return None
        lst_data = parameters_string.split(parameter_name)
        if len(lst_data) == 2:
            return lst_data[1][1:]
        return None

    #####
    #####        END -- FREQ VCF ONT
    #####
    #####################################################

    def set_default_software(
        self, software, user, type_of_use, project, project_sample, sample
    ):
        """
        Set default parameters for a software
        """
        vect_parameters = self._get_default_parameters(
            software.name,
            user,
            type_of_use,
            project,
            project_sample,
            sample,
            software.technology.name,
        )
        parameters = Parameter.objects.filter(
            software=software,
            project=project,
            project_sample=project_sample,
            sample=sample,
        )
        key_value = "{}_{}_{}".format(
            software.name,
            ConstantsSettings.TECHNOLOGY_illumina
            if software.technology is None
            else software.technology.name,
            user.username,
        )
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

    def is_change_values_for_software(self, software_name, technology_name, user_name):
        """Return if the software has a value changed"""
        key_value = "{}_{}_{}".format(software_name, technology_name, user_name)
        return self.change_values_software.get(key_value, False)

    def can_change_values_for_this_software(
        self, software, project, project_sample, sample, dataset=None
    ):
        """Return True if some of can_change is True"""
        parameters = Parameter.objects.filter(
            software=software,
            project=project,
            project_sample=project_sample,
            sample=sample,
            dataset=dataset,
        )
        for parameter in parameters:
            if parameter.can_change:
                return True
        return False

    def get_parameters(
        self,
        software_name,
        user,
        type_of_use,
        project,
        project_sample,
        sample,
        technology_name=ConstantsSettings.TECHNOLOGY_illumina,
        dataset=None,
        televir_project=None,
    ):
        """ """
        self.test_default_db(
            software_name,
            user,
            type_of_use,
            project,
            project_sample,
            sample,
            technology_name,
            dataset,
        )
        return self.default_parameters.get_parameters(
            software_name,
            user,
            type_of_use,
            project,
            project_sample,
            sample,
            technology_name,
            dataset,
            televir_project=televir_project,
        )

    def get_all_software(self):
        """
        get all softwares available by this class
        """
        vect_software = []
        vect_software.append(self.software_names.get_snippy_name())
        vect_software.append(self.software_names.get_medaka_name_consensus())
        vect_software.append(self.software_names.get_samtools_name_depth_ONT())
        vect_software.append(self.software_names.get_NanoFilt_name())
        vect_software.append(self.software_names.get_trimmomatic_name())
        vect_software.append(
            self.software_names.get_insaflu_parameter_mask_consensus_name()
        )
        vect_software.append(
            self.software_names.get_insaflu_parameter_limit_coverage_name()
        )
        vect_software.append(self.software_names.get_insaflu_parameter_freq_vcf_name())
        vect_software.append(self.software_names.get_abricate_name())
        vect_software.append(self.software_names.get_freebayes_name())
        return vect_software

    def _get_default_project(
        self, user, software_name, project, vect_parameters, technology_name
    ):
        """
        :param software_name name of the software
        :param project is None pass to global
        try to get project parameters
        """
        type_of_use = Software.TYPE_OF_USE_global
        if project is None:
            type_of_use = Software.TYPE_OF_USE_global
            if(software_name==self.software_names.get_abricate_name()):
                type_of_use = Software.TYPE_OF_USE_qc
            if(software_name==self.software_names.get_trimmomatic_name()):
                type_of_use = Software.TYPE_OF_USE_qc 
            if(software_name==self.software_names.get_NanoFilt_name()):
                type_of_use = Software.TYPE_OF_USE_qc   
        elif type(project) is Project:
            type_of_use = Software.TYPE_OF_USE_project
        elif type(project) is ProjectSample:
            type_of_use = Software.TYPE_OF_USE_project
        elif type(project) is Sample:
            type_of_use = Software.TYPE_OF_USE_sample

        try:
            software = Software.objects.get(
                name=software_name,
                owner=user,
                type_of_use=type_of_use,
                technology__name=technology_name,
                version_parameters=self.default_parameters.get_software_parameters_version(
                    software_name
                ),
            )
        except Software.DoesNotExist:
            return vect_parameters

        ## get parameters for a specific user and software
        if project is None:
            parameters = Parameter.objects.filter(software=software)
        elif type(project) is Project:
            parameters = Parameter.objects.filter(software=software, project=project)
        elif type(project) is ProjectSample:
            parameters = Parameter.objects.filter(
                software=software, project_sample=project
            )
        elif type(project) is Sample:
            parameters = Parameter.objects.filter(software=software, sample=project)

        ### parse them
        for parameter in parameters:
            for previous_parameter in vect_parameters:
                if previous_parameter.sequence_out == parameter.sequence_out:
                    previous_parameter.is_to_run = parameter.is_to_run
                    previous_parameter.software.is_to_run = parameter.is_to_run

                    ### don't set the not set parameters
                    if (
                        not parameter.not_set_value is None
                        and parameter.parameter == parameter.not_set_value
                    ):
                        break
                    previous_parameter.parameter = parameter.parameter
                    break
        return vect_parameters
