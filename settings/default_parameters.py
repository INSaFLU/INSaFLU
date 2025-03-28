"""
Created on 10/04/2021

@author: mmp
"""

import logging
import os
from typing import List, Optional

from django.conf import settings
from django.db import DatabaseError, transaction

from constants.meta_key_and_values import MetaKeyAndValue
from constants.software_names import SoftwareNames
from pathogen_identification.constants_settings import (
    ConstantsSettings as PI_ConstantsSettings,
)
from pathogen_identification.utilities.utilities_pipeline import (
    Utility_Pipeline_Manager,
)
from settings.constants_settings import ConstantsSettings
from settings.models import Parameter, PipelineStep, Project, Software, Technology
from utils.lock_atomic_transaction import LockedAtomicTransaction


class DefaultParameters(object):
    """
    Define default values for all softwares
    """

    software_names = SoftwareNames()

    ### used in snippy
    SNIPPY_COVERAGE_NAME = "--mincov"
    SNIPPY_MAPQUAL_NAME = "--mapqual"
    SNIPPY_PRIMER_NAME = "--primer"

    MEDAKA_PRIMER_NAME = "-p"

    ### used in NANOfilt
    NANOfilt_quality_read = "-q"

    ### used in Medaka
    MEDAKA_model = "-m"

    ### used in mask consensus
    MASK_CONSENSUS_threshold = "Threshold"

    ### used when the parameters are passed in other way
    ### parameters values has the META_KEY of the metakey_projectsample/metakey_project
    MASK_DONT_care = "dont_care"
    MASK_not_applicable = "Not applicable"

    ### clean human reads
    MASK_CLEAN_HUMAN_READS = software_names.SOFTWARE_CLEAN_HUMAN_READS_name

    ### MINIMUN of MAX of NanoFilt
    NANOFILT_MINIMUN_MAX = 100

    def __init__(self, prep_televir_dbs: bool = True):
        """
        Constructor
        """
        # televir_util = Parameter_DB_Utility()
        # software_list = televir_util.get_software_list()

        self.televir_db_manager = Utility_Pipeline_Manager()

        # self.televir_db_manager.set_software_list(software_list)

        if prep_televir_dbs:
            self.televir_db_manager.get_software_db_dict()
            self.televir_db_manager.get_host_dbs()

    def get_software_parameters_version(self, software_name):
        """
        With time it is possible to add or remove parameters available in softwares
        Add one every time that is necessary to add or remove parameters in software
        """
        if software_name == SoftwareNames.SOFTWARE_TRIMMOMATIC_name:
            return 1

        if software_name == SoftwareNames.SOFTWARE_FREEBAYES_name:
            return 1
        ## all others remain zero, didn't change anything
        return 0

    def _get_software_obsolete(self):
        """
        Return all softwares/parameters version obsolete
        """
        vect_return = []  ## [[name, version], [name, version], ...]
        vect_return.append([SoftwareNames.SOFTWARE_TRIMMOMATIC_name, 0])
        vect_return.append([SoftwareNames.SOFTWARE_FREEBAYES_name, 0])
        return vect_return

    @transaction.atomic
    def set_software_obsolete(self):
        """
        set software obsolete
        """
        vect_to_obsolete = self._get_software_obsolete()

        for data_obsolete in vect_to_obsolete:
            query_set = Software.objects.filter(
                name=data_obsolete[0],
                version_parameters=data_obsolete[1],
                is_obsolete=False,
            )
            for software in query_set:
                try:
                    with transaction.atomic():
                        software.is_obsolete = True
                except DatabaseError:
                    software.is_obsolete = True
                software.save()

    @transaction.atomic
    def persist_parameters_create(self, vect_parameters: List[Parameter]):
        software = None
        dt_out_sequential = {}
        # with LockedAtomicTransaction(Software), LockedAtomicTransaction(Parameter):
        for parameter in vect_parameters:
            assert parameter.sequence_out not in dt_out_sequential
            if software is None:
                try:
                    software = parameter.software
<<<<<<< HEAD
=======

>>>>>>> 6a32c5a4 (remove some locks, update request ajax identification.)
                    software.save()

                except Exception as e:
                    logging.error("Error persisting software: {}".format(e))
                    continue

            parameter.software = software

            try:

                parameter.save()
                dt_out_sequential[parameter.sequence_out] = 1

            except Exception as e:
                logging.error("Error persisting parameter: {}".format(e))
                continue

    @transaction.atomic
    def persist_parameters(self, vect_parameters: List[Parameter], type_of_use: int):
        """
        persist a specific software by default
        param: type_of_use Can by Software.TYPE_OF_USE_project; Software.TYPE_OF_USE_project_sample
        """
        software = None
        dt_out_sequential = {}

        for parameter in vect_parameters:
            assert parameter.sequence_out not in dt_out_sequential
            if software is None:
                # with LockedAtomicTransaction(Software), LockedAtomicTransaction(
                #    Parameter
                # ):
                try:
                    software = Software.objects.get(
                        name=parameter.software.name,
                        name_extended=parameter.software.name_extended,
                        owner=parameter.software.owner,
                        version=parameter.software.version,
                        type_of_software=parameter.software.type_of_software,
                        type_of_use=parameter.software.type_of_use,
                        technology=parameter.software.technology,
                    )
                except Software.DoesNotExist:
                    software = parameter.software
                    try:
                        # with LockedAtomicTransaction(Software):
                        software = Software.objects.create(
                            name=parameter.software.name,
                            name_extended=parameter.software.name_extended,
                            owner=parameter.software.owner,
                            version=parameter.software.version,
                            type_of_use=parameter.software.type_of_use,
                            technology=parameter.software.technology,
                            version_parameters=parameter.software.version_parameters,
                            pipeline_step=parameter.software.pipeline_step,
                            can_be_on_off_in_pipeline=parameter.software.can_be_on_off_in_pipeline,
                        )
                    except Exception as e:
                        logging.error("Error persisting software: {}".format(e))
                        continue

                except Software.MultipleObjectsReturned:
                    sof = Software.objects.filter(
                        name=parameter.software.name,
                        name_extended=parameter.software.name_extended,
                        owner=parameter.software.owner,
                        type_of_use=parameter.software.type_of_use,
                        technology=parameter.software.technology,
                        version_parameters=parameter.software.version_parameters,
                        pipeline_step=parameter.software.pipeline_step,
                    )

                    # keep last one
                    software = sof.last()

                    print("MULTIPLE SOFTWARES: ", sof.count(), software.name)
                    if sof.count() > 1:
                        sof_delete = sof.exclude(pk=software.pk)
                        Parameter.objects.filter(software__in=sof_delete).delete()
                        sof_delete.delete()

                        # raise Exception("MultipleObjectsReturned")

            parameter.software = software
            try:

                parameter.save()
                dt_out_sequential[parameter.sequence_out] = 1

            except Exception as e:
                logging.error("Error persisting parameter: {}".format(e))
                continue

            ## set sequential number

    def persist_parameters_update(
        self, vect_parameters: List[Parameter], software: Software, range_update=False
    ):
        """
        persist a specific software by default
        param: type_of_use Can by Software.TYPE_OF_USE_project; Software.TYPE_OF_USE_project_sample
        """
        # software = None
        dt_out_sequential = {}
        for parameter in vect_parameters:
            assert parameter.sequence_out not in dt_out_sequential

            try:
                parameter_existing = Parameter.objects.get(
                    software=software,
                    name=parameter.name,
                    sequence_out=parameter.sequence_out,
                )

                if range_update:
                    if parameter.type_data == Parameter.PARAMETER_int:

                        if (
                            parameter.range_available
                            != parameter_existing.range_available
                        ):
                            parameter_existing.range_available = (
                                parameter.range_available
                            )
                            parameter_existing.save()

                        if parameter.range_max != parameter_existing.range_max:
                            parameter_existing.range_max = parameter.range_max
                            parameter_existing.save()

                        if parameter.range_min != parameter_existing.range_min:
                            parameter_existing.range_min = parameter.range_min
                            parameter_existing.save()

                        if parameter.parameter is not None:

                            if int(parameter_existing.parameter) > int(
                                parameter_existing.range_max
                            ):
                                parameter_existing.parameter = (
                                    parameter_existing.range_max
                                )
                                parameter_existing.save()

            except Parameter.DoesNotExist:

                parameter.software = software
                try:

                    parameter.save()
                    dt_out_sequential[parameter.sequence_out] = 1

                except Exception as e:
                    logging.error("Error persisting parameter: {}".format(e))
                    continue

            ## set sequential number

        return len(dt_out_sequential) > 0

    def get_software_global_with_step(
        self,
        user,
        software_name,
        technology_name,
        type_of_use,
        pipeline_step,
        televir_project=None,
    ):
        """
        Get software global using pipeline step."""

        software_list = Software.objects.filter(
            name=software_name,
            owner=user,
            type_of_use=type_of_use,
            technology__name=technology_name,
            version_parameters=self.get_software_parameters_version(software_name),
            pipeline_step__name=pipeline_step,
            parameter__televir_project=televir_project,
        )

        if len(software_list) == 0:
            software_list = Software.objects.filter(
                name=software_name,
                owner=user,
                type_of_use=type_of_use,
                version_parameters=self.get_software_parameters_version(software_name),
                pipeline_step__name=pipeline_step,
                parameter__televir_project=televir_project,
            )

        if len(software_list) == 0:
            return None

        return software_list.first()

    def get_software_project_sample_mdcg_illumina(
        self,
        user,
        # technology_name,
        project=None,
        project_sample=None,
        name_extended=None,
    ) -> Optional[Software]:
        """
        Get software global
        """
        software_list = Software.objects.filter(
            owner=user,
            technology__name__in=[
                ConstantsSettings.TECHNOLOGY_illumina,
                ConstantsSettings.TECHNOLOGY_illumina_old,
            ],
            parameter__project=project,
            parameter__project_sample=project_sample,
            pipeline_step__name=ConstantsSettings.PIPELINE_NAME_variant_detection,
            is_to_run=True,
        ).distinct()

        if len(software_list) == 0:
            software_list = Software.objects.filter(
                owner=user,
                technology__name__in=[
                    ConstantsSettings.TECHNOLOGY_illumina,
                    ConstantsSettings.TECHNOLOGY_illumina_old,
                ],
                parameter__project=project,
                pipeline_step__name=ConstantsSettings.PIPELINE_NAME_variant_detection,
                is_to_run=True,
            ).distinct()

        if len(software_list) == 0:
            software_list = Software.objects.filter(
                owner=user,
                technology__name__in=[
                    ConstantsSettings.TECHNOLOGY_illumina,
                    ConstantsSettings.TECHNOLOGY_illumina_old,
                ],
                parameter__project_sample=project_sample,
                pipeline_step__name=ConstantsSettings.PIPELINE_NAME_variant_detection,
                is_to_run=True,
            ).distinct()

        if name_extended is not None:
            software_list = software_list.filter(name_extended=name_extended)

        if len(software_list) == 0:
            return None

        return software_list.first()

    def get_software_global(
        self,
        user,
        software_name,
        technology_name,
        type_of_use,
        project=None,
        televir_project=None,
        is_to_run=True,
        name_extended=None,
        project_sample=None,
    ) -> Optional[Software]:
        """
        Get software global
        """
        software_list = Software.objects.filter(
            name=software_name,
            owner=user,
            type_of_use=type_of_use,
            technology__name=technology_name,
            version_parameters=self.get_software_parameters_version(software_name),
            parameter__televir_project=televir_project,
            parameter__project=project,
        ).distinct()

        if len(software_list) == 0:
            software_list = Software.objects.filter(
                name=software_name,
                owner=user,
                type_of_use=type_of_use,
                version_parameters=self.get_software_parameters_version(software_name),
                parameter__televir_project=televir_project,
            ).distinct()

        if project_sample is not None:
            software_list = software_list.filter(
                parameter__project_sample=project_sample
            )

        if is_to_run == True:
            software_list = software_list.filter(is_to_run=True)

        if name_extended is not None:
            software_list = software_list.filter(name_extended=name_extended)

        if len(software_list) == 0:
            return None

        return software_list[0]

    def hide_parameter_name_check(self, par_name: str):
        if par_name.startswith("[") and par_name.endswith("]"):
            return True
        return False

    def get_software(
        self,
        software_name,
        user,
        type_of_use,
        project,
        project_sample,
        technology_name=ConstantsSettings.TECHNOLOGY_illumina,
        televir_project=None,
        pipeline_step=None,
        software_name_extended=None,
        is_to_run=False,
    ) -> Optional[Software]:
        if self.check_software_is_polyvalent(software_name):
            if pipeline_step is None:
                prefered_pipeline = self.get_polyvalent_software_pipeline(software_name)
            else:
                prefered_pipeline = pipeline_step

            software = self.get_software_global_with_step(
                user,
                software_name,
                technology_name,
                type_of_use,
                prefered_pipeline,
                televir_project=televir_project,
            )

        elif self.check_software_with_duplicates(software_name):
            software = self.get_software_global(
                user,
                software_name,
                technology_name,
                type_of_use,
                televir_project=televir_project,
                is_to_run=is_to_run,
                name_extended=software_name_extended,
                project=project,
                project_sample=project_sample,
            )
        else:
            software = self.get_software_global(
                user,
                software_name,
                technology_name,
                type_of_use,
                televir_project=televir_project,
                is_to_run=is_to_run,
                name_extended=software_name_extended,
                project=project,
            )

        return software

    def get_parameters_parsed(
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
        pipeline_step=None,
        software_name_extended=None,
        is_to_run=False,
    ):

        parameters = self.get_parameters(
            software_name,
            user,
            type_of_use,
            project,
            project_sample,
            sample,
            technology_name,
            dataset,
            televir_project,
            pipeline_step,
            software_name_extended,
            is_to_run,
        )

        if parameters is None:
            return None

        return self.parse_parameters(parameters, software_name)

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
        pipeline_step=None,
        software_name_extended=None,
        is_to_run=False,
    ):

        software = self.get_software(
            software_name,
            user,
            type_of_use,
            project,
            project_sample,
            technology_name,
            televir_project,
            pipeline_step,
            software_name_extended,
            is_to_run,
        )

        if software is None:
            return software

        if not project_sample is None:
            project = None

        ## logger.debug("Get parameters: software-{} user-{} typeofuse-{} project-{} psample-{} sample-{} tec-{} dataset-{}",software, user, type_of_use, project, project_sample, sample, technology_name, dataset)
        ## get parameters for a specific user  #
        parameters = Parameter.objects.filter(
            software=software,
            project_sample=project_sample,
            sample=sample,
            dataset=dataset,
            televir_project=televir_project,
        ).distinct()

        if not project is None:
            parameters = parameters.filter(project=project)

        # logger.debug("Get parameters: {}".format(parameters))
        ### if only one parameter and it is don't care, return dont_care

        return parameters

    def get_parameters_specific(
        self,
        software_name,
        user,
        type_of_use,
        project: Project,
        project_sample: ProjectSample,
        sample,
        technology_name=ConstantsSettings.TECHNOLOGY_illumina,
        dataset=None,
        televir_project=None,
        pipeline_step=None,
        software_name_extended=None,
        is_to_run=False,
    ):
        """
        get software_name parameters, if it saved in database...
        """
        # logger = logging.getLogger("fluWebVirus.debug")
        # logger.debug("Get parameters: software-{} user-{} typeofuse-{} project-{} psample-{} sample-{} tec-{} dataset-{}",software_name, user, type_of_use, project, project_sample, sample, technology_name, dataset)
        if not project_sample is None:
            project = project_sample.project

        if self.check_software_is_polyvalent(software_name):
            if pipeline_step is None:
                prefered_pipeline = self.get_polyvalent_software_pipeline(software_name)
            else:
                prefered_pipeline = pipeline_step

            software = self.get_software_global_with_step(
                user,
                software_name,
                technology_name,
                type_of_use,
                prefered_pipeline,
                televir_project=televir_project,
            )

        else:
            software = self.get_software_global(
                user,
                software_name,
                technology_name,
                type_of_use,
                televir_project=televir_project,
                is_to_run=is_to_run,
                name_extended=software_name_extended,
                project=project,
                project_sample=project_sample,
            )

        if software is None:
            return software

        if not project_sample is None:
            project = None

        ## logger.debug("Get parameters: software-{} user-{} typeofuse-{} project-{} psample-{} sample-{} tec-{} dataset-{}",software, user, type_of_use, project, project_sample, sample, technology_name, dataset)
        ## get parameters for a specific user  #
        parameters = Parameter.objects.filter(
            software=software,
            project_sample=project_sample,
            sample=sample,
            dataset=dataset,
            televir_project=televir_project,
        ).distinct()

        if not project is None:
            parameters = parameters.filter(project=project)

        # logger.debug("Get parameters: {}".format(parameters))
        ### if only one parameter and it is don't care, return dont_care

        if len(list(parameters)) == 1 and list(parameters)[0].name in [
            DefaultParameters.MASK_not_applicable,
            DefaultParameters.MASK_DONT_care,
        ]:
            return DefaultParameters.MASK_not_applicable

        return self.parse_parameters(parameters, software_name)

    def get_parameters_classic(
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
        pipeline_step=None,
        software_name_extended=None,
        is_to_run=False,
    ):
        """
        get software_name parameters, if it saved in database...
        """
        # logger = logging.getLogger("fluWebVirus.debug")
        # logger.debug("Get parameters: software-{} user-{} typeofuse-{} project-{} psample-{} sample-{} tec-{} dataset-{}",software_name, user, type_of_use, project, project_sample, sample, technology_name, dataset)

        if self.check_software_is_polyvalent(software_name):
            if pipeline_step is None:
                prefered_pipeline = self.get_polyvalent_software_pipeline(software_name)
            else:
                prefered_pipeline = pipeline_step

            software = self.get_software_global_with_step(
                user,
                software_name,
                technology_name,
                type_of_use,
                prefered_pipeline,
                televir_project=televir_project,
            )

        else:

            software = self.get_software_global(
                user,
                software_name,
                technology_name,
                type_of_use,
                televir_project=televir_project,
                is_to_run=is_to_run,
                name_extended=software_name_extended,
                project=project,
            )

        if software is None:
            return software

        if not project_sample is None:
            project = None

        ## logger.debug("Get parameters: software-{} user-{} typeofuse-{} project-{} psample-{} sample-{} tec-{} dataset-{}",software, user, type_of_use, project, project_sample, sample, technology_name, dataset)
        ## get parameters for a specific user  #
        parameters = Parameter.objects.filter(
            software=software,
            project=project,
            project_sample=project_sample,
            sample=sample,
            dataset=dataset,
            televir_project=televir_project,
        )

        if not project is None:
            parameters = parameters.filter(project=project)

        # logger.debug("Get parameters: {}".format(parameters))
        ### if only one parameter and it is don't care, return dont_care

        if len(list(parameters)) == 1 and list(parameters)[0].name in [
            DefaultParameters.MASK_not_applicable,
            DefaultParameters.MASK_DONT_care,
        ]:
            return DefaultParameters.MASK_not_applicable

        return self.parse_parameters(parameters, software_name)

    def edit_paramaters_show(self, parameters: str, software_name) -> str:
        """
        edits for frontend"""

        parameters = str(parameters).strip()

        if software_name == SoftwareNames.SOFTWARE_IRMA_name:
            # find and remove --mincov and value
            parameter_list = parameters.split(" ")
            if "--mincov" in parameter_list:
                index = parameter_list.index("--mincov")
                parameter_list.pop(index)
                parameter_list.pop(index)
            parameters = " ".join(parameter_list)
        return parameters

    def parse_parameters(self, parameters: List[Parameter], software_name) -> str:

        if len(list(parameters)) == 1 and list(parameters)[0].name in [
            DefaultParameters.MASK_not_applicable,
            DefaultParameters.MASK_DONT_care,
        ]:
            return DefaultParameters.MASK_not_applicable

        ### parse them
        dict_out = {}
        vect_order_ouput = []
        for parameter in parameters:
            ### don't set the not set parameters

            if (
                not parameter.not_set_value is None
                and parameter.parameter == parameter.not_set_value
            ):
                continue

            ### create a dict with parameters
            if parameter.name in dict_out:
                dict_out[parameter.name][1].append(parameter.parameter)
                dict_out[parameter.name][0].append(parameter.union_char)
            else:
                vect_order_ouput.append(parameter.name)
                dict_out[parameter.name] = [
                    [parameter.union_char],
                    [parameter.parameter],
                ]

        return_parameter = ""
        for par_name in vect_order_ouput:

            if self.hide_parameter_name_check(par_name) is False:
                return_parameter += " {}".format(par_name)

            if len(dict_out[par_name][1]) == 1 and len(dict_out[par_name][1][0]) == 0:
                continue
            else:
                # return_parameter += " {}".format(par_name)
                ### exception SOFTWARE_TRIMMOMATIC_illuminaclip SOFTWARE_TRIMMOMATIC_name
                if (
                    software_name == SoftwareNames.SOFTWARE_TRIMMOMATIC_name
                    and par_name == SoftwareNames.SOFTWARE_TRIMMOMATIC_illuminaclip
                ):
                    return_parameter += "{}{}{}".format(
                        dict_out[par_name][0][0],
                        dict_out[par_name][1][0],
                        # os.path.join(
                        #    settings.DIR_SOFTWARE,
                        #    "trimmomatic/adapters",
                        #    dict_out[par_name][1][0],
                        # ),
                        SoftwareNames.SOFTWARE_TRIMMOMATIC_addapter_trim_used_to_assemble,
                    )
                elif (
                    software_name == SoftwareNames.SOFTWARE_DIAMOND_name
                    and par_name
                    == SoftwareNames.SOFTWARE_DIAMOND_PARAMETER_SENSITIVITY_name
                ):
                    return_parameter += " {}".format(dict_out[par_name][1][0])
                elif (
                    par_name == DefaultParameters.SNIPPY_PRIMER_NAME
                    and dict_out[par_name][1][0]
                    == SoftwareNames.SOFTWARE_SNIPPY_no_primer
                ):
                    return_parameter += " {}".format(dict_out[par_name][1][0])
                elif par_name == DefaultParameters.SNIPPY_PRIMER_NAME:

                    return_parameter += " {}".format(
                        dict_out[par_name][1][0]
                        # os.path.join(
                        #    settings.DIR_SOFTWARE,
                        #    "trimmomatic/adapters",
                        #    dict_out[par_name][1][0],
                        # )
                    )
                elif (
                    par_name == DefaultParameters.MEDAKA_PRIMER_NAME
                    and software_name == SoftwareNames.SOFTWARE_Medaka_name_consensus
                ):

                    if (
                        dict_out[par_name][1][0]
                        == SoftwareNames.SOFTWARE_SNIPPY_no_primer
                    ):
                        return_parameter += " {}".format(dict_out[par_name][1][0])
                    else:
                        return_parameter += " {}".format(
                            dict_out[par_name][1][0]
                            # os.path.join(
                            #    settings.DIR_SOFTWARE,
                            #    "trimmomatic/adapters",
                            #    dict_out[par_name][1][0],
                            # )
                        )

                elif par_name == "--db":
                    return_parameter += "{}{}".format(
                        dict_out[par_name][0][0],
                        os.path.basename(dict_out[par_name][1][0]),
                    )
                else:
                    for _ in range(len(dict_out[par_name][0])):
                        return_parameter += "{}{}".format(
                            dict_out[par_name][0][_], dict_out[par_name][1][_]
                        )
        # logger.debug("Get parameters return output: {}".format(return_parameter))
        #### This is the case where all the options can be "not to set"

        if len(return_parameter.strip()) == 0 and len(parameters) == 0:
            return None
        return return_parameter.strip()

    def get_list_parameters(
        self,
        software_name,
        user,
        type_of_use,
        project,
        project_sample,
        sample,
        technology_name=ConstantsSettings.TECHNOLOGY_illumina,
        dataset=None,
    ):
        """
        get software_name parameters, if it saved in database...
        """
        try:
            software = Software.objects.get(
                name=software_name,
                owner=user,
                type_of_use=type_of_use,
                technology__name=technology_name,
                version_parameters=self.get_software_parameters_version(software_name),
            )
        except Software.DoesNotExist:
            if type_of_use == Software.TYPE_OF_USE_global:
                try:
                    software = Software.objects.get(
                        name=software_name,
                        owner=user,
                        type_of_use=type_of_use,
                        version_parameters=self.get_software_parameters_version(
                            software_name
                        ),
                    )
                except Software.DoesNotExist:
                    return None
            else:
                return None

        ## get parameters for a specific user
        parameters = Parameter.objects.filter(
            software=software,
            project=project,
            project_sample=project_sample,
            sample=sample,
            dataset=dataset,
        )

        ### if only one parameter and it is don't care, return dont_care
        return list(parameters)

    def is_software_to_run(
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
        """Test if it is necessary to run this software, By default return True"""
        try:
            software = Software.objects.get(
                name=software_name,
                owner=user,
                type_of_use=type_of_use,
                technology__name=technology_name,
                version_parameters=self.get_software_parameters_version(software_name),
            )
        except Software.DoesNotExist:
            if type_of_use == Software.TYPE_OF_USE_global:
                try:
                    software = Software.objects.get(
                        name=software_name,
                        owner=user,
                        type_of_use=type_of_use,
                        version_parameters=self.get_software_parameters_version(
                            software_name
                        ),
                    )
                except Software.DoesNotExist:
                    return True
            else:
                if type_of_use == Software.TYPE_OF_USE_qc:
                    try:
                        software = Software.objects.get(
                            name=software_name,
                            owner=user,
                            type_of_use=type_of_use,
                            version_parameters=self.get_software_parameters_version(
                                software_name
                            ),
                        )
                    except Software.DoesNotExist:
                        return True
                else:
                    return True

        ### if it is Global it is software that is mandatory
        if type_of_use == Software.TYPE_OF_USE_global:
            return software.is_to_run

        ### if it is Global it is software that is mandatory
        if type_of_use == Software.TYPE_OF_USE_qc:
            return software.is_to_run

        ## get parameters for a specific sample, project or project_sample
        parameters = Parameter.objects.filter(
            software=software,
            project=project,
            project_sample=project_sample,
            sample=sample,
            dataset=dataset,
        )

        ### Try to find the parameter of sequence_out == 1. It is the one that has the flag to run or not.
        for parameter in parameters:
            if parameter.sequence_out == 1:
                return parameter.is_to_run
        return True

    def set_software_to_run(
        self,
        software_name,
        user,
        type_of_use,
        project,
        project_sample,
        sample,
        technology_name,
        is_to_run,
        dataset=None,
    ):
        """set software to run ON/OFF
        :output True if the is_to_run is changed"""

        try:
            software = Software.objects.get(
                name=software_name,
                owner=user,
                type_of_use=type_of_use,
                technology__name=technology_name,
                version_parameters=self.get_software_parameters_version(software_name),
            )
        except Software.DoesNotExist:
            if type_of_use == Software.TYPE_OF_USE_global:
                try:
                    software = Software.objects.get(
                        name=software_name,
                        owner=user,
                        type_of_use=type_of_use,
                        version_parameters=self.get_software_parameters_version(
                            software_name
                        ),
                    )
                except Software.DoesNotExist:
                    return False
            else:
                return False

        ## if the software can not be change return False
        if not software.can_be_on_off_in_pipeline:
            return False

        self.set_software_to_run_by_software(
            software,
            project,
            None,
            project_sample,
            sample,
            is_to_run=is_to_run,
            dataset=dataset,
        )
        return True

    def set_software_to_run_by_software(
        self,
        software: Software,
        project,
        televir_project,
        project_sample,
        sample,
        is_to_run=None,
        dataset=None,
    ):
        """set software to run ON/OFF
        :output True if the is_to_run is changed"""

        # with LockedAtomicTransaction(Software), LockedAtomicTransaction(Parameter):
        ## get parameters for a specific sample, project or project_sample
        parameters = Parameter.objects.filter(
            software=software,
            project=project,
            project_sample=project_sample,
            televir_project=televir_project,
            sample=sample,
            dataset=dataset,
        )

        ## if None need to take the value from database
        if is_to_run is None:
            if software.type_of_use in [
                Software.TYPE_OF_USE_qc,
                Software.TYPE_OF_USE_global,
                Software.TYPE_OF_USE_televir_project,
                Software.TYPE_OF_USE_televir_global,
                Software.TYPE_OF_USE_televir_project,
                Software.TYPE_OF_USE_televir_settings,
                Software.TYPE_OF_USE_televir_project_settings,
            ]:
                is_to_run = not software.is_to_run
            elif len(parameters) > 0:
                is_to_run = not parameters[0].is_to_run
            else:
                is_to_run = not software.is_to_run

<<<<<<< HEAD
            ## if the software can not be change return False
            if not software.can_be_on_off_in_pipeline:
                if software.type_of_use in [
                    Software.TYPE_OF_USE_qc,
                    Software.TYPE_OF_USE_global,
                    Software.TYPE_OF_USE_televir_global,
                    Software.TYPE_OF_USE_televir_project,
                    Software.TYPE_OF_USE_televir_settings,
                    Software.TYPE_OF_USE_televir_project_settings,
                ]:
                    return software.is_to_run
                elif len(parameters) > 0:
                    return parameters[0].is_to_run
                return True

        # if it is Global it is software that is mandatory
        # only can change if TYPE_OF_USE_global, other type_of_use is not be tested
        if software.type_of_use in [
            Software.TYPE_OF_USE_qc,
            Software.TYPE_OF_USE_project,
            Software.TYPE_OF_USE_global,
            Software.TYPE_OF_USE_televir_global,
            Software.TYPE_OF_USE_televir_project,
            Software.TYPE_OF_USE_televir_settings,
            Software.TYPE_OF_USE_televir_project_settings,
        ]:
            software.is_to_run = is_to_run
            software.save()

        if (
            software.pipeline_step.name
            == ConstantsSettings.PIPELINE_NAME_variant_detection
        ):
            software.is_to_run = is_to_run
            software.save()

        ## get parameters for a specific sample, project or project_sample

        parameters = Parameter.objects.filter(
            software=software,
            project=project,
            televir_project=televir_project,
            project_sample=project_sample,
            sample=sample,
            dataset=dataset,
        )

<<<<<<< HEAD
        ### Try to find the parameter of sequence_out == 1. It is the one that has the flag to run or not.
        for parameter in parameters:
            parameter.is_to_run = is_to_run
            parameter.save()
=======
            ### Try to find the parameter of sequence_out == 1. It is the one that has the flag to run or not.
            for parameter in parameters:
                parameter.is_to_run = is_to_run
                parameter.save()
>>>>>>> cb39ec32 (merge flumut / insaflu pipelins into ubuntu development branch)
=======
        ## if the software can not be change return False
        if not software.can_be_on_off_in_pipeline:
            if software.type_of_use in [
                Software.TYPE_OF_USE_qc,
                Software.TYPE_OF_USE_global,
                Software.TYPE_OF_USE_televir_global,
                Software.TYPE_OF_USE_televir_project,
                Software.TYPE_OF_USE_televir_settings,
                Software.TYPE_OF_USE_televir_project_settings,
            ]:
                return software.is_to_run
            elif len(parameters) > 0:
                return parameters[0].is_to_run
            return True

        # if it is Global it is software that is mandatory
        # only can change if TYPE_OF_USE_global, other type_of_use is not be tested
        if software.type_of_use in [
            Software.TYPE_OF_USE_qc,
            Software.TYPE_OF_USE_project,
            Software.TYPE_OF_USE_global,
            Software.TYPE_OF_USE_televir_global,
            Software.TYPE_OF_USE_televir_project,
            Software.TYPE_OF_USE_televir_settings,
            Software.TYPE_OF_USE_televir_project_settings,
        ]:
            software.is_to_run = is_to_run
            software.save()

        if (
            software.pipeline_step.name
            == ConstantsSettings.PIPELINE_NAME_variant_detection
        ):
            software.is_to_run = is_to_run
            software.save()

        ## get parameters for a specific sample, project or project_sample

        parameters = Parameter.objects.filter(
            software=software,
            project=project,
            televir_project=televir_project,
            project_sample=project_sample,
            sample=sample,
            dataset=dataset,
        )
>>>>>>> 6a32c5a4 (remove some locks, update request ajax identification.)

        ### Try to find the parameter of sequence_out == 1. It is the one that has the flag to run or not.
        for parameter in parameters:
            parameter.is_to_run = is_to_run
            parameter.save()

        return is_to_run

    def get_vect_parameters(self, software: Software):
        """return all parameters, by software instance"""

        if software.name == SoftwareNames.SOFTWARE_IVAR_name:
            return self.get_ivar_default(
                software.owner,
                Software.TYPE_OF_USE_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            )

        elif software.name == SoftwareNames.SOFTWARE_IRMA_name:
            return self.get_irma_default(
                software.owner,
                Software.TYPE_OF_USE_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            )

        if software.name == SoftwareNames.SOFTWARE_SNIPPY_name:

            return self.get_snippy_default(
                software.owner,
                Software.TYPE_OF_USE_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            )
        elif software.name == SoftwareNames.SOFTWARE_TRIMMOMATIC_name:
            return self.get_trimmomatic_default(
                software.owner,
                Software.TYPE_OF_USE_qc,
                ConstantsSettings.TECHNOLOGY_illumina,
            )
        elif software.name == SoftwareNames.SOFTWARE_NanoFilt_name:
            return self.get_nanofilt_default(
                software.owner,
                Software.TYPE_OF_USE_qc,
                ConstantsSettings.TECHNOLOGY_minion,
            )
        elif software.name == SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name:
            return self.get_mask_consensus_threshold_default(
                software.owner,
                Software.TYPE_OF_USE_global,
                (
                    ConstantsSettings.TECHNOLOGY_illumina
                    if software.technology is None
                    else software.technology.name
                ),
            )
        elif software.name == SoftwareNames.SOFTWARE_CLEAN_HUMAN_READS_name:
            return self.get_clean_human_reads_default(
                software.owner,
                Software.TYPE_OF_USE_global,
                (
                    ConstantsSettings.TECHNOLOGY_illumina
                    if software.technology is None
                    else software.technology.name
                ),
            )
        elif software.name == SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name:
            return self.get_limit_coverage_ONT_threshold_default(
                software.owner,
                Software.TYPE_OF_USE_global,
                ConstantsSettings.TECHNOLOGY_minion,
            )
        elif software.name == SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name:
            return self.get_vcf_freq_ONT_threshold_default(
                software.owner,
                Software.TYPE_OF_USE_global,
                ConstantsSettings.TECHNOLOGY_minion,
            )
        elif software.name == SoftwareNames.SOFTWARE_Medaka_name_consensus:
            return self.get_medaka_model_default(
                software.owner,
                Software.TYPE_OF_USE_global,
                ConstantsSettings.TECHNOLOGY_minion,
            )
        elif software.name == SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT:
            return self.get_samtools_depth_default_ONT(
                software.owner,
                Software.TYPE_OF_USE_global,
                ConstantsSettings.TECHNOLOGY_minion,
            )
        elif software.name == SoftwareNames.SOFTWARE_ABRICATE_name:
            return self.get_abricate_default(
                software.owner,
                Software.TYPE_OF_USE_qc,
                (
                    ConstantsSettings.TECHNOLOGY_illumina
                    if software.technology is None
                    else software.technology.name
                ),
            )
        elif software.name == SoftwareNames.SOFTWARE_MASK_CONSENSUS_BY_SITE_name:
            return self.get_mask_consensus_by_site_default(
                software.owner,
                Software.TYPE_OF_USE_global,
                (
                    ConstantsSettings.TECHNOLOGY_illumina
                    if software.technology is None
                    else software.technology.name
                ),
            )
        elif software.name == SoftwareNames.SOFTWARE_GENERATE_CONSENSUS_name:
            return self.get_generate_consensus_default(software.owner)
        elif software.name == SoftwareNames.SOFTWARE_FREEBAYES_name:
            return self.get_freebayes_default(
                software.owner,
                Software.TYPE_OF_USE_global,
                (
                    ConstantsSettings.TECHNOLOGY_illumina
                    if software.technology is None
                    else software.technology.name
                ),
            )
        elif software.name == SoftwareNames.SOFTWARE_NEXTSTRAIN_name:
            return self.get_nextstrain_default(software.owner)

        ####
        #### PATHOGEN IDENTIFICATION
        ####

        elif software.name == SoftwareNames.SOFTWARE_REMAP_PARAMS_name:
            return self.get_remap_defaults(
                software.owner,
                Software.TYPE_OF_USE_televir_settings,
                software.technology.name,
            )

        elif software.name == SoftwareNames.SOFTWARE_BAMUTIL_name:
            return self.get_bamutil_defaults(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                software.technology.name,
            )

        elif software.name == SoftwareNames.SOFTWARE_DUSTMASKER_name:
            return self.get_dustmasker_defaults(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                software.technology.name,
            )

        elif software.name == SoftwareNames.SOFTWARE_MSAMTOOLS_name:
            return self.get_msamtools_defaults(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                software.technology.name,
            )

        elif software.name == SoftwareNames.SOFTWARE_PRINSEQ_name:
            return self.get_prinseq_defaults(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                software.technology.name,
            )

        elif software.name == SoftwareNames.SOFTWARE_televir_report_layout_name:
            return self.get_televir_report_defaults(
                software.owner,
                Software.TYPE_OF_USE_televir_settings,
                software.technology.name,
            )

        elif software.name == SoftwareNames.SOFTWARE_METAGENOMICS_SETTINGS_name:
            return self.get_metagenomics_settings_defaults(
                software.owner,
                Software.TYPE_OF_USE_televir_settings,
                software.technology.name,
            )

        elif software.name == SoftwareNames.SOFTWARE_CENTRIFUGE_name:
            return self.get_centrifuge_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            )

        elif software.name == SoftwareNames.SOFTWARE_BWA_name:
            return self.get_bwa_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            )

        elif software.name == SoftwareNames.SOFTWARE_KRAKEN2_name:
            return self.get_kraken2_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            )
        elif software.name == SoftwareNames.SOFTWARE_KRAKENUNIQ_name:
            return self.get_krakenuniq_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            )
        elif software.name == SoftwareNames.SOFTWARE_DIAMOND_name:
            return self.get_diamond_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            )
        elif software.name == SoftwareNames.SOFTWARE_KAIJU_name:
            return self.get_kaiju_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_read_classification,
            )

        elif software.name == SoftwareNames.SOFTWARE_DESAMBA_name:
            return self.get_desamba_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            )
        elif software.name == SoftwareNames.SOFTWARE_BLAST_name:
            return self.get_blast_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            )
        elif software.name == SoftwareNames.SOFTWARE_RAVEN_name:
            return self.get_raven_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            )
        elif software.name == SoftwareNames.SOFTWARE_SPAdes_name:
            return self.get_spades_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            )
        elif software.name == SoftwareNames.SOFTWARE_FASTVIROMEEXPLORER_name:
            return self.get_fastviromeexplorer_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            )
        elif software.name == SoftwareNames.SOFTWARE_SNIPPY_PI_name:
            return self.get_snippy_pi_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
                pipeline_step=software.pipeline_step.name,
            )

        elif software.name == SoftwareNames.SOFTWARE_BOWTIE2_DEPLETE_name:
            return self.get_bowtie2_deplete_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            )

        elif software.name == SoftwareNames.SOFTWARE_BOWTIE2_REMAP_name:
            return self.get_bowtie2_remap_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            )

        elif software.name == SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ILLU_name:

            if (
                software.name_extended
                == SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ILLU_name_extended_screening
            ):
                return self.get_minimap2_remap_illumina_default(
                    software.owner,
                    Software.TYPE_OF_USE_televir_global,
                    ConstantsSettings.TECHNOLOGY_illumina,
                    pipeline_step=ConstantsSettings.PIPELINE_NAME_metagenomics_screening,
                )

            return self.get_minimap2_remap_illumina_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
                pipeline_step=ConstantsSettings.PIPELINE_NAME_metagenomics_screening,
            )
        elif software.name == SoftwareNames.SOFTWARE_MINIMAP2_MAP_ASSEMBLY_name:
            return self.get_minimap2_map_assembly_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            )
        elif software.name == SoftwareNames.SOFTWARE_MINIMAP2_DEPLETE_ILLU_name:
            return self.get_minimap2_depletion_illumina_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_illumina,
            )

        elif software.name == SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ONT_name:

            if (
                software.name_extended
                == SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ONT_name_extended_screening
            ):
                return self.get_minimap2_remap_ONT_default(
                    software.owner,
                    Software.TYPE_OF_USE_televir_global,
                    ConstantsSettings.TECHNOLOGY_minion,
                    pipeline_step=ConstantsSettings.PIPELINE_NAME_metagenomics_screening,
                    job="screening",
                )
            elif (
                software.name_extended
                == SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ONT_name_extended_request_mapping
            ):
                return self.get_minimap2_remap_ONT_default(
                    software.owner,
                    Software.TYPE_OF_USE_televir_global,
                    ConstantsSettings.TECHNOLOGY_minion,
                    pipeline_step=ConstantsSettings.PIPELINE_NAME_request_mapping,
                    job="request_mapping",
                )
            else:
                return self.get_minimap2_remap_ONT_default(
                    software.owner,
                    Software.TYPE_OF_USE_televir_global,
                    ConstantsSettings.TECHNOLOGY_minion,
                )

        elif software.name == SoftwareNames.SOFTWARE_MINIMAP2_DEPLETE_ONT_name:
            return self.get_minimap2_depletion_ONT_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            )

        elif software.name == SoftwareNames.SOFTWARE_EMPTY_name:
            return self.get_empty_software_default(
                software.owner,
                Software.TYPE_OF_USE_televir_global,
                ConstantsSettings.TECHNOLOGY_minion,
            )

        else:
            return None

    def check_software_with_duplicates(self, software_name):
        """return true if software is in duplocate_list"""
        return software_name in SoftwareNames.duplicate_softwares

    def check_software_is_polyvalent(self, software_name):
        """return True if the software is polyvalent"""
        if software_name in SoftwareNames.polyvalent_software:
            return True
        else:
            return False

    def get_polyvalent_software_pipeline(self, software_name):
        """return first pipeline where software is used"""
        if software_name in SoftwareNames.polyvalent_software:
            return SoftwareNames.polyvalent_software_pipelines[software_name][0]
        else:
            return None

    def _get_pipeline(self, pipeline_name):
        """return a record for a pipeline step name"""
        if pipeline_name is None:
            return None

        try:
            pipeline_step = PipelineStep.objects.get(name=pipeline_name)
        except PipelineStep.DoesNotExist as e:
            pipeline_step = PipelineStep()
            pipeline_step.name = pipeline_name
            pipeline_step.save()
        return pipeline_step

    def get_technology(self, technology_name) -> Technology:
        """return a record for a pipeline step name"""
        if technology_name is None:
            return None

        try:
            technology = Technology.objects.get(name=technology_name)
        except Technology.DoesNotExist:
            technology = Technology()
            technology.name = technology_name
            technology.save()
        return technology

    def get_snippy_default(
        self,
        user,
        type_of_use,
        technology_name,
        project=None,
        project_sample=None,
        pipeline_step=None,
    ):
        """
        –mapqual: minimum mapping quality to allow (–mapqual 20)
        —mincov: minimum coverage of variant site (–mincov 10)
        –minfrac: minumum proportion for variant evidence (–minfrac 0.51)
        primer: Fasta of amplicon scheme primers for filtering ("")
        """
        if not pipeline_step:
            pipeline_step = ConstantsSettings.PIPELINE_NAME_variant_detection

        software = Software()
        software.name = SoftwareNames.SOFTWARE_SNIPPY_name
        software.name_extended = SoftwareNames.SOFTWARE_SNIPPY_name_extended
        software.version = SoftwareNames.SOFTWARE_SNIPPY_VERSION
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True  ## set to True if it is going to run, for example Trimmomatic can run or not

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(pipeline_step)
        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = DefaultParameters.SNIPPY_MAPQUAL_NAME
        parameter.parameter = "20"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.range_available = "[0:100]"
        parameter.range_max = "100"
        parameter.range_min = "0"
        parameter.description = "MAPQUAL: is the minimum mapping quality to accept in variant calling (–mapqual 20)."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = DefaultParameters.SNIPPY_COVERAGE_NAME
        parameter.parameter = "10"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.sequence_out = 2
        parameter.range_available = "[4:100]"
        parameter.range_max = "100"
        parameter.range_min = "4"
        parameter.description = "MINCOV: the minimum number of reads covering a site to be considered (–mincov 10)."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--minfrac"
        parameter.parameter = "0.51"
        parameter.type_data = Parameter.PARAMETER_float
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.sequence_out = 3
        parameter.range_available = "[0.5:1.0]"
        parameter.range_max = "1.0"
        parameter.range_min = "0.5"
        parameter.description = (
            "MINFRAC: minimum proportion for variant evidence (–minfrac 0.51)"
        )
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = DefaultParameters.SNIPPY_PRIMER_NAME
        parameter.parameter = SoftwareNames.SOFTWARE_SNIPPY_no_primer
        # parameter.not_set_value = SoftwareNames.SOFTWARE_SNIPPY_no_primer
        parameter.not_set_value = "NA"
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.sequence_out = 4
        parameter.description = "PRIMER: fasta of primers used for amplicon sequencing"
        vect_parameters.append(parameter)

        return vect_parameters

    def get_ivar_default(
        self,
        user,
        type_of_use,
        technology_name,
        project=None,
        project_sample=None,
        pipeline_step=None,
    ):
        """
        –mapqual: minimum mapping quality to allow (–mapqual 20)
        —mincov: minimum coverage of variant site (–mincov 10)
        –minfrac: minumum proportion for variant evidence (–minfrac 0.51)
        primer: Fasta of amplicon scheme primers for filtering ("")
        """
        if not pipeline_step:
            pipeline_step = ConstantsSettings.PIPELINE_NAME_variant_detection

        software = Software()
        software.name = SoftwareNames.SOFTWARE_IVAR_name
        software.name_extended = SoftwareNames.SOFTWARE_IVAR_name_extended
        software.version = SoftwareNames.SOFTWARE_IVAR_VERSION
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = False  ## set to True if it is going to run, for example Trimmomatic can run or not

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(pipeline_step)
        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = DefaultParameters.SNIPPY_MAPQUAL_NAME
        parameter.parameter = "20"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.range_available = "[0:100]"
        parameter.range_max = "100"
        parameter.range_min = "0"
        parameter.description = "MAPQUAL: is the minimum mapping quality to accept in variant calling (–mapqual 20)."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = DefaultParameters.SNIPPY_COVERAGE_NAME
        parameter.parameter = "10"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.sequence_out = 2
        parameter.range_available = "[4:100]"
        parameter.range_max = "100"
        parameter.range_min = "4"
        parameter.description = "MINCOV: the minimum number of reads covering a site to be considered (–mincov 10)."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--minfrac"
        parameter.parameter = "0.51"
        parameter.type_data = Parameter.PARAMETER_float
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.sequence_out = 3
        parameter.range_available = "[0.5:1.0]"
        parameter.range_max = "1.0"
        parameter.range_min = "0.5"
        parameter.description = (
            "MINFRAC: minimum proportion for variant evidence (–minfrac 0.51)"
        )
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = DefaultParameters.SNIPPY_PRIMER_NAME
        parameter.parameter = SoftwareNames.SOFTWARE_SNIPPY_no_primer
        # parameter.not_set_value = SoftwareNames.SOFTWARE_SNIPPY_no_primer
        parameter.not_set_value = "NA"
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.sequence_out = 4
        parameter.description = "PRIMER: fasta of primers used for amplicon sequencing"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--ivar"
        parameter.parameter = ""
        parameter.type_data = Parameter.PARAMETER_char
        parameter.software = software
        parameter.project_sample = project_sample
        parameter.project = project
        parameter.union_char = " "
        parameter.can_change = False
        parameter.sequence_out = 5
        parameter.description = "iVar software"
        vect_parameters.append(parameter)

        return vect_parameters

    def get_irma_default(
        self, user, type_of_use, technology_name, project=None, project_sample=None
    ):
        """
        –mapqual: minimum mapping quality to allow (–mapqual 20)
        —mincov: minimum coverage of variant site (–mincov 10)
        –minfrac: minumum proportion for variant evidence (–minfrac 0.51)
        primer: Fasta of amplicon scheme primers for filtering ("")
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_IRMA_name
        software.name_extended = SoftwareNames.SOFTWARE_IRMA_name_extended
        software.version = SoftwareNames.SOFTWARE_IRMA_VERSION
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = False

        ###  small description of software
        software.help_text = "iterative refinement meta-assembler"

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_variant_detection
        )

        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = "--module"
        parameter.parameter = SoftwareNames.SOFTWARE_IRMA_PARAMETER_model_options[0]
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.sequence_out = 1
        parameter.description = "Organism specific assembly configuration"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = DefaultParameters.SNIPPY_COVERAGE_NAME
        parameter.parameter = "1"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.sequence_out = 2
        parameter.range_available = "[1:1]"
        parameter.range_max = "1"
        parameter.range_min = "1"
        parameter.description = "MINCOV: the minimum number of reads covering a site to be considered (–mincov 10)."
        vect_parameters.append(parameter)

        return vect_parameters

    def get_freebayes_default(
        self, user, type_of_use, technology_name, project=None, project_sample=None
    ):
        """
        –min-mapping-quality: excludes read alignments from analysis if they have a mapping quality less than Q (–min-mapping-quality 20)
        —min-base-quality: excludes alleles from iSNV analysis if their supporting base quality is less than Q (–min-base-quality 20)
        –min-coverage: requires at least 100-fold of coverage to process a site (–min-coverage 100)
        —min-alternate-count: require at least 10 reads supporting an alternate allele within a single individual in order to evaluate the position (–min-alternate-count 10)
        –min-alternate-fraction: defines the minimum intra-host frequency of the alternate allele to be assumed (–min-alternate-fraction 0.01). This frequency is contingent on the depth of coverage of each processed site since min-alternate-count is set to 10, i.e., the identification of iSNV sites at frequencies of 10%, 2% and 1% is only allowed for sites with depth of coverage of at least 100-fold, 500-fold and 1000-fold, respectively.
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_FREEBAYES_name
        software.name_extended = SoftwareNames.SOFTWARE_FREEBAYES_name_extended
        software.version = SoftwareNames.SOFTWARE_FREEBAYES_VERSION
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True  ## set to True if it is going to run, for example Trimmomatic can run or not

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_intra_host_minor_variant_detection
        )
        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = "--min-mapping-quality"
        parameter.parameter = "20"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.range_available = "[10:50]"
        parameter.range_max = "50"
        parameter.range_min = "10"
        parameter.not_set_value = "0"
        parameter.description = "min-mapping-quality: excludes read alignments from analysis if they have a mapping quality less than Q (–min-mapping-quality 20)."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--min-base-quality"
        parameter.parameter = "20"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.sequence_out = 2
        parameter.range_available = "[10:50]"
        parameter.range_max = "50"
        parameter.range_min = "10"
        parameter.not_set_value = "0"
        parameter.description = "min-base-quality: excludes alleles from iSNV analysis if their supporting base quality is less than Q (–min-base-quality 20)."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--min-coverage"
        parameter.parameter = "100"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.sequence_out = 3
        parameter.range_available = "[20:500]"
        parameter.range_max = "500"
        parameter.range_min = "20"
        parameter.description = "min-coverage: requires at least 100-fold of coverage to process a site (–min-coverage 100)."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--min-alternate-count"
        parameter.parameter = "10"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.sequence_out = 4
        parameter.range_available = "[5:100]"
        parameter.range_max = "100"
        parameter.range_min = "5"
        parameter.description = (
            "min-alternate-count: require at least 10 reads supporting an alternate allele within a single "
            + "individual in order to evaluate the position (–min-alternate-count 10)."
        )
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--min-alternate-fraction"
        parameter.parameter = "0.01"
        parameter.type_data = Parameter.PARAMETER_float
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.sequence_out = 5
        parameter.range_available = "[0.01:0.5]"
        parameter.range_max = "0.5"
        parameter.range_min = "0.01"
        parameter.description = "min-alternate-fraction: defines the minimum intra-host frequency of the alternate allele to be assumed (–min-alternate-fraction 0.01)."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--ploidy"
        parameter.parameter = "2"
        parameter.type_data = Parameter.PARAMETER_float
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.sequence_out = 6
        parameter.range_available = "[1:5]"
        parameter.range_max = "5"
        parameter.range_min = "1"
        parameter.description = (
            "ploidy: Sets the default ploidy for the analysis to N. (--ploidy 2)."
        )
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "-V"
        parameter.parameter = ""
        parameter.type_data = Parameter.PARAMETER_null
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = ""
        parameter.can_change = False
        parameter.sequence_out = 7
        parameter.description = "binomial-obs-priors-off: Disable incorporation of prior expectations about observations."
        vect_parameters.append(parameter)

        return vect_parameters

    def get_nextstrain_default(self, user, dataset=None):
        """
        build: excludes read alignments from analysis if they have a mapping quality less than Q (–min-mapping-quality 20)
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_NEXTSTRAIN_name
        software.name_extended = SoftwareNames.SOFTWARE_NEXTSTRAIN_name_extended
        # TODO Add version
        software.version = SoftwareNames.SOFTWARE_NEXTSTRAIN_VERSION
        software.type_of_use = Software.TYPE_OF_USE_dataset
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(ConstantsSettings.TECHNOLOGY_generic)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True  ## set to True if it is going to run, for example Trimmomatic can run or not

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = None
        software.owner = user

        vect_parameters = []

        # For the moment only has one parameter...
        parameter = Parameter()
        parameter.name = "build"
        parameter.parameter = SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_parameter
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.dataset = dataset
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.not_set_value = SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_parameter
        parameter.description = "Define the build to be used"

        vect_parameters.append(parameter)

        return vect_parameters

    def get_mask_consensus_threshold_default(
        self, user, type_of_use, technology_name, project=None, project_sample=None
    ):
        """
        Threshold of mask not consensus coverage
        """
        software = Software()
        software.name = SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name
        software.name_extended = (
            SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name_extended
        )
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
        software.version = "1.0"
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            False  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True  ## set to True if it is going to run, for example Trimmomatic can run or not

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_coverage_analysis
        )
        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = DefaultParameters.MASK_CONSENSUS_threshold
        parameter.parameter = "70"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = ":"
        parameter.can_change = True
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.range_available = "[50:100]"
        parameter.range_max = "100"
        parameter.range_min = "50"
        parameter.description = "Minimum percentage of locus horizontal coverage with depth of coverage equal or above –mincov (see Snippy) to generate consensus sequence."
        vect_parameters.append(parameter)
        return vect_parameters

    def get_clean_human_reads_default(self, user, type_of_use, technology_name):
        """
        Threshold of mask not consensus coverage
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_CLEAN_HUMAN_READS_name
        software.name_extended = SoftwareNames.SOFTWARE_CLEAN_HUMAN_READS_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
        software.version = "1.0"
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True  ## set to True if it is going to run, for example Trimmomatic can run or not

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_read_quality_analysis
        )
        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = DefaultParameters.MASK_CLEAN_HUMAN_READS
        parameter.parameter = SoftwareNames.SOFTWARE_TAG_no
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.union_char = ": "
        parameter.can_change = True
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.description = "Clean human reads from fastq files."
        vect_parameters.append(parameter)
        return vect_parameters

    def get_limit_coverage_ONT_threshold_default(
        self, user, type_of_use, technology_name, project=None, project_sample=None
    ):
        """
        Minimum depth of coverage per site to validate the sequence (default: –mincov 30)
        Where to use this cut-off:
                This cut-off is used to exclude from vcf files sites with DEPTH <=30
                This cut-off is used to mask consensus sequences with DEPTH <=30 in msa_masker (-c: 30-1 = 29)
        """
        software = Software()
        software.name = SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name
        software.name_extended = (
            SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name_extended
        )
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
        software.version = "1.0"
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            False  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True  ## set to True if it is going to run, for example Trimmomatic can run or not

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_variant_detection
        )
        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = "Threshold"
        parameter.parameter = "30"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = ":"
        parameter.can_change = True
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.range_available = "[4:100]"
        parameter.range_max = "100"
        parameter.range_min = "4"
        parameter.description = "This cut-off is used to exclude from vcf files sites and to mask consensus sequences. Value in percentage"
        vect_parameters.append(parameter)
        return vect_parameters

    def get_vcf_freq_ONT_threshold_default(
        self, user, type_of_use, technology_name, project=None, project_sample=None
    ):
        """
        MINFRAC: minumum proportion for variant evidence (–minfrac 51) Range: [10:100]

        """
        software = Software()
        software.name = SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name
        software.name_extended = (
            SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name_extended
        )
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
        software.version = "1.0"
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            False  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True  ## set to True if it is going to run, for example Trimmomatic can run or not

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_variant_detection
        )
        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = "Threshold"
        parameter.parameter = "0.80"
        parameter.type_data = Parameter.PARAMETER_float
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = ":"
        parameter.can_change = True
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.range_available = "[0.10:1.0]"
        parameter.range_max = "1.0"
        parameter.range_min = "0.10"
        parameter.description = "MINFRAC: minimum proportion for variant evidence (–minfrac) Range: [0.1:1.0]"
        vect_parameters.append(parameter)
        return vect_parameters

    def get_medaka_model_default(
        self, user, type_of_use, technology_name, project=None, project_sample=None
    ):
        """
        Minimum depth of coverage per site to validate the sequence (default: –mincov 30)
        Where to use this cut-off:
                This cut-off is used to exclude from vcf files sites with DEPTH <=30
                This cut-off is used to mask consensus sequences with DEPTH <=30 in msa_masker (-c: 30-1 = 29)
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_Medaka_name_consensus
        software.name_extended = SoftwareNames.SOFTWARE_Medaka_name_extended_consensus
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_Medaka_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            False  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True  ## set to True if it is going to run, for example Trimmomatic can run or not

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_variant_detection
        )
        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = "-m"
        parameter.parameter = (
            SoftwareNames.SOFTWARE_Medaka_default_model
        )  ## default model
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.description = (
            "Medaka models are named to indicate: "
            + "i) the pore type; "
            + "ii) the sequencing device (min -> MinION, prom -> PromethION); "
            + "iii) the basecaller variant (only high and variant available in INSAFlu); "
            + "iv) the Guppy basecaller version. "
            + "Complete format: "
            + "{pore}_{device}_{caller variant}_{caller version}"
        )
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = DefaultParameters.MEDAKA_PRIMER_NAME
        parameter.parameter = SoftwareNames.SOFTWARE_SNIPPY_no_primer
        # parameter.not_set_value = SoftwareNames.SOFTWARE_SNIPPY_no_primer
        parameter.not_set_value = "NA"
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.sequence_out = 2
        parameter.description = "PRIMER: fasta of primers used for amplicon sequencing"
        vect_parameters.append(parameter)

        return vect_parameters

    def get_bamutil_defaults(
        self,
        user,
        type_of_use,
        technology_name,
        sample=None,
        is_to_run=False,
    ):
        """
        remapping parameters, namely:
            max number of taxids to map against.
            max number of acccids to map for each taxid.
            minimum coverage?
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_BAMUTIL_name
        software.name_extended = SoftwareNames.SOFTWARE_BAMUTIL_name_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
        software.version = SoftwareNames.SOFTWARE_BAMUTIL_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = is_to_run  ## set to True if it is going to run, for example Trimmomatic can run or not

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_remap_filtering
        )

        software.owner = user
        vect_parameters = []

        parameter = Parameter()
        parameter.name = SoftwareNames.SOFTWARE_REMAP_PARAMS_min_quality
        parameter.parameter = "120"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 3
        parameter.range_available = "[40:640]"
        parameter.range_max = "640"
        parameter.range_min = "40"
        parameter.description = "Filter: maximum sum of the mismatch qualities before marking a read unmapped. (Defaults to  120)"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = SoftwareNames.SOFTWARE_REMAP_PARAMS_max_mismatch
        parameter.parameter = ".1"
        parameter.type_data = Parameter.PARAMETER_float
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 4
        parameter.range_available = "[0:1]"
        parameter.range_max = "1"
        parameter.range_min = "0"
        parameter.description = "Soft clipping: maximum fraction of mismatches allowed before clipping from the ends. (Defaults to 0.1)"
        vect_parameters.append(parameter)

        return vect_parameters

    def get_msamtools_defaults(
        self,
        user,
        type_of_use,
        technology_name,
        sample=None,
        is_to_run=False,
        pipeline_step=None,
    ):
        """
        remapping parameters, namely:
            max number of taxids to map against.
            max number of acccids to map for each taxid.
            minimum coverage?
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_MSAMTOOLS_name
        software.name_extended = SoftwareNames.SOFTWARE_MSAMTOOLS_name_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
        software.version = SoftwareNames.SOFTWARE_MSAMTOOLS_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = is_to_run

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        if not pipeline_step:
            pipeline_step = ConstantsSettings.PIPELINE_NAME_remap_filtering
        software.pipeline_step = self._get_pipeline(pipeline_step)

        software.owner = user
        vect_parameters = []

        parameter = Parameter()
        parameter.name = SoftwareNames.SOFTWARE_REMAP_PARAMS_min_length
        parameter.parameter = "100"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.range_available = "[50:500]"
        parameter.range_max = "500"
        parameter.range_min = "50"
        parameter.description = (
            "filter: minimum length of the read to be kept. (Defaults to 100)"
        )
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = SoftwareNames.SOFTWARE_REMAP_PARAMS_min_identity
        parameter.parameter = "0.95"
        parameter.type_data = Parameter.PARAMETER_float
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 2
        parameter.range_available = "[0.5:1.0]"
        parameter.range_max = "1.0"
        parameter.range_min = "0.5"
        parameter.description = "filter: minimum percent identity of the read to be kept. (Defaults to 0.95)"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = SoftwareNames.SOFTWARE_REMAP_PARAMS_min_cover
        parameter.parameter = "0.5"
        parameter.type_data = Parameter.PARAMETER_float
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 3
        parameter.range_available = "[0.5:1.0]"
        parameter.range_max = "1.0"
        parameter.range_min = "0.5"
        parameter.description = "filter: minimum percent coverage of the read to be kept. (Defaults to 0.95)"
        vect_parameters.append(parameter)

        return vect_parameters

    def get_dustmasker_defaults(
        self,
        user,
        type_of_use,
        technology_name,
        sample=None,
        is_to_run=False,
        pipeline_step=None,
    ):

        software = Software()
        software.name = SoftwareNames.SOFTWARE_DUSTMASKER_name
        software.name_extended = SoftwareNames.SOFTWARE_DUSTMASKER_name_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
        software.version = SoftwareNames.SOFTWARE_DUSTMASKER_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = True
        software.is_to_run = is_to_run

        #### small description of software
        software.help_text = "identify and mask low complexity regions in the sequence"

        ### which part of pipeline is going to run
        if pipeline_step == None:
            software.pipeline_step = self._get_pipeline(
                ConstantsSettings.PIPELINE_NAME_map_filtering
            )

        else:
            software.pipeline_step = self._get_pipeline(pipeline_step)

        software.owner = user
        vect_parameters = []

        parameter = Parameter()
        parameter.name = SoftwareNames.SOFTWARE_DUSTMASKER_PARAM_MASK_name
        parameter.parameter = SoftwareNames.SOFTWARE_DUSTMASKER_PARAM_MASK_HARD
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.description = (
            "mask: mask low complexity regions with hard masking"
            "soft: mask low complexity regions with soft masking"
            "both: mask low complexity regions with both hard and soft masking"
        )

        vect_parameters.append(parameter)

        return vect_parameters

    def get_metagenomics_settings_defaults(
        self,
        user,
        type_of_use,
        technology_name,
        sample=None,
        is_to_run=True,
    ):
        """
        remapping parameters, namely:
            max number of taxids to map against.
            max number of acccids to map for each taxid.
            minimum coverage?
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_METAGENOMICS_SETTINGS_name
        software.name_extended = (
            SoftwareNames.SOFTWARE_METAGENOMICS_SETTINGS_name_extended
        )
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
        software.version = SoftwareNames.SOFTWARE_METAGENOMICS_SETTINGS_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = False
        software.is_to_run = is_to_run

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_metagenomics_settings
        )

        software.owner = user
        vect_parameters = []

        parameter = Parameter()
        parameter.name = SoftwareNames.SOFTWARE_COMBINED_min_score  # "min_score"
        parameter.parameter = "0"
        parameter.type_data = Parameter.PARAMETER_float
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = ":"
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.range_available = "[0:1]"
        parameter.range_max = "1"
        parameter.range_min = "0"
        parameter.description = "Minimum standardized score to consider a reference for confirmatory mapping. (Defaults to 0)"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = (
            SoftwareNames.SOFTWARE_COMBINED_include_screening
        )  # "min_length"
        parameter.parameter = "OFF"
        parameter.type_data = Parameter.PARAMETER_char_list

        parameter.software = software
        parameter.sample = sample
        parameter.union_char = ":"
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 2
        parameter.description = (
            "Include screening step in the pipeline. (Defaults to OFF)"
        )
        vect_parameters.append(parameter)

        return vect_parameters

    def get_remap_defaults(
        self,
        user,
        type_of_use,
        technology_name,
        sample=None,
        is_to_run=True,
    ):
        """
        remapping parameters, namely:
            max number of taxids to map against.
            max number of acccids to map for each taxid.
            minimum coverage?
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_REMAP_PARAMS_name
        software.name_extended = SoftwareNames.SOFTWARE_REMAP_PARAMS_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
        software.version = SoftwareNames.SOFTWARE_REMAP_PARAMS_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            False  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = is_to_run  ## set to True if it is going to run, for example Trimmomatic can run or not

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_reporting
        )
        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = SoftwareNames.SOFTWARE_REMAP_PARAMS_max_taxids
        parameter.parameter = "12"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.range_available = "[5:30]"
        parameter.range_max = "30"
        parameter.range_min = "5"
        parameter.description = "Maximum number of taxids to map against"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = SoftwareNames.SOFTWARE_REMAP_PARAMS_max_accids
        parameter.parameter = "4"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 2
        parameter.range_available = "[1:30]"
        parameter.range_max = "30"
        parameter.range_min = "1"
        parameter.description = "Number of accession IDs to map against."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = SoftwareNames.SOFTWARE_REMAP_PARAMS_include_manual
        parameter.parameter = "ON"
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 3
        parameter.description = (
            "Include manual curation in the pipeline. (Defaults to OFF)"
        )
        vect_parameters.append(parameter)

        return vect_parameters

    def get_prinseq_defaults(
        self, user, type_of_use, technology_name, sample=None, is_to_run=False
    ):
        """
        -lc_entropy <float>  Filter on a minimum entropy score. Range: [0.0:1.0].
        -lc_dust <float>     Filter on a maximum DUST score. Range: [0.0:1.0].
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_PRINSEQ_name
        software.name_extended = SoftwareNames.SOFTWARE_PRINSEQ_name_extended
        software.version = SoftwareNames.SOFTWARE_PRINSEQ_VERSION
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = is_to_run  ## set to True if it is going to run, for example Trimmomatic can run or not

        ###  small description of software
        software.help_text = (
            "Prinseq is a program for processing and filtering sequencing data."
        )

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_extra_qc
        )
        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = SoftwareNames.SOFTWARE_PRINSEQ_lc_entropy
        parameter.parameter = "0.5"
        parameter.type_data = Parameter.PARAMETER_float
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.range_available = "[0.0:1.0]"
        parameter.range_max = "1.0"
        parameter.range_min = "0.0"
        parameter.description = "Filter on a minimum entropy score."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = SoftwareNames.SOFTWARE_PRINSEQ_lc_dust
        parameter.parameter = "0.7"
        parameter.type_data = Parameter.PARAMETER_float
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 2
        parameter.range_available = "[0.0:1.0]"
        parameter.range_max = "1.0"
        parameter.range_min = "0.0"
        parameter.description = "Filter on a maximum DUST score."
        vect_parameters.append(parameter)

        return vect_parameters

    def get_televir_report_defaults(
        self, user, type_of_use, technology_name, sample=None, is_to_run=True
    ):
        """
        flag calculations to use for mapping. defaults for viruses, bacteria, probes.
        """

        software = Software()
        software.name = SoftwareNames.SOFTWARE_televir_report_layout_name
        software.name_extended = (
            SoftwareNames.SOFTWARE_televir_report_layout_name_extended
        )
        software.version = SoftwareNames.SOFTWARE_televir_report_layout_version
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_INSAFLU_PARAMETER
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )

        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = False
        software.is_to_run = is_to_run

        ###  small description of software
        software.help_text = "flag calculations to use for mapping. defaults for viruses, bacteria, probes."

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_reporting
        )

        software.owner = user

        vect_parameters = []

        parameter = Parameter()

        parameter.name = SoftwareNames.SOFTWARE_televir_report_layout_flag_name
        parameter.parameter = "viruses"
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.description = "flag calculations to use for mapping. defaults for viruses, bacteria, probes."

        vect_parameters.append(parameter)

        parameter = Parameter()

        parameter.name = SoftwareNames.SOFTWARE_televir_report_layout_threshold_name
        parameter.parameter = "0.5"
        parameter.type_data = Parameter.PARAMETER_float
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.range_available = "[0.0:1.0]"
        parameter.range_max = "1.0"
        parameter.range_min = "0.0"
        parameter.sequence_out = 2
        parameter.description = (
            "read sharing threshold used to group references into groups"
        )
        vect_parameters.append(parameter)

        return vect_parameters

    def get_nanofilt_default(self, user, type_of_use, technology_name, sample=None):
        """
        -l <LENGTH>, Filter on a minimum read length. Range: [50:1000].
        --maxlength Filter on a maximum read length
        dá para colocar outro parametro, por exemplo: -ml <MAX_LENGTH>, Set the maximum read length.
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_NanoFilt_name
        software.name_extended = SoftwareNames.SOFTWARE_NanoFilt_name_extended
        software.version = SoftwareNames.SOFTWARE_NanoFilt_VERSION
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True  ## set to True if it is going to run, for example Trimmomatic can run or not

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_read_quality_analysis
        )
        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = "-q"
        parameter.parameter = "8"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.range_available = "[5:30]"
        parameter.range_max = "30"
        parameter.range_min = "5"
        parameter.description = (
            "-q <QUALITY>, Filter on a minimum average read quality score."
        )
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "-l"
        parameter.parameter = "50"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.sequence_out = 2
        parameter.range_available = "[50:1000]"
        parameter.range_max = "1000"
        parameter.range_min = "50"
        parameter.description = "-l <LENGTH>, Filter on a minimum read length."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--headcrop"
        parameter.parameter = "30"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.sequence_out = 3
        parameter.range_available = "[1:1000]"
        parameter.range_max = "1000"
        parameter.range_min = "1"
        parameter.not_set_value = "0"
        parameter.description = (
            "--headcrop <HEADCROP>, Trim n nucleotides from start of read."
        )
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--tailcrop"
        parameter.parameter = "30"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.sequence_out = 4
        parameter.range_available = "[1:1000]"
        parameter.range_max = "1000"
        parameter.range_min = "1"
        parameter.not_set_value = "0"
        parameter.description = (
            "--tailcrop <TAILCROP>, Trim n nucleotides from end of read."
        )
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--maxlength"
        parameter.parameter = "50000"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.sequence_out = 5
        parameter.range_available = "[{}:50000]".format(
            DefaultParameters.NANOFILT_MINIMUN_MAX
        )
        parameter.range_max = "50000"
        parameter.range_min = "0"
        parameter.not_set_value = "0"
        parameter.description = "--maxlength <LENGTH>, Set a maximum read length."
        vect_parameters.append(parameter)

        return vect_parameters

    def get_samtools_depth_default_ONT(
        self, user, type_of_use, technology_name, project=None, project_sample=None
    ):
        """
        samtools depth for ONT
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_SAMTOOLS_name_depth_ONT
        software.name_extended = SoftwareNames.SOFTWARE_SAMTOOLS_name_extended_depth_ONT
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_SAMTOOLS_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            False  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True  ## set to True if it is going to run, for example Trimmomatic can run or not

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_coverage_analysis
        )

        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = "-q"
        parameter.parameter = "0"  ## default model
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.range_available = "[0:100]"
        parameter.range_max = "100"
        parameter.range_min = "0"
        parameter.not_set_value = "0"
        parameter.description = "-q <Quality> base quality threshold."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "-Q"
        parameter.parameter = "0"  ## default model
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.sequence_out = 2
        parameter.range_available = "[0:100]"
        parameter.range_max = "100"
        parameter.range_min = "0"
        parameter.not_set_value = "0"
        parameter.description = "-Q <Quality> mapping quality threshold."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "-aa"
        parameter.parameter = ""  ## default model
        parameter.type_data = Parameter.PARAMETER_null
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = ""
        parameter.can_change = False
        parameter.sequence_out = 3
        parameter.description = (
            "Output absolutely all positions, including unused ref. sequences."
        )
        vect_parameters.append(parameter)

        return vect_parameters

    def get_trimmomatic_default(self, user, type_of_use, technology_name, sample=None):
        software = Software()
        software.name = SoftwareNames.SOFTWARE_TRIMMOMATIC_name
        software.name_extended = SoftwareNames.SOFTWARE_TRIMMOMATIC_name_extended
        software.version = SoftwareNames.SOFTWARE_TRIMMOMATIC_VERSION
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True  ## set to True if it is going to run, for example Trimmomatic can run or not

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_read_quality_analysis
        )
        software.owner = user

        vect_parameters = []

        # ILLUMINACLIP:/usr/local/software/bioinformatics/trimmomatic/adapters/adapters.fa:3:30:10:6:true
        parameter = Parameter()
        parameter.name = SoftwareNames.SOFTWARE_TRIMMOMATIC_illuminaclip
        parameter.parameter = SoftwareNames.SOFTWARE_TRIMMOMATIC_addapter_not_apply
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.not_set_value = SoftwareNames.SOFTWARE_TRIMMOMATIC_addapter_not_apply
        parameter.union_char = ":"
        parameter.can_change = True
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.description = (
            "To clip the Illumina adapters or PCR primers from the input file using the adapter / primer sequences.\n"
            + "ILLUMINACLIP:<ADAPTER_FILE>:3:30:10:6:true"
        )
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "HEADCROP"
        parameter.parameter = "0"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = ":"
        parameter.can_change = True
        parameter.sequence_out = 2
        parameter.range_available = "[0:100]"
        parameter.range_max = "100"
        parameter.range_min = "0"
        parameter.not_set_value = "0"
        parameter.description = "HEADCROP:<length> Cut the specified number of bases from the start of the read."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "CROP"
        parameter.parameter = "0"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = ":"
        parameter.can_change = True
        parameter.sequence_out = 3
        parameter.range_available = "[0:400]"
        parameter.range_max = "400"
        parameter.range_min = "0"
        parameter.not_set_value = "0"
        parameter.description = "CROP:<length> Cut the read to a specified length."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "SLIDINGWINDOW"
        parameter.parameter = "5"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = ":"
        parameter.can_change = True
        parameter.sequence_out = 4
        parameter.range_available = "[3:50]"
        parameter.range_max = "50"
        parameter.range_min = "3"
        parameter.description = (
            "SLIDINGWINDOW:<windowSize> specifies the number of bases to average across"
        )
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "SLIDINGWINDOW"
        parameter.parameter = "20"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = ":"
        parameter.can_change = True
        parameter.sequence_out = 5
        parameter.range_available = "[10:100]"
        parameter.range_max = "100"
        parameter.range_min = "10"
        parameter.description = (
            "SLIDINGWINDOW:<requiredQuality> specifies the average quality required"
        )
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "LEADING"
        parameter.parameter = "3"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = ":"
        parameter.can_change = True
        parameter.sequence_out = 6
        parameter.range_available = "[0:100]"
        parameter.range_max = "100"
        parameter.range_min = "0"
        parameter.not_set_value = "0"
        parameter.description = (
            "LEADING:<quality> Remove low quality bases from the beginning."
        )
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "TRAILING"
        parameter.parameter = "3"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = ":"
        parameter.can_change = True
        parameter.sequence_out = 7
        parameter.range_available = "[0:100]"
        parameter.range_max = "100"
        parameter.range_min = "0"
        parameter.not_set_value = "0"
        parameter.description = (
            "TRAILING:<quality> Remove low quality bases from the end."
        )
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "MINLEN"
        parameter.parameter = "35"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = ":"
        parameter.can_change = True
        parameter.sequence_out = 8
        parameter.range_available = "[5:500]"
        parameter.range_max = "500"
        parameter.range_min = "5"
        parameter.description = "MINLEN:<length> This module removes reads that fall below the specified minimal length."
        vect_parameters.append(parameter)

        ##        Only available in 0.30 version
        #
        #         parameter = Parameter()
        #         parameter.name = "AVGQUAL"
        #         parameter.parameter = "0"
        #         parameter.type_data = Parameter.PARAMETER_int
        #         parameter.software = software
        #         parameter.sample = sample
        #         parameter.union_char = ":"
        #         parameter.can_change = True
        #         parameter.sequence_out = 8
        #         parameter.range_available = "[0:100]"
        #         parameter.range_max = "100"
        #         parameter.range_min = "0"
        #         parameter.not_set_value = "0"
        #         parameter.description = "AVGQUAL:<quality> Drop the read if the average quality is below the specified level."
        #         vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "TOPHRED33"
        parameter.parameter = ""
        parameter.type_data = Parameter.PARAMETER_null
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = ""
        parameter.can_change = False
        parameter.sequence_out = 9
        parameter.description = (
            "This (re)encodes the quality part of the FASTQ file to base 33."
        )
        vect_parameters.append(parameter)
        return vect_parameters

    def get_centrifuge_default(
        self,
        user,
        type_of_use,
        technology_name,
        sample=None,
        pipeline_step="",
        is_to_run=True,
    ):
        """
        centrifuge default illumina
        """
        if not pipeline_step:
            pipeline_step = ConstantsSettings.PIPELINE_NAME_viral_enrichment

        software = Software()
        software.name = SoftwareNames.SOFTWARE_CENTRIFUGE_name
        software.name_extended = SoftwareNames.SOFTWARE_CENTRIFUGE_name_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_CENTRIFUGE_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = is_to_run

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(pipeline_step)

        software.owner = user

        ### software db
        dbs_available = self.televir_db_manager.software_dbs_dict.get(
            software.name.lower(), ["None"]
        )

        vect_parameters = []

        parameter = Parameter()

        parameter = Parameter()
        parameter.name = "--min-hitlen"
        parameter.parameter = "20"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True if PI_ConstantsSettings.METAGENOMICS else False
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.range_available = "[16:30]"
        parameter.range_max = "30"
        parameter.range_min = "16"
        parameter.description = "Minimum length of partial hits."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "-k"
        parameter.parameter = "5"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True if PI_ConstantsSettings.METAGENOMICS else False
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 2
        parameter.range_available = "[1:5]"
        parameter.range_max = "5"
        parameter.range_min = "1"
        parameter.description = "report up to k distinc assignments per read or pair."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--db"
        parameter.parameter = dbs_available[0]
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 4
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "Database to use"

        vect_parameters.append(parameter)

        return vect_parameters

    def get_bowtie2_deplete_default(
        self, user, type_of_use, technology_name, sample=None, pipeline_step=""
    ):
        """
        bowtie2 default illumina
        """
        if not pipeline_step:
            pipeline_step = ConstantsSettings.PIPELINE_NAME_host_depletion

        software = Software()
        software.name = SoftwareNames.SOFTWARE_BOWTIE2_DEPLETE_name
        software.name_extended = SoftwareNames.SOFTWARE_BOWTIE2_DEPLETE_name_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_BOWTIE2_DEPLETE_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = False

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(pipeline_step)

        software.owner = user

        ### software db
        dbs_available = self.televir_db_manager.get_from_host_db(
            software.name.lower(), ["None"]
        )
        vect_parameters = []

        parameter = Parameter()
        parameter.name = "--db"
        parameter.parameter = dbs_available[0][0]
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "Database to use"
        vect_parameters.append(parameter)

        return vect_parameters

    def get_bowtie2_remap_default(
        self,
        user,
        type_of_use,
        technology_name,
        sample=None,
        pipeline_step="",
        job="confirmatory",
    ):
        if not pipeline_step:
            pipeline_step = ConstantsSettings.PIPELINE_NAME_remapping

        software = Software()
        software.name = SoftwareNames.SOFTWARE_BOWTIE2_REMAP_name

        if job == "confirmatory":
            software.name_extended = SoftwareNames.SOFTWARE_BOWTIE2_REMAP_name_extended
        if job == "screening":
            software.name_extended = (
                SoftwareNames.SOFTWARE_BOWTIE2_REMAP_name_extended_screening
            )
        elif job == "request_mapping":
            software.name_extended = (
                SoftwareNames.SOFTWARE_BOWTIE2_REMAP_name_extended_request_mapping
            )

        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_BOWTIE2_REMAP_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)

        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )

        software.is_to_run = False  ## set to True if it is going to run, for example Trimmomatic can run or not

        # if job == "request_mapping":
        #    software.can_be_on_off_in_pipeline = False
        #    software.is_to_run = True
        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(pipeline_step)

        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = "--no-unal"
        parameter.parameter = ""
        parameter.type_data = Parameter.PARAMETER_null
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = ""
        parameter.can_change = False
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.description = "suppress SAM records for unaligned reads."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "[mode]"
        parameter.parameter = "--end-to-end"
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 2
        parameter.description = "end-to-end alignment (default) or local alignment."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "[preset]"
        parameter.parameter = "--sensitive"
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 3
        parameter.description = "preset to use (default: --sensitive)."
        vect_parameters.append(parameter)

        return vect_parameters

    def get_minimap2_remap_ONT_default(
        self,
        user,
        type_of_use,
        technology_name,
        sample=None,
        pipeline_step="",
        is_to_run=True,
        job="confirmatory",
    ):
        """
        minimap remap ONT default
        """
        if not pipeline_step:
            pipeline_step = ConstantsSettings.PIPELINE_NAME_remapping

        software = Software()
        software.name = SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ONT_name

        if job == "confirmatory":
            software.name_extended = (
                SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ONT_name_extended
            )
        elif job == "screening":
            software.name_extended = (
                SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ONT_name_extended_screening
            )

        elif job == "request_mapping":
            software.name_extended = (
                SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ONT_name_extended_request_mapping
            )

        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ONT_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )

        software.is_to_run = is_to_run

        if job in ["request_mapping", "screening"]:
            software.can_be_on_off_in_pipeline = False
            software.is_to_run = True

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(pipeline_step)

        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = "-ax map-ont"
        parameter.parameter = ""
        parameter.type_data = Parameter.PARAMETER_char
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.description = "preset for ONT data"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--secondary=no"
        parameter.parameter = ""
        parameter.type_data = Parameter.PARAMETER_char
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 2
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "do not output secondary alignments, default no"
        vect_parameters.append(parameter)

        return vect_parameters

    def get_minimap2_remap_illumina_default(
        self,
        user,
        type_of_use,
        technology_name,
        sample=None,
        pipeline_step="",
        is_to_run=False,
    ):
        """
        minimap remap illlumina default
        """

        if not pipeline_step:
            pipeline_step = ConstantsSettings.PIPELINE_NAME_remapping

        software = Software()
        software.name = SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ILLU_name

        if pipeline_step == ConstantsSettings.PIPELINE_NAME_metagenomics_screening:
            software.name_extended = (
                SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ILLU_name_extended_screening
            )
        else:
            software.name_extended = (
                SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ILLU_name_extended
            )

        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_MINIMAP2_REMAP_ILLU_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = is_to_run

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(pipeline_step)

        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = "-ax sr"
        parameter.parameter = ""
        parameter.type_data = Parameter.PARAMETER_char
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.description = "preset for Illumina data"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--secondary=no"
        parameter.parameter = ""
        parameter.type_data = Parameter.PARAMETER_char
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 2
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "do not output secondary alignments, default no"
        vect_parameters.append(parameter)

        return vect_parameters

    def get_minimap2_depletion_ONT_default(
        self, user, type_of_use, technology_name, sample=None, pipeline_step=""
    ):
        """
        minimap remap ONT default
        """
        if not pipeline_step:
            pipeline_step = ConstantsSettings.PIPELINE_NAME_host_depletion

        software = Software()
        software.name = SoftwareNames.SOFTWARE_MINIMAP2_DEPLETE_ONT_name
        software.name_extended = (
            SoftwareNames.SOFTWARE_MINIMAP2_DEPLETE_ONT_name_extended
        )
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_MINIMAP2_DEPLETE_ONT_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = False

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(pipeline_step)

        software.owner = user

        vect_parameters = []

        ### software db
        dbs_available = self.televir_db_manager.get_from_host_db(
            software.name.lower(), [["None", "None"]]
        )

        parameter = Parameter()
        parameter.name = "--db"
        parameter.parameter = dbs_available[0][0]
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "Database to use"

        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "-x map-ont"
        parameter.parameter = ""
        parameter.type_data = Parameter.PARAMETER_char
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 2
        parameter.description = "preset for ONT data"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--secondary=no"
        parameter.parameter = ""
        parameter.type_data = Parameter.PARAMETER_char
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 3
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "do not output secondary alignments, default no"
        vect_parameters.append(parameter)

        return vect_parameters

    def get_minimap2_depletion_illumina_default(
        self, user, type_of_use, technology_name, sample=None, pipeline_step=""
    ):
        """
        minimap remap illlumina default
        """

        if not pipeline_step:
            pipeline_step = ConstantsSettings.PIPELINE_NAME_host_depletion

        software = Software()
        software.name = SoftwareNames.SOFTWARE_MINIMAP2_DEPLETE_ILLU_name
        software.name_extended = (
            SoftwareNames.SOFTWARE_MINIMAP2_DEPLETE_ILLU_name_extended
        )
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_MINIMAP2_DEPLETE_ILLU_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = False

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(pipeline_step)

        software.owner = user

        vect_parameters = []

        ### software db
        dbs_available = self.televir_db_manager.get_from_host_db(
            software.name.lower(), ["None"]
        )

        parameter = Parameter()
        parameter.name = "--db"
        parameter.parameter = dbs_available[0][0]
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "Database to use"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "-x sr"
        parameter.parameter = ""
        parameter.type_data = Parameter.PARAMETER_char
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 2
        parameter.description = "preset for Illumina data"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--secondary=no"
        parameter.parameter = ""
        parameter.type_data = Parameter.PARAMETER_char
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 3
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "do not output secondary alignments, default no"
        vect_parameters.append(parameter)

        return vect_parameters

    def get_empty_software_default(
        self, user, type_of_use, technology_name, sample=None, pipeline_step=""
    ):
        """
        empty software default
        """
        if not pipeline_step:
            pipeline_step = ConstantsSettings.PIPELINE_NAME_extra_qc

        software = Software()
        software.name = SoftwareNames.SOFTWARE_EMPTY_name
        software.name_extended = SoftwareNames.SOFTWARE_EMPTY_name_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_EMPTY_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)

        software.can_be_on_off_in_pipeline = (
            False  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )

        software.is_to_run = False

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(pipeline_step)
        software.owner = user

        vect_parameters = []

        return vect_parameters

    def get_minimap2_map_assembly_default(
        self, user, type_of_use, technology_name, sample=None, pipeline_step=""
    ):
        """
        minimap assembly default
        """
        if not pipeline_step:
            pipeline_step = ConstantsSettings.PIPELINE_NAME_contig_classification

        software = Software()
        software.name = SoftwareNames.SOFTWARE_MINIMAP2_MAP_ASSEMBLY_name
        software.name_extended = (
            SoftwareNames.SOFTWARE_MINIMAP2_MAP_ASSEMBLY_name_extended
        )
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_MINIMAP2_MAP_ASSEMBLY_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)

        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )

        software.is_to_run = False

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(pipeline_step)

        software.owner = user

        vect_parameters = []

        dbs_available = self.televir_db_manager.get_from_software_db_dict(
            software_name=software.name, empty=["None"]
        )

        parameter = Parameter()
        parameter.name = "--db"
        parameter.parameter = dbs_available[0]
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "Database to use"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "-x"
        parameter.parameter = "asm10"
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 2
        parameter.description = (
            "preset for assembly to ref mapping, sequence divergence."
        )
        vect_parameters.append(parameter)

        return vect_parameters

    def get_kraken2_default(
        self,
        user,
        type_of_use,
        technology_name,
        sample=None,
        pipeline_step="",
        is_to_run=True,
    ):
        """
        kraken2 default
        """

        if not pipeline_step:
            pipeline_step = ConstantsSettings.PIPELINE_NAME_read_classification

        software = Software()
        software.name = SoftwareNames.SOFTWARE_KRAKEN2_name
        software.name_extended = SoftwareNames.SOFTWARE_KRAKEN2_name_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_KRAKEN2_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = is_to_run

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(pipeline_step)

        software.owner = user

        ### software db
        dbs_available = self.televir_db_manager.software_dbs_dict.get(
            software.name.lower(), ["None"]
        )

        vect_parameters = []

        parameter = Parameter()

        parameter = Parameter()
        parameter.name = "--quick"
        parameter.parameter = "ON"
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True if PI_ConstantsSettings.METAGENOMICS else False
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.description = "operation mode"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--confidence"
        parameter.parameter = "0.5"
        parameter.type_data = Parameter.PARAMETER_float
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True if PI_ConstantsSettings.METAGENOMICS else False
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 2
        parameter.range_available = "[0.4:1.0]"
        parameter.range_max = "1.0"
        parameter.range_min = "0.4"
        parameter.description = "confidence threshold for reporting a taxon"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--db"
        parameter.parameter = dbs_available[0]
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 4
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "Database to use"

        vect_parameters.append(parameter)

        return vect_parameters

    def get_kaiju_default(
        self, user, type_of_use, technology_name, sample=None, pipeline_step=""
    ):
        """
        kaiju default
        """
        if not pipeline_step:
            pipeline_step = ConstantsSettings.PIPELINE_NAME_viral_enrichment

        software = Software()
        software.name = SoftwareNames.SOFTWARE_KAIJU_name
        software.name_extended = SoftwareNames.SOFTWARE_KAIJU_name_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_KAIJU_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(pipeline_step)

        software.owner = user

        ### software db
        dbs_available = self.televir_db_manager.software_dbs_dict.get(
            software.name.lower(), ["None"]
        )

        vect_parameters = []

        parameter = Parameter()

        parameter = Parameter()
        parameter.name = "-e"
        parameter.parameter = "3"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True if PI_ConstantsSettings.METAGENOMICS else False
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.range_available = "[1:5]"
        parameter.range_max = "5"
        parameter.range_min = "1"
        parameter.description = "Number of mismatches allowed"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "-s"
        parameter.parameter = "65"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True if PI_ConstantsSettings.METAGENOMICS else False
        parameter.is_to_run = True
        parameter.sequence_out = 2
        parameter.range_available = "[60:80]"
        parameter.range_max = "80"
        parameter.range_min = "60"
        parameter.description = "Minimum match score in Greedy mode"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "-X"
        parameter.parameter = ""
        parameter.type_data = Parameter.PARAMETER_null
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 3
        parameter.description = "Disable SEG low complexity filter"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--db"
        parameter.parameter = dbs_available[0]
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 4
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "Database to use"

        vect_parameters.append(parameter)

        return vect_parameters

    def get_diamond_default(self, user, type_of_use, technology_name, sample=None):
        """
        diamond default
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_DIAMOND_name
        software.name_extended = SoftwareNames.SOFTWARE_DIAMOND_name_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_DIAMOND_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_read_classification
        )

        software.owner = user

        vect_parameters = []

        dbs_available = self.televir_db_manager.software_dbs_dict.get(
            software.name.lower(), ["None"]
        )
        vect_parameters = []

        parameter = Parameter()
        parameter.name = SoftwareNames.SOFTWARE_DIAMOND_PARAMETER_SENSITIVITY_name
        parameter.parameter = "--fast"
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.description = "sensitivity mode"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--id"
        parameter.parameter = "40"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 3
        parameter.range_available = "[30:80]"
        parameter.range_max = "80"
        parameter.range_min = "40"
        parameter.description = "Minimum percent identity"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--query-cover"
        parameter.parameter = "40"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 4
        parameter.range_available = "[20:60]"
        parameter.range_max = "60"
        parameter.range_min = "20"
        parameter.description = "Minimum query coverage"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--db"
        parameter.parameter = dbs_available[0]
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 5
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "Database to use"

        vect_parameters.append(parameter)

        return vect_parameters

    def get_krakenuniq_default(
        self, user, type_of_use, technology_name, sample=None, pipeline_step=""
    ):
        """
        krakenuniq default
        """

        if not pipeline_step:
            pipeline_step = ConstantsSettings.PIPELINE_NAME_read_classification

        software = Software()
        software.name = SoftwareNames.SOFTWARE_KRAKENUNIQ_name
        software.name_extended = SoftwareNames.SOFTWARE_KRAKENUNIQ_name_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_KRAKENUNIQ_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(pipeline_step)

        software.owner = user

        ### dbs
        dbs_available = self.televir_db_manager.software_dbs_dict.get(
            software.name.lower(), ["None"]
        )

        vect_parameters = []

        parameter = Parameter()

        parameter = Parameter()
        parameter.name = "--hll-precision"
        parameter.parameter = "14"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 2
        parameter.range_available = "[10:18]"
        parameter.range_max = "18"
        parameter.range_min = "13"
        parameter.description = (
            "Precision for HyperLogLog k-mer cardinality estimation."
        )
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--exact"
        parameter.parameter = ""
        parameter.type_data = Parameter.PARAMETER_char
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 3
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "Use exact k-mer counting (slower, but more accurate)."
        vect_parameters.append(parameter)

        vect_parameters = []

        parameter = Parameter()
        parameter.name = "--db"
        parameter.parameter = dbs_available[0]
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 5
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "Database to use"

        vect_parameters.append(parameter)

        return vect_parameters

    def get_bwa_default(
        self, user, type_of_use, technology_name, sample=None, pipeline_step=""
    ):
        if not pipeline_step:
            pipeline_step = ConstantsSettings.PIPELINE_NAME_host_depletion

        software = Software()
        software.name = SoftwareNames.SOFTWARE_BWA_name
        software.name_extended = SoftwareNames.SOFTWARE_BWA_name_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_BWA_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = False

        ###  small description of software
        software.help_text = ""
        software.pipeline_step = self._get_pipeline(pipeline_step)

        software.owner = user

        ### software db
        dbs_available = self.televir_db_manager.get_from_host_db(
            software.name.lower(), ["None"]
        )

        vect_parameters = []

        parameter = Parameter()
        parameter.name = "--db"
        parameter.parameter = dbs_available[0][0]
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "Database to use"

        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "-M"
        parameter.parameter = ""
        parameter.type_data = Parameter.PARAMETER_null
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = ""
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 2
        parameter.description = (
            "Mark shorter split hits as secondary (for Picard compatibility)."
        )
        vect_parameters.append(parameter)

        return vect_parameters

    def get_blast_default(self, user, type_of_use, technology_name, sample=None):
        """
        blast default
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_BLAST_name
        software.name_extended = SoftwareNames.SOFTWARE_BLAST_name_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_BLAST_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_contig_classification
        )

        ### dbs available
        dbs_available = self.televir_db_manager.software_dbs_dict.get(
            software.name.lower(), ["None"]
        )

        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = "-max_target_seqs"
        parameter.parameter = "5"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.range_available = "[1:10]"
        parameter.range_max = "10"
        parameter.range_min = "1"
        parameter.description = (
            "Maximum number of target sequences to report alignments for"
        )
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "-evalue"
        parameter.parameter = "1e-5"
        parameter.type_data = Parameter.PARAMETER_char
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 2
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "E-value threshold for inclusion in the output"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--db"
        parameter.parameter = dbs_available[0]
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 3
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "Database to use"

        vect_parameters.append(parameter)

        return vect_parameters

    def get_fastviromeexplorer_default(
        self, user, type_of_use, technology_name, sample=None
    ):
        """
        fastviromeExplorer default
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_FASTVIROMEEXPLORER_name
        software.name_extended = SoftwareNames.SOFTWARE_FASTVIROMEEXPLORER_name_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_FASTVIROMEEXPLORER_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = False

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_read_classification
        )

        software.owner = user

        vect_parameters = []

        dbs_available = self.televir_db_manager.software_dbs_dict.get(
            software.name.lower(), ["None"]
        )

        prefered = "virosaurus.idx"
        db_prefered = [x for x in dbs_available if prefered in x]
        if db_prefered:
            dbs_available.remove(db_prefered[0])
            dbs_available.insert(0, db_prefered[0])

        vect_parameters = []

        parameter = Parameter()
        parameter.name = "--db"
        parameter.parameter = dbs_available[0]
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "Database to use"

        vect_parameters.append(parameter)

        parameter = Parameter()

        parameter = Parameter()
        parameter.name = "-cr"
        parameter.parameter = "0.3"
        parameter.type_data = Parameter.PARAMETER_float
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 2
        parameter.range_available = "[0.2:0.4]"
        parameter.range_max = "0.4"
        parameter.range_min = "0.2"
        parameter.description = "value of ratio criteria, default 0.3"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "-co"
        parameter.parameter = "0.1"
        parameter.type_data = Parameter.PARAMETER_float
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 3
        parameter.range_available = "[0.1:0.4]"
        parameter.range_max = "0.4"
        parameter.range_min = "0.1"
        parameter.description = "value of overlap criteria, default 0.1"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "-cn"
        parameter.parameter = "10"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 4
        parameter.range_available = "[8:14]"
        parameter.range_max = "14"
        parameter.range_min = "8"
        parameter.description = "value of number of contigs criteria, default 10"
        vect_parameters.append(parameter)

        return vect_parameters

    def get_desamba_default(self, user, type_of_use, technology_name, sample=None):
        """
        desamba default
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_DESAMBA_name
        software.name_extended = SoftwareNames.SOFTWARE_DESAMBA_name_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_DESAMBA_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_read_classification
        )

        software.owner = user

        vect_parameters = []

        dbs_available = self.televir_db_manager.software_dbs_dict.get(
            software.name.lower(), ["None"]
        )
        vect_parameters = []

        parameter = Parameter()
        parameter.name = "--db"
        parameter.parameter = dbs_available[0]
        parameter.type_data = Parameter.PARAMETER_char_list
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "Database to use"

        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = ""
        parameter.parameter = ""
        parameter.type_data = Parameter.PARAMETER_char
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = False
        parameter.sequence_out = 2
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = ""
        vect_parameters.append(parameter)

        return vect_parameters

    def get_spades_default(self, user, type_of_use, technology_name, sample=None):
        """
        spades default
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_SPAdes_name
        software.name_extended = SoftwareNames.SOFTWARE_SPAdes_name_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_SPAdes_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_assembly
        )

        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = "-k"
        parameter.parameter = "31"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.range_available = "[31:41]"
        parameter.range_max = "41"
        parameter.range_min = "31"
        parameter.description = "value of k-mer size, default 31"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--phred-offset"
        parameter.parameter = "33"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 3
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "value of phred-offset, default 33"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--only-assembler"
        parameter.parameter = ""
        parameter.type_data = Parameter.PARAMETER_char
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 4
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "skip read error correction step"
        vect_parameters.append(parameter)

        return vect_parameters

    def get_raven_default(self, user, type_of_use, technology_name, sample=None):
        """
        raven default
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_RAVEN_name
        software.name_extended = SoftwareNames.SOFTWARE_RAVEN_name_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_RAVEN_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_assembly
        )

        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = "-f"
        parameter.parameter = "0.001"
        parameter.type_data = Parameter.PARAMETER_float
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = (
            "minimum fraction of reads required to support a variant"
        )
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "-w"
        parameter.parameter = "5"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 2
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "minimum number of reads required to support a variant"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "-k"
        parameter.parameter = "21"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 3
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = "k-mer size"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "-p"
        parameter.parameter = "2"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 4
        parameter.range_available = "[1:4]"
        parameter.range_max = "4"
        parameter.range_min = "1"
        parameter.description = "number of times racon is invoked, default 2"
        vect_parameters.append(parameter)

        return vect_parameters

    def get_snippy_pi_default(
        self, user, type_of_use, technology_name, sample=None, pipeline_step=None
    ):
        """
        snippy for televir mapping
        –mapqual: minimum mapping quality to allow (–mapqual 20)
        —mincov: minimum coverage of variant site (–mincov 10)
        –minfrac: minumum proportion for variant evidence (–minfrac 0.51)
        """
        if not pipeline_step:
            pipeline_step = ConstantsSettings.PIPELINE_NAME_remapping

        else:
            pipeline_step = pipeline_step

        software = Software()
        software.name = SoftwareNames.SOFTWARE_SNIPPY_PI_name
        software.name_extended = SoftwareNames.SOFTWARE_SNIPPY_PI_name_extended
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = SoftwareNames.SOFTWARE_SNIPPY_VERSION
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run; NEED TO CHECK
        software.pipeline_step = self._get_pipeline(pipeline_step)

        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = "--mapqual"
        parameter.parameter = "20"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = True
        parameter.is_to_run = True
        parameter.sequence_out = 1
        parameter.range_available = "[20:60]"
        parameter.range_max = "60"
        parameter.range_min = "20"
        parameter.description = "minimum mapping quality, default 20"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--mincov"
        parameter.parameter = "2"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True
        parameter.sequence_out = 2
        parameter.range_available = "[2:10]"
        parameter.range_max = "10"
        parameter.range_min = "2"
        parameter.description = "minimum coverage, default 10"
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--minfrac"
        parameter.parameter = "0.51"
        parameter.type_data = Parameter.PARAMETER_float
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.sequence_out = 3
        parameter.range_available = "[0.5:1.0]"
        parameter.range_max = "1.0"
        parameter.range_min = "0.5"
        parameter.description = (
            "MINFRAC: minimum proportion for variant evidence (–minfrac 0.51)"
        )
        vect_parameters.append(parameter)
        return vect_parameters

    ##############################
    ############################## END OF PI SOFTWARE DEFAULTS.

    def get_abricate_default(self, user, type_of_use, technology_name, sample=None):
        software = Software()
        software.name = SoftwareNames.SOFTWARE_ABRICATE_name
        software.name_extended = SoftwareNames.SOFTWARE_ABRICATE_name_extended
        software.version = SoftwareNames.SOFTWARE_ABRICATE_VERSION
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_type_and_subtype_analysis
        )
        software.owner = user

        vect_parameters = []
        parameter = Parameter()
        parameter.name = "--minid"
        parameter.parameter = "70"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.range_available = "[1:100]"
        parameter.range_max = "100"
        parameter.range_min = "1"
        parameter.description = "Minimum DNA %identity."
        vect_parameters.append(parameter)

        parameter = Parameter()
        parameter.name = "--mincov"
        parameter.parameter = "40"
        parameter.type_data = Parameter.PARAMETER_int
        parameter.software = software
        parameter.sample = sample
        parameter.union_char = " "
        parameter.can_change = False
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 2
        parameter.range_available = "[0:100]"
        parameter.range_max = "100"
        parameter.range_min = "0"
        parameter.description = "Minimum DNA %coverage."
        vect_parameters.append(parameter)
        return vect_parameters

    def get_mask_consensus_by_site_default(
        self, user, type_of_use, technology_name, project=None, project_sample=None
    ):
        """
        Mask consensus by site
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_MASK_CONSENSUS_BY_SITE_name
        software.name_extended = (
            SoftwareNames.SOFTWARE_MASK_CONSENSUS_BY_SITE_name_extended
        )
        software.type_of_use = type_of_use
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = "1.0"
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(technology_name)
        software.can_be_on_off_in_pipeline = (
            False  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True  ## set to True if it is going to run, for example Trimmomatic can run or not

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_variant_detection
        )
        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = DefaultParameters.MASK_DONT_care
        parameter.parameter = MetaKeyAndValue.META_KEY_Masking_consensus
        parameter.type_data = Parameter.PARAMETER_char
        parameter.software = software
        parameter.project = project
        parameter.project_sample = project_sample
        parameter.union_char = ""
        parameter.can_change = False
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = ""
        vect_parameters.append(parameter)
        return vect_parameters

    def get_generate_consensus_default(self, user, project_sample=None):
        """
        Mask consensus by site
        """
        software = Software()
        software.name = SoftwareNames.SOFTWARE_GENERATE_CONSENSUS_name
        software.name_extended = SoftwareNames.SOFTWARE_GENERATE_CONSENSUS_name_extended
        software.type_of_use = Software.TYPE_OF_USE_project_sample
        software.type_of_software = Software.TYPE_SOFTWARE
        software.version = "1.0"
        software.version_parameters = self.get_software_parameters_version(
            software.name
        )
        software.technology = self.get_technology(ConstantsSettings.TECHNOLOGY_generic)
        software.can_be_on_off_in_pipeline = (
            True  ## set to True if can be ON/OFF in pipeline, otherwise always ON
        )
        software.is_to_run = True  ## set to True if it is going to run, for example Trimmomatic can run or not

        ###  small description of software
        software.help_text = ""

        ###  which part of pipeline is going to run
        software.pipeline_step = self._get_pipeline(
            ConstantsSettings.PIPELINE_NAME_variant_detection
        )
        software.owner = user

        vect_parameters = []

        parameter = Parameter()
        parameter.name = DefaultParameters.MASK_DONT_care
        parameter.parameter = DefaultParameters.MASK_DONT_care
        parameter.type_data = Parameter.PARAMETER_none
        parameter.software = software
        parameter.project_sample = project_sample
        parameter.union_char = ""
        parameter.can_change = False
        parameter.is_to_run = True  ### by default it's True
        parameter.sequence_out = 1
        parameter.range_available = ""
        parameter.range_max = ""
        parameter.range_min = ""
        parameter.description = ""
        vect_parameters.append(parameter)
        return vect_parameters


### "Generate consensus" -> it is used for set ON/OFF consensus in the AllConsensus File


### "Generate consensus" -> it is used for set ON/OFF consensus in the AllConsensus File
