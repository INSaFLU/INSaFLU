"""
Created on Dec 6, 2017

@author: mmp
"""
from constants.meta_key_and_values import MetaKeyAndValue
from constants.software_names import SoftwareNames
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt, csrf_protect
from extend_user.models import Profile
from managing_files.manage_database import ManageDatabase
from managing_files.models import Project, ProjectSample, Sample
from pathogen_identification.constants_settings import Pipeline_Makeup
from pathogen_identification.models import Projects as PIProjects
from utils.process_SGE import ProcessSGE
from utils.result import DecodeObjects, MaskingConsensus

from settings.constants_settings import ConstantsSettings
from settings.default_parameters import DefaultParameters
from settings.default_software import DefaultSoftware
from settings.default_software_project_sample import DefaultProjectSoftware
from settings.models import Parameter, Software


@csrf_protect
def set_default_parameters(request):
    """
    remove a reference. It can only be removed if not belongs to any deleted project
    """
    if request.is_ajax():
        data = {"is_ok": False}
        software_id_a = "software_id"
        project_id_a = "project_id"
        project_sample_id_a = "project_sample_id"
        sample_id_a = "sample_id"

        ## some pre-requisites
        if not request.user.is_active or not request.user.is_authenticated:
            return JsonResponse(data)
        try:
            profile = Profile.objects.get(user__pk=request.user.pk)
        except Profile.DoesNotExist:
            return JsonResponse(data)
        if profile.only_view_project:
            return JsonResponse(data)
        project, project_sample = None, None

        if software_id_a in request.GET:
            software_id = request.GET[software_id_a]
            b_change = False
            try:
                software = Software.objects.get(pk=software_id)
                if project_id_a in request.GET:
                    project_id = request.GET[project_id_a]
                    project = Project.objects.get(pk=project_id)

                    if (
                        software.name
                        == SoftwareNames.SOFTWARE_MASK_CONSENSUS_BY_SITE_name
                    ):
                        masking_consensus_original = get_fresh_masking_consensus(
                            project.reference
                        )
                        if change_mask_consensus_in_project(
                            project, masking_consensus_original
                        ):
                            b_change = True
                        else:
                            b_change = False
                        data["default"] = "Has no values, yet"
                    else:
                        default_project_software = DefaultProjectSoftware()
                        default_project_software.set_default_software(
                            software,
                            request.user,
                            Software.TYPE_OF_USE_project,
                            project,
                            None,
                            None,
                        )
                        ## set a new default
                        data["default"] = default_project_software.get_parameters(
                            software.name,
                            request.user,
                            Software.TYPE_OF_USE_project,
                            project,
                            None,
                            None,
                            software.technology.name,
                        )

                        b_change = (
                            default_project_software.is_change_values_for_software(
                                software.name,
                                software.technology.name,
                                request.user.username,
                            )
                        )
                elif project_sample_id_a in request.GET:
                    project_sample_id = request.GET[project_sample_id_a]
                    project_sample = ProjectSample.objects.get(pk=project_sample_id)

                    if (
                        software.name
                        == SoftwareNames.SOFTWARE_MASK_CONSENSUS_BY_SITE_name
                    ):
                        masking_consensus_original = get_fresh_masking_consensus(
                            project_sample.project.reference
                        )
                        if change_mask_consensus_in_project_sample(
                            project_sample, masking_consensus_original
                        ):
                            b_change = True
                        else:
                            b_change = False
                        data["default"] = "Has no values, yet"
                    else:
                        default_project_software = DefaultProjectSoftware()
                        default_project_software.set_default_software(
                            software,
                            request.user,
                            Software.TYPE_OF_USE_project_sample,
                            None,
                            project_sample,
                            None,
                        )
                        ## set a new default
                        data["default"] = default_project_software.get_parameters(
                            software.name,
                            request.user,
                            Software.TYPE_OF_USE_project_sample,
                            None,
                            project_sample,
                            None,
                            software.technology.name,
                        )

                        ### need to re-run this sample with snippy if the values change
                        if default_project_software.is_change_values_for_software(
                            software.name,
                            ConstantsSettings.TECHNOLOGY_illumina
                            if project_sample.is_sample_illumina()
                            else ConstantsSettings.TECHNOLOGY_minion,
                            request.user.username,
                        ):
                            b_change = True

                            ### re-run data
                            metaKeyAndValue = MetaKeyAndValue()
                            manageDatabase = ManageDatabase()
                            process_SGE = ProcessSGE()

                            ### change flag to not finished
                            project_sample.is_finished = False
                            project_sample.save()

                            ### create a task to perform the analysis of snippy and freebayes
                            try:
                                (
                                    job_name_wait,
                                    job_name,
                                ) = request.user.profile.get_name_sge_seq(
                                    Profile.SGE_PROCESS_dont_care, Profile.SGE_GLOBAL
                                )
                                if project_sample.is_sample_illumina():
                                    taskID = process_SGE.set_second_stage_snippy(
                                        project_sample,
                                        request.user,
                                        job_name,
                                        [job_name_wait],
                                    )
                                else:
                                    taskID = process_SGE.set_second_stage_medaka(
                                        project_sample,
                                        request.user,
                                        job_name,
                                        [job_name_wait],
                                    )

                                ### set project sample queue ID
                                manageDatabase.set_project_sample_metakey(
                                    project_sample,
                                    request.user,
                                    metaKeyAndValue.get_meta_key_queue_by_project_sample_id(
                                        project_sample.id
                                    ),
                                    MetaKeyAndValue.META_VALUE_Queue,
                                    taskID,
                                )

                                ### need to collect global files again
                                taskID = process_SGE.set_collect_global_files(
                                    project, request.user
                                )
                                manageDatabase.set_project_metakey(
                                    project,
                                    request.user,
                                    metaKeyAndValue.get_meta_key(
                                        MetaKeyAndValue.META_KEY_Queue_TaskID_Project,
                                        project.id,
                                    ),
                                    MetaKeyAndValue.META_VALUE_Queue,
                                    taskID,
                                )
                            except:
                                pass

                elif sample_id_a in request.GET:
                    sample_id = request.GET[sample_id_a]
                    sample = Sample.objects.get(pk=sample_id)

                    ### can not do anything because the sample is running
                    if sample.is_sample_in_the_queue:
                        data = {"is_ok": False}
                        data[
                            "message"
                        ] = "You cannot do this operation. The sample '{}' is in pipeline to run.".format(
                            sample.name
                        )
                        return JsonResponse(data)

                    default_project_software = DefaultProjectSoftware()
                    default_project_software.set_default_software(
                        software,
                        request.user,
                        Software.TYPE_OF_USE_sample,
                        None,
                        None,
                        sample,
                    )
                    ## set a new default
                    data["default"] = default_project_software.get_parameters(
                        software.name,
                        request.user,
                        Software.TYPE_OF_USE_sample,
                        None,
                        None,
                        sample,
                        software.technology.name,
                    )

                    ### need to re-run this sample with NanoFilt if the values change
                    if default_project_software.is_change_values_for_software(
                        software.name,
                        ConstantsSettings.TECHNOLOGY_illumina
                        if sample.is_type_fastq_gz_sequencing()
                        else ConstantsSettings.TECHNOLOGY_minion,
                        request.user.username,
                    ):
                        b_change = True

                        ### re-run data
                        manageDatabase = ManageDatabase()
                        process_SGE = ProcessSGE()

                        ### create a task to perform the analysis of NanoFilt
                        try:
                            (
                                job_name_wait,
                                job_name,
                            ) = request.user.profile.get_name_sge_seq(
                                Profile.SGE_PROCESS_clean_sample, Profile.SGE_SAMPLE
                            )
                            if sample.is_type_fastq_gz_sequencing():
                                taskID = process_SGE.set_run_trimmomatic_species(
                                    sample, request.user, job_name
                                )
                            else:
                                taskID = process_SGE.set_run_clean_minion(
                                    sample, request.user, job_name
                                )

                            ### set sample queue ID
                            manageDatabase.set_sample_metakey(
                                sample,
                                sample.owner,
                                MetaKeyAndValue.META_KEY_Queue_TaskID,
                                MetaKeyAndValue.META_VALUE_Queue,
                                taskID,
                            )
                        except:
                            sample.is_sample_in_the_queue = False
                            sample.save()
                            data["message"] = "Error in the queue system."
                            return JsonResponse(data)

                        ## refresh sample list for this user
                        if not job_name is None:
                            process_SGE.set_create_sample_list_by_user(
                                request.user, [job_name]
                            )
                else:
                    default_software = DefaultSoftware()
                    default_software.set_default_software(software)

                    ## set a new default
                    data["default"] = default_software.get_parameters(
                        software.name, request.user, software.technology.name
                    )
                    b_change = default_software.is_change_values_for_software(
                        software, request.user.username
                    )

                ### clean values in mask site consensus
                if software.name == SoftwareNames.SOFTWARE_MASK_CONSENSUS_BY_SITE_name:
                    if project is None and not project_sample is None:
                        project = project_sample.project

                    if not project is None and b_change:
                        ## check if they have projects
                        count = ProjectSample.objects.filter(
                            project=project,
                            is_deleted=False,
                            is_error=False,
                            is_finished=True,
                        ).count()
                        if (
                            count > 0
                        ):  ### need to send a message to recalculate the global files
                            metaKeyAndValue = MetaKeyAndValue()
                            manageDatabase = ManageDatabase()
                            try:
                                process_SGE = ProcessSGE()
                                taskID = process_SGE.set_collect_global_files(
                                    project, request.user
                                )
                                manageDatabase.set_project_metakey(
                                    project,
                                    request.user,
                                    metaKeyAndValue.get_meta_key(
                                        MetaKeyAndValue.META_KEY_Queue_TaskID_Project,
                                        project.id,
                                    ),
                                    MetaKeyAndValue.META_VALUE_Queue,
                                    taskID,
                                )
                                data["is_ok"] = True
                                data["message"] = " clean all sites."
                            except:
                                data = {"is_ok": False}
                    else:
                        data["message"] = " already has no sites to mask."
                else:
                    ### message to show
                    if b_change:
                        data["message"] = "were set to default values."
                    else:
                        data["message"] = "already had the default values."

            except Software.DoesNotExist:
                return JsonResponse(data)
            except Project.DoesNotExist:
                return JsonResponse(data)
            except ProjectSample.DoesNotExist:
                return JsonResponse(data)
            except Sample.DoesNotExist:
                return JsonResponse(data)
            data["is_ok"] = True
        return JsonResponse(data)


## @csrf_protect
@csrf_exempt
def mask_consensus(request):
    """
    mask consensus
    """
    if request.is_ajax():
        data = {"is_ok": False}

        ## some pre-requisites
        if not request.user.is_active or not request.user.is_authenticated:
            return JsonResponse(data)
        try:
            profile = Profile.objects.get(user__pk=request.user.pk)
        except Profile.DoesNotExist:
            return JsonResponse(data)
        if profile.only_view_project:
            return JsonResponse(data)

        all_data_a = "all_data"
        project_id_a = "project_id"
        project_sample_id_a = "project_sample_id"

        project_id = None
        project_sample_id = None
        if project_id_a in request.POST:
            project_id = request.POST[project_id_a]
        elif project_sample_id_a in request.POST:
            project_sample_id = request.POST[project_sample_id_a]
        project, project_sample = None, None
        b_change_data = False

        manageDatabase = ManageDatabase()
        genetic_element = DecodeObjects()
        ## only for project
        if not project_id is None:
            try:
                project = Project.objects.get(pk=project_id)
            except ProjectSample.DoesNotExist:
                return JsonResponse(data)

            masking_consensus_proposed = genetic_element.decode_result(
                request.POST[all_data_a]
            )
            masking_consensus_proposed.cleaning_mask_results()  ## clean data

            if change_mask_consensus_in_project(project, masking_consensus_proposed):
                b_change_data = True
                data[
                    "message"
                ] = "The project '{}' is going to mask/unmask consensus. ".format(
                    project.name
                )
                data["is_going_to_mask"] = True
            else:
                data[
                    "message"
                ] = "Masking regions are the same, nothing to do for project '{}'".format(
                    project.name
                )
                data["is_going_to_mask"] = False

            data[
                "new_title_i"
            ] = masking_consensus_proposed.get_message_mask_to_show_in_web_site()
            data[
                "new_class_i"
            ] = "padding-button-table {} fa fa-2x fa-pencil padding-button-table tip".format(
                "warning_fa_icon"
                if masking_consensus_proposed.has_masking_data()
                else ""
            )
            data["default"] = (
                "Has positions masked"
                if masking_consensus_proposed.has_masking_data()
                else "Has no values, yet"
            )
        ## only for project sample
        if not project_sample_id is None:  ## for project sample
            try:
                project_sample = ProjectSample.objects.get(pk=project_sample_id)
            except ProjectSample.DoesNotExist:
                return JsonResponse(data)

            masking_consensus_proposed = genetic_element.decode_result(
                request.POST[all_data_a]
            )
            masking_consensus_proposed.cleaning_mask_results()  ## clean data

            if change_mask_consensus_in_project_sample(
                project_sample, masking_consensus_proposed
            ):
                b_change_data = True
                data[
                    "message"
                ] = "The project sample is going to mask/unmask consensus. "
                data["is_going_to_mask"] = True
            else:
                data[
                    "message"
                ] = "Masking regions are the same, nothing to do for project sample"
                data["is_going_to_mask"] = False

            data[
                "new_title_i"
            ] = masking_consensus_proposed.get_message_mask_to_show_in_web_site()
            data[
                "new_class_i"
            ] = "padding-button-table {} fa fa-2x fa-pencil padding-button-table tip".format(
                "warning_fa_icon"
                if masking_consensus_proposed.has_masking_data()
                else ""
            )
            data["default"] = (
                "Has positions masked"
                if masking_consensus_proposed.has_masking_data()
                else "Has no values, yet"
            )
        ## define project if not set
        if project is None and not project_sample is None:
            project = project_sample.project

        if not project is None and b_change_data:
            ## check if they have projects
            count = ProjectSample.objects.filter(
                project=project, is_deleted=False, is_error=False, is_finished=True
            ).count()
            if count > 0:  ### need to send a message to recalculate the global files
                metaKeyAndValue = MetaKeyAndValue()
                manageDatabase = ManageDatabase()
                try:
                    process_SGE = ProcessSGE()
                    taskID = process_SGE.set_collect_global_files(project, request.user)
                    manageDatabase.set_project_metakey(
                        project,
                        request.user,
                        metaKeyAndValue.get_meta_key(
                            MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id
                        ),
                        MetaKeyAndValue.META_VALUE_Queue,
                        taskID,
                    )
                    data["is_ok"] = True
                except:
                    data = {"is_ok": False}
            else:  ### there's no projects to apply
                data["is_ok"] = True
                data[
                    "message"
                ] = "Mask/unmask consensus for the project '{}' are set.".format(
                    project.name
                )

        return JsonResponse(data)


def change_mask_consensus_in_project(project, masking_consensus_proposed):
    """change mask in project"""
    manageDatabase = ManageDatabase()
    genetic_element = DecodeObjects()

    meta_value = manageDatabase.get_project_metakey_last(
        project,
        MetaKeyAndValue.META_KEY_Masking_consensus,
        MetaKeyAndValue.META_VALUE_Success,
    )
    masking_consensus_original = None
    if not meta_value is None:
        masking_consensus_original = genetic_element.decode_result(
            meta_value.description
        )

    if (
        masking_consensus_original is None
        or masking_consensus_original != masking_consensus_proposed
    ):
        manageDatabase.set_project_metakey(
            project,
            project.owner,
            MetaKeyAndValue.META_KEY_Masking_consensus,
            MetaKeyAndValue.META_VALUE_Success,
            masking_consensus_proposed.to_json(),
        )
        ## need to mask all the project_sample if exist
        for project_sample in project.project_samples.all():
            if project_sample.is_deleted:
                continue
            meta_value = manageDatabase.get_project_sample_metakey_last(
                project_sample,
                MetaKeyAndValue.META_KEY_Masking_consensus,
                MetaKeyAndValue.META_VALUE_Success,
            )
            if not meta_value is None:
                manageDatabase.set_project_sample_metakey(
                    project_sample,
                    project.owner,
                    MetaKeyAndValue.META_KEY_Masking_consensus,
                    MetaKeyAndValue.META_VALUE_Success,
                    masking_consensus_proposed.to_json(),
                )
        return True
    return False


def change_mask_consensus_in_project_sample(project_sample, masking_consensus_proposed):
    """change mask in project sample"""

    manageDatabase = ManageDatabase()
    genetic_element = DecodeObjects()
    meta_value = manageDatabase.get_project_sample_metakey_last(
        project_sample,
        MetaKeyAndValue.META_KEY_Masking_consensus,
        MetaKeyAndValue.META_VALUE_Success,
    )
    masking_consensus_original = None
    if not meta_value is None:
        masking_consensus_original = genetic_element.decode_result(
            meta_value.description
        )

    if (
        masking_consensus_original is None
        or masking_consensus_original != masking_consensus_proposed
    ):
        manageDatabase.set_project_sample_metakey(
            project_sample,
            project_sample.project.owner,
            MetaKeyAndValue.META_KEY_Masking_consensus,
            MetaKeyAndValue.META_VALUE_Success,
            masking_consensus_proposed.to_json(),
        )
        return True
    return False


def get_fresh_masking_consensus(reference):
    """get fresh masking consensus"""
    manageDatabase = ManageDatabase()
    genetic_element = DecodeObjects()
    meta_value = manageDatabase.get_reference_metakey_last(
        reference,
        MetaKeyAndValue.META_KEY_Elements_And_CDS_Reference,
        MetaKeyAndValue.META_VALUE_Success,
    )
    masking_consensus_original = genetic_element.decode_result(meta_value.description)
    for element in masking_consensus_original.get_sorted_elements():
        masking_consensus_original.dt_elements_mask[element] = MaskingConsensus()
    return masking_consensus_original


@csrf_protect
def get_mask_consensus_actual_values(request):
    """
    return mask consensus of actual values
    """
    data = {"is_ok": False}
    if request.is_ajax():

        ## some pre-requisites
        if not request.user.is_active or not request.user.is_authenticated:
            return JsonResponse(data)
        try:
            profile = Profile.objects.get(user__pk=request.user.pk)
        except Profile.DoesNotExist:
            return JsonResponse(data)
        if profile.only_view_project:
            return JsonResponse(data)

        project_id_a = "project_id"
        project_sample_id_a = "project_sample_id"

        project_id = None
        project_sample_id = None
        if project_id_a in request.GET:
            project_id = request.GET[project_id_a]
        elif project_sample_id_a in request.GET:
            project_sample_id = request.GET[project_sample_id_a]

        manageDatabase = ManageDatabase()
        genetic_element = DecodeObjects()
        if not project_id is None:
            try:
                project = Project.objects.get(pk=project_id)
            except Project.DoesNotExist:
                return JsonResponse(data)

            meta_value = manageDatabase.get_project_metakey_last(
                project,
                MetaKeyAndValue.META_KEY_Masking_consensus,
                MetaKeyAndValue.META_VALUE_Success,
            )
            if meta_value is None:
                masking_consensus_original = get_fresh_masking_consensus(
                    project.reference
                )
            else:
                masking_consensus_original = genetic_element.decode_result(
                    meta_value.description
                )

            ### passing data
            data["all_data"] = masking_consensus_original.to_json()
            if (
                ProjectSample.objects.filter(
                    project=project, is_deleted=False, is_finished=True
                ).count()
                > 0
            ):
                data[
                    "warning_project"
                ] = "If you modify mask values it will change all the masks in the samples associated to the this project."
            data["is_ok"] = True
        elif not project_sample_id is None:
            try:
                project_sample = ProjectSample.objects.get(pk=project_sample_id)
            except ProjectSample.DoesNotExist:
                return JsonResponse(data)

            meta_value = manageDatabase.get_project_sample_metakey_last(
                project_sample,
                MetaKeyAndValue.META_KEY_Masking_consensus,
                MetaKeyAndValue.META_VALUE_Success,
            )
            if meta_value is None:
                ### test project
                meta_value = manageDatabase.get_project_metakey_last(
                    project_sample.project,
                    MetaKeyAndValue.META_KEY_Masking_consensus,
                    MetaKeyAndValue.META_VALUE_Success,
                )
                if meta_value is None:
                    masking_consensus_original = get_fresh_masking_consensus(
                        project_sample.project.reference
                    )
                else:
                    masking_consensus_original = genetic_element.decode_result(
                        meta_value.description
                    )
            else:
                masking_consensus_original = genetic_element.decode_result(
                    meta_value.description
                )

            ### passing data
            data["all_data"] = masking_consensus_original.to_json()
            data["is_ok"] = True
        return JsonResponse(data)


@csrf_protect
def turn_on_off_software(request):
    """
    Denies is_to_run in main software description.
    Don't do this if the software already run
    """
    if request.is_ajax():
        data = {"is_ok": False}
        software_id_a = "software_id"
        sample_id_a = "sample_id"
        project_id_a = "project_id"
        project_sample_id_a = "project_sample_id"
        type_of_use_id_a = "type_of_use_id"
        televir_project_id_a = "televir_project_id"

        ## some pre-requisites
        if not request.user.is_active or not request.user.is_authenticated:
            return JsonResponse(data)
        try:
            profile = Profile.objects.get(user__pk=request.user.pk)
        except Profile.DoesNotExist:
            return JsonResponse(data)
        if profile.only_view_project:
            return JsonResponse(data)

        pipeline_makeup = Pipeline_Makeup()

        sample_id = None
        project_id = None
        televir_project_id = None
        project_sample_id = None
        type_of_use_id = None

        if type_of_use_id_a in request.GET:
            type_of_use_id = int(request.GET[type_of_use_id_a])
        if televir_project_id_a in request.GET:
            televir_project_id = int(request.GET[televir_project_id_a])
        elif sample_id_a in request.GET:
            sample_id = request.GET[sample_id_a]
        elif project_id_a in request.GET:
            project_id = request.GET[project_id_a]
        elif project_sample_id_a in request.GET:
            project_sample_id = request.GET[project_sample_id_a]

        default_parameters = DefaultParameters()
        if software_id_a in request.GET:
            software_id = request.GET[software_id_a]

            try:
                project, televir_project, project_sample, sample = (
                    None,
                    None,
                    None,
                    None,
                )
                software = Software.objects.get(pk=software_id)
                current_is_to_run = software.is_to_run
                if not televir_project_id is None:

                    televir_project = PIProjects.objects.get(pk=televir_project_id)
                    pipeline_steps_project = (
                        pipeline_makeup.get_pipeline_makeup_result_of_operation(
                            software,
                            turn_off=current_is_to_run,
                            televir_project=televir_project,
                        )
                    )

                    makeup = pipeline_makeup.match_makeup_name_from_list(
                        pipeline_steps_project
                    )
                    if makeup is None:
                        if current_is_to_run:
                            data[
                                "message"
                            ] = f"You cannot perform this operation. Project '{televir_project}' would not meet minimum pipeline step requirements."

                            return JsonResponse(data)

                if not type_of_use_id is None:
                    if type_of_use_id == Software.TYPE_OF_USE_televir_global:

                        pipeline_steps_project = (
                            pipeline_makeup.get_pipeline_makeup_result_of_operation(
                                software, turn_off=current_is_to_run
                            )
                        )

                        makeup = pipeline_makeup.match_makeup_name_from_list(
                            pipeline_steps_project
                        )

                        if makeup is None:
                            if current_is_to_run:

                                data[
                                    "message"
                                ] = "You cannot perform this operation. Deployment would not meet minimum pipeline step requirements."

                                return JsonResponse(data)

                if not project_id is None:  ##    project
                    project = Project.objects.get(pk=project_id)
                elif not project_sample_id is None:  ##    project sample
                    project_sample = ProjectSample.objects.get(pk=project_sample_id)

                    ## sample is not ready for run again, To prevent to many ON|OFF
                    if not project_sample.is_finished:
                        data[
                            "message"
                        ] = "You cannot do this operation. The project '{}' with sample '{}' is in pipeline to run.".format(
                            project_sample.project.name, project_sample.sample.name
                        )
                        return JsonResponse(data)

                    ### create a task to perform the analysis of snippy and freebayes
                    manageDatabase = ManageDatabase()
                    process_SGE = ProcessSGE()
                    metaKeyAndValue = MetaKeyAndValue()

                    try:
                        ### the unique can be ON|OFF
                        if software.name == SoftwareNames.SOFTWARE_FREEBAYES_name:
                            ### change flag to not finished
                            project_sample.is_finished = False
                            project_sample.save()

                            (
                                job_name_wait,
                                job_name,
                            ) = request.user.profile.get_name_sge_seq(
                                Profile.SGE_PROCESS_dont_care, Profile.SGE_GLOBAL
                            )

                            if project_sample.is_sample_illumina():
                                taskID = process_SGE.set_second_stage_snippy(
                                    project_sample,
                                    request.user,
                                    job_name,
                                    [job_name_wait],
                                )
                            else:
                                taskID = process_SGE.set_second_stage_medaka(
                                    project_sample,
                                    request.user,
                                    job_name,
                                    [job_name_wait],
                                )

                            ### set project sample queue ID
                            manageDatabase.set_project_sample_metakey(
                                project_sample,
                                request.user,
                                metaKeyAndValue.get_meta_key_queue_by_project_sample_id(
                                    project_sample.id
                                ),
                                MetaKeyAndValue.META_VALUE_Queue,
                                taskID,
                            )

                        ## Generate Consensus
                        ### need to collect global files again
                        taskID = process_SGE.set_collect_global_files(
                            project_sample.project, request.user
                        )
                        manageDatabase.set_project_metakey(
                            project,
                            request.user,
                            metaKeyAndValue.get_meta_key(
                                MetaKeyAndValue.META_KEY_Queue_TaskID_Project,
                                project.id,
                            ),
                            MetaKeyAndValue.META_VALUE_Queue,
                            taskID,
                        )
                    except:
                        pass

                elif not sample_id is None:  ## for Sample
                    sample = Sample.objects.get(pk=sample_id)

                    ## sample is not ready for run again
                    if sample.is_sample_in_the_queue:
                        data[
                            "message"
                        ] = "You cannot do this operation. The sample '{}' is in pipeline to run.".format(
                            sample.name
                        )
                        return JsonResponse(data)

                    ### re-run data
                    manageDatabase = ManageDatabase()
                    process_SGE = ProcessSGE()

                    ### create a task to perform the analysis again
                    try:
                        (
                            job_name_wait,
                            job_name,
                        ) = request.user.profile.get_name_sge_seq(
                            Profile.SGE_PROCESS_clean_sample, Profile.SGE_SAMPLE
                        )
                        if sample.is_type_fastq_gz_sequencing():
                            taskID = process_SGE.set_run_trimmomatic_species(
                                sample, request.user, job_name
                            )
                        else:
                            taskID = process_SGE.set_run_clean_minion(
                                sample, request.user, job_name
                            )

                        ### set sample queue ID
                        manageDatabase.set_sample_metakey(
                            sample,
                            sample.owner,
                            MetaKeyAndValue.META_KEY_Queue_TaskID,
                            MetaKeyAndValue.META_VALUE_Queue,
                            taskID,
                        )
                    except:
                        sample.is_sample_in_the_queue = False
                        sample.save()
                        data["message"] = "Error in the queue system."
                        return JsonResponse(data)

                    ## refresh sample list for this user
                    if not job_name is None:
                        process_SGE.set_create_sample_list_by_user(
                            request.user, [job_name]
                        )

                ## set ON|OFF software
                is_to_run = default_parameters.set_software_to_run_by_software(
                    software, project, televir_project, project_sample, sample
                )

                ## set a new default
                data["is_to_run"] = is_to_run
                data["message"] = "The '{}' in '{}' technology was turned '{}'.".format(
                    software.name_extended,
                    software.technology.name,
                    "ON" if is_to_run else "OFF",
                )

            except Software.DoesNotExist:
                return JsonResponse(data)
            except Project.DoesNotExist:
                return JsonResponse(data)
            except ProjectSample.DoesNotExist:
                return JsonResponse(data)
            except Sample.DoesNotExist:
                return JsonResponse(data)
            data["is_ok"] = True
        return JsonResponse(data)


@csrf_protect
def get_software_name_to_turn_on_off(request):

    data = {"is_ok": False, "message": "You are not allow to do this operation."}

    if request.is_ajax():

        software_id_a = "software_id"
        sample_id_a = "sample_id"
        project_id_a = "project_id"
        project_sample_id_a = "project_sample_id"
        type_of_use_id_a = "type_of_use_id"
        televir_project_id_a = "televir_project_id"

        ## some pre-requisites
        if not request.user.is_active or not request.user.is_authenticated:
            return JsonResponse(data)
        try:
            profile = Profile.objects.get(user__pk=request.user.pk)
        except Profile.DoesNotExist:
            return JsonResponse(data)
        if profile.only_view_project:
            return JsonResponse(data)

        sample_id = None
        project_id = None
        project_sample_id = None
        type_of_use_id = None
        televir_project_id = None

        if type_of_use_id_a in request.GET:
            type_of_use_id = request.GET[type_of_use_id_a]
        if televir_project_id_a in request.GET:
            televir_project_id = request.GET[televir_project_id_a]
        if sample_id_a in request.GET:
            sample_id = request.GET[sample_id_a]
        elif project_id_a in request.GET:
            project_id = request.GET[project_id_a]
        elif project_sample_id_a in request.GET:
            project_sample_id = request.GET[project_sample_id_a]

        if software_id_a in request.GET:
            software_id = request.GET[software_id_a]
            try:
                project, project_sample, sample = None, None, None
                software = Software.objects.get(pk=software_id)
                if not project_id is None:
                    project = Project.objects.get(pk=project_id)
                elif not project_sample_id is None:
                    project_sample = ProjectSample.objects.get(pk=project_sample_id)
                elif not sample_id is None:
                    sample = Sample.objects.get(pk=sample_id)

                if software.can_be_on_off_in_pipeline:
                    ## get parameters for a specific sample, project or project_sample
                    parameters = Parameter.objects.filter(
                        software=software,
                        project=project,
                        project_sample=project_sample,
                        sample=sample,
                    )

                    if software.type_of_use == Software.TYPE_OF_USE_global:
                        is_to_run = software.is_to_run
                    elif len(parameters) > 0:
                        is_to_run = parameters[0].is_to_run
                    else:
                        is_to_run = software.is_to_run

                    data["is_ok"] = True
                    data[
                        "message"
                    ] = "Do you want to turn {} '{}' in '{}' technology?.".format(
                        "OFF" if is_to_run else "ON",
                        software.name_extended,
                        software.technology.name,
                    )
                else:
                    data[
                        "message"
                    ] = "'{}' in '{}' technology can not be ON/OFF.".format(
                        software.name_extended, software.technology.name
                    )

            except Software.DoesNotExist:
                return JsonResponse(data)
        return JsonResponse(data)


@csrf_protect
def reset_project_settings(request):
    data = {"is_ok": False, "message": "You are not allow to do this operation."}
    if request.is_ajax():

        televir_project_id_a = "project_id"
        televir_id = int(request.GET[televir_project_id_a])

        try:

            project_parameters = Parameter.objects.filter(
                televir_project__pk=televir_id
            )
            project_software = Software.objects.filter(
                parameter__televir_project__pk=televir_id
            )

        except Exception as e:
            print(e)
            data["message"] = "Error in the database."
            return JsonResponse(data)

        for parameter in project_parameters:
            parameter.delete()

        for software in project_software:
            software.delete()

        data["is_ok"] = True
        data["message"] = "Project settings were reset."

    return JsonResponse(data)
