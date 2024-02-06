import mimetypes
import os
from datetime import datetime

import pandas as pd
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.http import HttpResponse, JsonResponse
from django.utils.safestring import mark_safe
from django.utils.translation import ugettext_lazy as _
from django.views.decorators.csrf import csrf_protect
from django.views.decorators.http import require_POST

from constants.constants import Constants
from constants.meta_key_and_values import MetaKeyAndValue
from fluwebvirus.settings import STATIC_ROOT, STATIC_URL
from managing_files.models import ProcessControler
from pathogen_identification.models import (
    FinalReport,
    MetaReference,
    ParameterSet,
    PIProject_Sample,
    Projects,
    ReferenceMap_Main,
    RunMain,
    TeleFluProject,
    TeleFluSample,
)
from pathogen_identification.tables import ReferenceSourceTable
from pathogen_identification.utilities.reference_utils import (
    check_metaReference_exists_from_ids,
    create_combined_reference,
)
from pathogen_identification.utilities.televir_parameters import TelevirParameters
from pathogen_identification.utilities.utilities_general import get_services_dir
from pathogen_identification.utilities.utilities_pipeline import SoftwareTreeUtils
from pathogen_identification.utilities.utilities_views import (
    ReportSorter,
    SampleReferenceManager,
    set_control_reports,
)
from utils.process_SGE import ProcessSGE
from utils.utils import Utils


def simplify_name(name):
    return (
        name.replace("_", "_")
        .replace("-", "_")
        .replace(" ", "_")
        .replace(".", "_")
        .lower()
    )


@login_required
@require_POST
def submit_sample_metagenomics_televir(request):
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()

        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(id=int(sample_id))

        user = sample.project.owner
        project = sample.project

        software_utils = SoftwareTreeUtils(user, project, sample=sample)
        runs_to_deploy = software_utils.check_runs_to_submit_metagenomics_sample(sample)
        reference_manager = SampleReferenceManager(sample)

        try:
            if len(runs_to_deploy) > 0:
                for sample, leaves_to_deploy in runs_to_deploy.items():
                    for leaf in leaves_to_deploy:
                        metagenomics_run = reference_manager.mapping_run_from_leaf(leaf)
                        taskID = process_SGE.set_submit_televir_sample_metagenomics(
                            user=request.user,
                            sample_pk=sample.pk,
                            leaf_pk=leaf.pk,
                            combined_analysis=True,
                            map_run_pk=metagenomics_run.pk,
                        )

                data["is_deployed"] = True

        except Exception as e:
            print(e)
            data["is_deployed"] = False

        data["is_ok"] = True
        print(data)
        return JsonResponse(data)


@login_required
@require_POST
def submit_sample_screening_televir(request):
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()

        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(id=int(sample_id))

        user = sample.project.owner
        project = sample.project
        reference_manager = SampleReferenceManager(sample)

        software_utils = SoftwareTreeUtils(user, project, sample=sample)
        runs_to_deploy = software_utils.check_runs_to_submit_screening_sample(sample)

        try:
            if len(runs_to_deploy) > 0:
                for sample, leaves_to_deploy in runs_to_deploy.items():
                    for leaf in leaves_to_deploy:
                        screening_run = reference_manager.screening_run_from_leaf(leaf)
                        taskID = process_SGE.set_submit_televir_sample_metagenomics(
                            user=request.user,
                            sample_pk=sample.pk,
                            leaf_pk=leaf.pk,
                            map_run_pk=screening_run.pk,
                        )

                data["is_deployed"] = True

        except Exception as e:
            print(e)
            data["is_deployed"] = False

        data["is_ok"] = True
        print(data)
        return JsonResponse(data)


from pathogen_identification.models import RawReference


@login_required
@require_POST
def submit_sample_mapping_televir(request):
    if request.is_ajax():
        data = {"is_ok": True, "is_deployed": False, "is_empty": False}

        process_SGE = ProcessSGE()

        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(id=int(sample_id))
        user = sample.project.owner
        project = sample.project
        reference_manager = SampleReferenceManager(sample)

        ### check if all references are already mapped
        reference_id_list = request.POST.getlist("reference_ids[]")
        added_references = RawReference.objects.filter(
            run__sample__pk=sample_id, run__run_type=RunMain.RUN_TYPE_STORAGE
        )

        if len(reference_id_list) == 0 and len(added_references) == 0:
            data["is_empty"] = True
            return JsonResponse(data)

        already_mapped = True

        #### check among reference id list

        for reference_id in reference_id_list:
            reference = RawReference.objects.get(pk=int(reference_id))
            if (
                RawReference.objects.filter(
                    accid=reference.accid,
                    status__in=[
                        RawReference.STATUS_MAPPED,
                        RawReference.STATUS_MAPPING,
                    ],
                    run__sample__pk=sample_id,
                ).exists()
                is False
            ):
                already_mapped = False

        ##### check among added references

        for added_reference in added_references:
            if (
                RawReference.objects.filter(
                    accid=added_reference.accid,
                    status__in=[
                        RawReference.STATUS_MAPPED,
                        RawReference.STATUS_MAPPING,
                    ],
                    run__sample__pk=sample_id,
                ).exists()
                is False
            ):
                already_mapped = False

        if already_mapped is True:
            data["is_ok"] = True
            data["is_deployed"] = False
            data["is_empty"] = False
            data["is_already_mapped"] = True
            return JsonResponse(data)

        ### runs to deploy
        software_utils = SoftwareTreeUtils(user, project, sample=sample)
        runs_to_deploy = software_utils.check_runs_to_submit_mapping_only(sample)

        if len(runs_to_deploy) == 0:
            return JsonResponse(data)

        try:
            if len(runs_to_deploy) > 0:
                for sample, leaves_to_deploy in runs_to_deploy.items():
                    for leaf in leaves_to_deploy:
                        mapping_run = reference_manager.mapping_request_run_from_leaf(
                            leaf
                        )

                        for reference_id in reference_id_list:
                            reference = RawReference.objects.get(pk=int(reference_id))
                            reference.pk = None
                            reference.run = mapping_run
                            reference.save()

                        for added_reference in added_references:
                            added_reference.pk = None
                            added_reference.run = mapping_run
                            added_reference.save()

                        taskID = process_SGE.set_submit_televir_sample_metagenomics(
                            user=request.user,
                            sample_pk=sample.pk,
                            leaf_pk=leaf.pk,
                            mapping_request=True,
                            map_run_pk=mapping_run.pk,
                        )

                data["is_deployed"] = True

        except Exception as e:
            print(e)
            data["is_ok"] = False

        print(data)
        return JsonResponse(data)


@login_required
@require_POST
def submit_project_samples_mapping_televir(request):
    if request.is_ajax():
        data = {
            "is_ok": True,
            "is_deployed": False,
            "is_empty": False,
            "samples_deployed": 0,
        }

        process_SGE = ProcessSGE()

        project_id = int(request.POST["project_id"])
        project = Projects.objects.get(id=int(project_id))
        user = project.owner

        project_samples = PIProject_Sample.objects.filter(project=project)

        try:
            samples_map_launched = []
            for sample in project_samples:
                sample_id = sample.pk
                reference_manager = SampleReferenceManager(sample)

                ### check if all references are already mapped
                added_references = RawReference.objects.filter(
                    run__sample__pk=sample_id, run__run_type=RunMain.RUN_TYPE_STORAGE
                )
                already_mapped = True

                for added_reference in added_references:
                    if (
                        RawReference.objects.filter(
                            accid=added_reference.accid,
                            status__in=[
                                RawReference.STATUS_MAPPED,
                                RawReference.STATUS_MAPPING,
                            ],
                            run__sample__pk=sample_id,
                        ).exists()
                        is False
                    ):
                        already_mapped = False

                if already_mapped is True:
                    continue
                else:
                    samples_map_launched.append(sample)

                ### runs to deploy
                software_utils = SoftwareTreeUtils(user, project, sample=sample)
                runs_to_deploy = software_utils.check_runs_to_submit_mapping_only(
                    sample
                )

                if len(runs_to_deploy) == 0:
                    continue

                for sample, leaves_to_deploy in runs_to_deploy.items():
                    for leaf in leaves_to_deploy:
                        mapping_run = reference_manager.mapping_request_run_from_leaf(
                            leaf
                        )

                        for added_reference in added_references:
                            added_reference.pk = None
                            added_reference.run = mapping_run
                            added_reference.save()

                        taskID = process_SGE.set_submit_televir_sample_metagenomics(
                            user=request.user,
                            sample_pk=sample.pk,
                            leaf_pk=leaf.pk,
                            mapping_request=True,
                            map_run_pk=mapping_run.pk,
                        )

            if len(samples_map_launched) > 0:
                data["is_deployed"] = True
                data["samples_deployed"] = len(samples_map_launched)

        except Exception as e:
            print(e)
            data["is_ok"] = False

        print(data)
        return JsonResponse(data)


@login_required
@require_POST
def deploy_ProjectPI(request):
    """
    prepare data for deployment of pathogen identification.
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()

        project_id = int(request.POST["project_id"])
        project = Projects.objects.get(id=int(project_id))

        user_id = int(request.POST["user_id"])
        user = User.objects.get(id=int(user_id))

        samples = PIProject_Sample.objects.filter(
            project=project, is_deleted_in_file_system=False
        )

        software_utils = SoftwareTreeUtils(user, project)

        try:
            for sample in samples:
                runs_to_deploy = software_utils.check_runs_to_deploy_sample(sample)

                if len(runs_to_deploy) > 0:
                    for sample, leafs_to_deploy in runs_to_deploy.items():
                        taskID = process_SGE.set_submit_televir_sample(
                            user=user,
                            project_pk=project.pk,
                            sample_pk=sample.pk,
                        )

                    data["is_deployed"] = True

        except Exception as e:
            print(e)
            data["is_deployed"] = False

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def deploy_ProjectPI_runs(request):
    """
    prepare data for deployment of pathogen identification.
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()

        project_id = int(request.POST["project_id"])
        project = Projects.objects.get(id=int(project_id))

        user_id = int(request.POST["user_id"])
        user = User.objects.get(id=int(user_id))

        software_utils = SoftwareTreeUtils(user, project)
        runs_to_deploy = software_utils.check_runs_to_deploy_project()

        try:
            if len(runs_to_deploy) > 0:
                for sample, leaves_to_deploy in runs_to_deploy.items():
                    for leaf in leaves_to_deploy:
                        taskID = process_SGE.set_submit_televir_run(
                            user=request.user,
                            project_pk=project.pk,
                            sample_pk=sample.pk,
                            leaf_pk=leaf.pk,
                        )

                data["is_deployed"] = True

        except Exception as e:
            print(e)
            data["is_deployed"] = False

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def deploy_ProjectPI_combined_runs(request):
    """
    prepare data for deployment of pathogen identification.
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()

        project_id = int(request.POST["project_id"])
        project = Projects.objects.get(id=int(project_id))

        user_id = int(request.POST["user_id"])
        user = User.objects.get(id=int(user_id))

        samples = PIProject_Sample.objects.filter(
            project=project, is_deleted_in_file_system=False
        )
        first_sample = samples.first()
        software_utils = SoftwareTreeUtils(user, project, sample=first_sample)
        runs_to_deploy = software_utils.check_runs_to_submit_metagenomics_sample(
            first_sample
        )

        if len(runs_to_deploy) == 0:
            data["is_ok"] = True
            return JsonResponse(data)

        try:
            if len(runs_to_deploy) > 0:
                for sample in samples:
                    software_utils = SoftwareTreeUtils(user, project, sample=sample)
                    runs_to_deploy = (
                        software_utils.check_runs_to_submit_metagenomics_sample(sample)
                    )
                    for sample, leaves_to_deploy in runs_to_deploy.items():
                        reference_manager = SampleReferenceManager(sample)
                        for leaf in leaves_to_deploy:
                            metagenomics_run = reference_manager.mapping_run_from_leaf(
                                leaf
                            )

                            taskID = process_SGE.set_submit_televir_sample_metagenomics(
                                user=request.user,
                                sample_pk=sample.pk,
                                leaf_pk=leaf.pk,
                                combined_analysis=True,
                                map_run_pk=metagenomics_run.pk,
                            )

                data["is_deployed"] = True

        except Exception as e:
            print(e)
            data["is_deployed"] = False

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def submit_televir_project_sample_runs(request):
    """
    submit a new sample to televir project
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user

        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(id=int(sample_id))
        project = Projects.objects.get(id=int(sample.project.pk))

        software_utils = SoftwareTreeUtils(user, project)
        runs_to_deploy = software_utils.check_runs_to_deploy_sample(sample)

        try:
            if len(runs_to_deploy) > 0:
                for sample, leafs_to_deploy in runs_to_deploy.items():
                    for leaf in leafs_to_deploy:
                        taskID = process_SGE.set_submit_televir_run(
                            user=request.user,
                            project_pk=project.pk,
                            sample_pk=sample.pk,
                            leaf_pk=leaf.pk,
                        )

                data["is_deployed"] = True

        except Exception as e:
            print(e)
            data["is_deployed"] = False

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def submit_televir_project_sample(request):
    """
    submit a new sample to televir project
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}
        process_SGE = ProcessSGE()
        user = request.user

        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(id=int(sample_id))
        project = Projects.objects.get(id=int(sample.project.pk))

        software_utils = SoftwareTreeUtils(user, project=project)
        runs_to_deploy = software_utils.check_runs_to_deploy_sample(sample)

        try:
            if len(runs_to_deploy) > 0:
                for sample, leafs_to_deploy in runs_to_deploy.items():
                    taskID = process_SGE.set_submit_televir_sample(
                        user=request.user,
                        project_pk=project.pk,
                        sample_pk=sample.pk,
                    )

                data["is_deployed"] = True

        except Exception as e:
            print(e)
            data["is_deployed"] = False

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def Project_explify_merge(request):
    """
    submit a new sample to televir project
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user
        utils: Utils = Utils()
        try:
            temp_directory = utils.get_temp_dir()

            project_id = int(request.POST["project_id"])
            project = Projects.objects.get(id=int(project_id))

            rpip_file = request.FILES["rpip_file"]
            upip_file = request.FILES["upip_file"]

            rpip_report_path = os.path.join(
                temp_directory,
                rpip_file.name.replace(" ", "_").replace("(", "_").replace(")", "_"),
            )
            upip_report_path = os.path.join(
                temp_directory,
                upip_file.name.replace(" ", "_").replace("(", "_").replace(")", "_"),
            )

            with open(rpip_report_path, "wb") as f:
                f.write(rpip_file.file.read())
            with open(upip_report_path, "wb") as f:
                f.write(upip_file.file.read())

        except Exception as e:
            print(e)
            return JsonResponse(data)

        ### check process not running for this project
        process_controler = ProcessControler()
        try:
            ProcessControler.objects.get(
                owner__id=user.pk,
                name=process_controler.get_name_televir_project_merge_explify(
                    project_pk=project.pk,
                ),
                is_running=True,
            )

        except ProcessControler.DoesNotExist:
            all_reports = FinalReport.objects.filter(
                run__project__pk=int(project.pk)
            ).order_by("-coverage")

            televir_reports = pd.DataFrame.from_records(all_reports.values())

            report_path = os.path.join(
                temp_directory, f"televir_project.{project.pk}.tsv"
            )
            televir_reports.to_csv(report_path, sep="\t", index=False)

            taskID = process_SGE.set_submit_televir_explify_merge(
                user=request.user,
                project_pk=project.pk,
                rpip_filepath=rpip_report_path,
                upip_filepath=upip_report_path,
                televir_report_filepath=report_path,
                out_dir=temp_directory,
            )

            data["is_deployed"] = True

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def Project_explify_merge_external(request):
    """
    merge explify rpip and upip reports to televir report, all provided by user.
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user
        utils: Utils = Utils()
        try:
            temp_directory = utils.get_temp_dir()

            rpip_file = request.FILES["rpip_file"]
            upip_file = request.FILES["upip_file"]
            project_file = request.FILES["project_file"]

            rpip_report_path = os.path.join(
                temp_directory,
                rpip_file.name.replace(" ", "_").replace("(", "_").replace(")", "_"),
            )
            upip_report_path = os.path.join(
                temp_directory,
                upip_file.name.replace(" ", "_").replace("(", "_").replace(")", "_"),
            )

            report_path = os.path.join(
                temp_directory,
                project_file.name.replace(" ", "_").replace("(", "_").replace(")", "_"),
            )

            with open(rpip_report_path, "wb") as f:
                f.write(rpip_file.file.read())
            with open(upip_report_path, "wb") as f:
                f.write(upip_file.file.read())
            with open(report_path, "wb") as f:
                f.write(project_file.file.read())

        except Exception as e:
            print(e)
            return JsonResponse(data)

        ### check process not running for this project
        process_controler = ProcessControler()
        try:
            ProcessControler.objects.get(
                owner__id=user.pk,
                name=process_controler.get_name_televir_project_merge_explify_external(
                    user_pk=user.pk,
                ),
                is_running=True,
            )

        except ProcessControler.DoesNotExist:
            taskID = process_SGE.set_submit_televir_explify_merge_external(
                user=request.user,
                rpip_filepath=rpip_report_path,
                upip_filepath=upip_report_path,
                televir_report_filepath=report_path,
                out_dir=temp_directory,
            )

            data["is_deployed"] = True

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def Project_explify_delete_external(request):
    """
    delete external televir report
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user
        utils: Utils = Utils()
        try:
            output_file_merged = os.path.join(
                get_services_dir(user), "merged_televir_explify.tsv"
            )
            os.remove(output_file_merged)
        except Exception as e:
            print(e)
            return JsonResponse(data)

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def kill_televir_project_sample(request):
    """
    kill all processes a sample, set queued to false
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user

        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(id=int(sample_id))
        project = Projects.objects.get(id=int(sample.project.pk))

        runs = ParameterSet.objects.filter(
            sample=sample,
            status__in=[
                ParameterSet.STATUS_RUNNING,
                ParameterSet.STATUS_QUEUED,
            ],
        )

        for run in runs:
            try:  # kill process
                process_SGE.kill_televir_process_controler_runs(
                    user.pk, project.pk, sample.pk, run.leaf.pk
                )

            except ProcessControler.DoesNotExist as e:
                print(e)
                print("ProcessControler.DoesNotExist")
                pass

            if run.status == ParameterSet.STATUS_RUNNING:
                run.delete_run_data()

            run.status = ParameterSet.STATUS_KILLED
            run.save()

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def kill_televir_project_tree_sample(request):
    """
    kill all processes a sample, set queued to false
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user

        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(id=int(sample_id))
        project = Projects.objects.get(id=int(sample.project.pk))

        try:  # kill process
            process_SGE.kill_televir_process_controler_samples(
                user.pk,
                project.pk,
                sample.pk,
            )

        except ProcessControler.DoesNotExist as e:
            print(e)
            print("ProcessControler.DoesNotExist")
            pass

        runs = ParameterSet.objects.filter(
            sample=sample,
            status__in=[
                ParameterSet.STATUS_RUNNING,
                ParameterSet.STATUS_QUEUED,
            ],
        )

        for run in runs:
            if run.status == ParameterSet.STATUS_RUNNING:
                run.delete_run_data()

            run.status = ParameterSet.STATUS_KILLED
            run.save()

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def sort_report_projects(request):
    """
    sort report projects
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}
        process_SGE = ProcessSGE()
        samples = PIProject_Sample.objects.filter(
            project__pk=int(request.POST["project_id"])
        )

        project = Projects.objects.get(id=int(request.POST["project_id"]))
        samples = PIProject_Sample.objects.filter(project=project)
        report_layout_params = TelevirParameters.get_report_layout_params(
            project_pk=project.pk
        )
        try:
            for sample in samples:
                final_reports = FinalReport.objects.filter(sample=sample)

                report_sorter = ReportSorter(final_reports, report_layout_params)

                if report_sorter.run is None:
                    pass
                elif report_sorter.check_analyzed():
                    pass
                else:
                    taskID = process_SGE.set_submit_televir_sort_pisample_reports(
                        user=request.user,
                        pisample_pk=sample.pk,
                    )
                    data["is_deployed"] = True

        except Exception as e:
            print(e)
            return JsonResponse(data)

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def sort_report_sample(request):
    """
    sort report projects
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}
        process_SGE = ProcessSGE()
        sample = PIProject_Sample.objects.get(pk=int(request.POST["sample_id"]))
        references = request.POST.getlist("references[]")

        project = sample.project
        report_layout_params = TelevirParameters.get_report_layout_params(
            project_pk=project.pk
        )
        try:
            final_reports = FinalReport.objects.filter(sample=sample)
            report_sorter = ReportSorter(final_reports, report_layout_params)

            if report_sorter.run is None:
                pass
            elif report_sorter.check_analyzed():
                pass
            else:
                taskID = process_SGE.set_submit_televir_sort_pisample_reports(
                    user=request.user,
                    pisample_pk=sample.pk,
                )
                data["is_deployed"] = True

        except Exception as e:
            print(e)
            return JsonResponse(data)

        data["is_ok"] = True
        return JsonResponse(data)


from constants.constants import Constants, FileExtensions, FileType, TypePath
from constants.software_names import SoftwareNames
from managing_files.models import ProjectSample as InsafluProjectSample
from pathogen_identification.models import RawReference, ReferenceSourceFileMap
from pathogen_identification.utilities.reference_utils import (
    check_reference_exists,
    check_reference_submitted,
)
from pathogen_identification.utilities.televir_bioinf import TelevirBioinf


@login_required
@require_POST
def teleflu_igv_create(request):
    print("teleflu_igv_create")
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        teleflu_project_pk = int(request.POST["pk"])
        teleflu_project = TeleFluProject.objects.get(pk=teleflu_project_pk)
        insaflu_project = teleflu_project.insaflu_project

        ### get reference
        reference = teleflu_project.reference
        if reference is None:
            return JsonResponse(data)

        reference_file = reference.get_reference_fasta(TypePath.MEDIA_ROOT)

        samples = InsafluProjectSample.objects.filter(project=insaflu_project)
        # samples= [sample.sample for sample in samples]
        sample_dict = {}

        ### get sample files
        software_names = SoftwareNames()

        for sample in samples:
            bam_file = sample.get_file_output(
                TypePath.MEDIA_ROOT, FileType.FILE_BAM, software_names.get_snippy_name()
            )
            bam_file_index = sample.get_file_output(
                TypePath.MEDIA_ROOT,
                FileType.FILE_BAM_BAI,
                software_names.get_snippy_name(),
            )
            vcf_file = sample.get_file_output(
                TypePath.MEDIA_ROOT, FileType.FILE_VCF, software_names.get_snippy_name()
            )

            if bam_file and bam_file_index and vcf_file:
                sample_dict[sample.sample.pk] = {
                    "name": sample.sample.name,
                    "bam_file": bam_file,
                    "bam_file_index": bam_file_index,
                    "vcf_file": vcf_file,
                }

        ### merge vcf files
        televir_bioinf = TelevirBioinf()
        vcf_files = [files["vcf_file"] for sample_pk, files in sample_dict.items()]
        group_vcf = teleflu_project.project_vcf
        stacked_html = teleflu_project.project_igv_report_media

        os.makedirs(teleflu_project.project_vcf_directory, exist_ok=True)

        merged_success = televir_bioinf.merge_vcf_files(vcf_files, group_vcf)

        try:

            televir_bioinf.create_igv_report(
                reference_file,
                vcf_file=group_vcf,
                tracks=sample_dict,
                output_html=stacked_html,
            )

            # for sample_pk, files in sample_dict.items():
            #    print(sample_pk, files)
        except Exception as e:
            print(e)
            return JsonResponse(data)

        return JsonResponse(data)


@login_required
@require_POST
def create_insaflu_reference(request):
    if request.is_ajax():
        data = {"is_ok": False, "exists": False}

        ref_id = int(request.POST["ref_id"])
        user_id = int(request.POST["user_id"])
        user = User.objects.get(id=user_id)
        process_SGE = ProcessSGE()

        try:
            if check_reference_exists(ref_id, user_id) or check_reference_submitted(
                ref_id=ref_id, user_id=user_id
            ):
                data["is_ok"] = True
                data["exists"] = True
                return JsonResponse(data)
            # success = create_reference(ref_id, user_id)
            taskID = process_SGE.set_submit_televir_teleflu_create(user, ref_id)

        except Exception as e:
            print(e)
            return JsonResponse(data)

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def add_references_to_sample(request):
    """
    add references to sample
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_error": False, "is_empty": False}
        temp_directory = Utils().get_temp_dir()
        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(pk=sample_id)

        reference_id_list = request.POST.getlist("reference_ids[]")

        if len(reference_id_list) == 0:
            data["is_empty"] = True
            return JsonResponse(data)

        reference_id_list = [int(x) for x in reference_id_list]

        references_existing = []

        sample_reference_manager = SampleReferenceManager(sample)

        try:
            ref_sources = ReferenceSourceFileMap.objects.filter(
                pk__in=reference_id_list
            )

            for reference in ref_sources:
                if RawReference.objects.filter(
                    accid=reference.accid,
                    run__sample__pk=sample_id,
                ).exists():
                    references_existing.append(reference.accid)
                    continue

                # truncate reference description: max 150 characters
                ref_description = reference.description
                if len(ref_description) > 150:
                    ref_description = ref_description[:150]

                new_reference = RawReference(
                    run=sample_reference_manager.storage_run,
                    accid=reference.accid,
                    taxid=reference.taxid,
                    description=ref_description,
                    status=RawReference.STATUS_UNMAPPED,
                    counts=0,
                    classification_source="none",
                )

                new_reference.save()

        except Exception as e:
            print(e)
            data["is_error"] = True
            return JsonResponse(data)

        data = {"is_ok": True}
        return JsonResponse(data)


from django.template.loader import render_to_string

from fluwebvirus.settings import BASE_DIR


def inject_references(references: list, request):
    context = {}
    data = {}

    context["references_table"] = ReferenceSourceTable(references)
    context["references_count"] = len(references)

    template_table_html = os.path.join(
        BASE_DIR,
        "templates",
        "pathogen_identification/references_table_table_only.html",
    )
    template_table_html = "pathogen_identification/references_table_table_only.html"

    # render tamplate using context
    rendered_table = render_to_string(template_table_html, context, request=request)
    data["my_content"] = rendered_table
    data["references_count"] = len(references)

    return data


@login_required
@csrf_protect
def create_teleflu_project(request):
    """
    create teleflu project
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_error": False, "exists": False, "is_empty": False}

        print(request.POST)
        ref_ids = request.POST.getlist("ref_ids[]")
        sample_ids = request.POST.getlist("sample_ids[]")

        def teleflu_project_name_from_refs(ref_ids):
            refs = [RawReference.objects.get(pk=int(x)) for x in ref_ids]
            date_now_str = datetime.now().strftime("%Y%m%d")
            if len(refs) == 0:
                return "teleflu_project"

            if len(refs) == 1:
                return f"televir_project_{refs[0].accid}_{date_now_str}"

            return f"televir_project_multiple_refs_{date_now_str}"

        def teleflu_project_description(ref_ids):
            if len(ref_ids) == 0:
                return "teleflu_project"
            if len(ref_ids) == 1:
                return "single reference project"
            return "multiple references project"

        project_name = teleflu_project_name_from_refs(ref_ids)
        first_ref = RawReference.objects.get(pk=int(ref_ids[0]))
        project = first_ref.run.project
        date = datetime.now()
        process_SGE = ProcessSGE()

        try:
            print("reference exists")
            if check_metaReference_exists_from_ids(ref_ids):
                data["exists"] = True
                return JsonResponse(data)

            metareference = create_combined_reference(ref_ids, project_name)
            print(metareference)

            if not metareference:
                data["is_error"] = True
                return JsonResponse(data)

            teleflu_project = TeleFluProject(
                televir_project=project,
                name=project_name,
                last_change_date=date,
                description=teleflu_project_description(ref_ids),
                raw_reference=metareference,
            )
            teleflu_project.save()
            print(teleflu_project)

            for sample_id in sample_ids:
                sample = PIProject_Sample.objects.get(pk=int(sample_id))
                TeleFluSample.objects.create(
                    teleflu_project=teleflu_project,
                    televir_sample=sample,
                )

            process_SGE.set_submit_televir_teleflu_project_create(
                user=request.user,
                project_pk=teleflu_project.pk,
            )

            data["is_ok"] = True
            data["project_id"] = teleflu_project.pk
            data["project_name"] = teleflu_project.name

        except Exception as e:
            print(e)
            data["is_error"] = True
            return JsonResponse(data)

        return JsonResponse(data)


@csrf_protect
def set_teleflu_check_box_values(request):
    """
    manage check boxes through ajax
    """
    if request.is_ajax():
        data = {"is_ok": False}
        utils = Utils()
        print(request.GET)
        if Constants.GET_CHECK_BOX_SINGLE in request.GET:
            data["is_ok"] = True
            for key in request.session.keys():
                if (
                    key.startswith(Constants.TELEFLU_CHECK_BOX)
                    and len(key.split("_")) == 4
                    and utils.is_integer(key.split("_")[3])
                ):
                    data[key] = request.session[key]
        ## change single status of a check_box_single
        elif Constants.GET_CHANGE_CHECK_BOX_SINGLE in request.GET:
            data["is_ok"] = True
            key_name = "{}_{}".format(
                Constants.TELEFLU_CHECK_BOX, request.GET.get(Constants.CHECK_BOX_VALUE)
            )
            print(request.session.keys())
            for key in request.session.keys():
                if (
                    key.startswith(Constants.TELEFLU_CHECK_BOX)
                    and len(key.split("_")) == 4
                    and utils.is_integer(key.split("_")[3])
                ):
                    if request.session[key]:
                        data[key] = False
                    if key == key_name:
                        request.session[key] = utils.str2bool(
                            request.GET.get(Constants.GET_CHANGE_CHECK_BOX_SINGLE)
                        )
                    else:
                        request.session[key] = False

        return JsonResponse(data)


@login_required
@require_POST
def add_references_all_samples(request):
    """
    add references to sample
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_error": False, "is_empty": False}
        temp_directory = Utils().get_temp_dir()
        project_id = int(request.POST["ref_id"])
        project = Projects.objects.get(pk=project_id)
        samples = PIProject_Sample.objects.filter(project=project)

        reference_id_list = request.POST.getlist("reference_ids[]")
        data["empty_content"] = inject_references([], request)["my_content"]
        if len(reference_id_list) == 0:
            data["is_empty"] = True
            return JsonResponse(data)

        reference_id_list = [int(x) for x in reference_id_list]

        references_existing = []

        try:
            ref_sources = ReferenceSourceFileMap.objects.filter(
                pk__in=reference_id_list
            )
            for sample in samples:
                sample_reference_manager = SampleReferenceManager(sample)
                sample_id = sample.pk

                for reference in ref_sources:
                    if RawReference.objects.filter(
                        accid=reference.accid,
                        run__sample__pk=sample_id,
                    ).exists():
                        references_existing.append(reference.accid)
                        continue

                    # truncate reference description: max 150 characters
                    ref_description = reference.description

                    if len(ref_description) > 150:
                        ref_description = ref_description[:150]

                    new_reference = RawReference(
                        run=sample_reference_manager.storage_run,
                        accid=reference.accid,
                        taxid=reference.taxid,
                        description=ref_description,
                        status=RawReference.STATUS_UNMAPPED,
                        counts=0,
                        classification_source="none",
                    )

                    new_reference.save()

        except Exception as e:
            print(e)
            data["is_error"] = True
            return JsonResponse(data)

        data["is_ok"] = True

        return JsonResponse(data)


from pathogen_identification.views import inject__added_references


@login_required
@require_POST
def remove_added_reference(request):
    """
    remove added reference
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_error": False}

        reference_id = int(request.POST["reference_id"])
        sample_id = int(request.POST["sample_id"])

        try:
            reference = RawReference.objects.get(pk=reference_id)
            reference.delete()
        except Exception as e:
            print(e)
            data["is_error"] = True
            return JsonResponse(data)

        query_set_added_manual = RawReference.objects.filter(
            run__sample__pk=sample_id, run__run_type=RunMain.RUN_TYPE_STORAGE
        )

        context = inject__added_references(query_set_added_manual, request)
        data["added_references"] = context["my_content"]

        data = {"is_ok": True}
        return JsonResponse(data)


@login_required
@require_POST
def deploy_televir_map(request):
    """
    prepare data for deployment of pathogen identification.
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user

        reference_id = int(request.POST["reference_id"])
        project_id = int(request.POST["project_id"])
        taskID = process_SGE.set_submit_televir_map(
            user, reference_pk=reference_id, project_pk=project_id
        )

        data["is_ok"] = True

        return JsonResponse(data)


@login_required
@require_POST
def set_sample_reports_control(request):
    """
    set sample reports control
    """
    if request.is_ajax():
        data = {"is_ok": False}
        data["set_control"] = False
        sample_id = int(request.POST["sample_id"])
        try:
            sample = PIProject_Sample.objects.get(pk=int(sample_id))

            sample_control_flag = False if sample.is_control else True

            sample_reports = FinalReport.objects.filter(sample=sample)
            for report in sample_reports:
                report.control_flag = FinalReport.CONTROL_FLAG_NONE

                report.save()

            sample.is_control = sample_control_flag
            sample.save()

            set_control_reports(sample.project.pk)

            data["is_ok"] = True
            data["set_control"] = sample_control_flag
            return JsonResponse(data)

        except Exception as e:
            print(e)
            data["is_ok"] = False
            return JsonResponse(data)


@csrf_protect
def validate_project_name(request):
    """
    test if exist this project name
    """
    if request.is_ajax():
        project_name = request.GET.get("project_name")

        data = {
            "is_taken": Projects.objects.filter(
                name__iexact=project_name,
                is_deleted=False,
                owner__username=request.user.username,
            ).exists()
        }

        ## check if name has spaces:
        if " " in project_name:
            data["has_spaces"] = True
            data["error_message"] = _("Spaces are not allowed in the project name.")

        ## check if name has special characters:
        if not project_name.replace("_", "").isalnum():
            data["has_special_characters"] = True
            data["error_message"] = _(
                "Special characters are not allowed in the project name."
            )

        if data["is_taken"]:
            data["error_message"] = _("Exists a project with this name.")

        return JsonResponse(data)


@csrf_protect
def IGV_display(request):
    """display python plotly app"""

    if request.is_ajax():
        data = {"is_ok": False}
        if request.method == "GET":
            sample_pk = request.GET.get("sample_pk")
            run_pk = request.GET.get("run_pk")
            reference = request.GET.get("accid")
            unique_id = request.GET.get("unique_id")

            sample = PIProject_Sample.objects.get(pk=int(sample_pk))
            sample_name = sample.name
            run = RunMain.objects.get(pk=int(run_pk))

            ref_map = ReferenceMap_Main.objects.get(
                reference=unique_id, sample=sample, run=run
            )

            def remove_pre_static(path: str, pattern: str) -> str:
                cwd = os.getcwd()
                if path.startswith(cwd):
                    path = path[len(cwd) :]

                path = path.replace(STATIC_ROOT, STATIC_URL)

                return path

            path_name_bam = remove_pre_static(
                ref_map.bam_file_path, "/insaflu_web/INSaFLU/"
            )
            path_name_bai = remove_pre_static(
                ref_map.bai_file_path, "/insaflu_web/INSaFLU/"
            )
            path_name_reference = remove_pre_static(
                ref_map.fasta_file_path, "/insaflu_web/INSaFLU/"
            )
            path_name_reference_index = remove_pre_static(
                ref_map.fai_file_path, "/insaflu_web/INSaFLU/"
            )
            path_name_vcf = remove_pre_static(ref_map.vcf, "/insaflu_web/INSaFLU/")

            data["is_ok"] = True
            data["path_bam"] = mark_safe(request.build_absolute_uri(path_name_bam))

            data["path_reference"] = mark_safe(
                request.build_absolute_uri(path_name_reference)
            )
            data["path_reference_index"] = mark_safe(
                request.build_absolute_uri(path_name_reference_index)
            )
            data["reference_name"] = reference

            #### other files
            data["bam_file_id"] = mark_safe(
                '<strong>Bam file:</strong> <a href="{}" download="{}"> {}</a>'.format(
                    path_name_bam,
                    os.path.basename(path_name_bam),
                    os.path.basename(path_name_bam),
                )
            )
            data["bai_file_id"] = mark_safe(
                '<strong>Bai file:</strong> <a href="{}" download="{}"> {}</a>'.format(
                    path_name_bai,
                    os.path.basename(path_name_bai),
                    os.path.basename(path_name_bai),
                )
            )
            data["vcf_file_id"] = mark_safe(
                '<strong>Vcf file:</strong> <a href="{}" download="{}"> {}</a>'.format(
                    path_name_vcf,
                    os.path.basename(path_name_vcf),
                    os.path.basename(path_name_vcf),
                )
            )
            data["reference_id"] = mark_safe(
                '<strong>Reference:</strong> <a href="{}" download="{}"> {}</a>'.format(
                    path_name_reference,
                    os.path.basename(path_name_reference),
                    os.path.basename(path_name_reference),
                )
            )
            data["reference_index_id"] = mark_safe(
                '<strong>Ref. index:</strong> <a href="{}" download="{}"> {}</a>'.format(
                    path_name_reference_index,
                    os.path.basename(path_name_reference_index),
                    os.path.basename(path_name_reference_index),
                )
            )

            data["static_dir"] = run.static_dir
            data["sample_name"] = sample_name

        return JsonResponse(data)
