import mimetypes
import os

from django.contrib.auth.decorators import login_required
from django.http import HttpResponse, JsonResponse
from django.utils.safestring import mark_safe
from django.utils.translation import ugettext_lazy as _
from django.views.decorators.csrf import csrf_protect
from django.views.decorators.http import require_POST

from constants.meta_key_and_values import MetaKeyAndValue
from fluwebvirus.settings import STATIC_ROOT, STATIC_URL
from managing_files.models import ProcessControler
from pathogen_identification.models import (FinalReport, ParameterSet,
                                            PIProject_Sample, Projects,
                                            ReferenceMap_Main, RunMain)
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager
from utils.process_SGE import ProcessSGE


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
def submit_televir_project_sample(request):
    """
    submit a new sample to televir project
    """
    print("heloo")
    print(request.is_ajax())
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user

        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(id=int(sample_id))
        project = Projects.objects.get(id=int(sample.project.pk))

        utils = Utils_Manager()

        runs_to_deploy = utils.check_runs_to_deploy_sample(user, project, sample)
        print("RUNS TO DEPLOY")

        try:
            if len(runs_to_deploy) > 0:
                for sample, leafs_to_deploy in runs_to_deploy.items():
                    for leaf in leafs_to_deploy:
                        print(request.user.username, project.pk, sample.pk, leaf.pk)

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
                process_SGE.kill_televir_process_controler(
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
def deploy_ProjectPI(request):
    """
    prepare data for deployment of pathogen identification.
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user

        project_id = int(request.POST["project_id"])
        project = Projects.objects.get(id=int(project_id))

        utils = Utils_Manager()
        runs_to_deploy = utils.check_runs_to_deploy_project(user, project)
        print(runs_to_deploy)

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
def deploy_televir_map(request):
    """
    prepare data for deployment of pathogen identification.
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user

        reference_id = int(request.POST["reference_id"])
        taskID = process_SGE.set_submit_televir_map(user, reference_pk=reference_id)

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
