import os

from constants.meta_key_and_values import MetaKeyAndValue
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse, JsonResponse
from django.utils.safestring import mark_safe
from django.views.decorators.csrf import csrf_protect
from django.views.decorators.http import require_POST
from fluwebvirus.settings import STATIC_ROOT, STATIC_URL
from utils.process_SGE import ProcessSGE

from pathogen_identification.models import (
    PIProject_Sample,
    Projects,
    ReferenceMap_Main,
    RunMain,
)
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager


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
        runs_to_deploy = utils.check_runs_to_deploy(user, project)

        try:
            if runs_to_deploy:

                taskID = process_SGE.set_submit_televir_job(
                    user=request.user,
                    project_pk=project.pk,
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


def validate_project_name(request):
    if request.is_ajax():
        data = {"is_taken": False}

        if request.method == "GET":
            user_obj = Projects.objects.filter(
                owner=request.user,
                name=request.GET.get("projectname"),
                is_deleted=False,
            ).exists()

            has_spaces = " " in request.GET.get("projectname")

            if user_obj:
                return HttpResponse("exists")

            if has_spaces:
                return HttpResponse("has_spaces")

            return HttpResponse(False)


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
