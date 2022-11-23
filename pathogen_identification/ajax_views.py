from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_protect
from django.views.decorators.http import require_POST
from utils.process_SGE import ProcessSGE

from pathogen_identification.deployment_main import Run_Main_from_Leaf
from pathogen_identification.models import (
    FinalReport,
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
        print("checking")

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
    print(request.is_ajax())
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user

        print(request.POST)

        reference_id = int(request.POST["reference_id"])
        taskID = process_SGE.set_submit_televir_map(user, reference_pk=reference_id)

        data["is_ok"] = True

        return JsonResponse(data)


@csrf_protect
def IGV_display(request):
    """display python plotly app"""
    if request.is_ajax():
        data = {"is_ok": False}
        if request.method == "GET":
            project_pk = request.GET.get("project_pk")
            sample_pk = request.GET.get("sample_pk")
            run_pk = request.GET.get("run_pk")
            reference = request.GET.get("reference")
            unique_id = request.GET.get("unique_id")

            sample = PIProject_Sample.objects.get(pk=int(sample_pk))
            sample_name = sample.name
            run = RunMain.objects.get(pk=int(run_pk))
            run_name = run.name

            try:
                ref_map = ReferenceMap_Main.objects.get(
                    reference=reference, sample=sample, run=run
                )
                final_report = FinalReport.objects.get(
                    sample=sample, run=run, unique_id=unique_id
                )

                def remove_pre_static(path, pattern):
                    # path = path.split(pattern)[1]
                    # path = f"/{pattern}{path}"
                    # path = \path.replace(pattern, "")
                    print(path)
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

                data["is_ok"] = True

                data["path_reference"] = path_name_reference
                data["path_reference_index"] = path_name_reference_index
                data["path_bam"] = path_name_bam
                data["path_bai"] = path_name_bai

                data["reference_name"] = sample_name
                data["sample_name"] = final_report.reference_contig_str

                #### other files
                data["bam_file_id"] = mark_safe(
                    '<strong>Bam file:</strong> <a href="{}" filename="{}">{}</a>'.format(
                        "download_file",
                        os.path.basename(path_name_bam),
                        os.path.basename(path_name_bam),
                    )
                )
                data["bai_file_id"] = mark_safe(
                    '<strong>Bai file:</strong> <a href="{}" filename="{}">{}</a>'.format(
                        path_name_bai,
                        os.path.basename(path_name_bai),
                        os.path.basename(path_name_bai),
                    )
                )
                data["reference_id"] = mark_safe(
                    '<strong>Reference:</strong> <a href="{}" filename="{}">{}</a>'.format(
                        path_name_reference,
                        os.path.basename(path_name_reference),
                        os.path.basename(path_name_reference),
                    )
                )
                data["reference_index_id"] = mark_safe(
                    '<strong>Ref. index:</strong> <a href="{}" filename="{}">{}</a>'.format(
                        path_name_reference_index,
                        os.path.basename(path_name_reference_index),
                        os.path.basename(path_name_reference_index),
                    )
                )

                data["static_dir"] = run.static_dir
                print(run.static_dir)
                data["sample_name"] = sample_name

            except ReferenceMap_Main.DoesNotExist as e:
                pass
        return JsonResponse(data)
