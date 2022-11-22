from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.views.decorators.http import require_POST
from utils.process_SGE import ProcessSGE

from pathogen_identification.deployment_main import Run_Main_from_Leaf
from pathogen_identification.models import Projects
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
