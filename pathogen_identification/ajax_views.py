from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.views.decorators.http import require_POST
from utils.process_SGE import ProcessSGE

from pathogen_identification.models import Projects


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
        data = {"is_ok": False}

        process_SGE = ProcessSGE()

        project_id = int(request.POST["project_id"])
        project = Projects.objects.get(id=int(project_id))

        taskID = process_SGE.set_submit_televir_job(
            user=request.user,
            project_pk=project.pk,
        )
        data = {"is_ok": True}

        return JsonResponse(data)
