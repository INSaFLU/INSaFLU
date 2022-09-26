import json

import pandas as pd
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.shortcuts import get_object_or_404, render
from django.views.decorators.http import require_POST

from pathogen_identification.models import PIProject_Sample, Projects
from pathogen_identification.utils import Utils


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
        utils = Utils()
        project_id = int(request.POST["project_id"])
        parameters_table = utils.get_parameters_available(project_id)

        print("parameters_table", parameters_table.columns)
