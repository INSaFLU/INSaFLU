import json

import pandas as pd
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.http import HttpResponse
from django.shortcuts import get_object_or_404, render
from django.views.decorators.http import require_POST

from pathogen_identification.constants_settings import ConstantsSettings
from pathogen_identification.models import PIProject_Sample, Projects, SoftwareTreeNode
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
        data = {"is_ok": False}
        project_id = int(request.POST["project_id"])
        project = Projects.objects.get(id=int(project_id))
        technology = project.technology
        samples = PIProject_Sample.objects.filter(project=project)

        utils = Utils_Manager(owner=request.user)
        combined_table = utils.parameter_util.generate_combined_parameters_table(
            utils.owner
        )
        all_paths = utils.get_all_technology_pipelines(technology)
        pipeline_tree = utils.generate_software_tree(technology)
        print(all_paths)
        print("user pk: ", request.user.pk)
        print("project pk: ", project.pk)
        print("sample pk: ", samples[0].pk)

        leaf_node_key = SoftwareTreeNode.LEAF_node
        leaf_node = SoftwareTreeNode.objects.filter(
            software_tree__technology=technology, node_place=leaf_node_key
        )
        print("leaf node pk: ", leaf_node[0].pk)
