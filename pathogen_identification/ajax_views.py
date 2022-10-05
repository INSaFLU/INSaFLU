import json

import pandas as pd
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.http import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404, render
from django.views.decorators.http import require_POST
from utils.process_SGE import ProcessSGE

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
        process_SGE = ProcessSGE()

        data = {"is_ok": False}
        project_id = int(request.POST["project_id"])
        project = Projects.objects.get(id=int(project_id))
        technology = project.technology
        print(technology)
        samples = PIProject_Sample.objects.filter(project=project)

        utils = Utils_Manager(owner=request.user)

        pipeline_tree = utils.generate_software_tree(technology)
        pipeline_tree_index = utils.get_software_tree_index(technology)

        print("user pk: ", request.user.pk)
        print("project pk: ", project.pk)
        print("sample pk: ", samples[0].pk)
        print("pipeline_tree_index: ", pipeline_tree_index)

        local_tree = utils.generate_current_tree(technology)

        local_paths = local_tree.get_all_graph_paths_explicit()
        print("local paths: ", local_paths.keys())

        try:
            for leaf, path in local_paths.items():
                print("leaf: ", leaf)

                matched_path = utils.utility_manager.match_path_to_tree(
                    path, pipeline_tree
                )

                matched_path_pk = SoftwareTreeNode.objects.get(
                    software_tree__pk=pipeline_tree_index, index=matched_path
                ).pk

                print("matched path: ", leaf, matched_path)
                if matched_path:
                    taskID = process_SGE.set_submit_televir_job(
                        user=request.user,
                        sample_pk=samples[1].pk,
                        project_pk=project.pk,
                        pipeline_pk=matched_path_pk,
                    )
            data = {"is_ok": True}
        except:
            data = {"is_ok": False}

        return JsonResponse(data)
