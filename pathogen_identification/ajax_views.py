import json

import pandas as pd
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.shortcuts import get_object_or_404, render
from django.views.decorators.http import require_POST

from pathogen_identification.constants_settings import ConstantsSettings
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

        utils.generate_software_parameter_dict(parameters_table)

        utils.create_pipe_tree()
        utils.generate_graph()
        all_paths = utils.get_all_graph_paths()  #
        print(all_paths)

        soft_unique = parameters_table.software_name.unique()
        for sof in soft_unique:

            exists_ok = utils.check_software_is_installed(sof)
            print(sof, exists_ok)

        # for path in all_paths:
        #    path_explicit = [utils.node_index[x][1] for x in path]
        #    path_df = utils.df_from_path(path_explicit)
        #    print(path_df)

        # print(parameters_table.software_id)
        # print("parameters_table", parameters_table.pipeline_step_id)
