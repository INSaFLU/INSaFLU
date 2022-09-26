import argparse
import itertools as it
import json
import os
import shutil

import numpy as np
import pandas as pd
from settings.models import Parameter, Software
from this import d


class Utils:
    def __init__():
        pass

    def get_parameters_available(project_id):
        """
        Get the parameters available for a project
        """
        project = Projects.objects.get(id=project_id)
        owner = project.owner
        software_available = Software.objects.filter(owner=owner)

        try:
            parameters_available = project.parameters.all()
        except:
            parameters_available = Parameter.objects.filter(
                software__in=software_available
            )

        return parameters_available
