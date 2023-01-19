from django.contrib.auth.models import User
from django.core.management.base import BaseCommand
from managing_files.models import ProcessControler
from pathogen_identification.constants_settings import ConstantsSettings
from pathogen_identification.deployment_main import Run_Main_from_Leaf
from pathogen_identification.models import (
    ParameterSet,
    PIProject_Sample,
    Projects,
    SoftwareTree,
    SoftwareTreeNode,
)
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager
from utils.process_SGE import ProcessSGE

import pandas as pd


class Run_Node:
    run: str
    step: tuple
    module: str
    params: pd.DataFrame

    def __init__(self, run: str, step: tuple, module: str):
        self.run = run
        self.step = step
        self.module = module

    def __repr__(self):
        return f"RunCapsule({self.run}, {self.step}"

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        return self.run == other.run and self.step == other.step


class Capsule_Manager:
    module_order: list = [
        "PREPROCESS",
        "ENRICHMENT",
        "ASSEMBLY",
        "READ_CLASSIFICATION",
        "CONTIG_CLASSIFICATION",
        "REMAPPING",
    ]
    tree: pipeline_tree
    capsule_dict: dict

    def __init__(self, tree: pipeline_tree):
        self.tree = tree
        self.capsule_dict = {}
