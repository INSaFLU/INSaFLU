from settings.models import Software, Parameter
from constants.software_names import SoftwareNames
from pathogen_identification.models import Projects
from dataclasses import dataclass

@dataclass
class RemapParams:
    max_taxids: int
    max_accids: int

def get_remap_software(username:str, project_name):
    """
    Get remap software
    """
    try: 
        project= Projects.objects.get(name=project_name)
        try: 

            remap = Software.objects.filter(
                name=SoftwareNames.SOFTWARE_REMAP_PARAMS_name, 
                owner__username= username, 
                technology__name=project.technology,
                type_of_use= Software.TYPE_OF_USE_televir_project)
        
        except Software.DoesNotExist:
            remap = Software.objects.filter(
                name=SoftwareNames.SOFTWARE_REMAP_PARAMS_name, 
                owner__username= username, 
                technology__name=project.technology,
                type_of_use= Software.TYPE_OF_USE_televir_global)
        
    except Software.DoesNotExist:
        raise Exception(f"Remap software not found for user {username}")    

    remap_params= Parameter.objects.filter(software=remap, televir_project__name=project_name)
    if remap_params.count() == 0:
        remap_params= Parameter.objects.filter(software=remap, televir_project__name=None)
    if remap_params.count() == 0:
        raise Exception(f"Remap software parameters not found for user {username} and project {project_name}")
    max_taxids= 0
    max_accids= 0
    for param in remap_params:

        if param.name == SoftwareNames.SOFTWARE_REMAP_PARAMS_max_taxids:
            max_taxids= int(param.parameter)
        elif param.name == SoftwareNames.SOFTWARE_REMAP_PARAMS_max_accids:
            max_accids= int(param.parameter)
    
    remap= RemapParams(max_taxids=max_taxids, max_accids=max_accids)

    
    return remap