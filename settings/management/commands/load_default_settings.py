'''
Created on 29/11/2021

@author: mmp
'''
from django.core.management import BaseCommand
from settings.default_parameters import DefaultParameters
from settings.constants_settings import ConstantsSettings
from settings.models import PipelineStep, Software, Technology
import logging

class Command(BaseCommand):
    '''
    classdocs
    '''
    help = "Test settings parameters."
    
    ## logging
    logger_debug = logging.getLogger("fluWebVirus.debug")
    logger_production = logging.getLogger("fluWebVirus.production")
    
    def __init__(self, *args, **kwargs):
        super(Command, self).__init__(*args, **kwargs)
    
    def load_default_settings(self):
        """
        Upload default files
        """
        ### Pipeline Name
        constants_settings = ConstantsSettings()
        for name in constants_settings.vect_pipeline_names:
            try:
                pipeline_step = PipelineStep.objects.get(name=name)
            except PipelineStep.DoesNotExist as e:
                pipeline_step = PipelineStep()
                pipeline_step.name = name
                pipeline_step.save()
        
        ### Technology Name
        for name in constants_settings.vect_technology:
            try:
                technology = Technology.objects.get(name=name)
            except Technology.DoesNotExist:
                technology = Technology()
                technology.name = name
                technology.save()
    
    def set_pipelines_in_previous_softwares(self):
        """
        set pipelines,this is done because of compatibilities between versions
        """
        default_parameters = DefaultParameters()
        for software in Software.objects.all():
            ## if obsolete continue
            if software.is_obsolete: continue
            
            ### need to set the pipeline...
            if software.pipeline_step is None:
                vect_parameters = default_parameters.get_vect_parameters(software)
                software.can_be_on_off_in_pipeline = vect_parameters[0].software.can_be_on_off_in_pipeline
                software.is_to_run = vect_parameters[0].software.is_to_run
                software.pipeline_step = vect_parameters[0].software.pipeline_step
                software.help_text = vect_parameters[0].software.help_text
                software.save()
            else:   ### if PipelineStep not none, test if it is correct
                vect_parameters = default_parameters.get_vect_parameters(software)
                if (software.pipeline_step.name != vect_parameters[0].software.pipeline_step.name):
                    software.pipeline_step = vect_parameters[0].software.pipeline_step
                    software.save()
                    
           

    # A command must define handle()
    def handle(self, *args, **options):
        
        #### set default fields
        self.stdout.write("Set default pipelines..")
        self.load_default_settings()
        
        #### set obsolete some softwares because of parameters
        self.stdout.write("Set obsolete softwares because of parameters...")
        default_parameters = DefaultParameters()
        default_parameters.set_software_obsolete()
        
        #### set obsolete some softwares because of parameters
        self.stdout.write("Set pipelines in previous softwares...")
        self.set_pipelines_in_previous_softwares()
        
        self.stdout.write("End")


