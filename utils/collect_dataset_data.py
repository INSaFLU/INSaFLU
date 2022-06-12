'''
Created on 12/06/2022

@author: mmp
'''
import os, csv, time, json, logging
from constants.constants import Constants, TypePath, FileExtensions
from utils.utils import Utils
from django.conf import settings
from utils.software import Software
from utils.process_SGE import ProcessSGE
from managing_files.models import ProcessControler
from datasets.models import Dataset
from constants.meta_key_and_values import MetaKeyAndValue
from datasets.manage_database import ManageDatabase
from utils.tree import CreateTree

class CollectExtraDatasetData(object):
    '''
    classdocs
    '''

    ## Type of sample list
    DATASET_LIST_simple = 0          ## used in the tree
    DATASET_LIST_list = 1            ## list of samples/references/consensus with metadata
    
    utils = Utils()
    software = Software()
    if settings.DEBUG: logger = logging.getLogger("fluWebVirus.debug")
    else: logger = logging.getLogger("fluWebVirus.production")

    def __init__(self):
        '''
        Constructor
        '''
        pass
    
    def collect_extra_data_for_dataset(self, dataset, user):
        """
        """
        ### make it running 
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        process_SGE.set_process_controler(user, process_controler.get_name_dataset(dataset), ProcessControler.FLAG_RUNNING)
        
        ## need to add a delay for the test in command line
        if (settings.RUN_TEST_IN_COMMAND_LINE): time.sleep(4)
        
        ## run collect data
        self.__collect_extra_data_for_dataset(dataset, user)
    
    def collect_update_extra_metadata_for_dataset(self, dataset, user):
        """
        Only for update metadata
        """
        ### make it running 
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        process_SGE.set_process_controler(user, process_controler.get_name_dataset(dataset), ProcessControler.FLAG_RUNNING)
        
        ## need to add a delay for the test in command line
        if (settings.RUN_TEST_IN_COMMAND_LINE): time.sleep(4)
        
        ## run collect data
        self.__collect_update_extra_metadata_for_dataset(dataset, user)
        
    def __collect_update_extra_metadata_for_dataset(self, dataset, user):
        """
        Only for update metadata
        It is not necessary 
        """
        ### get the taskID and seal it
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        manage_database = ManageDatabase()
        metaKeyAndValue = MetaKeyAndValue()
        
        try:
            ## collect sample table with plus type and subtype, mixed infection, equal to upload table
            self.calculate_global_files(Dataset.DATASET_FILE_NAME_SAMPLE_RESULT_CSV, dataset, user)
            self.calculate_global_files(Dataset.DATASET_FILE_NAME_SAMPLE_RESULT_TSV, dataset, user)
            ## IMPORTANT -> this need to be after of Project.DATASET_FILE_NAME_SAMPLE_RESULT_CSV
            #self.calculate_global_files(Project.DATASET_FILE_NAME_SAMPLE_RESULT_json, dataset, user)
            
            ### zip several files to download 
            self.zip_several_files(dataset)
        except:
            ## finished with error
            process_SGE.set_process_controler(user, process_controler.get_name_dataset(dataset), ProcessControler.FLAG_ERROR)
            return
        
        ## seal the tag        
        meta_dataset = manage_database.get_dataset_metakey_last(dataset, metaKeyAndValue.get_meta_key(\
                        MetaKeyAndValue.META_KEY_Queue_TaskID_Project, dataset.id), MetaKeyAndValue.META_VALUE_Queue)
        if (meta_dataset != None):
            manage_database.set_dataset_metakey(dataset, user, metaKeyAndValue.get_meta_key(\
                    MetaKeyAndValue.META_KEY_Queue_TaskID_Project, dataset.id),
                    MetaKeyAndValue.META_VALUE_Success, meta_dataset.description)
                
        ### finished
        process_SGE.set_process_controler(user, process_controler.get_name_dataset(dataset), ProcessControler.FLAG_FINISHED)
        
    def __collect_extra_data_for_dataset(self, dataset, user):
        """
        Everything that is necessary to do in the dataset
        Collect all extra data after all samples are finished
        only run after all the vect_taskID are finished
        """
        ### get the taskID and seal it
        metaKeyAndValue = MetaKeyAndValue()
        manage_database = ManageDatabase()
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        
        try:
            count = 0
            start = time.time()
            self.logger.info("COLLECT_EXTRA_FILES: Start")
            
            ## calculate the max sample label size of the samples that belong to this dataset
            ## used in MSA viewer 
            b_calculate_again = True
            manage_database.get_max_length_label(dataset, user, b_calculate_again)
            self.logger.info("COLLECT_EXTRA_FILES: Step {}  diff_time:{}".format(count, time.time() - start))
            count += 1
            start = time.time()
            
            ## collect all consensus files for a dataset
            self.calculate_global_files(Dataset.DATASET_FILE_NAME_RESULT_all_consensus, dataset, user)
            self.logger.info("COLLECT_EXTRA_FILES: Step {}  diff_time:{}".format(count, time.time() - start))
            count += 1
            start = time.time()
            
            ## collect sample table with plus type and subtype, mixed infection, equal to upload table
            self.calculate_global_files(Dataset.DATASET_FILE_NAME_SAMPLE_RESULT_CSV, dataset, user)
            self.calculate_global_files(Dataset.DATASET_FILE_NAME_SAMPLE_RESULT_TSV, dataset, user)
            self.logger.info("COLLECT_EXTRA_FILES: Step {}  diff_time:{}".format(count, time.time() - start))
            count += 1
            
            ### create trees
            createTree = CreateTree()
            createTree.create_tree_and_alignments_dataset(dataset, user)
            self.logger.info("COLLECT_EXTRA_FILES: Step {}  diff_time:{}".format(count, time.time() - start))
            count += 1
            start = time.time()
            
            meta_dataset = manage_database.get_dataset_metakey_last(dataset, metaKeyAndValue.get_meta_key(\
                        MetaKeyAndValue.META_KEY_Queue_TaskID_Project, dataset.id), MetaKeyAndValue.META_VALUE_Queue)
            if (meta_dataset != None):
                manage_database.set_dataset_metakey(dataset, user, metaKeyAndValue.get_meta_key(\
                        MetaKeyAndValue.META_KEY_Queue_TaskID_Project, dataset.id),
                        MetaKeyAndValue.META_VALUE_Success, meta_dataset.description)
                
            ### zip several files to download 
            self.zip_several_files(dataset)
            self.logger.info("COLLECT_EXTRA_FILES: Step {}  diff_time:{}".format(count, time.time() - start))
            count += 1
            start = time.time()
        except:
            ## finished with error
            process_SGE.set_process_controler(user, process_controler.get_name_dataset(dataset), ProcessControler.FLAG_ERROR)
            return
        
        ### finished
        process_SGE.set_process_controler(user, process_controler.get_name_dataset(dataset), ProcessControler.FLAG_FINISHED)

    def calculate_global_files(self, type_file, dataset, user):
        """
        Collect extra files
         """
        out_file = None
        out_file_file_system = None
        if (type_file == Dataset.DATASET_FILE_NAME_SAMPLE_RESULT_CSV):
            ## samples csv
            out_file = self.collect_sample_table(dataset, Constants.SEPARATOR_COMMA, CollectExtraDatasetData.SAMPLE_LIST_list)
            out_file_file_system = dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, type_file)
        elif (type_file == Dataset.DATASET_FILE_NAME_SAMPLE_RESULT_TSV):
            ## samples tsv
            out_file = self.collect_sample_table(dataset, Constants.SEPARATOR_TAB, CollectExtraDatasetData.SAMPLE_LIST_list)
            out_file_file_system = dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, type_file)
        elif (type_file == Dataset.DATASET_FILE_NAME_RESULT_all_consensus):
            out_file = self.merge_all_consensus_files(dataset)
            out_file_file_system = dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, type_file)

        ## copy file
        if (not out_file is None):
            self.utils.copy_file(out_file, out_file_file_system)
            self.utils.remove_file(out_file)
        elif (not out_file_file_system is None and os.path.exists(out_file_file_system)): self.utils.remove_file(out_file_file_system)

    def merge_all_consensus_files(self, dataset):
        """
        merge all consensus files
        """
        out_file = self.utils.get_temp_file('all_consensus', FileExtensions.FILE_FASTA)
        vect_to_process = []
        for dataset_consensus in dataset.dataset_consensus.all():
            if (dataset_consensus.is_deleted): continue
            if (dataset_consensus.is_error): continue
            
            if (not dataset_consensus.is_ready_to_proccess()): continue
            if not os.path.exists(dataset_consensus.get_consensus_file(TypePath.MEDIA_ROOT)): continue

            vect_to_process.append([
                dataset_consensus.get_consensus_file(TypePath.MEDIA_ROOT),\
                dataset_consensus.get_name(), -1])
        
        ### set the number sequences that passed     
        dataset.number_passed_sequences = len(vect_to_process)
        dataset.save()
        
        self.utils.merge_fasta_files(vect_to_process, out_file)
        return out_file
