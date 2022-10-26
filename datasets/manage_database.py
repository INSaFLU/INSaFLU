'''
Created on 12/06/2022

@author: mmp
'''

from datasets.models import MetaKey, MetaKeyDataset
from constants.meta_key_and_values import MetaKeyAndValue
from utils.lock_atomic_transaction import LockedAtomicTransaction

class ManageDatabase(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        pass
    
    def _get_metakey(self, meta_key_name):
        """
        get metakey with lobk table
        """
        with LockedAtomicTransaction(MetaKey):
            try:
                metaKey = MetaKey.objects.get(name=meta_key_name)
            except MetaKey.DoesNotExist:
                metaKey = MetaKey()
                metaKey.name = meta_key_name
                metaKey.save()
        return metaKey
    
    def set_dataset_metakey(self, dataset, owner, meta_key_name, value, description):
        """
        save a meta key
        """
        metaKey = self._get_metakey(meta_key_name)
        
        metaKeyDataset = MetaKeyDataset()
        metaKeyDataset.dataset = dataset
        metaKeyDataset.meta_tag = metaKey
        metaKeyDataset.owner = owner
        metaKeyDataset.value = value
        metaKeyDataset.description = description
        metaKeyDataset.save()
        return MetaKeyDataset
    
    def update_dataset_metakey(self, dataset, owner, meta_key_name, value, description):
        """
        update the meta_key, if not exist create other
        """
        metaKey = self._get_metakey(meta_key_name)
        
        metaKeyDataset_list = MetaKeyDataset.objects.filter(dataset=dataset, meta_tag=metaKey, owner=owner)
        if (metaKeyDataset_list.count() == 0):
            return self.set_dataset_metakey(dataset, owner, meta_key_name, value, description)
        else:
            meta_key_dataset = metaKeyDataset_list[0]
            meta_key_dataset.value = value
            meta_key_dataset.description = description
            meta_key_dataset.save()
        return meta_key_dataset

    
    def get_dataset_metakey(self, dataset, meta_key_name, value):
        """
        value = None, return a list
        """
        try:
            if (value == None): return MetaKeyDataset.objects.filter(dataset__id=dataset.id, meta_tag__name=meta_key_name).order_by('-creation_date')
            return MetaKeyDataset.objects.get(dataset__id=dataset.id, meta_tag__name=meta_key_name, value=value)
        except MetaKeyDataset.DoesNotExist:
            return None
        
    def get_dataset_metakey_last(self, dataset, meta_key_name, value):
        """
        value = None, return a list
        """
        try:
            if (value == None): query_set = MetaKeyDataset.objects.filter(dataset__id=dataset.id, meta_tag__name=meta_key_name).order_by('-creation_date')
            else: query_set = MetaKeyDataset.objects.filter(dataset__id=dataset.id, meta_tag__name=meta_key_name, value=value).order_by('-creation_date')
            if (query_set.count() > 0 ): return query_set[0]
            return None
        except MetaKeyDataset.DoesNotExist:
            return None    


    def is_dataset_processing_step(self, dataset, tag_key):
        """
            META_VALUE_Error = "Error"
            META_VALUE_Success = "Success"
            META_VALUE_Queue = "Queue"
            the end of a process is 
        """ 
    
        meta_dataset_queue = self.get_dataset_metakey_last(dataset, tag_key, MetaKeyAndValue.META_VALUE_Queue)
        if (meta_dataset_queue is None): return True
        
        ### try to find queue_success 
        meta_dataset_queue_success = self.get_dataset_metakey_last(dataset, tag_key, MetaKeyAndValue.META_VALUE_Success)        
        if (meta_dataset_queue_success is None): return True
        if (meta_dataset_queue_success.creation_date > meta_dataset_queue.creation_date and \
            meta_dataset_queue_success.description == meta_dataset_queue.description): return False
        return True
    
    #######################################
    ###
    ###        Other methods
    ###
    
    def get_max_length_label(self, dataset, user, b_calculate_again):
        """
        Get the max length of the samples in a specific dataset in chars
        b_calculate_again = True calculate again
        """
        
        if (not b_calculate_again):
            meta_data = self.get_dataset_metakey_last(dataset, MetaKeyAndValue.META_KEY_Dataset_max_name_length,\
                                        MetaKeyAndValue.META_VALUE_Success)
            if (meta_data != None):
                return int(meta_data.description)
            
        ### calculate and save
        n_max_value = 0
        for dataset_consensus in dataset.dataset_consensus.all():
            if (dataset_consensus.is_deleted): continue
            if (dataset_consensus.is_error): continue
            if (len(dataset_consensus.get_name()) > n_max_value): n_max_value = len(dataset_consensus.get_name())
            
        meta_data = self.update_dataset_metakey(dataset, user, MetaKeyAndValue.META_KEY_Dataset_max_name_length,\
                            MetaKeyAndValue.META_VALUE_Success, str(n_max_value))
        return n_max_value
    ###
    ###
    ###
    #######################################
