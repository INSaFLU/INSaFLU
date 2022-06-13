'''
Created on 13/06/2022

@author: mmp
'''

from utils.collect_extra_data import CollectExtraData

class MetaRow(object):
    
    def __init__(self, seq_name_consensus, row):
        self.seq_name_consensus = seq_name_consensus
        self.row = row

class Metadata(object):
    
    def __init__(self, header):
        '''
        Constructor
        '''
        self.header = header
        self.dt_header = dict(zip(header, list(range(0, len(header)))))
        self.vect_rows_id = []
        self.dt_rows_id = {}
        
    def add_metadata(self, project_sample_pk, seq_name_consensus, row):
        
        if not project_sample_pk in self.dt_rows_id:
            self.dt_rows_id[project_sample_pk] = MetaRow(seq_name_consensus, row)
    
    def get_vect_out(self, vect_header_out, csv_writer):
    
        dt_out_id_project_sample = {}
        count = 0
        for project_sample_pk in self.dt_rows_id:
            if project_sample_pk in dt_out_id_project_sample: continue
            dt_out_id_project_sample[project_sample_pk] = 1
            vect_out = [self.dt_rows_id[project_sample_pk].seq_name_consensus]
                
            for column in vect_header_out:
                if column == CollectExtraData.HEADER_SAMPLE_OUT_ID: continue 
                if column in self.dt_header: vect_out.append(self.dt_rows_id[project_sample_pk].row[self.dt_header[column]])
                else: vect_out.append("")
            count += 1 
            csv_writer.writerow(vect_out)
        return count


class DataColumns(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.dt_project = {}
    
    def add_header(self, project_pk, header):
        """
        add header for this project
        """
    
        if not project_pk in self.dt_project: 
            self.dt_project[project_pk] = Metadata(header)
        
    def add_metadata(self, project_pk, project_sample_pk, seq_name_consensus, row):
        """
        add metadata for a specific project,project_sample
        """
        if project_pk in self.dt_project:
            self.dt_project[project_pk].add_metadata(project_sample_pk, seq_name_consensus, row)
        
        
    def _get_header(self):
        """
        return header
        ## CollectExtraData.HEADER_SAMPLE_OUT_CSV_simple
        """
        
        self.vect_header_out = []
        dt_header_out = {}
        
        if len(self.dt_project) > 0:
            for key_metadata in self.dt_project:
                if len(self.vect_header_out) == 0: 
                    self.vect_header_out = self.dt_project[key_metadata].header
                    ## set CollectExtraData.HEADER_SAMPLE_OUT_ID be the first
                    try:
                        self.vect_header_out.pop(self.vect_header_out.index(CollectExtraData.HEADER_SAMPLE_OUT_ID))
                    except ValueError as e:
                        pass
                    self.vect_header_out = [CollectExtraData.HEADER_SAMPLE_OUT_ID] + self.vect_header_out
                    dt_header_out = dict(zip(self.vect_header_out, [1] * len(self.vect_header_out)))
                else:
                    for column in self.dt_project[key_metadata].header:
                        if not column in dt_header_out:
                            self.vect_header_out.append(column)
                            dt_header_out[column] = 1
        return self.vect_header_out

    
    def save_rows(self, csv_writer):
        """
        save all data
        """
        count = 0 ## number of rows saved
        
        ## header save
        csv_writer.writerow(self._get_header())
        
        ## save data
        for key_metadata in self.dt_project:
            count += self.dt_project[key_metadata].get_vect_out(
                self.vect_header_out, csv_writer) 
        return count
            