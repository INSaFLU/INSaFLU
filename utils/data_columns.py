'''
Created on 13/06/2022

@author: mmp
'''

import csv
from utils.collect_extra_data import CollectExtraData
from constants.constants import Constants

"""
#*Mandatory metadata fields:
#strain    Sample ID (Characters “()[]{}|#><” are disallowed)
#date    YEAR-MONTH-DAY ex: 2021-02-19
#virus    ncov
#region    Africa, Asia, Europe, North America, Oceania or South America
#gisaid_epi_isl    GISAID ID; if not available needs to be “?” 
#genbank_accession    Genbank accession #; if not available needs to be “?” 
#length    Genome length; can be filled with “?” 
#segment        Filled with “genome”
#sex    host sex; if not available needs to be “?” 
#age    host age; if not available needs to be “?” 
#host    host; if not available needs to be “?”  - from ncov apparently it is not mandatory??

# id in the Sample_list.tsv corresponds to strain in nextstrain metadata
# fastq1 and fastq2 in Sample_list is to ignore
# "data set" in Sample_list is to ignore

# Lineage Pangolin (for ncov) corresponds to lineage in nextstrain metadata
"""


NEXTSTRAIN_strain = "strain"    ##    Sample ID (Characters “()[]{}|#><” are disallowed)
NEXTSTRAIN_date = "date"        ##    YEAR-MONTH-DAY ex: 2021-02-19
NEXTSTRAIN_virus = "virus"      ## ncov
NEXTSTRAIN_region = "region"    ##   Africa, Asia, Europe, North America, Oceania or South America
NEXTSTRAIN_gisaid_epi_isl = "gisaid_epi_isl" #    GISAID ID; if not available needs to be “?” 
NEXTSTRAIN_genbank_accession = "genbank_accession"  #  Genbank accession #; if not available needs to be “?” 
NEXTSTRAIN_length = "length"    ##   Genome length; can be filled with “?” 
NEXTSTRAIN_segment = "segment"  ##   Filled with “genome”
NEXTSTRAIN_sex = "sex"          ## host sex; if not available needs to be “?” 
NEXTSTRAIN_age = "age"          ## host age; if not available needs to be “?” 
NEXTSTRAIN_host = "host"        ## host; if not available needs to be “?”  - from ncov apparently it is not mandatory??
VECT_NEXTSTRAIN_manatory = [
        NEXTSTRAIN_strain,
        NEXTSTRAIN_date,
        NEXTSTRAIN_virus,
        NEXTSTRAIN_region,
        NEXTSTRAIN_gisaid_epi_isl, 
        NEXTSTRAIN_genbank_accession, 
        NEXTSTRAIN_length, 
        NEXTSTRAIN_segment,
        NEXTSTRAIN_sex, 
        NEXTSTRAIN_age, 
        NEXTSTRAIN_host, 
    ]
     
DICT_NEXTSTRAIN_default = {
        NEXTSTRAIN_virus : "ncov",
        NEXTSTRAIN_region : "Europe",
        NEXTSTRAIN_gisaid_epi_isl : "?",
        NEXTSTRAIN_genbank_accession : "?",
        NEXTSTRAIN_segment : "genome",
        NEXTSTRAIN_sex : "?",
        NEXTSTRAIN_age : "?",
        NEXTSTRAIN_host : "?",
    }

### if None pass
DICT_INSAFLU_to_NEXTSTRAIN = {
        'id' : NEXTSTRAIN_strain,
        'collection date': NEXTSTRAIN_date,
        'onset date': NEXTSTRAIN_date,
        'lab reception date': NEXTSTRAIN_date,
        'data set' : None,
        'fastq2' : None,
        'fastq1' : None,
    }


DICT_NEXTSTRAIN_to_INSAFLU = {
        NEXTSTRAIN_strain : ['id'],
        NEXTSTRAIN_date : ['collection date', 'onset date', 'lab reception date'],
}

class MetaRow(object):
    
    def __init__(self, seq_name_consensus, row, consensus_length):
        self.seq_name_consensus = seq_name_consensus
        self.row = row
        self.consensus_length = consensus_length

class Metadata(object):
    
    def __init__(self, header):
        '''
        Constructor
        '''
        self.header = header
        self.dt_header = dict(zip(header, list(range(0, len(header)))))
        self.vect_rows_id = []
        self.dt_rows_id = {}
        
    def add_metadata(self, project_sample_pk, seq_name_consensus, row, consensus_length):
        
        if not project_sample_pk in self.dt_rows_id:
            self.dt_rows_id[project_sample_pk] = MetaRow(seq_name_consensus, row, consensus_length)
    
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
    
    def get_vect_out_nextstrain(self, vect_header_out, dt_header_normal_out, csv_writer):
        """
        :param vect_header_out all headers out
        :param dt_header_normal_out keys that are present only in INSAFLu list files
        :param csv_writer 
        """
    
        dt_out_id_project_sample = {}
        count = 0
        for project_sample_pk in self.dt_rows_id:
            if project_sample_pk in dt_out_id_project_sample: continue
            dt_out_id_project_sample[project_sample_pk] = 1
            ## NEXTSTRAIN_strain
            vect_out = [self.dt_rows_id[project_sample_pk].seq_name_consensus]

            dt_out_header = {}
            for column in vect_header_out:
                if column == NEXTSTRAIN_strain: continue    ## already out
                if column in dt_out_header: continue        ## already out, can be synonymous
                
                ## exception
                if column == NEXTSTRAIN_length:
                    if self.dt_header.get(column, -1) > 0 and \
                        len(self.dt_rows_id[project_sample_pk].row[self.dt_header[column]]) > 0:
                        vect_out.append(self.dt_rows_id[project_sample_pk].row[self.dt_header[column]])
                    elif not self.dt_rows_id[project_sample_pk].consensus_length is None:
                        vect_out.append(str(self.dt_rows_id[project_sample_pk].consensus_length))
                    else: vect_out.append('?')
                    continue
                
                ### test synonymous
                for column_insaflu in DICT_NEXTSTRAIN_to_INSAFLU.get(column, []):     
                    if self.dt_header.get(column, -1) > 0 and \
                        len(self.dt_rows_id[project_sample_pk].row[self.dt_header[column_insaflu]]) > 0:
                        vect_out.append(self.dt_rows_id[project_sample_pk].row[self.dt_header[column_insaflu]])
                        for column_insaflu in DICT_NEXTSTRAIN_to_INSAFLU[column]: dt_out_header[column_insaflu] = 1
                        break
                
                ## test default NEXTstrain columns names
                if column in DICT_NEXTSTRAIN_default:
                    if self.dt_header.get(column, -1) > 0 and \
                        len(self.dt_rows_id[project_sample_pk].row[self.dt_header[column]]) > 0:
                        vect_out.append(self.dt_rows_id[project_sample_pk].row[self.dt_header[column]])
                    else:
                        vect_out.append(DICT_NEXTSTRAIN_default[column])
                    dt_out_header[column] = 1
                    continue
                
                ### regular INSAFLU
                if column in self.dt_header: vect_out.append(self.dt_rows_id[project_sample_pk].row[self.dt_header[column]])
                else: vect_out.append("?")
            count += 1
            csv_writer.writerow(['?' if len(_) == 0 else _ for _ in vect_out])
        return count
    
    

class Reference(object):
    
    def __init__(self, file_name):
        self.file_name = file_name
        self._read_file()

    def _read_file(self):

        self.dt_out_rows = {}
        with open(self.file_name) as handle_in: 
            csv_reader = csv.reader(handle_in, delimiter=Constants.SEPARATOR_TAB)
            for line, row in enumerate(csv_reader):
                if line == 0: self.dt_header = dict(zip(row, list(range(0, len(row)))))
                else: self.dt_out_rows[row[0]] = row
        ## done
               
    def save_out_nextstrain(self, csv_writer, vect_header_out):
        """
        """
        
        count = 0
        for ref_id in self.dt_out_rows:
            vect_out = []
            for index, column in enumerate(VECT_NEXTSTRAIN_manatory):
                if (self.dt_header.get(column, -1) < 0): vect_out.append(DICT_NEXTSTRAIN_default.get(column, '?'))
                else: vect_out.append(self.dt_out_rows[ref_id][self.dt_header[column]])
            vect_out += ['?'] * (len(vect_header_out) - index - 1)   
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
        
    def add_metadata(self, project_pk, project_sample_pk, seq_name_consensus, row, consensus_length = None):
        """
        add metadata for a specific project,project_sample
        """
        if project_pk in self.dt_project:
            self.dt_project[project_pk].add_metadata(project_sample_pk, seq_name_consensus, row, consensus_length)
        
        
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


    def _get_header_nextstrain(self):
        """
        return header
        ## CollectExtraData.HEADER_SAMPLE_OUT_CSV_simple
        """
        regular_header = self._get_header()
        self.vect_header_out = []
        self.dt_header_normal_out = {}
        
        ## mandatory fields
        self.vect_header_out = VECT_NEXTSTRAIN_manatory.copy()
        
        ## try to find others in the other file
        for regular_name in regular_header:
            if regular_name in DICT_INSAFLU_to_NEXTSTRAIN: continue
            if regular_name in self.vect_header_out: continue
            self.vect_header_out.append(regular_name)
            self.dt_header_normal_out[regular_name] = 1
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


    def save_rows_nextstrain(self, csv_writer, reference_tsv):
        """
        create file for nextstrain
        """
        count = 0 ## number of rows saved
        
        ## save header
        csv_writer.writerow(self._get_header_nextstrain())
        
        ## read NextStrain reference 
        reference = Reference(reference_tsv)
        
        ### save reference
        count += reference.save_out_nextstrain(csv_writer, self.vect_header_out)
        
        ## save data
        for key_metadata in self.dt_project:
            count += self.dt_project[key_metadata].get_vect_out_nextstrain(
                self.vect_header_out, self.dt_header_normal_out, csv_writer) 
        return count
        
            