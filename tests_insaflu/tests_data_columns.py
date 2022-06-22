'''
Created on 13/06/2022

@author: mmp
'''
import unittest, os, filecmp, csv
from datetime import date
from django.conf import settings
from utils.data_columns import DataColumns
from utils.utils import Utils, FileExtensions
from utils.software import Software
from constants.constantsTestsCase import ConstantsTestsCase
from constants.constants import Constants

class Test(unittest.TestCase):

    utils = Utils()
    software = Software()

    def setUp(self):
        self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
        pass
    
    def test_data_columns(self):
        
        data_columns = DataColumns()
        
        data_columns.add_header(1, ['id', '1', '2'])
        data_columns.add_metadata(1, 1, "xpto", ['id', 'a', 'b'])
        data_columns.add_metadata(1, 1, "xpto", ['id', 'a', 'b'])
        data_columns.add_metadata(1, 3, "xpto3", ['id', 'a3', 'b3'])       
        
        data_columns.add_header(2, ['id', '2', '3'])
        data_columns.add_metadata(2, 4, "xpto2", ['id', 'a1', 'b1'])  

        expected_file_coverage = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_DATASET_FILES, "insa_flu_dataset_list_output.csv")
        
        out_file = self.utils.get_temp_file('dataset_list', FileExtensions.FILE_CSV)
        with open(out_file, 'w', newline='') as handle_out:
            csv_writer = csv.writer(handle_out, delimiter=Constants.SEPARATOR_COMMA, quotechar='"',
                        quoting=csv.QUOTE_MINIMAL)
            data_columns.save_rows(csv_writer)
        self.assertTrue(os.path.exists(out_file))
        self.assertTrue(filecmp.cmp(out_file, expected_file_coverage))
        if (os.path.exists(out_file)): os.unlink(out_file)


    def test_data_columns_nextstrain(self):
        data_columns = DataColumns()
        
        data_columns.add_header(1, ['id', '1', '2', 'collection date'])
        data_columns.add_metadata(1, 1, "xpto", ['id', 'a', 'b', '2019-12-26'])
        data_columns.add_metadata(1, 1, "xpto", ['id', 'a', 'b', '2019-12-25'])
        data_columns.add_metadata(1, 3, "xpto3", ['id', 'a3', 'b3', '2019-12-21'])       
        data_columns.add_header(2, ['id', '2', '3', 'onset date'])
        data_columns.add_metadata(2, 4, "xpto2", ['id', 'a1', 'b1', '2018-12-21'])
        data_columns.add_header(3, ['id', '2', '3', 'date'])
        data_columns.add_metadata(3, 5, "xptow2", ['id', 'a1', 'b1', '2015-12-21'])
        data_columns.add_metadata(3, 6, "xptow3", ['id', 'a1', 'b1', ''])

        reference_nextstrain = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_DATASET_FILES, "references_metadata.tsv")
        expected_file_nextstrain = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_DATASET_FILES, "insa_flu_dataset_output_nextstrain.csv")
        
        out_file = self.utils.get_temp_file('dataset_list', FileExtensions.FILE_CSV)
        with open(out_file, 'w', newline='') as handle_out:
            csv_writer = csv.writer(handle_out, delimiter=Constants.SEPARATOR_COMMA, quotechar='"',
                        quoting=csv.QUOTE_MINIMAL)
            data_columns.save_rows_nextstrain(csv_writer, reference_nextstrain)
        self.assertTrue(os.path.exists(out_file))
        self.software.dos_2_unix(out_file)
        
        ## update date
        temp_file_compare = self.utils.get_temp_file_from_dir(os.path.dirname(out_file), "nextstrain_compare", ".txt")
        vect_data = self.utils.read_text_file(expected_file_nextstrain)
        with open(temp_file_compare, 'w') as handle_write:
            for line in vect_data:
                handle_write.write(line.strip().replace('DATE_TO_REPLACE', date.today().strftime(settings.DATE_FORMAT_FOR_SHOW)) + "\n")
            
        self.assertTrue(filecmp.cmp(out_file, temp_file_compare))
        if (os.path.exists(out_file)): os.unlink(out_file)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_data_columns']
    unittest.main()