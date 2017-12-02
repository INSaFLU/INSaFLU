'''
Created on Nov 1, 2017

@author: mmp
'''

from utils.utils import Utils
from django.utils.translation import ugettext_lazy as _

class ParseOutFiles(object):
	'''
	classdocs
	'''

	GENE = 'Gene'
	COVERAGE = 'Coverage'
	IDENTITY = 'Identity'
	TYPE = 'Type'
	ACCESSION = 'Accession'

	utils = Utils()
	
	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	def parse_abricate_file(self, file_name):
		"""
		#FILE	 SEQUENCE	 START   END	 GENE	 COVERAGE	 COVERAGE_MAP	 GAPS  %COVERAGE  %IDENTITY  DATABASE   ACCESSION  PRODUCT
		6159.fna  NC_017338.1  39177   41186   mecA_15  1-2010/2010  ===============  0/0   100.00	 100.000	resfinder  AB505628   n/a
		6159.fna  NC_017338.1  727191  728356  norA_1   1-1166/1167  ===============  0/0   99.91	  92.367	 resfinder  M97169	 n/a
		6159.fna  NC_017339.1  10150   10995   blaZ_32  1-846/846	===============  0/0   100.00	 100.000	resfinder  AP004832   betalactamase
		
		Parse out abricate files 
		"""
		handle = open(file_name)
		b_fisrt_line_found = False 
		vect_out = []
		for line in handle:
			sz_temp = line.strip()
			if (len(sz_temp) == 0): continue
			if (sz_temp.find('#FILE') == 0): b_fisrt_line_found = True
			elif (b_fisrt_line_found):
				lst_split = line.split('\t')
				if (len(lst_split) != 13): continue
				dt_data = {}
				dt_data[self.GENE] = lst_split[4]
				if self.utils.is_float(lst_split[8]): dt_data[self.COVERAGE] = float(lst_split[8])
				else: raise ValueError(_("Value must be float '" + lst_split[8] + "'"))
				if self.utils.is_float(lst_split[8]): dt_data[self.IDENTITY] = float(lst_split[9])
				else: raise ValueError(_("Value must be float '" + lst_split[9] + "'"))
				dt_data[self.TYPE] = lst_split[10]
				dt_data[self.ACCESSION] = lst_split[11]
				vect_out.append(dt_data)
		handle.close()
		
		## order by coverage and identity
		return sorted(vect_out, reverse=True, key=lambda k: "%03d %03d" % (int(k[self.COVERAGE]) , int(k[self.IDENTITY])))

