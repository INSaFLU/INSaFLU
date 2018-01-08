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

	## header tab file
	HEADER_TAB_FILE = "CHROM	POS	TYPE	REF	ALT	FREQ	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	LOCUS_TAG	GENE	PRODUCT"
	HEADER_TAB_FILE_WRITE = "CHROM	POS	TYPE	REF	ALT	FREQ	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	EFFECT.c	EFFECT.p	LOCUS_TAG	GENE	PRODUCT"
	
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


	def parse_tab_files(self, file_to_parse, handle_out, vect_type_out, limit_freq, b_add_header):
		"""
		limit_freq -> the tre max freq to accept the variation
		process tab files
		possible in variation type
				snp	Single Nucleotide Polymorphism	A => T
				mnp	Multiple Nuclotide Polymorphism	GC => AT
				ins	Insertion	ATT => AGTT
				del	Deletion	ACGG => ACG
				complex	Combination of snp/mnp
		vect_count_type = ['snp', 'ins']
		"""
		
		if (b_add_header): handle_out.write(self.HEADER_TAB_FILE_WRITE + "\n")
		with open(file_to_parse) as handle_to_process:
			b_header = False
			for line in handle_to_process:
				sz_temp = line.strip()
				if (len(sz_temp) == 0): continue
				if (sz_temp == self.HEADER_TAB_FILE): b_header = True
				elif (b_header): 	## start processing data
					lst_data = sz_temp.split('\t')
					if (len(lst_data) > 5):
						lst_freq_data = lst_data[5].split(',')
						lst_type_var = lst_data[2].split(',')
						for i in range(0, len(lst_type_var)):
							if (lst_type_var[i] in vect_type_out and self.utils.is_float(lst_freq_data[i]) and float(lst_freq_data[i]) < limit_freq):
								if (len(lst_data) > 10):
									handle_out.write('\t'.join(lst_data[:10]) + '\t' + lst_data[10].replace(' ', '\t') + '\t' + '\t'.join(lst_data[11:]) + "\n")
								else: handle_out.write(sz_temp + "\n")
								break

		return True


