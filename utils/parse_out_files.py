'''
Created on Nov 1, 2017

@author: mmp
'''

from utils.utils import Utils
from django.utils.translation import ugettext_lazy as _
from constants.constants import Constants, FileExtensions
from constants.software_names import SoftwareNames
import csv, os

class ParseOutFiles(object):
	'''
	classdocs
	'''

	## header tab file
	HEADER_TAB_FILE = "CHROM	POS	TYPE	REF	ALT	FREQ	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	LOCUS_TAG	GENE	PRODUCT"
	HEADER_TAB_FILE_WRITE = "ID	CHROM	POS	TYPE	REF	ALT	FREQ	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	NT CHANGE	AA CHANGE	LOCUS_TAG	GENE	PRODUCT"
	HEADER_TAB_FILE_WRITE_WITHOUT_SAMPLE_ID = "CHROM	POS	TYPE	REF	ALT	FREQ	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	NT CHANGE	AA CHANGE	LOCUS_TAG	GENE	PRODUCT"
	
	HEADER_TAB_FILE_WRITE_snippy_expanded = "ID	CHROM	POS	TYPE	REF	ALT	EVIDENCE	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	NT CHANGE	AA CHANGE	LOCUS_TAG	GENE	PRODUCT	VARIANTS IN INCOMPLETE LOCUS"
	HEADER_TAB_FILE_snippy_changed = "CHROM	POS	TYPE	REF	ALT	EVIDENCE	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	NT CHANGE	AA CHANGE	LOCUS_TAG	GENE	PRODUCT	VARIANTS IN INCOMPLETE LOCUS"
	HEADER_TAB_FILE_snippy_original = "CHROM	POS	TYPE	REF	ALT	EVIDENCE	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	LOCUS_TAG	GENE	PRODUCT"
	
	GENE = 'Gene'
	COVERAGE = 'Coverage'
	IDENTITY = 'Identity'
	TYPE = 'Type'
	ACCESSION = 'Accession'
	SEQ_NAME = 'Seq_Name'

	utils = Utils()
	
	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	def parse_abricate_file(self, file_name, file_name_out, clean_hits_below_value):
		"""
		#FILE	 SEQUENCE	 START   END	 GENE	 COVERAGE	 COVERAGE_MAP	 GAPS  %COVERAGE  %IDENTITY  DATABASE   ACCESSION  PRODUCT
		6159.fna  NC_017338.1  39177   41186   mecA_15  1-2010/2010  ===============  0/0   100.00	 100.000	resfinder  AB505628   n/a
		6159.fna  NC_017338.1  727191  728356  norA_1   1-1166/1167  ===============  0/0   99.91	  92.367	 resfinder  M97169	 n/a
		6159.fna  NC_017339.1  10150   10995   blaZ_32  1-846/846	===============  0/0   100.00	 100.000	resfinder  AP004832   betalactamase
		
		#FILE	SEQUENCE	START	END	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT
		/tmp/insaFlu/insa_flu_85008740/contigs.fasta	NODE_1_length_3374_cov_697.364266	2372	3353	A	1-982/982	===============	0/0	100.00	99.69	influenza_type	XXXX	N/A
		/tmp/insaFlu/insa_flu_85008740/contigs.fasta	NODE_27_length_226_cov_1.035088	1	226	N2	41-266/1410	===............	0/0	16.03	96.46	influenza_subtype	XXXXXX	N/A

		Parse out abricate files
		return also a file with coverage below clean_hits_below_value removed
		Can be empty
		"""
		clean_abricate_file = self.utils.get_temp_file('clean_abricate', FileExtensions.FILE_TXT)
		with open(file_name) as handle, open(clean_abricate_file, 'w') as handle_new:
			b_fisrt_line_found = False 
			vect_out = []
			for line in handle:
				sz_temp = line.strip()
				if (len(sz_temp) == 0): continue
				if (sz_temp.find('#FILE') == 0):
					handle_new.write(line) 
					b_fisrt_line_found = True
				elif (b_fisrt_line_found):
					lst_split = line.split('\t')
					if (len(lst_split) != 13): continue
					
					### clean SoftwareNames.SOFTWARE_SPAdes_CLEAN_HITS_BELLOW_VALUE
					lst_coverage = lst_split[1].split('_')
					if (len(lst_coverage) > 4 and self.utils.is_float(lst_coverage[-1]) and\
							(float(lst_coverage[-1]) < clean_hits_below_value)):
						continue
					
					## create new file without the coverage below clean_hits_below_value
					handle_new.write('{}\t{}'.format(file_name_out, "\t".join(lst_split[1:])))
					
					dt_data = {}
					dt_data[self.GENE] = lst_split[4]
					if self.utils.is_float(lst_split[8]): dt_data[self.COVERAGE] = float(lst_split[8])
					else: raise ValueError(_("Value must be float '" + lst_split[8] + "'"))
					if self.utils.is_float(lst_split[8]): dt_data[self.IDENTITY] = float(lst_split[9])
					else: raise ValueError(_("Value must be float '" + lst_split[9] + "'"))
					dt_data[self.TYPE] = lst_split[10]
					dt_data[self.ACCESSION] = lst_split[11]
					dt_data[self.SEQ_NAME] = lst_split[1]
					vect_out.append(dt_data)
		
		## order by coverage and identity
		return sorted(vect_out, reverse=True, key=lambda k: "%03d %03d" % (int(k[self.COVERAGE]) , int(k[self.IDENTITY]))), clean_abricate_file


	def parse_tab_files(self, sample_name, file_to_parse, csv_writer, vect_type_out, vect_type_remove, limit_freq, b_add_header):
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
		vect_type_remove = ['ins', 'del'] , for example
		"""
		
		count_print = 0
		if (b_add_header): 
			if (sample_name == None): csv_writer.writerow(self.HEADER_TAB_FILE_WRITE_WITHOUT_SAMPLE_ID.split('\t'))
			else: csv_writer.writerow(self.HEADER_TAB_FILE_WRITE.split('\t'))
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
						b_exist = False
						for i in range(0, len(lst_type_var)):
							for value in vect_type_remove: 
								if (value in lst_type_var):
									b_exist = True
									break
							if (b_exist): break	### exist, don't print
							
							if ((lst_type_var[i] in vect_type_out or len(vect_type_out) == 0) and self.utils.is_float(lst_freq_data[i])\
										and float(lst_freq_data[i]) < limit_freq):
								vect_to_write = []
								if (sample_name != None): vect_to_write = [sample_name]
								if (len(lst_data) > 10):
									vect_to_write.extend(lst_data[:10])
									vect_to_write.extend(lst_data[10].split(' '))
									vect_to_write.extend(lst_data[11:])
								else: vect_to_write.extend(lst_data)
								csv_writer.writerow(vect_to_write)
								count_print += 1
								break

		return count_print

	def parse_tab_files_snippy(self, sample_name, file_to_parse, csv_writer, vect_type_out, b_add_header):
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
		
		if (b_add_header): csv_writer.writerow(self.HEADER_TAB_FILE_WRITE_snippy_expanded.split('\t'))
		with open(file_to_parse) as handle_to_process:
			b_header = False
			for line in handle_to_process:
				sz_temp = line.strip()
				if (len(sz_temp) == 0): continue
				if (sz_temp == self.HEADER_TAB_FILE_snippy_changed): b_header = True
				elif (b_header): 	## start processing data
					lst_data = sz_temp.split('\t')
					if (len(lst_data) > 5):
						lst_type_var = lst_data[2].split(',')
						for i in range(0, len(lst_type_var)):
							if (len(vect_type_out) == 0 or lst_type_var[i] in vect_type_out):
								vect_to_write = [sample_name]
								vect_to_write.extend(lst_data)
								csv_writer.writerow(vect_to_write)
								break
		return True


	def add_variants_in_incomplete_locus(self, file_to_process, coverage):
		"""
		add variants_in_incomplete_locus in low coverage locus
		"""
		out_file = self.utils.get_temp_file('variants_in_incomplete_locus', FileExtensions.FILE_TSV)
		with open(out_file, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=Constants.SEPARATOR_TAB, quotechar='"', quoting=csv.QUOTE_NONE)
			csv_writer.writerow(self.HEADER_TAB_FILE_snippy_changed.split('\t'))
			n_size_array = len(self.HEADER_TAB_FILE_snippy_changed.split('\t'))
			with open(file_to_process) as handle_to_process:
				b_header = False
				for line in handle_to_process:
					sz_temp = line.strip()
					if (len(sz_temp) == 0): continue
					if (sz_temp == self.HEADER_TAB_FILE_snippy_original): b_header = True
					elif (b_header): 	## start processing data
						lst_data = sz_temp.split('\t')
						vect_to_write = []
						if (len(lst_data) > 10):
							vect_to_write.extend(lst_data[:10])
							vect_to_write.extend(lst_data[10].split(' '))
							vect_to_write.extend(lst_data[11:])
							if (len(lst_data) > 1 and coverage.exist_this_element(lst_data[0]) and\
									(not coverage.is_100_more_0(lst_data[0]) or not coverage.is_100_more_9(lst_data[0]))):
								vect_to_write.append('yes')
						else:
							vect_to_write.extend(lst_data)
							if (len(lst_data) > 1 and coverage.exist_this_element(lst_data[0]) and\
									(not coverage.is_100_more_0(lst_data[0]) or not coverage.is_100_more_9(lst_data[0]))):
								vect_to_write.extend([''] * (n_size_array - len(lst_data) - 1))
								vect_to_write.append('yes')
						csv_writer.writerow(vect_to_write)

		### copy to the origin to the result file
		if (os.path.exists(out_file)):
			if (os.path.getsize(out_file) > 0): self.utils.copy_file(out_file, file_to_process)
			os.unlink(out_file)



