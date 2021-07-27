'''
Created on Nov 1, 2017

@author: mmp
'''

from utils.utils import Utils
from django.utils.translation import ugettext_lazy as _
from constants.constants import Constants, FileExtensions
import csv, os

class ParseOutFiles(object):
	'''
	classdocs
	'''

	## header tab file
	HEADER_TAB_FILE 		= "CHROM	POS	TYPE	REF	ALT	FREQ	COVERAGE	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	LOCUS_TAG	GENE	PRODUCT"
	HEADER_TAB_FILE_old 	= "CHROM	POS	TYPE	REF	ALT	FREQ	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	LOCUS_TAG	GENE	PRODUCT"
	### Add COVERAGE after FREQ
	HEADER_TAB_FILE_after_change = "FREQ"
	HEADER_TAB_FILE_after_change_add_fields = 1
	
	HEADER_TAB_FILE_WRITE = "ID	CHROM	POS	TYPE	REF	ALT	FREQ	COVERAGE	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	NT CHANGE	AA CHANGE	AA CHANGE ALT	LOCUS_TAG	GENE	PRODUCT"
	HEADER_TAB_FILE_WRITE_WITHOUT_SAMPLE_ID = "CHROM	POS	TYPE	REF	ALT	FREQ	COVERAGE	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	NT CHANGE	AA CHANGE	AA CHANGE ALT	LOCUS_TAG	GENE	PRODUCT"
	
	HEADER_TAB_FILE_WRITE_snippy_expanded 	= "ID	CHROM	POS	TYPE	REF	ALT	FREQ	COVERAGE	EVIDENCE	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	NT CHANGE	AA CHANGE	AA CHANGE ALT	LOCUS_TAG	GENE	PRODUCT	VARIANTS IN INCOMPLETE LOCUS"
	### three types of header
	HEADER_TAB_FILE_snippy_changed 			= 		"CHROM	POS	TYPE	REF	ALT	FREQ	COVERAGE	EVIDENCE	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	NT CHANGE	AA CHANGE	AA CHANGE ALT	LOCUS_TAG	GENE	PRODUCT	VARIANTS IN INCOMPLETE LOCUS"
	HEADER_TAB_FILE_snippy_changed_old_1	= 		"CHROM	POS	TYPE	REF	ALT	FREQ	COVERAGE	EVIDENCE	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	NT CHANGE	AA CHANGE	LOCUS_TAG	GENE	PRODUCT	VARIANTS IN INCOMPLETE LOCUS"
	HEADER_TAB_FILE_snippy_changed_old_2	= 		"CHROM	POS	TYPE	REF	ALT	EVIDENCE	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	NT CHANGE	AA CHANGE	LOCUS_TAG	GENE	PRODUCT	VARIANTS IN INCOMPLETE LOCUS"
	### Add FREQ and COVERAGE after ALT, IMPORTANT, first the most distants
	HEADER_TAB_FILE_snippy_after_change_old_2 = "ALT"
	HEADER_TAB_FILE_snippy_after_change_add_fields_old_2 = 2
	
	### original file header that came from snippy
	HEADER_TAB_FILE_snippy_original = "CHROM	POS	TYPE	REF	ALT	FREQ	COVERAGE	EVIDENCE	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	LOCUS_TAG	GENE	PRODUCT"

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
	
	def get_pos_from_header(self, vect_to_search, value_to_search):
		"""
		"""
		count = 0
		for data_ in vect_to_search.split():
			if (value_to_search == data_): return count + 1
			count += 1
		return 0

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
		:out dictonary { 'Seq.id' : [dt_data, dt_data2, ...], 'Seq.id1' : [dt_data5, dt_data6, ...]
		"""
		clean_abricate_file = self.utils.get_temp_file('clean_abricate', FileExtensions.FILE_TXT)
		with open(file_name) as handle, open(clean_abricate_file, 'w') as handle_new:
			b_fisrt_line_found = False 
			dict_data_out = {}
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
					if (dt_data[self.SEQ_NAME] in dict_data_out): dict_data_out[dt_data[self.SEQ_NAME]].append(dt_data)
					else: dict_data_out[dt_data[self.SEQ_NAME]] = [dt_data]
		
		### order the geneId by coverage
		for key in dict_data_out:
			dict_data_out[key] = sorted(dict_data_out[key], reverse=True, key=lambda k: (k[self.COVERAGE],
															k[self.IDENTITY]) )

		
		### order the key 
		## order by coverage and identity, bigger better...
		return dict_data_out, clean_abricate_file


	def parse_tab_files(self, sample_name, file_to_parse, csv_writer, vect_type_out, vect_type_remove, limit_freq, b_add_header):
		"""
		Only used in FreeBayes
		limit_freq -> the max freq to accept the variation 
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
			b_header_old = False
			for line in handle_to_process:
				sz_temp = line.strip()
				if (len(sz_temp) == 0): continue
				if (sz_temp == self.HEADER_TAB_FILE): b_header = True
				elif (sz_temp == self.HEADER_TAB_FILE_old): b_header_old = True
				elif (b_header or b_header_old): 	## start processing data
					lst_data = sz_temp.split('\t')
					if (len(lst_data) > 5):
						lst_freq_data = lst_data[5].split(',')
						lst_type_var = lst_data[2].split(',')
						b_exist = False
						for i in range(0, len(lst_type_var)):
							for value in vect_type_remove: 	## don't print specific var types
								if (value in lst_type_var):
									b_exist = True
									break
							if (b_exist): break	### exist, don't print
							
							if ((lst_type_var[i] in vect_type_out or len(vect_type_out) == 0) and self.utils.is_float(lst_freq_data[i])\
										and float(lst_freq_data[i]) < limit_freq):
								vect_to_write = []
								if (sample_name != None): vect_to_write = [sample_name]
								## transform 'synonymous_variant c.981A>G p.Glu327Glu' to ["synonymous_variant", "c.981A>G", "p.Glu327Glu"]
								if (b_header_old):
									if (len(lst_data) > 10):
										vect_to_write.extend(lst_data[:10])
										vect_to_write.extend(lst_data[10].split(' '))
										vect_to_write.extend([self.utils.parse_amino_HGVS_code(vect_to_write[-1])])
										vect_to_write.extend(lst_data[11:])
									else: vect_to_write.extend(lst_data)
									
									### need to add empty fields
									for _ in range(self.HEADER_TAB_FILE_after_change_add_fields):
										vect_to_write.insert(self.get_pos_from_header(self.HEADER_TAB_FILE_old, self.HEADER_TAB_FILE_after_change) + 1, "")
								else: 
									if (len(lst_data) > 11):
										vect_to_write.extend(lst_data[:11])
										vect_to_write.extend(lst_data[11].split(' '))
										vect_to_write.extend([self.utils.parse_amino_HGVS_code(vect_to_write[-1])])
										vect_to_write.extend(lst_data[12:])
									else: vect_to_write.extend(lst_data)
								csv_writer.writerow(vect_to_write)
								count_print += 1
								break
		return count_print

	def parse_tab_files_snippy(self, sample_name, file_to_parse, csv_writer, vect_type_out, b_add_header):
		"""
		Join and Snippy/Medaka variant files
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
			b_header_old_1 = False
			b_header_old_2 = False
			for line in handle_to_process:
				sz_temp = line.strip()
				if (len(sz_temp) == 0): continue
				if (sz_temp == self.HEADER_TAB_FILE_snippy_changed): b_header = True
				elif (sz_temp == self.HEADER_TAB_FILE_snippy_changed_old_1): b_header_old_1 = True
				elif (sz_temp == self.HEADER_TAB_FILE_snippy_changed_old_2): b_header_old_2 = True
				elif (b_header or b_header_old_1 or b_header_old_2): 	## start processing data
					lst_data = sz_temp.split('\t')
					if (len(lst_data) > 5):
						lst_type_var = lst_data[2].split(',')
						for i in range(0, len(lst_type_var)):
							if (len(vect_type_out) == 0 or lst_type_var[i] in vect_type_out):
								vect_to_write = [sample_name]
								### add HGMD amino p.V323S 
								if b_header_old_2: position = 12	## oldest versions
								elif b_header_old_1: position = 14	## old versions
								if (not b_header and len(lst_data) > position):
									amino_p = lst_data[position]
									lst_transformed = []
									for amino in amino_p.split(','):
										lst_transformed.append(self.utils.parse_amino_HGVS_code(amino))
									lst_data = lst_data[:position+1] + [",".join(lst_transformed)] + lst_data[position+1:]
								vect_to_write.extend(lst_data)
								break
						
						### test header old
						if b_header_old_2:
							for _ in range(self.HEADER_TAB_FILE_snippy_after_change_add_fields_old_2):
								vect_to_write.insert(self.get_pos_from_header(self.HEADER_TAB_FILE_snippy_changed_old_2,
									self.HEADER_TAB_FILE_snippy_after_change_old_2) + 1, "")
						csv_writer.writerow(vect_to_write)
		return True


	def add_variants_in_incomplete_locus(self, file_to_process, coverage):
		"""
		Used in snippy and Medaka
		add variants_in_incomplete_locus in low coverage locus
		AND -> transform 'synonymous_variant c.981A>G p.Glu327Glu' to ["synonymous_variant", "c.981A>G", "p.Glu327Glu"] 
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
						if (len(lst_data) > 12):
							## transform 'synonymous_variant c.981A>G p.Glu327Glu' to ["synonymous_variant", "c.981A>G", "p.Glu327Glu"]
							vect_to_write.extend(lst_data[:12])
							vect_to_write.extend(lst_data[12].split(' '))
							vect_to_write.extend([self.utils.parse_amino_HGVS_code(vect_to_write[-1])])
							vect_to_write.extend(lst_data[13:])
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


	def add_amino_single_letter_code(self, vcf_file):
		"""
		add single letter Amino change to VCF files 
		:param	vcf_file in
		:out	file with amino transformation  
		"""
		### AB=0;ABP=0;AC=2;AF=1;AN=2;AO=289;CIGAR=1X;DP=289;DPB=289;DPRA=0;EPP=7.70639;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=0;NS=1;NUMALT=1;ODDS=204.925;PAIRED=0;PAIREDR=0;
		###		PAO=0;PQA=0;PQR=0;PRO=0;QA=10861;QR=0;RO=0;RPL=151;RPP=4.28012;RPPR=0;RPR=138;RUN=1;SAF=133;SAP=6.98507;SAR=156;SRF=0;SRP=0;SRR=0;TYPE=snp;
		###		ANN=T|synonymous_variant|LOW|Exon_PB1_1_2274|Gene_Exon_PB1_1_2274|transcript|Transcript_Exon_PB1_1_2274|Coding|1/1|c.876C>T|p.Asn292Asn|876/2274|876/2274|292/757||;FREQ=100
		temp_file = self.utils.get_temp_file("add_amino_", ".vcf")
		with open(vcf_file) as handle_in, open(temp_file, 'w') as handle_out:
			for line in handle_in:
				sz_temp = line.strip()
				if (len(sz_temp) == 0 or sz_temp[0] == '#'): 
					
					### check annotation info, because now they have two HGVS.p
					if line.startswith("##INFO=<ID=ANN,Number="): line = line.replace('| HGVS.p', '| HGVS.p | HGVS.p') 
					handle_out.write(line)
					continue
				lst_data = sz_temp.split()
				if len(lst_data) > 7:
					lst_tag = lst_data[7].split(";")
					for pos_1, tag in enumerate(lst_tag):
						if (tag.startswith('ANN=')):
							lst_value = tag.split("|")
							for pos_2, value_ in enumerate(lst_value):
								if value_.startswith("p."):
									parse_amino = self.utils.parse_amino_HGVS_code(value_)
									if (not parse_amino is None): lst_value.insert(pos_2 + 1, parse_amino)
									break
							lst_tag[pos_1] = "|".join(lst_value)
							lst_data[7] = ";".join(lst_tag)
							break
					handle_out.write("\t".join(lst_data) + "\n")
				else: handle_out.write(line)
		return temp_file
		
	
		
		
		