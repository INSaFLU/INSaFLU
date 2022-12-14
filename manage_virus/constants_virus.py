'''
Created on Jan 21, 2018

@author: mmp
'''

class ConstantsVirus(object):
	'''
	Identification from file
		files_helpful/Criteria_flag_mixed_infections_alerts_v2.xls
	Source
		type_identification/db_influenza_typing_vXXXX.fasta
	'''
	
	
	##############################33
	####
	####	TYPE/GENUS/etc... identification
	####
	####	first level
	####	influenza_type~~~B~~~AF100378,
	####	coronavirus_genus~~~BetaCoV~~~M_gene_MN908947
	####
	SEQ_VIRUS_TYPE = "Type"				### A
	SEQ_VIRUS_GENUS = "Genus"
	SEQ_VIRUS_GENOTYPE = "Genotype"
	VECT_SEQ_VIRUS_TYPE = [SEQ_VIRUS_GENOTYPE.lower(), SEQ_VIRUS_TYPE.lower(),
				SEQ_VIRUS_GENUS.lower()]
	
	### sub type identification
	##############################33
	####
	####	SUBTYPE/LINEAGE/HUMAN/etc... identification
	####
	####	second level
	####	influenza_subtype~~~N11~~~CY125947
	####	influenza_subtype~~~H2~~~L11142
	####	influenza_lineage~~~Victoria~~~M58428
	####	coronavirus_human~~~2019_nCoV~~~S_gene_MN908947
	####	coronavirus_bat~~~bat_RaTG13~~~S_gene_MN996532
	####	monkeypox_species~~~MPVX~~~rpo18_gene_ON649720
	####
	SEQ_VIRUS_SUB_TYPE = "Subtype"
	SEQ_VIRUS_LINEAGE = "Lineage"		### B
	SEQ_VIRUS_HUMAN = "Human"
	SEQ_VIRUS_BAT = "Bat"
	SEQ_VIRUS_SPECIES = "Species"

	VECT_SEQ_VIRUS_SUB_TYPE = [SEQ_VIRUS_SUB_TYPE.lower(), SEQ_VIRUS_LINEAGE.lower(),
					SEQ_VIRUS_HUMAN.lower(), SEQ_VIRUS_BAT.lower(),
					SEQ_VIRUS_SPECIES.lower()]
	
	## used to set classification in Sample.type_subtype
	VECT_SUB_TYPE_other_than_influenza = [SEQ_VIRUS_HUMAN, SEQ_VIRUS_BAT,
					SEQ_VIRUS_SPECIES]
	
	### all possible tags available
	### this tags came from the file db/type_identification/db_influenza_typing_v2.fasta
	VECT_ALL_POSSIBLE_TAGS_lower = VECT_SEQ_VIRUS_TYPE + VECT_SEQ_VIRUS_SUB_TYPE
	VECT_ALL_POSSIBLE_TAGS = [SEQ_VIRUS_TYPE, SEQ_VIRUS_GENUS, SEQ_VIRUS_SUB_TYPE, SEQ_VIRUS_LINEAGE,
					SEQ_VIRUS_HUMAN, SEQ_VIRUS_BAT, SEQ_VIRUS_SPECIES, SEQ_VIRUS_GENOTYPE]
	
	TYPE_A = 'A'
	TYPE_B = 'B'
	TYPE_BetaCoV = 'BetaCoV'
	
	### SEQ_VIRUS_TYPE A starts
	SUB_TYPE_STARTS_N = 'N'
	SUB_TYPE_STARTS_H = 'H'


