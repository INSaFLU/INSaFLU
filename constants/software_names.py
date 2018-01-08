'''
Created on Nov 26, 2017

@author: mmp
'''

import os
from constants.constants import Constants


class SoftwareNames(object):
	'''
	classdocs
	'''

	## dir with software
	DIR_SOFTWARE = "/usr/local/software/insaflu"
	SOFTWARE_SAMTOOLS = os.path.join(DIR_SOFTWARE, "snippy/bin/samtools")
	SOFTWARE_SAMTOOLS_name = "Samtools"
	SOFTWARE_SAMTOOLS_VERSION = "1.3"
	SOFTWARE_SAMTOOLS_PARAMETERS = ""
	SOFTWARE_BGZIP = os.path.join(DIR_SOFTWARE, "snippy/bin/bgzip")
	SOFTWARE_BGZIP_name = "bgzip"
	SOFTWARE_BGZIP_VERSION = "1.3"
	SOFTWARE_BGZIP_PARAMETERS = ""
	SOFTWARE_GZIP_name = "gzip"
	SOFTWARE_GZIP_VERSION = "1.6"
	SOFTWARE_GZIP_PARAMETERS = ""
	SOFTWARE_TABIX = os.path.join(DIR_SOFTWARE, "snippy/bin/tabix")
	SOFTWARE_TABIX_name = "tabix"
	SOFTWARE_TABIX_VERSION = "1.3"
	SOFTWARE_TABIX_PARAMETERS = ""
	SOFTWARE_IGVTOOLS = os.path.join(DIR_SOFTWARE, "IGVTools/igvtools.jar")
	SOFTWARE_IGVTOOLS_name = "igvtools"
	SOFTWARE_IGVTOOLS_VERSION = "2.3.98"
	SOFTWARE_IGVTOOLS_PARAMETERS = ""
	SOFTWARE_SPAdes = os.path.join(DIR_SOFTWARE, "SPAdes-3.11.1-Linux/bin/spades.py")
	SOFTWARE_SPAdes_name = "SPAdes" 
	SOFTWARE_SPAdes_VERSION = "3.11.1"
	SOFTWARE_SPAdes_PARAMETERS = "--meta --only-assembler"
	SOFTWARE_SPAdes_PARAMETERS_single = "--only-assembler"
	SOFTWARE_ABYSS = os.path.join(DIR_SOFTWARE, "abyss/bin/abyss-pe")
	SOFTWARE_ABYSS_name = "Abyss" 
	SOFTWARE_ABYSS_VERSION = "2.0"
	SOFTWARE_ABYSS_PARAMETERS = 'k=63'
	SOFTWARE_ABRICATE = os.path.join(DIR_SOFTWARE, "abricate/bin/abricate")
	SOFTWARE_ABRICATE_name = "Abricate"
	SOFTWARE_ABRICATE_DB = os.path.join(DIR_SOFTWARE, "abricate/db")
	SOFTWARE_ABRICATE_VERSION = "0.8-dev"
	SOFTWARE_ABRICATE_PARAMETERS = "--minid 70"
	SOFTWARE_FASTQ = os.path.join(DIR_SOFTWARE, "FastQC/fastqc")
	SOFTWARE_FASTQ_name = "FastQC"
	SOFTWARE_FASTQ_VERSION = "0.11.5"
	SOFTWARE_FASTQ_PARAMETERS = ""
	SOFTWARE_TRIMMOMATIC = os.path.join(DIR_SOFTWARE, "trimmomatic/classes/trimmomatic.jar")
	SOFTWARE_TRIMMOMATIC_name = "Trimmomatic"
	SOFTWARE_TRIMMOMATIC_VERSION = "0.27"
	SOFTWARE_TRIMMOMATIC_PARAMETERS = "SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33"
	SOFTWARE_SNIPPY = os.path.join(DIR_SOFTWARE, "snippy/bin/snippy")
	SOFTWARE_SNIPPY_name = "Snippy"
	SOFTWARE_SNIPPY_VERSION = "3.2-dev"
	SOFTWARE_SNIPPY_PARAMETERS = "--mapqual 20 --mincov 10 --minfrac 0.51"
	SOFTWARE_SNIPPY_VCF_TO_TAB = os.path.join(DIR_SOFTWARE, "snippy/bin/snippy-vcf_to_tab_add_freq")
	SOFTWARE_SNIPPY_VCF_TO_TAB_name = "Snippy-vcf_to_tab_add_freq"
	SOFTWARE_SNIPPY_VCF_TO_TAB_VERSION = "3.2-dev"
	SOFTWARE_SNIPPY_VCF_TO_TAB_PARAMETERS = ""
	SOFTWARE_SNP_EFF = os.path.join(DIR_SOFTWARE, "snippy/bin/snpEff")
	SOFTWARE_SNP_EFF_config = os.path.join(DIR_SOFTWARE, "snippy/etc/snpeff.config")
	SOFTWARE_SNP_EFF_name = "snpEff"
	SOFTWARE_SNP_EFF_VERSION = "4.3p"
	SOFTWARE_SNP_EFF_PARAMETERS = "-no-downstream -no-upstream -no-intergenic -no-utr -noStats"
	SOFTWARE_GENBANK2GFF3 = 'bp_genbank2gff3'
	SOFTWARE_GENBANK2GFF3_name = 'Genbank2gff3'
	SOFTWARE_GENBANK2GFF3_VERSION = 'Unknown'
	SOFTWARE_GENBANK2GFF3_PARAMETERS = ''
	SOFTWARE_FREEBAYES = os.path.join(DIR_SOFTWARE, "snippy/bin/freebayes")
	SOFTWARE_FREEBAYES_PARALLEL = os.path.join(DIR_SOFTWARE, "snippy/bin/freebayes-parallel")
	SOFTWARE_FREEBAYES_name = "Freebayes"
	SOFTWARE_FREEBAYES_VERSION = "v1.1.0-54-g49413aa"
	SOFTWARE_FREEBAYES_PARAMETERS = "-p 2 -q 20 -m 20 --min-coverage 100 --min-alternate-fraction 0.01 --min-alternate-count 10 -V"
	
	SOFTWARE_FASTA_GENERATE_REGIONS = os.path.join(DIR_SOFTWARE, "snippy/bin/fasta_generate_regions.py")
	SOFTWARE_FASTA_GENERATE_REGIONS_name = "Fasta Generate Regions"
	SOFTWARE_FASTA_GENERATE_REGIONS_VERSION = "1.0"
	SOFTWARE_FASTA_GENERATE_REGIONS_PARAMETERS = ""
	SOFTWARE_COVERAGE_TO_REGIONS = os.path.join(DIR_SOFTWARE, "freebayes/scripts/coverage_to_regions.py")
	SOFTWARE_COVERAGE_TO_REGIONS_name = "Coverage to Regions"
	SOFTWARE_COVERAGE_TO_REGIONS_VERSION = "1.0"
	SOFTWARE_COVERAGE_TO_REGIONS_PARAMETERS = ""
	
	SOFTWARE_BAMTOOLS = os.path.join(DIR_SOFTWARE, "bamtools/build/src/toolkit/bamtools")
	SOFTWARE_BAMTOOLS_name = "Bamtools"
	SOFTWARE_BAMTOOLS_VERSION = "2.5"
	SOFTWARE_BAMTOOLS_PARAMETERS = ""
	
	SOFTWARE_COVERAGE = "Coverage, in-house script"
	SOFTWARE_COVERAGE_name = "Coverage"
	SOFTWARE_COVERAGE_VERSION = "v1.1"
	SOFTWARE_COVERAGE_PARAMETERS = ""
	
	SOFTWARE_PROKKA = os.path.join(DIR_SOFTWARE, "prokka/bin/prokka")
	SOFTWARE_PROKKA_name = "Prokka"
	SOFTWARE_PROKKA_VERSION = "1.2"
	SOFTWARE_PROKKA_PARAMETERS = "--kingdom Viruses --locustag locus --genus Influenzavirus --species Influenzavirus --strain "\
					"ref_PREFIX_FILES_OUT --gcode " + str(Constants.TRANSLATE_TABLE_NUMBER)
	
	SOFTWARE_MAUVE = os.path.join(DIR_SOFTWARE, "mauve/progressiveMauve")
	SOFTWARE_MAUVE_name = "Mauve"
	SOFTWARE_MAUVE_VERSION = "2.4.0, Feb 13 2015"
	SOFTWARE_MAUVE_PARAMETERS = ""
	
	SOFTWARE_CONVERT = os.path.join(DIR_SOFTWARE, "scripts/convert.pl")
	SOFTWARE_CONVERT_name = "Convert"
	SOFTWARE_CONVERT_VERSION = ""
	SOFTWARE_CONVERT_PARAMETERS = ""

	SOFTWARE_MAFFT = os.path.join(DIR_SOFTWARE, "mafft-7.313-without-extensions/scripts/mafft")
	SOFTWARE_SET_ENV_MAFFT = "export MAFFT_BINARIES={}".format(os.path.join(DIR_SOFTWARE, "mafft-7.313-without-extensions/binaries"))
	SOFTWARE_MAFFT_name = "Mafft"
	SOFTWARE_MAFFT_VERSION = "7.313"
	SOFTWARE_MAFFT_PARAMETERS_TWO_SEQUENCES = "--maxiterate 1000 --localpair --preservecase --leavegappyregion"
	SOFTWARE_MAFFT_PARAMETERS_PROTEIN = "--maxiterate 1000 --localpair --preservecase --amino"
	SOFTWARE_MAFFT_PARAMETERS = "--preservecase"
	
	SOFTWARE_SEQRET = "seqret"
	SOFTWARE_SEQRET_name = "seqret (EMBOSS)"
	SOFTWARE_SEQRET_VERSION = "6.6.0.0"
	SOFTWARE_SEQRET_NEX_PARAMETERS = "-sformat fasta -osformat2 nexusnon"
	
#	SOFTWARE_FASTTREE = os.path.join(DIR_SOFTWARE, "fasttree/FastTree")
#	SOFTWARE_FASTTREE_name = "FastTree"
#	SOFTWARE_FASTTREE_VERSION = "2.1.10 SSE3"
	SOFTWARE_FASTTREE = os.path.join(DIR_SOFTWARE, "fasttree/FastTreeDbl")
	SOFTWARE_FASTTREE_name = "FastTreeDbl"
	SOFTWARE_FASTTREE_VERSION = "2.1.10 Double precision"
	SOFTWARE_FASTTREE_PARAMETERS = "-gtr -boot 1000 -nt"
	SOFTWARE_FASTTREE_PARAMETERS_PROTEIN = "-gtr -boot 1000"
	
	def __init__(self):
		'''
		Constructor
		'''
		pass

	"""
	return samtools software
	"""
	def get_samtools(self): return self.SOFTWARE_SAMTOOLS
	def get_samtools_version(self): return self.SOFTWARE_SAMTOOLS_VERSION

	"""
	return spades software
	"""
	def get_spades(self): return self.SOFTWARE_SPAdes
	def get_spades_name(self): return self.SOFTWARE_SPAdes_name
	def get_spades_version(self): return self.SOFTWARE_SPAdes_VERSION
	def get_spades_parameters(self): return self.SOFTWARE_SPAdes_PARAMETERS
	def get_spades_parameters_single(self): return self.SOFTWARE_SPAdes_PARAMETERS_single


	"""
	return abyss software
	"""
	def get_abyss(self): return self.SOFTWARE_ABYSS
	def get_abyss_name(self): return self.SOFTWARE_ABYSS_name
	def get_abyss_version(self): return self.SOFTWARE_ABYSS_VERSION
	def get_abyss_parameters(self): return self.SOFTWARE_ABYSS_PARAMETERS

	"""
	return abricate software
	"""
	def get_abricate(self): return self.SOFTWARE_ABRICATE
	def get_abricate_name(self): return self.SOFTWARE_ABRICATE_name
	def get_abricate_version(self): return self.SOFTWARE_ABRICATE_VERSION
	def get_abricate_parameters(self): return self.SOFTWARE_ABRICATE_PARAMETERS

	"""
	return FASTq software
	"""
	def get_fastq(self): return self.SOFTWARE_FASTQ
	def get_fastq_name(self): return self.SOFTWARE_FASTQ_name
	def get_fastq_version(self): return self.SOFTWARE_FASTQ_VERSION
	def get_fastq_parameters(self): return self.SOFTWARE_FASTQ_PARAMETERS
	
	"""
	return trimmomatic software
	"""
	def get_trimmomatic(self): return self.SOFTWARE_TRIMMOMATIC
	def get_trimmomatic_name(self): return self.SOFTWARE_TRIMMOMATIC_name
	def get_trimmomatic_version(self): return self.SOFTWARE_TRIMMOMATIC_VERSION
	def get_trimmomatic_parameters(self): return self.SOFTWARE_TRIMMOMATIC_PARAMETERS
	
	"""
	return snippy software
	"""
	def get_snippy(self): return self.SOFTWARE_SNIPPY
	def get_snippy_name(self): return self.SOFTWARE_SNIPPY_name
	def get_snippy_version(self): return self.SOFTWARE_SNIPPY_VERSION
	def get_snippy_parameters(self): return self.SOFTWARE_SNIPPY_PARAMETERS

	"""
	return snippy-vcf_to_tab software
	"""
	def get_snippy_vcf_to_tab(self): return self.SOFTWARE_SNIPPY_VCF_TO_TAB
	def get_snippy_vcf_to_tab_name(self): return self.SOFTWARE_SNIPPY_VCF_TO_TAB_name
	def get_snippy_vcf_to_tab_version(self): return self.SOFTWARE_SNIPPY_VCF_TO_TAB_VERSION
	def get_snippy_vcf_to_tab_parameters(self): return self.SOFTWARE_SNIPPY_VCF_TO_TAB_PARAMETERS
	
	"""
	return snpEff software
	"""
	def get_snp_eff(self): return self.SOFTWARE_SNP_EFF
	def get_snp_eff_name(self): return self.SOFTWARE_SNP_EFF_name
	def get_snp_eff_config(self): return self.SOFTWARE_SNP_EFF_config
	def get_snp_eff_version(self): return self.SOFTWARE_SNP_EFF_VERSION
	def get_snp_eff_parameters(self): return self.SOFTWARE_SNP_EFF_PARAMETERS
	
	"""
	return genbank2gff3 software
	"""
	def get_genbank2gff3(self): return self.SOFTWARE_GENBANK2GFF3
	def get_genbank2gff3_name(self): return self.SOFTWARE_GENBANK2GFF3_name
	def get_genbank2gff3_version(self): return self.SOFTWARE_GENBANK2GFF3_VERSION
	def get_genbank2gff3_parameters(self): return self.SOFTWARE_GENBANK2GFF3_PARAMETERS
	
	
	"""
	return freebayes software
	"""
	def get_freebayes(self): return self.SOFTWARE_FREEBAYES
	def get_freebayes_parallel(self): return self.SOFTWARE_FREEBAYES_PARALLEL
	def get_freebayes_name(self): return self.SOFTWARE_FREEBAYES_name
	def get_freebayes_version(self): return self.SOFTWARE_FREEBAYES_VERSION
	def get_freebayes_parameters(self): return self.SOFTWARE_FREEBAYES_PARAMETERS

	"""
	get fasta generate regions
	"""
	def get_fasta_generate_regions(self): return self.SOFTWARE_FASTA_GENERATE_REGIONS
	def get_fasta_generate_regions_name(self): return self.SOFTWARE_FASTA_GENERATE_REGIONS_name
	def get_fasta_generate_regions_version(self): return self.SOFTWARE_FASTA_GENERATE_REGIONS_VERSION
	def get_fasta_generate_regions_parameters(self): return self.SOFTWARE_FASTA_GENERATE_REGIONS_PARAMETERS
	
	"""
	get coverge to regions
	"""
	def get_coverage_to_regions(self): return self.SOFTWARE_COVERAGE_TO_REGIONS
	def get_coverage_to_regions_name(self): return self.SOFTWARE_COVERAGE_TO_REGIONS_name
	def get_coverage_to_regions_version(self): return self.SOFTWARE_COVERAGE_TO_REGIONS_VERSION
	def get_coverage_to_regions_parameters(self): return self.SOFTWARE_COVERAGE_TO_REGIONS_PARAMETERS
	
	"""
	get bamtools
	"""
	def get_bamtools(self): return self.SOFTWARE_BAMTOOLS
	def get_bamtools_name(self): return self.SOFTWARE_BAMTOOLS_name
	def get_bamtools_version(self): return self.SOFTWARE_BAMTOOLS_VERSION
	def get_bamtools_parameters(self): return self.SOFTWARE_BAMTOOLS_PARAMETERS
	
	"""
	return bgzip software
	"""
	def get_bgzip(self): return self.SOFTWARE_BGZIP
	def get_bgzip_name(self): return self.SOFTWARE_BGZIP_name
	def get_bgzip_version(self): return self.SOFTWARE_BGZIP_VERSION
	def get_bgzip_parameters(self): return self.SOFTWARE_BGZIP_PARAMETERS

	"""
	return gzip software
	"""
	def get_gzip(self): return self.SOFTWARE_GZIP
	def get_gzip_name(self): return self.SOFTWARE_GZIP_name
	def get_gzip_version(self): return self.SOFTWARE_GZIP_VERSION
	def get_gzip_parameters(self): return self.SOFTWARE_GZIP_PARAMETERS
	
	"""
	return tabix software
	"""
	def get_tabix(self): return self.SOFTWARE_TABIX
	def get_tabix_name(self): return self.SOFTWARE_TABIX_name
	def get_tabix_version(self): return self.SOFTWARE_TABIX_VERSION
	def get_tabix_parameters(self): return self.SOFTWARE_TABIX_PARAMETERS
	
	"""
	return igvtools software
	"""
	def get_igvtools(self): return self.SOFTWARE_IGVTOOLS
	def get_igvtools_name(self): return self.SOFTWARE_IGVTOOLS_name
	def get_igvtools_version(self): return self.SOFTWARE_IGVTOOLS_VERSION
	def get_igvtools_parameters(self): return self.SOFTWARE_IGVTOOLS_PARAMETERS
	
	"""
	return Coverage software
	"""
	def get_coverage(self): return self.SOFTWARE_COVERAGE
	def get_coverage_name(self): return self.SOFTWARE_COVERAGE_name
	def get_coverage_version(self): return self.SOFTWARE_COVERAGEVERSION
	def get_coverage_parameters(self): return self.SOFTWARE_COVERAGE_PARAMETERS
	
	"""
	return Prokka software
	"""
	def get_prokka(self): return self.SOFTWARE_PROKKA
	def get_prokka_name(self): return self.SOFTWARE_PROKKA_name
	def get_prokka_version(self): return self.SOFTWARE_PROKKA_VERSION
	def get_prokka_parameters(self): return self.SOFTWARE_PROKKA_PARAMETERS

	"""
	return Mauve software
	"""
	def get_mauve(self): return self.SOFTWARE_MAUVE
	def get_mauve_name(self): return self.SOFTWARE_MAUVE_name
	def get_mauve_version(self): return self.SOFTWARE_MAUVE_VERSION
	def get_mauve_parameters(self): return self.SOFTWARE_MAUVE_PARAMETERS

	"""
	return Convert Mauve software
	"""
	def get_convert_mauve(self): return self.SOFTWARE_CONVERT
	def get_convert_mauve_name(self): return self.SOFTWARE_CONVERT_name
	def get_convert_mauve_version(self): return self.SOFTWARE_CONVERT_VERSION
	def get_convert_mauve_parameters(self): return self.SOFTWARE_CONVERT_PARAMETERS
	
	"""
	return Mafft software
	"""
	def get_mafft(self): return self.SOFTWARE_MAFFT
	def get_mafft_set_env_variable(self): return self.SOFTWARE_SET_ENV_MAFFT
	def get_mafft_name(self): return self.SOFTWARE_MAFFT_name
	def get_mafft_version(self): return self.SOFTWARE_MAFFT_VERSION
	def get_mafft_parameters(self): return self.SOFTWARE_MAFFT_PARAMETERS
	def get_mafft_parameters_two_parameters(self): return self.SOFTWARE_MAFFT_PARAMETERS_TWO_SEQUENCES
	
	
	"""
	return seqret software EMBOSS
 	"""
	def get_seqret(self): return self.SOFTWARE_SEQRET
	def get_seqret_name(self): return self.SOFTWARE_SEQRET_name
	def get_seqret_version(self): return self.SOFTWARE_SEQRET_VERSION
	def get_seqret_nex_parameters(self): return self.SOFTWARE_SEQRET_NEX_PARAMETERS
	
	"""
	return FastTree software
	"""
	def get_fasttree(self): return self.SOFTWARE_FASTTREE
	def get_fasttree_name(self): return self.SOFTWARE_FASTTREE_name
	def get_fasttree_version(self): return self.SOFTWARE_FASTTREE_VERSION
	def get_fasttree_parameters(self): return self.SOFTWARE_FASTTREE_PARAMETERS
	def get_fasttree_parameters_protein(self): return self.SOFTWARE_FASTTREE_PARAMETERS_PROTEIN


