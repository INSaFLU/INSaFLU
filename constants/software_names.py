"""
Created on Nov 26, 2017

@author: mmp
"""

import os

from django.conf import settings
from settings.constants_settings import ConstantsSettings

from constants.constants import Constants


class SoftwareNames(object):
    """
    classdocs

    """

    ## Simple tags, Yes, No
    SOFTWARE_TAG_yes = "Yes"
    SOFTWARE_TAG_no = "No"

    ## some software is distributed by snippy
    DIR_SOFTWARE_SNIPPY = os.path.join(settings.DIR_SOFTWARE, "snippy")
    SOFTWARE_DEPTH_SAMTOOLS_file_flag = "DEPTH_SAMTOOLS"
    SOFTWARE_SAMTOOLS = os.path.join(DIR_SOFTWARE_SNIPPY, "bin/samtools")
    SOFTWARE_SAMTOOLS_name = "Samtools"
    SOFTWARE_SAMTOOLS_name_depth_ONT = "Samtools depth ONT"
    SOFTWARE_SAMTOOLS_name_extended_depth_ONT = "Coverage calculation (Samtools)"
    SOFTWARE_SAMTOOLS_VERSION = "1.3"
    SOFTWARE_SAMTOOLS_PARAMETERS = ""
    SOFTWARE_SAMTOOLS_parameters_depth_ONT = ""

    SOFTWARE_BGZIP = os.path.join(DIR_SOFTWARE_SNIPPY, "bin/bgzip")
    SOFTWARE_BGZIP_name = "bgzip"
    SOFTWARE_BGZIP_VERSION = "1.3"
    SOFTWARE_BGZIP_PARAMETERS = ""
    SOFTWARE_GZIP = "gzip"
    SOFTWARE_GZIP_name = "gzip"
    SOFTWARE_GZIP_VERSION = "1.6"
    SOFTWARE_GZIP_PARAMETERS = ""
    SOFTWARE_TABIX = os.path.join(DIR_SOFTWARE_SNIPPY, "bin/tabix")
    SOFTWARE_TABIX_name = "tabix"
    SOFTWARE_TABIX_VERSION = "1.3"
    SOFTWARE_TABIX_PARAMETERS = ""
    SOFTWARE_IGVTOOLS = os.path.join(settings.DIR_SOFTWARE, "IGVTools/igvtools.jar")
    SOFTWARE_IGVTOOLS_name = "igvtools"
    SOFTWARE_IGVTOOLS_VERSION = "2.3.98"
    SOFTWARE_IGVTOOLS_PARAMETERS = ""

    SOFTWARE_SPAdes = os.path.join(
        settings.DIR_SOFTWARE, "SPAdes-3.11.1-Linux/bin/spades.py"
    )
    #     SOFTWARE_SPAdes = os.path.join(settings.DIR_SOFTWARE, "SPAdes-3.13.0-Linux/bin/spades.py")
    SOFTWARE_SPAdes_name = "SPAdes"
    SOFTWARE_SPAdes_name_extended = "SPAdes"
    SOFTWARE_SPAdes_VERSION = "3.11.1"  ### older version change at 25/11/2109 to 3.
    #     SOFTWARE_SPAdes_VERSION = "3.13.0"
    SOFTWARE_SPAdes_PARAMETERS = "--only-assembler"
    SOFTWARE_SPAdes_PARAMETERS_single = "--only-assembler"  ## same at this point
    SOFTWARE_SPAdes_CLEAN_HITS_BELLOW_VALUE = (
        3  ## clean the values bellow of this value "NODE_128_length_572_cov_3.682785"
    )
    SOFTWARE_ABRICATE = os.path.join(settings.DIR_SOFTWARE, "abricate/bin/abricate")
    SOFTWARE_ABRICATE_name = "Abricate"
    SOFTWARE_ABRICATE_name_extended = "Abricate"
    SOFTWARE_ABRICATE_DB = os.path.join(settings.DIR_SOFTWARE, "abricate/db")
    SOFTWARE_ABRICATE_VERSION = "0.8-dev"
    SOFTWARE_ABRICATE_PARAMETERS = "--minid 70 --mincov 60"
    SOFTWARE_ABRICATE_PARAMETERS_mincov_30 = "--minid 70 --mincov 30"
    SOFTWARE_FASTQ_name = "FastQC"
    SOFTWARE_FASTQ_VERSION = "0.11.9"
    SOFTWARE_FASTQ_PARAMETERS = ""
    SOFTWARE_FASTQ = os.path.join(
        settings.DIR_SOFTWARE, "FastQC/{}/FastQC/fastqc".format(SOFTWARE_FASTQ_VERSION)
    )
    SOFTWARE_TRIMMOMATIC = os.path.join(
        settings.DIR_SOFTWARE, "trimmomatic/classes/trimmomatic.jar"
    )
    SOFTWARE_TRIMMOMATIC_name = "Trimmomatic"
    SOFTWARE_TRIMMOMATIC_name_extended = "Quality analysis and control (Trimmomatic)"
    SOFTWARE_TRIMMOMATIC_VERSION = "0.39"
    SOFTWARE_TRIMMOMATIC_PARAMETERS = (
        "SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:35 TOPHRED33"
    )
    SOFTWARE_TRIMMOMATIC_vect_info_to_collect = [
        "Input Read Pairs:",
        "Both Surviving:",
        "Forward Only Surviving:",
        "Reverse Only Surviving:",
        "Dropped:",
    ]

    SOFTWARE_TRIMMOMATIC_illuminaclip = "ILLUMINACLIP"
    SOFTWARE_TRIMMOMATIC_addapter_to_replace = "ADAPTER_FILE"
    SOFTWARE_TRIMMOMATIC_addapter_trim_used_to_assemble = ":3:30:10:6:true"
    SOFTWARE_TRIMMOMATIC_addapter_trim = "ILLUMINACLIP:{}:3:30:10:6:true".format(
        SOFTWARE_TRIMMOMATIC_addapter_to_replace
    )
    SOFTWARE_TRIMMOMATIC_addapter_not_apply = "Not apply"
    SOFTWARE_TRIMMOMATIC_addapter_apply_all = "All-adapters.fa"
    SOFTWARE_TRIMMOMATIC_addapter_vect_available = [
        SOFTWARE_TRIMMOMATIC_addapter_not_apply,
        SOFTWARE_TRIMMOMATIC_addapter_apply_all,
        "NexteraPE-PE.fa",
        "TruSeq2-PE.fa",
        "TruSeq2-SE.fa",
        "TruSeq3-PE-2.fa",
        "TruSeq3-PE.fa",
        "TruSeq3-SE.fa",
        "SARS_CoV_2_MN908947_artic_3.fa",
        "SARS_CoV_2_MN908947_artic_4.1.fa",
        "MPXV_MT903345_Yale_PrimalSeq_v.1.fa",
        "MPXV_comb-ccc7sszn.fa",
    ]
    ### collect stat data for ILLUMINA, in form of key value
    SOFTWARE_ILLUMINA_stat = "illumina_stat"
    SOFTWARE_ILLUMINA_stat_collect = [
        "Number of reads R1",
        "Number of reads R2",
        "Average read length R1",
        "Average read length R2",
        "STDEV read length R1",
        "STDEV read length R2",
        "Total reads",
    ]
    SOFTWARE_ILLUMINA_stat_collect_show_percentage = [
        "Number of reads R1",
        "Number of reads R2",
        "Total reads",
    ]

    SOFTWARE_RabbitQC = os.path.join(settings.DIR_SOFTWARE, "RabbitQC/rabbit_qc")
    SOFTWARE_RabbitQC_name = "RabbitQC"
    SOFTWARE_RabbitQC_name_extended = "Quality analysis and control (RabbitQC)"
    SOFTWARE_RabbitQC_VERSION = "0.0.1"
    SOFTWARE_RabbitQC_PARAMETERS = "-w 3 -D"  ## for long reads
    SOFTWARE_NanoStat = "NanoStat"
    SOFTWARE_NanoStat_name = "NanoStat"
    SOFTWARE_NanoStat_name_extended = (
        "Calculate various statistics from a long read sequencing (NanoStat)"
    )
    SOFTWARE_NanoStat_VERSION = "1.4.0"
    SOFTWARE_NanoStat_PARAMETERS = (
        ""  ## -o temp -n temp.txt --fastq ERR4082025_1.fastq.gz
    )
    SOFTWARE_NANOSTAT_vect_info_to_collect = [
        "Mean read length",
        "Mean read quality",
        "Median read length",
        "Median read quality",
        "Number of reads",
        "Read length N50",
        "STDEV read length",
        "Total bases",
    ]
    SOFTWARE_NANOSTAT_vect_info_to_collect_show_percentage = [
        "Mean read length",
        "Median read length",
        "Number of reads",
        "Total bases",
    ]

    SOFTWARE_NanoFilt = "NanoFilt"
    SOFTWARE_NanoFilt_name = "NanoFilt"
    SOFTWARE_NanoFilt_name_extended = (
        "Filtering and trimming of ONT sequencing data (NanoFilt)"
    )
    SOFTWARE_NanoFilt_VERSION = "2.6.0"
    SOFTWARE_NanoFilt_PARAMETERS = "-l 50 -q 10 --headcrop 70 --tailcrop 70"  ## gzip -cd ERR4082025_1.fastq.gz | NanoFilt -q 10 --headcrop 40 --tailcrop 50 | gzip > trimmed-reads.fastq.gz
    SOFTWARE_Medaka_Env = ". {};".format(
        os.path.join(settings.DIR_SOFTWARE, "medaka/bin/activate")
    )
    SOFTWARE_Medaka = "medaka"
    SOFTWARE_Medaka_name = "Medaka"
    SOFTWARE_Medaka_name_consensus = "Medaka consensus"
    SOFTWARE_Medaka_name_variants = "Medaka variants"
    SOFTWARE_Medaka_default_model = "r941_min_high_g360"
    SOFTWARE_Medaka_remove_tags_model = ["_snp_", "_fast_"]
    SOFTWARE_Medaka_name_extended_consensus = "Consensus Generation (Medaka)"
    SOFTWARE_Medaka_name_extended_variant = "Call Variants (Medaka)"
    SOFTWARE_Medaka_PARAMETERS_variant = "--verbose"
    SOFTWARE_Medaka_PARAMETERS_consensus = "-m {}".format(SOFTWARE_Medaka_default_model)
    SOFTWARE_Medaka_VERSION = "1.2.1"

    SOFTWARE_Pangolin_Env = ". {};".format(
        os.path.join(settings.DIR_SOFTWARE, "pangolin/bin/activate")
    )
    SOFTWARE_Pangolin = "pangolin"
    # SOFTWARE_Pangolin_VERSION = "v3.1.14"  ## Version Name: pangolin
    SOFTWARE_Pangolin_VERSION = "v4.0.6"  ## Version Name: pangolin
    SOFTWARE_Pangolin_name = "Pangolin"  ## Pangolin
    SOFTWARE_Pangolin_analysis_mode = "Analysis mode (Pangolin)"  ## Pangolin
    SOFTWARE_Pangolin_name_search_name = (
        "Pango"  ## only to help on the search of Pango Name in output
    )

    #    SOFTWARE_Pangolin_name_search_name = "Pango"        ## only to help on the search of Pango Name in output
    #    SOFTWARE_Pangolin_learn_name_old = "PangolinLearn"         ## was "PangolinLearn", now PangoLearn
    #    SOFTWARE_Pangolin_learn_name = "PangoLearn"         ## was "PangolinLearn", now PangoLearn
    #    SOFTWARE_Pangolin_designation_name = "Pango-designation"        ## Pango Designation
    #    SOFTWARE_Pangolin_constellations_name = "Constellations"    ## Constelations
    #    SOFTWARE_Pangolin_scorpio_name = "scorpio"                    ## Scorpio
    #    SOFTWARE_Pangolin_name_extended = "Phylogenetic Assignment of Named Global Outbreak LINeages (Pangolin)"
    #    SOFTWARE_Pangolin_learn_VERSION = "2021-09-28"        ## Version Name: pangoLearn
    #    SOFTWARE_Pangolin_designation_VERSION = "v1.2.86"    ## Version Name: pango
    #    SOFTWARE_Pangolin_scorpio_VERSION = "v0.3.12"    ## Version Name: pango
    #    SOFTWARE_Pangolin_constellations_VERSION = "v0.0.16"    ## Version Name: pango

    #     VECT_PANGOLIN_TO_TEST = [
    #         SOFTWARE_Pangolin_name,
    #         SOFTWARE_Pangolin_learn_name,
    #         SOFTWARE_Pangolin_designation_name,
    #         SOFTWARE_Pangolin_constellations_name,
    #         SOFTWARE_Pangolin_scorpio_name,
    #         ]

    SOFTWARE_BCFTOOLS = os.path.join(settings.DIR_SOFTWARE, "medaka/bin/bcftools")
    SOFTWARE_BCFTOOLS_name = "bcftools"
    SOFTWARE_BCFTOOLS_VERSION = "1.9"
    SOFTWARE_BCFTOOLS_NEX_PARAMETERS = ""

    ## --nano-raw ../flu_minion/test_minion_seq.fastq.gz --no-alt-contigs -o temp -t 4
    ## ll temp/assembly.fasta
    SOFTWARE_FLYE = os.path.join(settings.DIR_SOFTWARE, "flye/flye.sh")
    SOFTWARE_FLYE_name = "flye"
    SOFTWARE_FLYE_name_extended = "Assembly ONT (Flye)"
    SOFTWARE_FLYE_VERSION = "2.9.1-b1780"
    SOFTWARE_FLYE_PARAMETERS = (
        "--no-alt-contigs "  ## genomeSize=<number>[g|m|k] file1 file2
    )

    ### used to create a file with variations table
    ## git: https://github.com/SantosJGND/INSA
    ## python algn2pheno.py --db DB_COG_UK_antigenic_mutations_2022-05-30.tsv -g S --algn Alignment_aa_SARS_CoV_2_S.fasta -r <referenceName inside Alignment> --odir algn2pheno --output COG_UK_antigenic_mutations
    ####
    ## outputs to provide to the users:
    ## ...._final_report.tsv
    ## ...._flagged_mutation_report.tsv
    ## (if possible, add the Lineage column to these two reports)
    ## ....log
    SOFTWARE_Alignment2phenotype_Env = ". {};".format(
        os.path.join(settings.DIR_SOFTWARE, "Alignment2phenotype/bin/activate")
    )
    SOFTWARE_Alignment2phenotype = "algn2pheno.py"
    SOFTWARE_Alignment2phenotype_VERSION = "v1.0"  ## Version Name: pangolin
    SOFTWARE_Alignment2phenotype_name = "Alignment2phenotype"  ## Alignment2phenotype
    SOFTWARE_Alignment2phenotype_analysis_mode = (
        "Analysis mode (Alignment2phenotype)"  ## Alignment2phenotype
    )
    SOFTWARE_Alignment2phenotype_parameters = "-g S"  ## Alignment2phenotype

    SOFTWARE_SNIPPY = os.path.join(DIR_SOFTWARE_SNIPPY, "bin/snippy")
    SOFTWARE_SNIPPY_name = "Snippy"
    SOFTWARE_SNIPPY_name_extended = "Mapping (Snippy)"
    SOFTWARE_SNIPPY_VERSION = "3.2-dev"
    SOFTWARE_SNIPPY_PARAMETERS = "--mapqual 20 --mincov 10 --minfrac 0.51"

    #### VERY important, change in snippy-vcf
    #     mmp@california:/usr/local/software/insaflu/snippy/bin$ diff snippy-vcf_to_tab_add_freq snippy-vcf_to_tab_add_freq~
    # 57c57
    # < print join("\t", qw(CHROM POS TYPE REF ALT FREQ), @ANNO), "\n";
    # ---
    # > print join("\t", qw(CHROM POS TYPE REF ALT EVIDENCE), @ANNO), "\n";
    SOFTWARE_SNIPPY_VCF_TO_TAB = os.path.join(
        DIR_SOFTWARE_SNIPPY, "bin/snippy-vcf_to_tab_add_freq"
    )
    SOFTWARE_SNIPPY_VCF_TO_TAB_name = "Snippy-vcf_to_tab_add_freq"
    SOFTWARE_SNIPPY_VCF_TO_TAB_VERSION = "3.2-dev"
    SOFTWARE_SNIPPY_VCF_TO_TAB_PARAMETERS = ""
    SOFTWARE_SNIPPY_VCF_TO_TAB_AND_EVIDENCE = os.path.join(
        DIR_SOFTWARE_SNIPPY, "bin/snippy-vcf_to_tab_add_freq_and_evidence"
    )
    SOFTWARE_SNIPPY_VCF_TO_TAB_AND_EVIDENCE_name = (
        "Snippy-vcf_to_tab_add_freq_and_evidence"
    )
    SOFTWARE_SNIPPY_VCF_TO_TAB_AND_EVIDENCE_VERSION = "3.2-dev"
    SOFTWARE_SNIPPY_VCF_TO_TAB_AND_EVIDENCE_PARAMETERS = ""
    SOFTWARE_SNP_EFF = os.path.join(DIR_SOFTWARE_SNIPPY, "bin/snpEff")
    SOFTWARE_SNP_EFF_config = os.path.join(DIR_SOFTWARE_SNIPPY, "etc/snpeff.config")
    SOFTWARE_SNP_EFF_name = "snpEff"
    SOFTWARE_SNP_EFF_VERSION = "4.3p"
    SOFTWARE_SNP_EFF_PARAMETERS = (
        "-no-downstream -no-upstream -no-intergenic -no-utr -noStats"
    )

    SOFTWARE_MSA_MASKER = os.path.join(DIR_SOFTWARE_SNIPPY, "bin/msa_masker.py")
    SOFTWARE_MSA_MASKER_name = "MSA Masker"
    SOFTWARE_MSA_MASKER_VERSION = "1.0"
    ##    SOFTWARE_MSA_MASKER_PARAMETERS = "--g --c"    ### other possibility (--g If the process should mask gaps)
    SOFTWARE_MSA_MASKER_PARAMETERS = "--c"

    SOFTWARE_FREEBAYES = os.path.join(DIR_SOFTWARE_SNIPPY, "bin/freebayes")
    SOFTWARE_FREEBAYES_PARALLEL = os.path.join(
        DIR_SOFTWARE_SNIPPY, "bin/freebayes-parallel"
    )
    SOFTWARE_FREEBAYES_name = "Freebayes"
    SOFTWARE_FREEBAYES_name_extended = "Mapping (Freebayes)"
    SOFTWARE_FREEBAYES_VERSION = "v1.1.0-54-g49413aa"
    SOFTWARE_FREEBAYES_PARAMETERS = "--min-mapping-quality 20 --min-base-quality 20 --min-coverage 100 --min-alternate-count 10 --min-alternate-fraction 0.01 --ploidy 2 -V"

    SOFTWARE_FASTA_GENERATE_REGIONS = os.path.join(
        DIR_SOFTWARE_SNIPPY, "bin/fasta_generate_regions.py"
    )
    SOFTWARE_FASTA_GENERATE_REGIONS_name = "Fasta Generate Regions"
    SOFTWARE_FASTA_GENERATE_REGIONS_VERSION = "1.0"
    SOFTWARE_FASTA_GENERATE_REGIONS_PARAMETERS = ""
    SOFTWARE_COVERAGE_TO_REGIONS = os.path.join(
        settings.DIR_SOFTWARE, "freebayes/scripts/coverage_to_regions.py"
    )
    SOFTWARE_COVERAGE_TO_REGIONS_name = "Coverage to Regions"
    SOFTWARE_COVERAGE_TO_REGIONS_VERSION = "1.0"
    SOFTWARE_COVERAGE_TO_REGIONS_PARAMETERS = ""

    SOFTWARE_BAMTOOLS = os.path.join(
        settings.DIR_SOFTWARE, "bamtools/build/src/toolkit/bamtools"
    )
    SOFTWARE_BAMTOOLS_name = "Bamtools"
    SOFTWARE_BAMTOOLS_VERSION = "2.5"
    SOFTWARE_BAMTOOLS_PARAMETERS = ""

    SOFTWARE_PROKKA = os.path.join(settings.DIR_SOFTWARE, "prokka/bin/prokka")
    SOFTWARE_PROKKA_name = "Prokka"
    SOFTWARE_PROKKA_VERSION = "1.2"
    SOFTWARE_PROKKA_PARAMETERS = (
        "--kingdom Viruses --locustag locus --genus Influenzavirus --species Influenzavirus --strain "
        "ref_PREFIX_FILES_OUT --gcode " + str(Constants.TRANSLATE_TABLE_NUMBER)
    )

    SOFTWARE_MAUVE = os.path.join(settings.DIR_SOFTWARE, "mauve/progressiveMauve")
    SOFTWARE_MAUVE_name = "Mauve"
    SOFTWARE_MAUVE_VERSION = "2.4.0, Feb 13 2015"
    SOFTWARE_MAUVE_PARAMETERS = ""

    SOFTWARE_CONVERT = os.path.join(
        settings.DIR_SOFTWARE, "scripts/convertAlignment.pl"
    )
    SOFTWARE_CONVERT_name = "Convert"
    SOFTWARE_CONVERT_VERSION = ""
    SOFTWARE_CONVERT_PARAMETERS = ""

    vect_versions_available = ["7.453", "7.313"]
    for version in vect_versions_available:
        SOFTWARE_MAFFT_VERSION = version
        SOFTWARE_MAFFT = os.path.join(
            settings.DIR_SOFTWARE,
            "mafft-{}-without-extensions/scripts/mafft".format(version),
        )
        if os.path.exists(SOFTWARE_MAFFT):
            break

    SOFTWARE_SET_ENV_MAFFT = "export MAFFT_BINARIES={}".format(
        os.path.join(
            settings.DIR_SOFTWARE,
            "mafft-{}-without-extensions/binaries".format(version),
        )
    )
    SOFTWARE_MAFFT_name = "Mafft"

    SOFTWARE_MAFFT_PARAMETERS_TWO_SEQUENCES = (
        "--maxiterate 1000 --localpair --preservecase --leavegappyregion --quiet"
    )
    SOFTWARE_MAFFT_PARAMETERS_PROTEIN = "--preservecase --amino --quiet"
    SOFTWARE_MAFFT_PARAMETERS = "--preservecase --quiet"

    SOFTWARE_CLUSTALO = os.path.join(settings.DIR_SOFTWARE, "clustalo/clustalo")
    SOFTWARE_CLUSTALO_name = "clustalo"
    SOFTWARE_CLUSTALO_VERSION = "1.2.4"
    SOFTWARE_CLUSTALO_PARAMETERS = ""

    SOFTWARE_SEQRET = "/usr/bin/seqret"
    SOFTWARE_SEQRET_name = "seqret (EMBOSS)"
    SOFTWARE_SEQRET_VERSION = "6.6.0.0"
    SOFTWARE_SEQRET_NEX_PARAMETERS = "-sformat fasta -osformat2 nexusnon"

    #     SOFTWARE_FASTTREE = os.path.join(DIR_SOFTWARE, "fasttree/FastTree")
    #     SOFTWARE_FASTTREE_name = "FastTree"
    #     SOFTWARE_FASTTREE_VERSION = "2.1.10 SSE3"
    SOFTWARE_FASTTREE = os.path.join(settings.DIR_SOFTWARE, "fasttree/FastTreeDbl")
    SOFTWARE_FASTTREE_name = "FastTreeDbl"
    SOFTWARE_FASTTREE_VERSION = "2.1.10 Double precision"
    SOFTWARE_FASTTREE_PARAMETERS = "-gtr -boot 1000 -nt"
    SOFTWARE_FASTTREE_PARAMETERS_PROTEIN = "-gtr -boot 1000"

    SOFTWARE_FASTQ_TOOLS_SAMPLE = os.path.join(
        settings.DIR_SOFTWARE, "fastq-tools/src/fastq-sample"
    )
    SOFTWARE_FASTQ_TOOLS_SAMPLE_name = "fastq-tools"
    SOFTWARE_FASTQ_TOOLS_SAMPLE_VERSION = "0.8"
    SOFTWARE_FASTQ_TOOLS_SAMPLE_PARAMETERS = ""

    ### not necessary to install, it's a caveat for older versions, before 1/07/2018
    SOFTWARE_CREATE_NEW_REFERENCE_TO_SNIPPY = os.path.join(
        settings.DIR_SOFTWARE, "scripts/create_new_reference_to_snippy.pl"
    )
    SOFTWARE_CREATE_NEW_REFERENCE_TO_SNIPPY_name = "create_new_reference_to_snippy"
    SOFTWARE_CREATE_NEW_REFERENCE_TO_SNIPPY_vesion = "1"
    SOFTWARE_CREATE_NEW_REFERENCE_TO_SNIPPY_parameters = ""

    ### software to produce alignments
    ## nextalign --sequences sequences.fasta --reference SARS_CoV_2_COVID_19_Wuhan_Hu_1_MN908947.fasta --genemap covid.gff3 --genes S --output-dir output
    ## output/sequences.aligned.fasta
    ## output/sequences.gene.S.fasta
    SOFTWARE_NEXTALIGN = os.path.join(
        settings.DIR_SOFTWARE, "nextalign/1.7.0/nextalign"
    )
    SOFTWARE_NEXTALIGN_name = "nextalign"
    SOFTWARE_NEXTALIGN_vesion = "1.7.0"
    SOFTWARE_NEXTALIGN_parameters = ""

    ### genebank to perl
    SOFTWARE_genbank_to_perl = os.path.join(
        settings.DIR_SOFTWARE, "scripts/bp_genbank2gff3.pl"
    )

    SOFTWARE_MASK_CONSENSUS_BY_SITE_name = "Mask consensus by sites"
    SOFTWARE_MASK_CONSENSUS_BY_SITE_name_extended = "Mask consensus by sites"
    SOFTWARE_MASK_CONSENSUS_BY_SITE_vesion = "1"
    SOFTWARE_MASK_CONSENSUS_BY_SITE_parameters = ""  ## don't have anything by default

    ### set consensus on AllConsensus file ON|OFF
    SOFTWARE_GENERATE_CONSENSUS_name = "Generate consensus"
    SOFTWARE_GENERATE_CONSENSUS_name_extended = "Generate consensus"
    SOFTWARE_GENERATE_CONSENSUS_vesion = "1"
    SOFTWARE_GENERATE_CONSENSUS_parameters = ""  ## don't have anything by default

    ### algn2pheno
    SOFTWARE_ALN2PHENO = os.path.join(settings.DIR_SOFTWARE, "aln2pheno/aln2pheno.sh")
    SOFTWARE_ALN2PHENO_README = os.path.join(
        settings.DIR_SOFTWARE, "aln2pheno/algn2pheno_readme.txt"
    )

    ### Nextstrain
    SOFTWARE_NEXTSTRAIN_VERSION = "0.1"
    SOFTWARE_NEXTSTRAIN = os.path.join(
        settings.DIR_SOFTWARE, "nextstrain/nextstrain.sh"
    )
    SOFTWARE_NEXTSTRAIN_MPX = os.path.join(
        settings.DIR_SOFTWARE, "nextstrain/nextstrain_mpx.sh"
    )

    ### Several types of builds for NEXTSTRAIN
    SOFTWARE_NEXTSTRAIN_BUILDS_BASE = os.path.join(
        settings.DIR_SOFTWARE, "nextstrain/nextstrain_builds/"
    )
    SOFTWARE_NEXTSTRAIN_BUILDS_generic = "generic"
    SOFTWARE_NEXTSTRAIN_BUILDS_ncov = "ncov"
    SOFTWARE_NEXTSTRAIN_BUILDS_mpx = "mpx"

    SOFTWARE_NEXTSTRAIN_BUILDS_flu_h3n2_12y = "flu_h3n2_12y"
    SOFTWARE_NEXTSTRAIN_BUILDS_flu_h1n1pdm_12y = "flu_h1n1pdm_12y"
    SOFTWARE_NEXTSTRAIN_BUILDS_flu_vic_12y = "flu_vic_12y"
    SOFTWARE_NEXTSTRAIN_BUILDS_flu_yam_12y = "flu_yam_12y"

    SOFTWARE_NEXTSTRAIN_BUILDS_flu = [
        SOFTWARE_NEXTSTRAIN_BUILDS_flu_h3n2_12y,
        SOFTWARE_NEXTSTRAIN_BUILDS_flu_h1n1pdm_12y,
        SOFTWARE_NEXTSTRAIN_BUILDS_flu_vic_12y,
        SOFTWARE_NEXTSTRAIN_BUILDS_flu_yam_12y,
    ]

    SOFTWARE_NEXTSTRAIN_BUILDS = [
        SOFTWARE_NEXTSTRAIN_BUILDS_generic,
        SOFTWARE_NEXTSTRAIN_BUILDS_ncov,
        SOFTWARE_NEXTSTRAIN_BUILDS_flu_h3n2_12y,
        SOFTWARE_NEXTSTRAIN_BUILDS_flu_h1n1pdm_12y,
        SOFTWARE_NEXTSTRAIN_BUILDS_flu_vic_12y,
        SOFTWARE_NEXTSTRAIN_BUILDS_flu_yam_12y,
        SOFTWARE_NEXTSTRAIN_BUILDS_mpx,
    ]
    SOFTWARE_NEXTSTRAIN_BUILDS_DESC = [
        [SOFTWARE_NEXTSTRAIN_BUILDS_generic, "Generic"],
        [SOFTWARE_NEXTSTRAIN_BUILDS_ncov, "SARS-CoV-2"],
        [SOFTWARE_NEXTSTRAIN_BUILDS_flu_h3n2_12y, "Influenza (H3N2 HA 12years)"],
        [SOFTWARE_NEXTSTRAIN_BUILDS_flu_h1n1pdm_12y, "Influenza (H1N1PDM HA 12years)"],
        [SOFTWARE_NEXTSTRAIN_BUILDS_flu_vic_12y, "Influenza (Vic HA 12years)"],
        [SOFTWARE_NEXTSTRAIN_BUILDS_flu_yam_12y, "Influenza (Yam HA 12years)"],
        [SOFTWARE_NEXTSTRAIN_BUILDS_mpx, "Monkeypox"],
    ]

    # default build
    SOFTWARE_NEXTSTRAIN_BUILDS_parameter = SOFTWARE_NEXTSTRAIN_BUILDS_generic
    SOFTWARE_NEXTSTRAIN_name = "Nextstrain"
    SOFTWARE_NEXTSTRAIN_name_extended = "Nextstrain pathogen evolution"

    ### PATHOGEN IDENTIFICATION SOFTWARE

    PATHOGEN_IDENTIFICATION_name = "Pathogen Identification"
    PATHOGEN_IDENTIFICATION_deplete = True
    PATHOGEN_IDENTIFICATION_enrich = True
    PATHOGEN_IDENTIFICATION_deplete_skip = False
    PATHOGEN_IDENTIFICATION_enrich_skip = False
    PATHOGEN_IDENTIFICATION_assemble = True
    PATHOGEN_IDENTIFICATION_remap = True

    ### CENTRIFUGE
    SOFTWARE_CENTRIFUGE = os.path.join(
        settings.DIR_SOFTWARE,
        "host_depletion/hostdep_env/bin/centrifuge",
    )
    SOFTWARE_CENTRIFUGE_name = "Centrifuge"
    SOFTWARE_CENTRIFUGE_name_extended = "Centrifuge"

    SOFTWARE_CENTRIFUGE_VERSION = "1.0.4_beta"
    SOFTWARE_CENTRIFUGE_parameters = "-p 4 --time -k 3 --out-fmt sam"
    SOFTWARE_CENTRIFUGE_PARAMETERS_TWO_SEQUENCES = "-p 4 --time -k 3 --out-fmt sam"

    ### BWA
    SOFTWARE_BWA = os.path.join(
        settings.DIR_SOFTWARE,
        "preprocess/preproc/bin/bwa",
    )

    SOFTWARE_BWA_name = "BWA"
    SOFTWARE_BWA_name_extended = "BWA"
    SOFTWARE_BWA_VERSION = "0.7.17"
    SOFTWARE_BWA_parameters = "mem -t 4 -M -R '@RG\tID:foo\tSM:bar\tPL:illumina'"
    SOFTWARE_BWA_PARAMETERS_TWO_SEQUENCES = (
        "mem -t 4 -M -R '@RG\tID:foo\tSM:bar\tPL:illumina'"
    )

    ### KRAKEN2
    SOFTWARE_KRAKEN2 = os.path.join(
        settings.DIR_SOFTWARE, "kraken2/kraken_env/bin/kraken2"
    )
    SOFTWARE_KRAKEN2_name = "Kraken2"
    SOFTWARE_KRAKEN2_name_extended = "Kraken2"
    SOFTWARE_KRAKEN2_VERSION = "2.1.2"
    SOFTWARE_KRAKEN2_parameters = "--threads 4 --gzip-compressed  --confidence .5"
    SOFTWARE_KRAKEN2_PARAMETERS_TWO_SEQUENCES = (
        "--threads 4 --gzip-compressed  --confidence .5"
    )
    ### KRAKENUNIQ
    SOFTARE_KRAKENUNIQ = os.path.join(
        settings.DIR_SOFTWARE,
        "host_depletion/hostdep_env/bin/krakenuniq",
    )
    SOFTWARE_KRAKENUNIQ_name = "Krakenuniq"
    SOFTWARE_KRAKENUNIQ_name_extended = "Krakenuniq"
    SOFTWARE_KRAKENUNIQ_VERSION = "0.7.3"
    SOFTWARE_KRAKENUNIQ_parameters = (
        "--threads 4 --gzip-compressed --fastq-input --hll-precision 12"
    )
    SOFTWARE_KRAKENUNIQ_PARAMETERS_TWO_SEQUENCES = (
        "--threads 4 --gzip-compressed --fastq-input --hll-precision 12"
    )

    ### FAST VIROME EXPLORER

    SOFTWARE_FASTVIROMEEXPLORER_name = "FastViromeExplorer"
    SOFTWARE_FASTVIROMEEXPLORER_name_extended = "FastViromeExplorer"
    SOFTWARE_FASTVIROMEEXPLORER = os.path.join(
        settings.DIR_SOFTWARE,
        "FastViromeExplorer/FastViromeExplorer/bin/FastViromeExplorer",
    )
    SOFTWARE_FASTVIROMEEXPLORER_VERSION = "1.3"
    SOFTWARE_FASTVIROMEEXPLORER_parameters = ""

    ### DIAMOND
    SOFTWARE_DIAMOND_name = "Diamond"
    SOFTWARE_DIAMOND_name_extended = "Diamond"
    SOFTWARE_DIAMOND = os.path.join(
        settings.DIR_SOFTWARE,
        "hostDepletion/hostdep_env/bin/diamond",
    )
    SOFTWARE_DIAMOND_VERSION = "2.0.15"
    SOFTWARE_DIAMOND_parameters = (
        "-p 5 --top 5 -e 0.01 --id 65 --query-cover 50 --sensitive"
    )

    ### BLASTN
    SOFTWARE_BLAST_name = "Blastn"
    SOFTWARE_BLAST_name_extended = "Blastn"
    SOFTWARE_BLAST = os.path.join(
        settings.DIR_SOFTWARE, "hostDepletion/hostdep_env/bin/blast"
    )
    SOFTWARE_BLAST_VERSION = "2.12.0"
    SOFTWARE_BLAST_parameters = "-num_threads 5 -outfmt 6"

    ### KAIJU
    SOFTWARE_KAIJU = os.path.join(
        settings.DIR_SOFTWARE,
        "hostDepletion/hostdep_env/bin/kaiju",
    )
    SOFTWARE_KAIJU_name = "Kaiju"
    SOFTWARE_KAIJU_name_extended = "Kaiju"
    SOFTWARE_KAIJU_VERSION = "1.8.2"
    SOFTWARE_KAIJU_parameters = (
        "-t 4 -z -f kaiju_db/nodes.dmp -i kaiju_db/kaiju_db.fmi -o -v -x -u -v"
    )

    ### DESAMBA

    SOFTWARE_DESAMBA = os.path.join(
        settings.DIR_SOFTWARE,
        "classm_lc/deSAMBA/bin/deSAMBA",
    )
    SOFTWARE_DESAMBA_name = "Desamba"
    SOFTWARE_DESAMBA_name_extended = "deSAMBA"
    SOFTWARE_DESAMBA_VERSION = "1.0.0"
    SOFTWARE_DESAMBA_parameters = ""
    ### RAVEN

    SOFTWARE_RAVEN_name = "Raven"
    SOFTWARE_RAVEN_name_extended = "Raven"
    SOFTWARE_RAVEN = os.path.join(settings.DIR_SOFTWARE, "assembly/assembly/bin/raven")
    SOFTWARE_RAVEN_VERSION = "1.8.1"
    SOFTWARE_RAVEN_parameters = "--threads 4 --gzip-compressed -p2"

    ### SNIPPY

    SOFTWARE_SNIPPY_PI_name = "Snippy_PI"
    SOFTWARE_SNIPPY_PI_name_extended = "Snippy"
    SOFTWARE_SNIPPY_PI = os.path.join(
        settings.DIR_SOFTWARE,
        "preprocess/preproc/bin/snippy",
    )
    SOFTWARE_SNIPPY_PI_VERSION = "3.2_dev"
    SOFTWARE_SNIPPY_PI_parameters = "--cpus 4 mapqual=40 --mincov 10 --minfrac 0.1"

    ### MINIMAP2 REMAP

    SOFTWARE_MINIMAP2_REMAP_ONT_name = "Minimap2"
    SOFTWARE_MINIMAP2_REMAP_ONT_name_extended = "Minimap2"
    SOFTWARE_MINIMAP2_REMAP_ONT = os.path.join(
        settings.DIR_SOFTWARE,
        "preprocess/preproc/bin/minimap2",
    )
    SOFTWARE_MINIMAP2_REMAP_ONT_VERSION = "2.24"
    SOFTWARE_MINIMAP2_REMAP_ONT_parameters = "-a -x map-ont -t 4"

    ### MINIMAP2 DEPLETE

    SOFTWARE_MINIMAP2_DEPLETE_ONT_name = "Minimap2_ONT"
    SOFTWARE_MINIMAP2_DEPLETE_ONT_name_extended = "Minimap2"
    SOFTWARE_MINIMAP2_DEPLETE_ONT = os.path.join(
        settings.DIR_SOFTWARE,
        "hostDepletion/hostdep_env/bin/minimap2",
    )
    SOFTWARE_MINIMAP2_DEPLETE_ONT_VERSION = "2.24"
    SOFTWARE_MINIMAP2_DEPLETE_ONT_parameters = "-a -x map-ont -t 4"

    ### Bowtie2 DEPLETE

    SOFTWARE_BOWTIE2_DEPLETE_name = "Bowtie2"
    SOFTWARE_BOWTIE2_DEPLETE_name_extended = "Bowtie2"
    SOFTWARE_BOWTIE2_DEPLETE = os.path.join(
        settings.DIR_SOFTWARE,
        "remap/remap/bin/bowtie2",
    )
    SOFTWARE_BOWTIE2_DEPLETE_VERSION = "2.4.5"
    SOFTWARE_BOWTIE2_DEPLETE_parameters = "-p 4"

    ###################################
    ###################################
    #####
    #####    Global parameters for INSaFLU
    #####

    INSAFLU_PARAMETER_MASK_CONSENSUS_name = "Threshold to mask consensus sequences"
    INSAFLU_PARAMETER_MASK_CONSENSUS_name_extended = (
        "Minimum percentage of horizontal coverage to generate consensus"
    )
    INSAFLU_PARAMETER_MASK_CONSENSUS_vesion = "1"
    INSAFLU_PARAMETER_MASK_CONSENSUS_parameters = "Threshold:70"

    INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name = (
        "Minimum depth of coverage per site to validate the sequence"
    )
    INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name_extended = (
        "Minimum depth of coverage per site to validate the sequence (mincov)"
    )
    INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_vesion = "1"
    INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_parameters = "Threshold:30"

    ### It also put Ns on positions where variations are below this threshold (masking consensus on ONT)
    INSAFLU_PARAMETER_VCF_FREQ_ONT_name = "Minimum proportion for variant evidence"
    INSAFLU_PARAMETER_VCF_FREQ_ONT_name_extended = (
        "Minimum proportion for variant evidence (minfrac)"
    )
    INSAFLU_PARAMETER_VCF_FREQ_ONT_vesion = "1"
    INSAFLU_PARAMETER_VCF_FREQ_ONT_parameters = "Threshold:0.51"

    ### clean human reads
    SOFTWARE_CLEAN_HUMAN_READS_name = "Clean human reads"
    SOFTWARE_CLEAN_HUMAN_READS_extended = "Remove human reads on fastq files"
    SOFTWARE_CLEAN_HUMAN_READS_VERSION = "1"
    SOFTWARE_CLEAN_HUMAN_READS_PARAMETERS = SOFTWARE_TAG_no
    SOFTWARE_CLEAN_HUMAN_READS_vect_available = [
        SOFTWARE_TAG_yes,
        SOFTWARE_TAG_no,
    ]

    ### Extend Parameters names in output tables, "csv" and "tst"
    VECT_INSAFLU_PARAMETER = [
        INSAFLU_PARAMETER_MASK_CONSENSUS_name,
        #                 SOFTWARE_CLEAN_HUMAN_READS_name,
    ]
    #         INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name]

    #####
    #####    END Global parameters for INSaFLU
    #####
    ###################################

    ###################################
    ###################################
    ###
    ###        Relation between Technology -> Software
    ###
    DICT_SOFTWARE_RELATION = {
        ConstantsSettings.TECHNOLOGY_illumina: [
            SOFTWARE_TRIMMOMATIC_name,
            SOFTWARE_SNIPPY_name,
            SOFTWARE_FREEBAYES_name,
            INSAFLU_PARAMETER_MASK_CONSENSUS_name,
            SOFTWARE_MASK_CONSENSUS_BY_SITE_name,
            SOFTWARE_CLEAN_HUMAN_READS_name,
            SOFTWARE_ABRICATE_name,
            SOFTWARE_CENTRIFUGE_name,
            SOFTWARE_KRAKEN2_name,
            SOFTWARE_KRAKENUNIQ_name,
            SOFTWARE_KAIJU_name,
        ],
        ConstantsSettings.TECHNOLOGY_minion: [
            SOFTWARE_NanoFilt_name,
            SOFTWARE_Medaka_name,
            INSAFLU_PARAMETER_MASK_CONSENSUS_name,
            SOFTWARE_CLEAN_HUMAN_READS_name,
            INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name,
            SOFTWARE_MASK_CONSENSUS_BY_SITE_name,
            INSAFLU_PARAMETER_VCF_FREQ_ONT_name,
            SOFTWARE_ABRICATE_name,
            SOFTWARE_DIAMOND_name,
            SOFTWARE_CENTRIFUGE_name,
            SOFTWARE_KRAKEN2_name,
            SOFTWARE_KRAKENUNIQ_name,
            SOFTWARE_KAIJU_name,
        ],
    }
    ###
    ###################################

    ### software with application in multiple pipeline_steps:
    polyvalent_software = [
        SOFTWARE_CENTRIFUGE_name,
        SOFTWARE_SNIPPY_name,
    ]
    # pipeline_steps per software, for software with multiple pipeline_steps.
    polyvalent_software_pipelines = {
        SOFTWARE_CENTRIFUGE_name: [
            ConstantsSettings.PIPELINE_NAME_viral_enrichment,
            ConstantsSettings.PIPELINE_NAME_read_classification,
        ],
        SOFTWARE_SNIPPY_name: [
            ConstantsSettings.PIPELINE_NAME_variant_detection,
            ConstantsSettings.PIPELINE_NAME_remapping,
        ],
        SOFTWARE_MINIMAP2_REMAP_ONT_name: [
            ConstantsSettings.PIPELINE_NAME_remapping,
            ConstantsSettings.PIPELINE_NAME_host_depletion,
        ],
    }

    ###################################

    def __init__(self):
        """
        Constructor
        """
        pass

    """
    return samtools software
    """

    def get_samtools(self):
        return self.SOFTWARE_SAMTOOLS

    def get_samtools_name(self):
        return self.SOFTWARE_SAMTOOLS_name

    def get_samtools_name_depth_ONT(self):
        return self.SOFTWARE_SAMTOOLS_name_depth_ONT

    def get_samtools_version(self):
        return self.SOFTWARE_SAMTOOLS_VERSION

    def get_samtools_parameters_depth_ONT(self):
        return self.SOFTWARE_SAMTOOLS_parameters_depth_ONT

    """
    return spades software
    """

    def get_spades(self):
        return self.SOFTWARE_SPAdes

    def get_spades_name(self):
        return self.SOFTWARE_SPAdes_name

    def get_spades_version(self):
        return self.SOFTWARE_SPAdes_VERSION

    def get_spades_parameters(self):
        return self.SOFTWARE_SPAdes_PARAMETERS

    def get_spades_parameters_single(self):
        return self.SOFTWARE_SPAdes_PARAMETERS_single

    """
    return abricate software
    """

    def get_abricate(self):
        return self.SOFTWARE_ABRICATE

    def get_abricate_name(self):
        return self.SOFTWARE_ABRICATE_name

    def get_abricate_version(self):
        return self.SOFTWARE_ABRICATE_VERSION

    def get_abricate_parameters(self):
        return self.SOFTWARE_ABRICATE_PARAMETERS

    def get_abricate_parameters_mincov_30(self):
        return self.SOFTWARE_ABRICATE_PARAMETERS_mincov_30

    """
    return FASTq software
    """

    def get_fastq(self):
        return self.SOFTWARE_FASTQ

    def get_fastq_name(self):
        return self.SOFTWARE_FASTQ_name

    def get_fastq_version(self):
        return self.SOFTWARE_FASTQ_VERSION

    def get_fastq_parameters(self):
        return self.SOFTWARE_FASTQ_PARAMETERS

    """
    return trimmomatic software
    """

    def get_trimmomatic(self):
        return self.SOFTWARE_TRIMMOMATIC

    def get_trimmomatic_name(self):
        return self.SOFTWARE_TRIMMOMATIC_name

    def get_trimmomatic_name_extended(self):
        return self.SOFTWARE_TRIMMOMATIC_name_extended

    def get_trimmomatic_version(self):
        return self.SOFTWARE_TRIMMOMATIC_VERSION

    def get_trimmomatic_parameters(self):
        return self.SOFTWARE_TRIMMOMATIC_PARAMETERS

    """
    return msa masker software
    """

    def get_msa_masker(self):
        return self.SOFTWARE_MSA_MASKER

    def get_msa_masker_name(self):
        return self.SOFTWARE_MSA_MASKER_name

    def get_msa_masker_version(self):
        return self.SOFTWARE_MSA_MASKER_VERSION

    def get_msa_masker_parameters(self):
        return self.SOFTWARE_MSA_MASKER_PARAMETERS

    """
    return snippy software
    """

    def get_snippy(self):
        return self.SOFTWARE_SNIPPY

    def get_snippy_name(self):
        return self.SOFTWARE_SNIPPY_name

    def get_snippy_name_extended(self):
        return self.SOFTWARE_SNIPPY_name_extended

    def get_snippy_version(self):
        return self.SOFTWARE_SNIPPY_VERSION

    def get_snippy_parameters(self):
        return self.SOFTWARE_SNIPPY_PARAMETERS

    """
    return snippy-vcf_to_tab software. Add FRED
    """

    def get_snippy_vcf_to_tab(self):
        return self.SOFTWARE_SNIPPY_VCF_TO_TAB

    def get_snippy_vcf_to_tab_name(self):
        return self.SOFTWARE_SNIPPY_VCF_TO_TAB_name

    def get_snippy_vcf_to_tab_version(self):
        return self.SOFTWARE_SNIPPY_VCF_TO_TAB_VERSION

    def get_snippy_vcf_to_tab_parameters(self):
        return self.SOFTWARE_SNIPPY_VCF_TO_TAB_PARAMETERS

    """
    return snippy-vcf_to_tab software. Add FRED and Evidence
    """

    def get_snippy_vcf_to_tab_freq_and_evidence(self):
        return self.SOFTWARE_SNIPPY_VCF_TO_TAB_AND_EVIDENCE

    def get_snippy_vcf_to_tab_freq_and_evidence_name(self):
        return self.SOFTWARE_SNIPPY_VCF_TO_TAB_AND_EVIDENCE_name

    def get_snippy_vcf_to_tab_freq_and_evidence_version(self):
        return self.SOFTWARE_SNIPPY_VCF_TO_TAB_AND_EVIDENCE_VERSION

    def get_snippy_vcf_to_tab_freq_and_evidence_parameters(self):
        return self.SOFTWARE_SNIPPY_VCF_TO_TAB_AND_EVIDENCE_PARAMETERS

    """
    return snpEff software
    """

    def get_snp_eff(self):
        return self.SOFTWARE_SNP_EFF

    def get_snp_eff_name(self):
        return self.SOFTWARE_SNP_EFF_name

    def get_snp_eff_config(self):
        return self.SOFTWARE_SNP_EFF_config

    def get_snp_eff_version(self):
        return self.SOFTWARE_SNP_EFF_VERSION

    def get_snp_eff_parameters(self):
        return self.SOFTWARE_SNP_EFF_PARAMETERS

    """
    return freebayes software
    """

    def get_freebayes(self):
        return self.SOFTWARE_FREEBAYES

    def get_freebayes_parallel(self):
        return self.SOFTWARE_FREEBAYES_PARALLEL

    def get_freebayes_name(self):
        return self.SOFTWARE_FREEBAYES_name

    def get_freebayes_version(self):
        return self.SOFTWARE_FREEBAYES_VERSION

    def get_freebayes_parameters(self):
        return self.SOFTWARE_FREEBAYES_PARAMETERS

    """
    get fasta generate regions
    """

    def get_fasta_generate_regions(self):
        return self.SOFTWARE_FASTA_GENERATE_REGIONS

    def get_fasta_generate_regions_name(self):
        return self.SOFTWARE_FASTA_GENERATE_REGIONS_name

    def get_fasta_generate_regions_version(self):
        return self.SOFTWARE_FASTA_GENERATE_REGIONS_VERSION

    def get_fasta_generate_regions_parameters(self):
        return self.SOFTWARE_FASTA_GENERATE_REGIONS_PARAMETERS

    """
    get coverge to regions
    """

    def get_coverage_to_regions(self):
        return self.SOFTWARE_COVERAGE_TO_REGIONS

    def get_coverage_to_regions_name(self):
        return self.SOFTWARE_COVERAGE_TO_REGIONS_name

    def get_coverage_to_regions_version(self):
        return self.SOFTWARE_COVERAGE_TO_REGIONS_VERSION

    def get_coverage_to_regions_parameters(self):
        return self.SOFTWARE_COVERAGE_TO_REGIONS_PARAMETERS

    """
    get bamtools
    """

    def get_bamtools(self):
        return self.SOFTWARE_BAMTOOLS

    def get_bamtools_name(self):
        return self.SOFTWARE_BAMTOOLS_name

    def get_bamtools_version(self):
        return self.SOFTWARE_BAMTOOLS_VERSION

    def get_bamtools_parameters(self):
        return self.SOFTWARE_BAMTOOLS_PARAMETERS

    """
    return bgzip software
    """

    def get_bgzip(self):
        return self.SOFTWARE_BGZIP

    def get_bgzip_name(self):
        return self.SOFTWARE_BGZIP_name

    def get_bgzip_version(self):
        return self.SOFTWARE_BGZIP_VERSION

    def get_bgzip_parameters(self):
        return self.SOFTWARE_BGZIP_PARAMETERS

    """
    return gzip software
    """

    def get_gzip(self):
        return self.SOFTWARE_GZIP

    def get_gzip_name(self):
        return self.SOFTWARE_GZIP_name

    def get_gzip_version(self):
        return self.SOFTWARE_GZIP_VERSION

    def get_gzip_parameters(self):
        return self.SOFTWARE_GZIP_PARAMETERS

    """
    return tabix software
    """

    def get_tabix(self):
        return self.SOFTWARE_TABIX

    def get_tabix_name(self):
        return self.SOFTWARE_TABIX_name

    def get_tabix_version(self):
        return self.SOFTWARE_TABIX_VERSION

    def get_tabix_parameters(self):
        return self.SOFTWARE_TABIX_PARAMETERS

    """
    return igvtools software
    """

    def get_igvtools(self):
        return self.SOFTWARE_IGVTOOLS

    def get_igvtools_name(self):
        return self.SOFTWARE_IGVTOOLS_name

    def get_igvtools_version(self):
        return self.SOFTWARE_IGVTOOLS_VERSION

    def get_igvtools_parameters(self):
        return self.SOFTWARE_IGVTOOLS_PARAMETERS

    """
    return Coverage software
    """
    #     def get_coverage(self): return self.SOFTWARE_COVERAGE
    #     def get_coverage_name(self): return self.SOFTWARE_COVERAGE_name
    #     def get_coverage_version(self): return self.SOFTWARE_COVERAGEVERSION
    #     def get_coverage_parameters(self): return self.SOFTWARE_COVERAGE_PARAMETERS

    """
    return Prokka software
    """

    def get_prokka(self):
        return self.SOFTWARE_PROKKA

    def get_prokka_name(self):
        return self.SOFTWARE_PROKKA_name

    def get_prokka_version(self):
        return self.SOFTWARE_PROKKA_VERSION

    def get_prokka_parameters(self):
        return self.SOFTWARE_PROKKA_PARAMETERS

    """
    return Mauve software
    """

    def get_mauve(self):
        return self.SOFTWARE_MAUVE

    def get_mauve_name(self):
        return self.SOFTWARE_MAUVE_name

    def get_mauve_version(self):
        return self.SOFTWARE_MAUVE_VERSION

    def get_mauve_parameters(self):
        return self.SOFTWARE_MAUVE_PARAMETERS

    """
    return Convert Mauve software
    """

    def get_convert_mauve(self):
        return self.SOFTWARE_CONVERT

    def get_convert_mauve_name(self):
        return self.SOFTWARE_CONVERT_name

    def get_convert_mauve_version(self):
        return self.SOFTWARE_CONVERT_VERSION

    def get_convert_mauve_parameters(self):
        return self.SOFTWARE_CONVERT_PARAMETERS

    """
    return Mafft software
    """

    def get_mafft(self):
        return self.SOFTWARE_MAFFT

    def get_mafft_set_env_variable(self):
        return self.SOFTWARE_SET_ENV_MAFFT

    def get_mafft_name(self):
        return self.SOFTWARE_MAFFT_name

    def get_mafft_version(self):
        return self.SOFTWARE_MAFFT_VERSION

    def get_mafft_parameters(self):
        return self.SOFTWARE_MAFFT_PARAMETERS

    def get_mafft_parameters_two_parameters(self):
        return self.SOFTWARE_MAFFT_PARAMETERS_TWO_SEQUENCES

    """
    return ClustalO software
    """

    def get_clustalo(self):
        return self.SOFTWARE_CLUSTALO

    def get_clustalo_name(self):
        return self.SOFTWARE_CLUSTALO_name

    def get_clustalo_version(self):
        return self.SOFTWARE_CLUSTALO_VERSION

    def get_clustalo_parameters(self):
        return self.SOFTWARE_CLUSTALO_PARAMETERS

    """
    return seqret software EMBOSS
    """

    def get_seqret(self):
        return self.SOFTWARE_SEQRET

    def get_seqret_name(self):
        return self.SOFTWARE_SEQRET_name

    def get_seqret_version(self):
        return self.SOFTWARE_SEQRET_VERSION

    def get_seqret_nex_parameters(self):
        return self.SOFTWARE_SEQRET_NEX_PARAMETERS

    """
    return FastTree software
    """

    def get_fasttree(self):
        return self.SOFTWARE_FASTTREE

    def get_fasttree_name(self):
        return self.SOFTWARE_FASTTREE_name

    def get_fasttree_version(self):
        return self.SOFTWARE_FASTTREE_VERSION

    def get_fasttree_parameters(self):
        return self.SOFTWARE_FASTTREE_PARAMETERS

    def get_fasttree_parameters_protein(self):
        return self.SOFTWARE_FASTTREE_PARAMETERS_PROTEIN

    """
    return FastqTools sample software
    """

    def get_fastqtools_sample(self):
        return self.SOFTWARE_FASTQ_TOOLS_SAMPLE

    def get_fastqtools_sample_name(self):
        return self.SOFTWARE_FASTQ_TOOLS_SAMPLE_name

    def get_fastqtools_sample_version(self):
        return self.SOFTWARE_FASTQ_TOOLS_SAMPLE_VERSION

    def get_fastqtools_sample_parameters(self):
        return self.SOFTWARE_FASTQ_TOOLS_SAMPLE_PARAMETERS

    """
    return BCFTOOLS sample software
    """

    def get_bcftools(self):
        return self.SOFTWARE_BCFTOOLS

    def get_bcftools_name(self):
        return self.SOFTWARE_BCFTOOLS_name

    def get_bcftools_version(self):
        return self.SOFTWARE_BCFTOOLS_VERSION

    def get_bcftools_parameters(self):
        return self.SOFTWARE_BCFTOOLS_PARAMETERS

    """
    return 
    ### not necessary to install, it's a caveat for older versions, before 1/07/2018
    """

    def get_create_new_reference_to_snippy(self):
        return self.SOFTWARE_CREATE_NEW_REFERENCE_TO_SNIPPY

    def get_create_new_reference_to_snippy_name(self):
        return self.SOFTWARE_CREATE_NEW_REFERENCE_TO_SNIPPY_name

    def get_create_new_reference_to_snippy_version(self):
        return self.SOFTWARE_CREATE_NEW_REFERENCE_TO_SNIPPY_vesion

    def get_create_new_reference_to_snippy_parameters(self):
        return self.SOFTWARE_CREATE_NEW_REFERENCE_TO_SNIPPY_parameters

    """
    return mask consensus by sites
    """

    def get_software_mask_consensus_by_site(self):
        return self.SOFTWARE_MASK_CONSENSUS_BY_SITE_name

    def get_software_mask_consensus_by_site_name(self):
        return self.SOFTWARE_MASK_CONSENSUS_BY_SITE_name

    def get_software_mask_consensus_by_site_name_extended(self):
        return self.SOFTWARE_MASK_CONSENSUS_BY_SITE_name_extended

    def get_software_mask_consensus_by_site_version(self):
        return self.SOFTWARE_MASK_CONSENSUS_BY_SITE_vesion

    def get_software_mask_consensus_by_site_parameters(self):
        return self.SOFTWARE_MASK_CONSENSUS_BY_SITE_parameters

    """
    return generate consensus
    """

    def get_software_generate_consensus_(self):
        return self.SOFTWARE_GENERATE_CONSENSUS_name

    def get_software_generate_consensus_name(self):
        return self.SOFTWARE_GENERATE_CONSENSUS_name

    def get_software_generate_consensus_name_extended(self):
        return self.SOFTWARE_GENERATE_CONSENSUS_name_extended

    def get_software_generate_consensus_version(self):
        return self.SOFTWARE_GENERATE_CONSENSUS_vesion

    def get_software_generate_consensus_parameters(self):
        return self.SOFTWARE_GENERATE_CONSENSUS_parameters

    """
    return 
    ### nextalign
    """

    def get_nextalign(self):
        return self.SOFTWARE_NEXTALIGN

    def get_nextalign_name(self):
        return self.SOFTWARE_NEXTALIGN_name

    def get_nextalign_version(self):
        return self.SOFTWARE_NEXTALIGN_vesion

    def get_nextalign_parameters(self):
        return self.SOFTWARE_NEXTALIGN_parameters

    """
    Used in Alignment2phenotype
    """

    def get_Alignment2phenotype_env(self):
        return self.SOFTWARE_Alignment2phenotype_Env

    def get_Alignment2phenotype(self):
        return self.SOFTWARE_Alignment2phenotype

    def get_Alignment2phenotype_name(self):
        return self.SOFTWARE_Alignment2phenotype_name

    def get_Alignment2phenotype_name_extended(self):
        return self.SOFTWARE_Alignment2phenotype_name_extended

    def get_Alignment2phenotype_version(self):
        return self.SOFTWARE_Alignment2phenotype_VERSION

    def get_Alignment2phenotype_parameters(self):
        return self.SOFTWARE_Alignment2phenotype_PARAMETERS

    ###### START minion software
    ###
    def get_rabbitQC(self):
        return self.SOFTWARE_RabbitQC

    def get_rabbitQC_name(self):
        return self.SOFTWARE_RabbitQC_name

    def get_rabbitQC_name_extended(self):
        return self.SOFTWARE_RabbitQC_name_extended

    def get_rabbitQC_version(self):
        return self.SOFTWARE_RabbitQC_VERSION

    def get_rabbitQC_parameters(self):
        return self.SOFTWARE_RabbitQC_PARAMETERS

    def get_NanoStat(self):
        return self.SOFTWARE_NanoStat

    def get_NanoStat_name(self):
        return self.SOFTWARE_NanoStat_name

    def get_NanoStat_name_extended(self):
        return self.SOFTWARE_NanoStat_name_extended

    def get_NanoStat_version(self):
        return self.SOFTWARE_NanoStat_VERSION

    def get_NanoStat_parameters(self):
        return self.SOFTWARE_NanoStat_PARAMETERS

    def get_NanoFilt(self):
        return self.SOFTWARE_NanoFilt

    def get_NanoFilt_name(self):
        return self.SOFTWARE_NanoFilt_name

    def get_NanoFilt_name_extended(self):
        return self.SOFTWARE_NanoFilt_name_extended

    def get_NanoFilt_version(self):
        return self.SOFTWARE_NanoFilt_VERSION

    def get_NanoFilt_parameters(self):
        return self.SOFTWARE_NanoFilt_PARAMETERS

    def get_medaka_env(self):
        return self.SOFTWARE_Medaka_Env

    def get_medaka(self):
        return self.SOFTWARE_Medaka

    def get_medaka_name(self):
        return self.SOFTWARE_Medaka_name

    def get_medaka_name_consensus(self):
        return self.SOFTWARE_Medaka_name_consensus

    def get_medaka_name_variants(self):
        return self.SOFTWARE_Medaka_name_variants

    def get_medaka_name_extended_consensus(self):
        return self.SOFTWARE_Medaka_name_extended_consensus

    def get_medaka_name_extended_variants(self):
        return self.SOFTWARE_Medaka_name_extended_variants

    def get_medaka_parameters_consensus(self):
        return self.SOFTWARE_Medaka_PARAMETERS_consensus

    def get_medaka_parameters_variant(self):
        return self.SOFTWARE_Medaka_PARAMETERS_variant

    def get_medaka_default_model(self):
        return self.SOFTWARE_Medaka_default_model

    def get_medaka_remove_tags_model(self):
        return self.SOFTWARE_Medaka_remove_tags_model

    def get_medaka_version(self):
        return self.SOFTWARE_Medaka_VERSION

    def get_flye(self):
        return self.SOFTWARE_FLYE

    def get_flye_name(self):
        return self.SOFTWARE_FLYE_name

    def get_flye_name_extended(self):
        return self.SOFTWARE_FLYE_name_extended

    def get_flye_version(self):
        return self.SOFTWARE_FLYE_VERSION

    def get_flye_parameters(self):
        return self.SOFTWARE_FLYE_PARAMETERS

    ###
    ###### END minion software

    ### START PATHOGEN IDENTIFIGATION SOFTWARE

    def get_centrifuge(self):
        return self.SOFTWARE_CENTRIFUGE

    def get_centrifuge_name(self):
        return self.SOFTWARE_CENTRIFUGE_name

    def get_centrifuge_parameters(self):
        return self.SOFTWARE_CENTRIFUGE_parameters

    def get_kraken2(self):
        return self.SOFTWARE_KRAKEN2

    def get_kraken2_name(self):
        return self.SOFTWARE_KRAKEN2_name

    def get_kraken2_parameters(self):
        return self.SOFTWARE_KRAKEN2_parameters

    def get_krakenuniq_parameters(self):
        return self.SOFTWARE_KRAKENUNIQ_parameters

    def get_krakenuniq_name(self):
        return self.SOFTWARE_KRAKENUNIQ_name

    def get_diamond(self):
        return self.SOFTWARE_DIAMOND

    def get_diamond_name(self):
        return self.SOFTWARE_DIAMOND_name

    def get_kaiju(self):
        return self.SOFTWARE_KAIJU

    def get_kaiju_name(self):
        return self.SOFTWARE_KAIJU_name

    def get_kaiju_parameters(self):
        return self.SOFTWARE_KAIJU_parameters

    def get_desamba_name(self):
        return self.SOFTWARE_DESAMBA_name

    def get_desamba_parameters(self):
        return self.SOFTWARE_DESAMBA_parameters

    def get_desamba_name(self):
        return self.SOFTWARE_DESAMBA_name

    def get_blast(self):
        return self.SOFTWARE_BLAST

    def get_blast_name(self):
        return self.SOFTWARE_BLAST_name

    def get_blast_parameters(self):
        return self.SOFTWARE_BLAST_parameters

    def get_raven(self):
        return self.SOFTWARE_RAVEN

    def get_raven_name(self):
        return self.SOFTWARE_RAVEN_name

    def get_raven_parameters(self):
        return self.SOFTWARE_RAVEN_parameters

    def get_fastviromeexplorer(self):
        return self.SOFTWARE_FASTVIROMEEXPLORER

    def get_fastviromeexplorer_name(self):
        return self.SOFTWARE_FASTVIROMEEXPLORER_name

    def get_fastviromeexplorer_parameters(self):
        return self.SOFTWARE_FASTVIROMEEXPLORER_parameters

    ###################################
    #####
    #####    Global parameters for INSaFLU
    #####
    """
    ### it's a way to define a threshold to mask consensus sequences, with regions with low coverage, obtained by snippy
    ### It's a global parameter for INSaFLU
    """
    ### Define the Minimum percentage of horizontal coverage to generate consensus, otherwise drop sequence
    ##def get_insaflu_parameter_mask_consensus(self): return ""
    def get_insaflu_parameter_mask_consensus_name(self):
        return self.INSAFLU_PARAMETER_MASK_CONSENSUS_name

    def get_insaflu_parameter_mask_consensus_name_extended(self):
        return self.INSAFLU_PARAMETER_MASK_CONSENSUS_name_extended

    def get_insaflu_parameter_mask_consensus_vesion(self):
        return self.INSAFLU_PARAMETER_MASK_CONSENSUS_vesion

    def get_insaflu_parameter_mask_consensus_parameters(self):
        return self.INSAFLU_PARAMETER_MASK_CONSENSUS_parameters

    ### Define the minimum coverage
    ##def get_insaflu_parameter_limit_coverage(self): return ""
    def get_insaflu_parameter_limit_coverage_name(self):
        return self.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name

    def get_insaflu_parameter_limit_coverage_name_extended(self):
        return self.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name_extended

    def get_insaflu_parameter_limit_coverage_vesion(self):
        return self.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_vesion

    def get_insaflu_parameter_limit_coverage_parameters(self):
        return self.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_parameters

    ### Define the minimum freq
    ##def get_insaflu_parameter_limit_freq_cov(self): return ""
    def get_insaflu_parameter_freq_vcf_name(self):
        return self.INSAFLU_PARAMETER_VCF_FREQ_ONT_name

    def get_insaflu_parameter_freq_vcf_name_extended(self):
        return self.INSAFLU_PARAMETER_VCF_FREQ_ONT_name_extended

    def get_insaflu_parameter_freq_vcf_vesion(self):
        return self.INSAFLU_PARAMETER_VCF_FREQ_ONT_vesion

    def get_insaflu_parameter_freq_vcf_parameters(self):
        return self.INSAFLU_PARAMETER_VCF_FREQ_ONT_parameters

    ### get pangolin
    def get_pangolin_env(self):
        return self.SOFTWARE_Pangolin_Env

    def get_pangolin(self):
        return self.SOFTWARE_Pangolin

    def get_pangolin_name(self):
        return self.SOFTWARE_Pangolin_name

    def get_pangolin_version(self, software):
        if software == SoftwareNames.SOFTWARE_Pangolin_name:
            return SoftwareNames.SOFTWARE_Pangolin_VERSION
        return ""

    ## version must be obtain from ManagingFiles.Software.Version

    def get_pangolin_all_names_version(self):
        dt_return_version = {}
        for name in self.VECT_PANGOLIN_TO_TEST:
            dt_return_version[name] = self.get_pangolin_version(name)
        return dt_return_version

    ## version must be obtain from ManagingFiles.Software.Version
