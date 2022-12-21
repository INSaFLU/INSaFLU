"""
Created on Oct 28, 2017

@author: mmp
"""
import cmd
import datetime
import gzip
import logging
import os
import re
import subprocess

from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord
from constants.constants import Constants, FileExtensions, FileType, TypePath
from constants.meta_key_and_values import MetaKeyAndValue
from constants.software_names import SoftwareNames
from django.conf import settings
from django.template.defaultfilters import filesizeformat
from ete3 import Tree
from manage_virus.models import UploadFile
from manage_virus.uploadFiles import UploadFiles
from managing_files.manage_database import ManageDatabase
from managing_files.models import (MixedInfectionsTag, ProcessControler,
                                   ProjectSample, Reference, Sample)
from settings.constants_settings import ConstantsSettings
from settings.default_parameters import DefaultParameters
from settings.default_software_project_sample import DefaultProjectSoftware

from utils.coverage import DrawAllCoverage
from utils.mixed_infections_management import MixedInfectionsManagement
from utils.parse_coverage_file import GetCoverage
from utils.parse_out_files import ParseOutFiles
from utils.process_SGE import ProcessSGE
from utils.result import (CountHits, DecodeObjects, KeyValue, MaskingConsensus,
                          Result, ResultAverageAndNumberReads, SoftwareDesc)
from utils.utils import Utils


class Software(object):
    """
    classdocs
    """

    utils = Utils()
    software_names = SoftwareNames()

    ## logging
    logger_debug = logging.getLogger("fluWebVirus.debug")
    logger_production = logging.getLogger("fluWebVirus.production")

    def test_bgzip_and_tbi_in_vcf(self, vcf_file):
        """
        test if a a vcf file has a gzip file, if not create it
        """
        ### create the gzip file
        self.utils.compress_files(self.software_names.get_bgzip(), vcf_file)
        ### create the tabix
        if vcf_file.endswith(".gz"):
            self.create_index_files(vcf_file)
        else:
            self.create_index_files(vcf_file + ".gz")

    def get_vect_type_files_to_copy(self, software):
        """
        get type of files to copy
        """
        if software == SoftwareNames.SOFTWARE_SNIPPY_name:
            return [
                FileType.FILE_BAM,
                FileType.FILE_BAM_BAI,
                FileType.FILE_CONSENSUS_FA,
                FileType.FILE_DEPTH_GZ,
                FileType.FILE_DEPTH_GZ_TBI,
                FileType.FILE_TAB,
                FileType.FILE_VCF_GZ,
                FileType.FILE_VCF,
                FileType.FILE_VCF_GZ_TBI,
                FileType.FILE_CSV,
                FileType.FILE_REF_FASTA,
            ]
        elif software == SoftwareNames.SOFTWARE_FREEBAYES_name:
            return [FileType.FILE_VCF, FileType.FILE_TAB]
        elif software == SoftwareNames.SOFTWARE_Medaka_name:
            return [
                FileType.FILE_BAM,
                FileType.FILE_BAM_BAI,
                FileType.FILE_CONSENSUS_FA,
                FileType.FILE_DEPTH_GZ,
                FileType.FILE_DEPTH_GZ_TBI,
                FileType.FILE_TAB,
                FileType.FILE_VCF_GZ,
                FileType.FILE_VCF,
                FileType.FILE_VCF_GZ_TBI,
            ]

    def copy_files_to_project(self, project_sample, software, path_from):
        """
        copy temp files to the project_sample software names
        software : SOFTWARE_SNIPPY_name, SOFTWARE_FREEBAYES_name
        """
        for type_file in self.get_vect_type_files_to_copy(software):
            if type_file == FileType.FILE_CONSENSUS_FA:  ## if .fa file pass to .fasta
                self.utils.copy_file(
                    os.path.join(
                        path_from,
                        os.path.basename(
                            project_sample.get_file_output(
                                TypePath.MEDIA_ROOT, type_file, software
                            )
                        ),
                    ),
                    project_sample.get_file_output(
                        TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, software
                    ),
                )
            elif (
                type_file == FileType.FILE_VCF
                and software == SoftwareNames.SOFTWARE_FREEBAYES_name
            ):  ## vcf file
                ### create the bgzip file
                self.utils.compress_files(
                    self.software_names.get_bgzip(),
                    os.path.join(
                        path_from,
                        os.path.basename(
                            project_sample.get_file_output(
                                TypePath.MEDIA_ROOT, type_file, software
                            )
                        ),
                    ),
                )
                ### create the tabix
                self.create_index_files(
                    os.path.join(
                        path_from,
                        os.path.basename(
                            project_sample.get_file_output(
                                TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, software
                            )
                        ),
                    )
                )

                ### copy both
                self.utils.copy_file(
                    os.path.join(
                        path_from,
                        os.path.basename(
                            project_sample.get_file_output(
                                TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, software
                            )
                        ),
                    ),
                    project_sample.get_file_output(
                        TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, software
                    ),
                )
                self.utils.copy_file(
                    os.path.join(
                        path_from,
                        os.path.basename(
                            project_sample.get_file_output(
                                TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ_TBI, software
                            )
                        ),
                    ),
                    project_sample.get_file_output(
                        TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ_TBI, software
                    ),
                )

                ### if snippy copy also the vcf
                self.utils.copy_file(
                    os.path.join(
                        path_from,
                        os.path.basename(
                            project_sample.get_file_output(
                                TypePath.MEDIA_ROOT, type_file, software
                            )
                        ),
                    ),
                    project_sample.get_file_output(
                        TypePath.MEDIA_ROOT, type_file, software
                    ),
                )
            elif type_file == FileType.FILE_REF_FASTA:  ## this is only work for Snippy
                self.utils.copy_file(
                    os.path.join(
                        path_from,
                        "reference",
                        os.path.basename(
                            project_sample.get_file_output(
                                TypePath.MEDIA_ROOT, type_file, software
                            )
                        ),
                    ),
                    project_sample.get_file_output(
                        TypePath.MEDIA_ROOT, type_file, software
                    ),
                )
                ## create the FAI index
                self.create_fai_fasta(
                    project_sample.get_file_output(
                        TypePath.MEDIA_ROOT, type_file, software
                    )
                )
            else:
                self.utils.copy_file(
                    os.path.join(
                        path_from,
                        os.path.basename(
                            project_sample.get_file_output(
                                TypePath.MEDIA_ROOT, type_file, software
                            )
                        ),
                    ),
                    project_sample.get_file_output(
                        TypePath.MEDIA_ROOT, type_file, software
                    ),
                )

    def create_fai_fasta(self, fileFastaName):
        """
        Create fai for a fasta file
        """
        cmd = "%s faidx %s" % (self.software_names.get_samtools(), fileFastaName)
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run samtools")
        return cmd

    def create_index_bam(self, file_bam):
        """
        Create bai for a bam file
        """
        cmd = "%s index %s" % (self.software_names.get_samtools(), file_bam)
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run samtools")
        return cmd

    def get_statistics_bam(self, bam_file):
        """
        Create stats from bam
        ### ONT
        33022 + 0 in total (QC-passed reads + QC-failed reads)
        0 + 0 secondary
        0 + 0 supplementary
        0 + 0 duplicates
        33022 + 0 mapped (100.00% : N/A)
        0 + 0 paired in sequencing
        0 + 0 read1
        0 + 0 read2
        0 + 0 properly paired (N/A : N/A)
        0 + 0 with itself and mate mapped
        0 + 0 singletons (N/A : N/A)
        0 + 0 with mate mapped to a different chr
        0 + 0 with mate mapped to a different chr (mapQ>=5)
        ### Illumina
        40714 + 0 in total (QC-passed reads + QC-failed reads)
        0 + 0 secondary
        0 + 0 supplementary
        0 + 0 duplicates
        40714 + 0 mapped (100.00% : N/A)
        40714 + 0 paired in sequencing
        20367 + 0 read1
        20347 + 0 read2
        36712 + 0 properly paired (90.17% : N/A)
        40676 + 0 with itself and mate mapped
        38 + 0 singletons (0.09% : N/A)
        0 + 0 with mate mapped to a different chr
        0 + 0 with mate mapped to a different chr (mapQ>=5)
        """
        temp_file = self.utils.get_temp_file("bam_stat", FileExtensions.FILE_TXT)
        cmd = "%s flagstat %s > %s" % (
            self.software_names.get_samtools(),
            bam_file,
            temp_file,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run samtools")

        ##
        result = Result()
        if os.path.exists(temp_file):
            vect_result = self.utils.read_text_file(temp_file)
            self.utils.remove_file(temp_file)

            ### parse result
            for line in vect_result:
                sz_temp = line.strip()
                if len(sz_temp) == 0:
                    continue
                if sz_temp.find("in total") != -1:
                    result.add_key_value(
                        [
                            MetaKeyAndValue.SAMTOOLS_flagstat_total_reads,
                            sz_temp.split("+")[0].strip(),
                        ]
                    )
                if sz_temp.find("mapped (") != -1:
                    result.add_key_value(
                        [
                            MetaKeyAndValue.SAMTOOLS_flagstat_mapped_reads,
                            sz_temp.split("+")[0].strip(),
                        ]
                    )
        return result

    def creat_new_reference_to_snippy(self, project_sample):

        ### get temp file
        temp_file = self.utils.get_temp_file("new_reference", FileExtensions.FILE_FASTA)
        cmd = "perl %s %s %s" % (
            self.software_names.get_create_new_reference_to_snippy(),
            project_sample.project.reference.get_reference_gbk(TypePath.MEDIA_ROOT),
            temp_file,
        )
        #         print(cmd)
        exist_status = os.system(cmd)
        if exist_status != 0:
            if os.path.exists(temp_file):
                os.unlink(temp_file)
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run create_new_reference_to_snippy")
        path_dest = project_sample.get_file_output(
            TypePath.MEDIA_ROOT,
            FileType.FILE_REF_FASTA,
            self.software_names.get_snippy_name(),
        )
        #         print("Copy {} to {}".format(temp_file, path_dest))
        self.utils.copy_file(temp_file, path_dest)
        self.create_fai_fasta(path_dest)
        if os.path.exists(temp_file):
            os.unlink(temp_file)

    def run_spades(self, fastq_1, fastq_2, out_dir):
        """
        Run spades
        IF you have problems running spades.py change the spades.py file from:
        #!/usr/bin/env python
        to
        #!/usr/bin/env python3
        """
        if not os.path.exists(fastq_1):
            self.logger_production.error("Fastq 1 not found: " + fastq_1)
            self.logger_debug.error("Fastq 1 not found: " + fastq_1)
            raise Exception("Fastq 1 not found: " + fastq_1)

        if fastq_2 is None or len(fastq_2) == 0 or not os.path.exists(fastq_2):
            cmd = "%s -s %s %s -t %d -o %s" % (
                self.software_names.get_spades(),
                fastq_1,
                self.software_names.get_spades_parameters_single(),
                settings.THREADS_TO_RUN_FAST,
                out_dir,
            )
        else:
            cmd = "%s --pe1-1 %s --pe1-2 %s %s -t %d -o %s" % (
                self.software_names.get_spades(),
                fastq_1,
                fastq_2,
                self.software_names.get_spades_parameters(),
                settings.THREADS_TO_RUN_FAST,
                out_dir,
            )

        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run spades")
        return cmd

    def run_flye(self, fastq_1, out_dir):
        """
        Run flye
        """
        if not os.path.exists(fastq_1):
            self.logger_production.error("Fastq 1 not found: " + fastq_1)
            self.logger_debug.error("Fastq 1 not found: " + fastq_1)
            raise Exception("Fastq 1 not found: " + fastq_1)

        cmd = "%s --nano-raw %s --threads %d --out-dir %s %s" % (
            self.software_names.SOFTWARE_FLYE,
            fastq_1,
            settings.THREADS_TO_RUN_FAST,
            out_dir,
            self.software_names.SOFTWARE_FLYE_PARAMETERS,
        )

        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run flye")

        orig = os.path.join(out_dir, "assembly.fasta")
        dest = os.path.join(out_dir, "contigs.fasta")
        self.utils.copy_file(orig, dest)

        return cmd

    def convert_fastq_to_fasta(self, fastq_1, fasta_out_file):
        """
        Convert fastq to Fasta
        """
        if not os.path.exists(fastq_1):
            self.logger_production.error("Fastq 1 not found: " + fastq_1)
            self.logger_debug.error("Fastq 1 not found: " + fastq_1)
            raise Exception("Fastq 1 not found: " + fastq_1)

        cmd = "gzip -cd {} | sed -n '1~4s/^@/>/p;2~4p' > {}".format(
            fastq_1, fasta_out_file
        )

        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run fastq to fasta")

        ### need create a fasta file
        return cmd

    def is_exist_database_abricate(self, database):
        """
        DATABASE    SEQUENCES    DATE
        argannot    1749    2017-Oct-30
        card    2158    2017-Oct-30
        """
        temp_file = self.utils.get_temp_file("list_abricate", FileExtensions.FILE_TXT)
        cmd = "{} --list > {}".format(self.software_names.get_abricate(), temp_file)
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to collect list from abricate")
        vect_out = self.utils.read_text_file(temp_file)
        if os.path.exists(temp_file):
            os.unlink(temp_file)
        for line in vect_out:
            for field in line.split("\t"):
                if field.lower() == database.lower():
                    return True
        return False

    def create_database_abricate(self, database, file_name):
        """
        create a database on abricate
        """
        if not os.path.isfile(file_name):
            raise IOError("File not found: " + file_name)
        cmd = "mkdir -p %s/%s" % (self.software_names.SOFTWARE_ABRICATE_DB, database)
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run make directory")

        ## copy the file
        self.utils.copy_file(
            file_name,
            os.path.join(
                self.software_names.SOFTWARE_ABRICATE_DB, database, "sequences"
            ),
        )

        cmd = "%s --setupdb" % (self.software_names.get_abricate())
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run abricate --setupdb")

    def run_abricate(self, database, file_name, parameters, out_file):
        """
        Run abricator
        """
        cmd = "%s --db %s %s --quiet %s > %s" % (
            self.software_names.get_abricate(),
            database,
            parameters,
            file_name,
            out_file,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run abricate")
        return cmd

    """
    Global processing
    """
    #     @transaction.atomic
    def identify_type_and_sub_type(
        self, sample, fastq1_1, fastq1_2, owner, b_run_tests=False
    ):
        """
        Identify type and sub_type
        Because of the tests need to pass the files also as parameters
        """

        manageDatabase = ManageDatabase()
        ### temp dir out spades
        out_dir_result = self.utils.get_temp_dir()
        result_all = Result()

        ### test files
        if fastq1_1 is None or not os.path.exists(fastq1_1):
            return False
        if not fastq1_2 is None and not os.path.exists(fastq1_2):
            return False

        print("identify_type_and_sub_type")
        if sample.is_type_fastq_gz_sequencing():  ## illumina
            try:
                cmd = self.run_spades(fastq1_1, fastq1_2, out_dir_result)
                parameters = self.software_names.get_spades_parameters()
                if fastq1_2 == None or len(fastq1_2) == 0:
                    parameters = self.software_names.get_spades_parameters_single()
                result_all.add_software(
                    SoftwareDesc(
                        self.software_names.get_spades_name(),
                        self.software_names.get_spades_version(),
                        parameters,
                    )
                )
            except Exception:
                result = Result()
                result.set_error(
                    "Spades (%s) fail to run"
                    % (self.software_names.get_spades_version())
                )
                result.add_software(
                    SoftwareDesc(
                        self.software_names.get_spades_name(),
                        self.software_names.get_spades_version(),
                        self.software_names.get_spades_parameters(),
                    )
                )
                manageDatabase.set_sample_metakey(
                    sample,
                    owner,
                    MetaKeyAndValue.META_KEY_Identify_Sample,
                    MetaKeyAndValue.META_VALUE_Error,
                    result.to_json(),
                )
                self.utils.remove_dir(out_dir_result)
                return False
        else:  ### for minion
            try:
                cmd = self.convert_fastq_to_fasta(
                    fastq1_1, os.path.join(out_dir_result, "contigs.fasta")
                )

            except Exception:
                result = Result()
                result.set_error("Fail to convert fastq to fasta.")
                result.add_software(SoftwareDesc("sed", "", ""))
                manageDatabase.set_sample_metakey(
                    sample,
                    owner,
                    MetaKeyAndValue.META_KEY_Identify_Sample,
                    MetaKeyAndValue.META_VALUE_Error,
                    result.to_json(),
                )
                self.utils.remove_dir(out_dir_result)
                return False

        file_out_contigs = os.path.join(out_dir_result, "contigs.fasta")
        if (
            not os.path.exists(file_out_contigs)
            or os.path.getsize(file_out_contigs) < 50
        ):
            ## save error in MetaKeySample
            result = Result()
            if sample.is_type_fastq_gz_sequencing():
                result.set_error(
                    "Spades (%s) fail to run, empty contigs.fasta file."
                    % (self.software_names.get_spades_version())
                )
                result.add_software(
                    SoftwareDesc(
                        self.software_names.get_spades_name(),
                        self.software_names.get_spades_version(),
                        self.software_names.get_spades_parameters(),
                    )
                )
            else:
                result.set_error(
                    "Low number of reads in fasta file. Came from fastq.gz"
                )
                result.add_software(SoftwareDesc("sed", "", ""))
            manageDatabase.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_Identify_Sample,
                MetaKeyAndValue.META_VALUE_Error,
                result.to_json(),
            )
            self.utils.remove_dir(out_dir_result)
            return False

        ### test id abricate has the database
        try:
            uploadFile = UploadFile.objects.order_by("-version")[0]
        except UploadFile.DoesNotExist:
            ## save error in MetaKeySample
            result = Result()
            result.set_error(
                "Abricate (%s) fail to run"
                % (self.software_names.get_abricate_version())
            )
            result.add_software(
                SoftwareDesc(
                    self.software_names.get_abricate_name(),
                    self.software_names.get_abricate_version(),
                    self.software_names.get_abricate_parameters(),
                )
            )
            manageDatabase.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_Identify_Sample,
                MetaKeyAndValue.META_VALUE_Error,
                result.to_json(),
            )
            self.utils.remove_dir(out_dir_result)
            return False

        if not self.is_exist_database_abricate(uploadFile.abricate_name):
            try:
                self.create_database_abricate(uploadFile.abricate_name, uploadFile.path)
            except Exception:
                result = Result()
                result.set_error(
                    "Abricate (%s) fail to run --setupdb"
                    % (self.software_names.get_abricate_version())
                )
                result.add_software(
                    SoftwareDesc(
                        self.softwafile_outre_names.get_abricate_name(),
                        self.software_names.get_abricate_version(),
                        self.software_names.get_abricate_parameters(),
                    )
                )
                manageDatabase.set_sample_metakey(
                    sample,
                    owner,
                    MetaKeyAndValue.META_KEY_Identify_Sample,
                    MetaKeyAndValue.META_VALUE_Error,
                    result.to_json(),
                )
                self.utils.remove_dir(out_dir_result)
                return False

        ## run abricate
        out_file_abricate = self.utils.get_temp_file("temp_abricate", ".txt")
        try:
            cmd = self.run_abricate(
                uploadFile.abricate_name,
                file_out_contigs,
                SoftwareNames.SOFTWARE_ABRICATE_PARAMETERS,
                out_file_abricate,
            )
            result_all.add_software(
                SoftwareDesc(
                    self.software_names.get_abricate_name(),
                    self.software_names.get_abricate_version(),
                    self.software_names.get_abricate_parameters()
                    + " for type/subtype identification",
                )
            )
        except Exception:
            result = Result()
            result.set_error(
                "Abricate (%s) fail to run"
                % (self.software_names.get_abricate_version())
            )
            result.add_software(
                SoftwareDesc(
                    self.software_names.get_abricate_name(),
                    self.software_names.get_abricate_version(),
                    self.software_names.get_abricate_parameters()
                    + " for type/subtype identification",
                )
            )
            manageDatabase.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_Identify_Sample,
                MetaKeyAndValue.META_VALUE_Error,
                result.to_json(),
            )
            self.utils.remove_dir(out_dir_result)
            return False

        if not os.path.exists(out_file_abricate):
            ## save error in MetaKeySample
            result = Result()
            result.set_error(
                "Abricate (%s)identify_contigs fail to run"
                % (self.software_names.get_abricate_version())
            )
            result.add_software(
                SoftwareDesc(
                    self.software_names.get_abricate(),
                    self.software_names.get_abricate_version(),
                    self.software_names.get_abricate_parameters(),
                )
            )
            manageDatabase.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_Identify_Sample,
                MetaKeyAndValue.META_VALUE_Error,
                result.to_json(),
            )
            self.utils.remove_dir(out_dir_result)
            return False

        parseOutFiles = ParseOutFiles()
        (dict_data_out, clean_abricate_file) = parseOutFiles.parse_abricate_file(
            out_file_abricate,
            os.path.basename(sample.get_abricate_output(TypePath.MEDIA_ROOT)),
            SoftwareNames.SOFTWARE_SPAdes_CLEAN_HITS_BELLOW_VALUE,
        )

        ### set the identification in database
        uploadFiles = UploadFiles()
        vect_data = uploadFiles.uploadIdentifyVirus(
            dict_data_out, uploadFile.abricate_name
        )
        if len(vect_data) == 0:
            ## save error in MetaKeySample
            result = Result()
            result.set_error("Fail to identify type and sub type")
            result.add_software(
                SoftwareDesc(
                    self.software_names.get_abricate(),
                    self.software_names.get_abricate_version(),
                    self.software_names.get_abricate_parameters(),
                )
            )
            manageDatabase.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_Identify_Sample,
                MetaKeyAndValue.META_VALUE_Error,
                result.to_json(),
            )
        else:
            for identify_virus in vect_data:
                sample.identify_virus.add(identify_virus)
            sample.save()

        ## copy the abricate output
        self.utils.copy_file(
            clean_abricate_file, sample.get_abricate_output(TypePath.MEDIA_ROOT)
        )

        ## Only identify Contigs for Illuminua, because Spades runs. In ONT doesn't run because it is identify in reads.
        try:
            contigs_2_sequences = Contigs2Sequences(b_run_tests)
            (
                out_file_clean,
                clean_abricate_file,
            ) = contigs_2_sequences.identify_contigs(
                file_out_contigs,
                os.path.basename(
                    sample.get_draft_contigs_abricate_output(TypePath.MEDIA_ROOT)
                )
                if sample.is_type_fastq_gz_sequencing()
                else os.path.basename(
                    sample.get_draft_reads_abricate_output(TypePath.MEDIA_ROOT)
                ),
                True if sample.is_type_fastq_gz_sequencing() else False,
            )
            ## copy the contigs from spades
            if sample.is_type_fastq_gz_sequencing():  ## illumina
                if os.path.exists(out_file_clean):
                    self.utils.copy_file(
                        out_file_clean,
                        sample.get_draft_contigs_output(TypePath.MEDIA_ROOT),
                    )
                if os.path.exists(clean_abricate_file):
                    self.utils.copy_file(
                        clean_abricate_file,
                        sample.get_draft_contigs_abricate_output(TypePath.MEDIA_ROOT),
                    )
            else:
                print(
                    "Should be writing abricate results to {}".format(
                        sample.get_draft_reads_abricate_output(TypePath.MEDIA_ROOT)
                    )
                )
                if os.path.exists(clean_abricate_file):
                    self.utils.copy_file(
                        clean_abricate_file,
                        sample.get_draft_reads_abricate_output(TypePath.MEDIA_ROOT),
                    )
            result_all.add_software(
                SoftwareDesc(
                    self.software_names.get_abricate_name(),
                    self.software_names.get_abricate_version(),
                    self.software_names.get_abricate_parameters_mincov_30()
                    + " for segments/references assignment",
                )
            )
            if not out_file_clean is None:
                self.utils.remove_file(out_file_clean)
        except Exception as e:
            result = Result()
            result.set_error(
                "Abricate (%s) fail to run"
                % (self.software_names.get_abricate_version())
            )
            result.add_software(
                SoftwareDesc(
                    self.software_names.get_abricate_name(),
                    self.software_names.get_abricate_version(),
                    self.software_names.get_abricate_parameters_mincov_30()
                    + " for segments/references assignment",
                )
            )
            manageDatabase.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_Identify_Sample,
                MetaKeyAndValue.META_VALUE_Error,
                result.to_json(),
            )
            return False

        ## save everything OK
        if sample.is_type_fastq_gz_sequencing():
            manageDatabase.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_Identify_Sample,
                MetaKeyAndValue.META_VALUE_Success,
                "Success, Spades(%s), Abricate(%s)"
                % (
                    self.software_names.get_spades_version(),
                    self.software_names.get_abricate_version(),
                ),
            )
        else:
            manageDatabase.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_Identify_Sample,
                MetaKeyAndValue.META_VALUE_Success,
                "Success, Abricate(%s)" % (self.software_names.get_abricate_version()),
            )
        manageDatabase.set_sample_metakey(
            sample,
            owner,
            MetaKeyAndValue.META_KEY_Identify_Sample_Software,
            MetaKeyAndValue.META_VALUE_Success,
            result_all.to_json(),
        )
        self.utils.remove_dir(out_dir_result)
        self.utils.remove_file(out_file_abricate)
        self.utils.remove_file(clean_abricate_file)
        return True

    def get_species_tag(self, reference):
        """
        :param reference instance, database
        :return specie TAGs: Reference.SPECIES_SARS_COV_2; Reference.SPECIES_MPXV; Reference.SPECIES_INFLUENZA, ...
        """

        ### Done
        if (
            len(reference.specie_tag) > 0
            and reference.specie_tag != Reference.SPECIES_NOT_SET
        ):
            return reference.specie_tag

        ### test id abricate has the database
        try:
            uploadFile = UploadFile.objects.order_by("-version")[0]
        except UploadFile.DoesNotExist:
            return Reference.SPECIES_NOT_SET

        ## if not exist try to create it
        if not self.is_exist_database_abricate(uploadFile.abricate_name):
            try:
                self.create_database_abricate(uploadFile.abricate_name, uploadFile.path)
            except Exception:
                return Reference.SPECIES_NOT_SET

        ## run abricate
        out_file_abricate = self.utils.get_temp_file("temp_abricate", ".txt")
        try:
            cmd = self.run_abricate(
                uploadFile.abricate_name,
                reference.get_reference_fasta(TypePath.MEDIA_ROOT),
                SoftwareNames.SOFTWARE_ABRICATE_PARAMETERS,
                out_file_abricate,
            )
        except Exception:
            return Reference.SPECIES_NOT_SET

        if not os.path.exists(out_file_abricate):
            return Reference.SPECIES_NOT_SET

        parseOutFiles = ParseOutFiles()
        (vect_data, clean_abricate_file) = parseOutFiles.parse_abricate_file(
            out_file_abricate,
            "doesnt_mather.txt",
            SoftwareNames.SOFTWARE_SPAdes_CLEAN_HITS_BELLOW_VALUE,
        )
        ## don't need this file, only need vect_data_file
        self.utils.remove_file(clean_abricate_file)
        self.utils.remove_file(out_file_abricate)

        uploadFiles = UploadFiles()
        try:
            ## could fail some identification, so, return False if it occurs
            vect_data = uploadFiles.uploadIdentifyVirus(
                vect_data, uploadFile.abricate_name
            )
            if len(vect_data) == 0:
                return Reference.SPECIES_NOT_SET
        except Exception as e:
            return Reference.SPECIES_NOT_SET

        ### test number of right segments
        number_right_beta_cov, number_right_mpxv, number_right_influenza = (0, 0, 0)
        for identify_virus in vect_data:
            if (
                identify_virus.seq_virus.name == "BetaCoV"
                and identify_virus.seq_virus.kind_type.name == "Genus"
            ):
                number_right_beta_cov += 1
            ### need to read from the file
            elif (
                identify_virus.seq_virus.name
                in (
                    "SARS_CoV_2",
                    "SARS_CoV",
                    "SCoV2_potential_Omicron",
                    "HCoV_OC43",
                    "HCoV_HKU1",
                    "MERS_CoV",
                )
                and identify_virus.seq_virus.kind_type.name == "Human"
            ):
                number_right_beta_cov += 1
            elif (
                identify_virus.seq_virus.name in ("A", "B")
                and identify_virus.seq_virus.kind_type.name == "Type"
            ):
                number_right_influenza += 1
            elif (
                identify_virus.seq_virus.name in ("Yamagata", "Victoria")
                and identify_virus.seq_virus.kind_type.name.lower() == "lineage"
            ):
                number_right_influenza += 1
            elif identify_virus.seq_virus.kind_type.name.lower() == "subtype" and (
                identify_virus.seq_virus.name.startswith("H")
                or identify_virus.seq_virus.name.startswith("N")
            ):
                number_right_influenza += 1
            elif (
                identify_virus.seq_virus.name == Reference.SPECIES_MPXV
                and identify_virus.seq_virus.kind_type.name.lower() == "species"
            ):
                number_right_mpxv += 1

        ## if right at least two
        if number_right_beta_cov > 0:
            reference.specie_tag = Reference.SPECIES_SARS_COV_2
            reference.save()
            return Reference.SPECIES_SARS_COV_2
        elif number_right_mpxv > 0:
            reference.specie_tag = Reference.SPECIES_MPXV
            reference.save()
            return Reference.SPECIES_MPXV
        elif number_right_influenza > 0:
            reference.specie_tag = Reference.SPECIES_INFLUENZA
            reference.save()
            return Reference.SPECIES_INFLUENZA

        ## Not ser
        return Reference.SPECIES_NOT_SET

    def run_fastq(self, file_name_1, file_name_2):
        """
        run fastQ, return output directory
        -o OUT_FOLDER --nogroup --format fastq --threads 10 --dir OUT_FOLDER FILE1 FILE2
        """
        temp_dir = self.utils.get_temp_dir()
        if not file_name_2 is None and len(file_name_2) > 0:
            cmd = "%s -o %s --nogroup --format fastq --threads %d --dir %s %s %s" % (
                self.software_names.get_fastq(),
                temp_dir,
                settings.THREADS_TO_RUN_FASTQC,
                temp_dir,
                file_name_1,
                file_name_2,
            )
        else:
            cmd = "%s -o %s --nogroup --format fastq --threads %d --dir %s %s" % (
                self.software_names.get_fastq(),
                temp_dir,
                settings.THREADS_TO_RUN_FASTQC,
                temp_dir,
                file_name_1,
            )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run run_fastq")
        return temp_dir

    def run_prokka(self, fasta_file_name, original_file_name):
        """
        run prokka, in fasta file
        out: genbank
        {bpipe_prokka} FILE1 --kingdom Viruses --locustag locus --kingdom Viruses --locustag locus --genus Influenzavirus
                --species Influenzavirus --strain ref_PREFIX_FILES_OUT --outdir OUT_FOLDER/PREFIX_FILES_OUT --prefix PREFIX_FILES_OUT
        """
        if not os.path.exists(fasta_file_name) or os.path.getsize(fasta_file_name) == 0:
            raise Exception("File doesn't exist")
        temp_dir = self.utils.get_temp_dir()
        name_strain = self.utils.clean_extension(os.path.basename(original_file_name))
        cmd = "{} {} {} --strain {} --force --outdir {} --prefix {}".format(
            self.software_names.get_prokka(),
            fasta_file_name,
            self.software_names.get_prokka_parameters(),
            name_strain,
            temp_dir,
            name_strain,
        )
        print(cmd)
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run prokka")

        ## clean /mol_type=
        gbk_file = os.path.join(
            temp_dir,
            self.utils.clean_extension(original_file_name) + FileExtensions.FILE_GBK,
        )
        temp_file = self.utils.get_temp_file_from_dir(
            temp_dir, "new_file", FileExtensions.FILE_GBK
        )
        cmd = "grep -E -v '/mol_type=' {} > {}".format(gbk_file, temp_file)
        os.system(cmd)
        self.utils.move_file(temp_file, gbk_file)
        return temp_dir

    def run_mauve(self, dir_to_process, out_file):
        """
        run mauve
        out: out_file
        --output=alignment_all_96_samples_H3.xmfa *fasta
        """
        dir_present = os.getcwd()
        os.chdir(dir_to_process)
        cmd = "{} --output={} *fasta".format(self.software_names.get_mauve(), out_file)
        exist_status = os.system(cmd)
        os.chdir(dir_present)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run progressive mauve")
        return out_file

    def run_convert_mauve(self, input_file, out_file):
        """
        run convert mauve
        out: out_file
        """
        temp_file = self.utils.get_temp_file(
            "convert_mauve_software", FileExtensions.FILE_FASTA
        )
        cmd = "perl {} -c -i {} -o {} -f fasta".format(
            self.software_names.get_convert_mauve(), input_file, temp_file
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run convert mauve")

        #### clean names
        ### important, must have this order
        vect_names_to_clean = [
            FileExtensions.FILE_CONSENSUS_FASTA,
            FileExtensions.FILE_FASTA,
            FileExtensions.FILE_FA,
        ]
        self.utils.clean_fasta_names(vect_names_to_clean, temp_file, out_file)
        os.unlink(temp_file)
        return out_file

    def run_mafft(self, input_file, out_file, parameters):
        """
        run mafft
        out: out_file
        """
        cmd = "{}; {} {} --thread {} {} > {}".format(
            self.software_names.get_mafft_set_env_variable(),
            self.software_names.get_mafft(),
            parameters,
            settings.THREADS_TO_RUN_SLOW,
            input_file,
            out_file,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run mafft")
        return out_file

    def run_clustalo(self, input_file, out_file, parameters=""):
        """
        run clustalo
        out: out_file
        """
        cmd = "{} --force --infmt=fa --outfmt=fa --seqtype dna --MAC-RAM 8000 {} --threads={} -i {} -o {}".format(
            self.software_names.get_clustalo(),
            parameters,
            settings.THREADS_TO_RUN_SLOW,
            input_file,
            out_file,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run mafft")
        return out_file

    def run_seqret_nex(self, input_file, out_file):
        """
        run mafft
        out: out_file
        """
        cmd = "{} -sequence {} {} -outseq {}".format(
            self.software_names.get_seqret(),
            input_file,
            self.software_names.get_seqret_nex_parameters(),
            out_file,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run seqret")
        return out_file

    def run_fasttree(self, input_file, out_file, parameters, root_leaf=None):
        """
        run fasttree, can pass a root leaf
        out: out_file
        """
        cmd = "{} {} {} > {}".format(
            self.software_names.get_fasttree(), parameters, input_file, out_file
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run fasttree")

        ## reroot the tree
        if not root_leaf is None:
            tree = Tree(out_file)
            if len(tree.get_leaves_by_name(root_leaf)) > 0:
                leaf = tree.get_leaves_by_name(root_leaf)[0]
                tree.set_outgroup(leaf)
                tree.write(outfile=out_file)
        return out_file

    def run_trimmomatic(self, file_name_1, file_name_2, sample, user=None):
        """
        run trimmomatic
        return output directory

        #${bpipe_trimmomatic} PE -threads 3 -basein FILE1 -baseout OUT_FOLDER/PREFIX_FILES_OUT.fastq.gz SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:55 TOPHRED33

        PE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
                or:
        SE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] <inputFile> <outputFile> <trimmer1>
        :out log.txt with data of the result
        """

        ## get dynamic parameters
        if user is None:
            parameters = self.software_names.get_trimmomatic_parameters()
        else:
            default_software_project = DefaultProjectSoftware()
            parameters = (
                default_software_project.get_trimmomatic_parameters_all_possibilities(
                    user, sample
                )
            )

        ### run software
        temp_dir = self.utils.get_temp_dir()
        if file_name_2 is None or len(file_name_2) == 0:
            cmd = "java -jar %s SE -threads %d %s %s_1P.fastq.gz %s" % (
                self.software_names.get_trimmomatic(),
                settings.THREADS_TO_RUN_FAST,
                file_name_1,
                os.path.join(temp_dir, sample.name),
                parameters,
            )
        else:
            ### need to make links the files to trimmomatic identify the _R1_ and _R2_
            new_file_name = os.path.join(temp_dir, "name_R1_001.fastq.gz")
            cmd = "ln -s {} {}".format(file_name_1, new_file_name)
            os.system(cmd)
            cmd = "ln -s {} {}".format(
                file_name_2, os.path.join(temp_dir, "name_R2_001.fastq.gz")
            )
            os.system(cmd)
            cmd = "java -jar %s PE -threads %d -basein %s -baseout %s.fastq.gz %s" % (
                self.software_names.get_trimmomatic(),
                settings.THREADS_TO_RUN_FAST,
                new_file_name,
                os.path.join(temp_dir, sample.name),
                parameters,
            )
        (exist_status, output) = subprocess.getstatusoutput(cmd)
        if exist_status != 0:
            self.utils.remove_dir(temp_dir)
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run trimmomatic")

        result = Result()
        for line in output.split("\n"):
            ### line of interested
            if line.startswith("Input Read Pairs:"):
                ### parse line
                ### Input Read Pairs: 44425 Both Surviving: 41254 (92,86%) Forward Only Surviving: 2306 (5,19%) Reverse Only Surviving: 431 (0,97%) Dropped: 434 (0,98%)
                for _ in range(
                    len(SoftwareNames.SOFTWARE_TRIMMOMATIC_vect_info_to_collect)
                ):
                    if _ + 1 == len(
                        SoftwareNames.SOFTWARE_TRIMMOMATIC_vect_info_to_collect
                    ):
                        result.add_key_value(
                            KeyValue(
                                SoftwareNames.SOFTWARE_TRIMMOMATIC_vect_info_to_collect[
                                    _
                                ],
                                line.split(
                                    SoftwareNames.SOFTWARE_TRIMMOMATIC_vect_info_to_collect[
                                        _
                                    ]
                                )[1].strip(),
                            )
                        )
                    else:
                        result.add_key_value(
                            KeyValue(
                                SoftwareNames.SOFTWARE_TRIMMOMATIC_vect_info_to_collect[
                                    _
                                ],
                                line.split(
                                    SoftwareNames.SOFTWARE_TRIMMOMATIC_vect_info_to_collect[
                                        _
                                    ]
                                )[1]
                                .split(
                                    SoftwareNames.SOFTWARE_TRIMMOMATIC_vect_info_to_collect[
                                        _ + 1
                                    ][
                                        0
                                    ]
                                )[0]
                                .strip(),
                            )
                        )
            elif line.startswith("Quality encoding detected as"):
                result.add_key_value(
                    KeyValue(
                        "Quality encoding detected:",
                        line.split("Quality encoding detected as")[1].strip(),
                    )
                )
        ### dt_info has the information
        return (temp_dir, result, parameters)

    def run_snippy(
        self,
        file_name_1,
        file_name_2,
        path_reference_fasta,
        path_reference_genbank,
        sample_name,
        snippy_parameters,
    ):
        """
        run snippy
        return output directory

        ## ./snippy --cpus 3 --outdir /tmp/insafli/xpto --prefix xpto --ref ~/fluWeb/ref_H3.fasta --mapqual 20 --mincov 10 --minfrac 0.51 --R1 ~/fluWeb/EVA001_S66_L001_R1_001.fastq.gz --R2 ~/fluWeb/EVA001_S66_L001_R2_001.fastq.gz SNIPPY

        Out file
        [06:29:16] * /tmp/insafli/xpto/xpto.aligned.fa
        [06:29:16] * /tmp/insafli/xpto/xpto.bam
        [06:29:16] * /tmp/insafli/xpto/xpto.bam.bai
        [06:29:16] * /tmp/insafli/xpto/xpto.bed
        [06:29:16] * /tmp/insafli/xpto/xpto.consensus.fa
        [06:29:16] * /tmp/insafli/xpto/xpto.consensus.subs.fa
        [06:29:16] * /tmp/insafli/xpto/xpto.csv
        [06:29:16] * /tmp/insafli/xpto/xpto.depth.gz
        [06:29:16] * /tmp/insafli/xpto/xpto.depth.gz.tbi
        [06:29:16] * /tmp/insafli/xpto/xpto.filt.subs.vcf
        [06:29:16] * /tmp/insafli/xpto/xpto.filt.subs.vcf.gz
        [06:29:16] * /tmp/insafli/xpto/xpto.filt.subs.vcf.gz.tbi
        [06:29:16] * /tmp/insafli/xpto/xpto.filt.vcf
        [06:29:16] * /tmp/insafli/xpto/xpto.gff
        [06:29:16] * /tmp/insafli/xpto/xpto.html
        [06:29:16] * /tmp/insafli/xpto/xpto.log
        [06:29:16] * /tmp/insafli/xpto/xpto.raw.vcf
        [06:29:16] * /tmp/insafli/xpto/xpto.tab
        [06:29:16] * /tmp/insafli/xpto/xpto.txt
        [06:29:16] * /tmp/insafli/xpto/xpto.vcf
        [06:29:16] * /tmp/insafli/xpto/xpto.vcf.gz
        [06:29:16] * /tmp/insafli/xpto/xpto.vcf.gz.tbi
        """
        temp_dir = os.path.join(self.utils.get_temp_dir(), sample_name)
        if file_name_2 is None or len(file_name_2) == 0:
            cmd = "%s --cpus %d --outdir %s --prefix %s --ref %s %s --se %s" % (
                self.software_names.get_snippy(),
                settings.THREADS_TO_RUN_SLOW,
                temp_dir,
                sample_name,
                path_reference_genbank,
                snippy_parameters,
                file_name_1,
            )
        else:
            cmd = "%s --cpus %d --outdir %s --prefix %s --ref %s %s --R1 %s --R2 %s" % (
                self.software_names.get_snippy(),
                settings.THREADS_TO_RUN_SLOW,
                temp_dir,
                sample_name,
                path_reference_genbank,
                snippy_parameters,
                file_name_1,
                file_name_2,
            )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            self.utils.remove_dir(temp_dir)
            raise Exception("Fail to run snippy")

        ### add FREQ to VCF file
        parse_out_files = ParseOutFiles()
        out_file_transformed_amino = parse_out_files.add_amino_single_letter_code(
            os.path.join(temp_dir, sample_name + ".vcf")
        )
        self.utils.add_freq_to_vcf(
            out_file_transformed_amino, os.path.join(temp_dir, sample_name + "_2.vcf")
        )
        self.utils.remove_file(out_file_transformed_amino)

        ### add FREQ and other things to TAB file
        self.run_snippy_vcf_to_tab_freq_and_evidence(
            path_reference_fasta,
            path_reference_genbank,
            os.path.join(temp_dir, sample_name + "_2.vcf"),
            os.path.join(temp_dir, sample_name + ".tab"),
        )
        return temp_dir

    def run_nextalign(
        self, reference_fasta, sequences, gff_file, genes_to_process, output_path
    ):
        """
        :reference_fasta only one sequence
        :param genes_to_process -> has the name of the genes to process, if has no values don't produce these files
        output path: has all files (align, genes.fasta)
        nextalign --sequences sequences.fasta --reference SARS_CoV_2_COVID_19_Wuhan_Hu_1_MN908947.fasta --genemap covid.gff3 --genes S --output-dir output
        :out output path, alignment file, [gene file names]
        """

        if len(gff_file) > 0 and len(genes_to_process) > 0:
            cmd = (
                "%s --reference %s --sequences %s --genemap %s --genes %s --output-dir %s"
                % (
                    self.software_names.get_nextalign(),
                    reference_fasta,
                    sequences,
                    gff_file,
                    genes_to_process,
                    output_path,
                )
            )
        else:
            cmd = "%s --reference %s --sequences %s --output-dir %s" % (
                self.software_names.get_nextalign(),
                reference_fasta,
                sequences,
                output_path,
            )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run nextalign")

        ### return file names
        vect_file_gene_out = []
        if len(genes_to_process) > 0:
            vect_file_gene_out = [
                sequences.replace(".fasta", ".gene.{}.fasta".format(gene_name))
                for gene_name in genes_to_process.split(",")
            ]
        return (
            sequences.replace(".fasta", ".aligned.fasta"),
            vect_file_gene_out,
            sequences.replace(".fasta", ".insertions.csv"),
        )

    def run_genbank2gff3(self, genbank, out_file, for_nextclade=False):
        """
        for_nextclade = True; need to add gene annotation and gene_name in INFO
        """
        temp_file = self.utils.get_temp_file("gbk_to_gff3", ".txt")

        ## set VERSION ID equal to ACCESSION
        out_file_gb = self.utils.get_temp_file("file_name", ".gb")
        self.utils.clean_genbank_version_name(genbank, out_file_gb)

        cmd = "perl {} ".format(
            SoftwareNames.SOFTWARE_genbank_to_perl
        ) + "-f GenBank {} -out stdout -x gene > {}".format(out_file_gb, temp_file)
        exist_status = os.system(cmd)
        if exist_status != 0:
            os.unlink(out_file_gb)
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run genbank2gff3")

        ##        with open(temp_file, "w") as out_handle, open(out_file_gb) as in_handle:
        ##            GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)

        ### remove clean version gb file
        os.unlink(out_file_gb)

        ### filter file
        #     vect_filter = ['remark', 'source', 'gene']
        vect_pass = ["CDS"]
        ##### CAVEAT...
        ### If you lead ID= in the gff the orf1ab will be continued (Transcript_Exon_MN908947_13468_21555|Coding|1/1|c.395C>T|p.Thr132Ile|p.T132I)
        ### if it is removed ID=  is going like (Transcript_orf1ab|Coding|1/1|c.13597C>T|p.His4533Tyr|p.H4533Y)
        ### the last one is equal to SNIPPY
        vect_remove_info = ["ID="]
        with open(temp_file) as handle, open(out_file, "w") as handle_write:
            for line in handle:
                sz_temp = line.strip()
                if (
                    len(sz_temp) == 0
                    or sz_temp.find("# Input") == 0
                    or sz_temp.find("# GFF3 saved") == 0
                ):
                    continue
                if sz_temp.find("##FASTA") == 0:
                    break
                if sz_temp[0] == "#":
                    handle_write.write(sz_temp + "\n")
                elif (
                    len(sz_temp.split("\t")) > 3 and sz_temp.split("\t")[2] in vect_pass
                ):
                    lst_data = sz_temp.split("\t")
                    if len(lst_data) > 8:
                        lst_data_info = lst_data[8].split(";")
                        has_gene, has_gene_name = False, False
                        gene, gene_name = "", ""
                        vect_remove_index = []
                        for _, data_ in enumerate(lst_data_info):
                            if data_.lower().startswith("gene="):
                                has_gene = True
                            if data_.lower().startswith("name="):
                                gene = "gene={}".format("=".join(data_.split("=")[1:]))

                            ### for next clade
                            if data_.lower().startswith("gene_name="):
                                has_gene_name = True
                            if data_.lower().startswith(
                                "gene="
                            ) or data_.lower().startswith("name="):
                                gene_name = "gene_name={}".format(
                                    "=".join(data_.split("=")[1:])
                                )
                            for to_remove in vect_remove_info:
                                if data_.lower().startswith(to_remove.lower()):
                                    vect_remove_index.append(_)
                                    break

                        ### remove index
                        vect_remove_index = sorted(vect_remove_index, reverse=True)
                        for remove_index in vect_remove_index:
                            lst_data_info.pop(remove_index)

                        ### zero base
                        if self.utils.is_integer(lst_data[7]) and int(lst_data[7]) > 0:
                            lst_data[7] = "{}".format(int(lst_data[7]) - 1)

                        ### add gene
                        if not has_gene and len(gene) > 0:
                            lst_data[8] = gene + ";" + ";".join(lst_data_info)
                            sz_temp = "\t".join(lst_data)

                        ### only for nextclade, that need to add gene
                        if for_nextclade:
                            if not has_gene_name and len(gene_name) > 0:
                                lst_data[8] = gene_name + ";" + ";".join(lst_data_info)
                            lst_data[2] = "gene"
                            handle_write.write("\t".join(lst_data) + "\n")
                    handle_write.write(sz_temp + "\n")
        os.unlink(temp_file)
        return out_file

    def run_genbank2gff3_positions_comulative(self, genbank, out_file):
        """ """
        temp_file = self.utils.get_temp_file("gbk_to_gff3", ".txt")

        ## set VERSION ID equal to ACCESSION
        out_file_gb = self.utils.get_temp_file("file_name", ".gb")
        self.utils.clean_genbank_version_name(genbank, out_file_gb)

        with open(temp_file, "w") as out_handle, open(out_file_gb) as in_handle:
            GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)

        ### remove clean version gb file
        os.unlink(out_file_gb)
        ### filter file
        #     vect_filter = ['remark', 'source', 'gene']
        vect_pass = ["CDS"]
        b_fail = False
        dt_data_sequence_length = {}
        sequence_name = ""
        actual_position = 0
        with open(temp_file) as handle, open(out_file, "w") as handle_write:
            for line in handle:
                sz_temp = line.strip()
                if (
                    len(sz_temp) == 0
                    or sz_temp.find("# Input") == 0
                    or sz_temp.find("# GFF3 saved") == 0
                ):
                    continue
                if sz_temp.find("##FASTA") == 0:
                    break

                lst_data = sz_temp.split("\t")
                if sz_temp[0] == "#":
                    lst_data = sz_temp.split()
                    if sz_temp.startswith("##sequence-region") and len(lst_data) == 4:
                        dt_data_sequence_length[lst_data[1]] = int(lst_data[3])
                    handle_write.write(sz_temp + "\n")
                elif len(lst_data) > 3 and lst_data[2] in vect_pass:
                    if lst_data[0] in dt_data_sequence_length:
                        if len(sequence_name) == 0:
                            sequence_name = lst_data[0]
                            actual_position = 0
                        elif sequence_name != lst_data[0]:
                            actual_position += dt_data_sequence_length[sequence_name]
                            sequence_name = lst_data[0]
                    else:
                        b_fail = True
                        break
                    lst_data[3] = str(int(lst_data[3]) + actual_position)
                    lst_data[4] = str(int(lst_data[4]) + actual_position)

                    vect_out = [
                        data
                        for data in lst_data[8].split(";")
                        if not data.startswith("translation=")
                    ]
                    lst_data[8] = ";".join(vect_out)
                    handle_write.write("\t".join(lst_data) + "\n")
        os.unlink(temp_file)
        return None if b_fail else out_file

    def get_snpeff_config(self, fasta_file):
        """
        return a file for snpEff config for the fasta
        genome_name: this name is going to set in properties in the file
        """
        temp_file = self.utils.get_temp_file("snpEff_config", ".config")
        if not os.path.exists(self.software_names.get_snp_eff_config()):
            raise IOError(
                "Error: file not found {}".format(
                    self.software_names.get_snp_eff_config()
                )
            )

        self.utils.copy_file(self.software_names.get_snp_eff_config(), temp_file)
        handle = open(temp_file, "a")
        base_file_name = os.path.basename(fasta_file)
        base_file_name = base_file_name[0 : base_file_name.rfind(".")]
        handle.write(
            "{}.genome : {} reference".format(base_file_name, Constants.INSAFLU_NAME)
        )

        ### open fasta
        if self.utils.is_gzip(fasta_file):
            handle_fasta = gzip.open(fasta_file, mode="rt")
        else:
            handle_fasta = open(fasta_file)
        vect_names = []
        for record in SeqIO.parse(handle_fasta, Constants.FORMAT_FASTA):
            vect_names.append(record.name)
        handle_fasta.close()

        handle.write(
            "\n{}.chromosome : {}".format(base_file_name, ", ".join(vect_names))
        )
        for chromosome in vect_names:
            handle.write(
                "\n{}.{}.codonTable : Bacterial_and_Plant_Plastid".format(
                    base_file_name, chromosome
                )
            )
        handle.close()
        return (base_file_name, temp_file)

    def run_snpEff(self, fasta_file, genbank, vcf_file, out_file):
        """
        ./snpEff ann -no-downstream -no-upstream -no-intergenic -no-utr -c ../path_to/reference/snpeff.config -dataDir . -noStats ref sample.vcf > sample_annot.vcf
        """
        temp_dir = self.utils.get_temp_dir()
        (fasta_file_name, snpeff_config) = self.get_snpeff_config(fasta_file)

        temp_vcf_file = os.path.join(temp_dir, os.path.basename(vcf_file))
        self.utils.copy_file(vcf_file, temp_vcf_file)

        ## create the database
        out_gff_file = self.utils.get_temp_file("temp_gff", ".gff")
        self.run_genbank2gff3(genbank, out_gff_file)

        ### count sequences, if none return None
        if self.utils.get_number_sequeces_in_gff_file(out_gff_file) == 0:
            os.unlink(out_gff_file)
            os.unlink(snpeff_config)
            self.utils.remove_dir(temp_dir)
            return None

        datase_dir = "{}".format(fasta_file_name)
        cmd = "mkdir -p {}".format(os.path.join(temp_dir, datase_dir))
        os.system(cmd)
        cmd = "mkdir -p {}".format(os.path.join(temp_dir, "genomes"))
        os.system(cmd)
        self.utils.copy_file(
            out_gff_file, os.path.join(temp_dir, datase_dir, "genes.gff")
        )
        temp_file = os.path.join(temp_dir, "genomes", fasta_file_name + ".fa")
        self.utils.copy_file(fasta_file, temp_file)
        os.unlink(out_gff_file)

        ## indexing database
        ## snpEff build -c reference/snpeff.config -dataDir . -gff3 ref 2>> run_snippy2_1.log
        cmd = "%s build -c %s -dataDir %s -gff3 %s" % (
            self.software_names.get_snp_eff(),
            snpeff_config,
            temp_dir,
            fasta_file_name,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            os.unlink(snpeff_config)
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to create snpEff database")

        ### create the annotation
        cmd = "%s ann %s -c %s -dataDir %s %s %s > %s" % (
            self.software_names.get_snp_eff(),
            self.software_names.get_snp_eff_parameters(),
            snpeff_config,
            temp_dir,
            fasta_file_name,
            temp_vcf_file,
            out_file,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            os.unlink(snpeff_config)
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run snpEff")

        #### add the transform p.Val423Glu to p.V423G
        parse_out_files = ParseOutFiles()
        out_file_transformed_amino = parse_out_files.add_amino_single_letter_code(
            out_file
        )
        self.utils.move_file(out_file_transformed_amino, out_file)

        os.unlink(snpeff_config)
        self.utils.remove_dir(temp_dir)
        return out_file

    def run_snippy_vcf_to_tab(self, fasta, genbank, vcf_file, out_file):
        """
        ./snippy-vcf_to_tab_add_freq [options] --ref ref.fa [--gff ref.gff] --vcf snps.vcf > snp.tab
        """

        temp_file = self.utils.get_temp_file("snippy_vcf_to_tab", ".gff")
        self.run_genbank2gff3(genbank, temp_file)

        cmd = "%s --ref %s --gff %s --vcf %s > %s" % (
            self.software_names.get_snippy_vcf_to_tab(),
            fasta,
            temp_file,
            vcf_file,
            out_file,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            os.unlink(temp_file)
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run snippy-vcf-to-tab")
        os.unlink(temp_file)
        return out_file

    def run_snippy_vcf_to_tab_freq_and_evidence(
        self, fasta, genbank, vcf_file, out_file
    ):
        """
        ./snippy-vcf_to_tab_add_freq_and_evidence [options] --ref ref.fa [--gff ref.gff] --vcf snps.vcf > snp.tab
        """

        temp_file = self.utils.get_temp_file("snippy_vcf_to_tab", ".gff")
        self.run_genbank2gff3(genbank, temp_file)

        cmd = "%s --ref %s --gff %s --vcf %s > %s" % (
            self.software_names.get_snippy_vcf_to_tab_freq_and_evidence(),
            fasta,
            temp_file,
            vcf_file,
            out_file,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            os.unlink(temp_file)
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run snippy-vcf-to-tab. Add freq and evidence")
        os.unlink(temp_file)
        return out_file

    def run_freebayes(self, bam_file, reference_fasta, genbank_file, sample_name):
        """
        run freebayes
        return output directory

        ## freebayes -p 1 -q 20 -m 20 --min-coverage 100 --min-alternate-fraction 0.01 --min-alternate-count 10 -V -f ../$input2 -b {} > {.}.vcf'
        """
        temp_dir = os.path.join(self.utils.get_temp_dir())

        file_to_process = os.path.join(temp_dir, sample_name + ".bam")
        cmd = "ln -s {} {}".format(
            bam_file, os.path.join(temp_dir, sample_name + ".bam")
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run freebayes")

        cmd = "ln -s {} {}.bai".format(bam_file + ".bai", file_to_process)
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run freebayes")

        reference_fasta_temp = os.path.join(temp_dir, os.path.basename(reference_fasta))
        cmd = "ln -s {} {}".format(reference_fasta, reference_fasta_temp)
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run freebayes")

        reference_fasta_fai = reference_fasta + FileExtensions.FILE_FAI
        if not os.path.exists(reference_fasta_fai):
            self.logger_production.error("Files doesnt exist: " + reference_fasta_fai)
            self.logger_debug.error("Files doesnt exist: " + reference_fasta_fai)
            raise Exception("Fail to run freebayes")

        reference_fasta_temp_fai = os.path.join(
            temp_dir, os.path.basename(reference_fasta) + FileExtensions.FILE_FAI
        )
        cmd = "ln -s {} {}".format(reference_fasta_fai, reference_fasta_temp_fai)
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run freebayes")

        temp_file = self.utils.get_temp_file("freebayes_temp", ".vcf")
        cmd = "%s %s -f %s -b %s > %s" % (
            self.software_names.get_freebayes(),
            self.software_names.get_freebayes_parameters(),
            reference_fasta_temp,
            file_to_process,
            temp_file,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run freebayes")

        if os.path.getsize(temp_file) > 0:
            ### run snpEff
            temp_file_2 = self.utils.get_temp_file("vcf_file", ".vcf")
            output_file = self.run_snpEff(
                reference_fasta,
                genbank_file,
                temp_file,
                os.path.join(temp_dir, os.path.basename(temp_file_2)),
            )

            if output_file is None:  ## sometimes the gff does not have amino sequences
                self.utils.copy_file(
                    temp_file, os.path.join(temp_dir, os.path.basename(temp_file_2))
                )

            self.test_bgzip_and_tbi_in_vcf(
                os.path.join(temp_dir, os.path.basename(temp_file_2))
            )

            ### add FREQ to vcf file
            vcf_file_out_temp = self.utils.add_freq_to_vcf(
                os.path.join(temp_dir, os.path.basename(temp_file_2)),
                os.path.join(temp_dir, sample_name + ".vcf"),
            )
            os.unlink(temp_file)
            if os.path.exists(temp_file_2):
                os.unlink(temp_file_2)

            ### pass vcf to tab
            self.run_snippy_vcf_to_tab(
                reference_fasta,
                genbank_file,
                vcf_file_out_temp,
                "{}.tab".format(os.path.join(temp_dir, sample_name)),
            )
        return temp_dir

    def run_freebayes_parallel(
        self, bam_file, reference_fasta, genbank_file, sample_name
    ):
        """
        run freebayes
        return output directory

        ## freebayes -p 1 -q 20 -m 20 --min-coverage 100 --min-alternate-fraction 0.01 --min-alternate-count 10 -V -f ../$input2 -b {} > {.}.vcf'
        """
        temp_dir = os.path.join(self.utils.get_temp_dir())

        file_to_process = os.path.join(temp_dir, sample_name + ".bam")
        cmd = "ln -s {} {}".format(
            bam_file, os.path.join(temp_dir, sample_name + ".bam")
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run freebayes parallel")

        cmd = "ln -s {} {}.bai".format(bam_file + ".bai", file_to_process)
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run freebayes parallel")

        reference_fasta_temp = os.path.join(temp_dir, os.path.basename(reference_fasta))
        cmd = "ln -s {} {}".format(reference_fasta, reference_fasta_temp)
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run freebayes parallel")

        reference_fasta_fai = reference_fasta + FileExtensions.FILE_FAI
        if not os.path.exists(reference_fasta_fai):
            self.logger_production.error("Files doesnt exist: " + reference_fasta_fai)
            self.logger_debug.error("Files doesnt exist: " + reference_fasta_fai)
            raise Exception("Fail to run freebayes parallel")

        reference_fasta_temp_fai = os.path.join(
            temp_dir, os.path.basename(reference_fasta) + FileExtensions.FILE_FAI
        )
        cmd = "ln -s {} {}".format(reference_fasta_fai, reference_fasta_temp_fai)
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run freebayes parallel")

        ### create regions
        temp_file_bam_coverage = self.utils.get_temp_file("bamtools_coverage", ".txt")

        # bamtools coverage -in aln.bam | coverage_to_regions.py ref.fa 500 >ref.fa.500.regions
        cmd = "%s coverage -in %s > %s" % (
            self.software_names.get_bamtools(),
            bam_file,
            temp_file_bam_coverage,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            os.unlink(temp_file_bam_coverage)
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run bamtools coverage")

        #### test if it has some data, if all zero does not have data
        if (os.path.getsize(temp_file_bam_coverage)) == 0:
            if os.path.exists(temp_file_bam_coverage):
                os.unlink(temp_file_bam_coverage)
            self.utils.remove_dir(temp_dir)
            return None

        temp_file_regions = self.utils.get_temp_file("freebayes_regions", ".txt")
        cmd = "cat %s | %s %s 500 > %s" % (
            temp_file_bam_coverage,
            self.software_names.get_coverage_to_regions(),
            reference_fasta_temp_fai,
            temp_file_regions,
        )

        exist_status = os.system(cmd)
        if exist_status != 0:
            os.unlink(temp_file_bam_coverage)
            os.unlink(temp_file_regions)
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run freebayes parallel")

        temp_file = self.utils.get_temp_file("freebayes_temp", ".vcf")
        cmd = "%s %s %s %s -f %s -b %s > %s" % (
            self.software_names.get_freebayes_parallel(),
            temp_file_regions,
            settings.THREADS_TO_RUN_SLOW,
            self.software_names.get_freebayes_parameters(),
            reference_fasta_temp,
            file_to_process,
            temp_file,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            os.unlink(temp_file_regions)
            os.unlink(temp_file_bam_coverage)
            os.unlink(temp_file)
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run freebayes parallel")

        ### run snpEff
        if os.path.exists(temp_file):
            temp_file_2 = self.utils.get_temp_file("vcf_file", ".vcf")
            output_file = self.run_snpEff(
                reference_fasta,
                genbank_file,
                temp_file,
                os.path.join(temp_dir, os.path.basename(temp_file_2)),
            )
            if output_file is None:  ## sometimes the gff does not have amino sequences
                self.utils.copy_file(
                    temp_file, os.path.join(temp_dir, os.path.basename(temp_file_2))
                )

            self.test_bgzip_and_tbi_in_vcf(
                os.path.join(temp_dir, os.path.basename(temp_file_2))
            )

            ### add FREQ to vcf file
            vcf_file_out_temp = self.utils.add_freq_to_vcf(
                os.path.join(temp_dir, os.path.basename(temp_file_2)),
                os.path.join(temp_dir, sample_name + ".vcf"),
            )
            ### pass vcf to tab
            self.run_snippy_vcf_to_tab(
                reference_fasta,
                genbank_file,
                vcf_file_out_temp,
                "{}.tab".format(os.path.join(temp_dir, sample_name)),
            )
        if os.path.exists(temp_file):
            os.unlink(temp_file)
        if os.path.exists(temp_file_2):
            os.unlink(temp_file_2)
        if os.path.exists(temp_file_regions):
            os.unlink(temp_file_regions)
        if os.path.exists(temp_file_bam_coverage):
            os.unlink(temp_file_bam_coverage)
        return temp_dir

    """
    Global processing
    """
    #     @transaction.atomic
    def run_fastq_and_trimmomatic(self, sample, owner):
        """
        run fastq and trimmomatic
        Upload average and sequence numbers
        :output (Has data?, Is Trimmomatic run?)
        """
        manage_database = ManageDatabase()
        result_all = Result()

        ### first try run downsize if necessary
        if settings.DOWN_SIZE_FASTQ_FILES:
            (is_downsized, file_name_1, file_name_2) = self.make_downsize(
                sample.get_fastq(TypePath.MEDIA_ROOT, True),
                sample.get_fastq(TypePath.MEDIA_ROOT, False),
                settings.MAX_FASTQ_FILE_UPLOAD,
            )
            if is_downsized:
                if os.path.exists(file_name_1) and os.path.getsize(file_name_1) > 100:
                    self.utils.move_file(
                        file_name_1, sample.get_fastq(TypePath.MEDIA_ROOT, True)
                    )
                if (
                    file_name_2 != None
                    and len(file_name_2) > 0
                    and os.path.exists(file_name_2)
                    and os.path.getsize(file_name_2) > 100
                ):
                    self.utils.move_file(
                        file_name_2, sample.get_fastq(TypePath.MEDIA_ROOT, False)
                    )

                ### set the downsize message
                manage_database.set_sample_metakey(
                    sample,
                    owner,
                    MetaKeyAndValue.META_KEY_ALERT_DOWNSIZE_OF_FASTQ_FILES,
                    MetaKeyAndValue.META_VALUE_Success,
                    "Fastq files were down sized to values ~{}.".format(
                        filesizeformat(int(settings.MAX_FASTQ_FILE_UPLOAD))
                    ),
                )

        ### first run fastqc
        try:
            temp_dir = self.run_fastq(
                sample.get_fastq(TypePath.MEDIA_ROOT, True),
                sample.get_fastq(TypePath.MEDIA_ROOT, False),
            )
            result_all.add_software(
                SoftwareDesc(
                    self.software_names.get_fastq_name(),
                    self.software_names.get_fastq_version(),
                    self.software_names.get_fastq_parameters(),
                )
            )

            ### need to copy the files to samples/user path
            self.utils.copy_file(
                os.path.join(
                    temp_dir,
                    os.path.basename(
                        sample.get_fastqc_output(TypePath.MEDIA_ROOT, True)
                    ),
                ),
                sample.get_fastqc_output(TypePath.MEDIA_ROOT, True),
            )
            if sample.exist_file_2():
                self.utils.copy_file(
                    os.path.join(
                        temp_dir,
                        os.path.basename(
                            sample.get_fastqc_output(TypePath.MEDIA_ROOT, False)
                        ),
                    ),
                    sample.get_fastqc_output(TypePath.MEDIA_ROOT, False),
                )
            self.utils.remove_dir(temp_dir)
        except Exception as e:
            result = Result()
            result.set_error("Fail to run fastq software: " + e.args[0])
            result.add_software(
                SoftwareDesc(
                    self.software_names.get_fastq_name(),
                    self.software_names.get_fastq(),
                    self.software_names.get_fastq_parameters(),
                )
            )
            manage_database.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_Fastq_Trimmomatic,
                MetaKeyAndValue.META_VALUE_Error,
                result.to_json(),
            )
            return (False, False)

        ### get number of sequences and average length
        (number_1, average_1, std1) = self.get_number_sequences_in_fastq(
            sample.get_fastq(TypePath.MEDIA_ROOT, True)
        )
        if sample.exist_file_2():
            (number_2, average_2, std_2) = self.get_number_sequences_in_fastq(
                sample.get_fastq(TypePath.MEDIA_ROOT, False)
            )
        else:
            (number_2, average_2, std_2) = (None, None, None)

        ### create key/values for stat in Illumina, before trimmomatic
        result_all.add_software(
            SoftwareDesc(
                SoftwareNames.SOFTWARE_ILLUMINA_stat,
                "",
                "",
                self.get_key_values_stats_illumina(
                    number_1, average_1, std1, number_2, average_2, std_2
                ),
            )
        )

        ### get software parameters
        default_software_project = DefaultProjectSoftware()
        default_software_project.test_all_defaults(sample.owner, None, None, sample)

        ### test if the software ran
        if not default_software_project.is_to_run_trimmomatic(sample.owner, sample):
            ## set software run
            manage_database.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_Fastq_Trimmomatic,
                MetaKeyAndValue.META_VALUE_Success,
                "Success, Fastq({})".format(self.software_names.get_fastq_version()),
            )
            manage_database.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software,
                MetaKeyAndValue.META_VALUE_Success,
                result_all.to_json(),
            )

            ### set number of reads to show on the table
            result_average = ResultAverageAndNumberReads(
                number_1, average_1, number_2, average_2
            )
            manage_database.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_Number_And_Average_Reads,
                MetaKeyAndValue.META_VALUE_Success,
                result_average.to_json(),
            )

            ## remove trimmomatic result files
            self.utils.remove_file(
                sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)
            )
            self.utils.remove_file(
                sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False)
            )
            self.utils.remove_file(
                sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)
            )
            self.utils.remove_file(
                sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)
            )
            return (result_average.has_reads(), False)

        ### run trimmomatic
        try:
            (temp_dir, filtering_result, parameters) = self.run_trimmomatic(
                sample.get_fastq(TypePath.MEDIA_ROOT, True),
                sample.get_fastq(TypePath.MEDIA_ROOT, False),
                sample,
                owner,
            )

            ### it run
            result_all.add_software(
                SoftwareDesc(
                    self.software_names.get_trimmomatic_name(),
                    self.software_names.get_trimmomatic_version(),
                    parameters,
                    filtering_result.key_values,
                )
            )

            ### need to copy the files to samples/user path
            self.utils.copy_file(
                os.path.join(
                    temp_dir,
                    os.path.basename(
                        sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)
                    ),
                ),
                sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True),
            )
            if sample.exist_file_2():
                self.utils.copy_file(
                    os.path.join(
                        temp_dir,
                        os.path.basename(
                            sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False)
                        ),
                    ),
                    sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False),
                )
            self.utils.remove_dir(temp_dir)
        except Exception as e:
            result = Result()
            result.set_error("Fail to run trimmomatic software: " + e.args[0])
            result.add_software(
                SoftwareDesc(
                    self.software_names.get_trimmomatic_name(),
                    self.software_names.get_trimmomatic_version(),
                    self.software_names.get_trimmomatic_parameters(),
                )
            )
            manage_database.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_Fastq_Trimmomatic,
                MetaKeyAndValue.META_VALUE_Error,
                result.to_json(),
            )
            self.utils.remove_dir(temp_dir)
            return (False, False)

        ### run fastq again
        try:
            temp_dir = self.run_fastq(
                sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True),
                sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False),
            )

            ### need to copy the files to samples/user path
            self.utils.copy_file(
                os.path.join(
                    temp_dir,
                    os.path.basename(
                        sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)
                    ),
                ),
                sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True),
            )
            if sample.exist_file_2():
                self.utils.copy_file(
                    os.path.join(
                        temp_dir,
                        os.path.basename(
                            sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)
                        ),
                    ),
                    sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False),
                )
            self.utils.remove_dir(temp_dir)
        except Exception as e:
            result = Result()
            result.set_error("Fail to run fastq software: " + e.args[0])
            result.add_software(
                SoftwareDesc(
                    self.software_names.get_fastq_name(),
                    self.software_names.get_fastq(),
                    self.software_names.get_fastq_parameters(),
                )
            )
            manage_database.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_Fastq_Trimmomatic,
                MetaKeyAndValue.META_VALUE_Error,
                result.to_json(),
            )
            return (False, False)

        ### collect numbers
        (number_1, average_1, std1) = self.get_number_sequences_in_fastq(
            sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)
        )
        if sample.exist_file_2():
            (number_2, average_2, std2) = self.get_number_sequences_in_fastq(
                sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False)
            )
        else:
            (number_2, average_2, std2) = (None, None, None)

        result_average = ResultAverageAndNumberReads(
            number_1, average_1, number_2, average_2
        )
        manage_database.set_sample_metakey(
            sample,
            owner,
            MetaKeyAndValue.META_KEY_Number_And_Average_Reads,
            MetaKeyAndValue.META_VALUE_Success,
            result_average.to_json(),
        )

        ### create key/values for stat in Illumina, after trimmomatic
        result_all.add_software(
            SoftwareDesc(
                SoftwareNames.SOFTWARE_ILLUMINA_stat,
                "",
                "",
                self.get_key_values_stats_illumina(
                    number_1, average_1, std1, number_2, average_2, std2
                ),
            )
        )

        ## save everything OK
        manage_database.set_sample_metakey(
            sample,
            owner,
            MetaKeyAndValue.META_KEY_Fastq_Trimmomatic,
            MetaKeyAndValue.META_VALUE_Success,
            "Success, Fastq(%s), Trimmomatic(%s)"
            % (
                self.software_names.get_fastq_version(),
                self.software_names.get_trimmomatic_version(),
            ),
        )
        manage_database.set_sample_metakey(
            sample,
            owner,
            MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software,
            MetaKeyAndValue.META_VALUE_Success,
            result_all.to_json(),
        )
        return result_average.has_reads(), True

    def get_number_sequences_in_fastq(self, fastq_file):
        """return number of sequences, average and std"""
        (number_sequences, average_length, std) = (0, 0, 0)
        try:
            (
                number_sequences,
                average_length,
                std,
            ) = self.utils.get_number_sequences_fastq(fastq_file)
        except Exception as e:
            pass
        return (number_sequences, average_length, std)

    def get_key_values_stats_illumina(
        self, number_1, average_1, std1, number_2, average_2, std2
    ):
        """
        return key values for illumina
        Check the keys in the SoftwareNames.SOFTWARE_ILLUMINA_stat_collect
                "Number of reads R1", "Number of reads R2",
                "Average read length R1", "Average read length R2",
                "STDEV read length R1", "STDEV read length R2",
        """
        result = Result()
        total_reads = 0
        if not number_1 is None:
            result.add_key_value(
                KeyValue("Number of reads R1", "{:,.0f}".format(number_1))
            )
            total_reads += number_1
        if not average_1 is None:
            result.add_key_value(KeyValue("Average read length R1", str(average_1)))
        if not std1 is None:
            result.add_key_value(KeyValue("STDEV read length R1", str(std1)))
        if not number_2 is None:
            result.add_key_value(
                KeyValue("Number of reads R2", "{:,.0f}".format(number_2))
            )
            total_reads += number_2
        if not average_2 is None:
            result.add_key_value(KeyValue("Average read length R2", str(average_2)))
        if not std2 is None:
            result.add_key_value(KeyValue("STDEV read length R2", str(std2)))
        result.add_key_value(KeyValue("Total reads", "{:,.0f}".format(total_reads)))
        return result.key_values

    def get_stats_from_sample_reads(self, sample):
        """
        Return stats to show on sample details
        :out data_stat -> [[key name, value before filtering, value before filtering, percentage],
                                [key name, value before filtering, value before filtering, percentage],
                                [key name, value before filtering, value before filtering, percentage],... ]
        :out total_reads_start -> number of reads before mapping
        """
        if sample.is_type_fastq_gz_sequencing():
            meta_key_software = MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software
            software_name = SoftwareNames.SOFTWARE_ILLUMINA_stat
            show_percentage = (
                SoftwareNames.SOFTWARE_ILLUMINA_stat_collect_show_percentage
            )
        else:
            meta_key_software = MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software
            software_name = SoftwareNames.SOFTWARE_NanoStat_name
            show_percentage = (
                SoftwareNames.SOFTWARE_NANOSTAT_vect_info_to_collect_show_percentage
            )

        ## return variables
        data_stat, total_reads_start = [], -1
        decode_nanostat = DecodeObjects()
        manageDatabase = ManageDatabase()
        meta_data = manageDatabase.get_sample_metakey_last(
            sample, meta_key_software, None
        )
        if meta_data is None:
            return data_stat, total_reads_start

        result_data = decode_nanostat.decode_result(meta_data.description)
        vect_soft = result_data.get_list_software_instance(software_name)
        if len(vect_soft) == 2:  ## has data
            key_data_0 = vect_soft[0].get_vect_key_values()
            key_data_1 = vect_soft[1].get_vect_key_values()
            for _ in range(len(key_data_0)):
                percentage = "--"
                value_0 = key_data_0[_].value.replace(",", "")
                value_1 = key_data_1[_].value.replace(",", "")
                if key_data_0[_].key in show_percentage and (
                    len(value_0) > 0
                    and self.utils.is_float(value_0)
                    and float(value_0) > 0
                    and len(value_1) > 0
                ):
                    percentage = "{:,.1f}".format(float(value_1) / float(value_0) * 100)
                data_stat.append(
                    [
                        key_data_0[_].key,
                        key_data_0[_].value,
                        key_data_1[_].value,
                        "{:,.1f}".format(float(value_1) - float(value_0)),
                        percentage,
                    ]
                )

                ### get total reads
                if key_data_1[_].key in (
                    "Number of reads R1",
                    "Number of reads R2",
                    "Number of reads",
                ) and self.utils.is_float(value_1):
                    if total_reads_start == -1:
                        total_reads_start = 0
                    total_reads_start += int(float(value_1))

        elif len(vect_soft) == 1:  ## only has the fist analysis
            key_data_0 = vect_soft[0].get_vect_key_values()
            for _ in range(len(key_data_0)):
                percentage = "--"
                value_0 = key_data_0[_].value.replace(",", "")
                if key_data_0[_].key in show_percentage and (
                    len(value_0) > 0
                    and self.utils.is_float(value_0)
                    and float(value_0) > 0
                ):
                    percentage = "100"
                data_stat.append(
                    [key_data_0[_].key, key_data_0[_].value, "--", "0.0", percentage]
                )

                ### get total reads
                if key_data_0[_].key in (
                    "Number of reads",
                    "Number of reads R1",
                    "Number of reads R2",
                ) and self.utils.is_float(value_0):
                    if total_reads_start == -1:
                        total_reads_start = 0
                    total_reads_start += int(float(value_0))
        return data_stat, None if total_reads_start == -1 else total_reads_start

    """
    Global processing, fastQ, trimmomatic and GetSpecies
    """

    def run_fastq_and_trimmomatic_and_identify_species(self, sample, user):

        print("Start ProcessControler")
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        process_SGE.set_process_controler(
            user,
            process_controler.get_name_sample(sample),
            ProcessControler.FLAG_RUNNING,
        )

        ### it can be deleted
        if sample.is_deleted or not sample.is_valid_1:
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_sample(sample),
                ProcessControler.FLAG_FINISHED,
            )
            return True

        ################################
        ##################################
        ### remove possible previous alerts from others run
        manage_database = ManageDatabase()
        for keys_to_remove in MetaKeyAndValue.VECT_TO_REMOVE_RUN_SAMPLE:
            manage_database.remove_sample_start_metakey(sample, keys_to_remove)
        
        ### remove some other 
        sample.identify_virus.all().delete()
        if (not sample.mixed_infections_tag is None): sample.mixed_infections_tag = None
        sample.number_alerts = 0
        sample.save()
        
        try:
            print("Start run_fastq_and_trimmomatic")
            ### run trimmomatics
            b_has_data, b_it_ran = self.run_fastq_and_trimmomatic(sample, user)
            
            print("Result run_fastq_and_trimmomatic: " + str(b_has_data))
            
            ### test Abricate ON/OFF
            default_software_project = DefaultProjectSoftware()
            b_make_identify_species = default_software_project.is_to_run_abricate(sample.owner, sample,
                                            ConstantsSettings.TECHNOLOGY_illumina)
            
            ### queue the quality check and
            if (b_has_data and b_make_identify_species):    ## don't run for single file because spades doesn't work for one single file
                self.identify_type_and_sub_type(sample, sample.get_fastq_available(TypePath.MEDIA_ROOT, True),\
                    sample.get_fastq_available(TypePath.MEDIA_ROOT, False), user)
    
            ## set the flag that is ready for process
            sample_to_update = Sample.objects.get(pk=sample.id)
            sample_to_update.is_sample_in_the_queue = False
            if (b_has_data):
                sample_to_update.is_ready_for_projects = True
                
                ### make identify species
                if (b_make_identify_species):
                    sample_to_update.type_subtype = sample_to_update.get_type_sub_type()[:Sample.TYPE_SUBTYPE_LENGTH-1]
                    
                    (tag_mixed_infection, alert, message) = sample_to_update.get_mixed_infection()
                    if (sample_to_update.number_alerts == None): sample_to_update.number_alerts = alert
                    else: sample_to_update.number_alerts += alert
                        
                    manage_database = ManageDatabase()
                    if (message != None and len(message) > 0):
                        manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE,\
                                            MetaKeyAndValue.META_VALUE_Success, message)
        
                    ### save tag mixed_infecion
                    manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_TAG_MIXED_INFECTION_TYPE_SUBTYPE,\
                                MetaKeyAndValue.META_VALUE_Success, tag_mixed_infection)
    
                    try:
                        mixed_infections_tag = MixedInfectionsTag.objects.get(name=tag_mixed_infection)
                    except MixedInfectionsTag.DoesNotExist as e:
                        mixed_infections_tag = MixedInfectionsTag()
                        mixed_infections_tag.name = tag_mixed_infection
                        mixed_infections_tag.save()
                    
                    sample_to_update.mixed_infections_tag = mixed_infections_tag
                else:
                    sample_to_update.type_subtype = Constants.EMPTY_VALUE_NA
                    tag_mixed_infection = Constants.EMPTY_VALUE_NA
                    try:
                        mixed_infections_tag = MixedInfectionsTag.objects.get(name=tag_mixed_infection)
                    except MixedInfectionsTag.DoesNotExist as e:
                        mixed_infections_tag = MixedInfectionsTag()
                        mixed_infections_tag.name = tag_mixed_infection
                        mixed_infections_tag.save()
                    
                    sample_to_update.mixed_infections_tag = mixed_infections_tag
                    
                    manage_database = ManageDatabase()
                    message = "Info: Abricate turned OFF by the user."
                    manage_database.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE,\
                                MetaKeyAndValue.META_VALUE_Success, message)
            else:
                manage_database = ManageDatabase()
                manage_database.set_sample_metakey(sample_to_update, user, MetaKeyAndValue.META_KEY_ALERT_NO_READS_AFTER_FILTERING,\
                                        MetaKeyAndValue.META_VALUE_Success, "Warning: no reads left after filtering.")
                
                if (sample_to_update.number_alerts == None): sample_to_update.number_alerts = 1
                else: sample_to_update.number_alerts += 1
                sample_to_update.is_ready_for_projects = False
                sample_to_update.type_subtype = Constants.EMPTY_VALUE_TYPE_SUBTYPE
            sample_to_update.save()
            
            ### set the flag of the end of the task        
            meta_sample = manage_database.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue)
            if (meta_sample != None):
                manage_database.set_sample_metakey(sample, sample.owner, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Success, meta_sample.description)

        except:
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_sample(sample),
                ProcessControler.FLAG_ERROR,
            )
            return False

        ### finished
        process_SGE.set_process_controler(
            user,
            process_controler.get_name_sample(sample),
            ProcessControler.FLAG_FINISHED,
        )
        return b_has_data

    def process_second_stage_snippy_coverage_freebayes(self, project_sample, user):
        """ """
        ### make it running
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        process_SGE.set_process_controler(
            user,
            process_controler.get_name_project_sample(project_sample),
            ProcessControler.FLAG_RUNNING,
        )

        ## run collect data
        return self.__process_second_stage_snippy_coverage_freebayes(
            project_sample, user
        )

    """
    Global processing, Snippy, Coverage, Freebayes and MixedInfections
    """
    #     @transaction.atomic
    def __process_second_stage_snippy_coverage_freebayes(self, project_sample, user):
        """
        Global processing, snippy, coverage,
        """

        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        manageDatabase = ManageDatabase()
        result_all = Result()
        ### metakey for this process
        metaKeyAndValue = MetaKeyAndValue()
        try:
            meta_key_project_sample = (
                metaKeyAndValue.get_meta_key_queue_by_project_sample_id(
                    project_sample.id
                )
            )

            ### Test if this sample already run
            meta_sample = manageDatabase.get_project_sample_metakey_last(
                project_sample,
                meta_key_project_sample,
                MetaKeyAndValue.META_VALUE_Queue,
            )
            if (
                meta_sample != None
                and meta_sample.value == MetaKeyAndValue.META_VALUE_Success
            ):
                return

            ## test software parameters for project_sample
            default_project_software = DefaultProjectSoftware()
            default_project_software.test_all_defaults(user, None, project_sample, None)

            ## process snippy
            try:
                ### get snippy parameters
                snippy_parameters = (
                    default_project_software.get_snippy_parameters_all_possibilities(
                        user, project_sample
                    )
                )
                out_put_path = self.run_snippy(
                    project_sample.sample.get_fastq_available(
                        TypePath.MEDIA_ROOT, True
                    ),
                    project_sample.sample.get_fastq_available(
                        TypePath.MEDIA_ROOT, False
                    ),
                    project_sample.project.reference.get_reference_fasta(
                        TypePath.MEDIA_ROOT
                    ),
                    project_sample.project.reference.get_reference_gbk(
                        TypePath.MEDIA_ROOT
                    ),
                    project_sample.sample.name,
                    snippy_parameters,
                )
                result_all.add_software(
                    SoftwareDesc(
                        self.software_names.get_snippy_name(),
                        self.software_names.get_snippy_version(),
                        snippy_parameters,
                    )
                )
            except Exception as e:
                result = Result()
                result.set_error(e.args[0])
                result.add_software(
                    SoftwareDesc(
                        self.software_names.get_snippy_name(),
                        self.software_names.get_snippy_version(),
                        snippy_parameters,
                    )
                )
                manageDatabase.set_project_sample_metakey(
                    project_sample,
                    user,
                    MetaKeyAndValue.META_KEY_Snippy,
                    MetaKeyAndValue.META_VALUE_Error,
                    result.to_json(),
                )

                ### get again and set error
                project_sample = ProjectSample.objects.get(pk=project_sample.id)
                project_sample.is_error = True
                project_sample.save()

                meta_sample = manageDatabase.get_project_sample_metakey_last(
                    project_sample,
                    meta_key_project_sample,
                    MetaKeyAndValue.META_VALUE_Queue,
                )
                if meta_sample != None:
                    manageDatabase.set_project_sample_metakey(
                        project_sample,
                        user,
                        meta_key_project_sample,
                        MetaKeyAndValue.META_VALUE_Error,
                        meta_sample.description,
                    )
                process_SGE.set_process_controler(
                    user,
                    process_controler.get_name_project_sample(project_sample),
                    ProcessControler.FLAG_ERROR,
                )
                return False

            ## copy the files to the project sample directories
            self.copy_files_to_project(
                project_sample, self.software_names.get_snippy_name(), out_put_path
            )
            self.utils.remove_dir(out_put_path)

            ### make the link for the new tab file name
            path_snippy_tab = project_sample.get_file_output(
                TypePath.MEDIA_ROOT,
                FileType.FILE_TAB,
                self.software_names.get_snippy_name(),
            )
            if os.path.exists(path_snippy_tab):
                sz_file_to = project_sample.get_file_output_human(
                    TypePath.MEDIA_ROOT,
                    FileType.FILE_TAB,
                    self.software_names.get_snippy_name(),
                )
                self.utils.link_file(path_snippy_tab, sz_file_to)

            ### get mapped stast reads
            bam_file = project_sample.get_file_output(
                TypePath.MEDIA_ROOT,
                FileType.FILE_BAM,
                self.software_names.get_snippy_name(),
            )
            result = Result()
            if os.path.exists(bam_file):
                result = self.get_statistics_bam(bam_file)
            manageDatabase.set_project_sample_metakey(
                project_sample,
                user,
                MetaKeyAndValue.META_KEY_bam_stats,
                MetaKeyAndValue.META_VALUE_Success,
                result.to_json(),
            )

            ## get coverage from deep file
            get_coverage = GetCoverage()
            try:

                ### limit of the coverage for a project, can be None, if not exist
                coverage_for_project = (
                    default_project_software.get_snippy_single_parameter_for_project(
                        project_sample.project, DefaultParameters.SNIPPY_COVERAGE_NAME
                    )
                )
                if not coverage_for_project is None:
                    coverage_for_project = int(coverage_for_project)

                b_coverage_default = True
                if default_project_software.is_snippy_single_parameter_default(
                    project_sample, DefaultParameters.SNIPPY_COVERAGE_NAME
                ):
                    coverage = get_coverage.get_coverage(
                        project_sample.get_file_output(
                            TypePath.MEDIA_ROOT,
                            FileType.FILE_DEPTH_GZ,
                            self.software_names.get_snippy_name(),
                        ),
                        project_sample.project.reference.get_reference_fasta(
                            TypePath.MEDIA_ROOT
                        ),
                        None,
                        coverage_for_project,
                    )
                else:
                    b_coverage_default = False
                    default_coverage_value = (
                        default_project_software.get_snippy_single_parameter(
                            project_sample, DefaultParameters.SNIPPY_COVERAGE_NAME
                        )
                    )
                    coverage = get_coverage.get_coverage(
                        project_sample.get_file_output(
                            TypePath.MEDIA_ROOT,
                            FileType.FILE_DEPTH_GZ,
                            self.software_names.get_snippy_name(),
                        ),
                        project_sample.project.reference.get_reference_fasta(
                            TypePath.MEDIA_ROOT
                        ),
                        int(default_coverage_value),
                        coverage_for_project,
                    )
                ################################
                ##################################
                ### set the alerts in the coverage
                ### remove possible previous alerts from others run
                for keys_to_remove in MetaKeyAndValue.VECT_TO_REMOVE_RUN_PROJECT_SAMPLE:
                    manageDatabase.remove_project_sample_start_metakey(
                        project_sample, keys_to_remove
                    )

                project_sample = ProjectSample.objects.get(pk=project_sample.id)
                project_sample.alert_second_level = 0
                project_sample.alert_first_level = 0
                for element in coverage.get_dict_data():
                    if not coverage.is_100_more_9(element) and b_coverage_default:
                        project_sample.alert_second_level += 1
                        meta_key = metaKeyAndValue.get_meta_key(
                            MetaKeyAndValue.META_KEY_ALERT_COVERAGE_9, element
                        )
                        manageDatabase.set_project_sample_metakey(
                            project_sample,
                            user,
                            meta_key,
                            MetaKeyAndValue.META_VALUE_Success,
                            coverage.get_fault_message_9(element),
                        )
                    elif (
                        not coverage.is_100_more_defined_by_user(element)
                        and not b_coverage_default
                    ):
                        project_sample.alert_second_level += 1
                        meta_key = metaKeyAndValue.get_meta_key(
                            MetaKeyAndValue.META_KEY_ALERT_COVERAGE_value_defined_by_user,
                            element,
                        )
                        manageDatabase.set_project_sample_metakey(
                            project_sample,
                            user,
                            meta_key,
                            MetaKeyAndValue.META_VALUE_Success,
                            coverage.get_fault_message_defined_by_user(
                                element, default_coverage_value
                            ),
                        )
                    elif not coverage.is_100_more_0(element):
                        project_sample.alert_first_level += 1
                        meta_key = metaKeyAndValue.get_meta_key(
                            MetaKeyAndValue.META_KEY_ALERT_COVERAGE_0, element
                        )
                        manageDatabase.set_project_sample_metakey(
                            project_sample,
                            user,
                            meta_key,
                            MetaKeyAndValue.META_VALUE_Success,
                            coverage.get_fault_message_0(element),
                        )
                project_sample.save()

                ## set the coverage in database
                meta_sample = manageDatabase.set_project_sample_metakey(
                    project_sample,
                    user,
                    MetaKeyAndValue.META_KEY_Coverage,
                    MetaKeyAndValue.META_VALUE_Success,
                    coverage.to_json(),
                )
            except Exception as e:
                result = Result()
                result.set_error("Fail to get coverage: " + e.args[0])
                result.add_software(
                    SoftwareDesc(
                        self.software_names.get_coverage_name(),
                        self.software_names.get_coverage_version(),
                        self.software_names.get_coverage_parameters(),
                    )
                )
                manageDatabase.set_project_sample_metakey(
                    project_sample,
                    user,
                    MetaKeyAndValue.META_KEY_Coverage,
                    MetaKeyAndValue.META_VALUE_Error,
                    result.to_json(),
                )

                ### get again and set error
                project_sample = ProjectSample.objects.get(pk=project_sample.id)
                project_sample.is_error = True
                project_sample.save()

                meta_sample = manageDatabase.get_project_sample_metakey_last(
                    project_sample,
                    meta_key_project_sample,
                    MetaKeyAndValue.META_VALUE_Queue,
                )
                if meta_sample != None:
                    manageDatabase.set_project_sample_metakey(
                        project_sample,
                        user,
                        meta_key_project_sample,
                        MetaKeyAndValue.META_VALUE_Error,
                        meta_sample.description,
                    )
                process_SGE.set_process_controler(
                    user,
                    process_controler.get_name_project_sample(project_sample),
                    ProcessControler.FLAG_ERROR,
                )
                return False

            #####################
            ###
            ### make mask the consensus SoftwareNames.SOFTWARE_MSA_MASKER
            limit_to_mask_consensus = int(
                default_project_software.get_mask_consensus_single_parameter(
                    project_sample,
                    DefaultParameters.MASK_CONSENSUS_threshold,
                    ConstantsSettings.TECHNOLOGY_illumina,
                )
            )
            msa_parameters = self.make_mask_consensus_by_deep(
                project_sample.get_file_output(
                    TypePath.MEDIA_ROOT,
                    FileType.FILE_CONSENSUS_FASTA,
                    self.software_names.get_snippy_name(),
                ),
                project_sample.project.reference.get_reference_fasta(
                    TypePath.MEDIA_ROOT
                ),
                project_sample.get_file_output(
                    TypePath.MEDIA_ROOT,
                    FileType.FILE_DEPTH_GZ,
                    self.software_names.get_snippy_name(),
                ),
                coverage,
                project_sample.sample.name,
                limit_to_mask_consensus,
            )
            ### add version of mask
            result_all.add_software(
                SoftwareDesc(
                    self.software_names.get_msa_masker_name(),
                    self.software_names.get_msa_masker_version(),
                    "{}; for coverages less than {} in {}% of the regions.".format(
                        msa_parameters,
                        default_project_software.get_snippy_single_parameter(
                            project_sample, DefaultParameters.SNIPPY_COVERAGE_NAME
                        ),
                        100 - limit_to_mask_consensus,
                    ),
                )
            )

            ## identify VARIANTS IN INCOMPLETE LOCUS in all locus, set yes in variants if are in areas with coverage problems
            ## transform 'synonymous_variant c.981A>G p.Glu327Glu' to ["synonymous_variant", "c.981A>G", "p.Glu327Glu"]
            parse_out_files = ParseOutFiles()
            parse_out_files.add_variants_in_incomplete_locus(
                project_sample.get_file_output(
                    TypePath.MEDIA_ROOT,
                    FileType.FILE_TAB,
                    self.software_names.get_snippy_name(),
                ),
                coverage,
            )

            ################
            ################
            ## run freebayes if at least one segment has some coverage
            ## test if it is necessary to run freebayes
            count_hits = CountHits()
            if default_project_software.is_to_run_freebayes(user, project_sample):
                try:
                    out_put_path = self.run_freebayes_parallel(
                        project_sample.get_file_output(
                            TypePath.MEDIA_ROOT,
                            FileType.FILE_BAM,
                            self.software_names.get_snippy_name(),
                        ),
                        project_sample.project.reference.get_reference_fasta(
                            TypePath.MEDIA_ROOT
                        ),
                        project_sample.project.reference.get_reference_gbk(
                            TypePath.MEDIA_ROOT
                        ),
                        project_sample.sample.name,
                    )
                    result_all.add_software(
                        SoftwareDesc(
                            self.software_names.get_freebayes_name(),
                            self.software_names.get_freebayes_version(),
                            self.software_names.get_freebayes_parameters(),
                        )
                    )
                except Exception as e:

                    ### can fail the freebayes parallel and try the regular one
                    try:
                        out_put_path = self.run_freebayes(
                            project_sample.get_file_output(
                                TypePath.MEDIA_ROOT,
                                FileType.FILE_BAM,
                                self.software_names.get_snippy_name(),
                            ),
                            project_sample.project.reference.get_reference_fasta(
                                TypePath.MEDIA_ROOT
                            ),
                            project_sample.project.reference.get_reference_gbk(
                                TypePath.MEDIA_ROOT
                            ),
                            project_sample.sample.name,
                        )
                        result_all.add_software(
                            SoftwareDesc(
                                self.software_names.get_freebayes_name(),
                                self.software_names.get_freebayes_version(),
                                self.software_names.get_freebayes_parameters(),
                            )
                        )
                    except Exception as e:
                        result = Result()
                        result.set_error(e.args[0])
                        result.add_software(
                            SoftwareDesc(
                                self.software_names.get_freebayes_name(),
                                self.software_names.get_freebayes_version(),
                                self.software_names.get_freebayes_parameters(),
                            )
                        )
                        manageDatabase.set_project_sample_metakey(
                            project_sample,
                            user,
                            MetaKeyAndValue.META_KEY_Freebayes,
                            MetaKeyAndValue.META_VALUE_Error,
                            result.to_json(),
                        )

                        ### get again and set error
                        project_sample = ProjectSample.objects.get(pk=project_sample.id)
                        project_sample.is_error = True
                        project_sample.save()

                        meta_sample = manageDatabase.get_project_sample_metakey_last(
                            project_sample,
                            meta_key_project_sample,
                            MetaKeyAndValue.META_VALUE_Queue,
                        )
                        if meta_sample != None:
                            manageDatabase.set_project_sample_metakey(
                                project_sample,
                                user,
                                meta_key_project_sample,
                                MetaKeyAndValue.META_VALUE_Error,
                                meta_sample.description,
                            )
                        process_SGE.set_process_controler(
                            user,
                            process_controler.get_name_project_sample(project_sample),
                            ProcessControler.FLAG_ERROR,
                        )
                        return False

                ## count hits from tab file
                count_hits = CountHits()
                if not out_put_path is None:
                    file_tab = os.path.join(
                        out_put_path, project_sample.sample.name + ".tab"
                    )
                    if os.path.exists(file_tab):
                        vect_count_type = ["snp"]  ## only detects snp
                        count_hits = self.utils.count_hits_from_tab(
                            file_tab, vect_count_type
                        )
                        ### set flag that is finished
                        manageDatabase.set_project_sample_metakey(
                            project_sample,
                            user,
                            MetaKeyAndValue.META_KEY_Count_Hits,
                            MetaKeyAndValue.META_VALUE_Success,
                            count_hits.to_json(),
                        )
                    else:
                        result = Result()
                        result.set_error("Fail to collect tab file from freebayes")
                        result.add_software(
                            SoftwareDesc(
                                self.software_names.get_freebayes_name(),
                                self.software_names.get_freebayes_version(),
                                self.software_names.get_freebayes_parameters(),
                            )
                        )
                        manageDatabase.set_project_sample_metakey(
                            project_sample,
                            user,
                            MetaKeyAndValue.META_KEY_Freebayes,
                            MetaKeyAndValue.META_VALUE_Error,
                            result.to_json(),
                        )

                        ### get again and set error
                        project_sample = ProjectSample.objects.get(pk=project_sample.id)
                        project_sample.is_error = True
                        project_sample.save()

                        meta_sample = manageDatabase.get_project_sample_metakey_last(
                            project_sample,
                            meta_key_project_sample,
                            MetaKeyAndValue.META_VALUE_Queue,
                        )
                        if meta_sample != None:
                            manageDatabase.set_project_sample_metakey(
                                project_sample,
                                user,
                                meta_key_project_sample,
                                MetaKeyAndValue.META_VALUE_Error,
                                meta_sample.description,
                            )
                        process_SGE.set_process_controler(
                            user,
                            process_controler.get_name_project_sample(project_sample),
                            ProcessControler.FLAG_ERROR,
                        )
                        return False

                    self.copy_files_to_project(
                        project_sample,
                        self.software_names.get_freebayes_name(),
                        out_put_path,
                    )
                    ## remove path dir if exist
                    self.utils.remove_dir(out_put_path)
                else:
                    ### set count hits to zero
                    manageDatabase.set_project_sample_metakey(
                        project_sample,
                        user,
                        MetaKeyAndValue.META_KEY_Count_Hits,
                        MetaKeyAndValue.META_VALUE_Success,
                        count_hits.to_json(),
                    )

                ### mixed infection
                try:
                    ## get instances
                    mixed_infections_management = MixedInfectionsManagement()

                    ## set the alert also
                    mixed_infection = mixed_infections_management.get_mixed_infections(
                        project_sample, user, count_hits
                    )
                except:
                    result = Result()
                    result.set_error("Fail to calculate mixed infextion")
                    result.add_software(SoftwareDesc("In house software", "1.0", ""))
                    manageDatabase.set_project_sample_metakey(
                        project_sample,
                        user,
                        MetaKeyAndValue.META_KEY_Mixed_Infection,
                        MetaKeyAndValue.META_VALUE_Error,
                        result.to_json(),
                    )

                    ### get again and set error
                    project_sample = ProjectSample.objects.get(pk=project_sample.id)
                    project_sample.is_error = True
                    project_sample.save()

                    meta_sample = manageDatabase.get_project_sample_metakey_last(
                        project_sample,
                        meta_key_project_sample,
                        MetaKeyAndValue.META_VALUE_Queue,
                    )
                    if meta_sample != None:
                        manageDatabase.set_project_sample_metakey(
                            project_sample,
                            user,
                            meta_key_project_sample,
                            MetaKeyAndValue.META_VALUE_Error,
                            meta_sample.description,
                        )
                    process_SGE.set_process_controler(
                        user,
                        process_controler.get_name_project_sample(project_sample),
                        ProcessControler.FLAG_ERROR,
                    )
                    return False
            else:
                ### set count hits to zero
                mixed_infections_management = MixedInfectionsManagement()
                mixed_infection = (
                    mixed_infections_management.get_mixed_infections_empty_value()
                )
                manageDatabase.set_project_sample_metakey(
                    project_sample,
                    user,
                    MetaKeyAndValue.META_KEY_Count_Hits,
                    MetaKeyAndValue.META_VALUE_Success,
                    count_hits.to_json(),
                )

                ### remove several files that can exist form previous interactions
                tab_freebayes_file = project_sample.get_file_output(
                    TypePath.MEDIA_ROOT,
                    FileType.FILE_TAB,
                    SoftwareNames.SOFTWARE_FREEBAYES_name,
                )
                self.utils.remove_file(tab_freebayes_file)
                file_out = project_sample.get_file_output_human(
                    TypePath.MEDIA_ROOT,
                    FileType.FILE_TAB,
                    SoftwareNames.SOFTWARE_FREEBAYES_name,
                )
                self.utils.remove_file(file_out)
                b_second_choice = True
                file_out = project_sample.get_file_output_human(
                    TypePath.MEDIA_ROOT,
                    FileType.FILE_TAB,
                    SoftwareNames.SOFTWARE_FREEBAYES_name,
                    b_second_choice,
                )
                self.utils.remove_file(file_out)

                file_out = project_sample.get_file_output(
                    TypePath.MEDIA_ROOT,
                    FileType.FILE_VCF,
                    SoftwareNames.SOFTWARE_FREEBAYES_name,
                )
                self.utils.remove_file(file_out)
                file_out = project_sample.get_file_output(
                    TypePath.MEDIA_ROOT,
                    FileType.FILE_VCF_GZ,
                    SoftwareNames.SOFTWARE_FREEBAYES_name,
                )
                self.utils.remove_file(file_out)

            ### draw coverage
            try:
                ### make the coverage images
                draw_all_coverage = DrawAllCoverage()
                draw_all_coverage.draw_all_coverages(project_sample)
            except:
                result = Result()
                result.set_error("Fail to draw coverage images")
                result.add_software(SoftwareDesc("In house software", "1.0", ""))
                manageDatabase.set_project_sample_metakey(
                    project_sample,
                    user,
                    MetaKeyAndValue.META_KEY_Coverage,
                    MetaKeyAndValue.META_VALUE_Error,
                    result.to_json(),
                )

                ### get again and set error
                project_sample = ProjectSample.objects.get(pk=project_sample.id)
                project_sample.is_error = True
                project_sample.save()

                meta_sample = manageDatabase.get_project_sample_metakey_last(
                    project_sample,
                    meta_key_project_sample,
                    MetaKeyAndValue.META_VALUE_Queue,
                )
                if meta_sample != None:
                    manageDatabase.set_project_sample_metakey(
                        project_sample,
                        user,
                        meta_key_project_sample,
                        MetaKeyAndValue.META_VALUE_Error,
                        meta_sample.description,
                    )
                process_SGE.set_process_controler(
                    user,
                    process_controler.get_name_project_sample(project_sample),
                    ProcessControler.FLAG_ERROR,
                )
                return False

            ### get again
            manage_database = ManageDatabase()
            project_sample = ProjectSample.objects.get(pk=project_sample.id)
            project_sample.is_finished = True
            project_sample.is_deleted = False
            project_sample.is_deleted_in_file_system = False
            project_sample.date_deleted = None
            project_sample.is_error = False
            project_sample.is_mask_consensus_sequences = True
            project_sample.count_variations = manage_database.get_variation_count(
                count_hits
            )
            project_sample.mixed_infections = mixed_infection
            project_sample.save()
            ### add today date, last change
            project = project_sample.project
            project.last_change_date = datetime.datetime.now()
            project.save()

            from utils.collect_extra_data import CollectExtraData

            collect_extra_data = CollectExtraData()

            ### get a clean freebayes file
            tab_freebayes_file = project_sample.get_file_output(
                TypePath.MEDIA_ROOT,
                FileType.FILE_TAB,
                SoftwareNames.SOFTWARE_FREEBAYES_name,
            )
            if os.path.exists(tab_freebayes_file):
                file_out = project_sample.get_file_output_human(
                    TypePath.MEDIA_ROOT,
                    FileType.FILE_TAB,
                    SoftwareNames.SOFTWARE_FREEBAYES_name,
                )
                vect_type_remove = ["ins", "del"]
                collect_extra_data.collect_variations_freebayes_only_one_file(
                    tab_freebayes_file, file_out, vect_type_remove
                )
                vect_type_remove = []
                b_second_choice = True
                file_out = project_sample.get_file_output_human(
                    TypePath.MEDIA_ROOT,
                    FileType.FILE_TAB,
                    SoftwareNames.SOFTWARE_FREEBAYES_name,
                    b_second_choice,
                )
                collect_extra_data.collect_variations_freebayes_only_one_file(
                    tab_freebayes_file, file_out, vect_type_remove
                )

            ### get clean consensus file
            consensus_fasta = project_sample.get_file_output(
                TypePath.MEDIA_ROOT,
                FileType.FILE_CONSENSUS_FASTA,
                SoftwareNames.SOFTWARE_SNIPPY_name,
            )
            if os.path.exists(consensus_fasta):
                file_out = project_sample.get_consensus_file(TypePath.MEDIA_ROOT)
                self.utils.filter_fasta_all_sequences_file(
                    consensus_fasta, coverage, file_out, limit_to_mask_consensus, False
                )

                ## make a backup of this file to use has a starter of second stage analysis
                self.utils.copy_file(
                    project_sample.get_consensus_file(TypePath.MEDIA_ROOT),
                    project_sample.get_backup_consensus_file(),
                )

            ### set the tag of result OK
            manageDatabase.set_project_sample_metakey(
                project_sample,
                user,
                MetaKeyAndValue.META_KEY_Snippy_Freebayes,
                MetaKeyAndValue.META_VALUE_Success,
                result_all.to_json(),
            )
            ### set the flag of the end of the task
            meta_sample = manageDatabase.get_project_sample_metakey_last(
                project_sample,
                meta_key_project_sample,
                MetaKeyAndValue.META_VALUE_Queue,
            )
            if not meta_sample is None:
                manageDatabase.set_project_sample_metakey(
                    project_sample,
                    user,
                    meta_key_project_sample,
                    MetaKeyAndValue.META_VALUE_Success,
                    meta_sample.description,
                )
        except:
            ## finished with error
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_project_sample(project_sample),
                ProcessControler.FLAG_ERROR,
            )
            return False

        ### finished
        process_SGE.set_process_controler(
            user,
            process_controler.get_name_project_sample(project_sample),
            ProcessControler.FLAG_FINISHED,
        )
        return True

    def create_index_files(self, file_name):
        """
        create index, need to be .bz
        """
        file_to_index = file_name
        if not file_to_index.endswith(FileExtensions.FILE_VCF_GZ):
            file_to_index += FileExtensions.FILE_VCF_GZ
        if not os.path.exists(file_to_index):
            self.logger_production.error("File doesn't exist: " + file_to_index)
            self.logger_debug.error("Fail doesn't exist: " + file_to_index)
            raise Exception("File doesn't exist")

        self.create_index_with_tabix(file_to_index)

    def create_index_with_tabix(self, file_name, index_type="vcf"):
        """
        :param index_type [gff, bed, sam, vcf]
        """
        ## test if tbi exists
        if os.path.exists(file_name + FileExtensions.FILE_TBI):
            return

        if index_type == SoftwareNames.SOFTWARE_DEPTH_SAMTOOLS_file_flag:
            ### index -s 1 -> chr name
            ### index -b 2 -> start position
            ### index -e 2 -> end position
            cmd = "{} -s 1 -b 2 -e 2 {}".format(
                self.software_names.get_tabix(), file_name
            )
        else:
            cmd = "{} -p {} {}".format(
                self.software_names.get_tabix(), index_type, file_name
            )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to create index")

    def create_index_files_from_igv_tools(self, file_name):
        """
        Create index from igvtools
        """
        cmd = "java -jar {} index {}".format(
            self.software_names.get_igvtools(), file_name
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to create index")

    def dos_2_unix(self, file_name):
        """
        convert dos 2 unix
        """
        if not os.path.exists(file_name):
            return
        cmd = "dos2unix {}".format(file_name)
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to create index")

    def zip_files_in_path(self, path_name):
        """
        zip all files inside a directory and return a file
        :out os.path.join(path_name, "zip_file_out.zip")
        """
        if not os.path.exists(path_name):
            return
        current_dir = os.getcwd()
        os.chdir(path_name)  ## change to path that will be compressed

        file_name_zip = "zip_file_out.zip"
        cmd = "zip -r {} *".format(file_name_zip)
        exist_status = os.system(cmd)
        if exist_status != 0:
            os.chdir(current_dir)
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to zip files in '{}'".format(path_name))
        os.chdir(current_dir)
        return os.path.join(path_name, file_name_zip)

    def fasta_2_upper(self, file_name):
        """
        covert fasta 2 upper
        """
        if not os.path.exists(file_name):
            return
        temp_file = self.utils.get_temp_file("fasta_2_upper", ".fasta")
        cmd = (
            "awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}'"
            + " {} > {}; mv {} {}".format(file_name, temp_file, temp_file, file_name)
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to create index")

    def set_first_sequence_fasta(self, file_name):
        """
        covert fasta 2 upper
        """
        vect_record = []
        with open(file_name) as handle_in:
            for record in SeqIO.parse(handle_in, "fasta"):
                vect_record.append(record)
                break

            if len(vect_record) > 0:
                with open(file_name, "w") as handle_fasta_out_align:
                    SeqIO.write(vect_record, handle_fasta_out_align, "fasta")

    def make_downsize(self, path_1, path_2, max_fastq_file):
        """
        Reduce file size from X to Constants.MAX_FASTQ_FILE
        return (TRUE|FALSE if was reduce or not, new_path_1_reduced, new_path_2_reduced)
        """
        file_size_max = 0
        if os.path.exists(path_1):
            file_size_max = os.path.getsize(path_1)
        if (
            not path_2 is None
            and os.path.exists(path_2)
            and file_size_max < os.path.getsize(path_2)
        ):
            file_size_max = os.path.getsize(path_2)

        ### need to make the down size
        if file_size_max > max_fastq_file:
            ratio = max_fastq_file / file_size_max

            ## small correction
            if ratio < 0.8:
                ratio += 0.1

            path_to_work = self.utils.get_temp_dir()
            path_1_temp = self.utils.get_temp_file_from_dir(
                path_to_work, "fastq_1", ".fastq"
            )
            path_2_temp = self.utils.get_temp_file_from_dir(
                path_to_work, "fastq_2", ".fastq"
            )

            self.utils.uncompress_files(
                self.software_names.get_gzip(), path_1, path_1_temp
            )
            file_names = path_1_temp
            if path_2 != None and os.path.exists(path_2):
                self.utils.uncompress_files(
                    self.software_names.get_gzip(), path_2, path_2_temp
                )
                file_names += " " + path_2_temp

            cmd = "{} -p {:.2f} -o {}/sample {}".format(
                self.software_names.get_fastqtools_sample(),
                ratio,
                path_to_work,
                file_names,
            )
            exist_status = os.system(cmd)
            if exist_status != 0:
                self.logger_production.error("Fail to run: " + cmd)
                self.logger_debug.error("Fail to run: " + cmd)
                raise Exception("Fail to downsize")

            if not path_2 is None and os.path.exists(path_2):
                self.utils.compress_files(
                    self.software_names.get_gzip(),
                    os.path.join(path_to_work, "sample.1.fastq"),
                )
                path_1_temp = os.path.join(
                    path_to_work, "sample.1.fastq" + FileExtensions.FILE_GZ
                )
                self.utils.compress_files(
                    self.software_names.get_gzip(),
                    os.path.join(path_to_work, "sample.2.fastq"),
                )
                path_2_temp = os.path.join(
                    path_to_work, "sample.2.fastq" + FileExtensions.FILE_GZ
                )
            else:
                self.utils.compress_files(
                    self.software_names.get_gzip(),
                    os.path.join(path_to_work, "sample.fastq"),
                )
                path_1_temp = os.path.join(
                    path_to_work, "sample.fastq" + FileExtensions.FILE_GZ
                )

            return (
                True,
                path_1_temp,
                path_2_temp if (path_2 != None and os.path.exists(path_2)) else None,
            )
        return (False, path_1, path_2)

    def make_mask_consensus_by_deep(
        self,
        consensus_file,
        reference_fasta,
        deep_file,
        coverage,
        sample_name,
        limit_make_mask,
    ):
        """
        :param limit_to_mask_consensus, default 70%
        /usr/local/software/insaflu/snippy/bin/msa_masker.py -i /tmp/insaFlu/insa_flu_path_86811930/run_snippy1_sdfs/run_snippy1_sdfs.consensus.fa
                -df /tmp/insaFlu/insa_flu_path_86811930/run_snippy1_sdfs/run_snippy1_sdfs.depth.gz
                -o /tmp/insaFlu/insa_flu_path_86811930/temp.fasta --c 200
        """
        ## run all elements in reference
        temp_masked = self.utils.get_temp_file("masked_file", ".fasta")
        temp_to_join = self.utils.get_temp_file("join_file", ".fasta")
        temp_mafft_align = self.utils.get_temp_file("mafft_to_align", ".fasta")
        temp_new_consensus = self.utils.get_temp_file("new_consensus", ".fasta")
        vect_out_fasta = []

        msa_parameters = ""
        with open(reference_fasta, "rU") as handle_fasta:
            dt_consensus = SeqIO.to_dict(SeqIO.parse(consensus_file, "fasta"))
            for record in SeqIO.parse(handle_fasta, "fasta"):
                if (
                    record.id in dt_consensus
                    and coverage.ratio_value_coverage_bigger_limit(
                        record.id, limit_make_mask
                    )
                ):  ### make mask

                    ### get sequences
                    vect_out_fasta_to_align = []
                    record_id = record.id
                    record.id = record.id + "_ref"
                    vect_out_fasta_to_align.append(record)
                    vect_out_fasta_to_align.append(dt_consensus[record_id])

                    with open(temp_to_join, "w") as handle_fasta_out_align:
                        SeqIO.write(
                            vect_out_fasta_to_align, handle_fasta_out_align, "fasta"
                        )

                    ### run maft
                    temp_mafft_align = self.run_mafft(
                        temp_to_join,
                        temp_mafft_align,
                        SoftwareNames.SOFTWARE_MAFFT_PARAMETERS,
                    )

                    ### run mask
                    msa_parameters = self.run_mask_app(
                        temp_mafft_align,
                        deep_file,
                        temp_masked,
                        coverage.get_middle_limit(),
                    )

                    ### read output file
                    dt_mask_consensus = SeqIO.to_dict(SeqIO.parse(temp_masked, "fasta"))
                    if record_id in dt_mask_consensus:
                        record_temp = dt_mask_consensus[record_id]
                    else:
                        record_temp = dt_consensus[record_id]
                    ## add sample name to the consensus sequences
                    record_temp.description = sample_name
                    vect_out_fasta.append(record_temp)
                else:  ## write as is
                    ## add sample name to the consensus sequences
                    record.description = sample_name
                    vect_out_fasta.append(record)

        ### write the output
        with open(temp_new_consensus, "w") as handle_fasta_out:
            if len(vect_out_fasta) > 0:
                SeqIO.write(vect_out_fasta, handle_fasta_out, "fasta")

        ### move temp consensus to original position, if has info
        if os.stat(temp_new_consensus).st_size > 0:
            self.utils.move_file(temp_new_consensus, consensus_file)

        self.utils.remove_file(temp_new_consensus)
        self.utils.remove_file(temp_to_join)
        self.utils.remove_file(temp_masked)
        self.utils.remove_file(temp_mafft_align)
        return msa_parameters

    def run_mask_app(self, input_fasta, deep_file, out_file, coverage_limit):
        """
        run msa_masker
        /usr/local/software/insaflu/snippy/bin/msa_masker.py -i /tmp/insaFlu/insa_flu_path_86811930/run_snippy1_sdfs/run_snippy1_sdfs.consensus.fa -df /tmp/insaFlu/insa_flu_path_86811930/run_snippy1_sdfs/run_snippy1_sdfs.depth.gz -o /tmp/insaFlu/insa_flu_path_86811930/temp.fasta --c 200
        """
        msa_parameters = "{} {}".format(
            self.software_names.get_msa_masker_parameters(), int(coverage_limit) - 1
        )
        cmd = "{} -i {} -df {} -o {} {}".format(
            self.software_names.get_msa_masker(),
            input_fasta,
            deep_file,
            out_file,
            msa_parameters,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run msa_masker")

        return msa_parameters

    def mask_sequence(
        self,
        sequence_ref,
        sequence_consensus,
        mask_sites,
        mask_from_beginning,
        mask_from_end,
        mask_range,
    ):
        """Mask characters at the given sites in a single sequence record, modifying the
        record in place.
        Parameters
        ----------
        sequence_ref : Bio.SeqIO.SeqRecord
            A sequence that act like a reference. All the positions are referent to this seq.
        sequence_consensus : Bio.SeqIO.SeqRecord
            A sequence to be masked
        mask_sites: list[int]
            A list of site indexes to exclude from the FASTA.
        mask_from_beginning: int
            Number of sites to mask from the beginning of each sequence (default 0)
        mask_from_end: int
            Number of sites to mask from the end of each sequence (default 0)
        mask_invalid: bool
            Mask invalid nucleotides (default False)
        Returns
        -------
        Bio.SeqIO.SeqRecord
            Masked sequence in its original record object
        """
        if len(sequence_ref.seq) == 0 or len(sequence_consensus.seq) == 0:
            return sequence_consensus
        # need to align
        sequence_length = len(sequence_ref.seq)
        seq_ref, seq_consensus = self.align_two_sequences(
            str(sequence_ref.seq), str(sequence_consensus.seq)
        )

        # Convert to a mutable sequence to enable masking with Ns.
        beginning = (
            int(mask_from_beginning)
            if not mask_from_beginning is None and len(mask_from_beginning) > 0
            else -1
        )
        end = (
            int(mask_from_end)
            if not mask_from_end is None and len(mask_from_end) > 0
            else -1
        )

        ## collecting all positions to maks
        dt_positions = {}
        if beginning != -1:
            for _ in range(0, beginning):
                dt_positions[_] = 1
        if end != -1:
            if (len(str(sequence_ref.seq)) - end) < 0:
                pos_from = 0
            else:
                pos_from = len(str(sequence_ref.seq)) - end
            for _ in range(pos_from, len(str(sequence_ref.seq))):
                dt_positions[_] = 1

        ## several sites
        if not mask_sites is None and len(mask_sites.split(",")[0]) > 0:
            for site in [int(_) - 1 for _ in mask_sites.split(",")]:
                if site < sequence_length:
                    dt_positions[site] = 1
        ## several ranges
        if not mask_range is None:
            for data_ in mask_range.split(","):
                if len(data_) > 0 and len(data_.split("-")) == 2:
                    for site in range(
                        int(data_.split("-")[0]) - 1, int(data_.split("-")[1])
                    ):
                        if site < sequence_length:
                            dt_positions[site] = 1

        ## mask positions
        masked_sequence = MutableSeq(seq_consensus)
        ref_insertions = 0
        ref_pos = 0
        for _ in range(len(seq_ref)):
            if seq_ref[_] == "-":
                ref_insertions += 1
                continue
            if ref_pos in dt_positions:
                masked_sequence[ref_pos + ref_insertions] = "N"
            ref_pos += 1
            if (ref_pos + ref_insertions) >= len(seq_consensus):
                break

        sequence_consensus.seq = Seq(str(masked_sequence).replace("-", ""))
        return sequence_consensus

    def mask_sequence_by_sites(
        self, reference_fasta_file, consensus_fasta_file, genetic_elemets
    ):
        """masking consensus file with positions related with reference elements"""
        vect_record_out = []
        ## always work with the backup
        with open(reference_fasta_file, "rU") as handle_ref, open(
            consensus_fasta_file, "rU"
        ) as handle_consensus:
            dict_record_ref = SeqIO.to_dict(SeqIO.parse(handle_ref, "fasta"))
            for record_consensus in SeqIO.parse(handle_consensus, "fasta"):
                masking_consensus = genetic_elemets.dt_elements_mask.get(
                    record_consensus.id, MaskingConsensus()
                )
                if (
                    masking_consensus.has_data()
                    and record_consensus.id in dict_record_ref
                ):
                    vect_record_out.append(
                        self.mask_sequence(
                            dict_record_ref[record_consensus.id],
                            record_consensus,
                            masking_consensus.mask_sites,
                            masking_consensus.mask_from_beginning,
                            masking_consensus.mask_from_ends,
                            masking_consensus.mask_regions,
                        )
                    )
                else:
                    vect_record_out.append(record_consensus)

        if len(vect_record_out) > 0:
            temp_file = self.utils.get_temp_file("masked_seq_", ".fasta")
            with open(temp_file, "w") as handle_fasta_out:
                SeqIO.write(vect_record_out, handle_fasta_out, "fasta")

            ### move temp consensus to original position, if has info
            if os.stat(temp_file).st_size > 0:
                self.utils.move_file(temp_file, consensus_fasta_file)
            else:
                os.unlink(temp_file)

    def align_two_sequences(self, ref_seq, consensus_seq, out_dir_temp=None):
        """
        :param  ref_seq: sequence with nucleotides
        :param  consensus_seq: sequence with nucleotides
        :out (ref_seq, consensus_seq)
        """

        ### set temp dir
        if out_dir_temp is None:
            out_dir = self.utils.get_temp_dir()
        else:
            out_dir = out_dir_temp

        temp_file_name = self.utils.get_temp_file_from_dir(
            out_dir, "to_align", ".fasta"
        )
        temp_file_name_out = self.utils.get_temp_file_from_dir(
            out_dir, "to_align_out", ".fasta"
        )
        records = []
        ref_seq_name = "ref"
        records.append(SeqRecord(Seq(ref_seq), id=ref_seq_name, description=""))
        records.append(SeqRecord(Seq(consensus_seq), id="consensus", description=""))

        ### save file
        with open(temp_file_name, "w") as handle_write:
            SeqIO.write(records, handle_write, "fasta")
        try:
            # self.software.run_mafft(temp_file_name, temp_file_name_out, SoftwareNames.SOFTWARE_MAFFT_PARAMETERS_TWO_SEQUENCES)
            self.run_mafft(
                temp_file_name,
                temp_file_name_out,
                SoftwareNames.SOFTWARE_MAFFT_PARAMETERS,
            )
        #             self.software.run_clustalo(temp_file_name, temp_file_name_out, SoftwareNames.SOFTWARE_CLUSTALO_PARAMETERS)
        except Exception as a:
            return "", ""

        ## test if the output is in fasta
        try:
            number_record = self.utils.is_fasta(temp_file_name_out)
            if number_record != 2:
                return "", ""
        except IOError as e:
            return "", ""

        ## read file
        with open(temp_file_name_out) as handle:
            ## get both sequences
            seq_ref = ""
            seq_other = ""
            for record_dict in SeqIO.parse(handle, "fasta"):
                if record_dict.id == ref_seq_name:  ## ref seq
                    seq_ref = str(record_dict.seq).upper()
                else:
                    seq_other = str(record_dict.seq).upper()

            if out_dir_temp is None:
                self.utils.remove_dir(out_dir)
        return seq_ref, seq_other

    # TODO remove after everything is settled with the specific builds...
    # Actually make it a generic function that calls specific subsections related to builds
    def run_nextstrain(
        self,
        alignments,
        metadata,
        build=SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_parameter,
        cores=1,
    ):
        """
        run nextstrain_ncov
        :param  alignments: sequence file with nucleotides
        :param  metadata: tabbed table file with properties
        :param  build: the specific nextstrain build to be used (defaults to a generic build)
        :param  cores: the number of cores to be used in nextstrain (defaults to 1)
        :out temp folder with all data (including results)
        """

        # Create a temp folder
        temp_dir = self.utils.get_temp_dir()

        # copy the base nexstrain folder to a temp folder
        # TODO Make a function copy_folder in utils
        cmd = (
            "cp -r "
            + SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_BASE
            + "/"
            + build
            + "/* "
            + temp_dir
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception(
                "Fail to copy nexstrain folder "
                + SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_BASE
                + "/"
                + build
                + "/* "
                + temp_dir
            )

        # Copy the setup folder and alignment and metadata files to the appropriate place in the temp folder
        self.utils.copy_file(
            alignments, os.path.join(temp_dir, "data", "sequences.fasta")
        )
        self.utils.copy_file(metadata, os.path.join(temp_dir, "data", "metadata.tsv"))

        # Copy the build-specific config file to the appropriate place in the temp folder
        config_file = os.path.join(
            getattr(settings, "STATIC_ROOT", None),
            Constants.DIR_NEXTSTRAIN_tables,
            build,
            "config.yaml",
        )
        self.utils.copy_file(
            config_file, os.path.join(temp_dir, "config", "config.yaml")
        )

        # to generate the include: may need to remove "" from the names...
        cmd = (
            "cat "
            + metadata
            + " | cut -f 1 | sed 's/\"//g' | tail -n +2 > "
            + os.path.join(temp_dir, "data", "include.txt")
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to generate include file in temp folder " + temp_dir)

        ### add reference sequence Fasta to global alignment
        reference_fasta = os.path.join(
            getattr(settings, "STATIC_ROOT", None),
            Constants.DIR_NEXTSTRAIN_tables,
            build,
            "references_sequences.fasta",
        )
        cmd = "cat {} >> {}".format(
            reference_fasta, os.path.join(temp_dir, "data", "sequences.fasta")
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.utils.remove_dir(temp_dir)
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run nextstrain. Command: " + cmd)

        # Run nextstrain
        # right now this is an exception... eventually make it generic, as the sofware swtup may be different according to the build?
        cmd = (
            SoftwareNames.SOFTWARE_NEXTSTRAIN
            + " build --native "
            + temp_dir
            + " --cores "
            + str(cores)
            + " --configfile "
            + temp_dir
            + "/config/config.yaml"
        )
        if build == "mpx":
            cmd = (
                SoftwareNames.SOFTWARE_NEXTSTRAIN_MPX
                + " build --native "
                + temp_dir
                + " --cores "
                + str(cores)
                + " --configfile "
                + temp_dir
                + "/config/config.yaml"
            )

        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise CmdException("Fail to run nextstrain.", cmd, temp_dir)

        # Collect results

        return temp_dir

    def run_nextstrain_ncov(self, alignments, metadata, cores=1):
        """
        run nextstrain_ncov
        :param  alignments: sequence file with nucleotides
        :param  metadata: tabbed table file with properties
        :param  cores: the number of cores to be used in nextstrain (defaults to 1)
        :out temp folder with all data (including results)
        """

        # Create a temp folder
        temp_dir = self.utils.get_temp_dir()

        build = "ncov"
        # copy the base nexstrain folder to a temp folder
        # TODO Make a function copy_folder in utils
        cmd = (
            "cp -r "
            + SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_BASE
            + "/"
            + build
            + "/* "
            + temp_dir
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception(
                "Fail to copy nexstrain folder "
                + SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_BASE
                + "/"
                + build
                + "/* "
                + temp_dir
            )

        # TODO?: Add this to the collect_consensus section...?
        # Add references from the build to the consensus sequences (TODO merge fasta with function instead with a cat...)
        cmd = "cat {} {} > {}".format(
            alignments,
            os.path.join(temp_dir, "data", "references_sequences.fasta"),
            os.path.join(temp_dir, "data", "sequences.fasta"),
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception(
                "Fail to concatenate ncov references with alignments in temp folder "
                + temp_dir
            )

        # Note this only works for now, IF the metadata entry is the last row (because the columns won't match...)
        self.utils.copy_file(metadata, os.path.join(temp_dir, "data", "metadata.tsv"))
        # cmd = "cat {} {} > {}".format(metadata, os.path.join(temp_dir, 'data', 'references_metadata.tsv'), os.path.join(temp_dir, 'data', 'metadata.tsv'))
        # exit_status = os.system(cmd)
        # if (exit_status != 0):
        #    self.logger_production.error('Fail to run: ' + cmd)
        #    self.logger_debug.error('Fail to run: ' + cmd)
        #    raise Exception("Fail to concatenate ncov reference metadata with metadata in temp folder " + temp_dir)

        cmd = (
            "cat "
            + metadata
            + " | cut -f 1 | sed 's/\"//g' | tail -n +2 > "
            + os.path.join(temp_dir, "data", "include.txt")
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to generate include file in temp folder " + temp_dir)

        cmd = (
            SoftwareNames.SOFTWARE_NEXTSTRAIN
            + " build --native "
            + temp_dir
            + " --cores "
            + str(cores)
            + " --configfile "
            + temp_dir
            + "/config/config.yaml"
        )

        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise CmdException(
                message="Fail to run nextstrain.", cmd=cmd, output_path=temp_dir
            )

        tree_file = self.utils.get_temp_file("treefile.nwk", sz_type="nwk")
        # Convert json to tree
        cmd = "{} --tree {} --output-tree {}".format(
            os.path.join(settings.DIR_SOFTWARE, "nextstrain/auspice_tree_to_table.sh"),
            os.path.join(temp_dir, "auspice", "ncov_current.json"),
            tree_file,
        )

        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise CmdException(
                message="Fail to run conversion of json to tree.",
                cmd=cmd,
                output_path=temp_dir,
            )

        # Copy log folder to auspice to be included in the zip
        cmd = "cp -r {} {}".format(
            os.path.join(temp_dir, ".snakemake", "log"),
            os.path.join(temp_dir, "auspice"),
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise CmdException(
                message="Fail to copy log to output folder.",
                cmd=cmd,
                output_path=temp_dir,
            )

        # Collect results
        zip_out = self.zip_files_in_path(os.path.join(temp_dir, "auspice"))
        auspice_zip = self.utils.get_temp_file("tempfile.zip", sz_type="zip")
        self.utils.move_file(zip_out, auspice_zip)

        # Put out the alignments too...
        # results/aligned_current.fasta.xz
        exit_status = os.system(
            "xz -d {}".format(
                os.path.join(temp_dir, "results", "aligned_current.fasta.xz")
            )
        )
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise CmdException(
                message="Fail to run unzip alignment file.",
                cmd=cmd,
                output_path=temp_dir,
            )

        alignment_file = self.utils.get_temp_file("aligned.fasta", sz_type="fasta")
        self.utils.move_file(
            os.path.join(temp_dir, "results", "aligned_current.fasta"), alignment_file
        )

        self.utils.remove_dir(temp_dir)

        return [tree_file, alignment_file, auspice_zip]

    def run_nextstrain_generic(
        self, alignments, metadata, ref_fasta, ref_genbank, cores=1
    ):
        """
        run nextstrain
        :param  alignments: sequence file with nucleotides
        :param  metadata: tabbed table file with properties
        :param  ref_fasta: the reference fasta file to be used
        :param  ref_genbank: the reference genbank file to be used (MUST correspond to fasta)
        :param  cores: the number of cores to be used in nextstrain (defaults to 1)
        :out temp folder with all data (including results)
        """

        # make sure fasta and genbank correspond...
        self.utils.compare_locus_fasta_gb(ref_fasta, ref_genbank)

        # Assume the reference as one of the segments: this is unlikely to work if there is more than one segment
        reference = self.utils.get_elements_and_genes(
            ref_genbank
        ).get_sorted_elements()[0]

        # Create a temp folder
        temp_dir = self.utils.get_temp_dir()

        # copy the base nexstrain folder to a temp folder
        # TODO Make a function copy_folder in utils
        build = "generic"
        cmd = (
            "cp -r "
            + SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_BASE
            + "/"
            + build
            + "/* "
            + temp_dir
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception(
                "Fail to copy nexstrain folder "
                + SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_BASE
                + "/"
                + build
                + "/* "
                + temp_dir
            )

        # add reference.gb and reference.fasta to data folder
        self.utils.copy_file(
            ref_fasta, os.path.join(temp_dir, "data", "reference.fasta")
        )
        self.utils.copy_file(
            ref_genbank, os.path.join(temp_dir, "data", "reference.gb")
        )

        # add sequences.fasta and metadata.tsv to data folder
        self.utils.copy_file(
            alignments, os.path.join(temp_dir, "data", "sequences.fasta")
        )
        self.utils.copy_file(metadata, os.path.join(temp_dir, "data", "metadata.tsv"))

        # add 'root = "reference_id"' at the top of Snakefile_base, creating the final Snakefile
        cmd = "cat {} | sed 's/REFID/{}/' > {}".format(
            os.path.join(temp_dir, "Snakefile_base"),
            reference,
            os.path.join(temp_dir, "Snakefile"),
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception(
                "Fail to change root in config file in temp folder " + temp_dir
            )

        # Now run Nextstrain
        cmd = (
            SoftwareNames.SOFTWARE_NEXTSTRAIN
            + " build --native "
            + temp_dir
            + " --cores "
            + str(cores)
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise CmdException(
                message="Fail to run nextstrain.", cmd=cmd, output_path=temp_dir
            )

        tree_file = self.utils.get_temp_file("treefile.nwk", sz_type="nwk")
        # Convert json to tree
        cmd = "{} --tree {} --output-tree {}".format(
            os.path.join(settings.DIR_SOFTWARE, "nextstrain/auspice_tree_to_table.sh"),
            os.path.join(temp_dir, "auspice", "generic.json"),
            tree_file,
        )

        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise CmdException(
                message="Fail to run conversion of json to tree.",
                cmd=cmd,
                output_path=temp_dir,
            )

        # Copy log folder to auspice to be included in the zip
        cmd = "cp -r {} {}".format(
            os.path.join(temp_dir, ".snakemake", "log"),
            os.path.join(temp_dir, "auspice"),
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise CmdException(
                message="Fail to copy log to output folder.",
                cmd=cmd,
                output_path=temp_dir,
            )

        # Collect results
        zip_out = self.zip_files_in_path(os.path.join(temp_dir, "auspice"))
        auspice_zip = self.utils.get_temp_file("tempfile.zip", sz_type="zip")
        self.utils.move_file(zip_out, auspice_zip)

        alignment_file = self.utils.get_temp_file("aligned.fasta", sz_type="fasta")
        self.utils.move_file(
            os.path.join(temp_dir, "results", "aligned.fasta"), alignment_file
        )

        self.utils.remove_dir(temp_dir)

        # results/aligned.fasta

        return [tree_file, alignment_file, auspice_zip]

    def run_nextstrain_flu(
        self, alignments, metadata, strain="h3n2", period="12y", cores=1
    ):
        """
        run nextstrain
        :param  alignments: sequence file with nucleotides
        :param  metadata: tabbed table file with properties
        :param  strain: flu strain (one of: h3n2, h1n1pdm, vic, yam)
        :param  period: the reference time period (one of: 6m, 2y, 3y, 6y, 12y, 60y)
        :param  cores: the number of cores to be used in nextstrain (defaults to 1)
        :out temp folder with all data (including results)
        """

        # Create a temp folder
        temp_dir = self.utils.get_temp_dir()

        # copy the base nexstrain folder to a temp folder
        # TODO Make a function copy_folder in utils
        build = "flu"
        cmd = (
            "cp -r "
            + SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_BASE
            + "/"
            + build
            + "/* "
            + temp_dir
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception(
                "Fail to copy nexstrain folder "
                + SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_BASE
                + "/"
                + build
                + "/* "
                + temp_dir
            )

        # add sequences.fasta and metadata.tsv to data folder
        self.utils.copy_file(
            alignments,
            os.path.join(temp_dir, "data", "sequences_" + strain + "_ha.fasta"),
        )
        self.utils.copy_file(
            metadata, os.path.join(temp_dir, "data", "metadata_" + strain + "_ha.tsv")
        )

        # Now run Nextstrain
        cmd = "{} build --native {} targets/flu_{}_ha_{} --cores {}".format(
            SoftwareNames.SOFTWARE_NEXTSTRAIN, temp_dir, strain, period, str(cores)
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise CmdException(
                message="Fail to run nextstrain.", cmd=cmd, output_path=temp_dir
            )

        tree_file = self.utils.get_temp_file("treefile.nwk", sz_type="nwk")
        # Convert json to tree
        cmd = "{} --tree {} --output-tree {}".format(
            os.path.join(settings.DIR_SOFTWARE, "nextstrain/auspice_tree_to_table.sh"),
            os.path.join(
                temp_dir, "auspice", "flu_" + strain + "_ha_" + period + ".json"
            ),
            tree_file,
        )

        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise CmdException(
                message="Fail to run conversion of json to tree.",
                cmd=cmd,
                output_path=temp_dir,
            )

        # Copy log folder to auspice to be included in the zip
        cmd = "cp -r {} {}".format(
            os.path.join(temp_dir, ".snakemake", "log"),
            os.path.join(temp_dir, "auspice"),
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise CmdException(
                message="Fail to copy log to output folder.",
                cmd=cmd,
                output_path=temp_dir,
            )

        # Collect results
        zip_out = self.zip_files_in_path(os.path.join(temp_dir, "auspice"))
        auspice_zip = self.utils.get_temp_file("tempfile.zip", sz_type="zip")
        self.utils.move_file(zip_out, auspice_zip)

        # results/aligned_h3n2_ha_12y.fasta
        alignment_file = self.utils.get_temp_file("aligned.fasta", sz_type="fasta")
        self.utils.move_file(
            os.path.join(
                temp_dir, "results", "aligned_{}_ha_{}.fasta".format(strain, period)
            ),
            alignment_file,
        )

        self.utils.remove_dir(temp_dir)

        return [tree_file, alignment_file, auspice_zip]

    def run_nextstrain_mpx(self, alignments, metadata, cores=1):
        """
        run nextstrain
        :param  alignments: sequence file with nucleotides
        :param  metadata: tabbed table file with properties
        :param  cores: the number of cores to be used in nextstrain (defaults to 1)
        :out temp folder with all data (including results)
        """

        # Create a temp folder
        temp_dir = self.utils.get_temp_dir()

        # copy the base nexstrain folder to a temp folder
        # TODO Make a function copy_folder in utils
        build = "mpx"
        cmd = (
            "cp -r "
            + SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_BASE
            + "/"
            + build
            + "/* "
            + temp_dir
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception(
                "Fail to copy nexstrain folder "
                + SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_BASE
                + "/"
                + build
                + "/* "
                + temp_dir
            )

        # add sequences.fasta and metadata.tsv to data folder
        # self.utils.copy_file(alignments,os.path.join(temp_dir, 'data', "sequences.fasta"))
        self.utils.copy_file(metadata, os.path.join(temp_dir, "data", "metadata.tsv"))

        cmd = "cat {} {} > {}".format(
            os.path.join(temp_dir, "data", "references_sequences.fasta"),
            alignments,
            os.path.join(temp_dir, "data", "sequences.fasta"),
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception(
                "Fail to concatenate references with consensus in temp folder "
            )

        # Now run Nextstrain
        cmd = "{} build --native {} --cores {}  --configfile config/config_hmpxv1_big.yaml".format(
            SoftwareNames.SOFTWARE_NEXTSTRAIN_MPX, temp_dir, str(cores)
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise CmdException(
                message="Fail to run nextstrain.", cmd=cmd, output_path=temp_dir
            )

        tree_file = self.utils.get_temp_file("treefile.nwk", sz_type="nwk")
        # Convert json to tree
        cmd = "{} --tree {} --output-tree {}".format(
            os.path.join(settings.DIR_SOFTWARE, "nextstrain/auspice_tree_to_table.sh"),
            os.path.join(temp_dir, "auspice", "monkeypox_hmpxv1_big.json"),
            tree_file,
        )

        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise CmdException(
                message="Fail to run conversion of json to tree.",
                cmd=cmd,
                output_path=temp_dir,
            )

        # Copy log folder to auspice to be included in the zip
        cmd = "cp -r {} {}".format(
            os.path.join(temp_dir, ".snakemake", "log"),
            os.path.join(temp_dir, "auspice"),
        )
        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise CmdException(
                message="Fail to copy log to output folder.",
                cmd=cmd,
                output_path=temp_dir,
            )

        # Collect results
        zip_out = self.zip_files_in_path(os.path.join(temp_dir, "auspice"))
        auspice_zip = self.utils.get_temp_file("tempfile.zip", sz_type="zip")
        self.utils.move_file(zip_out, auspice_zip)

        # results/hmpxv1_big/aligned.fasta
        alignment_file = self.utils.get_temp_file("aligned.fasta", sz_type="fasta")
        self.utils.move_file(
            os.path.join(temp_dir, "results", "hmpxv1_big", "aligned.fasta"),
            alignment_file,
        )

        self.utils.remove_dir(temp_dir)

        return [tree_file, alignment_file, auspice_zip]

    def run_aln2pheno(
        self,
        sequences,
        reference,
        gene,
        report,
        flagged,
        db="DB_COG_UK_antigenic_mutations_2022-05-30.tsv",
    ):
        """
        run aln2pheno
        :param sequences: sequence file with aminoacids from the SARS-CoV-2 S protein
        :param reference: name of the reference (must be one of the sequences)
        :param reference: name of the gene in use
        :param report: output file with final report
        :param flagged: output file with flagged mutations
        :param db: database for aln2pheno
        :out exit status
        """

        # Create a temp folder
        temp_dir = self.utils.get_temp_dir()

        # Add as parameter...
        db_file = os.path.join(settings.STATIC_ROOT, Constants.DIR_TYPE_ALN2PHENO, db)

        # Run aln2pheno
        cmd = "{} --db {} -g {} --algn {} -r {} --odir {} --output prefix -f".format(
            SoftwareNames.SOFTWARE_ALN2PHENO,
            db_file,
            gene,
            sequences,
            reference,
            temp_dir,
        )

        exit_status = os.system(cmd)
        if exit_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run aln2pheno in temp folder " + temp_dir)

        # copy results to output
        self.utils.copy_file(temp_dir + "/prefix_final_report.tsv", report)
        self.utils.copy_file(temp_dir + "/prefix_flagged_mutation_report.tsv", flagged)

        ###
        self.utils.remove_dir(temp_dir)
        return exit_status


class Contigs2Sequences(object):
    """
    classdocs
    """

    DATABASE_NAME = "influenza_assign_segments2contigs"

    utils = Utils()

    def __init__(self, b_testing):
        """
        Constructor
        """
        self.b_testing = b_testing
        self.root_path = os.path.join(
            getattr(settings, "STATIC_ROOT", None), "tests" if b_testing else ""
        )

    def get_most_recent_database(self):
        """ """
        version = 0
        path_to_return = None
        path_to_find = os.path.join(
            self.root_path, Constants.DIR_TYPE_CONTIGS_2_SEQUENCES
        )
        for file in self.utils.get_all_files(path_to_find):
            base_name = os.path.basename(file)
            match = re.search("\w+(_[v|V]\d+)\.\w+", base_name)
            if match == None:
                continue
            if len(match.regs) == 2:
                try:
                    self.utils.is_fasta(os.path.join(path_to_find, file))
                except IOError as e:
                    continue

                (temp, path) = (
                    int(match.group(1).lower().replace("_v", "")),
                    os.path.join(path_to_find, file),
                )
                if temp > version:
                    version = temp
                    path_to_return = path
        return (str(version), path_to_return)

    def get_database_name(self):
        ### get database file name
        (version, database_file_name) = self.get_most_recent_database()
        return self.utils.clean_extension(os.path.basename(database_file_name))

    def identify_contigs(self, file_name, file_name_out, b_create_fasta=True):
        """
        params database_name: if not a database_name is going to test the last one
        identify contigs
        params in: fasta file from spades
        out: fasta file with low coverage removed and Elements ID in description

        Change line 110 in abricate from (because of ONT approach, identify in reads instead of contigs):
        . " blastn -db \Q$db_path\E -outfmt '$format'"
        to:
        . " blastn -db \Q$db_path\E -outfmt '$format' -num_threads 3"
        """
        software = Software()

        ### get database file name, if it is not passed
        (version, database_file_name) = self.get_most_recent_database()
        database_name = self.get_database_name()

        ### first create database
        if not software.is_exist_database_abricate(database_name):
            software.create_database_abricate(database_name, database_file_name)

        out_file = self.utils.get_temp_file(
            "abricate_contig2seq", FileExtensions.FILE_TXT
        )
        ### run abricate
        software.run_abricate(
            database_name,
            file_name,
            SoftwareNames.SOFTWARE_ABRICATE_PARAMETERS_mincov_30,
            out_file,
        )

        parseOutFiles = ParseOutFiles()
        (dict_data_out, clean_abricate_file) = parseOutFiles.parse_abricate_file(
            out_file,
            file_name_out,
            SoftwareNames.SOFTWARE_SPAdes_CLEAN_HITS_BELLOW_VALUE,
        )

        out_file_fasta = None
        if b_create_fasta:
            vect_out_fasta = []
            vect_out_fasta_without_id = []
            out_file_fasta = self.utils.get_temp_file(
                "abricate_out_identified", FileExtensions.FILE_FASTA
            )
            with open(file_name) as handle_in, open(out_file_fasta, "w") as handle_out:
                for record in SeqIO.parse(handle_in, Constants.FORMAT_FASTA):
                    vect_possible_id = []
                    #                 for dict_data in vect_data:
                    #                     if (dict_data['Seq_Name'] == record.name): vect_possible_id.append(dict_data['Gene'])
                    for dict_data in dict_data_out.get(record.name, []):
                        vect_possible_id.append(dict_data["Gene"])
                    if len(vect_possible_id) > 0:
                        vect_out_fasta.append(
                            SeqRecord(
                                Seq(str(record.seq)),
                                id="_".join(record.id.split(".")[0].split("_")[:4]),
                                description=";".join(vect_possible_id),
                            )
                        )
                    ## NEED to check coverage for CANU
                    elif (
                        record.id.find("_") != -1
                        and float(record.id.split("_")[-1])
                        > SoftwareNames.SOFTWARE_SPAdes_CLEAN_HITS_BELLOW_VALUE
                    ):
                        vect_out_fasta_without_id.append(
                            SeqRecord(
                                Seq(str(record.seq)), id=record.id, description=""
                            )
                        )

                if len(vect_out_fasta) > 0 or len(vect_out_fasta_without_id) > 0:
                    vect_out_fasta.extend(vect_out_fasta_without_id)
                    SeqIO.write(vect_out_fasta, handle_out, "fasta")

        if os.path.exists(out_file):
            os.unlink(out_file)
        return (out_file_fasta, clean_abricate_file)
