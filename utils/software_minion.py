"""
Created on 01/01/2021

@author: mmp
"""
import datetime
import logging
import os

from constants.constants import Constants, FileType, TypePath
from constants.meta_key_and_values import MetaKeyAndValue
from constants.software_names import SoftwareNames
from django.conf import settings
from django.template.defaultfilters import filesizeformat
from managing_files.manage_database import ManageDatabase
from managing_files.models import (
    MixedInfectionsTag,
    ProcessControler,
    ProjectSample,
    Sample,
)
from pysam import pysam
from settings.constants_settings import ConstantsSettings
from settings.default_parameters import DefaultParameters
from settings.default_software_project_sample import DefaultProjectSoftware
from settings.models import Software as SoftwareSettings

from utils.coverage import DrawAllCoverage
from utils.mixed_infections_management import MixedInfectionsManagement
from utils.parse_coverage_file import GetCoverage
from utils.parse_out_files import ParseOutFiles
from utils.process_SGE import ProcessSGE
from utils.result import (
    CountHits,
    DecodeObjects,
    KeyValue,
    MaskingConsensus,
    Result,
    ResultAverageAndNumberReads,
    SoftwareDesc,
)
from utils.software import Software
from utils.utils import Utils


######################################
####   Minion methods
class SoftwareMinion(object):
    """
    classdocs
    """

    utils = Utils()
    software_names = SoftwareNames()
    software = Software()

    ## logging
    logger_debug = logging.getLogger("fluWebVirus.debug")
    logger_production = logging.getLogger("fluWebVirus.production")

    def __init__(self):
        """
        Constructor
        """
        pass

    def run_clean_minion(self, sample, user, b_make_identify_species=False):
        """
        Global processing, RabbitQC, NanoStat, NanoFilt and GetSpecies
        """
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
        if not sample.mixed_infections_tag is None:
            sample.mixed_infections_tag = None
        sample.number_alerts = 0
        sample.save()

        try:

            ### run stat and rabbit for Images
            b_has_data, b_it_ran = self.run_nanofilt_and_stat(sample, user)

            ### test Abricate ON/OFF
            default_software_project = DefaultProjectSoftware()
            b_make_identify_species = default_software_project.is_to_run_abricate(
                sample.owner, sample, ConstantsSettings.TECHNOLOGY_minion
            )

            ### queue the quality check and
            if (
                b_has_data and b_make_identify_species
            ):  ## don't run for single file because spades doesn't work for one single file
                self.software.identify_type_and_sub_type(
                    sample, sample.get_fastq_available(TypePath.MEDIA_ROOT), None, user
                )

            ## set the flag that is ready for process
            sample_to_update = Sample.objects.get(pk=sample.id)
            sample_to_update.is_sample_in_the_queue = False
            if b_has_data:
                sample_to_update.is_ready_for_projects = True

                ### make identify species
                if b_make_identify_species:
                    sample_to_update.type_subtype = (
                        sample_to_update.get_type_sub_type()[
                            : Sample.TYPE_SUBTYPE_LENGTH - 1
                        ]
                    )
                    (
                        tag_mixed_infection,
                        alert,
                        message,
                    ) = sample_to_update.get_mixed_infection()
                    if sample_to_update.number_alerts == None:
                        sample_to_update.number_alerts = alert
                    else:
                        sample_to_update.number_alerts += alert

                    manage_database = ManageDatabase()
                    if message != None and len(message) > 0:
                        manage_database.set_sample_metakey(
                            sample,
                            user,
                            MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE,
                            MetaKeyAndValue.META_VALUE_Success,
                            message,
                        )

                    ### save tag mixed_infecion
                    manage_database.set_sample_metakey(
                        sample,
                        user,
                        MetaKeyAndValue.META_KEY_TAG_MIXED_INFECTION_TYPE_SUBTYPE,
                        MetaKeyAndValue.META_VALUE_Success,
                        tag_mixed_infection,
                    )

                    try:
                        mixed_infections_tag = MixedInfectionsTag.objects.get(
                            name=tag_mixed_infection
                        )
                    except MixedInfectionsTag.DoesNotExist as e:
                        mixed_infections_tag = MixedInfectionsTag()
                        mixed_infections_tag.name = tag_mixed_infection
                        mixed_infections_tag.save()

                    sample_to_update.mixed_infections_tag = mixed_infections_tag
                else:
                    sample_to_update.type_subtype = Constants.EMPTY_VALUE_NA
                    tag_mixed_infection = Constants.EMPTY_VALUE_NA
                    try:
                        mixed_infections_tag = MixedInfectionsTag.objects.get(
                            name=tag_mixed_infection
                        )
                    except MixedInfectionsTag.DoesNotExist as e:
                        mixed_infections_tag = MixedInfectionsTag()
                        mixed_infections_tag.name = tag_mixed_infection
                        mixed_infections_tag.save()

                    sample_to_update.mixed_infections_tag = mixed_infections_tag

                    manage_database = ManageDatabase()
                    message = "Info: Abricate turned OFF by the user."
                    manage_database.set_sample_metakey(
                        sample,
                        user,
                        MetaKeyAndValue.META_KEY_ALERT_MIXED_INFECTION_TYPE_SUBTYPE,
                        MetaKeyAndValue.META_VALUE_Success,
                        message,
                    )
            else:
                manage_database = ManageDatabase()
                manage_database.set_sample_metakey(
                    sample_to_update,
                    user,
                    MetaKeyAndValue.META_KEY_ALERT_NO_READS_AFTER_FILTERING,
                    MetaKeyAndValue.META_VALUE_Success,
                    "Warning: no reads left after filtering.",
                )

                if sample_to_update.number_alerts == None:
                    sample_to_update.number_alerts = 1
                else:
                    sample_to_update.number_alerts += 1
                sample_to_update.is_ready_for_projects = False
                sample_to_update.type_subtype = Constants.EMPTY_VALUE_TYPE_SUBTYPE
            sample_to_update.save()

            ### set the flag of the end of the task
            meta_sample = manage_database.get_sample_metakey_last(
                sample,
                MetaKeyAndValue.META_KEY_Queue_TaskID,
                MetaKeyAndValue.META_VALUE_Queue,
            )
            if meta_sample != None:
                manage_database.set_sample_metakey(
                    sample,
                    sample.owner,
                    MetaKeyAndValue.META_KEY_Queue_TaskID,
                    MetaKeyAndValue.META_VALUE_Success,
                    meta_sample.description,
                )

        except Exception as e:
            print("Exception in run_clean_minion")
            print(e)
            import traceback

            traceback.print_exc()
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

    def run_nanofilt_and_stat(self, sample, owner):
        """
        run clean and stat before and after
        :output (Has data?, Is NanoFilt run?)
        """
        manage_database = ManageDatabase()
        result_all = Result()

        ### first try run down size if necessary
        if settings.DOWN_SIZE_FASTQ_FILES:
            (is_downsized, file_name_1, file_name_2) = self.software.make_downsize(
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

        ### get number of sequences in fastq
        (number_of_sequences, average_1, std1) = self.utils.get_number_sequences_fastq(
            sample.get_fastq(TypePath.MEDIA_ROOT, True)
        )

        ### first run stat
        try:
            result_nano_stat = self.run_nanostat(
                sample.get_fastq(TypePath.MEDIA_ROOT, True), number_of_sequences
            )
            result_all.add_software(
                SoftwareDesc(
                    self.software_names.get_NanoStat_name(),
                    self.software_names.get_NanoStat_version(),
                    self.software_names.get_NanoStat_parameters(),
                    result_nano_stat.key_values,
                )
            )

        except Exception as e:
            result = Result()
            result.set_error("Fail to run nanostat software: " + e.args[0])
            result.add_software(
                SoftwareDesc(
                    self.software_names.get_NanoStat_name(),
                    self.software_names.get_NanoStat_version(),
                    self.software_names.get_NanoStat_parameters(),
                )
            )
            manage_database.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_NanoStat_NanoFilt,
                MetaKeyAndValue.META_VALUE_Error,
                result.to_json(),
            )
            return False, False

        ## if they are sequences
        if number_of_sequences > 0:
            ### run rabbitQC, only to have a image of the data
            try:
                out_file_html = self.run_rabbitQC(
                    sample.get_fastq(TypePath.MEDIA_ROOT, True)
                )
                result_all.add_software(
                    SoftwareDesc(
                        self.software_names.get_rabbitQC_name(),
                        self.software_names.get_rabbitQC_version(),
                        self.software_names.get_rabbitQC_parameters(),
                    )
                )
                self.utils.copy_file(
                    out_file_html, sample.get_rabbitQC_output(TypePath.MEDIA_ROOT)
                )
                self.utils.remove_file(out_file_html)
            except Exception as e:
                result = Result()
                result.set_error("Fail to run rabbitQC software: " + e.args[0])
                result.add_software(
                    SoftwareDesc(
                        self.software_names.get_rabbitQC_name(),
                        self.software_names.get_rabbitQC_version(),
                        self.software_names.get_rabbitQC_parameters(),
                    )
                )
                manage_database.set_sample_metakey(
                    sample,
                    owner,
                    MetaKeyAndValue.META_KEY_NanoStat_NanoFilt,
                    MetaKeyAndValue.META_VALUE_Error,
                    result.to_json(),
                )
                return False, False

        default_software_project = DefaultProjectSoftware()
        default_software_project.test_all_defaults(sample.owner, None, None, sample)

        ### test if the software ran
        if number_of_sequences == 0 or not default_software_project.is_to_run_nanofilt(
            sample.owner, sample
        ):
            ### collect numbers to show on the sample table
            result_average = ResultAverageAndNumberReads(
                number_of_sequences, average_1, None, None
            )
            manage_database.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_Number_And_Average_Reads,
                MetaKeyAndValue.META_VALUE_Success,
                result_average.to_json(),
            )

            ## save everything OK
            manage_database.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_NanoStat_NanoFilt,
                MetaKeyAndValue.META_VALUE_Success,
                "Success, NanoStat(%s)" % (self.software_names.get_NanoStat_version()),
            )
            manage_database.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software,
                MetaKeyAndValue.META_VALUE_Success,
                result_all.to_json(),
            )

            ### remove nanofilt and rabiQC
            self.utils.remove_file(sample.get_nanofilt_file(TypePath.MEDIA_ROOT))
            self.utils.remove_file(sample.get_rabbitQC_nanofilt(TypePath.MEDIA_ROOT))
            return (result_average.has_reads(), False)

        ### run nanoflit
        try:
            ## get dynamic parameters
            if owner is None:
                parameters = self.software_names.get_NanoFilt_parameters()
            else:
                parameters = (
                    default_software_project.get_nanofilt_parameters_all_possibilities(
                        owner, sample
                    )
                )
            result_file = self.run_nanofilt(
                sample.get_fastq(TypePath.MEDIA_ROOT, True), parameters, owner
            )
            result_all.add_software(
                SoftwareDesc(
                    self.software_names.get_NanoFilt_name(),
                    self.software_names.get_NanoFilt_version(),
                    parameters,
                )
            )

            ### need to copy the files to samples/user path
            self.utils.copy_file(
                result_file, sample.get_nanofilt_file(TypePath.MEDIA_ROOT)
            )
            self.utils.remove_file(result_file)
        except Exception as e:
            result = Result()
            result.set_error("Fail to run NanoFilt software: " + e.args[0])
            result.add_software(
                SoftwareDesc(
                    self.software_names.get_NanoFilt_name(),
                    self.software_names.get_NanoFilt_version(),
                    self.software_names.get_NanoFilt_parameters(),
                )
            )
            manage_database.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_NanoStat_NanoFilt,
                MetaKeyAndValue.META_VALUE_Error,
                result.to_json(),
            )
            return False, False

        ### collect numbers
        (number_1, average_1, std1) = self.utils.get_number_sequences_fastq(
            sample.get_nanofilt_file(TypePath.MEDIA_ROOT)
        )
        result_average = ResultAverageAndNumberReads(number_1, average_1, None, None)
        manage_database.set_sample_metakey(
            sample,
            owner,
            MetaKeyAndValue.META_KEY_Number_And_Average_Reads,
            MetaKeyAndValue.META_VALUE_Success,
            result_average.to_json(),
        )

        ### run start again
        try:
            ## if the user set a high coverage it can get no reads to process
            result_nano_stat = self.run_nanostat(
                sample.get_nanofilt_file(TypePath.MEDIA_ROOT), number_1
            )
            result_all.add_software(
                SoftwareDesc(
                    self.software_names.get_NanoStat_name(),
                    self.software_names.get_NanoStat_version(),
                    self.software_names.get_NanoStat_parameters(),
                    result_nano_stat.key_values,
                )
            )

        except Exception as e:
            result = Result()
            result.set_error("Fail to run nanostat software: " + e.args[0])
            result.add_software(
                SoftwareDesc(
                    self.software_names.get_NanoStat_name(),
                    self.software_names.get_NanoStat_version(),
                    self.software_names.get_NanoStat_parameters(),
                )
            )
            manage_database.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_NanoStat_NanoFilt,
                MetaKeyAndValue.META_VALUE_Error,
                result.to_json(),
            )
            return False, False

        ### run rabbitQC, only to have a image of the data
        try:
            out_file_html = self.run_rabbitQC(
                sample.get_nanofilt_file(TypePath.MEDIA_ROOT)
            )
            result_all.add_software(
                SoftwareDesc(
                    self.software_names.get_rabbitQC_name(),
                    self.software_names.get_rabbitQC_version(),
                    self.software_names.get_rabbitQC_parameters(),
                )
            )
            self.utils.copy_file(
                out_file_html, sample.get_rabbitQC_nanofilt(TypePath.MEDIA_ROOT)
            )
            self.utils.remove_file(out_file_html)
        except Exception as e:
            result = Result()
            result.set_error("Fail to run rabbitQC software: " + e.args[0])
            result.add_software(
                SoftwareDesc(
                    self.software_names.get_rabbitQC_name(),
                    self.software_names.get_rabbitQC_version(),
                    self.software_names.get_rabbitQC_parameters(),
                )
            )
            manage_database.set_sample_metakey(
                sample,
                owner,
                MetaKeyAndValue.META_KEY_NanoStat_NanoFilt,
                MetaKeyAndValue.META_VALUE_Error,
                result.to_json(),
            )
            return False, False

        ## save everything OK
        manage_database.set_sample_metakey(
            sample,
            owner,
            MetaKeyAndValue.META_KEY_NanoStat_NanoFilt,
            MetaKeyAndValue.META_VALUE_Success,
            "Success, NanoStat(%s), NanoFilt(%s)"
            % (
                self.software_names.get_NanoStat_version(),
                self.software_names.get_NanoFilt_version(),
            ),
        )
        manage_database.set_sample_metakey(
            sample,
            owner,
            MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software,
            MetaKeyAndValue.META_VALUE_Success,
            result_all.to_json(),
        )
        return result_average.has_reads(), True

    def run_nanostat(self, file_name, number_sequences):
        """
        run nanoStat, return result with keys
        crate a text file with output
                General summary:
                Mean read length:                  490.4
                Mean read quality:                  12.6
                Median read length:                488.0
                Median read quality:                12.7
                Number of reads:               211,388.0
                Read length N50:                   489.0
                STDEV read length:                  17.2
                Total bases:               103,672,824.0
                Number, percentage and megabases of reads above quality cutoffs
                >Q5:    211388 (100.0%) 103.7Mb
                >Q7:    211388 (100.0%) 103.7Mb
                >Q10:    191968 (90.8%) 93.9Mb
                >Q12:    136331 (64.5%) 66.5Mb
                >Q15:    16830 (8.0%) 8.2Mb
        """

        if number_sequences == 0:
            result = Result()
            for key in SoftwareNames.SOFTWARE_NANOSTAT_vect_info_to_collect:
                result.add_key_value(KeyValue(key, "0"))
        else:
            temp_dir = self.utils.get_temp_dir()
            out_path_file_name = os.path.join(temp_dir, "temp.txt")
            if settings.RUN_NANOFILT_AND_NANOSTAT_IN_MEDAKA_ENV:
                cmd = "{} {} -o {} -n {} --fastq {}".format(
                    self.software_names.get_medaka_env(),
                    self.software_names.get_NanoStat(),
                    temp_dir,
                    out_path_file_name,
                    file_name,
                )
            else:
                cmd = "{} -o {} -n {} --fastq {}".format(
                    self.software_names.get_NanoStat(),
                    temp_dir,
                    out_path_file_name,
                    file_name,
                )

            exist_status = os.system(cmd)
            if exist_status != 0 or not os.path.exists(out_path_file_name):
                self.utils.remove_dir(temp_dir)
                self.logger_production.error("Fail to run: " + cmd)
                self.logger_debug.error("Fail to run: " + cmd)
                raise Exception("Fail to run run_fastq")

            ### process out STATs
            result = Result()
            with open(out_path_file_name) as handle_in:
                for line in handle_in:
                    sz_temp = line.strip()
                    if len(sz_temp) == 0 or sz_temp[0] == "#":
                        continue

                    ### add tags
                    lst_data = sz_temp.split(":")
                    if (
                        len(lst_data) == 2
                        and lst_data[0].strip()
                        in SoftwareNames.SOFTWARE_NANOSTAT_vect_info_to_collect
                    ):
                        result.add_key_value(
                            KeyValue(lst_data[0].strip(), lst_data[1].strip())
                        )

            ### remove dir
            self.utils.remove_dir(temp_dir)
        return result

    def run_rabbitQC(self, file_name):
        """
        run rabbitQC, return output directory
        only for visual data
        return a html file
        """
        temp_file = self.utils.get_temp_file("rabbit", ".html")
        cmd = "{} {} -i {} -h {}".format(
            self.software_names.get_rabbitQC(),
            self.software_names.get_rabbitQC_parameters(),
            file_name,
            temp_file,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.utils.remove_file(temp_file)
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run rabbitQC")
        return temp_file

    def run_nanofilt(self, file_name, parameters, user=None):
        """
        run nanofilt
        :out clean file and parameters
        """

        ### run software
        temp_file_name = self.utils.get_temp_file("nanofilt", "fastq.fz")
        if settings.RUN_NANOFILT_AND_NANOSTAT_IN_MEDAKA_ENV:
            cmd = "{} gzip -cd {} | {} {} | gzip > {}".format(
                self.software_names.get_medaka_env(),
                file_name,
                self.software_names.get_NanoFilt(),
                parameters,
                temp_file_name,
            )
        else:
            cmd = "gzip -cd {} | {} {} | gzip > {}".format(
                file_name,
                self.software_names.get_NanoFilt(),
                parameters,
                temp_file_name,
            )

        exist_status = os.system(cmd)
        if exist_status != 0:
            self.utils.remove_file(temp_file_name)
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            raise Exception("Fail to run NanoFilt")
        return temp_file_name

    """
    Global processing, Medaka, Coverage, and MixedInfections
    """

    def process_second_stage_medaka(self, project_sample, user):
        """
        Global processing, medaka and coverage
        """
        ### make it running
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        manageDatabase = ManageDatabase()
        result_all = Result()

        process_SGE.set_process_controler(
            user,
            process_controler.get_name_project_sample(project_sample),
            ProcessControler.FLAG_RUNNING,
        )

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

            ## process medaka
            default_project_software = DefaultProjectSoftware()
            default_project_software.test_all_defaults(user, None, project_sample, None)
            try:
                #### need to check the models
                parameters_medaka_consensus = (
                    default_project_software.get_medaka_parameters_all_possibilities(
                        user, project_sample
                    )
                )
                coverage_limit = int(
                    default_project_software.get_limit_coverage_ONT_single_parameter(
                        project_sample, DefaultParameters.MASK_CONSENSUS_threshold
                    )
                )
                freq_vcf_limit = float(
                    default_project_software.get_freq_vcf_ONT_single_parameter(
                        project_sample, DefaultParameters.MASK_CONSENSUS_threshold
                    )
                )

                ## deactivate, rigth now
                parameters_depth = default_project_software.get_samtools_parameters_all_possibilities_ONT(
                    user, project_sample
                )

                out_put_path = self.run_medaka(
                    project_sample.sample.get_fastq_available(TypePath.MEDIA_ROOT),
                    project_sample.project.reference.get_reference_fasta(
                        TypePath.MEDIA_ROOT
                    ),
                    project_sample.project.reference.get_reference_gbk(
                        TypePath.MEDIA_ROOT
                    ),
                    project_sample.sample.name,
                    parameters_medaka_consensus,
                    parameters_depth,
                    coverage_limit,
                    freq_vcf_limit,
                    project_sample,
                )
                result_all.add_software(
                    SoftwareDesc(
                        self.software_names.get_medaka_name(),
                        self.software_names.get_medaka_version(),
                        "consensus " + parameters_medaka_consensus,
                    )
                )
                result_all.add_software(
                    SoftwareDesc(
                        self.software_names.get_samtools_name(),
                        self.software_names.get_samtools_version(),
                        "depth {}".format(parameters_depth),
                    )
                )
                result_all.add_software(
                    SoftwareDesc(
                        self.software_names.get_medaka_name(),
                        self.software_names.get_medaka_version(),
                        "variant --verbose",
                    )
                )
                result_all.add_software(
                    SoftwareDesc(
                        self.software_names.get_insaflu_parameter_limit_coverage_name(),
                        "",
                        "Threshold:{}".format(coverage_limit),
                    )
                )
                result_all.add_software(
                    SoftwareDesc(
                        self.software_names.get_insaflu_parameter_freq_vcf_name(),
                        "",
                        "Threshold:{}".format(freq_vcf_limit),
                    )
                )
                result_all.add_software(
                    SoftwareDesc(
                        self.software_names.get_bcftools_name(),
                        self.software_names.get_bamtools_version(),
                        "consensus generation from medaka filtered vcf",
                    )
                )
            except Exception as e:
                result = Result()
                result.set_error(e.args[0])
                result.add_software(
                    SoftwareDesc(
                        self.software_names.get_medaka_name(),
                        self.software_names.get_medaka_version(),
                        self.software_names.get_medaka_parameters_consensus(),
                    )
                )
                manageDatabase.set_project_sample_metakey(
                    project_sample,
                    user,
                    MetaKeyAndValue.META_KEY_Medaka,
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
                if not meta_sample is None:
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

            ### Also copy the file "{}.consensus.fa".format(sample_name)
            ## copy the files to the project sample directories
            self.software.copy_files_to_project(
                project_sample, self.software_names.get_medaka_name(), out_put_path
            )
            self.utils.remove_dir(out_put_path)

            ### make the link for the new tab file name
            path_medaka_tab = project_sample.get_file_output(
                TypePath.MEDIA_ROOT,
                FileType.FILE_TAB,
                self.software_names.get_medaka_name(),
            )
            if os.path.exists(path_medaka_tab):
                sz_file_to = project_sample.get_file_output_human(
                    TypePath.MEDIA_ROOT,
                    FileType.FILE_TAB,
                    self.software_names.get_medaka_name(),
                )
                self.utils.link_file(path_medaka_tab, sz_file_to)

            ## get coverage from deep file
            get_coverage = GetCoverage()
            try:

                ### limit of the coverage for a project, can be None, if not exist
                b_coverage_default = (
                    False  ## because in ONT the limit is different than10
                )
                coverage_for_project = int(
                    default_project_software.get_limit_coverage_ONT_single_parameter(
                        project_sample, DefaultParameters.MASK_CONSENSUS_threshold
                    )
                )
                default_coverage_value = int(
                    default_project_software.get_limit_coverage_ONT_single_parameter(
                        project_sample, DefaultParameters.MASK_CONSENSUS_threshold
                    )
                )

                coverage = get_coverage.get_coverage(
                    project_sample.get_file_output(
                        TypePath.MEDIA_ROOT,
                        FileType.FILE_DEPTH_GZ,
                        self.software_names.get_medaka_name(),
                    ),
                    project_sample.project.reference.get_reference_fasta(
                        TypePath.MEDIA_ROOT
                    ),
                    default_coverage_value,
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
                    ConstantsSettings.TECHNOLOGY_minion,
                )
            )
            msa_parameters = self.software.make_mask_consensus_by_deep(
                project_sample.get_file_output(
                    TypePath.MEDIA_ROOT,
                    FileType.FILE_CONSENSUS_FASTA,
                    self.software_names.get_medaka_name(),
                ),
                project_sample.project.reference.get_reference_fasta(
                    TypePath.MEDIA_ROOT
                ),
                project_sample.get_file_output(
                    TypePath.MEDIA_ROOT,
                    FileType.FILE_DEPTH_GZ,
                    self.software_names.get_medaka_name(),
                ),
                coverage,
                project_sample.sample.name,
                limit_to_mask_consensus,
            )
            result_all.add_software(
                SoftwareDesc(
                    self.software_names.get_msa_masker_name(),
                    self.software_names.get_msa_masker_version(),
                    "{}; for coverages less than {} in {}% of the regions.".format(
                        msa_parameters, coverage_limit, 100 - limit_to_mask_consensus
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
                    self.software_names.get_medaka_name(),
                ),
                coverage,
            )

            ### draw coverage
            try:
                ### make the coverage images
                draw_all_coverage = DrawAllCoverage()
                draw_all_coverage.draw_all_coverages(
                    project_sample, SoftwareNames.SOFTWARE_Medaka_name
                )
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

            ## count hits from tab file
            count_hits = CountHits()
            if not out_put_path is None:
                file_tab = project_sample.get_file_output(
                    TypePath.MEDIA_ROOT,
                    FileType.FILE_TAB,
                    self.software_names.get_medaka_name(),
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
                result.set_error("Fail to calculate mixed infection")
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

            ### get clean consensus file
            consensus_fasta = project_sample.get_file_output(
                TypePath.MEDIA_ROOT,
                FileType.FILE_CONSENSUS_FASTA,
                SoftwareNames.SOFTWARE_Medaka_name,
            )
            if os.path.exists(consensus_fasta):
                file_out = project_sample.get_consensus_file(
                    TypePath.MEDIA_ROOT
                )  ### this is going to main path
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
                MetaKeyAndValue.META_KEY_Medaka,
                MetaKeyAndValue.META_VALUE_Success,
                result_all.to_json(),
            )

            ### set the flag of the end of the task
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
                    MetaKeyAndValue.META_VALUE_Success,
                    meta_sample.description,
                )
        except Exception as e:
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

    def run_medaka(
        self,
        file_fastq,
        reference_fasta,
        reference_gbk,
        sample_name,
        parameters_consensus,
        parameters_depth,
        coverage_limit,
        freq_vcf_limit,
        project_sample,
    ):
        """
                run medaka
                return output directory of snippy, try to do something most close possible


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

                :param
                (default: r941_min_high_g360).
        Available: r103_min_high_g345, r103_min_high_g360, r103_prom_high_g360, r103_prom_snp_g3210, r103_prom_variant_g3210,
                r10_min_high_g303, r10_min_high_g340, r941_min_fast_g303, r941_min_high_g303, r941_min_high_g330,
                r941_min_high_g340_rle, r941_min_high_g344, r941_min_high_g351, r941_min_high_g360, r941_prom_fast_g303,
                r941_prom_high_g303, r941_prom_high_g330, r941_prom_high_g344, r941_prom_high_g360, r941_prom_high_g4011,
                r941_prom_snp_g303, r941_prom_snp_g322, r941_prom_snp_g360, r941_prom_variant_g303, r941_prom_variant_g322,
                r941_prom_variant_g360.
                ### output medaka files

                bcftools   1.9
                bgzip      1.9
                minimap2   2.11
                samtools   1.9
                tabix      1.9

                consensus.fasta
                calls_to_draft.bam
                :param coverage_limit cut all variants less than this coverage
        """
        ### make medaka consensus and bam
        temp_dir = os.path.join(self.utils.get_temp_dir(), sample_name)

        ###need to make a link to the reference files
        ### Always need to run consensus because do HDF file for variants
        reference_fasta_medaka = self.utils.get_temp_file_from_dir(
            temp_dir, "medaka_ref", ".fasta"
        )
        self.utils.copy_file(reference_fasta, reference_fasta_medaka)

        cmd = "{} {}_consensus -i {} -d {} -o {} -t {} {}".format(
            self.software_names.get_medaka_env(),
            self.software_names.get_medaka(),
            file_fastq,
            reference_fasta_medaka,
            temp_dir,
            settings.THREADS_TO_RUN_SLOW,
            parameters_consensus,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            self.utils.remove_dir(temp_dir)
            raise Exception("Fail to run medaka_consensus")

        ### test output files
        hdf_file = os.path.join(temp_dir, "consensus_probs.hdf")
        bam_file = os.path.join(temp_dir, "calls_to_draft.bam")
        consensus_file = os.path.join(temp_dir, "consensus.fasta")
        for file_to_test in [hdf_file, bam_file, consensus_file]:
            if not os.path.exists(file_to_test):
                message = (
                    "File '{}' not found after medaka_consensus process: ".format(
                        file_to_test
                    )
                    + cmd
                )
                self.logger_production.error(message)
                self.logger_debug.error(message)
                self.utils.remove_dir(temp_dir)
                raise Exception(message)

        ### change bam file names
        self.utils.move_file(
            bam_file, os.path.join(os.path.dirname(bam_file), sample_name + ".bam")
        )
        self.utils.move_file(
            bam_file + ".bai",
            os.path.join(os.path.dirname(bam_file), sample_name + ".bam.bai"),
        )
        ## get from BCFtools
        ## self.utils.move_file(consensus_file, os.path.join(os.path.dirname(bam_file), sample_name + ".consensus.fa"))
        bam_file = os.path.join(os.path.dirname(bam_file), sample_name + ".bam")

        ### get mapped stast reads
        result = Result()
        if os.path.exists(bam_file):
            result = self.software.get_statistics_bam(bam_file)
        manageDatabase = ManageDatabase()
        manageDatabase.set_project_sample_metakey(
            project_sample,
            project_sample.sample.owner,
            MetaKeyAndValue.META_KEY_bam_stats,
            MetaKeyAndValue.META_VALUE_Success,
            result.to_json(),
        )

        ### create depth
        depth_file = os.path.join(temp_dir, sample_name + ".depth.gz")
        ##        cmd =  "{} depth -aa -q 10 {} | {} -c > {}".format(   ### with quality
        cmd = "{} depth {} {} | {} -c > {}".format(
            self.software_names.get_samtools(),
            parameters_depth,
            bam_file,
            self.software_names.get_bgzip(),
            depth_file,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            self.utils.remove_dir(temp_dir)
            raise Exception("Fail to run samtools depth in nanopore")

        ### create tbi
        self.software.create_index_with_tabix(
            depth_file, SoftwareNames.SOFTWARE_DEPTH_SAMTOOLS_file_flag
        )

        ### vcf
        vcf_before_file = os.path.join(temp_dir, sample_name + "_before_annotation.vcf")
        cmd = "{} {} variant --verbose {} {} {};".format(
            #         cmd =  "{} {} snp --verbose {} {} {}".format(
            self.software_names.get_medaka_env(),
            self.software_names.get_medaka(),
            reference_fasta_medaka,
            hdf_file,
            vcf_before_file,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            self.utils.remove_dir(temp_dir)
            raise Exception("Fail to run medaka variant")

        ### annotate vcf
        vcf_file = os.path.join(temp_dir, sample_name + ".vcf")
        cmd = "{} {} tools annotate {} {} {} {}".format(
            self.software_names.get_medaka_env(),
            self.software_names.get_medaka(),
            vcf_before_file,
            reference_fasta_medaka,
            bam_file,
            vcf_file,
        )
        exist_status = os.system(cmd)
        if exist_status != 0:
            self.logger_production.error("Fail to run: " + cmd)
            self.logger_debug.error("Fail to run: " + cmd)
            self.utils.remove_dir(temp_dir)
            raise Exception("Fail to run medaka annotate")

        ### Add
        ##INFO=<ID=SR,Number=.,Type=Integer,Description="Depth of spanning reads by strand which best align to each allele (ref fwd, ref rev, alt1 fwd, alt1 rev, etc.)">
        ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
        ##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
        ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Number of observation for each allele">

        if os.path.getsize(vcf_file) > 0:
            ### Add ANN anotation with snpEff
            temp_file_2 = os.path.join(temp_dir, sample_name + "_ann.vcf")
            output_file = self.software.run_snpEff(
                reference_fasta_medaka, reference_gbk, vcf_file, temp_file_2
            )

            if (
                not output_file is None
            ):  ## sometimes the gff does not have amino sequences
                self.utils.copy_file(output_file, vcf_file)

            ### add FREQ to VCF file
            final_vcf = os.path.join(temp_dir, sample_name + "_2.vcf")
            final_vcf_with_removed_variants = os.path.join(
                temp_dir, sample_name + "_removed.vcf"
            )
            self.utils.add_freq_ao_ad_and_type_to_vcf(
                vcf_file,
                depth_file,
                final_vcf,
                final_vcf_with_removed_variants,
                coverage_limit,
                freq_vcf_limit,
            )
            self.utils.move_file(final_vcf, vcf_file)
            self.software.test_bgzip_and_tbi_in_vcf(vcf_file)

            ### IF YOU mask CONSENSUS the positions on VCF are not real for consensus
            ### mask REFERENCE for variants below minfrac, with vcf_removed_variants,
            manageDatabase = ManageDatabase()
            decode_results = DecodeObjects()
            meta_value = manageDatabase.get_reference_metakey_last(
                project_sample.project.reference,
                MetaKeyAndValue.META_KEY_Elements_And_CDS_Reference,
                MetaKeyAndValue.META_VALUE_Success,
            )
            if not meta_value is None:
                genetic_element = decode_results.decode_result(meta_value.description)

                vcf_hanlder = pysam.VariantFile(final_vcf_with_removed_variants, "r")
                mask_consensus = MaskingConsensus()
                element_name_old = ""
                vect_sites = []
                vect_ranges = []
                for variant in vcf_hanlder:
                    if element_name_old != variant.chrom:
                        if len(element_name_old) > 0:
                            mask_consensus = MaskingConsensus()
                            mask_consensus.set_mask_sites(",".join(vect_sites))
                            mask_consensus.set_mask_regions(",".join(vect_ranges))
                            genetic_element.set_mask_consensus_element(
                                element_name_old, mask_consensus
                            )

                        ## new one
                        element_name_old = variant.chrom
                        vect_sites = []
                        vect_ranges = []
                    ### MEDAKA output must have "TYPE" in info
                    if variant.info["TYPE"][0] == "snp":
                        vect_sites.append(str(variant.pos))
                    elif variant.info["TYPE"][0] == "ins":
                        vect_sites.append(str(variant.pos))
                    elif variant.info["TYPE"][0] == "del":
                        vect_ranges.append(
                            "{}-{}".format(
                                variant.pos + len(variant.alts[0]),
                                variant.pos - len(variant.alts[0]) + len(variant.ref),
                            )
                        )
                    else:
                        vect_ranges.append(
                            "{}-{}".format(
                                variant.pos, variant.pos + len(variant.ref) - 1
                            )
                        )

                ## last one
                if len(element_name_old) > 0:
                    mask_consensus = MaskingConsensus()
                    mask_consensus.set_mask_sites(",".join(vect_sites))
                    mask_consensus.set_mask_regions(",".join(vect_ranges))
                    mask_consensus.cleaning_mask_results()
                    genetic_element.set_mask_consensus_element(
                        element_name_old, mask_consensus
                    )
                ## mask
                temp_reference = os.path.join(temp_dir, sample_name + "_consensus.vcf")
                self.utils.copy_file(reference_fasta_medaka, temp_reference)
                self.utils.mask_sequence_by_sites(
                    temp_reference, temp_reference, genetic_element
                )

                ### save positions that are going to be masked by MinFrac, even if there's not...
                manageDatabase.set_project_sample_metakey(
                    project_sample,
                    project_sample.project.owner,
                    MetaKeyAndValue.META_KEY_Masking_consensus_by_minfrac_VCF_medaka,
                    MetaKeyAndValue.META_VALUE_Success,
                    genetic_element.to_json(),
                )

            ####################
            ### run BCF tools to get the consensus
            if not os.path.exists(vcf_file + ".gz"):
                self.logger_production.error("Fail to find: " + vcf_file + ".gz")
                self.logger_debug.error("Fail to find: " + vcf_file + ".gz")
                self.utils.remove_dir(temp_dir)
                raise Exception("Fail to find VCF file: " + vcf_file + ".gz")

            ## self.utils.move_file(consensus_file, os.path.join(os.path.dirname(bam_file), sample_name + ".consensus.fa"))
            cmd = "{} consensus -s SAMPLE -f {} {} -o {}".format(
                self.software_names.get_bcftools(),
                temp_reference,
                vcf_file + ".gz",
                os.path.join(os.path.dirname(bam_file), sample_name + ".consensus.fa"),
            )
            exist_status = os.system(cmd)
            if exist_status != 0:
                self.logger_production.error("Fail to run: " + cmd)
                self.logger_debug.error("Fail to run: " + cmd)
                self.utils.remove_dir(temp_dir)
                raise Exception("Fail to run bcftools consensus")

            ########################
            ### create TAB file
            self.software.run_snippy_vcf_to_tab_freq_and_evidence(
                reference_fasta_medaka,
                reference_gbk,
                vcf_file,
                os.path.join(temp_dir, sample_name + ".tab"),
            )
        return temp_dir
