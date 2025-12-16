"""
Created on 12/06/2022

@author: mmp
"""

import csv
import glob
import json
import logging
import os
import time

from django.conf import settings

from constants.constants import Constants, FileExtensions, TypeFile, TypePath
from constants.meta_key_and_values import MetaKeyAndValue
from constants.software_names import SoftwareNames
from datasets.manage_database import ManageDatabase
from datasets.models import Dataset, DatasetConsensus, UploadFiles
from managing_files.models import ProcessControler, Project, Reference
from settings.models import Parameter
from utils.data_columns import DataColumns
from utils.exceptions import CmdException
from utils.parse_in_files_nextstrain import ParseNextStrainFiles
from utils.parse_out_files import ParseOutFiles
from utils.process_SGE import ProcessSGE
from utils.software import Software
from utils.utils import Utils


class CollectExtraDatasetData(object):
    """
    classdocs
    """

    ## Type of sample list
    DATASET_LIST_simple = 0  ## used in the tree
    DATASET_LIST_list = 1  ## list of samples/references/consensus with metadata

    NEXTSTRAIN_need_key = "length"  ## Consensus length of the sample

    utils = Utils()
    software = Software()
    if settings.DEBUG:
        logger = logging.getLogger("fluWebVirus.debug")
    else:
        logger = logging.getLogger("fluWebVirus.production")

    def __init__(self):
        """
        Constructor
        """
        pass

    def collect_extra_data_for_dataset(self, dataset, user):
        """ """

        dataset.is_processed = False
        dataset.save()

        ### make it running
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        process_SGE.set_process_controler(
            user,
            process_controler.get_name_dataset(dataset),
            ProcessControler.FLAG_RUNNING,
        )

        ### set user globally
        self.user = user

        ## need to add a delay for the test in command line
        if settings.RUN_TEST_IN_COMMAND_LINE:
            time.sleep(4)

        ## run collect data
        self.__collect_extra_data_for_dataset(dataset, user)

        dataset.is_processed = True
        dataset.save()

    def collect_update_extra_metadata_for_dataset(self, dataset, user):
        """
        Only for update metadata
        """

        dataset.is_processed = False
        dataset.save()

        ### make it running
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        process_SGE.set_process_controler(
            user,
            process_controler.get_name_dataset(dataset),
            ProcessControler.FLAG_RUNNING,
        )

        ## need to add a delay for the test in command line
        if settings.RUN_TEST_IN_COMMAND_LINE:
            time.sleep(4)

        ## run collect data
        self.__collect_update_extra_metadata_for_dataset(dataset, user)

        dataset.is_processed = True
        dataset.save()

    def __collect_update_extra_metadata_for_dataset(self, dataset, user):
        """
        Only for update metadata
        It is not necessary
        """
        ### get the taskID and seal it
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        manage_database = ManageDatabase()
        metaKeyAndValue = MetaKeyAndValue()

        try:

            manage_database.get_max_length_label(dataset, user, True)
            self.calculate_global_files(
                Dataset.DATASET_FILE_NAME_RESULT_all_consensus, dataset
            )

            ## collect sample table with plus type and subtype, mixed infection, equal to upload table
            self.calculate_global_files(Dataset.DATASET_FILE_NAME_RESULT_CSV, dataset)
            self.calculate_global_files(Dataset.DATASET_FILE_NAME_RESULT_TSV, dataset)

            self.calculate_global_files(
                Dataset.DATASET_FILE_NAME_RESULT_NEXTSTRAIN_TSV, dataset
            )
            self.calculate_global_files(
                Dataset.DATASET_FILE_NAME_RESULT_NEXTSTRAIN_CSV, dataset
            )

            ## Important, this need to be after DATASET_FILE_NAME_RESULT_NEXTSTRAIN_CSV
            self.calculate_global_files(Dataset.DATASET_FILE_NAME_RESULT_json, dataset)

            ### zip several files to download
            self.zip_several_files(dataset)

        except Exception as e:
            ## finished with error
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_dataset(dataset),
                ProcessControler.FLAG_ERROR,
            )
            self.logger.error(e)
            return

        ## seal the tag
        meta_dataset = manage_database.get_dataset_metakey_last(
            dataset,
            metaKeyAndValue.get_meta_key(
                MetaKeyAndValue.META_KEY_Queue_TaskID_Project, dataset.id
            ),
            MetaKeyAndValue.META_VALUE_Queue,
        )
        if meta_dataset != None:
            manage_database.set_dataset_metakey(
                dataset,
                user,
                metaKeyAndValue.get_meta_key(
                    MetaKeyAndValue.META_KEY_Queue_TaskID_Project, dataset.id
                ),
                MetaKeyAndValue.META_VALUE_Success,
                meta_dataset.description,
            )

        ### finished
        process_SGE.set_process_controler(
            user,
            process_controler.get_name_dataset(dataset),
            ProcessControler.FLAG_FINISHED,
        )

    def __collect_extra_data_for_dataset(self, dataset, user):
        """
        Everything that is necessary to do in the dataset
        Collect all extra data after all samples are finished
        only run after all the vect_taskID are finished
        """
        ### get the taskID and seal it
        metaKeyAndValue = MetaKeyAndValue()
        manage_database = ManageDatabase()
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        try:
            count = 0
            start = time.time()
            self.logger.info("COLLECT_EXTRA_FILES: Start")
            ## calculate the max sample label size of the samples that belong to this dataset
            ## used in MSA viewer
            # b_calculate_again = True
            # manage_database.get_max_length_label(dataset, user, b_calculate_again)
            # self.logger.info("COLLECT_EXTRA_FILES: Step {}  diff_time:{}".format(count, time.time() - start))
            # count += 1
            # start = time.time()

            ## collect all consensus files for a dataset
            self.calculate_global_files(
                Dataset.DATASET_FILE_NAME_RESULT_all_consensus, dataset
            )
            self.logger.info(
                "COLLECT_EXTRA_FILES: Step {}  diff_time:{}".format(
                    count, time.time() - start
                )
            )
            count += 1
            start = time.time()
            ## collect sample table with plus type and subtype, mixed infection, equal to upload table
            self.calculate_global_files(Dataset.DATASET_FILE_NAME_RESULT_CSV, dataset)
            self.calculate_global_files(Dataset.DATASET_FILE_NAME_RESULT_TSV, dataset)            
            self.calculate_global_files(
                Dataset.DATASET_FILE_NAME_RESULT_NEXTSTRAIN_TSV, dataset
            )          
            self.calculate_global_files(
                Dataset.DATASET_FILE_NAME_RESULT_NEXTSTRAIN_CSV, dataset
            )
            ## Important, this need to be after DATASET_FILE_NAME_RESULT_NEXTSTRAIN_CSV
            self.calculate_global_files(Dataset.DATASET_FILE_NAME_RESULT_json, dataset)
            self.logger.info(
                "COLLECT_EXTRA_FILES: Step {}  diff_time:{}".format(
                    count, time.time() - start
                )
            )
            count += 1
            start = time.time()
            # self.calculate_global_files(Dataset.DATASET_FILE_NAME_nextstrain_default_build, dataset)
            self.calculate_global_files(
                Dataset.DATASET_FILE_NAME_nextstrain_auspice_zip, dataset
            )
            self.logger.info(
                "RUN nextStrain: Step {}  diff_time:{}".format(
                    count, time.time() - start
                )
            )
            count += 1

            ### create trees; this is now replaced by nextstrain results
            # createTree = CreateTree()
            # createTree.create_tree_and_alignments_dataset(dataset, user)
            # self.logger.info("Trees and alignments: Step {}  diff_time:{}".format(count, time.time() - start))
            # count += 1
            # start = time.time()

            meta_dataset = manage_database.get_dataset_metakey_last(
                dataset,
                metaKeyAndValue.get_meta_key(
                    MetaKeyAndValue.META_KEY_Queue_TaskID_Project, dataset.id
                ),
                MetaKeyAndValue.META_VALUE_Queue,
            )
            if meta_dataset != None:
                manage_database.set_dataset_metakey(
                    dataset,
                    user,
                    metaKeyAndValue.get_meta_key(
                        MetaKeyAndValue.META_KEY_Queue_TaskID_Project, dataset.id
                    ),
                    MetaKeyAndValue.META_VALUE_Success,
                    meta_dataset.description,
                )

            ### zip several files to download
            self.zip_several_files(dataset)
            self.logger.info(
                "COLLECT_EXTRA_FILES: Step {}  diff_time:{}".format(
                    count, time.time() - start
                )
            )
            count += 1
            start = time.time()

        except Exception as e:
            ## finished with error
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_dataset(dataset),
                ProcessControler.FLAG_ERROR,
            )
            return

        ### finished
        process_SGE.set_process_controler(
            user,
            process_controler.get_name_dataset(dataset),
            ProcessControler.FLAG_FINISHED,
        )

    def calculate_global_files(self, type_file, dataset):
        """
        Collect extra files
        """

        # check what is the build that is configured, otherwise use the default
        build = SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_parameter

        # See if there is a build parameter specific for this dataset, in which case use it
        parameters_list = Parameter.objects.filter(dataset=dataset)
        if len(list(parameters_list)) == 1:
            build = list(parameters_list)[0].parameter
        out_file = None
        out_file_file_system = None
        if type_file == Dataset.DATASET_FILE_NAME_RESULT_CSV:
            ## samples csv
            out_file = self.collect_sample_table(
                dataset,
                Constants.SEPARATOR_COMMA,
                CollectExtraDatasetData.DATASET_LIST_list,
                build,
            )
            out_file_file_system = dataset.get_global_file_by_dataset(
                TypePath.MEDIA_ROOT, type_file
            )
        elif type_file == Dataset.DATASET_FILE_NAME_RESULT_TSV:
            ## samples tsv
            out_file = self.collect_sample_table(
                dataset,
                Constants.SEPARATOR_TAB,
                CollectExtraDatasetData.DATASET_LIST_list,
                build,
            )
            out_file_file_system = dataset.get_global_file_by_dataset(
                TypePath.MEDIA_ROOT, type_file
            )
        elif type_file == Dataset.DATASET_FILE_NAME_RESULT_NEXTSTRAIN_TSV:
            ## samples tsv
            out_file = self.collect_nextstrain_table(
                dataset, Constants.SEPARATOR_TAB, build
            )
            out_file_file_system = dataset.get_global_file_by_dataset(
                TypePath.MEDIA_ROOT, type_file
            )
        elif type_file == Dataset.DATASET_FILE_NAME_RESULT_NEXTSTRAIN_CSV:
            ## samples csv
            out_file = self.collect_nextstrain_table(
                dataset, Constants.SEPARATOR_COMMA, build
            )
            out_file_file_system = dataset.get_global_file_by_dataset(
                TypePath.MEDIA_ROOT, type_file
            )
        elif type_file == Dataset.DATASET_FILE_NAME_RESULT_all_consensus:
            out_file = self.merge_all_consensus_files(dataset)
            out_file_file_system = dataset.get_global_file_by_dataset(
                TypePath.MEDIA_ROOT, type_file
            )
        elif type_file == Dataset.DATASET_FILE_NAME_RESULT_json:
            ## tree json
            out_file = self.create_json_file_from_sample_csv(dataset)
            out_file_file_system = dataset.get_global_file_by_dataset(
                TypePath.MEDIA_ROOT, type_file
            )
        # elif (type_file == Dataset.DATASET_FILE_NAME_nextstrain_default_build):
        elif type_file == Dataset.DATASET_FILE_NAME_nextstrain_auspice_zip:

            ## remove previous nextStrain
            # for type_file in Dataset.VECT_files_next_strain + [Dataset.DATASET_FILE_NAME_nextstrain_error]:
            #    out_file_file_system = dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, type_file)
            #    self.utils.remove_file(out_file_file_system)
            self.utils.remove_file(
                dataset.get_global_file_by_dataset(
                    TypePath.MEDIA_ROOT,
                    Dataset.DATASET_FILE_NAME_nextstrain_auspice_zip,
                )
            )

            auspice_zip_file = None
            try:

                tree_file, alignment_file, auspice_zip_file = self.run_nextstrain(
                    dataset, build
                )

                ## copy files if they exist, try to remove in destination
                # for type_file in Dataset.VECT_files_next_strain:
                #    out_file_file_system = dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, type_file)
                #    if os.path.exists(auspice_zip_file): self.utils.move_file(auspice_zip_file, out_file_file_system)
                #    elif (not out_file_file_system is None and os.path.exists(out_file_file_system)):
                #        self.utils.remove_file(out_file_file_system)

                out_file_file_system_tree1 = dataset.get_global_file_by_dataset(
                    TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_FASTTREE
                )
                if os.path.exists(tree_file):
                    self.utils.move_file(tree_file, out_file_file_system_tree1)

                out_file_file_system_tree2 = dataset.get_global_file_by_dataset(
                    TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_FASTTREE_tree
                )
                if os.path.exists(out_file_file_system_tree1):
                    self.utils.copy_file(
                        out_file_file_system_tree1, out_file_file_system_tree2
                    )

                out_file_file_system_alignment = dataset.get_global_file_by_dataset(
                    TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_MAFFT
                )
                if os.path.exists(alignment_file):
                    self.utils.move_file(alignment_file, out_file_file_system_alignment)

                out_file_file_system_auspice = dataset.get_global_file_by_dataset(
                    TypePath.MEDIA_ROOT,
                    Dataset.DATASET_FILE_NAME_nextstrain_auspice_zip,
                )
                if os.path.exists(auspice_zip_file):
                    self.utils.move_file(auspice_zip_file, out_file_file_system_auspice)

                ### remove possible error of previous run
                out_file_file_system = dataset.get_global_file_by_dataset(
                    TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_nextstrain_error
                )
                self.utils.remove_file(out_file_file_system)

            except CmdException as e:  ## copy the snakeMake file
                # this means there is potentially a lot of stuff left in the tmp directory that needs to be manually removed
                files = []
                if e.exist_path():
                    files = glob.glob(
                        os.path.join(e.output_path, ".snakemake/log/*snakemake.log")
                    )
                out_file_file_system = dataset.get_global_file_by_dataset(
                    TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_nextstrain_error
                )
                if len(files) > 0:
                    self.utils.copy_file(str(files[0]), out_file_file_system)
                else:
                    with open(out_file_file_system, "w") as handle_write:
                        handle_write.write(str(e))
            except Exception as e:  ## fail to run nextStrain;
                out_file_file_system = dataset.get_global_file_by_dataset(
                    TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_nextstrain_error
                )
                with open(out_file_file_system, "w") as handle_write:
                    handle_write.write(str(e))

            out_file = None  ## not copy anything with this variable
            out_file_file_system = None
            # if not temp_dir is None: self.utils.remove_dir(temp_dir)

        ## copy file
        if not out_file is None:
            self.utils.copy_file(out_file, out_file_file_system)
            self.utils.remove_file(out_file)
        elif not out_file_file_system is None and os.path.exists(out_file_file_system):
            self.utils.remove_file(out_file_file_system)

    def create_json_file_from_sample_csv(self, dataset):
        """
        Create JSON file to insaPhylo
        """
        vect_remove_keys = [
            "strain",
            "fastq1",
            "fastq2",
            "data set",
            "latitude",
            "longitude",
        ]
        ## create it from DATASET_FILE_NAME_RESULT_NEXTSTRAIN_CSV instead of DATASET_FILE_NAME_RESULT_CSV
        file_name_root_sample = dataset.get_global_file_by_dataset(
            TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_RESULT_NEXTSTRAIN_CSV
        )
        if os.path.exists(file_name_root_sample):
            out_file = self.utils.get_temp_file(
                "json_sample_file", FileExtensions.FILE_JSON
            )
            with open(out_file, "w", encoding="utf-8") as handle_write, open(
                file_name_root_sample
            ) as handle_in_csv:
                reader = csv.DictReader(handle_in_csv)
                all_data = json.loads(json.dumps(list(reader)))
                dt_result = {}
                for dict_data in all_data:
                    ##if ('id' in dict_data):
                    if "strain" in dict_data:  ## because of
                        dt_out = dict_data.copy()
                        for key_to_remove in vect_remove_keys:
                            try:
                                del dt_out[key_to_remove]
                            except KeyError:
                                pass
                        dt_result[dict_data["strain"]] = {
                            key: "" if dt_out[key] == "?" else dt_out[key]
                            for key in dt_out
                        }
                if len(dt_result) == len(all_data):
                    handle_write.write(json.dumps(dt_result))
                else:
                    os.unlink(out_file)
                    self.logger.error(
                        "ProjectID: {}  different number of lines processing Sample {} -> JSON {}".format(
                            dataset.id, len(dt_result), len(all_data)
                        )
                    )
                    return None
            return out_file
        else:
            self.logger.error(
                "Sample csv file does not exist: {}".format(file_name_root_sample)
            )
        return None

    def merge_all_consensus_files(self, dataset):
        """
        merge all consensus files
        """

        out_file = self.utils.get_temp_file("all_consensus", FileExtensions.FILE_FASTA)
        vect_to_process = []
        for dataset_consensus in dataset.dataset_consensus.all():
            if dataset_consensus.is_deleted:
                continue
            if dataset_consensus.is_error:
                continue

            if not dataset_consensus.is_ready_to_proccess():
                continue
            if not os.path.exists(
                dataset_consensus.get_consensus_file(TypePath.MEDIA_ROOT)
            ):
                continue

            vect_to_process.append(
                [
                    dataset_consensus.get_consensus_file(TypePath.MEDIA_ROOT),
                    dataset_consensus.get_name(),
                    dataset_consensus.pk,
                ]
            )  ## need to pass ID to save possible new name

        ### set the number sequences that passed
        dataset.number_passed_sequences = len(vect_to_process)
        dataset.save()

        # check what is the build that is configured, otherwise use the default
        build = SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_parameter
        # See if there is a build parameter specific for this dataset, in which case use it
        parameters_list = Parameter.objects.filter(dataset=dataset)
        if len(list(parameters_list)) == 1:
            build = list(parameters_list)[0].parameter

        temp_out_fasta = self.utils.get_temp_file(
            "temp_consensus_fasta", FileExtensions.FILE_FASTA
        )

        # These builds are not to be filtered by Abricate (for now)
        builds_to_ignore = [
            SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_generic,
            SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_generic_time,
            SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_mpx,
            SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_mpox_clade_i,
            SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_mpox_mpxv,            
            SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_ncov,
            SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_rsv_a,
            SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_rsv_b,
            SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_dengue_all,
            SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_dengue_denv1,
            SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_dengue_denv2,
            SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_dengue_denv3,
            SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_dengue_denv4,
        ]

        build_to_ignore = False
        if build in builds_to_ignore:
            build_to_ignore = True

        # Concatenates all segments together to potentially filter by abricate
        # obtains a dictionary linking fasta record.id to dataset_consensus id
        # Can have more than one record.id to the same dataset_consensus id
        possiblename_to_id = self.utils.merge_fasta_files_dataset(
            vect_to_process, temp_out_fasta
        )

        # List of final segments that will go to nextstrain
        # Maximum of one per dataset_consensus id (but can be 0)
        sample_list = []

        if build_to_ignore:

            # just make sure there is only one per dataset consensus
            # for now randomly chose one (just the first one will do)
            already_done_dataset_consensus = []
            for sample in possiblename_to_id.keys():
                if possiblename_to_id[sample] not in already_done_dataset_consensus:
                    sample_list.append(sample)
                    already_done_dataset_consensus.append(possiblename_to_id[sample])

        else:

            # Filter by abricate
            temp_out_abricate = self.utils.get_temp_file(
                "temp_abricate", FileExtensions.FILE_TXT
            )

            try:
                _ = self.software.run_abricate(
                    "nextstrain",  # Need to change this
                    temp_out_fasta,
                    SoftwareNames.SOFTWARE_ABRICATE_PARAMETERS,
                    temp_out_abricate,
                )
                parseOutFiles = ParseOutFiles()
                dict_data_out = parseOutFiles.parse_abricate_file_simple(
                    temp_out_abricate
                )
                # This doesn't really matter
                self.utils.remove_temp_file(temp_out_abricate)

                # for each dataset consensus, chose segment with best match to desired target
                best_matches = {}
                for sample in dict_data_out.keys():
                    for entry_dic in dict_data_out[sample]:
                        if entry_dic["Gene"] == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_REF_DICT[build]:
                            if sample in best_matches.keys():
                                if float(entry_dic["Identity"]) > best_matches[sample]:
                                    best_matches[sample] = float(entry_dic["Identity"])
                            else:
                                best_matches[sample] = float(entry_dic["Identity"])
                sample_list = best_matches.keys()
            except Exception as e:
                self.logger.error("Welp... Could not run abricate")
                self.logger.error(e)

        # Since there can be only one sequence per sample, get back sample names
        # unless there are repeated names, in which case, add extra number...
        # TODO find another way to avoid name conflicts...
        rename = {}
        used_names = {}
        for sample in sample_list:
            try:
                dataset_consensus = DatasetConsensus.objects.get(
                    id=possiblename_to_id[sample]
                )
                # Sample names can only contain letters, numbers and underscores so it should be ok
                # possible_name = "{}_{}".format(dataset_consensus.get_name(),dataset_consensus.get_project_name())
                possible_name = dataset_consensus.get_name()
                count = 1
                while True:
                    if possible_name in used_names:
                        possible_name = "{}_{}".format(
                            dataset_consensus.get_name(), count
                        )
                        count += 1
                    else:
                        break
                rename[sample] = possible_name
                used_names[possible_name] = 1
            except DatasetConsensus.DoesNotExist:  ## need to create with last version
                possible_name = sample
                count = 1
                while True:
                    if possible_name in used_names:
                        possible_name = "{}_{}".format(sample, count)
                        count += 1
                    else:
                        break
                rename[sample] = possible_name
                used_names[possible_name] = 1
                continue

        self.utils.filter_fasta_by_ids(temp_out_fasta, rename, out_file)

        self.utils.remove_temp_file(temp_out_fasta)

        # First mark all as not going (so they can be removed from metadata table)
        for sample in possiblename_to_id.keys():
            try:
                dataset_consensus = DatasetConsensus.objects.get(
                    id=possiblename_to_id[sample]
                )
                dataset_consensus.seq_name_all_consensus = (
                    SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_NOBUILD
                )
                dataset_consensus.save()
            except DatasetConsensus.DoesNotExist:  ## need to create with last version
                continue

        #  Now register the ones that are going to nextstrain
        for sample in rename.keys():
            ### set all_consensus name in table dataset_consensus
            if (sample in possiblename_to_id) and (possiblename_to_id[sample] != -1):
                try:
                    dataset_consensus = DatasetConsensus.objects.get(
                        id=possiblename_to_id[sample]
                    )
                    dataset_consensus.seq_name_all_consensus = rename[sample]
                    dataset_consensus.save()
                except (
                    DatasetConsensus.DoesNotExist
                ):  ## need to create with last version
                    continue

        return out_file

    def run_nextstrain(self, dataset, build):
        """
        Runs nextStrain
        File expected: DATASET_FILE_NAME_auspice_zip = "auspice.zip"
        """

        # check what is the build that is configured, otherwise use the default
        # build = SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_parameter

        # See if there is a build parameter specific for this dataset, in which case use it
        # parameters_list = Parameter.objects.filter(dataset=dataset)
        # if len(list(parameters_list)) == 1:
        #    build = list(parameters_list)[0].parameter

        sequences_file = dataset.get_global_file_by_dataset(
            TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_RESULT_all_consensus
        )
        metadata_file = dataset.get_global_file_by_dataset(
            TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_RESULT_NEXTSTRAIN_TSV
        )

        tree_file = None
        alignment_file = None
        auspice_zip = None
        # TODO Make this more generic...
        if build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_ncov:
            tree_file, alignment_file, auspice_zip = self.software.run_nextstrain_ncov(
                alignments=sequences_file, metadata=metadata_file
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_mpx:
            tree_file, alignment_file, auspice_zip = self.software.run_nextstrain_mpx(
                alignments=sequences_file, metadata=metadata_file
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_mpox_clade_i:
            tree_file, alignment_file, auspice_zip = self.software.run_nextstrain_mpox(
                alignments=sequences_file, metadata=metadata_file, type="clade-i"
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_mpox_mpxv:
            tree_file, alignment_file, auspice_zip = self.software.run_nextstrain_mpox(
                alignments=sequences_file, metadata=metadata_file, type="all-clades"
            )            
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_flu_h3n2_12y:
            # This one can have extra parameters such as strain (default: h3n2, h1n1, etc...) and time period (default: 12y)
            tree_file, alignment_file, auspice_zip = self.software.run_nextstrain_flu(
                alignments=sequences_file, metadata=metadata_file, strain="h3n2"
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_flu_h3n2_051225:
            # This one can have extra parameters such as strain (default: h3n2, h1n1, etc...) and segment (default: ha)
            tree_file, alignment_file, auspice_zip = self.software.run_nextstrain_flu_051225(
                alignments=sequences_file, metadata=metadata_file, strain="h3n2"
            )            
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_flu_h1n1pdm_12y:
            tree_file, alignment_file, auspice_zip = self.software.run_nextstrain_flu(
                alignments=sequences_file, metadata=metadata_file, strain="h1n1pdm"
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_flu_h1n1pdm_051225:
            # This one can have extra parameters such as strain (default: h3n2, h1n1, etc...) and segment (default: ha)
            tree_file, alignment_file, auspice_zip = self.software.run_nextstrain_flu_051225(
                alignments=sequences_file, metadata=metadata_file, strain="h1n1pdm"
            )               
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_flu_vic_12y:
            tree_file, alignment_file, auspice_zip = self.software.run_nextstrain_flu(
                alignments=sequences_file, metadata=metadata_file, strain="vic"
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_flu_vic_051225:
            tree_file, alignment_file, auspice_zip = self.software.run_nextstrain_flu_051225(
                alignments=sequences_file, metadata=metadata_file, strain="vic"
            )            
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_flu_yam_12y:
            tree_file, alignment_file, auspice_zip = self.software.run_nextstrain_flu(
                alignments=sequences_file, metadata=metadata_file, strain="yam"
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_avianflu_h5n1_ha:
            tree_file, alignment_file, auspice_zip = (
                self.software.run_nextstrain_avianflu(
                    alignments=sequences_file,
                    metadata=metadata_file,
                    strain="h5n1",
                    gene="ha",
                )
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_avianflu_h5n1_na:
            tree_file, alignment_file, auspice_zip = (
                self.software.run_nextstrain_avianflu(
                    alignments=sequences_file,
                    metadata=metadata_file,
                    strain="h5n1",
                    gene="na",
                )
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_avianflu_h5n1_pb2:
            tree_file, alignment_file, auspice_zip = (
                self.software.run_nextstrain_avianflu(
                    alignments=sequences_file,
                    metadata=metadata_file,
                    strain="h5n1",
                    gene="pb2",
                )
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_avianflu_h5n1_pb1:
            tree_file, alignment_file, auspice_zip = (
                self.software.run_nextstrain_avianflu(
                    alignments=sequences_file,
                    metadata=metadata_file,
                    strain="h5n1",
                    gene="pb1",
                )
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_avianflu_h5n1_pa:
            tree_file, alignment_file, auspice_zip = (
                self.software.run_nextstrain_avianflu(
                    alignments=sequences_file,
                    metadata=metadata_file,
                    strain="h5n1",
                    gene="pa",
                )
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_avianflu_h5n1_np:
            tree_file, alignment_file, auspice_zip = (
                self.software.run_nextstrain_avianflu(
                    alignments=sequences_file,
                    metadata=metadata_file,
                    strain="h5n1",
                    gene="np",
                )
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_avianflu_h5n1_mp:
            tree_file, alignment_file, auspice_zip = (
                self.software.run_nextstrain_avianflu(
                    alignments=sequences_file,
                    metadata=metadata_file,
                    strain="h5n1",
                    gene="mp",
                )
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_avianflu_h5n1_ns:
            tree_file, alignment_file, auspice_zip = (
                self.software.run_nextstrain_avianflu(
                    alignments=sequences_file,
                    metadata=metadata_file,
                    strain="h5n1",
                    gene="ns",
                )
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_rsv_a:
            tree_file, alignment_file, auspice_zip = self.software.run_nextstrain_rsv(
                alignments=sequences_file, metadata=metadata_file, type="a"
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_rsv_b:
            tree_file, alignment_file, auspice_zip = self.software.run_nextstrain_rsv(
                alignments=sequences_file, metadata=metadata_file, type="b"
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_dengue_all:
            tree_file, alignment_file, auspice_zip = (
                self.software.run_nextstrain_dengue(
                    alignments=sequences_file, metadata=metadata_file, type="all"
                )
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_dengue_denv1:
            tree_file, alignment_file, auspice_zip = (
                self.software.run_nextstrain_dengue(
                    alignments=sequences_file, metadata=metadata_file, type="denv1"
                )
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_dengue_denv2:
            tree_file, alignment_file, auspice_zip = (
                self.software.run_nextstrain_dengue(
                    alignments=sequences_file, metadata=metadata_file, type="denv2"
                )
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_dengue_denv3:
            tree_file, alignment_file, auspice_zip = (
                self.software.run_nextstrain_dengue(
                    alignments=sequences_file, metadata=metadata_file, type="denv3"
                )
            )
        elif build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_dengue_denv4:
            tree_file, alignment_file, auspice_zip = (
                self.software.run_nextstrain_dengue(
                    alignments=sequences_file, metadata=metadata_file, type="denv4"
                )
            )
        elif build in [
            SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_generic,
            SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_generic_time,
        ]:
            # Need to get the reference fasta and genbank (if there is more than one reference, get the first one??)
            time = False
            if build == SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_generic_time:
                time = True
            reference = dataset.get_first_reference()
            if reference is None or reference == "":
                raise Exception("No Reference was found. The generic build needs at least one reference")
            try:
                # Check for user?
                tree_file, alignment_file, auspice_zip = (
                    self.software.run_nextstrain_generic(
                        alignments=sequences_file,
                        metadata=metadata_file,
                        ref_fasta=reference.get_reference_fasta(TypePath.MEDIA_ROOT),
                        ref_genbank=reference.get_reference_gbk(TypePath.MEDIA_ROOT),
                        time=time,
                    )
                )
            except Reference.DoesNotExist:
                raise Exception("No Reference was found. The generic build needs at least one reference")
            except Exception as e:
                raise e
        else:
            # It is not supposed to arrive here
            raise Exception("Unknown error. Please contact the administrators")

        # temp_dir = self.software.run_nextstrain(Dataset.REFERENCE_NAME, sequences_file, metadata_file)
        return tree_file, alignment_file, auspice_zip

    def collect_sample_table(self, dataset, column_separator, type_list, build):
        """
        collect sample table
        column_separator : COMMA or TAB
        id,fastq1,fastq2,data set,vaccine status,week,onset date,collection date,lab reception date,latitude,longitude
        :param type_list
        2) sample list, used to upload in the tree
        "onset date" is the date

        """

        data_columns = DataColumns(build)

        ### join all
        dt_out_id_project = {}
        for dataset_consensus in dataset.dataset_consensus.all():
            if dataset_consensus.is_deleted:
                continue
            if dataset_consensus.is_error:
                continue

            ## add metadata to reference
            if not dataset_consensus.reference is None:
                # At the moment there is no upload of metadata so we just go with 'id'
                if dataset_consensus.reference.name not in dt_out_id_project:
                    data_columns.add_header(dataset_consensus.reference.name, ["id"])
                    data_columns.add_metadata(
                        dataset_consensus.reference.name,
                        dataset_consensus.reference.name,
                        dataset_consensus.reference.name,
                        dataset_consensus.seq_name_all_consensus,
                        [dataset_consensus.seq_name_all_consensus],
                    )
                continue

            ## add metadata to consensus
            if not dataset_consensus.consensus is None:
                # At the moment there is no upload of metadata so we just go with 'id'
                if dataset_consensus.consensus.name not in dt_out_id_project:
                    data_columns.add_header(dataset_consensus.consensus.name, ["id"])
                    data_columns.add_metadata(
                        dataset_consensus.consensus.name,
                        dataset_consensus.consensus.name,
                        dataset_consensus.consensus.name,
                        dataset_consensus.seq_name_all_consensus,
                        [dataset_consensus.seq_name_all_consensus],
                    )

                continue

            ## read metadata from file
            if not dataset_consensus.project_sample is None:
                ## test if already processed
                if dataset_consensus.project_sample.project.pk in dt_out_id_project:
                    row = dt_out_id_project[
                        dataset_consensus.project_sample.project.pk
                    ].get(dataset_consensus.project_sample.sample.name, None)
                    if not row is None:
                        data_columns.add_metadata(
                            dataset_consensus.project_sample.project.pk,
                            dataset_consensus.project_sample.pk,
                            dataset_consensus.project_sample.project.name,
                            dataset_consensus.seq_name_all_consensus,
                            row,
                        )
                    continue

                dt_out_temp = {}
                ## start looking for simple list
                file_csv = (
                    dataset_consensus.project_sample.project.get_global_file_by_project(
                        TypePath.MEDIA_ROOT,
                        Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV_simple,
                    )
                )
                if not os.path.exists(file_csv):
                    file_csv = dataset_consensus.project_sample.project.get_global_file_by_project(
                        TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV
                    )
                if os.path.exists(file_csv):

                    with open(file_csv) as handle_in:
                        csv_reader = csv.reader(
                            handle_in, delimiter=Constants.SEPARATOR_COMMA
                        )
                        for line, row in enumerate(csv_reader):
                            if line == 0:
                                data_columns.add_header(
                                    dataset_consensus.project_sample.project.pk, row
                                )
                            else:
                                dt_out_temp[row[0]] = row
                ### set processed
                dt_out_id_project[dataset_consensus.project_sample.project.pk] = (
                    dt_out_temp
                )
                row = dt_out_id_project[
                    dataset_consensus.project_sample.project.pk
                ].get(dataset_consensus.project_sample.sample.name, None)
                if not row is None:
                    data_columns.add_metadata(
                        dataset_consensus.project_sample.project.pk,
                        dataset_consensus.project_sample.pk,
                        dataset_consensus.project_sample.project.name,
                        dataset_consensus.seq_name_all_consensus,
                        row,
                    )

        ## save file
        out_file = self.utils.get_temp_file(
            "dataset_out",
            (
                FileExtensions.FILE_CSV
                if column_separator == Constants.SEPARATOR_COMMA
                else FileExtensions.FILE_TSV
            ),
        )
        with open(out_file, "w", newline="") as handle_out:
            csv_writer = csv.writer(
                handle_out,
                delimiter=column_separator,
                quotechar='"',
                quoting=(
                    csv.QUOTE_MINIMAL
                    if column_separator == Constants.SEPARATOR_COMMA
                    else csv.QUOTE_ALL
                ),
            )

            ### save metadata
            n_count = data_columns.save_rows(csv_writer)

        if n_count == 0:
            os.unlink(out_file)
            return None
        return out_file

    def collect_nextstrain_table(self, dataset, column_separator, build):
        """
        collect sample table
        column_separator : COMMA or TAB
        id,fastq1,fastq2,data set,vaccine status,week,onset date,collection date,lab reception date,latitude,longitude
        :param type_list
        2) sample list, used to upload in the tree

        """
        # May change depending on the build

        data_columns = DataColumns(build)
        ### join all
        dt_out_id_project = {}
        for dataset_consensus in dataset.dataset_consensus.all():
            if dataset_consensus.is_deleted:
                continue
            if dataset_consensus.is_error:
                continue

            ## add metadata to reference
            if not dataset_consensus.reference is None:
                # print("Need to do metadata for reference: {}".format(dataset_consensus.reference))
                if dataset_consensus.reference.name not in dt_out_id_project:
                    if (
                        dataset_consensus.seq_name_all_consensus
                        != SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_NOBUILD
                    ):
                        data_columns.add_header(
                            dataset_consensus.reference.name, ["id"]
                        )
                        consensus_length = self.utils.get_total_length_fasta(
                            dataset_consensus.get_consensus_file(TypePath.MEDIA_ROOT)
                        )
                        data_columns.add_metadata(
                            dataset_consensus.reference.name,
                            dataset_consensus.reference.name,
                            dataset_consensus.reference.name,
                            dataset_consensus.seq_name_all_consensus,
                            [dataset_consensus.seq_name_all_consensus],
                            consensus_length,
                        )
                continue
            ## add metadata to consensus
            if not dataset_consensus.consensus is None:
                # print("Need to do metadata for user-provided consensus: {}".format(dataset_consensus.consensus))
                if dataset_consensus.consensus.name not in dt_out_id_project:
                    if (
                        dataset_consensus.seq_name_all_consensus
                        != SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_NOBUILD
                    ):
                        data_columns.add_header(
                            dataset_consensus.consensus.name, ["id"]
                        )
                        consensus_length = self.utils.get_total_length_fasta(
                            dataset_consensus.get_consensus_file(TypePath.MEDIA_ROOT)
                        )
                        data_columns.add_metadata(
                            dataset_consensus.consensus.name,
                            dataset_consensus.consensus.name,
                            dataset_consensus.consensus.name,
                            dataset_consensus.seq_name_all_consensus,
                            [dataset_consensus.seq_name_all_consensus],
                            consensus_length,
                        )
                continue

            ## read metadata from file
            if not dataset_consensus.project_sample is None:
                # print("Need to do metadata for project-base consensus: {}".format(dataset_consensus.project_sample))
                ## test if already processed
                if dataset_consensus.project_sample.project.pk in dt_out_id_project:
                    row = dt_out_id_project[
                        dataset_consensus.project_sample.project.pk
                    ].get(dataset_consensus.project_sample.sample.name, None)

                    if not row is None:
                        ### get consensus length
                        if (
                            dataset_consensus.seq_name_all_consensus
                            != SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_NOBUILD
                        ):
                            consensus_length = self.utils.get_total_length_fasta(
                                dataset_consensus.get_consensus_file(
                                    TypePath.MEDIA_ROOT
                                )
                            )
                            data_columns.add_metadata(
                                dataset_consensus.project_sample.project.pk,
                                dataset_consensus.project_sample.pk,
                                dataset_consensus.project_sample.project.name,
                                dataset_consensus.seq_name_all_consensus,
                                row,
                                consensus_length,
                            )
                    continue

                dt_out_temp = {}
                ## start looking for simple list
                file_csv = (
                    dataset_consensus.project_sample.project.get_global_file_by_project(
                        TypePath.MEDIA_ROOT,
                        Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV_simple,
                    )
                )
                if not os.path.exists(file_csv):
                    file_csv = dataset_consensus.project_sample.project.get_global_file_by_project(
                        TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV
                    )
                if os.path.exists(file_csv):

                    with open(file_csv) as handle_in:
                        csv_reader = csv.reader(
                            handle_in, delimiter=Constants.SEPARATOR_COMMA
                        )
                        for line, row in enumerate(csv_reader):
                            if line == 0:
                                data_columns.add_header(
                                    dataset_consensus.project_sample.project.pk, row
                                )
                            else:
                                dt_out_temp[row[0]] = row

                ### set processed
                dt_out_id_project[dataset_consensus.project_sample.project.pk] = (
                    dt_out_temp
                )
                row = dt_out_id_project[
                    dataset_consensus.project_sample.project.pk
                ].get(dataset_consensus.project_sample.sample.name, None)
                if not row is None:
                    if (
                        dataset_consensus.seq_name_all_consensus
                        != SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_NOBUILD
                    ):
                        consensus_length = self.utils.get_total_length_fasta(
                            dataset_consensus.get_consensus_file(TypePath.MEDIA_ROOT)
                        )
                        data_columns.add_metadata(
                            dataset_consensus.project_sample.project.pk,
                            dataset_consensus.project_sample.pk,
                            dataset_consensus.project_sample.project.name,
                            dataset_consensus.seq_name_all_consensus,
                            row,
                            consensus_length,
                        )

        ## get reference file (if it exists)
        reference_tsv = os.path.join(
            getattr(settings, "STATIC_ROOT", None),
            SoftwareNames.SOFTWARE_NEXTSTRAIN_BUILDS_BASE,
            build,
            "data",
            "references_metadata.tsv",
        )

        if not os.path.exists(reference_tsv):
            reference_tsv = None

        ## read last metadata nextstrain file, can exist from external upload
        upload_metadata_file = (
            UploadFiles.objects.filter(
                owner__id=dataset.owner.id,
                is_deleted=False,
                type_file__name=TypeFile.TYPE_FILE_dataset_file_metadata,
                is_valid=True,
                dataset=dataset,
            )
            .order_by("-creation_date")
            .first()
        )
        parse_in_files = ParseNextStrainFiles()
        if not upload_metadata_file is None:
            b_test_char_encoding = True
            parse_in_files.parse_nextstrain_files(
                upload_metadata_file.get_path_to_file(TypePath.MEDIA_ROOT),
                None,
                b_test_char_encoding,
                ParseNextStrainFiles.STATE_READ_metadata_dont_detect_errors_and_chech_nexttrain,
            )

        ## save file
        out_file = self.utils.get_temp_file(
            "dataset_out",
            (
                FileExtensions.FILE_CSV
                if column_separator == Constants.SEPARATOR_COMMA
                else FileExtensions.FILE_TSV
            ),
        )

        with open(out_file, "w", newline="") as handle_out:
            csv_writer = csv.writer(
                handle_out,
                delimiter=column_separator,
                quotechar='"',
                quoting=(
                    csv.QUOTE_MINIMAL
                    if column_separator == Constants.SEPARATOR_COMMA
                    else csv.QUOTE_ALL
                ),
            )

            ### save metadata
            n_count = data_columns.save_rows_nextstrain(
                csv_writer, reference_tsv, parse_in_files
            )

        if n_count == 0:
            os.unlink(out_file)
            return None

        ### update the metadata file status, only process the last one
        if not upload_metadata_file is None:
            upload_metadata_file.is_processed = True
            upload_metadata_file.save()
        return out_file

    def zip_several_files(self, dataset):

        temp_dir = self.utils.get_temp_dir()

        ## sample file result
        for type_file in Dataset.VECT_files_to_zip:
            if os.path.exists(
                dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, type_file)
            ):
                self.utils.link_file(
                    dataset.get_global_file_by_dataset(TypePath.MEDIA_ROOT, type_file),
                    os.path.join(temp_dir, type_file),
                )

        ## all files zipped
        zip_out = self.software.zip_files_in_path(temp_dir)
        if os.path.exists(zip_out):
            self.utils.move_file(
                zip_out,
                dataset.get_global_file_by_dataset(
                    TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_all_files_zipped
                ),
            )
        else:
            self.utils.remove_file(
                dataset.get_global_file_by_dataset(
                    TypePath.MEDIA_ROOT, Dataset.DATASET_FILE_NAME_all_files_zipped
                )
            )


