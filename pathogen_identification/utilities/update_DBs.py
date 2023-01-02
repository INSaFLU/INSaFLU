import os
import sys
from ast import Param
from typing import Type

from django.contrib.auth.models import User
from django.core.files import File
from django.db import IntegrityError, transaction
from pathogen_identification.models import (
    QC_REPORT,
    ContigClassification,
    FinalReport,
    ParameterSet,
    PIProject_Sample,
    Projects,
    RawReference,
    ReadClassification,
    ReferenceContigs,
    ReferenceMap_Main,
    RunAssembly,
    RunDetail,
    RunIndex,
    RunMain,
    RunRemapMain,
    SampleQC,
)
from pathogen_identification.modules.object_classes import Sample_runClass
from pathogen_identification.modules.remap_class import Mapping_Instance
from pathogen_identification.modules.run_main import RunMain_class


####################################################################################################################
####################################################################################################################
def Update_project(project_directory_path, user: str = "admin"):
    """Updates the project"""
    project_directory_path = os.path.dirname(project_directory_path)
    project_name = os.path.basename(project_directory_path)

    try:
        user = User.objects.get(username=user)
    except User.DoesNotExist:
        print("User does not exist")
        sys.exit(1)

    try:
        project = Projects.objects.get(
            name=project_name, created_by=user, is_deleted=False
        )

    except Projects.DoesNotExist:
        print("project_name: ", project_name)
        project = Projects(
            name=project_name,
            full_path=project_directory_path,
            project_type=Projects.INHOUSE,
            created_by=user,
        )
        project.save()


def Update_Sample(sample_class: Sample_runClass):
    """
    Update Sample class.

    :param sample_class:
    :return: None
    """

    user = User.objects.get(username=sample_class.user_name)

    project = Projects.objects.get(
        name=sample_class.project_name, owner=user, is_deleted=False
    )

    try:
        PIProject_Sample.objects.get(
            name=sample_class.sample_name,
            project=project,
        )
    except PIProject_Sample.DoesNotExist:
        Update_sample(sample_class)

    sample = PIProject_Sample.objects.get(
        name=sample_class.sample_name,
        project=project,
    )

    try:
        SampleQC.objects.get(sample=sample)
    except SampleQC.DoesNotExist:
        Update_sample_qc(sample_class)


def Update_sample(sample_class: Sample_runClass):
    """update sample_class.
    :param sample_class:
    :return: None
    """
    user = User.objects.get(username=sample_class.user_name)
    project = Projects.objects.get(
        name=sample_class.project_name, owner=user, is_deleted=False
    )

    try:
        sample = PIProject_Sample.objects.get(
            name=sample_class.sample_name,
            project=project,
        )

    except PIProject_Sample.DoesNotExist:
        #
        project = Projects.objects.get(
            name=sample_class.project_name,
            owner=user,
            is_deleted=False,
        )

        sample = PIProject_Sample(
            project=project,
            name_extended=sample_class.sample_name,
            name=os.path.splitext(sample_class.sample_name)[0],
            technology=sample_class.technology,
            type=sample_class.type,  # SE or PE
            combinations=sample_class.combinations,
            input=sample_class.input,  # input files,
            report="report",
        )
        sample.save()


def Update_sample_qc(sample_class: Sample_runClass):
    """update sample_class.qc_data.
    :param sample_class:
    :return: None
    """

    user = User.objects.get(username=sample_class.user_name)
    project = Projects.objects.get(
        name=sample_class.project_name, owner=user, is_deleted=False
    )

    sample = PIProject_Sample.objects.get(
        project=project, name=sample_class.sample_name
    )

    percent_passed = (
        int(sample_class.reads_after_processing)
        / int(sample_class.reads_before_processing)
    ) * 100

    input_report = open(sample_class.input_fastqc_report, "r")
    processed_report = open(sample_class.processed_fastqc_report, "r")

    try:
        sampleqc = SampleQC(
            sample=sample,
            software=sample_class.qc_soft,
            encoding=sample_class.technology,
            input_reads=f"{int(sample_class.reads_before_processing):,}",
            processed_reads=f"{int(sample_class.reads_after_processing):,}",
            percent_passed=round(percent_passed, 2),
            sequence_length=sample_class.qcdata["processed"].loc["Sequence length"][
                "value"
            ],
            percent_gc=sample_class.qcdata["processed"].loc["%GC"]["value"],
            input_fastqc_report=File(
                input_report, name=os.path.basename(sample_class.input_fastqc_report)
            ),
            processed_fastqc_report=File(
                processed_report,
                name=os.path.basename(sample_class.processed_fastqc_report),
            ),
        )

        sampleqc.save()

    except:
        print(f"failed to input sample {sample_class.sample_name}")
    finally:
        input_report.close()
        processed_report.close()


def Update_QC_report(sample_class: Sample_runClass, parameter_set: ParameterSet):
    """
    Update QC data for sample_class.

    :param sample_class:
    :return: None
    """
    user = User.objects.get(username=sample_class.user_name)
    project = Projects.objects.get(
        name=sample_class.project_name, owner=user, is_deleted=False
    )

    sample = PIProject_Sample.objects.get(
        project=project, name=sample_class.sample_name
    )

    try:
        qc_report = QC_REPORT.objects.get(sample=sample, report_source=QC_REPORT.RAW)
    except QC_REPORT.DoesNotExist:
        qc_report = QC_REPORT(
            sample=sample,
            report_source=QC_REPORT.RAW,
            QC_report=sample_class.input_fastqc_report,
        )
        qc_report.save()

    try:
        qc_report = QC_REPORT.objects.get(
            sample=sample, report_source=QC_REPORT.PROCESSED
        )
    except QC_REPORT.DoesNotExist:
        qc_report = QC_REPORT(
            sample=sample,
            report_source=QC_REPORT.PROCESSED,
            QC_report=sample_class.processed_fastqc_report,
        )
        qc_report.save()


@transaction.atomic
def Update_Sample_Runs(run_class: RunMain_class, parameter_set: ParameterSet):
    """get run data
    Update ALL run TABLES:
    - RunMain,
    - RunDetail,
    - RunAssembly,
    - ReadClassification,
    - ContigClassification,
    - RunRemapMain,
    - ReferenceMap_Main
    - ReferenceContigs
    - FinalReport,

    :param sample_class:
    :return: run_data
    """

    try:
        with transaction.atomic():
            Update_RunMain(run_class, parameter_set)
            Update_Run_Detail(run_class, parameter_set)
            Update_Sample_Runs_DB(run_class, parameter_set)
            Update_RefMap_DB(run_class, parameter_set)

        return True

    except IntegrityError as e:
        print(f"failed to update sample {run_class.sample_name} {e}")
        return False


@transaction.atomic
def Update_RunMain_Initial(run_class: RunMain_class, parameter_set: ParameterSet):
    """get run data
    Update ALL run TABLES:
    - RunMain,
    :param sample_class:
    :return: run_data
    """

    try:
        with transaction.atomic():
            Update_RunMain(run_class, parameter_set)

        return True

    except IntegrityError as e:
        print(f"failed to update sample {run_class.sample_name}")
        return False


@transaction.atomic
def Update_RunMain_Secondary(run_class: RunMain_class, parameter_set: ParameterSet):
    """get run data
    Update ALL run TABLES:
    - RunMain,

    :param sample_class:
    :return: run_data
    """

    try:
        with transaction.atomic():
            Update_RunMain_noCheck(run_class, parameter_set)
            Update_Run_Detail_noCheck(run_class, parameter_set)
        return True

    except IntegrityError as e:
        print(f"failed to update sample {run_class.sample_name}")
        return False


@transaction.atomic
def Update_Assembly(run_class: RunMain_class, parameter_set: ParameterSet):
    """get run data
    Update TABLES:
    - RunMain,
    - RunAssembly,

    :param sample_class:
    :return: run_data
    """

    try:
        with transaction.atomic():
            Update_RunMain_noCheck(run_class, parameter_set)
            Update_Run_Assembly(run_class, parameter_set)

        return True

    except IntegrityError as e:
        print(f"failed to update sample {run_class.sample_name}")
        return False


@transaction.atomic
def Update_Classification(run_class: RunMain_class, parameter_set: ParameterSet):
    """get run data
    Update TABLES:
    - RunMain,
    - ReadClassification,
    - ContigClassification,

    :param sample_class:
    :return: run_data
    """

    try:
        with transaction.atomic():
            Update_RunMain_noCheck(run_class, parameter_set)
            Update_Run_Classification(run_class, parameter_set)

        return True

    except IntegrityError as e:
        print(f"failed to update sample {run_class.sample_name}")
        return False


@transaction.atomic
def Update_Remap(run_class: RunMain_class, parameter_set: ParameterSet):
    """get run data
    Update TABLES:
    - RunMain,
    - RunRemapMain,

    :param sample_class:
    :return: run_data
    """

    sample, runmain, _ = get_run_parents(run_class, parameter_set)
    try:
        with transaction.atomic():
            Update_FinalReport(run_class, runmain, sample)
            Update_RefMap_DB(run_class, parameter_set)
            Update_Run_Detail_noCheck(run_class, parameter_set)
            Update_RunMain_noCheck(run_class, parameter_set, tag="finished")
        return True

    except IntegrityError as e:
        print(f"failed to update sample {run_class.sample_name}")
        return False


def retrieve_number_of_runs(project_name, sample_name, username):
    """
    retrieve number of runs for a given project.

    :param run_class:
    :return: number of runs
    """

    user = User.objects.get(username=username)

    try:
        project = Projects.objects.get(name=project_name, owner=user, is_deleted=False)
    except Projects.DoesNotExist:
        print(f"project {project_name} does not exist")
        return 0

    try:
        sample = PIProject_Sample.objects.get(name=sample_name, project=project)
    except PIProject_Sample.DoesNotExist:
        return 0

    return RunMain.objects.filter(project=project, sample=sample).count() + 1


def RunIndex_Update_Retrieve_Key(project_name, sample_name):

    run_index = retrieve_number_of_runs(project_name, sample_name)

    new_name = f"run_{run_index}"

    project = Projects.objects.get(name=project_name, is_deleted=False)
    sample = PIProject_Sample.objects.get(
        name=sample_name,
        project__name=project_name,
    )

    try:
        run_index = RunIndex.objects.get(
            project=project, sample=sample, name=new_name
        )  # check if run_index already exists
    except RunIndex.DoesNotExist:
        run_index = RunIndex(project=project, sample=sample, name=new_name)
        run_index.save()

    return new_name


def Update_RunMain(run_class: RunMain_class, parameter_set: ParameterSet):
    """update run data for run_class. Update run_class.run_data.

    :param run_class:
    :return: None
    """
    user = User.objects.get(username=run_class.username)
    project = Projects.objects.get(
        name=run_class.project_name, owner=user, is_deleted=False
    )

    sample = PIProject_Sample.objects.get(
        name=run_class.sample.sample_name,
        project=project,
    )

    reads_after_processing = run_class.sample.reads_after_processing

    if run_class.sample.reads_before_processing > 0:
        reads_proc_percent = (
            reads_after_processing / run_class.sample.reads_before_processing
        ) * 100

    else:
        reads_proc_percent = 0

    enrichment_method = run_class.enrichment_drone.classifier_method.name
    enrichment = run_class.enrichment_drone.deployed

    host_depletion_method = run_class.depletion_drone.classifier_method.name
    host_depletion = run_class.depletion_drone.deployed
    try:
        runmain = RunMain.objects.get(
            project__name=run_class.sample.project_name,
            suprun=run_class.suprun,
            sample=sample,
            name=run_class.prefix,
            parameter_set=parameter_set,
        )
    except RunMain.DoesNotExist:

        runmain = RunMain(
            parameter_set=parameter_set,
            suprun=run_class.suprun,
            project=project,
            sample=sample,
            name=run_class.prefix,
            params_file_path=run_class.params_file_path,
            processed_reads_r1=run_class.sample.r1.current,
            processed_reads_r2=run_class.sample.r2.current,
            assembly_performed=run_class.assembly_drone.assembly_exists,
            assembly_method=run_class.assembly_drone.assembly_method.name,
            reads_after_processing=f"{reads_after_processing:,}",
            reads_proc_percent=round(reads_proc_percent, 2),
            host_depletion=host_depletion_method,
            host_depletion_performed=host_depletion,
            host_depletion_args=run_class.depletion_drone.classifier_method.args,
            host_depletion_db=run_class.depletion_drone.classifier_method.db_name,
            enrichment_performed=enrichment,
            enrichment=enrichment_method,
            enrichment_args=run_class.enrichment_drone.classifier_method.args,
            enrichment_db=run_class.enrichment_drone.classifier_method.db_name,
            assembly_max=f"{run_class.assembly_drone.contig_summary['contig_length'].max():,}",
            remap=run_class.remapping_method.name,
            remap_args=run_class.remapping_method.args,
            read_classification=run_class.read_classification_drone.classifier_method.name,
            contig_classification=run_class.contig_classification_drone.classifier_method.name,
            runtime=f"{run_class.exec_time / 60:.2f} m",
            report="initial",
            # static_dir=run_class.static_dir,
        )

        runmain.save()


def Sample_update_combinations(run_class: Type[RunMain_class]):

    user = User.objects.get(username=run_class.username)
    project = Projects.objects.get(
        name=run_class.project_name, owner=user, is_deleted=False
    )

    sample = PIProject_Sample.objects.get(
        project=project,
        name=run_class.sample.sample_name,
    )

    sample.combinations = sample.combinations + 1

    sample.save()


def get_run_parents(run_class: RunMain_class, parameter_set: ParameterSet):
    """get run parents for run_class. Update run_class.run_data."""
    user = User.objects.get(username=run_class.username)
    project = Projects.objects.get(
        name=run_class.project_name, owner=user, is_deleted=False
    )

    sample = PIProject_Sample.objects.get(
        project=project,
        name=run_class.sample.sample_name,
    )

    try:
        runmain = RunMain.objects.get(
            project=project,
            suprun=run_class.suprun,
            sample=sample,
            name=run_class.prefix,
            parameter_set=parameter_set,
        )

    except RunMain.DoesNotExist:
        return None, None, None

    return sample, runmain, project


def Update_RunMain_noCheck(
    run_class: RunMain_class, parameter_set: ParameterSet, tag="secondary"
):
    """update run data for run_class. Update run_class.run_data.

    :param run_class:
    :return: None
    """
    sample, runmain, project = get_run_parents(run_class, parameter_set)

    if sample is None or runmain is None:
        return

    reads_after_processing = run_class.sample.reads_after_processing

    if run_class.sample.reads_before_processing > 0:
        reads_proc_percent = (
            reads_after_processing / run_class.sample.reads_before_processing
        ) * 100

    else:
        reads_proc_percent = 0

    enrichment_method = run_class.enrichment_drone.classifier_method.name
    enrichment = run_class.enrichment_drone.deployed

    host_depletion_method = run_class.depletion_drone.classifier_method.name
    host_depletion = run_class.depletion_drone.deployed

    runmain.parameter_set = parameter_set
    runmain.suprun = run_class.suprun
    runmain.project = project
    runmain.sample = sample
    runmain.name = run_class.prefix
    runmain.params_file_path = run_class.params_file_path
    runmain.processed_reads_r1 = run_class.sample.r1.current
    runmain.processed_reads_r2 = run_class.sample.r2.current
    runmain.assembly_performed = run_class.assembly_drone.assembly_exists
    runmain.assembly_method = run_class.assembly_drone.assembly_method.name
    runmain.reads_after_processing = f"{reads_after_processing:,}"
    runmain.reads_proc_percent = round(reads_proc_percent, 2)
    runmain.host_depletion = host_depletion_method
    runmain.host_depletion_performed = host_depletion
    runmain.host_depletion_args = run_class.depletion_drone.classifier_method.args
    runmain.host_depletion_db = run_class.depletion_drone.classifier_method.db_name
    runmain.enrichment_performed = enrichment
    runmain.enrichment = enrichment_method
    runmain.enrichment_args = run_class.enrichment_drone.classifier_method.args
    runmain.enrichment_db = run_class.enrichment_drone.classifier_method.db_name
    runmain.assembly_max = (
        f"{run_class.assembly_drone.contig_summary['contig_length'].max():,}"
    )
    runmain.remap = run_class.remapping_method.name
    runmain.remap_args = run_class.remapping_method.args
    runmain.read_classification = (
        run_class.read_classification_drone.classifier_method.name
    )
    runmain.contig_classification = (
        run_class.contig_classification_drone.classifier_method.name
    )
    runmain.runtime = f"{run_class.exec_time / 60:.2f} m"
    runmain.report = tag
    # static_dir=run_class.static_dir,

    runmain.save()


def Update_Run_Detail(run_class: RunMain_class, parameter_set: ParameterSet):
    """
    Update ALL run TABLES for one run_class.:
    - RunMain,
    - RunDetail,
    - RunAssembly,
    - ReadClassification,
    - ContigClassification,
    - RunRemapMain,
    - ReferenceMap_Main
    - ReferenceContigs
    - FinalReport,

    :param run_class:
    :return: run_data
    """
    # Sample_update_combinations(run_class)

    sample, runmain, _ = get_run_parents(run_class, parameter_set)

    if sample is None or runmain is None:
        return

    try:
        run_detail = RunDetail.objects.get(run=runmain)
    except RunDetail.DoesNotExist:
        run_detail = RunDetail(
            run=runmain,
            sample=sample,
            max_depth=run_class.run_detail_report.max_depth,  #
            max_depthR=run_class.run_detail_report.max_depthR,  #
            max_gaps=run_class.run_detail_report.max_gaps,  #
            max_prop=run_class.run_detail_report.max_prop,  #
            max_mapped=run_class.run_detail_report.max_mapped,  #
            input=run_class.run_detail_report.input,  #
            enriched_reads=run_class.run_detail_report.enriched_reads,  #
            enriched_reads_percent=run_class.run_detail_report.enriched_reads_percent,  #
            depleted_reads=run_class.run_detail_report.depleted_reads,  #
            depleted_reads_percent=run_class.run_detail_report.depleted_reads_percent,  #
            processed=run_class.run_detail_report.processed,  #
            processed_percent=run_class.run_detail_report.processed_percent,  #
            sift_preproc=run_class.run_detail_report.sift_preproc,  #
            sift_remap=run_class.run_detail_report.sift_remap,  #
            sift_removed_pprc=run_class.run_detail_report.sift_removed_pprc,
            processing_final=run_class.run_detail_report.processing_final,  #
            processing_final_percent=run_class.run_detail_report.processing_final_percent,  #
            merged=run_class.run_detail_report.merged,  #
            merged_number=run_class.run_detail_report.merged_number,  #
            merged_files=run_class.run_detail_report.merged_files,  #
        )

        run_detail.save()


def Update_Run_Detail_noCheck(run_class: RunMain_class, parameter_set: ParameterSet):
    """
    Update ALL run TABLES for one run_class.:
    - RunMain,
    - RunDetail,
    - RunAssembly,
    - ReadClassification,
    - ContigClassification,
    - RunRemapMain,
    - ReferenceMap_Main
    - ReferenceContigs
    - FinalReport,

    :param run_class:
    :return: run_data
    """
    # Sample_update_combinations(run_class)

    sample, runmain, _ = get_run_parents(run_class, parameter_set)

    if sample is None or runmain is None:
        return

    run_detail_exists = RunDetail.objects.filter(run=runmain, sample=sample).exists()

    if run_detail_exists:
        run_detail = RunDetail.objects.get(run=runmain, sample=sample)

        run_detail.run = runmain
        run_detail.sample = sample
        run_detail.max_depth = run_class.run_detail_report.max_depth  #
        run_detail.max_depthR = run_class.run_detail_report.max_depthR  #
        run_detail.max_gaps = run_class.run_detail_report.max_gaps  #
        run_detail.max_prop = run_class.run_detail_report.max_prop  #
        run_detail.max_mapped = run_class.run_detail_report.max_mapped  #
        run_detail.input = run_class.run_detail_report.input  #
        run_detail.enriched_reads = run_class.run_detail_report.enriched_reads  #
        run_detail.enriched_reads_percent = (
            run_class.run_detail_report.enriched_reads_percent
        )  #
        run_detail.depleted_reads = run_class.run_detail_report.depleted_reads  #
        run_detail.depleted_reads_percent = (
            run_class.run_detail_report.depleted_reads_percent
        )  #
        run_detail.processed = run_class.run_detail_report.processed  #
        run_detail.processed_percent = run_class.run_detail_report.processed_percent  #
        run_detail.sift_preproc = run_class.run_detail_report.sift_preproc  #
        run_detail.sift_remap = run_class.run_detail_report.sift_remap  #
        run_detail.sift_removed_pprc = run_class.run_detail_report.sift_removed_pprc
        run_detail.processing_final = run_class.run_detail_report.processing_final  #
        run_detail.processing_final_percent = (
            run_class.run_detail_report.processing_final_percent
        )  #
        run_detail.merged = run_class.run_detail_report.merged  #
        run_detail.merged_number = run_class.run_detail_report.merged_number  #
        run_detail.merged_files = run_class.run_detail_report.merged_files  #

        run_detail.save()

    else:
        run_detail = RunDetail(
            run=runmain,
            sample=sample,
            max_depth=run_class.run_detail_report.max_depth,  #
            max_depthR=run_class.run_detail_report.max_depthR,  #
            max_gaps=run_class.run_detail_report.max_gaps,  #
            max_prop=run_class.run_detail_report.max_prop,  #
            max_mapped=run_class.run_detail_report.max_mapped,  #
            input=run_class.run_detail_report.input,  #
            enriched_reads=run_class.run_detail_report.enriched_reads,  #
            enriched_reads_percent=run_class.run_detail_report.enriched_reads_percent,  #
            depleted_reads=run_class.run_detail_report.depleted_reads,  #
            depleted_reads_percent=run_class.run_detail_report.depleted_reads_percent,  #
            processed=run_class.run_detail_report.processed,  #
            processed_percent=run_class.run_detail_report.processed_percent,  #
            sift_preproc=run_class.run_detail_report.sift_preproc,  #
            sift_remap=run_class.run_detail_report.sift_remap,  #
            sift_removed_pprc=run_class.run_detail_report.sift_removed_pprc,
            processing_final=run_class.run_detail_report.processing_final,  #
            processing_final_percent=run_class.run_detail_report.processing_final_percent,  #
            merged=run_class.run_detail_report.merged,  #
            merged_number=run_class.run_detail_report.merged_number,  #
            merged_files=run_class.run_detail_report.merged_files,  #
        )

        run_detail.save()


def Update_Run_Assembly(run_class: RunMain_class, parameter_set: ParameterSet):
    """
    Update ALL run TABLES for one run_class.:
    - RunMain,
    - RunDetail,
    - RunAssembly,
    - ReadClassification,
    - ContigClassification,
    - RunRemapMain,
    - ReferenceMap_Main
    - ReferenceContigs
    - FinalReport,

    :param run_class:
    :return: run_data
    """
    # Sample_update_combinations(run_class)

    sample, runmain, _ = get_run_parents(run_class, parameter_set)

    if sample is None or runmain is None:
        return

    try:
        run_assembly = RunAssembly.objects.get(run=runmain, sample=sample)
    except RunAssembly.DoesNotExist:

        run_assembly = RunAssembly(
            run=runmain,
            sample=sample,
            performed=run_class.assembly_report.performed,
            method=run_class.assembly_report.assembly_soft,
            args=run_class.assembly_report.assembly_args,  #
            contig_number=run_class.assembly_report.assembly_number,
            contig_max=run_class.assembly_report.assembly_max,
            contig_min=run_class.assembly_report.assembly_min,
            contig_mean=run_class.assembly_report.assembly_mean,
            contig_trim=run_class.assembly_report.assembly_trim,
            assembly_contigs=run_class.assembly_drone.assembly_file_fasta_gz,
        )
        run_assembly.save()


def Update_Run_Classification(run_class: RunMain_class, parameter_set: ParameterSet):
    """
    Update ALL run TABLES for one run_class.:
    - RunMain,
    - RunDetail,
    - RunAssembly,
    - ReadClassification,
    - ContigClassification,
    - RunRemapMain,
    - ReferenceMap_Main
    - ReferenceContigs
    - FinalReport,

    :param run_class:
    :return: run_data
    """
    # Sample_update_combinations(run_class)

    sample, runmain, _ = get_run_parents(run_class, parameter_set)

    try:
        read_classification = ReadClassification.objects.get(run=runmain, sample=sample)
    except ReadClassification.DoesNotExist:

        read_classification = ReadClassification(
            run=runmain,
            sample=sample,
            read_classification_report=run_class.read_classification_summary,
            performed=run_class.read_classification_results.performed,
            method=run_class.read_classification_results.method,
            args=run_class.read_classification_results.args,
            db=run_class.read_classification_results.db,
            classification_number=run_class.read_classification_results.classification_number,
            classification_minhit=run_class.read_classification_results.classification_minhit,
            success=run_class.read_classification_results.success,
        )
        read_classification.save()

    try:
        contig_classification = ContigClassification.objects.get(
            run=runmain, sample=sample
        )

    except ContigClassification.DoesNotExist:

        contig_classification = ContigClassification(
            run=runmain,
            sample=sample,
            contig_classification_report=run_class.assembly_classification_summary,
            performed=run_class.contig_classification_results.performed,
            method=run_class.contig_classification_results.method,
            args=run_class.contig_classification_results.args,
            db=run_class.contig_classification_results.db,
            classification_number=run_class.contig_classification_results.classification_number,
            classification_minhit=run_class.contig_classification_results.classification_minhit,
            # success=run_class.contig_classification_results.success,
        )
        contig_classification.save()

    try:
        remap_main = RunRemapMain.objects.get(run=runmain, sample=sample)

    except RunRemapMain.DoesNotExist:
        remap_main = RunRemapMain(
            run=runmain,
            sample=sample,
            merged_log=run_class.merged_classification_summary,
            remap_plan=run_class.remap_plan_path,
            performed=run_class.remap_main.performed,
            method=run_class.remap_main.method,
            found_total=run_class.remap_main.found_total,
            coverage_maximum=run_class.remap_main.coverage_max,
            coverage_minimum=run_class.remap_main.coverage_min,
            success=run_class.remap_main.success,
        )
        remap_main.save()

    for ref, row in run_class.raw_targets.iterrows():

        if row.status:
            status = RawReference.STATUS_MAPPED
        else:
            status = RawReference.STATUS_UNMAPPED
        try:
            remap_target = RawReference.objects.get(
                run=runmain,
                taxid=row.taxid,
                accid=row.accid,
            )
        except RawReference.DoesNotExist:
            remap_target = RawReference(
                run=runmain,
                taxid=row.taxid,
                accid=row.accid,
                status=status,
                description=row.description,
                counts=row.counts,
                classification_source=row.source,
            )

            remap_target.save()


def Update_Sample_Runs_DB(run_class: RunMain_class, parameter_set: ParameterSet):
    """
    Update ALL run TABLES for one run_class.:
    - RunMain,
    - RunDetail,
    - RunAssembly,
    - ReadClassification,
    - ContigClassification,
    - RunRemapMain,
    - ReferenceMap_Main
    - ReferenceContigs
    - FinalReport,

    :param run_class:
    :return: run_data
    """
    # Sample_update_combinations(run_class)

    user = User.objects.get(username=run_class.username)
    project = Projects.objects.get(
        name=run_class.project_name, owner=user, is_deleted=False
    )

    sample = PIProject_Sample.objects.get(
        project=project,
        name=run_class.sample.sample_name,
    )

    try:
        runmain = RunMain.objects.get(
            project=project,
            suprun=run_class.suprun,
            sample=sample,
            name=run_class.prefix,
            parameter_set=parameter_set,
        )

    except RunMain.DoesNotExist:
        return

    try:
        run_detail = RunDetail.objects.get(run=runmain, sample=sample)
    except RunDetail.DoesNotExist:

        run_detail = RunDetail(
            run=runmain,
            sample=sample,
            max_depth=run_class.run_detail_report.max_depth,  #
            max_depthR=run_class.run_detail_report.max_depthR,  #
            max_gaps=run_class.run_detail_report.max_gaps,  #
            max_prop=run_class.run_detail_report.max_prop,  #
            max_mapped=run_class.run_detail_report.max_mapped,  #
            input=run_class.run_detail_report.input,  #
            enriched_reads=run_class.run_detail_report.enriched_reads,  #
            enriched_reads_percent=run_class.run_detail_report.enriched_reads_percent,  #
            depleted_reads=run_class.run_detail_report.depleted_reads,  #
            depleted_reads_percent=run_class.run_detail_report.depleted_reads_percent,  #
            processed=run_class.run_detail_report.processed,  #
            processed_percent=run_class.run_detail_report.processed_percent,  #
            sift_preproc=run_class.run_detail_report.sift_preproc,  #
            sift_remap=run_class.run_detail_report.sift_remap,  #
            sift_removed_pprc=run_class.run_detail_report.sift_removed_pprc,
            processing_final=run_class.run_detail_report.processing_final,  #
            processing_final_percent=run_class.run_detail_report.processing_final_percent,  #
            merged=run_class.run_detail_report.merged,  #
            merged_number=run_class.run_detail_report.merged_number,  #
            merged_files=run_class.run_detail_report.merged_files,  #
        )

        run_detail.save()

    try:
        run_assembly = RunAssembly.objects.get(run=runmain, sample=sample)
    except RunAssembly.DoesNotExist:

        run_assembly = RunAssembly(
            run=runmain,
            sample=sample,
            performed=run_class.assembly_report.performed,
            method=run_class.assembly_report.assembly_soft,
            args=run_class.assembly_report.assembly_args,  #
            contig_number=run_class.assembly_report.assembly_number,
            contig_max=run_class.assembly_report.assembly_max,
            contig_min=run_class.assembly_report.assembly_min,
            contig_mean=run_class.assembly_report.assembly_mean,
            contig_trim=run_class.assembly_report.assembly_trim,
            assembly_contigs=run_class.assembly_drone.assembly_file_fasta_gz,
        )
        run_assembly.save()

    try:
        read_classification = ReadClassification.objects.get(run=runmain, sample=sample)
    except ReadClassification.DoesNotExist:

        read_classification = ReadClassification(
            run=runmain,
            sample=sample,
            read_classification_report=run_class.read_classification_summary,
            performed=run_class.read_classification_results.performed,
            method=run_class.read_classification_results.method,
            args=run_class.read_classification_results.args,
            db=run_class.read_classification_results.db,
            classification_number=run_class.read_classification_results.classification_number,
            classification_minhit=run_class.read_classification_results.classification_minhit,
            success=run_class.read_classification_results.success,
        )
        read_classification.save()

    try:
        contig_classification = ContigClassification.objects.get(
            run=runmain, sample=sample
        )

    except ContigClassification.DoesNotExist:

        contig_classification = ContigClassification(
            run=runmain,
            sample=sample,
            contig_classification_report=run_class.assembly_classification_summary,
            performed=run_class.contig_classification_results.performed,
            method=run_class.contig_classification_results.method,
            args=run_class.contig_classification_results.args,
            db=run_class.contig_classification_results.db,
            classification_number=run_class.contig_classification_results.classification_number,
            classification_minhit=run_class.contig_classification_results.classification_minhit,
            # success=run_class.contig_classification_results.success,
        )
        contig_classification.save()

    try:
        remap_main = RunRemapMain.objects.get(run=runmain, sample=sample)

    except RunRemapMain.DoesNotExist:
        remap_main = RunRemapMain(
            run=runmain,
            sample=sample,
            merged_log=run_class.merged_classification_summary,
            remap_plan=run_class.remap_plan_path,
            performed=run_class.remap_main.performed,
            method=run_class.remap_main.method,
            found_total=run_class.remap_main.found_total,
            coverage_maximum=run_class.remap_main.coverage_max,
            coverage_minimum=run_class.remap_main.coverage_min,
            success=run_class.remap_main.success,
        )
        remap_main.save()

    for ref, row in run_class.raw_targets.iterrows():

        if row.status:
            status = RawReference.STATUS_MAPPED
        else:
            status = RawReference.STATUS_UNMAPPED
        try:
            remap_target = RawReference.objects.get(
                run=runmain,
                taxid=row.taxid,
                accid=row.accid,
            )
        except RawReference.DoesNotExist:
            remap_target = RawReference(
                run=runmain,
                taxid=row.taxid,
                accid=row.accid,
                status=status,
                description=row.description,
                counts=row.counts,
                classification_source=row.source,
            )

            remap_target.save()

    Update_FinalReport(run_class, runmain, sample)


def Update_FinalReport(run_class, runmain, sample):

    for i, row in run_class.report.iterrows():
        if row["ID"] == "None":
            continue

        try:
            report_row = FinalReport.objects.get(
                run=runmain,
                sample=sample,
                unique_id=row["unique_id"],
            )
        except FinalReport.DoesNotExist:

            report_row = FinalReport(
                run=runmain,
                sample=sample,
                unique_id=row["unique_id"],
                reference_length=row["contig_length"],
                taxid=row["taxid"],
                accid=row["ID"],
                reference_contig_str=row["contig_string"],
                simple_id=row["simple_id"],
                description=row["description"],
                ref_db=row["refdb"],
                coverage=row["coverage"],
                windows_covered=row["windows_covered"],
                depth=row["Hdepth"],
                depthR=row["HdepthR"],
                ngaps=row["ngaps"],
                mapped_reads=row["mapped"],
                ref_proportion=row["ref_prop"],
                mapped_proportion=row["mapped_prop"],
                mapping_success=row["mapping_success"],
                classification_success=row["classification_success"],
                refa_dotplot_exists=row["refa_dotplot_exists"],
                covplot_exists=row["covplot_exists"],
                refa_dotplot=row["refa_dotplot_path"],
                covplot=row["covplot_path"],
                bam_path=row["bam_path"],
                bai_path=row["bam_index_path"],
                reference_path=row["reference_path"],
                reference_index_path=row["reference_index_path"],
                reference_assembly_paf=row["reference_assembly_paf"],
                mapped_scaffolds_path=row["mapped_scaffolds_path"],
                mapped_scaffolds_index_path=row["mapped_scaffolds_index_path"],
            )

            report_row.save()


def Update_RefMap_DB(run_class: RunMain_class, parameter_set: ParameterSet):
    """
    Update Remap TABLES with info on this run.

    :param run_class:
    :return: run_data
    """
    print(f"updating refmap_dbs run {run_class.prefix}")

    user = User.objects.get(username=run_class.username)
    project = Projects.objects.get(
        name=run_class.project_name, owner=user, is_deleted=False
    )

    sample = PIProject_Sample.objects.get(
        name=run_class.sample.sample_name,
        project=project,
    )

    run = RunMain.objects.get(
        project=project,
        suprun=run_class.suprun,
        name=run_class.prefix,
        sample=sample,
        parameter_set=parameter_set,
    )

    for ref_map in run_class.remap_manager.mapped_instances:

        Update_ReferenceMap(ref_map, run, sample)


def Update_ReferenceMap(
    ref_map: Mapping_Instance,
    run: RunMain,
    sample: PIProject_Sample,
):
    """
    Updates the reference map data to TABLES.
    - ReferenceMap_Main,
    - ReferenceContigs
    """

    try:
        map_db = ReferenceMap_Main.objects.get(
            reference=ref_map.reference.target.acc_simple,
            taxid=ref_map.reference.target.taxid,
            sample=sample,
            run=run,
        )
    except ReferenceMap_Main.DoesNotExist:
        map_db = ReferenceMap_Main(
            reference=ref_map.reference.target.acc_simple,
            sample=sample,
            run=run,
            taxid=ref_map.reference.target.taxid,
            bam_file_path=ref_map.reference.read_map_sorted_bam,
            bai_file_path=ref_map.reference.read_map_sorted_bam_index,
            fasta_file_path=ref_map.reference.reference_file,
            fai_file_path=ref_map.reference.reference_fasta_index,
            mapped_subset_r1=ref_map.reference.mapped_subset_r1,
            mapped_subset_r2=ref_map.reference.mapped_subset_r2,
            mapped_subset_r1_fasta=ref_map.reference.mapped_subset_r1_fasta,
            mapped_subset_r2_fasta=ref_map.reference.mapped_subset_r2_fasta,
            vcf=ref_map.reference.vcf,
        )
        map_db.save()

    if ref_map.assembly is not None:

        remap_stats = ref_map.assembly.report.set_index("ID")

        for seqid, row in remap_stats.iterrows():
            try:
                map_db_seq = ReferenceContigs.objects.get(
                    reference=map_db, run=run, contig=seqid
                )
            except ReferenceContigs.DoesNotExist:
                map_db_seq = ReferenceContigs(
                    contig=seqid,
                    reference=map_db,
                    run=run,
                    depth=row["Hdepth"],
                    depthr=row["HdepthR"],
                    coverage=row["coverage"],
                )
                map_db_seq.save()

            # map_db_seq.report = run_class.report
            # map_db_seq.save()
