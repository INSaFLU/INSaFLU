import os
from typing import Type

from django.core.files import File

from pathogen_identification.models import (
    QC_REPORT,
    ContigClassification,
    FinalReport,
    Projects,
    ReadClassification,
    ReferenceContigs,
    ReferenceMap_Main,
    RunAssembly,
    RunDetail,
    RunIndex,
    RunMain,
    RunRemapMain,
    Sample,
    SampleQC,
)
from pathogen_identification.object_classes import Sample_runClass
from pathogen_identification.remap_class import Mapping_Instance
from pathogen_identification.run_main import RunMain_class


####################################################################################################################
####################################################################################################################
def Update_project(project_directory_path):
    """Updates the project"""
    project_directory_path = os.path.dirname(project_directory_path)
    project_name = os.path.basename(project_directory_path)
    project_name_simple = project_name.replace(".", "_").replace(":", "_")

    try:
        project = Projects.objects.get(name=project_name)

    except Projects.DoesNotExist:
        print("project_name: ", project_name)
        project = Projects(name=project_name, full_path=project_directory_path)
        project.save()


def Update_Sample(sample_class: Type[Sample_runClass]):
    """
    Update Sample class.

    :param sample_class:
    :return: None
    """

    try:
        Sample.objects.get(
            name_extended=sample_class.sample_name,
            project__name=sample_class.project_name,
        )
    except Sample.DoesNotExist:

        Update_sample(sample_class)

    sample = Sample.objects.get(
        name_extended=sample_class.sample_name,
        project__name=sample_class.project_name,
    )

    try:
        SampleQC.objects.get(sample=sample)
    except SampleQC.DoesNotExist:
        Update_sample_qc(sample_class)


def Update_sample(sample_class: Type[Sample_runClass]):
    """update sample_class.
    :param sample_class:
    :return: None
    """
    try:
        sample = Sample.objects.get(
            name_extended=sample_class.sample_name,
            project__name=sample_class.project_name,
        )

    except Sample.DoesNotExist:
        #
        project = Projects.objects.get(name=sample_class.project_name)

        sample = Sample(
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


def Update_sample_qc(sample_class: Type[Sample_runClass]):
    """update sample_class.qc_data.
    :param sample_class:
    :return: None
    """

    sample = Sample.objects.get(
        project__name=sample_class.project_name, name_extended=sample_class.sample_name
    )

    percent_passed = (
        int(sample_class.reads_after_processing)
        / int(sample_class.reads_before_processing)
    ) * 100

    input_report = open(sample_class.input_fastqc_report, "r")
    processed_report = open(sample_class.processed_fastqc_report, "r")

    print(sample_class.qc_soft)
    print(sample_class.technology)
    print(sample_class.reads_before_processing)
    print(sample_class.reads_after_processing)
    print(percent_passed)
    print(sample_class.qcdata)

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


def Update_QC_report(sample_class: Type[Sample_runClass]):
    """
    Update QC data for sample_class.

    :param sample_class:
    :return: None
    """
    sample = Sample.objects.get(
        project__name=sample_class.project_name, name_extended=sample_class.sample_name
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


def Update_Sample_Runs(run_class: Type[RunMain_class]):
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

    Update_RunMain(run_class)
    Update_Sample_Runs_DB(run_class)
    Update_RefMap_DB(run_class)


def retrieve_number_of_runs(project_name, sample_name):
    """
    retrieve number of runs for a given project.

    :param run_class:
    :return: number of runs
    """

    try:
        project = Projects.objects.get(name=project_name)
    except Projects.DoesNotExist:
        print(f"project {project_name} does not exist")
        return 0

    try:
        sample = Sample.objects.get(name_extended=sample_name, project=project)
    except Sample.DoesNotExist:
        return 0

    return RunMain.objects.filter(project=project, sample=sample).count() + 1


def RunIndex_Update_Retrieve_Key(project_name, sample_name):

    run_index = retrieve_number_of_runs(project_name, sample_name)

    new_name = f"run_{run_index}"

    project = Projects.objects.get(name=project_name)
    sample = Sample.objects.get(
        name_extended=sample_name,
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


def Update_RunMain(run_class: Type[RunMain_class]):
    """update run data for run_class. Update run_class.run_data.

    :param run_class:
    :return: None
    """
    project = Projects.objects.get(name=run_class.sample.project_name)
    sample = Sample.objects.get(
        name_extended=run_class.sample.sample_name,
        project__name=run_class.sample.project_name,
    )

    reads_after_processing = run_class.sample.reads_after_processing
    reads_proc_percent = (
        reads_after_processing / run_class.sample.reads_before_processing
    ) * 100

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
        )
    except RunMain.DoesNotExist:

        runmain = RunMain(
            suprun=run_class.suprun,
            project=project,
            sample=sample,
            name=run_class.prefix,
            params_file_path=run_class.params_file_path,
            processed_reads_r1=os.path.basename(run_class.sample.r1.current),
            processed_reads_r2=os.path.basename(run_class.sample.r2.current),
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
            # finished=str(run_class.finished),
            report="report",
            static_dir=run_class.static_dir,
        )

        runmain.save()


def Sample_update_combinations(run_class: Type[RunMain_class]):
    sample = Sample.objects.get(
        project__name=run_class.sample.project_name,
        name_extended=run_class.sample.sample_name,
    )

    sample.combinations = sample.combinations + 1

    sample.save()


def Update_Sample_Runs_DB(run_class: Type[RunMain_class]):
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
    Sample_update_combinations(run_class)

    sample = Sample.objects.get(
        project__name=run_class.sample.project_name,
        name_extended=run_class.sample.sample_name,
    )

    try:
        runmain = RunMain.objects.get(
            project__name=run_class.sample.project_name,
            suprun=run_class.suprun,
            sample=sample,
            name=run_class.prefix,
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


def Update_RefMap_DB(run_class: Type[RunMain_class]):
    """
    Update Remap TABLES with info on this run.

    :param run_class:
    :return: run_data
    """
    print(f"updating refmap_dbs run {run_class.prefix}")

    for ref_map in run_class.remap_manager.mapped_instances:

        Update_ReferenceMap(
            ref_map,
            run_class,
        )


def Update_ReferenceMap(
    ref_map: Type[Mapping_Instance],
    run_class: Type[RunMain_class],
):
    """
    Updates the reference map data to TABLES.
    - ReferenceMap_Main,
    - ReferenceContigs
    """
    sample = Sample.objects.get(
        name_extended=run_class.sample.sample_name,
        project__name=run_class.sample.project_name,
    )

    run = RunMain.objects.get(
        project__name=run_class.sample.project_name,
        suprun=run_class.suprun,
        name=run_class.prefix,
        sample=sample,
    )

    try:
        map_db = ReferenceMap_Main.objects.get(
            reference=ref_map.reference.target.acc_simple,
            sample=sample,
            run=run,
        )
    except ReferenceMap_Main.DoesNotExist:
        map_db = ReferenceMap_Main(
            reference=ref_map.reference.target.acc_simple,
            sample=sample,
            run=run,
            taxid=ref_map.reference.target.taxid,
            # report=ref_map.report,
            bam_file_path=ref_map.reference.read_map_sorted_bam,
            bai_file_path=ref_map.reference.read_map_sorted_bam_index,
            fasta_file_path=ref_map.reference.reference_file,
            fai_file_path=ref_map.reference.reference_fasta_index,
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
