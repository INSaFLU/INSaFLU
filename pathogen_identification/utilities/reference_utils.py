import ntpath
import os
from typing import List, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from django.conf import settings
from django.contrib.auth.models import User
from django.core.files.temp import NamedTemporaryFile
from django.db.models import Q

from constants.constants import Constants, FileExtensions, FileType, TypeFile, TypePath
from constants.software_names import SoftwareNames
from constants.televir_directories import Televir_Directory_Constants
from managing_files.models import ProcessControler
from managing_files.models import ProjectSample as InsafluProjectSample
from managing_files.models import Reference
from pathogen_identification.models import (
    RawReference,
    ReferenceSourceFileMap,
    TelefluMapping,
)
from pathogen_identification.utilities.televir_bioinf import TelevirBioinf
from pathogen_identification.utilities.utilities_general import simplify_name
from utils.software import Software
from utils.utils import Utils


def description_to_name(description, keep: int = 4):
    """
    This function takes the description of the reference and returns a name"""

    utils = Utils()

    description = description.split(" ")[:keep]
    description = "_".join(description)

    description = utils.clean_name(description)
    description = description.replace(",", "")
    return description


def fasta_from_reference(reference_id):
    """
    This function takes the raw reference id and returns the fasta name"""
    reference = RawReference.objects.get(id=reference_id)

    description_clean = description_to_name(reference.description)
    fasta_name = f"{reference.accid}_{description_clean}.fasta"

    return fasta_name


def replace_degenerate_bases(sequence):
    degenerate_bases = {
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
        # Exclude 'N' from replacement
    }

    for degenerate, options in degenerate_bases.items():
        if degenerate != "N":
            sequence = sequence.replace(degenerate, options[0])

    return sequence


def add_prefix(description, output_prefix):
    return f"{output_prefix}_{description}"


def add_suffix(description, output_suffix):
    return f"{description}_{output_suffix}"


def process_fasta(
    input_file,
    output_prefix: Optional[str] = None,
    output_suffix: Optional[str] = None,
):
    """
    This function takes a fasta file and replaces degenerate bases with the first base in the list of options.
    """
    records = list(SeqIO.parse(input_file, "fasta"))
    modified_records = []

    for record in records:
        original_sequence = str(record.seq)
        modified_sequence = replace_degenerate_bases(original_sequence)

        # Check if any replacements were made
        if modified_sequence != original_sequence:
            # Create a new SeqRecord with the modified identifier and sequence
            modified_id = record.description

            if output_prefix:
                modified_id = add_prefix(modified_id, output_prefix)
            if output_suffix:
                # Apply the suffix to the modified identifier
                modified_id = add_suffix(modified_id, output_suffix)

            modified_record = SeqRecord(
                seq=Seq(modified_sequence), id=modified_id, description=""
            )
            modified_records.append(modified_record)
        else:
            # If no replacements, keep the original record
            modified_records.append(record)

    with open(input_file, "w") as output_handle:
        SeqIO.write(modified_records, output_handle, "fasta")


def extract_file(accid):
    """ "
    This function takes the accid and returns the fasta file"""

    utils = Utils()
    televir_bioinf = TelevirBioinf()

    fasta_directory = Televir_Directory_Constants.ref_fasta_directory

    references = ReferenceSourceFileMap.objects.filter(reference_source__accid=accid)

    for reference in references:
        description = reference.reference_source.description
        description_simple = description_to_name(description)
        tmp_fasta = utils.get_temp_file(description_simple, ".fasta")

        source_file = os.path.join(
            fasta_directory, reference.reference_source_file.file
        )
        extracted = televir_bioinf.extract_reference(source_file, accid, tmp_fasta)
        if extracted:
            return tmp_fasta
        else:
            if os.path.exists(tmp_fasta):
                os.remove(tmp_fasta)


def merge_multiple_refs(references: List[RawReference], output_prefix: str):
    """
    This function takes a list of references and creates a merged fasta file
    """
    merged_fasta = NamedTemporaryFile(
        prefix=output_prefix, suffix=".fasta", delete=False
    )

    for reference in references:
        fasta_file = extract_file(reference.accid)
        process_fasta(fasta_file)
        with open(fasta_file, "r") as reference_fasta:
            merged_fasta.write(reference_fasta.read().encode())

    merged_fasta.flush()
    merged_fasta.close()

    return merged_fasta.name


from pathogen_identification.models import (
    MetaReference,
    RawReferenceMap,
    TeleFluProject,
)


def check_metaReference_exists(references: List[RawReference]):
    """
    This function takes a list of references and checks if a meta reference exists
    """

    meta_references = MetaReference.objects.filter(
        project__id=references[0].run.project.id
    )

    metaref_ids = [metaref.metaid for metaref in meta_references]
    reference_ids = [reference.id for reference in references]
    reference_ids_sorted = sorted(reference_ids)
    refid_code = "_".join([str(refid) for refid in reference_ids_sorted])

    if refid_code in metaref_ids:
        return True

    return False


def retrieve_metaReference(references: List[RawReference]):
    """
    This function takes a list of references and checks if a meta reference exists
    """
    meta_references = MetaReference.objects.filter(
        project__id=references[0].run.project.id
    )

    reference_ids = [reference.id for reference in references]
    reference_ids_sorted = sorted(reference_ids)
    refid_code = "_".join([str(refid) for refid in reference_ids_sorted])

    for metaref in meta_references:
        if metaref.metaid == refid_code:
            return metaref

    return None


def create_metaReference(references: List[RawReference]):
    if len(references) == 0:
        return None

    if check_metaReference_exists(references):
        return retrieve_metaReference(references)

    description = (
        f"{references[0].description} {references[0].accid}"
        if len(references) == 1
        else f"combined reference n{len(references)}"
    )
    metaref = MetaReference(
        description=description,
        project=references[0].run.project,
    )
    metaref.save()

    for reference in references:
        raw_ref_map = RawReferenceMap(
            reference=metaref,
            raw_reference=reference,
        )
        raw_ref_map.save()

    return metaref


import os


def check_metaReference_exists_from_ids(reference_ids: List[int]):
    references = RawReference.objects.filter(id__in=reference_ids)
    references = [reference for reference in references]

    return check_metaReference_exists(references)


def create_combined_reference(
    reference_ids: List[RawReference], output_prefix: str
) -> Optional[MetaReference]:
    """
    This function takes a list of references and creates a combined fasta file
    """
    utils = Utils()
    references = RawReference.objects.filter(id__in=reference_ids)
    references = [reference for reference in references]

    if check_metaReference_exists(references):
        return retrieve_metaReference(references)

    ### create combined fasta
    combined_fasta_name = merge_multiple_refs(references, output_prefix)

    ### create meta reference
    metaref = create_metaReference(references)

    if metaref is None:
        return None

    user_id = references[0].run.project.owner.id

    ### Create reference filename
    combined_fasta_name_final = f"{output_prefix}_{metaref.pk}.fasta"

    ### move the files to the right place
    final_data_path = os.path.join(
        settings.MEDIA_ROOT,
        utils.get_path_to_teleflu_reference_file(user_id, metaref.id),
    )

    sz_file_to = os.path.join(
        final_data_path,
        combined_fasta_name_final,
    )
    os.makedirs(final_data_path, exist_ok=True)

    utils.copy_file(
        combined_fasta_name,
        sz_file_to,
    )

    metaref.file_path = sz_file_to
    metaref.save()

    return metaref


def temp_fasta_copy(fasta_filepath: str):
    reference_fasta_temp_file_name = NamedTemporaryFile(prefix="flu_fa_", delete=False)
    with open(fasta_filepath, "r") as reference_fasta:
        reference_fasta_temp_file_name.write(reference_fasta.read().encode())

    reference_fasta_temp_file_name.flush()
    reference_fasta_temp_file_name.close()

    return reference_fasta_temp_file_name.name


def create_genbank_for_fasta(fasta_filepath: str):
    software = Software()
    utils = Utils()

    reference_fasta_temp_file_name = temp_fasta_copy(fasta_filepath)
    reference_fasta_temp_file_name = open(reference_fasta_temp_file_name, "r")

    software.dos_2_unix(reference_fasta_temp_file_name.name)
    software.fasta_2_upper(reference_fasta_temp_file_name.name)
    process_fasta(reference_fasta_temp_file_name.name)

    file_name_cleaned = utils.clean_name(
        ntpath.basename(reference_fasta_temp_file_name.name)
    )
    try:
        temp_genbank_dir = software.run_prokka(
            reference_fasta_temp_file_name.name, file_name_cleaned
        )
    except Exception as e:
        os.unlink(reference_fasta_temp_file_name.name)

        return None

    os.unlink(reference_fasta_temp_file_name.name)

    temp_genbank_file = os.path.join(
        temp_genbank_dir,
        utils.clean_extension(file_name_cleaned) + FileExtensions.FILE_GBK,
    )

    if not os.path.exists(temp_genbank_file):
        utils.remove_dir(temp_genbank_dir)

        return None

    return temp_genbank_file


def check_reference_exists(raw_reference_id, user_id):
    raw_ref = RawReference.objects.get(id=raw_reference_id)

    description = raw_ref.description
    description_clean = description_to_name(description)
    accid = raw_ref.accid

    query_set = Reference.objects.filter(
        owner__id=user_id, is_obsolete=False, is_deleted=False
    ).order_by("-name")

    if query_set.filter(
        Q(name__icontains=description_clean)
        | Q(reference_genbank_name__icontains=accid)
        | Q(reference_fasta_name__icontains=accid)
    ).exists():
        return True

    return False


def delete_reference(raw_reference_id, user_id):
    raw_ref = RawReference.objects.get(id=raw_reference_id)

    description = raw_ref.description
    description_clean = description_to_name(description)
    accid = raw_ref.accid

    query_set = Reference.objects.filter(
        owner__id=user_id, is_obsolete=False, is_deleted=False
    ).order_by("-name")

    existing = query_set.filter(
        Q(name__icontains=description_clean)
        | Q(reference_genbank_name__icontains=accid)
        | Q(reference_fasta_name__icontains=accid)
    )

    if existing.exists():
        existing.delete()


def check_reference_submitted(ref_id, user_id):
    user = User.objects.get(pk=user_id)
    process_controler = ProcessControler()

    process = ProcessControler.objects.filter(
        owner__id=user.pk,
        name=process_controler.get_name_televir_teleflu_ref_create(
            ref_id=ref_id,
        ),
        is_finished=False,
        is_error=False,
    )
    return process.exists()


def raw_reference_to_insaflu(raw_reference_id: int, user_id: int):
    user = User.objects.get(id=user_id)
    raw_reference = RawReference.objects.get(id=raw_reference_id)
    accid = raw_reference.accid

    if check_reference_exists(raw_reference_id, user_id):
        return False, None

    reference_fasta = extract_file(accid)

    name = description_to_name(raw_reference.description)
    final_fasta_name = fasta_from_reference(raw_reference_id)

    if reference_fasta is None:
        return None, None

    success, reference_id = generate_insaflu_reference(
        reference_fasta, name, final_fasta_name, user
    )

    return success, reference_id


def teleflu_to_insaflu_reference(project_id: int, user_id: int):
    user = User.objects.get(id=user_id)
    teleflu_project = TeleFluProject.objects.get(id=project_id)
    metareference = teleflu_project.raw_reference
    reference_fasta = teleflu_project.raw_reference.file_path
    name = metareference.description

    success, reference_id = generate_insaflu_reference(
        reference_fasta, name, os.path.basename(reference_fasta), user
    )

    insaflu_reference = Reference.objects.get(id=reference_id)
    teleflu_project.reference = insaflu_reference
    teleflu_project.save()

    return success, insaflu_reference


def generate_insaflu_reference(
    reference_fasta: str, name: str, final_fasta_name: str, user: User
):
    utils = Utils()
    software = Software()
    final_gb_name = final_fasta_name.replace(".fasta", ".gbk")

    temp_genbank_file = create_genbank_for_fasta(reference_fasta)

    if reference_fasta is None:
        return None, None

    ### Create reference

    number_locus = utils.is_fasta(reference_fasta)

    reference = Reference()
    reference.name = name
    reference.owner = user
    reference.is_obsolete = False
    reference.number_of_locus = number_locus
    reference.reference_fasta_name = final_fasta_name
    reference.reference_genbank_name = final_gb_name
    reference.save()

    ## move the files to the right place
    final_data_path = os.path.join(
        settings.MEDIA_ROOT, utils.get_path_to_reference_file(user.id, reference.id)
    )
    os.makedirs(final_data_path, exist_ok=True)

    ## fasta file
    sz_file_to = os.path.join(
        final_data_path,
        reference.reference_fasta_name,
    )
    cwd = os.getcwd()
    os.chdir(os.path.join(Constants.TEMP_DIRECTORY, Constants.COUNT_DNA_TEMP_DIRECTORY))

    software.dos_2_unix(reference_fasta)
    ## test if bases all lower
    software.fasta_2_upper(reference_fasta)
    utils.copy_file(
        reference_fasta,
        sz_file_to,
    )

    reference.hash_reference_fasta = utils.md5sum(sz_file_to)
    # fill reference_fasta ContentTypeRestrictedFileField
    with open(sz_file_to, "rb") as f:
        reference.reference_fasta.save(os.path.basename(sz_file_to), f, save=True)
    reference.reference_fasta.name = os.path.join(
        utils.get_path_to_reference_file(user.id, reference.id),
        reference.reference_fasta_name,
    )

    ###
    ## genbank file
    sz_file_to = os.path.join(
        final_data_path,
        reference.reference_genbank_name,
    )

    utils.move_file(temp_genbank_file, sz_file_to)

    software.dos_2_unix(sz_file_to)
    reference.hash_reference_genbank = utils.md5sum(sz_file_to)
    # fill reference_genbank ContentTypeRestrictedFileField
    with open(sz_file_to, "rb") as f:
        reference.reference_genbank.save(os.path.basename(sz_file_to), f, save=True)

    reference.reference_genbank.name = os.path.join(
        utils.get_path_to_reference_file(user.id, reference.id),
        reference.reference_genbank_name,
    )
    reference.save()

    ### create bed and index for genbank
    utils.from_genbank_to_bed(
        sz_file_to, reference.get_reference_bed(TypePath.MEDIA_ROOT)
    )
    software.create_index_files_from_igv_tools(
        reference.get_reference_bed(TypePath.MEDIA_ROOT)
    )

    ### create some gff3  essential to run other tools

    software.run_genbank2gff3(sz_file_to, reference.get_gff3(TypePath.MEDIA_ROOT))
    software.run_genbank2gff3(
        sz_file_to,
        reference.get_gff3_with_gene_annotation(TypePath.MEDIA_ROOT),
        True,
    )
    software.run_genbank2gff3_positions_comulative(
        sz_file_to, reference.get_gff3_comulative_positions(TypePath.MEDIA_ROOT)
    )

    os.chdir(cwd)
    ### save in database the elements and coordinates
    utils.get_elements_from_db(reference, user)
    utils.get_elements_and_cds_from_db(reference, user)

    ## create the index before commit in database, throw exception if something goes wrong
    software.create_fai_fasta(
        os.path.join(
            getattr(settings, "MEDIA_ROOT", None), reference.reference_fasta.name
        )
    )

    ### set specie tag
    software.get_species_tag(reference)

    return True, reference.id


def create_teleflu_igv_report(teleflu_project_pk: int) -> bool:

    teleflu_project = TeleFluProject.objects.get(pk=teleflu_project_pk)
    insaflu_project = teleflu_project.insaflu_project

    ### get reference
    reference = teleflu_project.reference
    if reference is None:
        return False

    reference_file = reference.get_reference_fasta(TypePath.MEDIA_ROOT)

    samples = InsafluProjectSample.objects.filter(project=insaflu_project)
    # samples= [sample.sample for sample in samples]

    sample_dict = {}

    ### get sample files
    software_names = SoftwareNames()

    for sample in samples:

        if sample.sample.type_of_fastq == 0:
            filename = software_names.get_snippy_name()
        else:
            filename = software_names.get_medaka_name()

        bam_file = sample.get_file_output(
            TypePath.MEDIA_ROOT, FileType.FILE_BAM, filename
        )
        bam_file_index = sample.get_file_output(
            TypePath.MEDIA_ROOT,
            FileType.FILE_BAM_BAI,
            filename,
        )
        vcf_file = sample.get_file_output(
            TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, filename
        )

        if bam_file and bam_file_index and os.path.exists(vcf_file):
            sample_dict[sample.sample.pk] = {
                "name": sample.sample.name,
                "bam_file": bam_file,
                "bam_file_index": bam_file_index,
                "vcf_file": vcf_file,
            }

    ### merge vcf files
    televir_bioinf = TelevirBioinf()
    vcf_files = [files["vcf_file"] for sample_pk, files in sample_dict.items()]
    group_vcf = teleflu_project.project_vcf
    stacked_html = teleflu_project.project_igv_report_media

    os.makedirs(teleflu_project.project_vcf_directory, exist_ok=True)

    merged_success = televir_bioinf.merge_vcf_files(vcf_files, group_vcf)

    try:

        televir_bioinf.create_igv_report(
            reference_file,
            vcf_file=group_vcf,
            tracks=sample_dict,
            output_html=stacked_html,
        )

        return True
    except Exception as e:
        print(e)
        return False


from pathogen_identification.models import (
    ParameterSet,
    PIProject_Sample,
    ReferenceMap_Main,
)


def filter_reference_maps_select(
    sample: PIProject_Sample, leaf_id: int, reference: List[str]
) -> Optional[ReferenceMap_Main]:

    ref_maps = ReferenceMap_Main.objects.filter(
        sample=sample,
        run__parameter_set__leaf__pk=leaf_id,
        reference__in=reference,
    )

    refs = RawReference.objects.filter(
        run__parameter_set__sample=sample,
        accid__in=reference,
        run__parameter_set__leaf__pk=leaf_id,
        run__parameter_set__status=ParameterSet.STATUS_FINISHED,
    )

    for ref in ref_maps:

        print("############### ref")
        print(ref)
        print(ref.bam_file_path)
        print(ref.bai_file_path)

        if not ref.bam_file_path:
            continue
        if not ref.bai_file_path:
            continue

        if os.path.exists(ref.bam_file_path) == False:
            continue

        if os.path.exists(ref.vcf) == False:
            continue

        if os.path.exists(ref.bai_file_path) == False:
            continue

        return ref

    return None


def create_televir_igv_report(teleflu_project_pk: int, leaf_index: int) -> bool:

    teleflu_project = TeleFluProject.objects.get(pk=teleflu_project_pk)
    teleflu_mapping = TelefluMapping.objects.get(
        teleflu_project=teleflu_project, leaf__pk=leaf_index
    )

    print("teleflu_mapping", teleflu_mapping)
    # reference_accid= teleflu_project.raw_reference.

    ### get reference insaflu
    insaflu_reference = teleflu_project.reference
    print("insaflu_reference", insaflu_reference)
    if insaflu_reference is None:
        return False

    reference_file = insaflu_reference.get_reference_fasta(TypePath.MEDIA_ROOT)

    # televir_reference
    teleflu_refs = teleflu_project.televir_references
    print("teleflu_refs", teleflu_refs)

    if teleflu_refs is None:
        return False

    accid_list = [ref.accid for ref in teleflu_refs if ref.accid]
    accid_list_simple = [simplify_name(accid) for accid in accid_list]

    # samples
    televir_project_samples = teleflu_mapping.mapped_samples
    sample_dict = {}

    ### get sample files

    print("televir_project_samples", televir_project_samples)

    for sample in televir_project_samples:

        ref_select = filter_reference_maps_select(sample, leaf_index, accid_list_simple)

        if ref_select is None:
            continue

        sample_dict[sample.pk] = {
            "name": sample.name,
            "bam_file": ref_select.bam_file_path,
            "bam_file_index": ref_select.bai_file_path,
            "vcf_file": ref_select.vcf,
            "sample": sample,
        }

    ### merge vcf files
    print("sample_dict", sample_dict)
    if len(sample_dict) == 0:
        return False
    else:
        os.makedirs(teleflu_mapping.mapping_directory, exist_ok=True)

    televir_bioinf = TelevirBioinf()
    # vcf_files = [files["vcf_file"] for sample_pk, files in sample_dict.items()]
    group_vcf = teleflu_mapping.mapping_vcf
    stacked_html = teleflu_mapping.mapping_igv_report

    os.makedirs(teleflu_project.project_vcf_directory, exist_ok=True)

    # merged_success = televir_bioinf.merge_vcf_files(vcf_files, group_vcf)

    try:
        merged_success = televir_bioinf.vcf_from_bam(
            [files["bam_file"] for sample_pk, files in sample_dict.items()],
            reference_file,
            group_vcf,
        )

        print("merged_success", merged_success)

        for sample_pk, sample_info in sample_dict.items():
            sample = PIProject_Sample.objects.get(pk=sample_pk)
            teleflu_mapping.stacked_samples.add(sample)
            teleflu_mapping.save()

    except Exception as e:
        print(e)
        return False

    try:

        televir_bioinf.create_igv_report(
            reference_file,
            vcf_file=group_vcf,
            tracks=sample_dict,
            output_html=stacked_html,
        )

        return True
    except Exception as e:
        print(e)
        return False
