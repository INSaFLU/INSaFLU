
from utils.utils import Utils
from pathogen_identification.models import ReferenceSourceFileMap, ReferenceSourceFile, ReferenceSource, RawReference
from django.contrib.auth.models import User
from pathogen_identification.utilities.televir_bioinf import TelevirBioinf
from constants.televir_directories import Televir_Directory_Constants
import os
from managing_files.models import Reference
from django.core.files.temp import NamedTemporaryFile
from utils.software import Software
import ntpath
from constants.constants import (Constants, FileExtensions, FileType, TypeFile,
                                 TypePath)
from django.conf import settings
from django.db.models import Q



def description_to_name(description):
    """
    This function takes the description of the reference and returns a name"""

    utils= Utils()

    description= description.split(" ")[:4]
    description= "_".join(description)

    description= utils.clean_name(description)
    description= description.replace(",", "")
    return description


def fasta_from_reference(reference_id):
    """
    This function takes the raw reference id and returns the fasta name"""
    reference= RawReference.objects.get(id=reference_id)

    description_clean= description_to_name(reference.description)
    fasta_name= f"{reference.accid}_{description_clean}.fasta"

    return fasta_name


def extract_file(accid):
    """"
    This function takes the accid and returns the fasta file"""

    utils= Utils()
    televir_bioinf= TelevirBioinf()
    
    fasta_directory= Televir_Directory_Constants.ref_fasta_directory
    
    references=ReferenceSourceFileMap.objects.filter(reference_source__accid= accid)

    for reference in references:

        description= reference.reference_source.description
        description_simple= description_to_name(description)
        tmp_fasta= utils.get_temp_file(description_simple,".fasta")


        source_file= os.path.join(
            fasta_directory,
            reference.reference_source_file.file
        )
        extracted= televir_bioinf.extract_reference(source_file, accid, tmp_fasta)
        if extracted:
            return tmp_fasta
        else:
            if os.path.exists(tmp_fasta):
                os.remove(tmp_fasta)
    

def temp_fasta_copy(fasta_filepath: str):
    reference_fasta_temp_file_name = NamedTemporaryFile(
        prefix="flu_fa_", delete=False
    )
    with open(fasta_filepath, "r") as reference_fasta:
            reference_fasta_temp_file_name.write(reference_fasta.read().encode())

    reference_fasta_temp_file_name.flush()
    reference_fasta_temp_file_name.close()
    

    return reference_fasta_temp_file_name.name

def create_genbank_for_fasta(fasta_filepath: str):
        
        software= Software()
        utils= Utils()

        reference_fasta_temp_file_name = temp_fasta_copy(fasta_filepath)
        reference_fasta_temp_file_name = open(reference_fasta_temp_file_name, "r")

        software.dos_2_unix(reference_fasta_temp_file_name.name)
        software.fasta_2_upper(reference_fasta_temp_file_name.name)

        file_name_cleaned = utils.clean_name(ntpath.basename(reference_fasta_temp_file_name.name))
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

    raw_ref= RawReference.objects.get(id=raw_reference_id)

    description= raw_ref.description
    description_clean= description_to_name(description)
    accid= raw_ref.accid

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



def reference_to_teleflu(raw_reference_id: int, user_id: int):

    utils= Utils()
    software= Software()

    user= User.objects.get(id=user_id)
    raw_reference= RawReference.objects.get(id=raw_reference_id)
    accid= raw_reference.accid
    
    

    if check_reference_exists(raw_reference_id, user_id):
        return False, None
    
    reference_fasta= extract_file(accid)

    name= description_to_name(raw_reference.description)
    final_fasta_name= fasta_from_reference(raw_reference_id)
    final_gb_name= final_fasta_name.replace(".fasta", ".gbk")

    if reference_fasta is None:
        return None, None

    temp_genbank_file= create_genbank_for_fasta(reference_fasta)

    if reference_fasta is None:
        return None, None

    ### Create reference

    number_locus = utils.is_fasta(reference_fasta)

    reference= Reference()
    reference.name= name
    reference.owner= user
    reference.is_obsolete= False
    reference.number_of_locus= number_locus
    reference.reference_fasta_name= final_fasta_name
    reference.reference_genbank_name= final_gb_name
    reference.save()

    ## move the files to the right place
    final_data_path=  os.path.join(
        settings.MEDIA_ROOT,
        utils.get_path_to_reference_file(user.id, reference.id)
    )

    ## fasta file
    sz_file_to = os.path.join(
        final_data_path,
        reference.reference_fasta_name,
    )
    
    software.dos_2_unix(
        reference_fasta
    )
    ## test if bases all lower
    software.fasta_2_upper(
        reference_fasta
    )
    utils.move_file(
        reference_fasta,
        sz_file_to,
    )

    reference.hash_reference_fasta = utils.md5sum(sz_file_to)
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
    


