"""
Created on Dec 6, 2017

@author: mmp
"""

import csv
import json
import logging
import os
from datetime import datetime

from constants.constants import (Constants, FileExtensions, FileType, TypeFile,
                                 TypePath)
from constants.meta_key_and_values import MetaKeyAndValue
from constants.software_names import SoftwareNames
from django.conf import settings
from django.db import transaction
from django.http import JsonResponse
from django.utils.safestring import mark_safe
from django.utils.translation import ugettext_lazy as _
from django.views.decorators.csrf import csrf_protect
from extend_user.models import Profile
from pathogen_identification.models import PIProject_Sample
from pathogen_identification.models import Projects as Televir_Project
from settings.constants_settings import ConstantsSettings
from settings.default_parameters import DefaultParameters
from settings.default_software_project_sample import DefaultProjectSoftware
from utils.collect_extra_data import CollectExtraData
from utils.process_SGE import ProcessSGE
from utils.result import Coverage, DecodeObjects
from utils.software import Software
from utils.utils import Utils

from managing_files.manage_database import ManageDatabase
from managing_files.models import (DataSet, MetaKey, ProcessControler, Project,
                                   ProjectSample, Reference, Sample,
                                   UploadFiles, VaccineStatus)

### Logger
logger_debug = logging.getLogger("fluWebVirus.debug")
logger_production = logging.getLogger("fluWebVirus.production")


######################################
###
###        AJAX methods for check box in session
###


@csrf_protect
def set_check_box_values(request):
    """
    manage check boxes through ajax
    """
    if request.is_ajax():
        data = {"is_ok": False}
        utils = Utils()
        if Constants.CHECK_BOX_ALL in request.GET:
            request.session[Constants.CHECK_BOX_ALL] = utils.str2bool(
                request.GET.get(Constants.CHECK_BOX_ALL)
            )
            ## check all unique
            for key in request.session.keys():
                if (
                    key.startswith(Constants.CHECK_BOX)
                    and len(key.split("_")) == 3
                    and utils.is_integer(key.split("_")[2])
                ):
                    request.session[key] = request.session[Constants.CHECK_BOX_ALL]
            data = {"is_ok": True}
        ### one check box is pressed
        elif Constants.CHECK_BOX in request.GET:
            request.session[Constants.CHECK_BOX_ALL] = False  ### set All false
            key_name = "{}_{}".format(
                Constants.CHECK_BOX, request.GET.get(Constants.CHECK_BOX_VALUE)
            )
            request.session[key_name] = utils.str2bool(
                request.GET.get(Constants.CHECK_BOX)
            )
            total_checked = 0
            for key in request.session.keys():
                if (
                    key.startswith(Constants.CHECK_BOX)
                    and len(key.split("_")) == 3
                    and utils.is_integer(key.split("_")[2])
                ):
                    if request.session[key]:
                        total_checked += 1
            data = {
                "is_ok": True,
                "total_checked": total_checked,
            }
        ## get the status of a check_box_single
        elif Constants.GET_CHECK_BOX_SINGLE in request.GET:
            data["is_ok"] = True
            for key in request.session.keys():
                if (
                    key.startswith(Constants.CHECK_BOX)
                    and len(key.split("_")) == 3
                    and utils.is_integer(key.split("_")[2])
                ):
                    data[key] = request.session[key]
        ## change single status of a check_box_single
        elif Constants.GET_CHANGE_CHECK_BOX_SINGLE in request.GET:
            data["is_ok"] = True
            key_name = "{}_{}".format(
                Constants.CHECK_BOX, request.GET.get(Constants.CHECK_BOX_VALUE)
            )
            for key in request.session.keys():
                if (
                    key.startswith(Constants.CHECK_BOX)
                    and len(key.split("_")) == 3
                    and utils.is_integer(key.split("_")[2])
                ):
                    if request.session[key]:
                        data[key] = False
                    if key == key_name:
                        request.session[key] = utils.str2bool(
                            request.GET.get(Constants.GET_CHANGE_CHECK_BOX_SINGLE)
                        )
                    else:
                        request.session[key] = False
        elif Constants.COUNT_CHECK_BOX in request.GET:
            count = 0
            for key in request.session.keys():
                if (
                    key.startswith(Constants.CHECK_BOX)
                    and len(key.split("_")) == 3
                    and utils.is_integer(key.split("_")[2])
                ):
                    if request.session[key]:
                        count += 1
                elif (
                    key == Constants.CHECK_BOX_ALL
                    and request.session[Constants.CHECK_BOX_ALL]
                ):
                    count += 1
            data = {
                "is_ok": True,
                Constants.COUNT_CHECK_BOX: count,
            }
        return JsonResponse(data)


###
###        END AJAX methods for check box in session
###
######################################


@csrf_protect
def show_phylo_canvas(request):
    """
    manage check boxes through ajax
    """

    if request.is_ajax():
        data = {"is_ok": False}
        utils = Utils()
        key_with_project_id = "project_id"
        if key_with_project_id in request.GET:
            project_id = int(request.GET.get(key_with_project_id))
            element_name = "all_together"
            key_element_name = "key_element_name"
            if key_element_name in request.GET:
                element_name = request.GET.get(key_element_name)
            try:
                project = Project.objects.get(id=project_id)
                file_name_root_json = project.get_global_file_by_project(
                    TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_SAMPLE_RESULT_json
                )
                file_name_url_json = project.get_global_file_by_project(
                    TypePath.MEDIA_URL, Project.PROJECT_FILE_NAME_SAMPLE_RESULT_json
                )
                ### this is a little version of PROJECT_FILE_NAME_SAMPLE_RESULT_CSV
                file_name_root_sample = project.get_global_file_by_project(
                    TypePath.MEDIA_ROOT,
                    Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV_simple,
                )
                if (
                    not os.path.exists(file_name_root_sample)
                    or os.path.getsize(file_name_root_sample) == 0
                ):  ## need to be created
                    collect_extra_data = CollectExtraData()
                    collect_extra_data.calculate_global_files(
                        Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV_simple,
                        project,
                        project.owner,
                    )

                if element_name == "all_together":
                    file_name_root_nwk = project.get_global_file_by_project(
                        TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_FASTTREE_tree
                    )
                    file_name_nwk = project.get_global_file_by_project(
                        TypePath.MEDIA_URL, Project.PROJECT_FILE_NAME_FASTTREE
                    )
                    file_name_tree = project.get_global_file_by_project(
                        TypePath.MEDIA_URL, Project.PROJECT_FILE_NAME_FASTTREE_tree
                    )
                else:
                    file_name_root_nwk = project.get_global_file_by_element(
                        TypePath.MEDIA_ROOT,
                        element_name,
                        Project.PROJECT_FILE_NAME_FASTTREE_tree,
                    )
                    file_name_nwk = project.get_global_file_by_element(
                        TypePath.MEDIA_URL,
                        element_name,
                        Project.PROJECT_FILE_NAME_FASTTREE,
                    )
                    file_name_tree = project.get_global_file_by_element(
                        TypePath.MEDIA_URL,
                        element_name,
                        Project.PROJECT_FILE_NAME_FASTTREE_tree,
                    )

                if os.path.exists(file_name_root_nwk) and os.path.exists(
                    file_name_root_sample
                ):
                    string_file_content = utils.read_file_to_string(
                        file_name_root_nwk
                    ).strip()

                    if (
                        not os.path.exists(file_name_root_json)
                        or os.path.getsize(file_name_root_json) == 0
                    ):
                        with open(
                            file_name_root_json, "w", encoding="utf-8"
                        ) as handle_write, open(file_name_root_sample) as handle_in_csv:
                            reader = csv.DictReader(handle_in_csv)
                            all_data = json.loads(json.dumps(list(reader)))
                            dt_result = {}
                            for dict_data in all_data:
                                if "id" in dict_data:
                                    dt_out = dict_data.copy()
                                    del dt_out["id"]
                                    dt_result[dict_data["id"]] = dt_out
                            if len(dt_result) == len(all_data):
                                handle_write.write(json.dumps(dt_result))
                            else:
                                logger_production.error(
                                    "ProjectID: {}  different number of lines processing Sample {} -> JSON {}".format(
                                        project_id, len(dt_result), len(all_data)
                                    )
                                )
                                logger_debug.error(
                                    "ProjectID: {}  different number of lines processing Sample {} -> JSON {}".format(
                                        project_id, len(dt_result), len(all_data)
                                    )
                                )
                                string_file_content = None  ## return error

                    if string_file_content != None and len(string_file_content) > 0:
                        data["is_ok"] = True
                        data["tree"] = string_file_content
                        data["root"] = project.reference.name
                        data["url_sample"] = file_name_url_json
                        data["tree_nwk_id"] = mark_safe(
                            '<strong>Tree (.nwk):</strong> <a href="{}" download="{}"> {}</a>'.format(
                                file_name_nwk,
                                os.path.basename(file_name_nwk),
                                os.path.basename(file_name_nwk),
                            )
                        )
                        data["tree_tree_id"] = mark_safe(
                            '<strong>Tree (.tree):</strong> <a href="{}" download="{}"> {}</a>'.format(
                                file_name_tree,
                                os.path.basename(file_name_tree),
                                os.path.basename(file_name_tree),
                            )
                        )
            except Project.DoesNotExist:
                pass
        return JsonResponse(data)


@csrf_protect
def show_variants_as_a_table(request):
    """
    return table with variants
    """
    if request.is_ajax():
        data = {"is_ok": False}
        key_with_project_id = "project_id"
        if key_with_project_id in request.GET:
            project_id = int(request.GET.get(key_with_project_id))
            try:
                project = Project.objects.get(id=project_id)
                out_file = project.get_global_file_by_project(
                    TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_TAB_VARIATIONS_SNIPPY
                )
                if os.path.exists(out_file) and os.stat(out_file).st_size > 0:
                    data["is_ok"] = True
                    data["url_path_variant_table"] = mark_safe(
                        request.build_absolute_uri(
                            project.get_global_file_by_project(
                                TypePath.MEDIA_URL,
                                Project.PROJECT_FILE_NAME_TAB_VARIATIONS_SNIPPY,
                            )
                        )
                    )
                    data["static_table_filter"] = mark_safe(
                        request.build_absolute_uri(
                            os.path.join(settings.STATIC_URL, "vendor/tablefilter")
                        )
                    )
            except Project.DoesNotExist:
                pass
        return JsonResponse(data)


@csrf_protect
def show_aln2pheno(request):
    """
    return table with variants
    """
    if request.is_ajax():
        data = {"is_ok": False}
        key_with_project_id = "project_id"
        if key_with_project_id in request.GET:
            project_id = int(request.GET.get(key_with_project_id))
            try:
                project = Project.objects.get(id=project_id)
                out_file = project.get_global_file_by_project(
                    TypePath.MEDIA_ROOT,
                    Project.PROJECT_FILE_NAME_Aln2pheno_report_COG_UK,
                )
                if os.path.exists(out_file) and os.stat(out_file).st_size > 0:
                    data["is_ok"] = True
                    data["url_path_aln2pheno"] = mark_safe(
                        request.build_absolute_uri(
                            project.get_global_file_by_project(
                                TypePath.MEDIA_URL,
                                Project.PROJECT_FILE_NAME_Aln2pheno_report_COG_UK,
                            )
                        )
                    )
                    data["static_table_filter"] = mark_safe(
                        request.build_absolute_uri(
                            os.path.join(settings.STATIC_URL, "vendor/tablefilter")
                        )
                    )
            except Project.DoesNotExist:
                pass
        return JsonResponse(data)


@csrf_protect
def show_coverage_as_a_table(request):
    """
    return table with coverage
    """

    if request.is_ajax():
        data = {"is_ok": False}
        key_with_project_id = "project_id"
        key_client_width = "client_width"
        if key_with_project_id in request.GET:
            project_id = int(request.GET.get(key_with_project_id))
            client_width = int(request.GET.get(key_client_width)) - 300
            try:
                manageDatabase = ManageDatabase()
                default_software = DefaultProjectSoftware()
                utils = Utils()
                project = Project.objects.get(id=project_id)

                ### get all elements and gene names
                geneticElement = utils.get_elements_and_genes(
                    project.reference.get_reference_gbk(TypePath.MEDIA_ROOT)
                )

                size_elements = client_width // len(
                    geneticElement.get_sorted_elements()
                )
                size_samples_max = 160
                if size_elements > size_samples_max:
                    size_elements = size_samples_max
                len(geneticElement.get_sorted_elements())
                ## this line is passed at the end
                ## content = '<table id="table_with_coverage_variants_id" number_coluns={} number_rows={}>\n<thead>\n<tr>'
                content = '<td id="id_table-coverage_0_0" class="header-table" size="{}">Samples</td>'.format(
                    size_elements
                )
                count_sequences = 1
                for sequence_name in geneticElement.get_sorted_elements():
                    content += '<td id="id_table-coverage_0_{}" class="header-table" size="{}">{}</td>'.format(
                        count_sequences, size_elements, sequence_name
                    )
                    count_sequences += 1

                content += "</thead></tr><tbody>"
                count_projects = 1
                for project_sample in project.project_samples.all():
                    if not project_sample.get_is_ready_to_proccess():
                        continue
                    content += (
                        '<tr class="coverage_image_tr_sample">'
                        + '<td id="id_table-coverage_{}_0" class="table-header-coverage">'.format(
                            count_projects
                        )
                        + project_sample.sample.name
                        + "</td>"
                    )

                    meta_value = manageDatabase.get_project_sample_metakey_last(
                        project_sample,
                        MetaKeyAndValue.META_KEY_Coverage,
                        MetaKeyAndValue.META_VALUE_Success,
                    )
                    decode_coverage = DecodeObjects()
                    coverage = decode_coverage.decode_result(meta_value.description)

                    ### default parameters
                    limit_to_mask_consensus = int(
                        default_software.get_mask_consensus_single_parameter(
                            project_sample,
                            DefaultParameters.MASK_CONSENSUS_threshold,
                            ConstantsSettings.TECHNOLOGY_illumina
                            if project_sample.is_sample_illumina()
                            else ConstantsSettings.TECHNOLOGY_minion,
                        )
                    )

                    count_sequences = 1
                    for sequence_name in geneticElement.get_sorted_elements():
                        coverage_value_average = coverage.get_coverage(
                            sequence_name, Coverage.COVERAGE_ALL
                        )
                        coverage_value = (
                            coverage.get_coverage(
                                sequence_name, Coverage.COVERAGE_MORE_DEFINED_BY_USER
                            )
                            if coverage.is_exist_limit_defined_by_user()
                            else coverage.get_coverage(
                                sequence_name, Coverage.COVERAGE_MORE_9
                            )
                        )
                        href_sample = '<a href="#coverageModal" id="id_table-coverage_{}_{}" data-toggle="modal" class="tip" project_sample_id="{}" sequence="{}" title="{}"></a>'.format(
                            count_projects,
                            count_sequences,
                            project_sample.id,
                            sequence_name,
                            coverage.get_message_to_show_in_web_site(
                                project_sample.sample.name, sequence_name
                            ),
                        )
                        content += '<td id="id_table-coverage_content_{}_{}" class="table-coverage-image" value_data="{}" '.format(
                            count_projects, count_sequences, coverage_value
                        ) + 'value_data_average="{}" color_graphic="{}" size="{}" value_limit_coverage="{}">{}</td>'.format(
                            coverage_value_average,
                            coverage.get_color(sequence_name, limit_to_mask_consensus),
                            size_elements,
                            coverage.get_middle_limit(),
                            href_sample,
                        )
                        count_sequences += 1

                    count_projects += 1
                    content += "</tr>"
                content += "</tbody></table>"

                data["is_ok"] = True
                data["content"] = (
                    '<table id="table_with_coverage_variants_id" number_coluns={} number_rows={}>\n<thead>\n<tr>'.format(
                        count_sequences - 1, count_projects - 1
                    )
                    + content
                )
            except Project.DoesNotExist:
                pass
        return JsonResponse(data)


@csrf_protect
def show_msa_nucleotide(request):
    """
    manage msa nucleotide alignments
    """
    if request.is_ajax():
        software = Software()
        utils = Utils()
        data = {"is_ok": False}
        key_with_project_id = "project_id"
        if key_with_project_id in request.GET:
            project_id = int(request.GET.get(key_with_project_id))
            element_name = "all_together"
            key_element_name = "key_element_name"
            if key_element_name in request.GET:
                element_name = request.GET.get(key_element_name)
            try:
                manage_database = ManageDatabase()
                project = Project.objects.get(id=project_id)

                ### set a GFF3 file in MSA file
                file_name_gff3 = project.reference.get_gff3_comulative_positions(
                    TypePath.MEDIA_ROOT
                )
                ### old versions doesn't have this file
                if not os.path.exists(file_name_gff3):
                    software.run_genbank2gff3_positions_comulative(
                        project.reference.get_reference_gbk(TypePath.MEDIA_ROOT),
                        file_name_gff3,
                    )
                if not os.path.exists(file_name_gff3):
                    file_name_gff3 = None
                file_name_gff3 = project.reference.get_gff3_comulative_positions(
                    TypePath.MEDIA_URL
                )

                if element_name == "all_together":
                    file_name_fasta = project.get_global_file_by_project(
                        TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_MAFFT
                    )
                    last_name_seq = None
                    if not file_name_gff3 is None:
                        last_name_seq = utils.get_last_name_from_fasta(file_name_fasta)

                    if os.path.exists(file_name_fasta):
                        file_name_fasta = project.get_global_file_by_project(
                            TypePath.MEDIA_URL, Project.PROJECT_FILE_NAME_MAFFT
                        )
                        data["alignment_fasta_show_id"] = mark_safe(
                            request.build_absolute_uri(file_name_fasta)
                        )
                        if not file_name_gff3 is None:
                            data["gff3_show_id"] = mark_safe(
                                request.build_absolute_uri(file_name_gff3)
                            )
                            data["last_name_seq"] = last_name_seq
                        url_file_name_fasta = '<a href="{}" download> {}</a>'.format(
                            file_name_fasta, os.path.basename(file_name_fasta)
                        )
                    else:
                        url_file_name_fasta = "File not available"
                        data["alignment_fasta_show_id"] = "#"

                    file_name_nex = project.get_global_file_by_project(
                        TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_nex
                    )
                    if os.path.exists(file_name_nex):
                        file_name_nex = project.get_global_file_by_project(
                            TypePath.MEDIA_URL, Project.PROJECT_FILE_NAME_nex
                        )
                        url_file_name_nex = '<a href="{}" download> {}</a>'.format(
                            file_name_nex, os.path.basename(file_name_nex)
                        )
                    else:
                        url_file_name_nex = "File not available"
                else:
                    file_name_fasta = project.get_global_file_by_element(
                        TypePath.MEDIA_ROOT,
                        element_name,
                        Project.PROJECT_FILE_NAME_MAFFT,
                    )
                    if os.path.exists(file_name_fasta):
                        file_name_fasta = project.get_global_file_by_element(
                            TypePath.MEDIA_URL,
                            element_name,
                            Project.PROJECT_FILE_NAME_MAFFT,
                        )
                        data["alignment_fasta_show_id"] = mark_safe(
                            request.build_absolute_uri(file_name_fasta)
                        )
                        url_file_name_fasta = '<a href="{}" download> {}</a>'.format(
                            file_name_fasta, os.path.basename(file_name_fasta)
                        )
                    else:
                        url_file_name_fasta = "File not available"
                        data["alignment_fasta_show_id"] = "#"

                    file_name_nex = project.get_global_file_by_element(
                        TypePath.MEDIA_ROOT, element_name, Project.PROJECT_FILE_NAME_nex
                    )
                    if os.path.exists(file_name_nex):
                        file_name_nex = project.get_global_file_by_element(
                            TypePath.MEDIA_URL,
                            element_name,
                            Project.PROJECT_FILE_NAME_nex,
                        )
                        url_file_name_nex = '<a href="{}" download> {}</a>'.format(
                            file_name_nex, os.path.basename(file_name_nex)
                        )
                    else:
                        url_file_name_nex = "File not available"

                data["is_ok"] = True
                data["alignment_fasta_id"] = mark_safe(
                    "<strong>Alignment (.fasta):</strong> {}".format(
                        url_file_name_fasta
                    )
                )
                data["alignment_nex_id"] = mark_safe(
                    "<strong>Alignment (.nex):</strong> {}".format(url_file_name_nex)
                )
                b_calculate_again = False
                data["max_length_label"] = manage_database.get_max_length_label(
                    project, request.user, b_calculate_again
                )
            except Project.DoesNotExist:
                pass
        return JsonResponse(data)


@csrf_protect
def show_msa_protein(request):
    """
    manage check boxes through ajax
    """
    if request.is_ajax():
        data = {"is_ok": False}
        key_with_project_id = "project_id"
        if key_with_project_id in request.GET:
            project_id = int(request.GET.get(key_with_project_id))
            gene_name = ""
            element_name = ""
            key_element_name = "key_element_name"
            key_gene_name = "key_gene_name"
            if key_element_name in request.GET:
                element_name = request.GET.get(key_element_name)
            if key_gene_name in request.GET:
                gene_name = request.GET.get(key_gene_name)
            if len(element_name) > 0 and len(gene_name) > 0:
                try:
                    manage_database = ManageDatabase()
                    project = Project.objects.get(id=project_id)
                    file_name_fasta = project.get_global_file_by_element_and_cds(
                        TypePath.MEDIA_ROOT,
                        element_name,
                        gene_name,
                        Project.PROJECT_FILE_NAME_MAFFT,
                    )
                    if os.path.exists(file_name_fasta):
                        file_name_fasta = project.get_global_file_by_element_and_cds(
                            TypePath.MEDIA_URL,
                            element_name,
                            gene_name,
                            Project.PROJECT_FILE_NAME_MAFFT,
                        )
                        data["alignment_amino_fasta_show_id"] = mark_safe(
                            request.build_absolute_uri(file_name_fasta)
                        )
                        url_file_name_fasta = '<a href="{}" download> {}</a>'.format(
                            file_name_fasta, os.path.basename(file_name_fasta)
                        )
                    else:
                        url_file_name_fasta = "File not available"
                        data["alignment_amino_fasta_show_id"] = "#"

                    file_name_nex = project.get_global_file_by_element_and_cds(
                        TypePath.MEDIA_ROOT,
                        element_name,
                        gene_name,
                        Project.PROJECT_FILE_NAME_nex,
                    )
                    if os.path.exists(file_name_nex):
                        file_name_nex = project.get_global_file_by_element_and_cds(
                            TypePath.MEDIA_URL,
                            element_name,
                            gene_name,
                            Project.PROJECT_FILE_NAME_nex,
                        )
                        url_file_name_nex = '<a href="{}" download> {}</a>'.format(
                            file_name_nex, os.path.basename(file_name_nex)
                        )
                    else:
                        url_file_name_nex = "File not available"

                    data["is_ok"] = True
                    data["alignment_amino_fasta_id"] = mark_safe(
                        "<strong>Alignment 'fasta':</strong> {}".format(
                            url_file_name_fasta
                        )
                    )
                    data["alignment_amino_nex_id"] = mark_safe(
                        "<strong>Alignment 'nex':</strong> {}".format(url_file_name_nex)
                    )
                    b_calculate_again = False
                    data["max_length_label"] = manage_database.get_max_length_label(
                        project, request.user, b_calculate_again
                    )
                except Project.DoesNotExist:
                    pass
        return JsonResponse(data)


@csrf_protect
def show_count_variations(request):
    """
    get chart information
    """
    if request.is_ajax():
        data = {"is_ok": False}
        key_with_project_id = "project_id"
        if key_with_project_id in request.GET:
            project_id = int(request.GET.get(key_with_project_id))
            try:
                project = Project.objects.get(id=project_id)
                ## need to be order by total
                data["labels"] = []  ## vector with all sample names
                data["data_less_50"] = []  ## number of variations less than 50
                data["data_50_var_90_50"] = []  ##number of variations 50<var<90

                data_out = (
                    []
                )  ## [[#<50, 50<var<90, sample name], [#<50, 50<var<90, sample name], ....]
                for project_sample in project.project_samples.all():
                    if project_sample.is_deleted:
                        continue
                    if project_sample.is_error:
                        continue
                    if not project_sample.is_finished:
                        continue
                    if not project_sample.is_sample_illumina():
                        continue
                    data_out.append(
                        [
                            project_sample.count_variations.var_less_50,
                            project_sample.count_variations.var_bigger_50_90,
                            project_sample.sample.name,
                        ]
                    )

                if len(data_out) > 0:
                    data["is_ok"] = True
                    data_out = sorted(
                        data_out,
                        key=lambda data_temp: data_temp[0] + data_temp[1],
                        reverse=True,
                    )
                    for data_temp in data_out:
                        data["labels"].append(data_temp[2])
                        data["data_less_50"].append(data_temp[0])
                        data["data_50_var_90_50"].append(data_temp[1])
                    data["sample_number"] = len(data_out)
                else:
                    data[
                        "error_message"
                    ] = "There's no samples to collect data to show."
            except Project.DoesNotExist:
                pass
        return JsonResponse(data)


@csrf_protect
def get_cds_from_element(request):
    """
    return the cds's for a specific element, can be more than one
    """
    if request.is_ajax():
        data = {"is_ok": False}
        key_with_project_id = "project_id"
        if key_with_project_id in request.GET:
            project_id = int(request.GET.get(key_with_project_id))
            element_name = ""
            key_element_name = "key_element_name"
            if key_element_name in request.GET:
                element_name = request.GET.get(key_element_name)
            if len(element_name) > 0:
                try:
                    utils = Utils()
                    project = Project.objects.get(id=project_id)
                    vect_genes = utils.get_vect_cds_from_element_from_db(
                        element_name, project.reference, project.owner
                    )
                    if vect_genes != None:
                        data["is_ok"] = True
                        data["vect_genes"] = vect_genes
                except Project.DoesNotExist:
                    pass
        return JsonResponse(data)


@csrf_protect
def get_image_coverage(request):
    """
    get image coverage
    """
    if request.is_ajax():
        data = {"is_ok": False}
        key_with_project_sample_id = "project_sample_id"
        key_element = "element"
        if key_with_project_sample_id in request.GET and key_element in request.GET:
            try:
                project_sample = ProjectSample.objects.get(
                    id=request.GET.get(key_with_project_sample_id)
                )
                path_name = project_sample.get_global_file_by_element(
                    TypePath.MEDIA_URL,
                    ProjectSample.PREFIX_FILE_COVERAGE,
                    request.GET.get(key_element).replace("/", "_"),
                    FileExtensions.FILE_PNG,
                )
                data["is_ok"] = True
                data["image"] = mark_safe(
                    '<img id="coverage_image_id" src="{}" style="width: 100%;">'.format(
                        path_name
                    )
                )
                data["image_download"] = path_name
                data["image_download_name"] = os.path.basename(path_name)
                data["text"] = _(
                    "Coverage for locus '{}'".format(request.GET.get(key_element))
                )
            except ProjectSample.DoesNotExist as e:
                pass
        return JsonResponse(data)


@csrf_protect
def update_project_pangolin(request):
    """
    get image coverage
    """
    if request.is_ajax():
        data = {"is_ok": False}
        key_with_project_id = "project_id"
        if key_with_project_id in request.GET:

            project_id = request.GET[key_with_project_id]
            try:
                project = Project.objects.get(pk=project_id)
            except Project.DoesNotExist:
                return JsonResponse(data)

            try:
                ### need to send a message to recalculate the global files
                metaKeyAndValue = MetaKeyAndValue()
                manageDatabase = ManageDatabase()
                try:
                    process_SGE = ProcessSGE()
                    project = Project.objects.get(id=project_id)
                    taskID = process_SGE.set_collect_update_pangolin_lineage(
                        project, request.user
                    )
                    manageDatabase.set_project_metakey(
                        project,
                        request.user,
                        metaKeyAndValue.get_meta_key(
                            MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id
                        ),
                        MetaKeyAndValue.META_VALUE_Queue,
                        taskID,
                    )
                    data = {"is_ok": True}
                except:
                    data = {"is_ok": False}
            except ProjectSample.DoesNotExist as e:
                pass
        return JsonResponse(data)


@csrf_protect
def show_igv(request):
    """
    get data for IGV
    """
    if request.is_ajax():
        data = {"is_ok": False}
        key_with_project_sample_id = "project_sample_id"
        if key_with_project_sample_id in request.GET:
            try:
                project_sample = ProjectSample.objects.get(
                    id=request.GET.get(key_with_project_sample_id)
                )
                if project_sample.is_sample_illumina():
                    path_name_bam = project_sample.get_file_output(
                        TypePath.MEDIA_URL,
                        FileType.FILE_BAM,
                        SoftwareNames.SOFTWARE_SNIPPY_name,
                    )
                    path_name_bai = project_sample.get_file_output(
                        TypePath.MEDIA_URL,
                        FileType.FILE_BAM_BAI,
                        SoftwareNames.SOFTWARE_SNIPPY_name,
                    )
                    path_name_vcf = project_sample.get_file_output(
                        TypePath.MEDIA_URL,
                        FileType.FILE_VCF,
                        SoftwareNames.SOFTWARE_SNIPPY_name,
                    )
                else:
                    path_name_bam = project_sample.get_file_output(
                        TypePath.MEDIA_URL,
                        FileType.FILE_BAM,
                        SoftwareNames.SOFTWARE_Medaka_name,
                    )
                    path_name_bai = project_sample.get_file_output(
                        TypePath.MEDIA_URL,
                        FileType.FILE_BAM_BAI,
                        SoftwareNames.SOFTWARE_Medaka_name,
                    )
                    path_name_vcf = project_sample.get_file_output(
                        TypePath.MEDIA_URL,
                        FileType.FILE_VCF,
                        SoftwareNames.SOFTWARE_Medaka_name,
                    )

                path_name_reference = (
                    project_sample.project.reference.get_reference_fasta(
                        TypePath.MEDIA_URL
                    )
                )
                path_name_reference_index = (
                    project_sample.project.reference.get_reference_fasta_index(
                        TypePath.MEDIA_URL
                    )
                )

                path_bed = project_sample.project.reference.get_reference_bed(
                    TypePath.MEDIA_URL
                )
                path_bed_idx = project_sample.project.reference.get_reference_bed_index(
                    TypePath.MEDIA_URL
                )
                data["is_ok"] = True
                data["path_bam"] = mark_safe(request.build_absolute_uri(path_name_bam))
                data["path_bed"] = mark_safe(request.build_absolute_uri(path_bed))
                data["path_bed_idx"] = mark_safe(
                    request.build_absolute_uri(path_bed_idx)
                )
                data["path_reference"] = mark_safe(
                    request.build_absolute_uri(path_name_reference)
                )
                data["path_reference_index"] = mark_safe(
                    request.build_absolute_uri(path_name_reference_index)
                )
                data["reference_name"] = project_sample.project.reference.display_name
                data["sample_name"] = project_sample.sample.name

                #### other files
                data["bam_file_id"] = mark_safe(
                    '<strong>Bam file:</strong> <a href="{}" download="{}"> {}</a>'.format(
                        path_name_bam,
                        os.path.basename(path_name_bam),
                        os.path.basename(path_name_bam),
                    )
                )
                data["bai_file_id"] = mark_safe(
                    '<strong>Bai file:</strong> <a href="{}" download="{}"> {}</a>'.format(
                        path_name_bai,
                        os.path.basename(path_name_bai),
                        os.path.basename(path_name_bai),
                    )
                )
                data["vcf_file_id"] = mark_safe(
                    '<strong>Vcf file:</strong> <a href="{}" download="{}"> {}</a>'.format(
                        path_name_vcf,
                        os.path.basename(path_name_vcf),
                        os.path.basename(path_name_vcf),
                    )
                )
                data["reference_id"] = mark_safe(
                    '<strong>Reference:</strong> <a href="{}" download="{}"> {}</a>'.format(
                        path_name_reference,
                        os.path.basename(path_name_reference),
                        os.path.basename(path_name_reference),
                    )
                )
                data["reference_index_id"] = mark_safe(
                    '<strong>Ref. index:</strong> <a href="{}" download="{}"> {}</a>'.format(
                        path_name_reference_index,
                        os.path.basename(path_name_reference_index),
                        os.path.basename(path_name_reference_index),
                    )
                )
            except ProjectSample.DoesNotExist as e:
                pass
        return JsonResponse(data)


@csrf_protect
def validate_project_reference_name(request):
    """
    test if exist this project name
    """
    if request.is_ajax():
        project_name = request.GET.get("project_name")
        request.session[Constants.PROJECT_NAME_SESSION] = project_name

        data = {
            "is_taken": Project.objects.filter(
                name__iexact=project_name,
                is_deleted=False,
                owner__username=request.user.username,
            ).exists()
        }
        if data["is_taken"]:
            data["error_message"] = _("Exists a project with this name.")
        return JsonResponse(data)


@csrf_protect
def validate_reference_name(request):
    """
    test if exist this reference name
    """
    if request.is_ajax():
        reference_name = request.GET.get("reference_name")

        data = {
            "is_taken": Reference.objects.filter(
                name__iexact=reference_name,
                is_deleted=False,
                owner__username=request.user.username,
            ).exists()
        }
        if data["is_taken"]:
            data["error_message"] = _("Exists a reference with this name.")
        return JsonResponse(data)


@csrf_protect
def add_single_value_database(request):
    """
    add a single value to a table in database
    possible tables to add: TagName, DataSet, VaccineStatus
    """
    if request.is_ajax():
        data = {"is_ok": False, "message": "Fail in the system..."}
        key_type_data = "type_data"
        key_value = "value"
        value = request.GET[key_value].strip()

        ### test empty values...
        if len(value) == 0:
            data["message"] = "Empty values are not accepted..."
            return JsonResponse(data)

        ## some pre-requisites
        if not request.user.is_active or not request.user.is_authenticated:
            return JsonResponse(data)
        try:
            profile = Profile.objects.get(user__pk=request.user.pk)
        except Profile.DoesNotExist:
            return JsonResponse(data)
        if profile.only_view_project:
            data["message"] = "'{}' can not add anything.".format(request.user.username)
            return JsonResponse(data)

        if key_value in request.GET and key_type_data in request.GET and len(value) > 0:
            ## add to data_set
            if request.GET[key_type_data] == "id_data_set_add_modal":
                try:
                    DataSet.objects.get(name__iexact=value)
                except DataSet.DoesNotExist as e:
                    dataset = DataSet()
                    dataset.name = value
                    dataset.owner = request.user
                    dataset.save()

                    ### set the data for jscript
                    data["is_ok"] = True
                    data["text"] = dataset.name
                    data["value"] = dataset.id
                    data["message"] = "'{}' is added.".format(value)
            ## add to vaccine status
            elif request.GET[key_type_data] == "id_vaccine_add_modal":
                try:
                    VaccineStatus.objects.get(name__iexact=value)
                except VaccineStatus.DoesNotExist as e:
                    vacine_status = VaccineStatus()
                    vacine_status.name = value
                    vacine_status.owner = request.user
                    vacine_status.save()

                    ### set the data for jscript
                    data["is_ok"] = True
                    data["text"] = vacine_status.name
                    data["value"] = vacine_status.id
                    data["message"] = "'{}' is added.".format(value)
        return JsonResponse(data)


@csrf_protect
def remove_single_value_database(request):
    """
    test if is prossible to remove
    """
    if request.is_ajax():
        data = {"is_ok": False, "message": "Fail in the system..."}
        key_type_data = "type_data"
        key_value = "value"
        key_is_to_test = "is_to_test"

        ## some pre-requisites
        if not request.user.is_active or not request.user.is_authenticated:
            return JsonResponse(data)
        try:
            profile = Profile.objects.get(user__pk=request.user.pk)
        except Profile.DoesNotExist:
            return JsonResponse(data)
        if profile.only_view_project:
            data["message"] = "'{}' can not remove anything.".format(
                request.user.username
            )
            return JsonResponse(data)

        if (
            key_is_to_test in request.GET
            and key_value in request.GET
            and key_type_data in request.GET
            and len(request.GET[key_value].strip()) > 0
        ):
            utils = Utils()
            value = request.GET[key_value].strip()
            is_to_test = utils.str2bool(request.GET[key_is_to_test])

            if len(value) == 0:
                data["message"] = "Empty values are not accepted..."
                return JsonResponse(data)

            ## add to data_set
            if request.GET[key_type_data] == "id_data_set_remove_modal":
                if value == Constants.DATA_SET_GENERIC:
                    data["is_ok"] = False
                    data["is_can_remove"] = False
                    data["message"] = "You can't remove 'Generic' data set"
                else:
                    try:
                        data_set = DataSet.objects.get(name__iexact=value)
                        if data_set.sample.count() > 0:
                            if is_to_test:
                                data["is_ok"] = False
                                data["is_can_remove"] = False
                                data[
                                    "message"
                                ] = "You can't remove '{}' name because has a relation in database.".format(
                                    value
                                )
                            else:
                                data["is_ok"] = False
                                data["is_remove"] = False
                                data[
                                    "message"
                                ] = "You can't remove '{}' name because has a relation in database.".format(
                                    value
                                )
                        else:
                            if is_to_test:
                                data["is_ok"] = True
                                data["is_can_remove"] = True
                            else:
                                data["is_ok"] = True
                                data["is_remove"] = True
                                data["value_to_remove"] = data_set.id
                                data["message"] = "'{}' was removed.".format(value)
                                data_set.delete()
                    except DataSet.DoesNotExist as e:
                        pass

            ## add to vaccine status
            elif request.GET[key_type_data] == "id_vaccine_remove_modal":
                try:
                    vacine_status = VaccineStatus.objects.get(name__iexact=value)
                    if vacine_status.sample.count() > 0:
                        if is_to_test:
                            data["is_ok"] = False
                            data["is_can_remove"] = False
                            data[
                                "message"
                            ] = "You can't remove '{}' name because has a relation in database.".format(
                                value
                            )
                        else:
                            data["is_ok"] = False
                            data["is_remove"] = False
                            data[
                                "message"
                            ] = "You can't remove '{}' name because has a relation in database.".format(
                                value
                            )
                    else:
                        if is_to_test:
                            data["is_ok"] = True
                            data["is_can_remove"] = True
                        else:
                            data["is_ok"] = True
                            data["is_remove"] = False
                            data["message"] = "'{}' was removed.".format(value)
                            data["value_to_remove"] = vacine_status.id
                            vacine_status.delete()
                except VaccineStatus.DoesNotExist as e:
                    data["is_ok"] = False
                    data["is_remove"] = False
                    data["message"] = "You can't remove '{}'.".format(value)
                    pass

        return JsonResponse(data)


@transaction.atomic
@csrf_protect
def remove_reference(request):
    """
    remove a reference. It can only be removed if not belongs to any deleted project
    """
    if request.is_ajax():
        data = {"is_ok": False}
        reference_id_a = "reference_id"

        if reference_id_a in request.GET:

            ## some pre-requisites
            if not request.user.is_active or not request.user.is_authenticated:
                return JsonResponse(data)
            try:
                profile = Profile.objects.get(user__pk=request.user.pk)
            except Profile.DoesNotExist:
                return JsonResponse(data)
            if profile.only_view_project:
                return JsonResponse(data)

            reference_id = request.GET[reference_id_a]
            try:
                reference = Reference.objects.get(pk=reference_id)
            except Profile.DoesNotExist:
                return JsonResponse(data)

            ## different owner or belong to a project not deleted
            if (
                reference.project.all().filter(is_deleted=False).count() != 0
                or reference.owner.pk != request.user.pk
            ):
                return JsonResponse(data)

            ### now you can remove
            reference.is_deleted = True
            reference.is_deleted_in_file_system = False
            reference.date_deleted = datetime.now()
            reference.save()
            data = {"is_ok": True}

        return JsonResponse(data)


@transaction.atomic
@csrf_protect
def remove_sample(request):
    """
    remove a sample, It can only be removed if not belongs to any deleted project
    """
    if request.is_ajax():
        data = {"is_ok": False}
        sample_id_a = "sample_id"
        if sample_id_a in request.GET:

            ## some pre-requisites
            if not request.user.is_active or not request.user.is_authenticated:
                return JsonResponse(data)
            try:
                profile = Profile.objects.get(user__pk=request.user.pk)
            except Profile.DoesNotExist:
                return JsonResponse(data)
            if profile.only_view_project:
                return JsonResponse(data)

            sample_id = request.GET[sample_id_a]
            try:
                sample = Sample.objects.get(pk=sample_id)
            except Sample.DoesNotExist:
                return JsonResponse(data)

            ## different owner or belong to a project not deleted
            if (
                sample.owner.pk != request.user.pk
                or sample.project_samples.all()
                .filter(is_deleted=False, is_error=False, project__is_deleted=False)
                .count()
                != 0
            ):
                return JsonResponse(data)

            ## it can have project samples not deleted but in projects deleted
            for project_samples in sample.project_samples.all().filter(
                is_deleted=False
            ):
                if (
                    not project_samples.is_deleted
                    and not project_samples.project.is_deleted
                ):
                    return JsonResponse(data)

            ### now you can remove
            sample.is_deleted = True
            sample.is_deleted_in_file_system = False
            sample.date_deleted = datetime.now()
            sample.save()

            ## if sample is not processed yet and is waiting for a fastq.gz it is necessary to decrease the value in
            if not sample.has_files:
                for upload_file in sample.uploadfiles_set.all():
                    if (
                        not upload_file.is_processed
                        and upload_file.type_file.name == TypeFile.TYPE_FILE_sample_file
                    ):
                        upload_file.number_files_processed += 1
                        if (
                            upload_file.number_files_to_process
                            == upload_file.number_files_processed
                        ):
                            upload_file.is_processed = True
                        upload_file.save()
                        break

            ## refresh sample list for this user
            process_SGE = ProcessSGE()
            process_SGE.set_create_sample_list_by_user(sample.owner, [])
            data = {"is_ok": True}
        return JsonResponse(data)


@transaction.atomic
@csrf_protect
def remove_project(request):
    """
    remove a project.
    """
    if request.is_ajax():
        data = {"is_ok": False}
        project_id_a = "project_id"

        if project_id_a in request.GET:

            ## some pre-requisites
            if not request.user.is_active or not request.user.is_authenticated:
                return JsonResponse(data)
            try:
                profile = Profile.objects.get(user__pk=request.user.pk)
            except Profile.DoesNotExist:
                return JsonResponse(data)
            if profile.only_view_project:
                return JsonResponse(data)

            project_id = request.GET[project_id_a]
            try:
                project = Project.objects.get(pk=project_id)
            except Project.DoesNotExist:
                return JsonResponse(data)

            ## different owner or belong to a project not deleted
            if project.owner.pk != request.user.pk:
                return JsonResponse(data)

            ### now you can remove
            project.is_deleted = True
            project.is_deleted_in_file_system = False
            project.date_deleted = datetime.now()
            project.save()

            ### delete all project samples
            ### this is only necessary for consistency
            for project_sample in project.project_samples.all():
                project_sample.is_deleted = True
                project_sample.is_deleted_in_file_system = False
                project_sample.date_deleted = datetime.now()
                project_sample.save()

            ## refresh sample and project list for this user
            process_SGE = ProcessSGE()
            process_SGE.set_create_sample_list_by_user(request.user, [])
            process_SGE.set_create_project_list_by_user(request.user)
            data = {"is_ok": True}
        return JsonResponse(data)


@transaction.atomic
@csrf_protect
def remove_televir_project(request):
    """
    remove a project.
    """
    if request.is_ajax():
        data = {"is_ok": False}
        project_id_a = "project_id"

        if project_id_a in request.GET:
            ## some pre-requisites
            if not request.user.is_active or not request.user.is_authenticated:
                return JsonResponse(data)
            try:
                profile = Profile.objects.get(user__pk=request.user.pk)
            except Profile.DoesNotExist:
                return JsonResponse(data)
            if profile.only_view_project:
                return JsonResponse(data)

            project_id = request.GET[project_id_a]
            try:
                project = Televir_Project.objects.get(pk=project_id)
            except Televir_Project.DoesNotExist:
                return JsonResponse(data)

            ## different owner or belong to a project not deleted
            if project.owner.pk != request.user.pk:
                return JsonResponse(data)

            ### now you can remove
            project.is_deleted = True
            project.is_deleted_in_file_system = False
            project.date_deleted = datetime.now()
            project.save()

            ### delete all project samples
            ### this is only necessary for consistency
            for project_sample in PIProject_Sample.objects.filter(project=project):
                project_sample.is_deleted = True
                project_sample.is_deleted_in_file_system = False
                project_sample.date_deleted = datetime.now()
                project_sample.save()

            ## refresh sample and project list for this user
            # process_SGE = ProcessSGE()
            # process_SGE.set_create_sample_list_by_user(request.user, [])
            # process_SGE.set_create_project_list_by_user(request.user)
            data = {"is_ok": True}
        return JsonResponse(data)


@transaction.atomic
@csrf_protect
def remove_televir_project_sample(request):
    """
    remove a project sample.
    """
    if request.is_ajax():
        data = {"is_ok": False}
        project_sample_id_a = "project_sample_id"

        if project_sample_id_a in request.GET:

            ## some pre-requisites
            if not request.user.is_active or not request.user.is_authenticated:
                return JsonResponse(data)
            try:
                profile = Profile.objects.get(user__pk=request.user.pk)
            except Profile.DoesNotExist:
                return JsonResponse(data)
            if profile.only_view_project:
                return JsonResponse(data)

            project_sample_id = request.GET[project_sample_id_a]
            try:
                project_sample = PIProject_Sample.objects.get(pk=project_sample_id)
            except PIProject_Sample.DoesNotExist:
                return JsonResponse(data)

            ## different owner or belong to a project not deleted
            if project_sample.project.owner.pk != request.user.pk:
                return JsonResponse(data)

            ### now you can remove
            project_sample.is_deleted = True
            project_sample.is_deleted_in_file_system = False
            project_sample.date_deleted = datetime.now()
            project_sample.save()
            data = {"is_ok": True}

        return JsonResponse(data)


@transaction.atomic
@csrf_protect
def remove_project_sample(request):
    """
    remove a project sample.
    """
    if request.is_ajax():
        data = {"is_ok": False}
        project_sample_id_a = "project_sample_id"

        if project_sample_id_a in request.GET:

            ## some pre-requisites
            if not request.user.is_active or not request.user.is_authenticated:
                return JsonResponse(data)
            try:
                profile = Profile.objects.get(user__pk=request.user.pk)
            except Profile.DoesNotExist:
                return JsonResponse(data)
            if profile.only_view_project:
                return JsonResponse(data)

            project_sample_id = request.GET[project_sample_id_a]
            try:
                project_sample = ProjectSample.objects.get(pk=project_sample_id)
            except ProjectSample.DoesNotExist:
                return JsonResponse(data)

            ## different owner or belong to a project not deleted
            if project_sample.project.owner.pk != request.user.pk:
                return JsonResponse(data)

            ### now you can remove
            project_sample.is_deleted = True
            project_sample.is_deleted_in_file_system = False
            project_sample.date_deleted = datetime.now()
            project_sample.save()

            ### need to send a message to recalculate the global files
            metaKeyAndValue = MetaKeyAndValue()
            manageDatabase = ManageDatabase()
            try:
                process_SGE = ProcessSGE()
                taskID = process_SGE.set_collect_global_files(
                    project_sample.project, request.user
                )
                manageDatabase.set_project_metakey(
                    project_sample.project,
                    request.user,
                    metaKeyAndValue.get_meta_key(
                        MetaKeyAndValue.META_KEY_Queue_TaskID_Project,
                        project_sample.project.id,
                    ),
                    MetaKeyAndValue.META_VALUE_Queue,
                    taskID,
                )

                ## refresh sample list for this user
                ## project list is updated in collect global files
                process_SGE.set_create_sample_list_by_user(request.user, [])
                data = {"is_ok": True}
            except:
                data = {"is_ok": False}

        return JsonResponse(data)


@transaction.atomic
@csrf_protect
def remove_uploaded_file(request):
    """
    remove a project.
    """
    if request.is_ajax():
        data = {"is_ok": False}
        uploaded_file_id_a = "uploaded_file_id"

        if uploaded_file_id_a in request.GET:

            ## some pre-requisites
            if not request.user.is_active or not request.user.is_authenticated:
                return JsonResponse(data)
            try:
                profile = Profile.objects.get(user__pk=request.user.pk)
            except Profile.DoesNotExist:
                return JsonResponse(data)
            if profile.only_view_project:
                return JsonResponse(data)

            uploaded_file_id = request.GET[uploaded_file_id_a]
            try:
                uploaded_file = UploadFiles.objects.get(pk=uploaded_file_id)
            except UploadFiles.DoesNotExist:
                return JsonResponse(data)

            ## different owner or belong to a project not deleted
            if uploaded_file.owner.pk != request.user.pk:
                return JsonResponse(data)

            ## if it's processed need to test samples deleted
            if uploaded_file.is_processed:
                for sample in uploaded_file.samples.all():
                    if not sample.is_deleted:
                        return JsonResponse(data)

            ### now you can remove
            uploaded_file.is_deleted = True
            uploaded_file.is_deleted_in_file_system = False
            uploaded_file.date_deleted = datetime.now()
            uploaded_file.save()
            data = {"is_ok": True}

        return JsonResponse(data)


@transaction.atomic
@csrf_protect
def remove_uploaded_files(request):
    """
    remove fastq files, all not processed
    """
    if request.is_ajax():
        number_files_removed = 0
        data = {"is_ok": False}
        data["number_files_removed"] = number_files_removed
        data["message_number_files_removed"] = "There's no files removed."

        ## some pre-requisites
        if not request.user.is_active or not request.user.is_authenticated:
            return JsonResponse(data)
        try:
            profile = Profile.objects.get(user__pk=request.user.pk)
        except Profile.DoesNotExist:
            return JsonResponse(data)
        if profile.only_view_project:
            return JsonResponse(data)

        ### get all files that can be deleted, only not processed
        query_set = UploadFiles.objects.filter(
            owner__id=request.user.id,
            is_deleted=False,
            is_processed=False,
            type_file__name=TypeFile.TYPE_FILE_fastq_gz,
        )
        for uploaded_file in query_set:

            ### now you can remove
            uploaded_file.is_deleted = True
            uploaded_file.is_deleted_in_file_system = False
            uploaded_file.date_deleted = datetime.now()
            uploaded_file.save()

            ### new removed file
            number_files_removed += 1
        data = {"is_ok": True}
        data["number_files_removed"] = number_files_removed
        if number_files_removed == 1:
            data["message_number_files_removed"] = "One file removed..."
        elif number_files_removed > 1:
            data["message_number_files_removed"] = "{} files removed...".format(
                number_files_removed
            )

        return JsonResponse(data)


####
def unlock_upload_file(upload_file):
    """
    unlock upload file
    """
    for sample in upload_file.samples.all():
        if not sample.is_ready_for_projects and not sample.is_deleted:
            sample.is_deleted = True
            sample.date_deleted = datetime.now()
            sample.save()
    upload_file.is_processed = True
    upload_file.save()


@transaction.atomic
@csrf_protect
def unlock_sample_file(request):
    """
    unlock sample list files, drop all samples not processed yet
    """
    if request.is_ajax():
        number_of_changes = 0
        data = {"is_ok": False}
        data["number_of_changes"] = number_of_changes
        data["message_number_files_removed"] = "There's no files unlocked."

        ## some pre-requisites
        if not request.user.is_active or not request.user.is_authenticated:
            return JsonResponse(data)
        try:
            profile = Profile.objects.get(user__pk=request.user.pk)
        except Profile.DoesNotExist:
            return JsonResponse(data)
        if profile.only_view_project:
            return JsonResponse(data)

        ### try to find sample_file imported not processed yet
        try:
            metaKey = MetaKey.objects.get(name=TypeFile.TYPE_FILE_sample_file)
        except MetaKey.DoesNotExist:
            return JsonResponse(data)

        number_of_changes = 0
        lst_files = UploadFiles.objects.all().filter(
            is_deleted=False, is_processed=False, owner=request.user, type_file=metaKey
        )
        for uploadfile in lst_files:

            ## teste the number files already processed
            if uploadfile.number_files_processed == uploadfile.number_files_to_process:
                continue

            unlock_upload_file(uploadfile)
            number_of_changes += 1

        data = {"is_ok": True}
        data["number_of_changes"] = number_of_changes
        if number_of_changes == 1:
            data["message_number_of_changes"] = "One sample file was unlocked..."
        elif number_of_changes > 1:
            data[
                "message_number_of_changes"
            ] = "{} sample files were unlocked...".format(number_of_changes)

        return JsonResponse(data)


@csrf_protect
def get_process_running(request):
    """
    get process running and to run for a specific user
    """
    if request.is_ajax():
        data = {"is_ok": False}

        ## some pre-requisites
        if not request.user.is_active or not request.user.is_authenticated:
            return JsonResponse(data)

        ## if it's processed need to test samples deleted
        data["process_running"] = str(
            ProcessControler.objects.filter(
                owner__id=request.user.pk,
                is_finished=False,
                is_error=False,
                is_running=True,
            ).count()
        )
        data["process_to_run"] = str(
            ProcessControler.objects.filter(
                owner__id=request.user.pk,
                is_finished=False,
                is_error=False,
                is_running=False,
            ).count()
        )
        data["process_to_run_total"] = str(
            ProcessControler.objects.filter(is_finished=False, is_error=False).count()
        )
        data["is_ok"] = True
        return JsonResponse(data)


@csrf_protect
def submit_sge(request):
    """
    get process running and to run for a specific user
    """
    if request.is_ajax():
        data = {"is_ok": False}
        process_SGE = ProcessSGE()
        process_SGE.submit_dummy_sge()
        return JsonResponse(data)
