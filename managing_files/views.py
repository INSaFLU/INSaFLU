# Create your views here.

import logging
import ntpath
import os
import sys
from datetime import datetime
from itertools import chain
from operator import attrgetter

from braces.views import FormValidMessageMixin, LoginRequiredMixin
from constants.constants import Constants, FileExtensions, FileType, TypeFile, TypePath
from constants.meta_key_and_values import MetaKeyAndValue
from constants.nextclade_links import get_constext_nextclade
from constants.software_names import SoftwareNames
from django.conf import settings
from django.contrib import messages
from django.contrib.gis.geos import Point
from django.contrib.sites.shortcuts import get_current_site
from django.core.files.temp import NamedTemporaryFile
from django.db import transaction
from django.db.models import Q
from django.http import HttpResponseRedirect, JsonResponse
from django.template.defaultfilters import filesizeformat, pluralize
from django.urls import reverse_lazy
from django.utils.safestring import mark_safe
from django.views import generic
from django.views.generic import DetailView, ListView, TemplateView
from django_tables2 import RequestConfig
from extend_user.models import Profile
from manage_virus.constants_virus import ConstantsVirus
from settings.constants_settings import ConstantsSettings
from settings.default_software_project_sample import DefaultProjectSoftware
from settings.models import Software as SoftwareSettings
from settings.tables import SoftwaresTable
from utils.collect_extra_data import CollectExtraData
from utils.process_SGE import ProcessSGE
from utils.result import DecodeObjects
from utils.session_variables import (
    clean_check_box_in_session,
    is_all_check_box_in_session,
)
from utils.software import Software
from utils.software_pangolin import SoftwarePangolin
from utils.support_django_template import get_link_for_dropdown_item
from utils.utils import ShowInfoMainPage, Utils

from managing_files.forms import (
    AddSampleProjectForm,
    ReferenceForm,
    ReferenceProjectFormSet,
    SampleForm,
    SamplesUploadDescriptionForm,
    SamplesUploadDescriptionMetadataForm,
    SamplesUploadMultipleFastqForm,
)
from managing_files.manage_database import ManageDatabase
from managing_files.models import (
    MetaKey,
    Project,
    ProjectSample,
    Reference,
    Sample,
    UploadFiles,
)
from managing_files.tables import (
    AddSamplesFromCvsFileTable,
    AddSamplesFromCvsFileTableMetadata,
    AddSamplesFromFastqFileTable,
    ProjectTable,
    ReferenceProjectTable,
    ReferenceTable,
    SampleTable,
    SampleToProjectsTable,
    ShowProjectSamplesResults,
)

# http://www.craigderington.me/generic-list-view-with-django-tables/


class ProjectIndex(TemplateView):
    template_name = "project/project_index.html"


# class ReferencesView(LoginRequiredMixin, GroupRequiredMixin, ListView):
class ReferenceView(LoginRequiredMixin, ListView):
    model = Reference
    template_name = "references/references.html"
    context_object_name = "reference"
    ordering = ["id"]

    def get_context_data(self, **kwargs):
        context = super(ReferenceView, self).get_context_data(**kwargs)

        tag_search = "search_references"
        query_set = Reference.objects.filter(
            owner__id=self.request.user.id, is_obsolete=False, is_deleted=False
        ).order_by("-name")
        if self.request.GET.get(tag_search) != None and self.request.GET.get(
            tag_search
        ):
            query_set = query_set.filter(
                Q(name__icontains=self.request.GET.get(tag_search))
                | Q(owner__username__icontains=self.request.GET.get(tag_search))
                | Q(reference_genbank_name__icontains=self.request.GET.get(tag_search))
                | Q(reference_fasta_name__icontains=self.request.GET.get(tag_search))
                | Q(isolate_name__icontains=self.request.GET.get(tag_search))
            )

        ### get the references from the system
        query_set_system = Reference.objects.filter(
            owner__username=Constants.DEFAULT_USER, is_obsolete=False, is_deleted=False
        ).order_by("-name")
        if self.request.GET.get(tag_search) != None and self.request.GET.get(
            tag_search
        ):
            query_set_system = query_set_system.filter(
                Q(name__icontains=self.request.GET.get(tag_search))
                | Q(owner__username__icontains=self.request.GET.get(tag_search))
                | Q(reference_genbank_name__icontains=self.request.GET.get(tag_search))
                | Q(reference_fasta_name__icontains=self.request.GET.get(tag_search))
                | Q(isolate_name__icontains=self.request.GET.get(tag_search))
            )
        query_set_result = sorted(
            chain(query_set, query_set_system), key=attrgetter("creation_date")
        )
        table = ReferenceTable(query_set_result)
        RequestConfig(
            self.request, paginate={"per_page": Constants.PAGINATE_NUMBER}
        ).configure(table)
        if self.request.GET.get(tag_search) != None:
            context[tag_search] = self.request.GET.get(tag_search)
        context["table"] = table
        context["nav_reference"] = True
        context["show_paginatior"] = len(query_set_result) > Constants.PAGINATE_NUMBER
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        return context


class ReferenceAddView(LoginRequiredMixin, FormValidMessageMixin, generic.FormView):
    """
    Create a new reference
    """

    form_class = ReferenceForm
    success_url = reverse_lazy("references")
    template_name = "references/reference_add.html"

    ## Other solution to get the reference
    ## https://pypi.python.org/pypi?%3aaction=display&name=django-contrib-requestprovider&version=1.0.1
    def get_form_kwargs(self):
        """
        Set the request to pass in the form
        """
        kw = super(ReferenceAddView, self).get_form_kwargs()
        kw["request"] = self.request  # the trick!
        return kw

    def get_context_data(self, **kwargs):
        context = super(ReferenceAddView, self).get_context_data(**kwargs)
        context["nav_reference"] = True
        context["nav_modal"] = True  ## short the size of modal window
        context["user_mmp"] = (
            self.request.user.username == "mmp"
        )  ## short the size of modal window
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        return context

    @transaction.atomic
    def form_valid(self, form):
        from utils.software import Software

        ### test anonymous account
        try:
            profile = Profile.objects.get(user=self.request.user)
            if profile.only_view_project:
                messages.warning(
                    self.request,
                    "'{}' account can not add references.".format(
                        self.request.user.username
                    ),
                    fail_silently=True,
                )
                return super(ReferenceAddView, self).form_invalid(form)
        except Profile.DoesNotExist:
            pass

        software = Software()
        utils = Utils()
        name = form.cleaned_data["name"]
        scentific_name = form.cleaned_data["isolate_name"]
        reference_fasta = form.cleaned_data["reference_fasta"]
        reference_genbank = form.cleaned_data["reference_genbank"]
        if scentific_name is None:
            scentific_name = ""

        ## create a genbank file
        temp_genbank_dir = None
        temp_genbank_file = None
        if reference_genbank is None:
            reference_fasta_temp_file_name = NamedTemporaryFile(
                prefix="flu_fa_", delete=False
            )
            reference_fasta.file.seek(0)
            reference_fasta_temp_file_name.write(reference_fasta.read())
            reference_fasta_temp_file_name.flush()
            reference_fasta_temp_file_name.close()
            software.dos_2_unix(reference_fasta_temp_file_name.name)
            software.fasta_2_upper(reference_fasta_temp_file_name.name)

            file_name_cleaned = utils.clean_name(ntpath.basename(reference_fasta.name))
            try:
                temp_genbank_dir = software.run_prokka(
                    reference_fasta_temp_file_name.name, file_name_cleaned
                )
            except Exception as e:
                os.unlink(reference_fasta_temp_file_name.name)
                messages.error(
                    self.request, "Error creating the genbank file", fail_silently=True
                )
                return super(ReferenceAddView, self).form_invalid(form)
            os.unlink(reference_fasta_temp_file_name.name)

            temp_genbank_file = os.path.join(
                temp_genbank_dir,
                utils.clean_extension(file_name_cleaned) + FileExtensions.FILE_GBK,
            )
            if not os.path.exists(temp_genbank_file):
                utils.remove_dir(temp_genbank_dir)
                messages.error(
                    self.request, "Error creating the genbank file", fail_silently=True
                )
                return super(ReferenceAddView, self).form_invalid(form)

        reference = form.save(commit=False)
        ## set other data
        reference.display_name = reference.name
        reference.owner = self.request.user
        reference.is_obsolete = False
        reference.number_of_locus = self.request.session[
            Constants.NUMBER_LOCUS_FASTA_FILE
        ]
        reference.reference_fasta_name = utils.clean_name(
            ntpath.basename(reference_fasta.name)
        )
        reference.scentific_name = scentific_name
        if reference_genbank == None:
            reference.reference_genbank_name = ntpath.basename(temp_genbank_file)
        else:
            reference.reference_genbank_name = utils.clean_name(
                ntpath.basename(reference_genbank.name)
            )
        reference.save()

        ## move the files to the right place
        sz_file_to = os.path.join(
            settings.MEDIA_ROOT,
            utils.get_path_to_reference_file(self.request.user.id, reference.id),
            reference.reference_fasta_name,
        )
        software.dos_2_unix(
            os.path.join(settings.MEDIA_ROOT, reference.reference_fasta.name)
        )
        ## test if bases all lower
        software.fasta_2_upper(
            os.path.join(settings.MEDIA_ROOT, reference.reference_fasta.name)
        )
        utils.move_file(
            os.path.join(settings.MEDIA_ROOT, reference.reference_fasta.name),
            sz_file_to,
        )
        reference.hash_reference_fasta = utils.md5sum(sz_file_to)
        reference.reference_fasta.name = os.path.join(
            utils.get_path_to_reference_file(self.request.user.id, reference.id),
            reference.reference_fasta_name,
        )

        ###
        ## genbank file
        sz_file_to = os.path.join(
            settings.MEDIA_ROOT,
            utils.get_path_to_reference_file(self.request.user.id, reference.id),
            reference.reference_genbank_name,
        )
        if reference_genbank is None:
            if not temp_genbank_file is None:
                utils.move_file(temp_genbank_file, sz_file_to)
        else:
            utils.move_file(
                os.path.join(settings.MEDIA_ROOT, reference.reference_genbank.name),
                sz_file_to,
            )
        software.dos_2_unix(sz_file_to)
        reference.hash_reference_genbank = utils.md5sum(sz_file_to)
        reference.reference_genbank.name = os.path.join(
            utils.get_path_to_reference_file(self.request.user.id, reference.id),
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
        utils.get_elements_from_db(reference, self.request.user)
        utils.get_elements_and_cds_from_db(reference, self.request.user)

        ## create the index before commit in database, throw exception if something goes wrong
        software.create_fai_fasta(
            os.path.join(
                getattr(settings, "MEDIA_ROOT", None), reference.reference_fasta.name
            )
        )

        ### set specie tag
        software.get_species_tag(reference)

        ## remove genbank temp dir if exist
        if temp_genbank_dir != None:
            utils.remove_dir(temp_genbank_dir)

        messages.success(
            self.request,
            "Reference '" + name + "' was created successfully",
            fail_silently=True,
        )
        return super(ReferenceAddView, self).form_valid(form)

    ## static method, not need for now.
    form_valid_message = ""  ## need to have this


class SamplesView(LoginRequiredMixin, ListView):
    model = Sample
    utils = Utils()
    template_name = "samples/samples.html"
    context_object_name = "samples"

    def get_context_data(self, **kwargs):
        context = super(SamplesView, self).get_context_data(**kwargs)
        search_key = "search_samples"
        tag_search = self.request.GET.get(search_key)
        query_set = Sample.objects.filter(
            owner__id=self.request.user.id, is_deleted=False
        ).order_by("-creation_date")
        if not tag_search is None and len(tag_search) > 0:
            filter_dict = {
                "name__icontains": tag_search,
                "type_subtype__icontains": tag_search,
                "data_set__name__icontains": tag_search,
            }  # Dict with fields
            or_condition = Q()
            for key, value in filter_dict.items():
                or_condition.add(Q(**{key: value}), Q.OR)

            if tag_search.lower().startswith(
                ConstantsSettings.TECHNOLOGY_illumina.lower()
            ):
                or_condition.add(Q(type_of_fastq=Sample.TYPE_OF_FASTQ_illumina), Q.OR)
            elif tag_search.lower().startswith(
                ConstantsSettings.TECHNOLOGY_minion.lower()
            ):
                or_condition.add(Q(type_of_fastq=Sample.TYPE_OF_FASTQ_minion), Q.OR)
            elif tag_search.lower().startswith(
                ConstantsSettings.TECHNOLOGY_Undefined.lower()
            ):
                or_condition.add(
                    Q(type_of_fastq=Sample.TYPE_OF_FASTQ_not_defined), Q.OR
                )

            ### filtering
            query_set = query_set.filter(or_condition)

        table = SampleTable(query_set)
        RequestConfig(
            self.request, paginate={"per_page": Constants.PAGINATE_NUMBER}
        ).configure(table)
        if not tag_search is None:
            context[search_key] = tag_search

        ## list of all samples in CSV and TSV
        csv_file = self.utils.get_sample_list_by_user(
            self.request.user.id, "MEDIA_ROOT", FileExtensions.FILE_CSV
        )
        if os.path.exists(csv_file):
            context["list_samples_file_csv"] = get_link_for_dropdown_item(
                self.utils.get_sample_list_by_user(
                    self.request.user.id, "MEDIA_URL", FileExtensions.FILE_CSV
                )
            )
        tsv_file = self.utils.get_sample_list_by_user(
            self.request.user.id, "MEDIA_ROOT", FileExtensions.FILE_TSV
        )
        if os.path.exists(tsv_file):
            context["list_samples_file_tsv"] = get_link_for_dropdown_item(
                self.utils.get_sample_list_by_user(
                    self.request.user.id, "MEDIA_URL", FileExtensions.FILE_TSV
                )
            )

        context["table"] = table
        context["nav_sample"] = True
        context["total_itens"] = query_set.count()
        context["show_paginatior"] = query_set.count() > Constants.PAGINATE_NUMBER
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        return context


class SamplesAddView(LoginRequiredMixin, FormValidMessageMixin, generic.FormView):
	"""
	Create a new reference
	"""
	form_class = SampleForm
	success_url = reverse_lazy('samples')
	template_name = 'samples/sample_add.html'

	if settings.DEBUG: logger = logging.getLogger("fluWebVirus.debug")
	else: logger = logging.getLogger("fluWebVirus.production")

	def get_form_kwargs(self):
		"""
		Set the request to pass in the form
		"""
		kw = super(SamplesAddView, self).get_form_kwargs()
		kw['request'] = self.request 	# the trick!
		return kw


	def get_context_data(self, **kwargs):
		context = super(SamplesAddView, self).get_context_data(**kwargs)
		context['nav_sample'] = True
		context['nav_modal'] = True	## short the size of modal window
		context['show_info_main_page'] = ShowInfoMainPage()		## show main information about the institute
		return context


	def form_valid(self, form):
		"""
		Validate the form
		"""

		with transaction.atomic():
			### test anonymous account
			try:
				profile = Profile.objects.get(user=self.request.user)
				if (profile.only_view_project):
					messages.warning(self.request, "'{}' account can not add samples.".format(self.request.user.username), fail_silently=True)
					return super(SamplesAddView, self).form_invalid(form)
			except Profile.DoesNotExist:
				pass
	
			utils = Utils()
			name = form.cleaned_data['name']
			lat = form.cleaned_data['lat']
			lng = form.cleaned_data['lng']
			like_dates = form.cleaned_data['like_dates']
				
			sample = form.save(commit=False)
			## set other data
			sample.owner = self.request.user
			sample.is_deleted = False
			sample.is_obsolete = False
			sample.file_name_1 = utils.clean_name(os.path.basename(sample.path_name_1.name))
			sample.is_valid_1 = True
			if (sample.exist_file_2()):
				sample.file_name_2 = utils.clean_name(os.path.basename(sample.path_name_2.name))
				sample.is_valid_2 = True 
			else: sample.is_valid_2 = False
			sample.has_files = True
			
			### set type of sequencing, illumina, minion...
			## here, the file is already tested for gastq.gz and illumina and minion
			sample.set_type_of_fastq_sequencing(form.cleaned_data['type_fastq'])
			
			if (like_dates == 'date_of_onset'):
				sample.day = int(sample.date_of_onset.strftime("%d"))
				sample.week = int(sample.date_of_onset.strftime("%W")) + 1
				sample.year = int(sample.date_of_onset.strftime("%Y"))
				sample.month = int(sample.date_of_onset.strftime("%m"))
			elif (like_dates == 'date_of_collection'):
				sample.day = int(sample.date_of_collection.strftime("%d"))
				sample.week = int(sample.date_of_collection.strftime("%W")) + 1
				sample.year = int(sample.date_of_collection.strftime("%Y"))
				sample.month = int(sample.date_of_collection.strftime("%m"))
			elif (like_dates == 'date_of_receipt_lab'):
				sample.day = int(sample.date_of_receipt_lab.strftime("%d"))
				sample.week = int(sample.date_of_receipt_lab.strftime("%W")) + 1
				sample.year = int(sample.date_of_receipt_lab.strftime("%Y"))
				sample.month = int(sample.date_of_receipt_lab.strftime("%m"))
			
			### test geo spacing
			if (lat != None and lng != None): sample.geo_local = Point(lat, lng)
			sample.save()
	
			## move the files to the right place
			sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), utils.get_path_to_fastq_file(self.request.user.id, sample.id), sample.file_name_1)
			utils.move_file(os.path.join(getattr(settings, "MEDIA_ROOT", None), sample.path_name_1.name), sz_file_to)
			sample.path_name_1.name = os.path.join(utils.get_path_to_fastq_file(self.request.user.id, sample.id), sample.file_name_1)
			
			if (sample.exist_file_2()):
				sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), utils.get_path_to_fastq_file(self.request.user.id, sample.id), sample.file_name_2)
				utils.move_file(os.path.join(getattr(settings, "MEDIA_ROOT", None), sample.path_name_2.name), sz_file_to)
				sample.path_name_2.name = os.path.join(utils.get_path_to_fastq_file(self.request.user.id, sample.id), sample.file_name_2)
			sample.save()

		### create a task to perform the analysis of fastq and trimmomatic
		try:
			process_SGE = ProcessSGE()
			(job_name_wait, job_name) = self.request.user.profile.get_name_sge_seq(Profile.SGE_PROCESS_clean_sample, Profile.SGE_SAMPLE)
			if sample.is_type_fastq_gz_sequencing():	### default is Illumina
				taskID = process_SGE.set_run_trimmomatic_species(sample, self.request.user, job_name)
			else:										### Minion, codify with other
				taskID = process_SGE.set_run_clean_minion(sample, self.request.user, job_name)
		except Exception as e:
			self.logger.error('Fail to run: ProcessSGE - ' + str(e))
			return super(SamplesAddView, self).form_invalid(form)
		
		## refresh sample list for this user
		if not job_name is None:
			process_SGE.set_create_sample_list_by_user(self.request.user, [job_name])
		### 
		manageDatabase = ManageDatabase()
		manageDatabase.set_sample_metakey(sample, self.request.user, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue, taskID)
		
		messages.success(self.request, "Sample '" + name + "' was created successfully", fail_silently=True)
		return super(SamplesAddView, self).form_valid(form)

	form_valid_message = ""		## need to have this, even empty


class SamplesAddDescriptionFileView(LoginRequiredMixin, FormValidMessageMixin, generic.CreateView):
	"""
	Create a new reference
	"""
	utils = Utils()
	success_url = reverse_lazy('samples')
	template_name = 'samples/sample_description_file.html'
	model = UploadFiles
	fields = ['file_name']
	
	def get_form_kwargs(self):
		"""
		Set the request to pass in the form
		"""
		kw = super(SamplesAddDescriptionFileView, self).get_form_kwargs()
		return kw
	
	def get_context_data(self, **kwargs):
		context = super(SamplesAddDescriptionFileView, self).get_context_data(**kwargs)
		
		### test anonymous account
		disable_upload_files = False
		try:
			profile = Profile.objects.get(user=self.request.user)
			if (profile.only_view_project): disable_upload_files = True
		except Profile.DoesNotExist:
			disable_upload_files = True
		
		tag_search = 'search_samples'
		query_set = UploadFiles.objects.filter(owner__id=self.request.user.id, is_deleted=False,\
				type_file__name=TypeFile.TYPE_FILE_sample_file).order_by('-creation_date')
		if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)): 
			query_set = query_set.filter(Q(file_name__icontains=self.request.GET.get(tag_search)) |\
										Q(owner__username__icontains=self.request.GET.get(tag_search)))
		table = AddSamplesFromCvsFileTable(query_set)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
		if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
		context['table'] = table
		context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
		context['nav_sample'] = True
		context['disable_upload_files'] = disable_upload_files
		
		### test if exists files to process to match with (csv/tsv) file
		context['does_not_exists_fastq_files_to_process'] = UploadFiles.objects.filter(owner__id=self.request.user.id, is_deleted=False,\
				type_file__name=TypeFile.TYPE_FILE_sample_file).order_by('-creation_date').count() == 0
				
		### test if can add other csv file
		count_not_complete = UploadFiles.objects.filter(owner__id=self.request.user.id, is_deleted=False,\
				type_file__name=TypeFile.TYPE_FILE_sample_file, is_processed=False).count()
		if (count_not_complete > 0): 
			context['can_add_other_file'] = "You cannot add a new file because you must first upload NGS data regarding another file."
			context['disable_upload_files'] = True
			context['message_unlock_file'] = "It will drop remain samples ({}) not processed in last file".format(count_not_complete)
			context['message_unlock_file_question'] = "Do you want drop remain samples ({}) not processed in last file?".format(count_not_complete)
		else:
			context['message_unlock_file'] = "No samples to drop for last sample file."
			context['message_unlock_file_question'] = "No samples to drop for last sample file."
		
		context['show_info_main_page'] = ShowInfoMainPage()		## show main information about the institute 
		return context

	def form_valid(self, form):
		"""
		Validate the form
		"""
		### test anonymous account
		try:
			profile = Profile.objects.get(user=self.request.user)
			if (profile.only_view_project):
				messages.warning(self.request, "'{}' account can not add description files.".format(self.request.user.username), fail_silently=True)
				return super(SamplesAddDescriptionFileView, self).form_invalid(form)
		except Profile.DoesNotExist:
			pass
		
		return super(SamplesAddDescriptionFileView, self).form_valid(form)

	form_valid_message = ""		## need to have this, even empty



class SamplesUpdateMetadata(LoginRequiredMixin, FormValidMessageMixin, generic.CreateView):
	"""
	Update metadata
	"""
	utils = Utils()
	success_url = reverse_lazy('samples')
	template_name = 'samples/sample_update_metadata.html'
	model = UploadFiles
	fields = ['file_name']
	
	def get_form_kwargs(self):
		"""
		Set the request to pass in the form
		"""
		kw = super(SamplesUpdateMetadata, self).get_form_kwargs()
		return kw
	
	def get_context_data(self, **kwargs):
		context = super(SamplesUpdateMetadata, self).get_context_data(**kwargs)
		
		### test anonymous account
		disable_upload_files = False
		try:
			profile = Profile.objects.get(user=self.request.user)
			if (profile.only_view_project): disable_upload_files = True
		except Profile.DoesNotExist:
			pass
		
		tag_search = 'search_samples'
		query_set = UploadFiles.objects.filter(owner__id=self.request.user.id, is_deleted=False,\
				type_file__name=TypeFile.TYPE_FILE_sample_file_metadata).order_by('-creation_date')
		if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)): 
			query_set = query_set.filter(Q(file_name__icontains=self.request.GET.get(tag_search)) |\
										Q(owner__username__icontains=self.request.GET.get(tag_search)))
		table = AddSamplesFromCvsFileTableMetadata(query_set)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
		if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
		context['table'] = table
		context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
		context['nav_sample'] = True
		context['disable_upload_files'] = disable_upload_files
		
		### test if exists files to process to match with (csv/tsv) file
		context['does_not_exists_fastq_files_to_process'] = UploadFiles.objects.filter(owner__id=self.request.user.id, is_deleted=False,\
				type_file__name=TypeFile.TYPE_FILE_sample_file_metadata).order_by('-creation_date').count() == 0
				
		### test if can add other csv file
		count_not_complete = UploadFiles.objects.filter(owner__id=self.request.user.id, is_deleted=False,\
				type_file__name=TypeFile.TYPE_FILE_sample_file_metadata, is_processed=False).count()
		if (count_not_complete > 0): 
			context['can_add_other_file'] = "You cannot add other file because there is a file in pipeline."
			context['disable_upload_files'] = True
			
		context['show_info_main_page'] = ShowInfoMainPage()		## show main information about the institute 
		return context

	def form_valid(self, form):
		"""
		Validate the form
		"""
		### test anonymous account
		try:
			profile = Profile.objects.get(user=self.request.user)
			if (profile.only_view_project):
				messages.warning(self.request, "'{}' account can not add metadata.".format(self.request.user.username), fail_silently=True)
				return super(SamplesUpdateMetadata, self).form_invalid(form)
		except Profile.DoesNotExist:
			pass
		
		return super(SamplesUpdateMetadata, self).form_valid(form)

	form_valid_message = ""		## need to have this, even empty


class SamplesUploadDescriptionFileView(LoginRequiredMixin, FormValidMessageMixin, generic.FormView):
	"""
	Set new samples
	"""
	form_class = SamplesUploadDescriptionForm
	success_url = reverse_lazy('sample-add-file')
	template_name = 'samples/samples_upload_description_file.html'

	def get_form_kwargs(self):
		"""
		Set the request to pass in the form
		"""
		kw = super(SamplesUploadDescriptionFileView, self).get_form_kwargs()
		kw['request'] = self.request 	# get error
		return kw
	
	def get_context_data(self, **kwargs):
		context = super(SamplesUploadDescriptionFileView, self).get_context_data(**kwargs)
		if ('form' in kwargs and hasattr(kwargs['form'], 'error_in_file')):
			context['error_in_file'] = mark_safe(kwargs['form'].error_in_file.replace('\n', "<br>")) ## pass a list
		context['nav_sample'] = True
		context['nav_modal'] = True	## short the size of modal window
		context['show_info_main_page'] = ShowInfoMainPage()		## show main information about the institute
		return context

	def form_valid(self, form):
		
		### test anonymous account
		try:
			profile = Profile.objects.get(user=self.request.user)
			if (profile.only_view_project):
				messages.warning(self.request, "'{}' account can not add file with samples.".format(self.request.user.username), fail_silently=True)
				return super(SamplesUploadDescriptionFileView, self).form_invalid(form)
		except Profile.DoesNotExist:
			pass

		utils = Utils()
		software = Software()
		path_name = form.cleaned_data['path_name']

		## create a genbank file
		if (not path_name is None):
			upload_files = form.save(commit=False)
			upload_files.is_valid = True
			upload_files.is_processed = False
			upload_files.is_deleted = False
			upload_files.number_errors = 0
			upload_files.number_files_processed = 0
			upload_files.number_files_to_process = form.number_files_to_process
			
			try:
				type_file = MetaKey.objects.get(name=TypeFile.TYPE_FILE_sample_file)
			except MetaKey.DoesNotExist:
				type_file = MetaKey()
				type_file.name = TypeFile.TYPE_FILE_sample_file
				type_file.save()
			
			upload_files.type_file = type_file
			upload_files.file_name = utils.clean_name(ntpath.basename(path_name.name))
			upload_files.owner = self.request.user
			
			upload_files.description = ""
			upload_files.save()
		
			## move the files to the right place
			sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), utils.get_path_upload_file(self.request.user.id,\
													TypeFile.TYPE_FILE_sample_file), upload_files.file_name)
			sz_file_to, path_added = utils.get_unique_file(sz_file_to)		## get unique file name, user can upload files with same name...
			utils.move_file(os.path.join(getattr(settings, "MEDIA_ROOT", None), upload_files.path_name.name), sz_file_to)
			software.dos_2_unix(sz_file_to)
			if path_added is None:
				upload_files.path_name.name = os.path.join(utils.get_path_upload_file(self.request.user.id,\
									TypeFile.TYPE_FILE_sample_file), ntpath.basename(sz_file_to))
			else:
				upload_files.path_name.name = os.path.join(utils.get_path_upload_file(self.request.user.id,\
									TypeFile.TYPE_FILE_sample_file), path_added, ntpath.basename(sz_file_to))
			upload_files.save()
			
			try:
				process_SGE = ProcessSGE()
				taskID =  process_SGE.set_read_sample_file(upload_files, self.request.user)
			except:
				return super(SamplesUploadDescriptionFileView, self).form_invalid(form)
			
			messages.success(self.request, "File '" + upload_files.file_name + "' with samples was uploaded successfully", fail_silently=True)
			return super(SamplesUploadDescriptionFileView, self).form_valid(form)
		return super(SamplesUploadDescriptionFileView, self).form_invalid(form)

	## static method, not need for now.
	form_valid_message = ""		## need to have this


class SamplesUploadDescriptionFileViewMetadata(LoginRequiredMixin, FormValidMessageMixin, generic.FormView):
	"""
	Create a new reference
	"""
	form_class = SamplesUploadDescriptionMetadataForm
	success_url = reverse_lazy('sample-update-metadata')
	template_name = 'samples/samples_upload_description_file_metadata.html'

	def get_form_kwargs(self):
		"""
		Set the request to pass in the form
		"""
		kw = super(SamplesUploadDescriptionFileViewMetadata, self).get_form_kwargs()
		kw['request'] = self.request 	# get error
		return kw
	
	def get_context_data(self, **kwargs):
		context = super(SamplesUploadDescriptionFileViewMetadata, self).get_context_data(**kwargs)
		if ('form' in kwargs and hasattr(kwargs['form'], 'error_in_file')):
			context['error_in_file'] = mark_safe(kwargs['form'].error_in_file.replace('\n', "<br>")) ## pass a list
		context['nav_sample'] = True
		context['nav_modal'] = True	## short the size of modal window
		context['show_info_main_page'] = ShowInfoMainPage()		## show main information about the institute
		return context

	def form_valid(self, form):
		
		### test anonymous account
		try:
			profile = Profile.objects.get(user=self.request.user)
			if (profile.only_view_project):
				messages.warning(self.request, "'{}' account can not add file with samples.".format(self.request.user.username), fail_silently=True)
				return super(SamplesUploadDescriptionFileViewMetadata, self).form_invalid(form)
		except Profile.DoesNotExist:
			pass

		utils = Utils()
		software = Software()
		path_name = form.cleaned_data['path_name']

		## create a genbank file
		if (not path_name is None):
			upload_files = form.save(commit=False)
			upload_files.is_valid = True
			upload_files.is_processed = False
			upload_files.is_deleted = False
			upload_files.number_errors = 0
			upload_files.number_files_processed = 0
			upload_files.number_files_to_process = form.number_files_to_process
			
			try:
				type_file = MetaKey.objects.get(name=TypeFile.TYPE_FILE_sample_file_metadata)
			except MetaKey.DoesNotExist:
				type_file = MetaKey()
				type_file.name = TypeFile.TYPE_FILE_sample_file_metadata
				type_file.save()
			
			upload_files.type_file = type_file
			upload_files.file_name = utils.clean_name(ntpath.basename(path_name.name))
			upload_files.owner = self.request.user
			
			upload_files.description = ""
			upload_files.save()				## need this save because of 
		
			## move the files to the right place
			sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), utils.get_path_upload_file(self.request.user.id,\
													TypeFile.TYPE_FILE_sample_file_metadata), upload_files.file_name)
			sz_file_to, path_added = utils.get_unique_file(sz_file_to)		## get unique file name, user can upload files with same name...
			utils.move_file(os.path.join(getattr(settings, "MEDIA_ROOT", None), upload_files.path_name.name), sz_file_to)
			software.dos_2_unix(sz_file_to)
			if path_added is None:
				upload_files.path_name.name = os.path.join(utils.get_path_upload_file(self.request.user.id,\
									TypeFile.TYPE_FILE_sample_file_metadata), ntpath.basename(sz_file_to))
			else:
				upload_files.path_name.name = os.path.join(utils.get_path_upload_file(self.request.user.id,\
									TypeFile.TYPE_FILE_sample_file_metadata), path_added, ntpath.basename(sz_file_to))
			upload_files.save()
			
			try:
				process_SGE = ProcessSGE()
				taskID =  process_SGE.set_read_sample_file_with_metadata(upload_files, self.request.user)
			except:
				return super(SamplesUploadDescriptionFileViewMetadata, self).form_invalid(form)
			
			messages.success(self.request, "File '" + upload_files.file_name + "' with metadata was uploaded successfully", fail_silently=True)
			return super(SamplesUploadDescriptionFileViewMetadata, self).form_valid(form)
		return super(SamplesUploadDescriptionFileViewMetadata, self).form_invalid(form)

	## static method, not need for now.
	form_valid_message = ""		## need to have this


class SamplesAddFastQView(LoginRequiredMixin, FormValidMessageMixin, generic.FormView):
	"""
	Add fastq files to system
	"""
	form_class = SampleForm
	success_url = reverse_lazy('samples')
	template_name = 'samples/sample_fastq_file.html'


	def get_form_kwargs(self):
		"""
		Set the request to pass in the form
		"""
		kw = super(SamplesAddFastQView, self).get_form_kwargs()
		kw['request'] = self.request 	# the trick!
		return kw


	def get_context_data(self, **kwargs):
		context = super(SamplesAddFastQView, self).get_context_data(**kwargs)
		
		### test anonymous account
		disable_upload_files = False
		try:
			profile = Profile.objects.get(user=self.request.user)
			if (profile.only_view_project): disable_upload_files = True
		except Profile.DoesNotExist:
			pass
		
		### test to show only the processed or all
		if ('show-not-only-checked' in self.request.GET):
			b_show_all = self.request.GET.get('show-not-only-checked') != 'on'
		else: b_show_all = True
		
		### get number of files that can be removed
		number_files_can_be_removed = UploadFiles.objects.filter(owner__id=self.request.user.id, is_deleted=False,\
				is_processed=False, type_file__name=TypeFile.TYPE_FILE_fastq_gz).count()

		### 
		tag_search = 'search_samples'
		query_set = UploadFiles.objects.filter(owner__id=self.request.user.id, is_deleted=False,\
				type_file__name=TypeFile.TYPE_FILE_fastq_gz).order_by('-creation_date')
		if (self.request.GET.get(tag_search) != None and self.request.GET.get(tag_search)):
			query_set = query_set.filter(Q(file_name__icontains=self.request.GET.get(tag_search)) |\
							Q(owner__username__icontains=self.request.GET.get(tag_search)))
		if (not b_show_all):
			query_set = query_set.filter(Q(is_processed=b_show_all))
			
		table = AddSamplesFromFastqFileTable(query_set)
		RequestConfig(self.request, paginate={'per_page': Constants.PAGINATE_NUMBER}).configure(table)
		if (self.request.GET.get(tag_search) != None): context[tag_search] = self.request.GET.get(tag_search)
		context['table'] = table
		context['show_paginatior'] = query_set.count() > Constants.PAGINATE_NUMBER
		context['total_itens'] = query_set.count() 
		context['nav_sample'] = True
		context['disable_upload_files'] = disable_upload_files
		context['check_box_not_show_processed_files'] = not b_show_all

		### number of files to remove
		context['disable_remove_all_files'] = number_files_can_be_removed == 0		
		if (number_files_can_be_removed == 0):
			context['message_remove_files'] = "There's no files to remove"
			context['message_remove_files_2'] = "There's no files to remove..."
		elif (number_files_can_be_removed == 1):
			context['message_remove_files'] = "It is going to remove one file not processed"
			context['message_remove_files_2'] = "Do you want to remove one file not processed?"
		else:
			context['message_remove_files'] = "It is going to remove {} files not processed".format(number_files_can_be_removed)
			context['message_remove_files_2'] = "Do you want to remove {} files not processed?".format(number_files_can_be_removed)
			
		context['show_info_main_page'] = ShowInfoMainPage()		## show main information about the institute
		return context


	def form_valid(self, form):
		"""
		Validate the form
		"""
		return super(SamplesAddFastQView, self).form_valid(form)

	form_valid_message = ""		## need to have this, even empty
	
class SamplesUploadFastQView(LoginRequiredMixin, FormValidMessageMixin, generic.FormView):
	"""
	Create a new reference
	"""
	form_class = SamplesUploadMultipleFastqForm
	success_url = reverse_lazy('sample-add-fastq')
	template_name = 'samples/samples_upload_fastq_files.html'
	utils = Utils()

	if settings.DEBUG: logger = logging.getLogger("fluWebVirus.debug")
	else: logger = logging.getLogger("fluWebVirus.production")
	
	def get_form_kwargs(self):
		"""
		Set the request to pass in the form
		"""
		kw = super(SamplesUploadFastQView, self).get_form_kwargs()
		kw['request'] = self.request 	# get error
		return kw
	
	def get_context_data(self, **kwargs):
		context = super(SamplesUploadFastQView, self).get_context_data(**kwargs)
		context['nav_sample'] = True
		context['nav_modal'] = True	## short the size of modal window
		context['show_info_main_page'] = ShowInfoMainPage()		## show main information about the institute
		context['show_note_message_down_size'] = False
		if (settings.DOWN_SIZE_FASTQ_FILES):
			context['show_note_message_down_size'] = True
			context['message_note_2'] = "Files between {}-{} will be downsized randomly to ~{} before analysis.".format(
				filesizeformat(int(settings.MAX_FASTQ_FILE_UPLOAD)),
				filesizeformat(int(settings.MAX_FASTQ_FILE_WITH_DOWNSIZE)),
				filesizeformat(int(settings.MAX_FASTQ_FILE_UPLOAD))
				)		## show main information about the institute

			context['message_note_1'] = "Maximum size per fastq.gz file is {}.".format(
				filesizeformat(int(settings.MAX_FASTQ_FILE_WITH_DOWNSIZE)))		## show main information about the institute
		else:
			context['message_note_1'] = "Maximum size per fastq.gz file is {}.".format(
				filesizeformat(int(settings.MAX_FASTQ_FILE_UPLOAD)))		## show main information about the institute
			
		### message_note_3, type of files that can be uploaded
		
		return context

	def post(self, request):
		form = SamplesUploadMultipleFastqForm(request.POST, request.FILES, request=request)
		
		data = {}	## return data
		try:
			if form.is_valid():
				
				utils = Utils()
				## doesn't work like that
				#upload_files = form.save()
				
				### get the temporary variable
				path_name = form.cleaned_data['path_name']
				if (path_name is None):
					data = {'is_valid': False, 'name': self.request.FILES['path_name'].name, 'message' : 'Internal server error, path not found.' }
					return JsonResponse(data)

				upload_files = UploadFiles()
				upload_files.file_name = utils.clean_name(ntpath.basename(path_name.name))
				## move the files to the right place
				sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), self.utils.get_path_upload_file(self.request.user.id,\
														TypeFile.TYPE_FILE_fastq_gz), upload_files.file_name)
				sz_file_to, path_added = self.utils.get_unique_file(sz_file_to)		## get unique file name, user can upload files with same name...
				
				## because sometimes has 
				if (str(type(path_name.file)) == "<class '_io.BytesIO'>"):
					temp_file = self.utils.get_temp_file("upload_file", ".dat")
					with open(temp_file, 'wb') as out: ## Open temporary file as bytes
						path_name.file.seek(0)
						out.write(path_name.file.read())                ## Read bytes into file
					self.utils.move_file(temp_file, sz_file_to)
				else: self.utils.copy_file(path_name.file.name, sz_file_to)
				self.logger.info("Starting for file: " + str(upload_files.file_name))
			
				## test if file exist
				if (not os.path.exists(sz_file_to) and os.path.getsize(sz_file_to) > 10):
					data = {'is_valid': False, 'name': self.request.FILES['path_name'].name, 'message' : 'Internal server error, fail to copy file.' }
					return JsonResponse(data)
				
				if path_added is None:
					upload_files.path_name.name = os.path.join(self.utils.get_path_upload_file(self.request.user.id,\
										TypeFile.TYPE_FILE_fastq_gz), ntpath.basename(sz_file_to))
				else:
					upload_files.path_name.name = os.path.join(self.utils.get_path_upload_file(self.request.user.id,\
										TypeFile.TYPE_FILE_fastq_gz), path_added, ntpath.basename(sz_file_to))
				try:
					type_file = MetaKey.objects.get(name=TypeFile.TYPE_FILE_fastq_gz)
				except MetaKey.DoesNotExist:
					type_file = MetaKey()
					type_file.name = TypeFile.TYPE_FILE_fastq_gz
					type_file.save()
	
				upload_files.is_valid = True
				upload_files.is_processed = False			## True when all samples are set
				upload_files.owner = request.user
				upload_files.type_file = type_file
				upload_files.number_files_to_process = 1
				upload_files.number_files_processed = 0
				upload_files.description = ""
				upload_files.save()
	
				data = {'is_valid': True, 'name': upload_files.file_name, 'url': mark_safe(upload_files.get_path_to_file(TypePath.MEDIA_URL)) }
			else:
				data = {'is_valid': False, 'name': self.request.FILES['path_name'].name, 'message' : str(form.errors['path_name'][0]) }
		except:
			self.logger.error(sys.exc_info())
			data = {'is_valid': False, 'name': self.request.FILES['path_name'].name, 'message' : 'Internal server error, unknown error.' }
			return JsonResponse(data)
	
		## if is last file send a message to link files with sample csv file
		if ('is_valid' in data and data['is_valid']): 
			try:
				process_SGE = ProcessSGE()
				taskID = process_SGE.set_link_files(self.request.user)
			except:
				data = {'is_valid': False, 'name': self.request.FILES['path_name'].name, 'message' : 'Fail to submit SGE job.' }
				return JsonResponse(data)
		return JsonResponse(data)
	
	form_valid_message = ""		## need to have this, even empty



class SamplesDetailView(LoginRequiredMixin, DetailView):
    """
    Sample detail view
    """

    utils = Utils()
    software = Software()
    model = Sample
    template_name = "samples/sample_detail.html"

    def get_context_data(self, **kwargs):
        context = super(SamplesDetailView, self).get_context_data(**kwargs)
        sample = kwargs["object"]
        context["nav_sample"] = True
        if sample.owner.id != self.request.user.id:
            context["error_cant_see"] = "1"
            return context

        manageDatabase = ManageDatabase()
        list_meta = manageDatabase.get_sample_metakey(
            sample,
            MetaKeyAndValue.META_KEY_Fastq_Trimmomatic
            if sample.is_type_fastq_gz_sequencing()
            else MetaKeyAndValue.META_KEY_NanoStat_NanoFilt,
            None,
        )
        ## the info can be called from two distinct sources
        ## main info
        if (
            list_meta.count() > 0
            and list_meta[0].value == MetaKeyAndValue.META_VALUE_Success
        ):
            context["virus_identify"] = sample.get_type_sub_type()

            ### fastq original
            file_name = sample.get_fastq(TypePath.MEDIA_ROOT, True)
            if os.path.exists(file_name):
                context["href_fastq_1"] = mark_safe(
                    '<a rel="nofollow" href="'
                    + sample.get_fastq(TypePath.MEDIA_URL, True)
                    + '" download="'
                    + sample.file_name_1
                    + '">'
                    + sample.file_name_1
                    + "</a>"
                )
            else:
                context["href_fastq_1"] = "Not available"

            ### quality first step
            if sample.is_type_fastq_gz_sequencing():
                file_name = sample.get_fastqc_output(
                    TypePath.MEDIA_ROOT, True
                )  ## illumina default reads
            else:
                file_name = sample.get_rabbitQC_output(TypePath.MEDIA_ROOT)
            if os.path.exists(file_name):
                if sample.is_type_fastq_gz_sequencing():
                    context["href_fastq_quality_1"] = mark_safe(
                        '<a rel="nofollow" target="_blank" href="'
                        + sample.get_fastqc_output(TypePath.MEDIA_URL, True)
                        + '">'
                        + sample.file_name_1
                        + ".html</a>"
                    )
                else:
                    context["href_fastq_quality_1"] = mark_safe(
                        '<a rel="nofollow" target="_blank" href="'
                        + sample.get_rabbitQC_output(TypePath.MEDIA_URL)
                        + '">'
                        + sample.file_name_1
                        + ".html</a>"
                    )
            else:
                context["href_fastq_quality_1"] = "Not available"

            ### files to show
            if sample.is_type_fastq_gz_sequencing():  ## illumina default reads
                ### trimmomatic 1
                file_name = sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)
                if os.path.exists(file_name):
                    context["href_trimmonatic_1"] = mark_safe(
                        '<a rel="nofollow" href="'
                        + sample.get_trimmomatic_file(TypePath.MEDIA_URL, True)
                        + '" download="'
                        + os.path.basename(
                            sample.get_trimmomatic_file(TypePath.MEDIA_URL, True)
                        )
                        + '">'
                        + os.path.basename(
                            sample.get_trimmomatic_file(TypePath.MEDIA_URL, True)
                        )
                        + "</a>"
                    )
                else:
                    context["href_trimmonatic_1"] = "Not available"

                ### trimmomatic quality 1
                file_name = sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)
                if os.path.exists(file_name):
                    context["href_trimmonatic_quality_1"] = mark_safe(
                        '<a rel="nofollow" target="_blank" href="'
                        + sample.get_fastq_trimmomatic(TypePath.MEDIA_URL, True)
                        + '">'
                        + os.path.basename(
                            sample.get_fastq_trimmomatic(TypePath.MEDIA_URL, True)
                        )
                        + ".html</a>"
                    )
                else:
                    context["href_trimmonatic_quality_1"] = "Not available"

                ### testing second file
                trimmomatic_file_name = sample.get_trimmomatic_file(
                    TypePath.MEDIA_URL, False
                )
                if not trimmomatic_file_name is None:
                    ### trimmomatic fastq second
                    file_name = sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False)
                    if os.path.exists(file_name):
                        context["href_trimmonatic_2"] = mark_safe(
                            '<a rel="nofollow" href="'
                            + sample.get_trimmomatic_file(TypePath.MEDIA_URL, False)
                            + '" download="'
                            + os.path.basename(
                                sample.get_trimmomatic_file(TypePath.MEDIA_URL, False)
                            )
                            + '">'
                            + os.path.basename(
                                sample.get_trimmomatic_file(TypePath.MEDIA_URL, False)
                            )
                            + "</a>"
                        )
                    else:
                        context["href_trimmonatic_2"] = "Not available"

                    ### trimmomatic fastqc second
                    file_name = sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)
                    if os.path.exists(file_name):
                        context["href_trimmonatic_quality_2"] = mark_safe(
                            '<a rel="nofollow" target="_blank" href="'
                            + sample.get_fastq_trimmomatic(TypePath.MEDIA_URL, False)
                            + '">'
                            + os.path.basename(
                                sample.get_fastq_trimmomatic(TypePath.MEDIA_URL, False)
                            )
                            + ".html</a>"
                        )
                    else:
                        context["href_trimmonatic_quality_2"] = "Not available"

                    ### fastq second
                    file_name = sample.get_fastq(TypePath.MEDIA_ROOT, False)
                    if os.path.exists(file_name):
                        context["href_fastq_2"] = mark_safe(
                            '<a rel="nofollow" href="'
                            + sample.get_fastq(TypePath.MEDIA_URL, False)
                            + '" download="'
                            + sample.file_name_2
                            + '">'
                            + sample.file_name_2
                            + "</a>"
                        )
                    else:
                        context["href_fastq_2"] = "Not available"
                    ### fastqc second
                    file_name = sample.get_fastqc_output(TypePath.MEDIA_ROOT, False)
                    if os.path.exists(file_name):
                        context["href_fastq_quality_2"] = mark_safe(
                            '<a rel="nofollow" target="_blank" href="'
                            + sample.get_fastqc_output(TypePath.MEDIA_URL, False)
                            + '"">'
                            + sample.file_name_2
                            + ".html</a>"
                        )
                    else:
                        context["href_fastq_quality_2"] = "Not available"
                else:  ## there's no second file
                    context["href_fastq_2"] = "Not available"
                    context["href_fastq_quality_2"] = "Not available"
                    context["href_trimmonatic_2"] = "Not available"
                    context["href_trimmonatic_quality_2"] = "Not available"

                #### data from illumina stat
                stat_data, total_reads = self.software.get_stats_from_sample_reads(
                    sample
                )
                if not stat_data is None:
                    context["data_illuminastat"] = stat_data

            else:  ### other like Minion
                file_name = sample.get_nanofilt_file(TypePath.MEDIA_ROOT)
                if os.path.exists(file_name):
                    context["href_trimmonatic_1"] = mark_safe(
                        '<a rel="nofollow" href="'
                        + sample.get_nanofilt_file(TypePath.MEDIA_URL)
                        + '" download="'
                        + os.path.basename(sample.get_nanofilt_file(TypePath.MEDIA_URL))
                        + '">'
                        + os.path.basename(sample.get_nanofilt_file(TypePath.MEDIA_URL))
                        + "</a>"
                    )
                else:
                    context["href_trimmonatic_1"] = "Not available"
                file_name = sample.get_rabbitQC_nanofilt(TypePath.MEDIA_ROOT)
                if os.path.exists(file_name):
                    context["href_trimmonatic_quality_1"] = mark_safe(
                        '<a rel="nofollow" target="_blank" href="'
                        + sample.get_rabbitQC_nanofilt(TypePath.MEDIA_URL)
                        + '">'
                        + os.path.basename(
                            sample.get_rabbitQC_nanofilt(TypePath.MEDIA_URL)
                        )
                        + "</a>"
                    )
                else:
                    context["href_trimmonatic_quality_1"] = "Not available"

                #### data from nanoStat
                stat_data, total_reads = self.software.get_stats_from_sample_reads(
                    sample
                )
                if not stat_data is None:
                    context["data_nanostat"] = stat_data

            ### software
            if sample.is_type_fastq_gz_sequencing():  ### for illumina

                meta_sample = manageDatabase.get_sample_metakey_last(
                    sample,
                    MetaKeyAndValue.META_KEY_Fastq_Trimmomatic_Software,
                    MetaKeyAndValue.META_VALUE_Success,
                )
                if meta_sample is None:
                    meta_sample = manageDatabase.get_sample_metakey_last(
                        sample,
                        MetaKeyAndValue.META_KEY_Fastq_Trimmomatic,
                        MetaKeyAndValue.META_VALUE_Success,
                    )
                    if not meta_sample is None:
                        lst_data = meta_sample.description.split(",")
                        context["fastq_software"] = lst_data[1].strip()
                        context["trimmomatic_software"] = lst_data[2].strip()
                else:
                    decodeResult = DecodeObjects()
                    result = decodeResult.decode_result(meta_sample.description)
                    context["fastq_software"] = result.get_software(
                        SoftwareNames.SOFTWARE_FASTQ_name
                    )
                    trimmomatic_software = result.get_software(
                        SoftwareNames.SOFTWARE_TRIMMOMATIC_name
                    )
                    context["trimmomatic_software"] = (
                        trimmomatic_software
                        if len(trimmomatic_software) > 0
                        else "Trimmomatic not ran."
                    )

                ### species identification, only in Illumina
                meta_sample = manageDatabase.get_sample_metakey_last(
                    sample,
                    MetaKeyAndValue.META_KEY_Identify_Sample_Software,
                    MetaKeyAndValue.META_VALUE_Success,
                )
                if meta_sample is None:
                    meta_sample = manageDatabase.get_sample_metakey_last(
                        sample,
                        MetaKeyAndValue.META_KEY_Identify_Sample,
                        MetaKeyAndValue.META_VALUE_Success,
                    )
                    if meta_sample != None:
                        lst_data = meta_sample.description.split(",")
                        context["spades_software"] = lst_data[1].strip()
                        context["abricate_software"] = lst_data[2].strip()
                else:
                    decodeResult = DecodeObjects()
                    result = decodeResult.decode_result(meta_sample.description)
                    context["spades_software"] = result.get_software(
                        SoftwareNames.SOFTWARE_SPAdes_name
                    )
                    context["abricate_software"] = result.get_software(
                        SoftwareNames.SOFTWARE_ABRICATE_name
                    )

            else:  ### Software Minion
                meta_sample = manageDatabase.get_sample_metakey_last(
                    sample,
                    MetaKeyAndValue.META_KEY_NanoStat_NanoFilt_Software,
                    MetaKeyAndValue.META_VALUE_Success,
                )
                if meta_sample == None:
                    meta_sample = manageDatabase.get_sample_metakey_last(
                        sample,
                        MetaKeyAndValue.META_KEY_NanoStat_NanoFilt,
                        MetaKeyAndValue.META_VALUE_Success,
                    )
                    if not meta_sample is None:
                        lst_data = meta_sample.description.split(",")
                        context["fastq_software"] = lst_data[1].strip()
                        context["trimmomatic_software"] = lst_data[2].strip()
                else:
                    decodeResult = DecodeObjects()
                    result = decodeResult.decode_result(meta_sample.description)
                    context["fastq_software"] = result.get_software(
                        SoftwareNames.SOFTWARE_NanoStat_name
                    )
                    nanofilt_software = result.get_software(
                        SoftwareNames.SOFTWARE_NanoFilt_name
                    )
                    context["trimmomatic_software"] = (
                        nanofilt_software
                        if len(nanofilt_software) > 0
                        else "Nanofilt not ran."
                    )

                ### abricate type/subtype
                # meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Identify_Sample, MetaKeyAndValue.META_VALUE_Success)
                # if (not meta_sample is None):
                #    lst_data = meta_sample.description.split(',')
                #    if len(lst_data) > 0: context['abricate_software'] = lst_data[1].strip()
                ### species identification, only in Illumina
                meta_sample = manageDatabase.get_sample_metakey_last(
                    sample,
                    MetaKeyAndValue.META_KEY_Identify_Sample_Software,
                    MetaKeyAndValue.META_VALUE_Success,
                )
                if meta_sample is None:
                    meta_sample = manageDatabase.get_sample_metakey_last(
                        sample,
                        MetaKeyAndValue.META_KEY_Identify_Sample,
                        MetaKeyAndValue.META_VALUE_Success,
                    )
                    if meta_sample != None:
                        lst_data = meta_sample.description.split(",")
                        context["spades_software"] = lst_data[1].strip()
                        context["abricate_software"] = lst_data[2].strip()
                else:
                    decodeResult = DecodeObjects()
                    result = decodeResult.decode_result(meta_sample.description)
                    context["spades_software"] = result.get_software(
                        SoftwareNames.SOFTWARE_FLYE_name
                    )
                    context["abricate_software"] = result.get_software(
                        SoftwareNames.SOFTWARE_ABRICATE_name
                    )

            ##### extra data sample, columns added by the user
            ## [[header1, value1], [header2, value2], [header3, value3], ...]
            ### if it's to big expand button is better
            tag_names = sample.get_tag_names()
            context[
                "extra_data_sample_expand"
            ] = tag_names != None and tag_names.count() > (
                Constants.START_EXPAND_SAMPLE_TAG_NAMES_ROWS
            )
            if tag_names != None:
                context["extra_data_sample"] = self.utils.grouped(tag_names, 4)

            ### get different types of alerts
            metaKeyAndValue = MetaKeyAndValue()
            alert_out = []
            for key in metaKeyAndValue.get_keys_show_alerts_in_sample_details_view():
                meta_data = manageDatabase.get_sample_metakey_last(
                    sample, key, MetaKeyAndValue.META_VALUE_Success
                )
                if meta_data != None:
                    alert_out.append(meta_data.description)
            context["alerts"] = alert_out
            context["has_type_subtype"] = sample.identify_virus.all().count() > 0
            if sample.identify_virus.all().count() > 0:
                vect_identify_virus = []
                for identify_virus in sample.identify_virus.all():
                    if not identify_virus in vect_identify_virus:
                        vect_identify_virus.append(identify_virus)
                context["vect_identify_virus"] = vect_identify_virus

            ## files with contigs
            if sample.is_type_fastq_gz_sequencing():
                if os.path.exists(
                    sample.get_draft_contigs_output(TypePath.MEDIA_ROOT)
                ) and os.path.exists(
                    sample.get_draft_contigs_abricate_output(TypePath.MEDIA_ROOT)
                ):
                    context["file_draft_contigs"] = mark_safe(
                        '<a rel="nofollow" href="'
                        + sample.get_draft_contigs_output(TypePath.MEDIA_URL)
                        + '" download="'
                        + os.path.basename(
                            sample.get_draft_contigs_output(TypePath.MEDIA_ROOT)
                        )
                        + '">'
                        + os.path.basename(
                            sample.get_draft_contigs_output(TypePath.MEDIA_ROOT)
                        )
                        + "</a>"
                    )
                    context["file_draft_contigs_abricate"] = mark_safe(
                        '<a rel="nofollow" href="'
                        + sample.get_draft_contigs_abricate_output(TypePath.MEDIA_URL)
                        + '" download="'
                        + os.path.basename(
                            sample.get_draft_contigs_abricate_output(
                                TypePath.MEDIA_ROOT
                            )
                        )
                        + '">'
                        + os.path.basename(
                            sample.get_draft_contigs_abricate_output(
                                TypePath.MEDIA_ROOT
                            )
                        )
                        + "</a>"
                    )
                    context["has_draft_contigs"] = True
            elif os.path.exists(
                sample.get_draft_reads_abricate_output(TypePath.MEDIA_ROOT)
            ):
                context["file_draft_reads_abricate"] = mark_safe(
                    '<a rel="nofollow" href="'
                    + sample.get_draft_reads_abricate_output(TypePath.MEDIA_URL)
                    + '" download="'
                    + os.path.basename(
                        sample.get_draft_reads_abricate_output(TypePath.MEDIA_ROOT)
                    )
                    + '">'
                    + os.path.basename(
                        sample.get_draft_reads_abricate_output(TypePath.MEDIA_ROOT)
                    )
                    + "</a>"
                )
                context["has_draft_reads"] = True

        elif (
            sample.candidate_file_name_1 != None
            and len(sample.candidate_file_name_1) > 0
        ):
            context["candidate_file_name_1"] = sample.candidate_file_name_1
            if (
                sample.candidate_file_name_2 != None
                and len(sample.candidate_file_name_2) > 0
            ):
                context["candidate_file_name_2"] = sample.candidate_file_name_2

        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        return context


class ProjectsView(LoginRequiredMixin, ListView):
    utils = Utils()
    model = Project
    template_name = "project/projects.html"
    context_object_name = "projects"
    ##    group_required = u'company-user' security related with GroupRequiredMixin

    def get_context_data(self, **kwargs):
        context = super(ProjectsView, self).get_context_data(**kwargs)
        tag_search = "search_projects"
        query_set = Project.objects.filter(
            owner__id=self.request.user.id, is_deleted=False
        ).order_by("-creation_date")
        if self.request.GET.get(tag_search) != None and self.request.GET.get(
            tag_search
        ):
            query_set = query_set.filter(
                Q(name__icontains=self.request.GET.get(tag_search))
                | Q(reference__name__icontains=self.request.GET.get(tag_search))
                | Q(
                    project_samples__sample__name__icontains=self.request.GET.get(
                        tag_search
                    )
                )
            ).distinct()

        table = ProjectTable(query_set)
        RequestConfig(
            self.request, paginate={"per_page": Constants.PAGINATE_NUMBER}
        ).configure(table)
        if self.request.GET.get(tag_search) != None:
            context[tag_search] = self.request.GET.get(tag_search)

        ### clean check box in the session
        clean_check_box_in_session(self.request)  ## for both samples and references

        ### clean project name session
        if Constants.PROJECT_NAME_SESSION in self.request.session:
            del self.request.session[Constants.PROJECT_NAME_SESSION]

        context["table"] = table
        context["nav_project"] = True
        context["show_paginatior"] = query_set.count() > Constants.PAGINATE_NUMBER
        context["query_set_count"] = query_set.count()

        ## project list
        ## list of all samples in CSV and TSV
        csv_file = self.utils.get_project_list_by_user(
            self.request.user.id, "MEDIA_ROOT", FileExtensions.FILE_CSV
        )
        if os.path.exists(csv_file):
            context["list_project_file_csv"] = mark_safe(
                '<a rel="nofollow" href="'
                + self.utils.get_project_list_by_user(
                    self.request.user.id, "MEDIA_URL", FileExtensions.FILE_CSV
                )
                + '" download="'
                + os.path.basename(csv_file)
                + '" class="dropdown-item"> Download - '
                + os.path.basename(csv_file)
                + "</a>"
            )
        tsv_file = self.utils.get_project_list_by_user(
            self.request.user.id, "MEDIA_ROOT", FileExtensions.FILE_TSV
        )
        if os.path.exists(tsv_file):
            context["list_project_file_tsv"] = mark_safe(
                '<a rel="nofollow" href="'
                + self.utils.get_project_list_by_user(
                    self.request.user.id, "MEDIA_URL", FileExtensions.FILE_TSV
                )
                + '" download="'
                + os.path.basename(tsv_file)
                + '" class="dropdown-item"> Download - '
                + os.path.basename(tsv_file)
                + "</a>"
            )

        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        return context


class ProjectCreateView(LoginRequiredMixin, FormValidMessageMixin, generic.CreateView):
    """
    Create a new reference
    """

    utils = Utils()
    model = Project
    fields = ["name"]
    success_url = reverse_lazy("projects")
    template_name = "project/project_add.html"

    def get_form_kwargs(self):
        """
        Set the request to pass in the form
        """
        kw = super(ProjectCreateView, self).get_form_kwargs()
        return kw

    def get_context_data(self, **kwargs):
        context = super(ProjectCreateView, self).get_context_data(**kwargs)

        ### get project name
        project_name = (
            self.request.session[Constants.PROJECT_NAME_SESSION]
            if Constants.PROJECT_NAME_SESSION in self.request.session
            else ""
        )
        tag_search = "search_references"
        query_set = Reference.objects.filter(
            owner__id=self.request.user.id, is_obsolete=False, is_deleted=False
        ).order_by("-name")
        if self.request.GET.get(tag_search) != None and self.request.GET.get(
            tag_search
        ):
            query_set = query_set.filter(
                Q(name__icontains=self.request.GET.get(tag_search))
                | Q(owner__username__icontains=self.request.GET.get(tag_search))
                | Q(isolate_name__icontains=self.request.GET.get(tag_search))
            )

        ### get the references from the system
        query_set_system = Reference.objects.filter(
            owner__username=Constants.DEFAULT_USER, is_obsolete=False, is_deleted=False
        ).order_by("-name")
        if self.request.GET.get(tag_search) != None and self.request.GET.get(
            tag_search
        ):
            query_set_system = query_set_system.filter(
                Q(name__icontains=self.request.GET.get(tag_search))
                | Q(owner__username__icontains=self.request.GET.get(tag_search))
                | Q(isolate_name__icontains=self.request.GET.get(tag_search))
            )
        query_set_result = sorted(
            chain(query_set, query_set_system), key=attrgetter("creation_date")
        )
        table = ReferenceProjectTable(query_set_result)
        RequestConfig(
            self.request, paginate={"per_page": Constants.PAGINATE_NUMBER}
        ).configure(table)

        ### set the check_box, check_box_all is only to control if is the first time or not
        if Constants.CHECK_BOX_ALL not in self.request.session:
            self.request.session[Constants.CHECK_BOX_ALL] = False
            is_all_check_box_in_session(
                [
                    "{}_{}".format(Constants.CHECK_BOX, key.id)
                    for key in query_set_result
                ],
                self.request,
            )
        elif "search_references" in self.request.GET:
            # clean check boxes in search
            dt_sample_id_add_temp = {}
            for reference in query_set_result:
                dt_sample_id_add_temp[
                    reference.id
                ] = 1  ## add the ids that are in the tables
            for key in self.request.session.keys():
                if (
                    key.startswith(Constants.CHECK_BOX)
                    and len(key.split("_")) == 3
                    and self.utils.is_integer(key.split("_")[2])
                ):
                    ### this is necessary because of the search. Can occur some checked box that are out of filter.
                    if int(key.split("_")[2]) not in dt_sample_id_add_temp:
                        self.request.session[key] = False

        ## clean the
        if self.request.GET.get(tag_search) != None:
            context[tag_search] = self.request.GET.get(tag_search)
        if Constants.PROJECT_NAME in self.request.session:
            context[Constants.PROJECT_NAME] = self.request.session[
                Constants.PROJECT_NAME
            ]
            del self.request.session[Constants.PROJECT_NAME]
        if Constants.ERROR_PROJECT_NAME in self.request.session:
            context[Constants.ERROR_PROJECT_NAME] = self.request.session[
                Constants.ERROR_PROJECT_NAME
            ]
            del self.request.session[Constants.ERROR_PROJECT_NAME]
        if Constants.ERROR_REFERENCE in self.request.session:
            context[Constants.ERROR_REFERENCE] = self.request.session[
                Constants.ERROR_REFERENCE
            ]
            del self.request.session[Constants.ERROR_REFERENCE]

        context["project_name"] = project_name
        context["table"] = table
        context["show_paginatior"] = len(query_set_result) > Constants.PAGINATE_NUMBER
        context["nav_project"] = True
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        if self.request.POST:
            context["project_references"] = ReferenceProjectFormSet(self.request.POST)
        else:
            context["project_references"] = ReferenceProjectFormSet()
        return context

    def form_valid(self, form):
        """
        Validate the form
        """
        ### test anonymous account
        try:
            profile = Profile.objects.get(user=self.request.user)
            if profile.only_view_project:
                messages.warning(
                    self.request,
                    "'{}' account can not create projects.".format(
                        self.request.user.username
                    ),
                    fail_silently=True,
                )
                return super(ProjectCreateView, self).form_invalid(form)
        except Profile.DoesNotExist:
            pass

        context = self.get_context_data()
        project_references = context["project_references"]
        name = form.cleaned_data["name"]
        b_error = False
        try:
            Project.objects.get(
                name__iexact=name,
                is_deleted=False,
                owner__username=self.request.user.username,
            )
            self.request.session[
                Constants.ERROR_PROJECT_NAME
            ] = "Exists a project with this name."
            self.request.session[Constants.PROJECT_NAME] = name
            b_error = True
        except Project.DoesNotExist:
            pass

        select_ref = None
        if "select_ref" in project_references.data:
            select_ref = project_references.data["select_ref"]
        else:
            ### test if it's in the session
            select_ref = get_unique_pk_from_session(self.request)
            if select_ref == None:
                self.request.session[
                    Constants.ERROR_REFERENCE
                ] = "You need to select a reference."
                b_error = True
        if name is None:
            name = ""

        if select_ref != None:
            try:
                reference = Reference.objects.get(pk=select_ref)
            except Reference.DoesNotExist:
                self.request.session[
                    Constants.ERROR_REFERENCE
                ] = "You need to select a reference."
                self.request.session[Constants.PROJECT_NAME] = name
                return super(ProjectCreateView, self).form_invalid(form)

            if reference.is_deleted:
                self.request.session[
                    Constants.ERROR_REFERENCE
                ] = "The reference '{}' was removed.".format(reference.name)
                self.request.session[Constants.PROJECT_NAME] = name
                return super(ProjectCreateView, self).form_invalid(form)
        else:
            self.request.session[
                Constants.ERROR_REFERENCE
            ] = "You need to select a reference."
            self.request.session[Constants.PROJECT_NAME] = name
            return super(ProjectCreateView, self).form_invalid(form)

        ### exists an error
        if b_error:
            return super(ProjectCreateView, self).form_invalid(form)

        with transaction.atomic():
            project = form.save()
            project.reference = reference
            project.owner = self.request.user
            project.save()

        ### test all defaults first, if exist in database
        default_software = DefaultProjectSoftware()
        default_software.test_all_defaults(
            self.request.user, project, None, None
        )  ## the user can have defaults yet

        process_SGE = ProcessSGE()
        process_SGE.set_create_project_list_by_user(self.request.user)

        messages.success(
            self.request,
            "Project '" + name + "' was created successfully",
            fail_silently=True,
        )
        return super(ProjectCreateView, self).form_valid(form)

    form_valid_message = ""  ## need to have this, even empty


class AddSamplesProjectsView(
    LoginRequiredMixin, FormValidMessageMixin, generic.CreateView
):
    """
    Create a new reference
    """

    utils = Utils()
    model = Sample
    fields = ["name"]
    success_url = reverse_lazy("projects")
    template_name = "project_sample/project_sample_add.html"

    if settings.DEBUG:
        logger = logging.getLogger("fluWebVirus.debug")
    else:
        logger = logging.getLogger("fluWebVirus.production")

    def get_context_data(self, **kwargs):
        context = super(AddSamplesProjectsView, self).get_context_data(**kwargs)

        ### test if the user is the same of the page
        project = Project.objects.get(pk=self.kwargs["pk"])
        context["nav_project"] = True
        if project.owner.id != self.request.user.id:
            context["error_cant_see"] = "1"
            return context

        ## catch everything that is not in connection with project
        count_active_projects = ProjectSample.objects.filter(
            project=project, is_deleted=False, is_error=False
        ).count()
        samples_out = ProjectSample.objects.filter(
            Q(project=project) & ~Q(is_deleted=True) & ~Q(is_error=True)
        ).values("sample__pk")
        query_set = Sample.objects.filter(
            owner__id=self.request.user.id,
            is_obsolete=False,
            is_deleted=False,
            is_deleted_processed_fastq=False,
            is_ready_for_projects=True,
        ).exclude(pk__in=samples_out)

        tag_search = "search_add_project_sample"
        if self.request.GET.get(tag_search) != None and self.request.GET.get(
            tag_search
        ):
            query_set = query_set.filter(
                Q(name__icontains=self.request.GET.get(tag_search))
                | Q(type_subtype__icontains=self.request.GET.get(tag_search))
                | Q(data_set__name__icontains=self.request.GET.get(tag_search))
                | Q(week__icontains=self.request.GET.get(tag_search))
            )
            tag_search
        table = SampleToProjectsTable(query_set)

        ### set the check_box
        if Constants.CHECK_BOX_ALL not in self.request.session:
            self.request.session[Constants.CHECK_BOX_ALL] = False
            is_all_check_box_in_session(
                ["{}_{}".format(Constants.CHECK_BOX, key.id) for key in query_set],
                self.request,
            )

        context[Constants.CHECK_BOX_ALL] = self.request.session[Constants.CHECK_BOX_ALL]
        ## need to clean all the others if are reject in filter
        dt_sample_id_add_temp = {}
        if context[Constants.CHECK_BOX_ALL]:
            for sample in query_set:
                dt_sample_id_add_temp[
                    sample.id
                ] = 1  ## add the ids that are in the tables
            for key in self.request.session.keys():
                if (
                    key.startswith(Constants.CHECK_BOX)
                    and len(key.split("_")) == 3
                    and self.utils.is_integer(key.split("_")[2])
                ):
                    ### this is necessary because of the search. Can occur some checked box that are out of filter.
                    if int(key.split("_")[2]) not in dt_sample_id_add_temp:
                        self.request.session[key] = False
                    else:
                        self.request.session[key] = True
        ## END need to clean all the others if are reject in filter

        ### check if it was shown the settings message
        key_session_name_project_settings = "project_settings_{}".format(project.name)
        if not key_session_name_project_settings in self.request.session:
            self.request.session[key_session_name_project_settings] = True
        else:
            self.request.session[key_session_name_project_settings] = False

        RequestConfig(
            self.request, paginate={"per_page": Constants.PAGINATE_NUMBER}
        ).configure(table)
        if self.request.GET.get(tag_search) != None:
            context[tag_search] = self.request.GET.get(tag_search)
        context["table"] = table
        context["show_paginatior"] = query_set.count() > Constants.PAGINATE_NUMBER
        context["query_set_count"] = query_set.count()
        context["project_name"] = project.name
        context["nav_modal"] = True  ## short the size of modal window
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        context["show_message_change_settings"] = (
            count_active_projects == 0
            and self.request.session[key_session_name_project_settings]
        )  ## Show message to change settings to the project
        context["add_all_samples_message"] = "Add {} sample{}".format(
            query_set.count(), pluralize(query_set.count(), ",s")
        )
        if self.request.POST:
            context["project_sample"] = AddSampleProjectForm(self.request.POST)
        else:
            context["project_sample"] = AddSampleProjectForm()
        return context

    def form_valid(self, form):
        """
        Validate the form
        """

        ### test anonymous account
        try:
            profile = Profile.objects.get(user=self.request.user)
            if profile.only_view_project:
                messages.warning(
                    self.request,
                    "'{}' account can not add samples to a project.".format(
                        self.request.user.username
                    ),
                    fail_silently=True,
                )
                return super(AddSamplesProjectsView, self).form_invalid(form)
        except Profile.DoesNotExist:
            pass

        if form.is_valid():
            metaKeyAndValue = MetaKeyAndValue()
            manageDatabase = ManageDatabase()
            process_SGE = ProcessSGE()

            ### get project sample..
            context = self.get_context_data()

            ### get project
            project = Project.objects.get(pk=self.kwargs["pk"])

            vect_sample_id_add_temp = []
            for sample in context["table"].data:
                vect_sample_id_add_temp.append(sample.id)

            vect_sample_id_add = []
            if "submit_checked" in self.request.POST:
                for key in self.request.session.keys():
                    if (
                        self.request.session[key]
                        and key.startswith(Constants.CHECK_BOX)
                        and len(key.split("_")) == 3
                        and self.utils.is_integer(key.split("_")[2])
                    ):
                        ### this is necessary because of the search. Can occur some checked box that are out of filter.
                        if int(key.split("_")[2]) in vect_sample_id_add_temp:
                            vect_sample_id_add.append(int(key.split("_")[2]))
            elif "submit_all" in self.request.POST:
                vect_sample_id_add = vect_sample_id_add_temp

            ### get samples already out
            dt_sample_out = {}
            for project_sample in ProjectSample.objects.filter(
                project__id=project.id,
                is_deleted=False,
                is_error=False,
                is_deleted_in_file_system=False,
            ):
                dt_sample_out[project_sample.sample.id] = 1

            ### start adding...
            (job_name_wait, job_name) = ("", "")
            project_sample_add = 0
            for id_sample in vect_sample_id_add:
                ## keep track of samples out
                if id_sample in dt_sample_out:
                    continue
                dt_sample_out[id_sample] = 1

                try:
                    sample = Sample.objects.get(pk=id_sample)
                except Sample.DoesNotExist:
                    ## log
                    self.logger.error(
                        "Fail to get sample_id {} in ProjectSample".format(
                            key.split("_")[2]
                        )
                    )
                    continue

                ## the sample can be deleted by other session
                if sample.is_deleted:
                    continue

                ## get project sample
                try:
                    project_sample = ProjectSample.objects.get(
                        project__id=project.id, sample__id=sample.id
                    )

                    ### if exist can be deleted, pass to active
                    if (
                        project_sample.is_deleted
                        and not project_sample.is_error
                        and not project_sample.is_deleted_in_file_system
                    ):
                        project_sample.is_deleted = False
                        project_sample.is_finished = False
                        project_sample.save()
                        project_sample_add += 1

                except ProjectSample.DoesNotExist:
                    project_sample = ProjectSample()
                    project_sample.project = project
                    project_sample.sample = sample
                    project_sample.save()
                    project_sample_add += 1

                ### create a task to perform the analysis of snippy and freebayes
                ### Important, it is necessary to run again because can have some changes in the parameters.
                try:
                    if len(job_name_wait) == 0:
                        (
                            job_name_wait,
                            job_name,
                        ) = self.request.user.profile.get_name_sge_seq(
                            Profile.SGE_PROCESS_projects, Profile.SGE_GLOBAL
                        )
                    if sample.is_type_fastq_gz_sequencing():
                        taskID = process_SGE.set_second_stage_snippy(
                            project_sample, self.request.user, job_name, [job_name_wait]
                        )
                    else:
                        taskID = process_SGE.set_second_stage_medaka(
                            project_sample, self.request.user, job_name, [job_name_wait]
                        )

                    ### set project sample queue ID
                    manageDatabase.set_project_sample_metakey(
                        project_sample,
                        self.request.user,
                        metaKeyAndValue.get_meta_key_queue_by_project_sample_id(
                            project_sample.id
                        ),
                        MetaKeyAndValue.META_VALUE_Queue,
                        taskID,
                    )
                except:
                    pass

            ### necessary to calculate the global results again
            if project_sample_add > 0:
                try:
                    taskID = process_SGE.set_collect_global_files(
                        project, self.request.user
                    )
                    manageDatabase.set_project_metakey(
                        project,
                        self.request.user,
                        metaKeyAndValue.get_meta_key(
                            MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id
                        ),
                        MetaKeyAndValue.META_VALUE_Queue,
                        taskID,
                    )
                except:
                    pass

            if project_sample_add == 0:
                messages.warning(
                    self.request,
                    "No sample was added to the project '{}'".format(project.name),
                )
            else:
                if project_sample_add > 1:
                    messages.success(
                        self.request,
                        "'{}' samples were added to your project.".format(
                            project_sample_add
                        ),
                        fail_silently=True,
                    )
                else:
                    messages.success(
                        self.request,
                        "One sample was added to your project.",
                        fail_silently=True,
                    )
            return HttpResponseRedirect(reverse_lazy("projects"))
        else:
            return super(AddSamplesProjectsView, self).form_invalid(form)

    form_valid_message = ""  ## need to have this, even empty


class ShowSampleProjectsView(LoginRequiredMixin, ListView):
    model = Project
    template_name = "project_sample/project_sample_show.html"
    context_object_name = "project_sample"

    def get_context_data(self, **kwargs):
        context = super(ShowSampleProjectsView, self).get_context_data(**kwargs)
        project = Project.objects.get(pk=self.kwargs["pk"])
        software_pangolin = SoftwarePangolin()
        software = Software()

        ### can't see this project
        context["nav_project"] = True
        if project.owner.id != self.request.user.id:
            context["error_cant_see"] = "1"
            return context

        query_set = ProjectSample.objects.filter(
            project__id=project.id, is_finished=True, is_deleted=False, is_error=False
        ).order_by("creation_date")
        tag_search = "search_add_project_sample"
        ### filter the search
        if self.request.GET.get(tag_search) != None and self.request.GET.get(
            tag_search
        ):
            query_set = query_set.filter(
                Q(sample__name__icontains=self.request.GET.get(tag_search))
                | Q(sample__data_set__name__icontains=self.request.GET.get(tag_search))
                | Q(
                    mixed_infections__tag__name__icontains=self.request.GET.get(
                        tag_search
                    )
                )
                | Q(sample__type_subtype__icontains=self.request.GET.get(tag_search))
            )
        table = ShowProjectSamplesResults(query_set)
        RequestConfig(
            self.request, paginate={"per_page": Constants.PAGINATE_NUMBER}
        ).configure(table)

        if self.request.GET.get(tag_search) != None:
            context[tag_search] = self.request.GET.get("search_add_project_sample")
        context["table"] = table
        context["show_paginatior"] = query_set.count() > Constants.PAGINATE_NUMBER
        context["query_set_count"] = query_set.count()
        context["project_id"] = project.id
        context["spinner_url"] = os.path.join(
            "/" + Constants.DIR_STATIC, Constants.DIR_ICONS, Constants.AJAX_LOADING_GIF
        )
        context["nav_project"] = True
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute

        context["project_name"] = project.name
        context["reference_name"] = project.reference.get_reference_fasta_web()
        context["number_of_samples"] = ProjectSample.objects.filter(
            project=project, is_deleted=False, is_error=False, is_finished=True
        ).count()
        context["number_of_alerts"] = ProjectSample.objects.filter(
            Q(project=project)
            & Q(is_deleted=False)
            & Q(is_error=False)
            & Q(is_finished=True)
            & (Q(alert_first_level__gte=1) | Q(alert_second_level__gte=1))
        ).count()
        context["samples_in_process"] = ProjectSample.objects.filter(
            project=project, is_deleted=False, is_error=False, is_finished=False
        ).count()
        context["samples_error"] = ProjectSample.objects.filter(
            project=project, is_deleted=False, is_error=True, is_finished=False
        ).count()

        ## Files
        context["coverage_file"] = project.get_global_file_by_project_web(
            Project.PROJECT_FILE_NAME_COVERAGE
        )
        context["main_variations_snippy_file"] = project.get_global_file_by_project_web(
            Project.PROJECT_FILE_NAME_TAB_VARIATIONS_SNIPPY
        )

        ## coverage
        if os.path.exists(
            project.get_global_file_by_project(
                TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_COVERAGE
            )
        ):
            context["samples_file_coverage"] = get_link_for_dropdown_item(
                project.get_global_file_by_project(
                    TypePath.MEDIA_URL, Project.PROJECT_FILE_NAME_COVERAGE
                )
            )
        ## variants
        if os.path.exists(
            project.get_global_file_by_project(
                TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_TAB_VARIATIONS_SNIPPY
            )
        ):
            context["samples_file_variants"] = get_link_for_dropdown_item(
                project.get_global_file_by_project(
                    TypePath.MEDIA_URL, Project.PROJECT_FILE_NAME_TAB_VARIATIONS_SNIPPY
                )
            )
        ## minor intra host
        if os.path.exists(
            project.get_global_file_by_project(
                TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_TAB_VARIATIONS_FREEBAYES
            )
        ):
            context["samples_file_minor_intra_host"] = get_link_for_dropdown_item(
                project.get_global_file_by_project(
                    TypePath.MEDIA_URL,
                    Project.PROJECT_FILE_NAME_TAB_VARIATIONS_FREEBAYES,
                )
            )

        ## aln2pheno result
        if os.path.exists(
            project.get_global_file_by_project(
                TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_Aln2pheno_report_COG_UK
            )
        ):
            context["aln2pheno_report_cog"] = get_link_for_dropdown_item(
                project.get_global_file_by_project(
                    TypePath.MEDIA_URL,
                    Project.PROJECT_FILE_NAME_Aln2pheno_report_COG_UK,
                )
            )
        # if os.path.exists(project.get_global_file_by_project(TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_Aln2pheno_report_pokay)):
        #    context['aln2pheno_report_pokay'] = get_link_for_dropdown_item(
        #        project.get_global_file_by_project(TypePath.MEDIA_URL, Project.PROJECT_FILE_NAME_Aln2pheno_report_pokay))
        if os.path.exists(
            project.get_global_file_by_project(
                TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_Aln2pheno_zip
            )
        ):
            context["aln2pheno_zip"] = get_link_for_dropdown_item(
                project.get_global_file_by_project(
                    TypePath.MEDIA_URL, Project.PROJECT_FILE_NAME_Aln2pheno_zip
                )
            )

        ## all files zipped
        if os.path.exists(
            project.get_global_file_by_project(
                TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_all_files_zipped
            )
        ):
            context["download_all_files"] = get_link_for_dropdown_item(
                project.get_global_file_by_project(
                    TypePath.MEDIA_URL, Project.PROJECT_FILE_NAME_all_files_zipped
                ),
                "{}_{}_{}".format(
                    os.path.splitext(Project.PROJECT_FILE_NAME_all_files_zipped)[0],
                    project.get_clean_project_name(),
                    datetime.now().strftime(settings.DATE_FORMAT_FOR_SHOW),
                ),
            )

        if os.path.exists(
            project.get_global_file_by_project(
                TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV
            )
        ):
            context["sample_file_result_csv"] = get_link_for_dropdown_item(
                project.get_global_file_by_project(
                    TypePath.MEDIA_URL, Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV
                )
            )
        if os.path.exists(
            project.get_global_file_by_project(
                TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_SAMPLE_RESULT_TSV
            )
        ):
            context["sample_file_result_tsv"] = get_link_for_dropdown_item(
                project.get_global_file_by_project(
                    TypePath.MEDIA_URL, Project.PROJECT_FILE_NAME_SAMPLE_RESULT_TSV
                )
            )
        if os.path.exists(
            project.get_global_file_by_project(
                TypePath.MEDIA_ROOT,
                Project.PROJECT_FILE_NAME_SAMPLE_RESULT_SETTINGS_CSV,
            )
        ):
            context[
                "samples_file_settings_statistics_csv"
            ] = get_link_for_dropdown_item(
                project.get_global_file_by_project(
                    TypePath.MEDIA_URL,
                    Project.PROJECT_FILE_NAME_SAMPLE_RESULT_SETTINGS_CSV,
                )
            )
        if os.path.exists(
            project.get_global_file_by_project(
                TypePath.MEDIA_ROOT,
                Project.PROJECT_FILE_NAME_SAMPLE_RESULT_SETTINGS_TSV,
            )
        ):
            context[
                "samples_file_settings_statistics_tsv"
            ] = get_link_for_dropdown_item(
                project.get_global_file_by_project(
                    TypePath.MEDIA_URL,
                    Project.PROJECT_FILE_NAME_SAMPLE_RESULT_SETTINGS_TSV,
                )
            )

        if os.path.exists(
            project.get_global_file_by_project(
                TypePath.MEDIA_ROOT,
                Project.PROJECT_FILE_NAME_SAMPLE_RESULT_all_consensus,
            )
        ):
            context["sample_file_all_consensus"] = get_link_for_dropdown_item(
                project.get_global_file_by_project(
                    TypePath.MEDIA_URL,
                    Project.PROJECT_FILE_NAME_SAMPLE_RESULT_all_consensus,
                )
            )

        ### need to test because in the past this file was not created
        context["freebays_variations_50_file"] = project.get_global_file_by_project_web(
            Project.PROJECT_FILE_NAME_TAB_VARIATIONS_FREEBAYES
        )
        file_name = project.get_global_file_by_project(
            TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_TOTAL_VARIATIONS
        )
        if not os.path.exists(file_name):
            collect_extra_data = CollectExtraData()
            collect_extra_data.calculate_count_variations(project)
        if os.path.exists(file_name):
            context[
                "variations_statistics_file"
            ] = project.get_global_file_by_project_web(
                Project.PROJECT_FILE_NAME_TOTAL_VARIATIONS
            )

        utils = Utils()
        vect_elements = utils.get_elements_from_db(project.reference, self.request.user)
        if vect_elements != None and len(vect_elements) > 0:
            context["elements"] = vect_elements
        vect_elements_protein = utils.get_elements_with_CDS_from_db(
            project.reference, self.request.user
        )
        if vect_elements_protein != None and len(vect_elements_protein) > 0:
            context[
                "elements_protein"
            ] = vect_elements_protein  ## because some does not have CDS
            ### get a vect of genes name
            context["genes"] = utils.get_vect_cds_from_element_from_db(
                vect_elements_protein[0], project.reference, self.request.user
            )

        ### pangolin data
        context["update_pangolin"] = False
        file_pangolin_result = project.get_global_file_by_project(
            TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_Pangolin_lineage
        )
        ## first condition is to the ones without pangolin lineage
        ### need at least one sequence fasta to run pangolin
        specie_tag = software.get_species_tag(project.reference)
        if (
            project.number_passed_sequences != 0
            and (
                not os.path.exists(file_pangolin_result)
                and specie_tag == Reference.SPECIES_SARS_COV_2
            )
            or (
                os.path.exists(file_pangolin_result)
                and software_pangolin.pangolin_results_out_date(project)
            )
        ):
            context["update_pangolin"] = True
            context["update_pangolin_message"] = mark_safe(
                software_pangolin.get_update_message(project)
            )

        ## pangolin file
        if project.number_passed_sequences > 0 and os.path.exists(file_pangolin_result):
            context["pangolin_lineage"] = get_link_for_dropdown_item(
                project.get_global_file_by_project(
                    TypePath.MEDIA_URL, Project.PROJECT_FILE_NAME_Pangolin_lineage
                )
            )

        #### nextclade link
        if (
            os.path.exists(
                project.get_global_file_by_project(
                    TypePath.MEDIA_ROOT,
                    Project.PROJECT_FILE_NAME_SAMPLE_RESULT_all_consensus,
                )
            )
            and settings.SHOW_NEXTCLADE_LINK
        ):  ## docker versions doesn't show NextClade link
            context = get_constext_nextclade(
                project.get_global_file_by_project(
                    TypePath.MEDIA_URL,
                    Project.PROJECT_FILE_NAME_SAMPLE_RESULT_all_consensus,
                ),
                context,
                get_current_site(self.request),
                specie_tag,
            )

        return context


class ProjectsSettingsView(LoginRequiredMixin, ListView):
    """
    can change settings in the projects
    """

    model = Project
    template_name = "settings/settings.html"
    context_object_name = "project"

    def get_context_data(self, **kwargs):

        context = super(ProjectsSettingsView, self).get_context_data(**kwargs)
        project = Project.objects.get(pk=self.kwargs["pk"])

        ### can't see this project
        context["nav_project"] = True
        if project.owner.id != self.request.user.id:
            context["error_cant_see"] = "1"
            return context

        ### test all defaults first, if exist in database
        default_software = DefaultProjectSoftware()
        default_software.test_all_defaults(
            self.request.user, project, None, None, None
        )  ## the user can have defaults yet

        all_tables = []  ## order by Technology, PipelineStep, table
        ## [ [unique_id, Technology, [ [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], ...],
        ##    [unique_id, Technology, [ [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], ...], etc
        ## Technology goes to NAV-container, PipelineStep goes to NAV-container, then table
        ## Mix parameters with software

        ### count samples already assigned to Project
        count_project_sample = ProjectSample.objects.filter(
            project=project, is_deleted=False
        ).count()
        ### IMPORTANT, must have technology__name, because old versions don't
        for technology in ConstantsSettings.vect_technology:  ## run over all technology
            vect_pipeline_step = []
            for pipeline_step in ConstantsSettings.vect_pipeline_names:
                query_set = SoftwareSettings.objects.filter(
                    owner=self.request.user,
                    type_of_use=SoftwareSettings.TYPE_OF_USE_project,
                    parameter__project=project,
                    parameter__project_sample=None,
                    type_of_software__in=[
                        SoftwareSettings.TYPE_SOFTWARE,
                        SoftwareSettings.TYPE_INSAFLU_PARAMETER,
                    ],
                    technology__name=technology,
                    pipeline_step__name=pipeline_step,
                    is_obsolete=False,
                ).distinct()

                ### if there are software
                if query_set.count() > 0:
                    vect_pipeline_step.append(
                        [
                            "{}_{}".format(
                                pipeline_step.replace(" ", "").replace("/", ""),
                                technology.replace(" ", "").replace("/", ""),
                            ),
                            pipeline_step,
                            SoftwaresTable(
                                query_set,
                                project=project,
                                project_sample=None,
                                sample=None,
                                b_enable_options=count_project_sample == 0,
                            ),
                        ]
                    )

            ## if there is software for the pipeline step
            if len(vect_pipeline_step) > 0:
                all_tables.append(
                    [
                        technology.replace(" ", "").replace("/", ""),
                        technology,
                        vect_pipeline_step,
                    ]
                )

        context["all_softwares"] = all_tables
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        context["project"] = project
        context["main_settings"] = False
        context["project_settings"] = True
        return context


class SampleProjectsSettingsView(LoginRequiredMixin, ListView):
    """
    can change settings in the projects
    """

    model = ProjectSample
    template_name = "settings/settings.html"
    context_object_name = "project_sample"

    def get_context_data(self, **kwargs):

        context = super(SampleProjectsSettingsView, self).get_context_data(**kwargs)
        project_sample = ProjectSample.objects.get(pk=self.kwargs["pk"])

        ### can't see this project
        context["nav_project"] = True
        if project_sample.project.owner.id != self.request.user.id:
            context["error_cant_see"] = "1"
            return context

        ### test all defaults first, if exist in database
        default_software = DefaultProjectSoftware()
        default_software.test_all_defaults(
            self.request.user, None, project_sample, None
        )  ## the user can have defaults yet

        all_tables = []  ## order by Technology, PipelineStep, table
        ## [ [unique_id, Technology, [ [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], ...],
        ##    [unique_id, Technology, [ [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], ...], etc
        ## Technology goes to NAV-container, PipelineStep goes to NAV-container, then table
        ## Mix parameters with software

        ### IMPORTANT, must have technology__name, because old versions don't
        for technology in ConstantsSettings.vect_technology:  ## run over all technology
            vect_pipeline_step = []
            for pipeline_step in ConstantsSettings.vect_pipeline_names:
                query_set = SoftwareSettings.objects.filter(
                    owner=self.request.user,
                    type_of_use=SoftwareSettings.TYPE_OF_USE_project_sample,
                    parameter__project=None,
                    parameter__project_sample=project_sample,
                    type_of_software__in=[
                        SoftwareSettings.TYPE_SOFTWARE,
                        SoftwareSettings.TYPE_INSAFLU_PARAMETER,
                    ],
                    technology__name=technology,
                    pipeline_step__name=pipeline_step,
                    is_obsolete=False,
                ).distinct()

                ### if there are software
                if query_set.count() > 0:
                    vect_pipeline_step.append(
                        [
                            "{}_{}".format(
                                pipeline_step.replace(" ", "").replace("/", ""),
                                technology.replace(" ", "").replace("/", ""),
                            ),
                            pipeline_step,
                            SoftwaresTable(
                                query_set,
                                project=None,
                                project_sample=project_sample,
                                sample=None,
                            ),
                        ]
                    )
            ## if there is software for the pipeline step
            if len(vect_pipeline_step) > 0:
                all_tables.append(
                    [
                        technology.replace(" ", "").replace("/", ""),
                        technology,
                        vect_pipeline_step,
                    ]
                )

        context["all_softwares"] = all_tables
        context["project_sample"] = project_sample
        context["project_sample_settings"] = True
        context["main_settings"] = False
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        return context


class SampleSettingsView(LoginRequiredMixin, ListView):
    """
    can change settings in the projects
    """

    model = Sample
    template_name = "settings/settings.html"
    context_object_name = "sample"

    def get_context_data(self, **kwargs):
        context = super(SampleSettingsView, self).get_context_data(**kwargs)
        sample = Sample.objects.get(pk=self.kwargs["pk"])

        ### can't see this sample
        context["nav_sample"] = True
        if sample.owner.id != self.request.user.id:
            context["error_cant_see"] = "1"
            return context

        ### test all defaults first, if exist in database
        default_software = DefaultProjectSoftware()
        default_software.test_all_defaults(
            user=self.request.user,
            project=None,
            project_sample=None,
            sample=sample,
            dataset=None,
        )  ## the user can have defaults yet

        all_tables = []  ## order by Technology, PipelineStep, table
        ## [ [unique_id, Technology, [ [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], ...],
        ##    [unique_id, Technology, [ [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], [unique_id, PipelineStep, table], ...], etc
        ## Technology goes to NAV-container, PipelineStep goes to NAV-container, then table
        ## Mix parameters with software

        ### IMPORTANT, must have technology__name, because old versions don't
        for technology in ConstantsSettings.vect_technology:  ## run over all technology
            vect_pipeline_step = []
            for pipeline_step in ConstantsSettings.vect_pipeline_names:
                query_set = SoftwareSettings.objects.filter(
                    owner=self.request.user,
                    type_of_use=SoftwareSettings.TYPE_OF_USE_sample,
                    parameter__sample=sample,
                    type_of_software__in=[SoftwareSettings.TYPE_SOFTWARE],
                    technology__name=technology,
                    pipeline_step__name=pipeline_step,
                    is_obsolete=False,
                ).distinct()

                ### if there are software
                if query_set.count() > 0:
                    vect_pipeline_step.append(
                        [
                            "{}_{}".format(
                                pipeline_step.replace(" ", "").replace("/", ""),
                                technology.replace(" ", "").replace("/", ""),
                            ),
                            pipeline_step,
                            SoftwaresTable(
                                query_set,
                                project=None,
                                project_sample=None,
                                sample=sample,
                            ),
                        ]
                    )

            ## if there is software for the pipeline step
            if len(vect_pipeline_step) > 0:
                all_tables.append(
                    [
                        technology.replace(" ", "").replace("/", ""),
                        technology,
                        vect_pipeline_step,
                    ]
                )

        context["all_softwares"] = all_tables
        context["sample"] = sample
        context["sample_settings"] = True
        context["main_settings"] = False
        context[
            "show_info_main_page"
        ] = ShowInfoMainPage()  ## show main information about the institute
        return context


class ShowSampleProjectsDetailsView(LoginRequiredMixin, ListView):
    """ """

    utils = Utils()
    software = Software()
    model = ProjectSample
    template_name = "project_sample/project_sample_single_detail.html"
    context_object_name = "project_sample"

    def get_context_data(self, **kwargs):
        context = super(ShowSampleProjectsDetailsView, self).get_context_data(**kwargs)
        try:
            ### can't see this project
            project_sample = ProjectSample.objects.get(pk=self.kwargs["pk"])
            context["nav_project"] = True
            if project_sample.project.owner.id != self.request.user.id:
                context["error_cant_see"] = "1"
                return context

            ## several data to show
            context["project_sample"] = project_sample
            context["num_alerts"] = (
                project_sample.alert_first_level + project_sample.alert_second_level
            )
            context["nav_project"] = True
            context[
                "show_info_main_page"
            ] = ShowInfoMainPage()  ## show main information about the institute

            ## collect alerts
            alert_out = []
            manageDatabase = ManageDatabase()
            metaKeyAndValue = MetaKeyAndValue()
            vect_elements = self.utils.get_elements_from_db(
                project_sample.project.reference, project_sample.project.owner
            )
            for element_temp in vect_elements:
                for key_message in metaKeyAndValue.VECT_MESSAGE_ALERT_COVERAGE:
                    meta_key = metaKeyAndValue.get_meta_key(key_message, element_temp)
                    meta_data = manageDatabase.get_project_sample_metakey_last(
                        project_sample, meta_key, MetaKeyAndValue.META_VALUE_Success
                    )
                    if meta_data != None:
                        alert_out.append(meta_data.description)

            ### get different types of alerts
            for (
                key
            ) in metaKeyAndValue.get_keys_show_alerts_in_sample_projects_details_view():
                meta_data = manageDatabase.get_project_sample_metakey_last(
                    project_sample, key, MetaKeyAndValue.META_VALUE_Success
                )
                if meta_data != None:
                    alert_out.append(meta_data.description)
            context["alerts"] = alert_out

            ##### extra data sample, columns added by the user
            ## [[header1, value1], [header2, value2], [header3, value3], ...]
            ### if it's to big expand button is better
            tag_names = project_sample.sample.get_tag_names()
            context["extra_data_sample_expand"] = (
                tag_names != None
            )  ## (tag_names != None and tag_names.count()  > (Constants.START_EXPAND_SAMPLE_TAG_NAMES_ROWS))
            if tag_names != None:
                context["extra_data_sample"] = self.utils.grouped(tag_names, 4)

            default_software = DefaultProjectSoftware()
            context["consensus_file"] = project_sample.get_consensus_file_web(
                not default_software.include_consensus(project_sample)
            )
            software_used = (
                []
            )  ### has a list with all software used... [name, parameters]
            ### only for illumina
            decode_result = DecodeObjects()
            if project_sample.is_sample_illumina():
                context["snippy_variants_file"] = project_sample.get_file_web(
                    FileType.FILE_TAB, SoftwareNames.SOFTWARE_SNIPPY_name
                )
                context["freebayes_variants_file"] = project_sample.get_file_web(
                    FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name
                )
                b_second_choice = True
                context[
                    "freebayes_variants_file_snp_indel"
                ] = project_sample.get_file_web(
                    FileType.FILE_TAB,
                    SoftwareNames.SOFTWARE_FREEBAYES_name,
                    b_second_choice,
                )
                context["depth_file"] = project_sample.get_file_web(
                    FileType.FILE_DEPTH_GZ, SoftwareNames.SOFTWARE_SNIPPY_name
                )
                context["depth_tbi_file"] = project_sample.get_file_web(
                    FileType.FILE_DEPTH_GZ_TBI, SoftwareNames.SOFTWARE_SNIPPY_name
                )

                #### software versions...
                software_used.append([SoftwareNames.SOFTWARE_SNIPPY_name, "Fail"])
                software_used.append([SoftwareNames.SOFTWARE_FREEBAYES_name, "Fail"])
                list_meta = manageDatabase.get_project_sample_metakey(
                    project_sample, MetaKeyAndValue.META_KEY_Snippy_Freebayes, None
                )
                if (
                    list_meta[0].value == MetaKeyAndValue.META_VALUE_Success
                    and MetaKeyAndValue.META_KEY_Snippy_Freebayes
                    == list_meta[0].meta_tag.name
                ):
                    result = decode_result.decode_result(list_meta[0].description)
                    if not result is None:
                        software_used[0][1] = result.get_software(
                            SoftwareNames.SOFTWARE_SNIPPY_name
                        )
                        software_used[1][1] = result.get_software(
                            SoftwareNames.SOFTWARE_FREEBAYES_name
                        )

                        ### could have or not, in older versions
                        msa_markers_software = result.get_software(
                            SoftwareNames.SOFTWARE_MSA_MASKER_name
                        )
                        if len(msa_markers_software) > 0:
                            software_used.append(
                                [
                                    SoftwareNames.SOFTWARE_MSA_MASKER_name,
                                    result.get_software(
                                        SoftwareNames.SOFTWARE_MSA_MASKER_name
                                    ),
                                ]
                            )

            else:
                context["medaka_variants_file"] = project_sample.get_file_web(
                    FileType.FILE_TAB, SoftwareNames.SOFTWARE_Medaka_name
                )
                context["depth_file"] = project_sample.get_file_web(
                    FileType.FILE_DEPTH_GZ, SoftwareNames.SOFTWARE_Medaka_name
                )
                context["depth_tbi_file"] = project_sample.get_file_web(
                    FileType.FILE_DEPTH_GZ_TBI, SoftwareNames.SOFTWARE_Medaka_name
                )

                #### software versions...
                software_used.append([SoftwareNames.SOFTWARE_Medaka_name, "Fail"])
                list_meta = manageDatabase.get_project_sample_metakey(
                    project_sample, MetaKeyAndValue.META_KEY_Medaka, None
                )
                if (
                    list_meta[0].value == MetaKeyAndValue.META_VALUE_Success
                    and MetaKeyAndValue.META_KEY_Medaka == list_meta[0].meta_tag.name
                ):
                    result = decode_result.decode_result(list_meta[0].description)
                    if not result is None:
                        software_used[0][1] = result.get_software(
                            SoftwareNames.SOFTWARE_Medaka_name
                        )

                        software_software = result.get_software(
                            SoftwareNames.SOFTWARE_SAMTOOLS_name
                        )
                        if len(software_software) > 0:
                            software_used.append(
                                [
                                    SoftwareNames.SOFTWARE_SAMTOOLS_name,
                                    result.get_software(
                                        SoftwareNames.SOFTWARE_SAMTOOLS_name
                                    ),
                                ]
                            )

                        ### could have or not, in older versions
                        msa_markers_software = result.get_software(
                            SoftwareNames.SOFTWARE_MSA_MASKER_name
                        )
                        if len(msa_markers_software) > 0:
                            software_used.append(
                                [
                                    SoftwareNames.SOFTWARE_MSA_MASKER_name,
                                    result.get_software(
                                        SoftwareNames.SOFTWARE_MSA_MASKER_name
                                    ),
                                ]
                            )

                        software_software = result.get_software(
                            SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name
                        )
                        if len(software_software) > 0:
                            software_used.append(
                                [
                                    SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name_extended,
                                    result.get_software(
                                        SoftwareNames.INSAFLU_PARAMETER_LIMIT_COVERAGE_ONT_name
                                    ),
                                ]
                            )

                        software_software = result.get_software(
                            SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name
                        )
                        if len(software_software) > 0:
                            software_used.append(
                                [
                                    SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name_extended,
                                    result.get_software(
                                        SoftwareNames.INSAFLU_PARAMETER_VCF_FREQ_ONT_name
                                    ),
                                ]
                            )

                        software_software = result.get_software(
                            SoftwareNames.SOFTWARE_BCFTOOLS_name
                        )
                        if len(software_software) > 0:
                            software_used.append(
                                [
                                    SoftwareNames.SOFTWARE_BCFTOOLS_name,
                                    result.get_software(
                                        SoftwareNames.SOFTWARE_BCFTOOLS_name
                                    ),
                                ]
                            )

            ## pangolin version, it is transversel (illumina e ONT)
            meta_key = manageDatabase.get_project_metakey_last(
                project_sample.project,
                MetaKeyAndValue.META_KEY_Identify_pangolin,
                MetaKeyAndValue.META_VALUE_Success,
            )
            if not meta_key is None:
                result = decode_result.decode_result(meta_key.description)
                ## only need to call first Pangolin, PangolinLearn is added automatically
                software_used.append(
                    [
                        SoftwareNames.SOFTWARE_Pangolin_name,
                        result.get_software(
                            SoftwareNames.SOFTWARE_Pangolin_name_search_name
                        ),
                    ]
                )

            ### list of software to used
            context["software_used"] = software_used

            #### nextclade link
            # is_sars_cov_2 = software_pangolin.is_ref_sars_cov_2(project_sample.project.reference.get_reference_fasta(TypePath.MEDIA_ROOT))
            specie_tag = self.software.get_species_tag(project_sample.project.reference)
            if (
                os.path.exists(project_sample.get_consensus_file(TypePath.MEDIA_ROOT))
                and settings.SHOW_NEXTCLADE_LINK
                and default_software.include_consensus(project_sample)
            ):  ## docker versions doesn't show NextClade link
                context = get_constext_nextclade(
                    project_sample.get_consensus_file(TypePath.MEDIA_URL),
                    context,
                    get_current_site(self.request),
                    specie_tag,
                )

        except ProjectSample.DoesNotExist:
            context["error_cant_see"] = 1
        return context


def is_all_check_box_in_session(vect_check_to_test, request):
    """
    test if all check boxes are in session
    If not remove the ones that are in and create the new ones all False
    """
    utils = Utils()
    dt_data = {}

    ## get the dictonary
    for key in request.session.keys():
        if (
            key.startswith(Constants.CHECK_BOX)
            and len(key.split("_")) == 3
            and utils.is_integer(key.split("_")[2])
        ):
            dt_data[key] = True

    b_different = False
    if len(vect_check_to_test) != len(dt_data):
        b_different = True

    ## test the vector
    if not b_different:
        for key in vect_check_to_test:
            if key not in dt_data:
                b_different = True
                break

    if b_different:
        ## remove all
        for key in dt_data:
            del request.session[key]

        ## create new
        for key in vect_check_to_test:
            request.session[key] = False
        return False
    return True


def clean_check_box_in_session(request):
    """
    check all check boxes on samples/references to add samples
    """
    utils = Utils()
    ## clean all check unique
    if Constants.CHECK_BOX_ALL in request.session:
        del request.session[Constants.CHECK_BOX_ALL]
    vect_keys_to_remove = []
    for key in request.session.keys():
        if (
            key.startswith(Constants.CHECK_BOX)
            and len(key.split("_")) == 3
            and utils.is_integer(key.split("_")[2])
        ):
            vect_keys_to_remove.append(key)
    for key in vect_keys_to_remove:
        del request.session[key]


def get_first_pk_from_session(request):
    """
    return the first pk selected, if none, return none
    """
    utils = Utils()
    for key in request.session.keys():
        if (
            key.startswith(Constants.CHECK_BOX)
            and len(key.split("_")) == 3
            and utils.is_integer(key.split("_")[2])
            and request.session[key]
        ):
            return key.split("_")[2]
    return None


def get_unique_pk_from_session(request):
    """
    return the unique pk selected. If exists more than one return none
    """
    utils = Utils()
    return_pk = None
    for key in request.session.keys():
        if (
            key.startswith(Constants.CHECK_BOX)
            and len(key.split("_")) == 3
            and utils.is_integer(key.split("_")[2])
            and request.session[key]
        ):
            if return_pk == None:
                return_pk = key.split("_")[2]
            else:
                return None
    return return_pk
