import django_tables2 as tables
from constants.constants import Constants, TypePath
from constants.meta_key_and_values import MetaKeyAndValue
from django.conf import settings
from django.db.models import F
from django.urls import reverse
from django.utils.safestring import mark_safe
from django.utils.translation import ugettext_lazy as _
from settings.constants_settings import ConstantsSettings
from settings.default_parameters import DefaultParameters
from settings.default_software_project_sample import DefaultProjectSoftware
from utils.result import DecodeObjects

from managing_files.manage_database import ManageDatabase
from managing_files.models import Project, ProjectSample, Reference, Sample


class CheckBoxColumnWithName(tables.CheckBoxColumn):
    @property
    def header(self):
        default = {"type": "checkbox"}
        return self.verbose_name


class ReferenceTable(tables.Table):
    #   Renders a normal value as an internal hyperlink to another page.
    #   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
    reference_fasta_name = tables.LinkColumn(
        "reference_fasta_name", args=[tables.A("pk")], verbose_name="Fasta file"
    )
    reference_genbank_name = tables.LinkColumn(
        "reference_genbank_name", args=[tables.A("pk")], verbose_name="GenBank file"
    )
    owner = tables.Column("Owner", orderable=True, empty_values=())
    constants = Constants()

    SHORT_NAME_LENGTH = 20

    class Meta:
        model = Reference
        fields = (
            "name",
            "reference_fasta_name",
            "reference_genbank_name",
            "creation_date",
            "owner",
            "number_of_locus",
        )
        attrs = {"class": "table-striped table-bordered"}
        empty_text = "There are no References to show..."

    def render_owner(self, record):
        return record.owner.username

    def render_name(self, record):
        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()
        user = current_request.user
        if user.username == Constants.USER_ANONYMOUS:
            return record.name
        if (
            user.username == record.owner.username
            and record.project.all().filter(is_deleted=False).count() == 0
        ):  ## it can't be in any active project
            return mark_safe(
                '<a href="#modal_remove_reference" id="id_remove_reference_modal" data-toggle="modal"'
                + ' ref_name="'
                + record.name
                + '" pk="'
                + str(record.pk)
                + '"><i class="fa fa-trash"></i></span> </a>'
                + record.name
            )
        return record.name

    def render_reference_fasta_name(self, **kwargs):
        record = kwargs.pop("record")
        return record.get_reference_fasta_web()

    def render_reference_genbank_name(self, **kwargs):
        record = kwargs.pop("record")
        return record.get_reference_gb_web()

    def render_creation_date(self, **kwargs):
        record = kwargs.pop("record")
        return record.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)

    def order_owner(self, queryset, is_descending):
        queryset = queryset.annotate(owner_name=F("owner__username")).order_by(
            ("-" if is_descending else "") + "owner_name"
        )
        return (queryset, True)


class ReferenceProjectTable(tables.Table):
    """
    create a new project with a reference
    """

    select_ref = CheckBoxColumnWithName(
        verbose_name=("Select One"), accessor="pk", orderable=False
    )
    owner = tables.Column("Owner", orderable=True, empty_values=())
    #     number_of_locus = tables.Column(orderable=False)

    class Meta:
        model = Reference
        fields = (
            "select_ref",
            "name",
            "isolate_name",
            "creation_date",
            "owner",
            "number_of_locus",
        )
        attrs = {"class": "table-striped table-bordered"}
        empty_text = "There are no References to add..."

    def render_creation_date(self, **kwargs):
        record = kwargs.pop("record")
        return record.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)

    def render_select_ref(self, value, record):
        return mark_safe(
            '<input name="select_ref" id="{}_{}" type="checkbox" value="{}"/>'.format(
                Constants.CHECK_BOX, record.id, value
            )
        )

    def order_owner(self, queryset, is_descending):
        queryset = queryset.annotate(owner_name=F("owner__username")).order_by(
            ("-" if is_descending else "") + "owner_name"
        )
        return (queryset, True)


class SampleToProjectsTable(tables.Table):
    """
    To add samples to projects
    """

    type_subtype = tables.Column("Type-Subtype", orderable=True, empty_values=())
    select_ref = tables.CheckBoxColumn(
        accessor="pk", attrs={"th__input": {"id": "checkBoxAll"}}, orderable=False
    )

    class Meta:
        model = Sample
        fields = (
            "select_ref",
            "name",
            "creation_date",
            "type_subtype",
            "week",
            "data_set",
        )
        attrs = {"class": "table-striped table-bordered"}
        empty_text = "There are no Samples to add..."

    def render_select_ref(self, value, record):
        return mark_safe(
            '<input name="select_ref" id="{}_{}" type="checkbox" value="{}"/>'.format(
                Constants.CHECK_BOX, record.id, value
            )
        )

    def render_creation_date(self, **kwargs):
        record = kwargs.pop("record")
        return record.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)


class SampleTable(tables.Table):

    ### manage database
    manage_database = ManageDatabase()

    #   Renders a normal value as an internal hyperlink to another page.
    #   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
    number_quality_sequences = tables.Column(
        "#Quality Seq. (Fastq1)-(Fastq2)",
        footer="#original number sequence",
        orderable=False,
        empty_values=(),
    )
    #     extra_info = tables.LinkColumn('sample-description', args=[tables.A('pk')], orderable=False, verbose_name='Extra Information', empty_values=())
    extra_info = tables.LinkColumn(
        "Extra Information", orderable=False, empty_values=()
    )
    technology = tables.Column("Technology", empty_values=())
    type_and_subtype = tables.Column("Classification", empty_values=())
    fastq_files = tables.Column("#Fastq Files", empty_values=())
    data_set = tables.Column("Data set", empty_values=())

    class Meta:
        model = Sample
        fields = (
            "name",
            "creation_date",
            "fastq_files",
            "technology",
            "type_and_subtype",
            "data_set",
            "number_alerts",
            "number_quality_sequences",
            "extra_info",
        )
        attrs = {"class": "table-striped table-bordered"}
        empty_text = "There are no Samples to show..."
        tooltips = (
            "Name for everyone",
            "creation_date",
            "fastq_files",
            "type_and_subtype",
            "data_set",
            "number_quality_sequences",
            "extra_info",
        )

    def render_fastq_files(self, record):
        """
        number of fastqFiles, to show if there is
        """
        if record.has_files:
            if record.file_name_2 is None or len(record.file_name_2) == 0:
                return "1"
            return "2"
        return "0"

    def render_name(self, record):
        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()
        user = current_request.user

        ### get the link for sample, to expand data
        manageDatabase = ManageDatabase()
        sample_name = record.name
        if not manageDatabase.is_sample_processing_step(record):
            list_meta = manageDatabase.get_sample_metakey(
                record,
                MetaKeyAndValue.META_KEY_Fastq_Trimmomatic
                if record.is_type_fastq_gz_sequencing()
                else MetaKeyAndValue.META_KEY_NanoStat_NanoFilt,
                None,
            )
            if (
                record.is_ready_for_projects
                and list_meta.count() > 0
                and list_meta[0].value == MetaKeyAndValue.META_VALUE_Success
            ):
                sample_name = (
                    "<a href="
                    + reverse("sample-description", args=[record.pk])
                    + "> {}</a>".format(record.name)
                )
            ### read not attached, yet
            elif (
                len(list_meta) == 0
                and not record.candidate_file_name_1 is None
                and len(record.candidate_file_name_1) > 0
            ):
                sample_name = (
                    "<a href="
                    + reverse("sample-description", args=[record.pk])
                    + "> {}</a>".format(record.name)
                )

        if user.username == Constants.USER_ANONYMOUS:
            return mark_safe(sample_name)
        if user.username == record.owner.username:  ## it can't be in any active project
            ## it can have project samples not deleted but in projects deleted
            for project_samples in record.project_samples.all().filter(
                is_deleted=False, is_error=False
            ):
                if (
                    not project_samples.is_deleted
                    and not project_samples.project.is_deleted
                ):
                    return mark_safe(sample_name)

            return mark_safe(
                '<a href="#id_remove_modal" id="id_remove_reference_modal" data-toggle="modal"'
                + ' ref_name="'
                + record.name
                + '" pk="'
                + str(record.pk)
                + '"><i class="fa fa-trash"></i></span> </a>'
                + sample_name
            )
        return mark_safe(sample_name)

    def render_technology(self, record):
        """shows if it is Illumina or Minion"""
        ### is not processed yet
        if record.type_of_fastq == Sample.TYPE_OF_FASTQ_not_defined:
            return "Undefined"
        return (
            ConstantsSettings.TECHNOLOGY_illumina
            if record.is_type_fastq_gz_sequencing()
            else ConstantsSettings.TECHNOLOGY_minion
        )

    def render_creation_date(self, **kwargs):
        record = kwargs.pop("record")
        return record.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)

    def render_type_and_subtype(self, record):
        """
        get type and sub type
        """
        return _("Not yet") if record.type_subtype is None else record.type_subtype

    def render_data_set(self, record):
        """
        return name of dataset
        """
        return record.data_set.name

    def render_number_quality_sequences(self, record):
        """
        number of quality sequences and average
        """
        manageDatabase = ManageDatabase()

        ## test ready to show
        if not record.get_is_ready_for_projects():
            return _("Not yet")

        ## get all stats
        list_meta = manageDatabase.get_sample_metakey(
            record, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None
        )
        if (
            not list_meta is None
            and list_meta.count() > 0
            and list_meta[0].value == MetaKeyAndValue.META_VALUE_Success
        ):
            decodeResultAverageAndNumberReads = DecodeObjects()
            result_average = decodeResultAverageAndNumberReads.decode_result(
                list_meta[0].description
            )
            if result_average.number_file_2 is None:
                if record.is_type_fastq_gz_sequencing():  ### default illumina
                    tip_info = (
                        '<span ><i class="tip fa fa-info-circle" title="#Reads from fastq1: {}\n'
                        'Length average from fastq1: {}\n"></i></span>'.format(
                            result_average.number_file_1, result_average.average_file_1
                        )
                    )
                else:
                    tip_info = (
                        '<span ><i class="tip fa fa-info-circle" title="#Reads from fastq: {}\n'
                        'Length average from fastq: {}\n"></i></span>'.format(
                            result_average.number_file_1, result_average.average_file_1
                        )
                    )
                return mark_safe(
                    tip_info
                    + " (%s/%s) (%s)"
                    % (
                        result_average.number_file_1,
                        result_average.average_file_1,
                        result_average.number_file_1,
                    )
                )

            tip_info = (
                '<span ><i class="tip fa fa-info-circle" title="#Reads from fastq1: {}\n'
                "Length average from fastq1: {}\n\n#Reads from fastq2: {}\n"
                'Length average from fastq2: {}\n\n#Total reads: {}"></i></span>'.format(
                    result_average.number_file_1,
                    result_average.average_file_1,
                    result_average.number_file_2,
                    result_average.average_file_2,
                    int(result_average.number_file_1)
                    + int(result_average.number_file_2),
                )
            )
            return mark_safe(
                tip_info
                + " (%s/%s)-(%s/%s) (%d)"
                % (
                    result_average.number_file_1,
                    result_average.average_file_1,
                    result_average.number_file_2,
                    result_average.average_file_2,
                    int(result_average.number_file_1)
                    + int(result_average.number_file_2),
                )
            )
        elif (
            not list_meta is None
            and list_meta.count() > 0
            and list_meta[0].value.equals(MetaKeyAndValue.META_VALUE_Error)
        ):
            return _("Error")
        return _("Not yet")

    def render_extra_info(self, record):
        """
        icon with link to extra info
        """
        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()
        user = current_request.user

        manageDatabase = ManageDatabase()
        if manageDatabase.is_sample_wating_fastq_file(record):
            return _("Waiting FASTQ files")
        if manageDatabase.is_sample_processing_step(record):
            return _("Processing")

        list_meta = manageDatabase.get_sample_metakey(
            record,
            MetaKeyAndValue.META_KEY_Fastq_Trimmomatic
            if record.is_type_fastq_gz_sequencing()
            else MetaKeyAndValue.META_KEY_NanoStat_NanoFilt,
            None,
        )

        if (
            record.is_ready_for_projects
            and list_meta.count() > 0
            and list_meta[0].value == MetaKeyAndValue.META_VALUE_Success
        ):
            if user.username != Constants.USER_ANONYMOUS:
                ## test if it has the original fastq files
                str_links = self._get_magic_handle(record)
            else:
                str_links = ""
            return mark_safe(
                str_links
                + "<a href="
                + reverse("sample-description", args=[record.pk])
                + '><span ><i class="fa fa-plus-square"></i></span> More Info</a>'
            )
        ### read not attached, yet
        elif (
            len(list_meta) == 0
            and not record.candidate_file_name_1 is None
            and len(record.candidate_file_name_1) > 0
        ):
            return mark_safe(
                "<a href="
                + reverse("sample-description", args=[record.pk])
                + '><span ><i class="fa fa-plus-square"></i></span> More Info</a>'
            )
        elif (
            list_meta.count() > 0
            and list_meta[0].value == MetaKeyAndValue.META_VALUE_Error
        ):
            return _("Error")
        ### this case is when the number of remain reads are zero
        ### it is success but there's no reads
        elif (
            list_meta.count() > 0
            and list_meta[0].value == MetaKeyAndValue.META_VALUE_Success
        ):
            str_links = self._get_magic_handle(record)
            if len(str_links) > 0:
                return mark_safe(str_links + " No reads left")
        return _("Not yet")

    def _get_magic_handle(self, record):
        """ """
        str_links = ""
        if record.is_original_fastq_removed():
            str_links = (
                '<a href=# data-toggle="tooltip" title="Software settings are not enable because original fastq files were removed">'
                + '<span ><i class="padding-button-table fa fa-magic padding-button-table" style="color: grey"></i></span></a>'
            )
        ## test if it's in a project
        elif (
            not record.project_samples is None
            and record.project_samples.filter(is_deleted=False, is_error=False).count()
            > 0
        ):
            str_links = (
                '<a href=# data-toggle="tooltip" title="Software settings are not enable because this sample is in at least one project.">'
                + '<span ><i class="padding-button-table fa fa-magic padding-button-table" style="color: grey"></i></span></a>'
            )
        else:
            str_links = (
                "<a href="
                + reverse("sample-settings", args=[record.pk])
                + ' data-toggle="tooltip" title="Software settings">'
                + '<span ><i class="padding-button-table fa fa-magic padding-button-table"></i></span></a>'
            )

        ### check if it was downsized
        if self.manage_database.is_sample_downsized(record):
            str_links += (
                '<a href="javascript:void(0);" data-toggle="tooltip" title="Sample was downsized">'
                + '<span ><i class="padding-button-table warning_fa_icon fa fa-level-down"></i></a>'
            )
        return str_links

    def order_fastq_files(self, queryset, is_descending):
        queryset = queryset.annotate(
            is_valid__1=F("is_valid_1"), is_valid__2=F("is_valid_2")
        ).order_by(
            ("-" if is_descending else "") + "is_valid__1",
            ("-" if is_descending else "") + "is_valid__2",
        )
        return (queryset, True)

    def order_type_and_subtype(self, queryset, is_descending):
        queryset = queryset.annotate(type_and_subtype=F("type_subtype")).order_by(
            ("-" if is_descending else "") + "type_subtype"
        )
        return (queryset, True)

    def order_alerts(self, queryset, is_descending):
        queryset = queryset.annotate(number_alerts_=F("number_alerts")).order_by(
            ("-" if is_descending else "") + "number_alerts_"
        )
        return (queryset, True)

    def order_technology(self, queryset, is_descending):
        """shows if it is Illumina or Minion"""
        queryset = queryset.annotate(technology=F("type_of_fastq")).order_by(
            ("-" if is_descending else "") + "technology"
        )
        return (queryset, True)


class ProjectTable(tables.Table):
    #   Renders a normal value as an internal hyperlink to another page.
    #   account_number = tables.LinkColumn('customer-detail', args=[A('pk')])
    reference = tables.Column("Reference", empty_values=())
    samples = tables.Column("#Samples (P/W/E)", orderable=False, empty_values=())
    last_change_date = tables.Column("Last Change date", empty_values=())
    creation_date = tables.Column("Creation date", empty_values=())
    results = tables.LinkColumn("Options", orderable=False, empty_values=())

    class Meta:
        model = Project
        fields = (
            "name",
            "reference",
            "last_change_date",
            "creation_date",
            "samples",
            "results",
        )
        attrs = {"class": "table-striped table-bordered"}
        empty_text = "There are no Projects to show..."

    def render_name(self, record):
        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()
        user = current_request.user

        ## there's nothing to show
        count = ProjectSample.objects.filter(
            project__id=record.id, is_deleted=False, is_error=False, is_finished=True
        ).count()
        project_sample = record.name
        if count > 0:
            project_sample = (
                "<a href="
                + reverse("show-sample-project-results", args=[record.pk])
                + ' data-toggle="tooltip" title="See Results">'
                + "{}</a>".format(record.name)
            )
        else:
            project_sample = record.name

        if user.username == Constants.USER_ANONYMOUS:
            return mark_safe(project_sample)
        if user.username == record.owner.username:
            return mark_safe(
                '<a href="#id_remove_modal" id="id_remove_reference_modal" data-toggle="modal" data-toggle="tooltip" title="Delete"'
                + ' ref_name="'
                + record.name
                + '" pk="'
                + str(record.pk)
                + '"><i class="fa fa-trash"></i></span> </a>'
                + project_sample
            )
        return project_sample

    def render_samples(self, record):
        """
        return a reference name
        """
        add_remove = ""
        # if (ProjectSample.objects.filter(project__id=record.id, is_deleted=False).count() > 0):
        #     TODO
        #     add_remove = ' <a href=' + reverse('remove-sample-project', args=[record.pk]) + '><span ><i class="fa fa-trash"></i></span> Remove</a>'
        #     add_remove = ' <a href="#"><span ><i class="fa fa-trash"></i></span> Remove</a>'

        n_processed = ProjectSample.objects.filter(
            project__id=record.id, is_deleted=False, is_error=False, is_finished=True
        ).count()
        n_error = ProjectSample.objects.filter(
            project__id=record.id, is_deleted=False, is_error=True, is_finished=False
        ).count()
        n_processing = ProjectSample.objects.filter(
            project__id=record.id, is_deleted=False, is_error=False, is_finished=False
        ).count()
        tip_info = '<span ><i class="tip fa fa-info-circle" title="Processed: {}\nWaiting: {}\nError: {}"></i></span>'.format(
            n_processed, n_processing, n_error
        )
        return mark_safe(
            tip_info
            + " ({}/{}/{}) ".format(n_processed, n_processing, n_error)
            + "<a href="
            + reverse("add-sample-project", args=[record.pk])
            + ' data-toggle="tooltip" title="Add samples" ><i class="fa fa-plus-square"></i> Add</a>'  #         return mark_safe(tip_info + " ({}/{}/{}) ".format(n_processed, n_processing, n_error) + '<a href=# id="id_add_sample_message"' +\
            + add_remove
        )

    def render_creation_date(self, **kwargs):
        record = kwargs.pop("record")
        return record.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)

    def render_last_change_date(self, **kwargs):
        record = kwargs.pop("record")
        if record.last_change_date is None:
            return "Not set yet"
        return record.last_change_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)

    def render_results(self, record):
        """
        icon with link to extra info
        """
        ## there's nothing to show
        count = ProjectSample.objects.filter(
            project__id=record.id, is_deleted=False, is_error=False, is_finished=True
        ).count()
        count_not_finished = ProjectSample.objects.filter(
            project__id=record.id, is_deleted=False, is_error=False, is_finished=False
        ).count()

        sz_project_sample = ""
        if count > 0:
            sz_project_sample = (
                "<a href="
                + reverse("show-sample-project-results", args=[record.pk])
                + ' data-toggle="tooltip" title="See Results"> '
                + '<span ><i class="padding-button-table fa fa-info-circle padding-button-table"></i></span></a> '
            )
            ## only can change settings when has projects finished
            # sz_project_sample += '<a href=' + reverse('project-settings', args=[record.pk]) + ' data-toggle="tooltip" title="Software settings">' +\
            #     '<span ><i class="fa fa-magic padding-button-table"></i></span></a>'
            sz_project_sample += (
                "<a href="
                + reverse("project-settings", args=[record.pk])
                + ' data-toggle="tooltip" title="Software settings">'
                + '<span ><i class="padding-button-table fa fa-magic padding-button-table"></i></span></a>'
            )
        elif count_not_finished > 0:
            sz_project_sample = _("{} processing ".format(count_not_finished))
        else:  ## can change settings
            sz_project_sample += (
                "<a href="
                + reverse("project-settings", args=[record.pk])
                + ' data-toggle="tooltip" title="Software settings">'
                + '<span ><i class="padding-button-table fa fa-magic padding-button-table"></i></span></a>'
            )

        return mark_safe(sz_project_sample)


class ShowProjectSamplesResults(tables.Table):
    """
    Has the results from the result samples processed against references
    """

    manage_database = ManageDatabase()

    sample_name = tables.Column("Sample name", empty_values=())
    coverage = tables.Column("Coverage", orderable=False, empty_values=())
    alerts = tables.Column("Alerts", empty_values=())
    type_and_subtype = tables.LinkColumn("Classification", empty_values=())
    putative_mixed_infection = tables.LinkColumn(
        "Putative Mixed-infection", empty_values=()
    )
    technology = tables.Column("Technology", empty_values=())
    dataset = tables.LinkColumn("Dataset", empty_values=())
    results = tables.LinkColumn("Options", orderable=False, empty_values=())
    consensus_file = tables.LinkColumn(
        "Consensus File", orderable=False, empty_values=()
    )

    class Meta:
        model = ProjectSample
        fields = (
            "sample_name",
            "type_and_subtype",
            "putative_mixed_infection",
            "technology",
            "dataset",
            "coverage",
            "consensus_file",
            "alerts",
            "results",
        )
        attrs = {"class": "table-striped table-bordered"}
        empty_text = "There are no samples processed to show..."

    def render_sample_name(self, record):
        """
        return sample name
        """
        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()
        user = current_request.user
        if user.username == Constants.USER_ANONYMOUS:
            return record.sample.name
        if user.username == record.project.owner.username:
            return mark_safe(
                '<a href="#id_remove_modal" id="id_remove_reference_modal" data-toggle="modal"'
                + ' ref_name="'
                + record.sample.name
                + '" pk="'
                + str(record.pk)
                + '" +\
                    " ref_project="'
                + record.project.name
                + '" data-toggle="tooltip" title="Remove sample">'
                + '<i class="fa fa-trash"></i></span> </a>'
                + "<a href="
                + reverse("show-sample-project-single-detail", args=[record.pk])
                + ' data-toggle="tooltip" title="Show more information">'
                + "{}</a>".format(record.sample.name)
            )
        return record.sample.name

    def render_coverage(self, record):
        """
        return icons about coverage
        """
        manageDatabase = ManageDatabase()
        meta_value = manageDatabase.get_project_sample_metakey_last(
            record,
            MetaKeyAndValue.META_KEY_Coverage,
            MetaKeyAndValue.META_VALUE_Success,
        )

        ### coverage
        decode_coverage = DecodeObjects()
        coverage = decode_coverage.decode_result(meta_value.description)

        ### default parameters
        default_software = DefaultProjectSoftware()
        limit_to_mask_consensus = int(
            default_software.get_mask_consensus_single_parameter(
                record,
                DefaultParameters.MASK_CONSENSUS_threshold,
                ConstantsSettings.TECHNOLOGY_illumina
                if record.is_sample_illumina()
                else ConstantsSettings.TECHNOLOGY_minion,
            )
        )
        return_html = ""
        for key in coverage.get_sorted_elements_name():
            return_html += '<a href="#coverageModal" id="showImageCoverage" data-toggle="modal" project_sample_id="{}" '.format(
                record.id
            ) + 'sequence="{}"><img title="{}" class="tip" src="{}"></a>'.format(
                key,
                coverage.get_message_to_show_in_web_site(record.sample.name, key),
                coverage.get_icon(key, limit_to_mask_consensus),
            )
        return mark_safe(return_html)

    def render_technology(self, record):
        """shows if it is Illumina or Minion"""
        return (
            ConstantsSettings.TECHNOLOGY_illumina
            if record.sample.is_type_fastq_gz_sequencing()
            else ConstantsSettings.TECHNOLOGY_minion
        )

    def render_alerts(self, record):
        """
        return number
        """
        return "{}".format(record.alert_first_level + record.alert_second_level)

    def render_type_and_subtype(self, record):
        """
        return number
        """
        return record.sample.type_subtype

    def render_putative_mixed_infection(self, record):
        """
        return number
        """
        return (
            "" if record.mixed_infections is None else record.mixed_infections.tag.name
        )

    def render_dataset(self, record):
        """
        return number
        """
        return record.sample.data_set.name

    def render_consensus_file(self, record):
        """
        return link to consensus file
        """
        default_software = DefaultProjectSoftware()
        return record.get_consensus_file_web(
            not default_software.include_consensus(record)
        )

    def render_results(self, record):
        """
        icon with link to extra info
        """
        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()
        user = current_request.user
        str_links = ""
        if user.username != Constants.USER_ANONYMOUS:
            str_links = (
                "<a href="
                + reverse("sample-project-settings", args=[record.pk])
                + ' data-toggle="tooltip" title="Software settings">'
                + '<span ><i class="padding-button-table fa fa-magic padding-button-table"></i></span></a>'
            )

        b_space_add = False
        if record.is_mask_consensus_sequences:
            str_links += (
                '<a href="javascript:void(0);" data-toggle="tooltip" title="Sample run with user-selected parameters (see update 30 Oct 2020)">'
                + '<span ><i class="padding-button-table fa fa-calendar-check-o"></i></a>'
            )
        else:
            b_space_add = True
            str_links += "  "

        ### check if it was downsized
        if self.manage_database.is_sample_downsized(record.sample):
            str_links += (
                '<a href="javascript:void(0);" data-toggle="tooltip" title="Sample was downsized">'
                + '<span ><i class="padding-button-table warning_fa_icon fa fa-level-down"></i></a>'
            )
        elif not b_space_add:
            str_links += "  "

        str_links += (
            "<a href="
            + reverse("show-sample-project-single-detail", args=[record.pk])
            + ' data-toggle="tooltip" title="Show more information" class="padding-button-table">'
            + '<span ><i class="fa fa-info-circle"></i> More info</a>'
        )
        return mark_safe(str_links)

    def order_sample_name(self, queryset, is_descending):
        queryset = queryset.annotate(sample_name=F("sample__name")).order_by(
            ("-" if is_descending else "") + "sample_name"
        )
        return (queryset, True)

    def order_type_and_subtype(self, queryset, is_descending):
        queryset = queryset.annotate(type_subtype=F("sample__type_subtype")).order_by(
            ("-" if is_descending else "") + "type_subtype"
        )
        return (queryset, True)

    def order_dataset(self, queryset, is_descending):
        queryset = queryset.annotate(dataset=F("sample__data_set__name")).order_by(
            ("-" if is_descending else "") + "dataset"
        )
        return (queryset, True)

    def order_alerts(self, queryset, is_descending):
        queryset = queryset.annotate(
            alerts=F("alert_first_level") + F("alert_second_level")
        ).order_by(("-" if is_descending else "") + "alerts")
        return (queryset, True)

    def order_putative_mixed_infection(self, queryset, is_descending):
        queryset = queryset.annotate(
            mixed_infection=F("mixed_infections__tag__name")
        ).order_by(("-" if is_descending else "") + "mixed_infection")
        return (queryset, True)

    def order_technology(self, queryset, is_descending):
        """shows if it is Illumina or Minion"""
        queryset = queryset.annotate(technology=F("sample__type_of_fastq")).order_by(
            ("-" if is_descending else "") + "technology"
        )
        return (queryset, True)


class AddSamplesFromCvsFileTable(tables.Table):
    """
    To add samples to projects
    """

    samples_processed = tables.Column("#Samples processed", empty_values=())
    number_samples = tables.Column("#Samples", empty_values=())
    is_completed = tables.Column("Load completed", empty_values=())

    class Meta:
        model = Sample
        fields = (
            "file_name",
            "creation_date",
            "owner",
            "number_samples",
            "samples_processed",
            "is_completed",
        )
        attrs = {"class": "table-striped table-bordered"}
        empty_text = "There are no (csv) or (tsv) files with samples to add..."

    def render_creation_date(self, value, record):
        return record.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)

    def render_number_samples(self, value, record):
        return record.number_files_to_process

    def render_samples_processed(self, value, record):
        return record.number_files_processed

    def render_is_completed(self, value, record):
        return "True" if record.is_processed else "False"

    def render_file_name(self, value, record):
        href = record.get_path_to_file(TypePath.MEDIA_URL)
        return mark_safe('<a href="' + href + '" download>' + record.file_name + "</a>")

    def order_is_completed(self, queryset, is_descending):
        queryset = queryset.annotate(is_completed=F("is_processed")).order_by(
            ("-" if is_descending else "") + "is_processed"
        )
        return (queryset, True)

    def order_samples_processed(self, queryset, is_descending):
        queryset = queryset.annotate(
            samples_processed=F("number_files_processed")
        ).order_by(("-" if is_descending else "") + "number_files_processed")
        return (queryset, True)

    def order_number_samples(self, queryset, is_descending):
        queryset = queryset.annotate(
            number_samples=F("number_files_to_process")
        ).order_by(("-" if is_descending else "") + "number_files_to_process")
        return (queryset, True)


class AddSamplesFromCvsFileTableMetadata(tables.Table):
    """
    To add samples to projects
    """

    samples_processed = tables.Column("#Samples linked", empty_values=())
    number_samples = tables.Column("#Samples", empty_values=())
    is_completed = tables.Column("Load completed", empty_values=())

    class Meta:
        model = Sample
        fields = (
            "file_name",
            "creation_date",
            "owner",
            "number_samples",
            "samples_processed",
            "is_completed",
        )
        attrs = {"class": "table-striped table-bordered"}
        empty_text = "There are no (csv) or (tsv) files with samples to add..."

    def render_creation_date(self, value, record):
        return record.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)

    def render_number_samples(self, value, record):
        return record.number_files_to_process

    def render_samples_processed(self, value, record):
        return record.number_files_processed

    def render_is_completed(self, value, record):
        return "True" if record.is_processed else "False"

    def render_file_name(self, value, record):
        href = record.get_path_to_file(TypePath.MEDIA_URL)
        return mark_safe('<a href="' + href + '" download>' + record.file_name + "</a>")

    def order_is_completed(self, queryset, is_descending):
        queryset = queryset.annotate(is_completed=F("is_processed")).order_by(
            ("-" if is_descending else "") + "is_processed"
        )
        return (queryset, True)

    def order_samples_processed(self, queryset, is_descending):
        queryset = queryset.annotate(
            samples_processed=F("number_files_processed")
        ).order_by(("-" if is_descending else "") + "number_files_processed")
        return (queryset, True)

    def order_number_samples(self, queryset, is_descending):
        queryset = queryset.annotate(
            number_samples=F("number_files_to_process")
        ).order_by(("-" if is_descending else "") + "number_files_to_process")
        return (queryset, True)


class AddSamplesFromFastqFileTable(tables.Table):
    """
    To add samples to projects
    """

    is_completed = tables.Column("Is it attached?", empty_values=())
    sample_attached = tables.Column("Sample attached", empty_values=())

    class Meta:
        model = Sample
        fields = (
            "file_name",
            "creation_date",
            "owner",
            "attached_date",
            "sample_attached",
            "is_completed",
        )
        attrs = {"class": "table-striped table-bordered"}
        empty_text = "There's no 'fastq' files with to show..."

    def render_creation_date(self, value, record):
        return record.creation_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)

    def render_attached_date(self, value, record):
        if record.attached_date == None:
            return "---"
        return record.attached_date.strftime(settings.DATETIME_FORMAT_FOR_TABLE)

    def render_sample_attached(self, value, record):
        samples = record.samples.all()
        if samples.count() > 0:
            return samples[0].name
        return "---"

    def render_is_completed(self, value, record):
        return "True" if record.is_processed else "False"

    def render_file_name(self, value, record):
        from crequest.middleware import CrequestMiddleware

        current_request = CrequestMiddleware.get_request()
        user = current_request.user
        remove_str = ""
        if user.username == Constants.USER_ANONYMOUS:
            return record.file_name
        if user.username == record.owner.username:
            if not record.is_processed:
                remove_str = mark_safe(
                    '<a href="#id_remove_modal" id="id_remove_reference_modal" data-toggle="modal"'
                    + ' ref_name="'
                    + record.file_name
                    + '" pk="'
                    + str(record.pk)
                    + '"><i class="fa fa-trash"></i></span> </a> '
                )
            else:  ## if the attached sample is deleted you can remove the fastq.gz
                b_all_deleted = True
                for sample in record.samples.all():
                    if not sample.is_deleted:
                        b_all_deleted = False
                        break

                if b_all_deleted:
                    remove_str = mark_safe(
                        '<a href="#id_remove_modal" id="id_remove_reference_modal" data-toggle="modal"'
                        + ' ref_name="'
                        + record.file_name
                        + '" pk="'
                        + str(record.pk)
                        + '"><i class="fa fa-trash"></i></span> </a> '
                    )

        href = record.get_path_to_file(TypePath.MEDIA_URL)
        return mark_safe(
            remove_str + ' <a href="' + href + '" download>' + record.file_name + "</a>"
        )

    def order_is_completed(self, queryset, is_descending):
        queryset = queryset.annotate(is_completed=F("is_processed")).order_by(
            ("-" if is_descending else "") + "is_processed"
        )
        return (queryset, True)

    def order_sample_attached(self, queryset, is_descending):
        queryset = queryset.annotate(sample_name=F("samples__name")).order_by(
            ("-" if is_descending else "") + "sample_name"
        )
        return (queryset, True)
