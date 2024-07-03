import os

from Bio import SeqIO
from crispy_forms.helper import FormHelper
from crispy_forms.layout import HTML, Button, ButtonHolder, Div, Layout, Submit
from django import forms
from django.conf import settings
from django.core.files.temp import NamedTemporaryFile
from django.template.defaultfilters import filesizeformat
from django.urls import reverse

from constants.constants import Constants
from managing_files.models import Reference
from pathogen_identification.models import Projects
from utils.software import Software
from utils.utils import Utils


class ReferenceForm(forms.Form):
    search = forms.CharField(
        label="Search",
        widget=forms.TextInput(
            attrs={
                "class": "form-control",
                "placeholder": "Search",
                "name": "search_add_reference",
            }
        ),
    )


class UploadFileForm(forms.Form):
    description_help = (
        "Enter a short description for the upload file. "
        "This description will be used to identify the file in the future."
    )
    description_fasta = (
        "Upload a FASTA file. The file must contain the sequences of the reference. "
        + "The sequences must be in FASTA format. "
        + "<br>Sequence IDs must be unique. Text after space will be ignored."
        + "<br>Max total sequence length: {} bp.<br>".format(
            settings.MAX_LENGTH_SEQUENCE_TOTAL_FROM_FASTA * 10
        )
        + "Max FASTA file size: {}".format(
            filesizeformat(settings.MAX_FASTQ_FILE_UPLOAD)
        )
    )
    description_metadata = (
        "Upload a metadata file. The metadata file must be in CSV or TSV format with appropriate extentions."
        + "<br>File must contain the following columns:"
        + "<br>      <bold>Accession ID</bold>: The accession identifier the sequence in the FASTA file."
        + "<br><bold>TaxID</bold>: The NCBI Taxonomy ID of the sequence. "
        + "<br><bold>Description (optional)</bold>: A description of the sequence."
        + "<br>Only record IDs that are present in both the FASTA and metadata files will be uploaded."
    )

    description = forms.CharField(widget=forms.Textarea, help_text=description_help)
    fasta_file = forms.FileField(help_text=description_fasta)
    metadata = forms.FileField(help_text=description_metadata)

    def __init__(self, *args, **kwargs):
        ## add ids to filefields
        print(args, kwargs)
        self.request = kwargs.pop("request")
        super(UploadFileForm, self).__init__(*args, **kwargs)

        self.fields["fasta_file"].widget.attrs["id"] = "fasta_file"
        self.fields["metadata"].widget.attrs["id"] = "metadata"

    def is_valid1(self) -> bool:
        """
        Check if the form is valid
        """
        valid = super(UploadFileForm, self).is_valid()
        if not valid:
            return valid

        fasta_file = self.cleaned_data.get("fasta_file")
        metadata = self.cleaned_data.get("metadata")

        print("#TESTING")
        print(fasta_file)
        print(metadata)

        if fasta_file is None:
            self.add_error("fasta_file", "Please upload a FASTA file.")
            return False

        if metadata is None:
            self.add_error("metadata", "Please upload a metadata file.")
            return False

        return True


## https://kuanyui.github.io/2015/04/13/django-crispy-inline-form-layout-with-bootstrap/
class PanelReferencesUploadForm(forms.ModelForm):
    """
    Reference form, name, isolate_name and others
    """

    software = Software()
    utils = Utils()
    error_css_class = "error"

    class Meta:
        model = Reference
        # specify what fields should be used in this form.
        fields = ("name", "display_name", "reference_fasta")

    def __init__(self, *args, **kwargs):
        self.request = kwargs.pop("request")
        self.pk = kwargs.pop("pk")
        super(PanelReferencesUploadForm, self).__init__(*args, **kwargs)

        ## can exclude explicitly
        ## exclude = ('md5',)
        field_text = [
            # (field_name, Field title label, Detailed field description, requiered)
            (
                "name",
                "Prefix Name",
                "Prefix name to attach to the sequences names in fasta file. "
                + "<p><b><i>The prefix can be useful to select specific group of sequences to be added to datasets.</i></b></p><p>If empty, "
                + "only the names of the sequences will taken in consideration.</p><p>It can be a multi-fasta sequence file. If the name "
                + "already exists in database the sequence will be rejected.</p>",
                False,
            ),
            (
                "display_name",
                "Sequence names to upload, can be separated by comma",
                "If you want to upload specific sequences from multifasta sequences, "
                "set the names where separated by comma, if more than one. If empty, upload all.",
                False,
            ),
            (
                "reference_fasta",
                "Reference (Fasta/Multi-Fasta)",
                "Reference file in fasta format.<br>"
                + "Max total sequence length: {} bp<br>".format(
                    settings.MAX_LENGTH_SEQUENCE_TOTAL_FROM_FASTA * 10
                )
                + "Max FASTA file size: {}".format(
                    filesizeformat(settings.MAX_FASTQ_FILE_UPLOAD)
                ),
                True,
            ),
        ]
        for x in field_text:
            self.fields[x[0]].label = x[1]
            self.fields[x[0]].help_text = x[2]
            self.fields[x[0]].required = x[3]

        ## in case you have to undo to specific ID
        ##		cancel_url = reverse('references')
        ##		if self.instance: cancel_url = reverse('references', kwargs={'pk': self.instance.id})

        self.helper = FormHelper()
        self.helper.form_method = "POST"
        self.helper.attrs["data-validate-consensus-url"] = reverse(
            "validate-consensus-name"
        )
        self.helper.attrs["id"] = "id_form_consensus"
        self.helper.layout = Layout(
            Div(Div("name", css_class="col-sm-7"), css_class="row"),
            Div(Div("display_name", css_class="col-sm-7"), css_class="row"),
            Div("reference_fasta", css_class="show-for-sr"),
            ButtonHolder(
                Submit(
                    "save",
                    "Save",
                    css_class="btn-primary",
                    onclick="this.disabled=true,this.form.submit();",
                ),
                Button(
                    "cancel",
                    "Cancel",
                    css_class="btn-secondary",
                    onclick='window.location.href="{}"'.format(
                        reverse(
                            "reference_panels",
                        )
                    ),
                ),
            ),
        )

    def clean(self):
        """
        Clean all together because it's necessary to compare the genbank and fasta files
        """
        cleaned_data = super(PanelReferencesUploadForm, self).clean()
        ## It will be perfix, can be empty
        name = self.cleaned_data.get("name", "").strip()
        vect_names_to_upload = (
            self.cleaned_data.get("display_name", "").strip().split(",")
            if len(self.cleaned_data.get("display_name", "").strip()) > 0
            else []
        )

        dict_names = dict(zip(vect_names_to_upload, [0] * len(vect_names_to_upload)))
        ## test reference_fasta
        if "reference_fasta" not in cleaned_data:
            self.add_error(
                "reference_fasta", "Error: Must have a Fasta/Multi-Fasta file."
            )
            return cleaned_data

        ### testing file names
        reference_fasta = cleaned_data["reference_fasta"]
        print("HIOUHIOUB")

        ## testing fasta
        some_error_in_files = False
        reference_fasta_temp_file_name = NamedTemporaryFile(
            prefix="flu_fa_", delete=False
        )
        reference_fasta_temp_file_name.write(reference_fasta.read())
        reference_fasta_temp_file_name.flush()
        reference_fasta_temp_file_name.close()
        self.software.dos_2_unix(reference_fasta_temp_file_name.name)
        try:
            number_locus = self.utils.is_fasta(reference_fasta_temp_file_name.name)
            self.request.session[Constants.NUMBER_LOCUS_FASTA_FILE] = number_locus
            self.request.session[Constants.SEQUENCES_TO_PASS] = ""

            ## test the max numbers
            if number_locus > Constants.MAX_SEQUENCES_FROM_CONTIGS_FASTA:
                self.add_error(
                    "reference_fasta",
                    "Max allow number of contigs in Multi-Fasta: {}".format(
                        Constants.MAX_SEQUENCES_FROM_CONTIGS_FASTA
                    ),
                )
                some_error_in_files = True
            total_length_fasta = self.utils.get_total_length_fasta(
                reference_fasta_temp_file_name.name
            )
            if (
                not some_error_in_files
                and total_length_fasta
                > settings.MAX_LENGTH_SEQUENCE_TOTAL_FROM_FASTA * 10
            ):
                some_error_in_files = True
                self.add_error(
                    "reference_fasta",
                    "The max sum length of the sequences in fasta: {}".format(
                        settings.MAX_LENGTH_SEQUENCE_TOTAL_FROM_FASTA * 10
                    ),
                )

            fasta_file_size = os.path.getsize(reference_fasta_temp_file_name.name)
            if (
                not some_error_in_files
                and fasta_file_size > settings.MAX_FASTQ_FILE_UPLOAD
            ):
                some_error_in_files = True
                self.add_error(
                    "reference_fasta",
                    "The max size of the fasta file: {}".format(
                        filesizeformat(settings.MAX_FASTQ_FILE_UPLOAD)
                    ),
                )

            n_seq_name_bigger_than = self.utils.get_number_seqs_names_bigger_than(
                reference_fasta_temp_file_name.name,
                Constants.MAX_LENGTH_CONTIGS_SEQ_NAME,
                len(name),
            )

            if not some_error_in_files and n_seq_name_bigger_than > 0:
                some_error_in_files = True
                if n_seq_name_bigger_than == 1:
                    self.add_error(
                        "reference_fasta",
                        "There is one sequence name length bigger than {0}. The max. length name is {0}.".format(
                            Constants.MAX_LENGTH_CONTIGS_SEQ_NAME
                        ),
                    )
                else:
                    self.add_error(
                        "reference_fasta",
                        "There are {0} sequences with name length bigger than {1}. The max. length name is {1}.".format(
                            n_seq_name_bigger_than,
                            Constants.MAX_LENGTH_CONTIGS_SEQ_NAME,
                        ),
                    )

            ## if some errors in the files, fasta or genBank, return
            if some_error_in_files:
                return cleaned_data

            ### check if there all seq names are present in the database yet
            b_pass = False
            vect_error, vect_fail_seqs, vect_pass_seqs = [], [], []
            with open(reference_fasta_temp_file_name.name) as handle_in:
                for record in SeqIO.parse(handle_in, "fasta"):
                    ## only these ones can get in
                    if record.id in dict_names:
                        dict_names[record.id] = 1
                    if (
                        len(vect_names_to_upload)
                    ) > 0 and not record.id in vect_names_to_upload:
                        vect_fail_seqs.append(record.id)
                        continue

                    ## try to upload
                    try:
                        seq_name = (
                            "{}_{}".format(name, record.id)
                            if len(name) > 0
                            else record.id
                        )
                        Reference.objects.get(
                            name__iexact=seq_name,
                            owner=self.request.user,
                            is_obsolete=False,
                            is_deleted=False,
                        )
                        vect_error.append(
                            "Seq. name: '" + seq_name + "' already exist in database."
                        )
                    except Reference.DoesNotExist:
                        vect_pass_seqs.append(record.id)
                        b_pass = True

            ## if none of them pass throw an error
            if not b_pass:
                some_error_in_files = True
                if len(vect_names_to_upload) > 0 and len(vect_pass_seqs) == 0:
                    self.add_error(
                        "display_name",
                        "None of these names '{}' match to the sequences names".format(
                            "', '".join(vect_names_to_upload)
                        ),
                    )
                for message in vect_error:
                    self.add_error("reference_fasta", message)
            else:
                ## if empty load all
                self.request.session[Constants.SEQUENCES_TO_PASS] = ",".join(
                    vect_pass_seqs
                )

            ### some sequences names suggested are not present in the file
            vect_fail_seqs = [key for key in dict_names if dict_names[key] == 0]
            if len(vect_fail_seqs) > 0:
                self.add_error(
                    "display_name",
                    "Sequences names '{}' does not have match in the file".format(
                        ", '".join(vect_fail_seqs)
                    ),
                )

        except IOError as e:  ## (e.errno, e.strerror)
            os.unlink(reference_fasta_temp_file_name.name)
            some_error_in_files = True
            self.add_error("reference_fasta", e.args[0])
        except ValueError as e:  ## (e.errno, e.strerror)
            os.unlink(reference_fasta_temp_file_name.name)
            some_error_in_files = True
            self.add_error("reference_fasta", e.args[0])
        except:
            os.unlink(reference_fasta_temp_file_name.name)
            some_error_in_files = True
            self.add_error("reference_fasta", "Not a valid 'fasta' file.")

        ## remove temp files
        os.unlink(reference_fasta_temp_file_name.name)
        return cleaned_data
