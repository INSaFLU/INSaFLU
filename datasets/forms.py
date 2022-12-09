'''
Created on 31/05/2022

@author: mmp
'''
import os
from django import forms
from managing_files.models import Reference, Project
from datasets.models import UploadFiles
from datasets.models import Consensus, Dataset
from utils.utils import Utils
from django.conf import settings
from django.template.defaultfilters import filesizeformat
from crispy_forms.helper import FormHelper
from django.core.files.temp import NamedTemporaryFile
from crispy_forms.layout import Div, Layout, ButtonHolder, Submit, Button, HTML
from django.urls import reverse
from constants.constants import Constants
from utils.software import Software
from Bio import SeqIO
from utils.parse_in_files_nextstrain import ParseNextStrainFiles

class AddReferencesDatasetForm(forms.ModelForm):
	"""
	References model to dataset 
	"""
	error_css_class = 'error'
	
	class Meta:
		model = Reference
		exclude = ()

	def __init__(self, *args, **kwargs):
		super(AddReferencesDatasetForm, self).__init__(*args, **kwargs)
		
	def clean(self):
		"""
		"""
		cleaned_data = super(AddReferencesDatasetForm, self).clean()
		return cleaned_data
	
	
class AddConsensusDatasetForm(forms.ModelForm):
	"""
	Add consensus to Insaflu 
	"""
	error_css_class = 'error'
	
	class Meta:
		model = Consensus
		exclude = ()

	def __init__(self, *args, **kwargs):
		super(AddConsensusDatasetForm, self).__init__(*args, **kwargs)
		
	def clean(self):
		"""
		"""
		cleaned_data = super(AddConsensusDatasetForm, self).clean()
		return cleaned_data

class AddProjectsDatasetForm(forms.ModelForm):
	"""
	Add project to Insaflu 
	"""
	error_css_class = 'error'
	
	class Meta:
		model = Project
		exclude = ()

	def __init__(self, *args, **kwargs):
		super(AddProjectsDatasetForm, self).__init__(*args, **kwargs)
		
	def clean(self):
		"""
		"""
		cleaned_data = super(AddProjectsDatasetForm, self).clean()
		return cleaned_data
	
## https://kuanyui.github.io/2015/04/13/django-crispy-inline-form-layout-with-bootstrap/
class ConsensusForm(forms.ModelForm):
	"""
	Reference form, name, isolate_name and others
	"""
	software = Software()
	utils = Utils()
	error_css_class = 'error'
	
	class Meta:
		model = Consensus
		# specify what fields should be used in this form.
		fields = ('name', 'display_name', 'consensus_fasta')
		
	def __init__(self, *args, **kwargs):
		self.request = kwargs.pop('request')
		self.pk = kwargs.pop('pk')
		super(ConsensusForm, self).__init__(*args, **kwargs)

		## can exclude explicitly
		## exclude = ('md5',)
		field_text= [
			# (field_name, Field title label, Detailed field description, requiered)
			('name', 'Prefix Name', 'Prefix name to attach to the sequences names in fasta file. ' +\
			'<p><b><i>The prefix can be useful to select specific group of sequences to be added to datasets.</i></b></p><p>If empty, ' +\
			'only the names of the sequences will taken in consideration.</p><p>It can be a multi-fasta sequence file. If the name ' +\
			'already exists in database the sequence will be rejected.</p>', False),
			('display_name', 'Sequence names to upload, can be separated by comma', 'If you want to upload specific sequences from multifasta sequences, '
			'set the names where separated by comma, if more than one. If empty, upload all.', False),
			('consensus_fasta', 'Consensus (Fasta/Multi-Fasta)', "Consensus file in fasta format.<br>" +\
				"Max total sequence length: {}<br>".format(filesizeformat(settings.MAX_LENGTH_SEQUENCE_TOTAL_FROM_CONSENSUS_FASTA))  +\
				"Max FASTA file size: {}".format(filesizeformat(settings.MAX_CONSENSUS_FASTA_FILE)), True),
		]
		for x in field_text:
			self.fields[x[0]].label = x[1]
			self.fields[x[0]].help_text = x[2]
			self.fields[x[0]].required = x[3]

## in case you have to undo to specific ID
##		cancel_url = reverse('references')
##		if self.instance: cancel_url = reverse('references', kwargs={'pk': self.instance.id})

		self.helper = FormHelper()
		self.helper.form_method = 'POST'
		self.helper.attrs["data-validate-consensus-url"] = reverse('validate-consensus-name')
		self.helper.attrs["id"] = "id_form_consensus"
		self.helper.layout = Layout(
			Div(
				Div('name', css_class="col-sm-7"),
				css_class = 'row'
			),
			Div(
				Div('display_name', css_class="col-sm-7"),
				css_class = 'row'
			),
			Div('consensus_fasta', css_class = 'show-for-sr'),
			ButtonHolder(
				Submit('save', 'Save', css_class='btn-primary', onclick="this.disabled=true,this.form.submit();"),
				Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(
					reverse('add-consensus-dataset', args=[self.pk])))
			)
		)

	def clean(self):
		"""
		Clean all together because it's necessary to compare the genbank and fasta files
		"""
		cleaned_data = super(ConsensusForm, self).clean()
		## It will be perfix, can be empty
		name = self.cleaned_data.get('name', '').strip()
		vect_names_to_upload = self.cleaned_data.get('display_name', '').strip().split(',') if \
			len(self.cleaned_data.get('display_name', '').strip()) > 0 else []
		dict_names = dict(zip(vect_names_to_upload, [0] * len(vect_names_to_upload)))
		
		# if (len(name) == 0):
		#	 self.add_error('name', "Error: You must give a unique name.")
		#	 return cleaned_data
		#
		# try:
		#	 Consensus.objects.get(name__iexact=name, owner=self.request.user, is_obsolete=False, is_deleted=False)
		#	 self.add_error('name', "This name '" + name +"' already exist in database, please choose other.")
		#	 return cleaned_data
		# except Consensus.DoesNotExist:
		#	 pass
		
		## test reference_fasta
		if ('consensus_fasta' not in cleaned_data):
			self.add_error('consensus_fasta', "Error: Must have a Fasta/Multi-Fasta file.")
			return cleaned_data
		
		### testing file names
		reference_fasta = cleaned_data['consensus_fasta']
		
		## testing fasta
		some_error_in_files = False
		reference_fasta_temp_file_name = NamedTemporaryFile(prefix='flu_fa_', delete=False)
		reference_fasta_temp_file_name.write(reference_fasta.read())
		reference_fasta_temp_file_name.flush()
		reference_fasta_temp_file_name.close()
		self.software.dos_2_unix(reference_fasta_temp_file_name.name)
		try:
			number_locus = self.utils.is_fasta(reference_fasta_temp_file_name.name)
			self.request.session[Constants.NUMBER_LOCUS_FASTA_FILE] = number_locus
			self.request.session[Constants.SEQUENCES_TO_PASS] = ""
			
			## test the max numbers
			if (number_locus > Constants.MAX_SEQUENCES_FROM_CONTIGS_FASTA):
				self.add_error('consensus_fasta', 'Max allow number of contigs in Multi-Fasta: {}'.format(Constants.MAX_SEQUENCES_FROM_CONTIGS_FASTA))
				some_error_in_files = True
			total_length_fasta = self.utils.get_total_length_fasta(reference_fasta_temp_file_name.name)
			if (not some_error_in_files and total_length_fasta > settings.MAX_LENGTH_SEQUENCE_TOTAL_FROM_CONSENSUS_FASTA):
				some_error_in_files = True
				self.add_error('consensus_fasta', 'The max sum length of the sequences in fasta: {}'.format(settings.MAX_LENGTH_SEQUENCE_TOTAL_FROM_CONSENSUS_FASTA))
			
			n_seq_name_bigger_than = self.utils.get_number_seqs_names_bigger_than(reference_fasta_temp_file_name.name,
							Constants.MAX_LENGTH_CONTIGS_SEQ_NAME, len(name))
			if (not some_error_in_files and n_seq_name_bigger_than > 0):
				some_error_in_files = True
				if (n_seq_name_bigger_than == 1):
					self.add_error('consensus_fasta', 'There is one sequence name length bigger than {0}. The max. length name is {0}.'.format(Constants.MAX_LENGTH_CONTIGS_SEQ_NAME))
				else:
					self.add_error('consensus_fasta', 'There are {0} sequences with name length bigger than {1}. The max. length name is {1}.'.format(n_seq_name_bigger_than, Constants.MAX_LENGTH_CONTIGS_SEQ_NAME))
					
			## if some errors in the files, fasta or genBank, return
			if (some_error_in_files): return cleaned_data
				
			### check if there all seq names are present in the database yet
			b_pass = False
			vect_error, vect_fail_seqs, vect_pass_seqs = [], [], []
			with open(reference_fasta_temp_file_name.name) as handle_in:
				for record in SeqIO.parse(handle_in, "fasta"):
					## only these ones can get in
					if record.id in dict_names: dict_names[record.id] = 1
					if (len(vect_names_to_upload)) > 0 and not record.id in vect_names_to_upload:
						vect_fail_seqs.append(record.id)
						continue
					
					## try to upload
					try:
						seq_name = "{}_{}".format(name, record.id) if len(name) > 0 else record.id
						Consensus.objects.get(name__iexact=seq_name, owner=self.request.user, is_obsolete=False, is_deleted=False)
						vect_error.append("Seq. name: '" + seq_name +"' already exist in database.")
					except Consensus.DoesNotExist:
						vect_pass_seqs.append(record.id)
						b_pass = True

			## if none of them pass throw an error
			if not b_pass: 
				some_error_in_files = True
				if len(vect_names_to_upload) > 0 and len(vect_pass_seqs) == 0:
					self.add_error('display_name', "None of these names '{}' match to the sequences names".format(
						"', '".join(vect_names_to_upload)))
				for message in vect_error: self.add_error('consensus_fasta', message)
			else:
				## if empty load all
				self.request.session[Constants.SEQUENCES_TO_PASS] = ",".join(vect_pass_seqs)
				
			### some sequences names suggested are not present in the file
			vect_fail_seqs = [key for key in dict_names if dict_names[key] == 0]
			if len(vect_fail_seqs) > 0:
				self.add_error('display_name', "Sequences names '{}' does not have match in the file".format(
						", '".join(vect_fail_seqs)))
				
		except IOError as e:	## (e.errno, e.strerror)
			os.unlink(reference_fasta_temp_file_name.name)
			some_error_in_files = True
			self.add_error('consensus_fasta', e.args[0])
		except ValueError as e:	## (e.errno, e.strerror)
			os.unlink(reference_fasta_temp_file_name.name)
			some_error_in_files = True
			self.add_error('consensus_fasta', e.args[0])
		except: 
			os.unlink(reference_fasta_temp_file_name.name)
			some_error_in_files = True
			self.add_error('consensus_fasta', "Not a valid 'fasta' file.")

		### test if it has degenerated bases
		# if (os.path.exists(reference_fasta_temp_file_name.name)):
		#	 try:
		#		 self.utils.has_degenerated_bases(reference_fasta_temp_file_name.name)
		#	 except Exception as e:
		#		 os.unlink(reference_fasta_temp_file_name.name)
		#		 some_error_in_files = True
		#		 self.add_error('consensus_fasta', e.args[0])
			
		## remove temp files
		os.unlink(reference_fasta_temp_file_name.name)
		return cleaned_data



class DatastesUploadDescriptionMetadataForm(forms.ModelForm):
	"""
	Upload csv file to update metadata 
	"""
	error_css_class = 'error'
	
	class Meta:
		model = UploadFiles
		exclude = ()
		fields = ('path_name', )

	def __init__(self, *args, **kwargs):
		self.request = kwargs.pop('request')
		self.pk = kwargs.pop('pk')
		super(DatastesUploadDescriptionMetadataForm, self).__init__(*args, **kwargs)
		
		dataset_name = Dataset.objects.get(pk=self.pk) 
		
		field_text= [
			('path_name', 'File name', '"csv" or "tsv" file with metadata to update existing samples.', True),
		]
		for x in field_text:
			self.fields[x[0]].label = x[1]
			self.fields[x[0]].help_text = x[2]
			self.fields[x[0]].required = x[3]
		
		self.helper = FormHelper()
		self.helper.form_method = 'POST'
		self.helper.layout = Layout(
#			HTML('<p> </p>'),
#			HTML('<div class="alert alert-dark"> <span><i class="fa fa-download"></i></span> ' + \
#					dataset_name.get_global_file_by_dataset_web(Dataset.DATASET_FILE_NAME_RESULT_NEXTSTRAIN_CSV) + \
#					" Last metadata file 'csv' for '{}' dataset</div>".format(
#					dataset_name.name)),
			HTML('<p> </p>'),
			HTML('<div class="alert alert-dark"> <span><i class="fa fa-download"></i></span> ' + \
					 dataset_name.get_global_file_by_dataset_web(Dataset.DATASET_FILE_NAME_RESULT_NEXTSTRAIN_TSV) + \
					" Last metadata file 'tsv' for '{}' dataset.</div>".format(
					dataset_name.name)),
			HTML('<p> </p>'),
			Div('path_name', css_class="col-lm-3"),
			HTML('<p> </p>'),
			ButtonHolder(
				Submit('save', 'Upload', css_class='btn-primary', onclick="this.disabled=true,this.form.submit();"),
				Button('cancel', 'Cancel', css_class='btn-secondary', 
					   onclick='window.location.href="{}"'.format(reverse('dataset-update-metadata', args=[self.pk])))
			),
		)

	def clean(self):
		"""
		Clean all 
		"""
		software = Software()
		cleaned_data = super(DatastesUploadDescriptionMetadataForm, self).clean()
		
		### get path name
		path_name = self.cleaned_data['path_name']
		
		## testing fastq
		temp_file_name = NamedTemporaryFile(prefix='flu_fq_', suffix='.csv', delete=False)
		temp_file_name.write(path_name.file.read())
		temp_file_name.flush()
		temp_file_name.close()
		software.dos_2_unix(temp_file_name.name)
		
		parse_in_files = ParseNextStrainFiles()
		b_test_char_encoding = True
		parse_in_files.parse_nextstrain_files(temp_file_name.name, self.request.user, b_test_char_encoding,\
									ParseNextStrainFiles.STATE_READ_metadata_only_detect_errors_and_chech_nexttrain)
		
		os.unlink(temp_file_name.name)
		if (parse_in_files.get_errors().has_errors()):
			self.add_error('path_name', 'There are errors in the file')
			self.error_in_file = str(parse_in_files.get_errors())		## pass error_file kwargs for context
		else:
			self.number_files_to_process = parse_in_files.get_number_samples()
		return cleaned_data



