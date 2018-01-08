'''
Created on Oct 31, 2017

@author: mmp
'''
from django.db import transaction
from constants.constants import Constants, FileExtensions, TypePath
from utils.parse_out_files import ParseOutFiles
from utils.utils import Utils
from django.conf import settings
from .models import UploadFile, Tags, SeqVirus, IdentifyVirus
from django.utils.translation import ugettext_lazy as _
import os, re
from Bio import SeqIO

class UploadFiles(object):
	'''
	classdocs
	'''

	constants = Constants()
	utils = Utils()
	
	ACCESSION = "Accesion"
	DESCRIPTION = "Description"
	
	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	def get_file_to_upload(self, b_testing):
		"""
		"""
		version = 0
		path_to_return = None
		path_to_find = os.path.join( getattr(settings, "STATIC_ROOT", None), "tests" if b_testing else "", self.constants.DIR_TYPE_IDENTIFICATION)
		for file in self.utils.get_all_files(path_to_find):
			base_name = os.path.basename(file)
			match = re.search('\w+(_[v|V]\d+)\.\w+', base_name)
			if (match == None): continue
			if (len(match.regs) == 2):
				try:
					self.utils.is_fasta(os.path.join(path_to_find, file))
				except IOError as e:
					continue
				
				(temp, path) = (int(match.group(1).lower().replace("_v", "")), os.path.join(path_to_find, file))
				if (temp > version):
					version = temp
					path_to_return = path
		return (str(version), path_to_return)


	def parse_title_name(self, header_name):
		"""
		"""
		dt_return = {}
		split_name = header_name.split(" ")
		if (len(split_name) > 1): dt_return[self.DESCRIPTION] = split_name[1]
		split_name_2 = re.split("~+", split_name[0])
		if (len(split_name_2) != 3): raise Exception("Error: this seq. name '" + split_name + "' has a wrong format. Must have 'influenza_type|subtype|lineage~~~<value>~~~<accession> description'")
		split_type = split_name_2[0].split("_")
		if (len(split_type) != 2 and (split_type[1].lower() != Constants.SEQ_VIRUS_TYPE.lower() or split_type[1].lower() != Constants.SEQ_VIRUS_SUB_TYPE.lower()
									 or split_type[1].lower() != Constants.SEQ_VIRUS_LINEAGE.lower()) ):
			raise Exception("Error: this seq. name '" + split_name + "' has a wrong format. Must have 'influenza_type|subtype|lineage~~~<value>~~~<accession> description'")
		dt_return[self.__get_type__(split_name_2[0])] = split_name_2[1]
		dt_return[self.ACCESSION] = split_name_2[2]
		return dt_return

	def __get_type__(self, split_type_name):
		split_type = split_type_name.split("_")
		if (len(split_type) < 2): return Constants.SEQ_VIRUS_SUB_TYPE
		if (split_type[1].lower() == Constants.SEQ_VIRUS_TYPE.lower()): return Constants.SEQ_VIRUS_TYPE
		elif (split_type[1].lower() == Constants.SEQ_VIRUS_LINEAGE.lower()): return Constants.SEQ_VIRUS_LINEAGE
		return Constants.SEQ_VIRUS_SUB_TYPE
		
	@transaction.atomic
	def upload_file(self, version, path):
		"""
		Upload fasta files with this
		influenza_type~~~B~~~AF100378
		atomic transaction for database 
		"""
		base_name = os.path.basename(path)
		try:
			uploadFile = UploadFile.objects.get(name=base_name)
			return uploadFile ### it's already there
		except UploadFile.DoesNotExist:	## not exist
			uploadFile = UploadFile()
			uploadFile.name = base_name
			uploadFile.path = path
			uploadFile.version = version
			uploadFile.abricate_name = self.utils.get_file_name_without_extension(base_name)
			uploadFile.save()

		## type_tag
		try:
			tag_type = Tags.objects.get(name=Constants.SEQ_VIRUS_TYPE)
		except Tags.DoesNotExist:	## not exist
			tag_type = Tags()
			tag_type.name = Constants.SEQ_VIRUS_TYPE
			tag_type.save()
		
		## sub_type_tag
		try:
			tag_sub_type = Tags.objects.get(name=Constants.SEQ_VIRUS_SUB_TYPE)
		except Tags.DoesNotExist:	## not exist
			tag_sub_type = Tags()
			tag_sub_type.name = Constants.SEQ_VIRUS_SUB_TYPE
			tag_sub_type.save()
			
		## sub_lineage
		try:
			tag_leneage = Tags.objects.get(name=Constants.SEQ_VIRUS_LINEAGE)
		except Tags.DoesNotExist:	## not exist
			tag_leneage = Tags()
			tag_leneage.name = Constants.SEQ_VIRUS_LINEAGE
			tag_leneage.save()
		
		## load the file
		record_dict = SeqIO.index(path, "fasta")
		for seq in record_dict:
			dt_data = self.parse_title_name(seq)
			seqVirus = SeqVirus()
			if (Constants.SEQ_VIRUS_TYPE in dt_data): 
				seqVirus.kind_type = tag_type
				seqVirus.name = dt_data[Constants.SEQ_VIRUS_TYPE]
			elif (Constants.SEQ_VIRUS_SUB_TYPE in dt_data):
				seqVirus.kind_type = tag_sub_type
				seqVirus.name = dt_data[Constants.SEQ_VIRUS_SUB_TYPE]
			elif (Constants.SEQ_VIRUS_LINEAGE in dt_data):
				seqVirus.kind_type = tag_leneage
				seqVirus.name = dt_data[Constants.SEQ_VIRUS_LINEAGE]

			seqVirus.accession = dt_data[self.ACCESSION]
			seqVirus.file = uploadFile
			seqVirus.save()
			if (self.DESCRIPTION in dt_data): seqVirus.description = dt_data[self.DESCRIPTION]

		return uploadFile
		
		
	@transaction.atomic
	def uploadIdentifyVirus(self, vect_data, database_name):
		
		### first look for the type
		dt_type = self.__get_type_data__(vect_data)
		vect_data_type_and_lineage = self.__get_sub_type_and_lineage__(vect_data)
		
		vect_return = []
		rank = 0
		if (len(dt_type) > 0):
			seqVirus = None
			try:
				seqVirus = SeqVirus.objects.get(name=dt_type[ParseOutFiles.GENE], accession=dt_type[ParseOutFiles.ACCESSION],\
						kind_type__name=Constants.SEQ_VIRUS_TYPE, file__abricate_name=database_name)
			except SeqVirus.DoesNotExist:
				raise Exception(_("Gene '%s' not found in database '%s'" % (dt_type[ParseOutFiles.GENE], database_name)))
			
			identifyVirus = IdentifyVirus()
			identifyVirus.coverage = "%.2f" % (dt_type[ParseOutFiles.COVERAGE])
			identifyVirus.identity = "%.2f" % (dt_type[ParseOutFiles.IDENTITY])
			identifyVirus.rank = rank
			identifyVirus.seq_virus = seqVirus
			identifyVirus.save()
			rank += 1
			vect_return.append(identifyVirus)
		
		if (len(vect_data_type_and_lineage) > 0):
			for dt_data in vect_data_type_and_lineage:
				seqVirus = None
				try:
					seqVirus = SeqVirus.objects.get(name=dt_data[ParseOutFiles.GENE], kind_type__name=self.__get_type__(dt_data[ParseOutFiles.TYPE]), file__abricate_name=database_name)
				except SeqVirus.DoesNotExist:
					raise Exception(_("Gene '%s' not found in database '%s'" % (dt_data[ParseOutFiles.GENE], database_name)))
				
				identifyVirus = IdentifyVirus()
				identifyVirus.coverage = "%.2f" % (dt_data[ParseOutFiles.COVERAGE])
				identifyVirus.identity = "%.2f" % (dt_data[ParseOutFiles.IDENTITY])
				identifyVirus.rank = rank
				identifyVirus.seq_virus = seqVirus
				identifyVirus.save()
				rank += 1
				vect_return.append(identifyVirus)
		return vect_return


	def __get_type_data__(self, vect_data):
		"""
		Get the best type
		"""
		for dt_data in vect_data:
			if (ParseOutFiles.TYPE in dt_data):
				split_type = dt_data[ParseOutFiles.TYPE].split("_")
				if (len(split_type) == 2 and split_type[1].lower() == Constants.SEQ_VIRUS_TYPE.lower()):
					return dt_data
		return {}
	
	def __get_sub_type_and_lineage__(self, vect_data):
		"""
		Get the all subtype and lineage
		"""
		vect_data_out = []
		vect_out_gene = []
		for dt_data in vect_data:
			if (ParseOutFiles.TYPE in dt_data):
				type_result = self.__get_type__(dt_data[ParseOutFiles.TYPE])
				if (type_result == Constants.SEQ_VIRUS_SUB_TYPE or type_result == Constants.SEQ_VIRUS_LINEAGE):
					if (dt_data[ParseOutFiles.GENE] in vect_out_gene): continue
					vect_data_out.append(dt_data)
					vect_out_gene.append(dt_data[ParseOutFiles.GENE])
		return vect_data_out
	
	
	@transaction.atomic
	def upload_default_references(self, user, b_test):
		"""
		upload default files for reference
		"""
		from managing_files.models import Reference
		from utils.software import Software
		
		software = Software()
		path_to_find = os.path.join(getattr(settings, "STATIC_ROOT", None), Constants.DIR_TYPE_REFERENCES)
		n_upload = 0
		for file in self.utils.get_all_files(path_to_find):
			file = os.path.join(path_to_find, file)
			try:
				number_of_elements = self.utils.is_fasta(file)
			except IOError as e:
				print(e.args[0])
				continue
			
			name = self.utils.clean_extension(os.path.basename(file))
			try:
				reference = Reference.objects.get(owner=user, is_obsolete=False, is_deleted=False, name__iexact=name)
			except Reference.DoesNotExist as e:
				reference = Reference()
				reference.display_name = name
				reference.isolate_name = name
				reference.name = reference.display_name
				reference.owner = user
				reference.is_obsolete = False
				reference.number_of_locus = number_of_elements
				reference.hash_reference_fasta = self.utils.md5sum(file)
				reference.reference_fasta_name = os.path.basename(file)
				reference.scentific_name = os.path.basename(file)
				reference.reference_genbank_name = name + FileExtensions.FILE_GBK
				reference.save()
				
				## move the files to the right place
				sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), self.utils.get_path_to_reference_file(user.id, reference.id), reference.reference_fasta_name)
				self.utils.copy_file(file, sz_file_to)
				reference.reference_fasta.name = os.path.join(self.utils.get_path_to_reference_file(user.id, reference.id), reference.reference_fasta_name)
				
				temp_dir = software.run_prokka(file, os.path.basename(file))
				sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), self.utils.get_path_to_reference_file(user.id, reference.id), reference.reference_genbank_name)
				self.utils.move_file(os.path.join(temp_dir, reference.reference_genbank_name), sz_file_to)
				reference.reference_genbank.name = os.path.join(self.utils.get_path_to_reference_file(user.id, reference.id), reference.reference_genbank_name)
				reference.hash_reference_genbank = self.utils.md5sum(sz_file_to)
				reference.save()
				
				### create bed and index for genbank
				self.utils.from_genbank_to_bed(sz_file_to, reference.get_reference_bed(TypePath.MEDIA_ROOT))
				software.create_index_files_from_igv_tools(reference.get_reference_bed(TypePath.MEDIA_ROOT))
				
				### save in database the elements and coordinates
				self.utils.get_elements_from_db(reference, user)
				self.utils.get_elements_and_cds_from_db(reference, user)

				## create the index before commit in database, throw exception if something goes wrong
				software.create_fai_fasta(os.path.join(getattr(settings, "MEDIA_ROOT", None), reference.reference_fasta.name))
		
				### remove dir
				self.utils.remove_dir(temp_dir)
				n_upload += 1
			
			if (b_test and n_upload > 2): break
