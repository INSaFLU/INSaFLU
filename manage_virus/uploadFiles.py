'''
Created on Oct 31, 2017

@author: mmp
'''
from django.db import transaction
from utils.Constants import Constants
from utils.ParseOutFiles import ParseOutFiles
from utils.Utils import Utils
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
				if (not self.utils.is_fasta(os.path.join(path_to_find, file))): continue
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
		if (len(split_name_2) != 3): raise Exception("Error: this seq. name '" + split_name + "' has a wrong format. Must have 'influenza_type|subtype~~~<value>~~~<accession> description'")
		split_type = split_name_2[0].split("_")
		if (len(split_type) != 2 and (split_type[1].lower() != Constants.SEQ_VIRUS_TYPE.lower() or split_type[1].lower() != Constants.SEQ_VIRUS_SUB_TYPE.lower()) ):
			raise Exception("Error: this seq. name '" + split_name + "' has a wrong format. Must have 'influenza_type|subtype~~~<value>~~~<accession> description'")
		if (split_type[1].lower() == Constants.SEQ_VIRUS_TYPE.lower()): dt_return[Constants.SEQ_VIRUS_TYPE] = split_name_2[1]
		else: dt_return[Constants.SEQ_VIRUS_SUB_TYPE] = split_name_2[1]
		dt_return[self.ACCESSION] = split_name_2[2]
		return dt_return

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

			seqVirus.accession = dt_data[self.ACCESSION]
			seqVirus.file = uploadFile
			seqVirus.save()
			if (self.DESCRIPTION in dt_data): seqVirus.description = dt_data[self.DESCRIPTION]

		return uploadFile
		
		
	@transaction.atomic
	def uploadIdentifyVirus(self, vect_data, database_name):
		
		### first look for the type
		dt_type = self.__get_type__(vect_data)
		dt_first_sub_type = self.__get_sub_type__(vect_data, None)
		dt_second_sub_type = self.__get_sub_type__(vect_data, dt_first_sub_type)
		
		vect_return = []
		rank = 0
		if (len(dt_type) > 0):
			seqVirus = None
			try:
				seqVirus = SeqVirus.objects.get(name=dt_type[ParseOutFiles.GENE], kind_type__name=Constants.SEQ_VIRUS_TYPE, file__abricate_name=database_name)
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
		
		if (len(dt_first_sub_type) > 0):
			seqVirus = None
			try:
				seqVirus = SeqVirus.objects.get(name=dt_first_sub_type[ParseOutFiles.GENE], kind_type__name=Constants.SEQ_VIRUS_SUB_TYPE, file__abricate_name=database_name)
			except SeqVirus.DoesNotExist:
				raise Exception(_("Gene '%s' not found in database '%s'" % (dt_first_sub_type[ParseOutFiles.GENE], database_name)))
			
			identifyVirus = IdentifyVirus()
			identifyVirus.coverage = "%.2f" % (dt_first_sub_type[ParseOutFiles.COVERAGE])
			identifyVirus.identity = "%.2f" % (dt_first_sub_type[ParseOutFiles.IDENTITY])
			identifyVirus.rank = rank
			identifyVirus.seq_virus = seqVirus
			identifyVirus.save()
			rank += 1
			vect_return.append(identifyVirus)
			
		if (len(dt_second_sub_type) > 0):
			seqVirus = None
			try:
				seqVirus = SeqVirus.objects.get(name=dt_second_sub_type[ParseOutFiles.GENE], kind_type__name=Constants.SEQ_VIRUS_SUB_TYPE, file__abricate_name=database_name)
			except SeqVirus.DoesNotExist:
				raise Exception(_("Gene '%s' not found in database '%s'" % (dt_second_sub_type[ParseOutFiles.GENE], database_name)))
			
			identifyVirus = IdentifyVirus()
			identifyVirus.coverage = "%.2f" % (dt_second_sub_type[ParseOutFiles.COVERAGE])
			identifyVirus.identity = "%.2f" % (dt_second_sub_type[ParseOutFiles.IDENTITY])
			identifyVirus.rank = rank
			identifyVirus.seq_virus = seqVirus
			identifyVirus.save()
			rank += 1
			vect_return.append(identifyVirus)
		return vect_return


	def __get_type__(self, vect_data):
		"""
		Get the best type
		"""
		for dt_data in vect_data:
			if (ParseOutFiles.TYPE in dt_data):
				split_type = dt_data[ParseOutFiles.TYPE].split("_")
				if (len(split_type) == 2 and split_type[1].lower() == Constants.SEQ_VIRUS_TYPE.lower()):
					return dt_data
		return {}
	
	def __get_sub_type__(self, vect_data, dt_already_returned):
		"""
		Get the best type
		"""
		for dt_data in vect_data:
			if (ParseOutFiles.TYPE in dt_data):
				split_type = dt_data[ParseOutFiles.TYPE].split("_")
				if (len(split_type) == 2 and split_type[1].lower() == Constants.SEQ_VIRUS_SUB_TYPE.lower()):
					if (dt_already_returned == None): return dt_data
					if (dt_already_returned[ParseOutFiles.GENE] != dt_data[ParseOutFiles.GENE]): return dt_data
		return {}
		