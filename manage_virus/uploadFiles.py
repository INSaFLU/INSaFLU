'''
Created on Oct 31, 2017

@author: mmp
'''
from django.db import transaction
from constants.constants import Constants, FileExtensions, TypePath
from manage_virus.constants_virus import ConstantsVirus
from utils.parse_out_files import ParseOutFiles
from utils.utils import Utils
from django.conf import settings
from .models import UploadFile, Tags, SeqVirus, IdentifyVirus
from django.utils.translation import ugettext_lazy as _
from collections import OrderedDict
from Bio import SeqIO
import os, re, logging


class UploadFiles(object):
	'''
	classdocs
	'''

	constants = Constants()
	utils = Utils()
	
	ACCESSION = "Accesion"
	DESCRIPTION = "Description"
	
	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
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
		if (len(split_type) != 2 and (not split_type[1].lower() in ConstantsVirus.VECT_ALL_POSSIBLE_TAGS) ):
			raise Exception("Error: this seq. name '" + split_name + "' has a wrong format. Must have 'influenza_type|subtype|lineage~~~<value>~~~<accession> description'")
		dt_return[self.__get_type__(split_name_2[0])] = split_name_2[1]
		dt_return[self.ACCESSION] = split_name_2[2]
		return dt_return

	def __get_type__(self, split_type_name):
		split_type = split_type_name.split("_")
		if (len(split_type) < 2): return ConstantsVirus.SEQ_VIRUS_SUB_TYPE
		for tag in ConstantsVirus.VECT_ALL_POSSIBLE_TAGS:
			if (split_type[1].lower() == tag.lower()): return tag
		return ConstantsVirus.SEQ_VIRUS_SUB_TYPE
		
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

		## load the file
		record_dict = SeqIO.index(path, "fasta")
		for seq in record_dict:
			dt_data = self.parse_title_name(seq)
			seqVirus = SeqVirus()
			for possible_tag in ConstantsVirus.VECT_ALL_POSSIBLE_TAGS:
				if (possible_tag in dt_data):
					try:
						tag_type = Tags.objects.get(name=possible_tag)
					except Tags.DoesNotExist:	## not exist
						tag_type = Tags()
						tag_type.name = possible_tag
						tag_type.save()
					seqVirus.kind_type = tag_type
					seqVirus.name = dt_data[possible_tag]
					break

			if not seqVirus.kind_type is None:
				seqVirus.accession = dt_data.get(self.ACCESSION)
				seqVirus.description = dt_data.get(self.DESCRIPTION)
				seqVirus.file = uploadFile
				seqVirus.save()
			else:
				### error
				self.logger_production.error('Fail to upload identification tag: ' + seq)
				self.logger_debug.error('Fail to upload identification tag: ' + seq)
		return uploadFile
		
		
	@transaction.atomic
	def uploadIdentifyVirus(self, dict_data, database_name):
		"""
		:in dict_data -> dictonary { 'Seq.id' : [dt_data, dt_data2, ...], 'Seq.id1' : [dt_data5, dt_data6, ...]
		OrderedDict(sorted(dict_data_out.items(), reverse=True,
#			key=lambda k: (k[1][0][self.COVERAGE],
#			k[1][0][self.IDENTITY]) )), clean_abricate_file
		"""
		### first look for the type
		### this is necessary because of the order
		dict_data_out = OrderedDict(sorted(dict_data.items(), reverse=True,
								key=lambda k: (k[1][0][ParseOutFiles.COVERAGE],
									k[1][0][ParseOutFiles.IDENTITY]) ))

		vect_data_out = self.__get_type_data__(dict_data_out)
		vect_data_type_and_lineage = self.__get_sub_type_and_lineage__(dict_data_out)
		
		vect_return = []
		rank = 0
		if (len(vect_data_out + vect_data_type_and_lineage) > 0):
			for dt_type in vect_data_out + vect_data_type_and_lineage:
				seqVirus = None
				try:
					seqVirus = SeqVirus.objects.get(name=dt_type.get(ParseOutFiles.GENE), accession=dt_type.get(ParseOutFiles.ACCESSION),\
									file__abricate_name=database_name)
				except SeqVirus.DoesNotExist:
					raise Exception(_("Gene '%s', Accession '%s' not found in database '%s'" % (dt_type.get(ParseOutFiles.GENE),\
								dt_type.get(ParseOutFiles.ACCESSION), database_name)))
				
				identifyVirus = IdentifyVirus()
				identifyVirus.coverage = "%.2f" % (dt_type.get(ParseOutFiles.COVERAGE))
				identifyVirus.identity = "%.2f" % (dt_type.get(ParseOutFiles.IDENTITY))
				identifyVirus.rank = rank
				identifyVirus.seq_virus = seqVirus
				identifyVirus.save()
				rank += 1
				vect_return.append(identifyVirus)
		
		return vect_return


	def __get_type_data__(self, dict_data_ordered):
		"""
		Get the best type
		"""
		vect_data_out = []
		dict_out_type = {}
		for key in dict_data_ordered:
			for dict_info in dict_data_ordered[key]:
				if (ParseOutFiles.TYPE in dict_info):
					split_type = dict_info[ParseOutFiles.TYPE].split("_")
					if (len(split_type) == 2 and not (dict_info[ParseOutFiles.GENE] in dict_out_type) and \
						split_type[1].lower() in ConstantsVirus.VECT_SEQ_VIRUS_TYPE):
						dict_out_type[dict_info[ParseOutFiles.GENE]] = 1
						vect_data_out.append(dict_info)
		return vect_data_out
	
	def __get_sub_type_and_lineage__(self, dict_data_ordered):
		"""
		Get the all subtype and lineage
		"""
		vect_data_out = []
		dict_out_gene = {}
		for key in dict_data_ordered:
			for dict_info in dict_data_ordered[key]:
				if (ParseOutFiles.TYPE in dict_info):
					type_result = self.__get_type__(dict_info[ParseOutFiles.TYPE])
					if (type_result.lower() in ConstantsVirus.VECT_SEQ_VIRUS_SUB_TYPE and \
						not dict_info[ParseOutFiles.GENE] in dict_out_gene):
							vect_data_out.append(dict_info)
							dict_out_gene[dict_info[ParseOutFiles.GENE]] = 1
		return vect_data_out
	
	
	@transaction.atomic
	def upload_default_references(self, user, b_test):
		"""
		upload default files for reference
		"""
		from managing_files.models import Reference
		from utils.software import Software
		
		software = Software()
		path_to_find = os.path.join(getattr(settings, "STATIC_ROOT", None), \
					Constants.DIR_TEST_TYPE_REFERENCES if b_test else Constants.DIR_TYPE_REFERENCES)
		n_upload = 0
		for file in self.utils.get_all_files(path_to_find):
			file = os.path.join(path_to_find, file)
			
			#### don't make genbank files
			try:
				number_of_elements = self.utils.is_genbank(file)
				continue
			except IOError as e:
				pass
			
			try:
				number_of_elements = self.utils.is_fasta(file)
			except IOError as e:
				print(e.args[0])
				continue
			except ValueError as e:    ## (e.errno, e.strerror)
				print(e.args[0])
				continue
			
			name = self.utils.clean_extension(os.path.basename(file))
			try:
				reference = Reference.objects.get(owner=user, is_obsolete=False, is_deleted=False, name__iexact=name)
			except Reference.DoesNotExist as e:
				
				### check if exist genbank, and valid
				gbk_file_name = None
				number_of_elements_gbk = -1
				for extension in FileExtensions.VECT_ALL_GBK_EXTENSIONS:
					gbk_file_name = os.path.join(os.path.dirname(file), name + extension)
					try:
						number_of_elements_gbk = self.utils.is_genbank(gbk_file_name)
						break
					except IOError as e:
						gbk_file_name = None
			
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
				if (not gbk_file_name is None and number_of_elements_gbk == number_of_elements):
					reference.reference_genbank_name = os.path.basename(gbk_file_name)
				else:
					gbk_file_name = None
					reference.reference_genbank_name = name + FileExtensions.FILE_GBK
				reference.save()
				
				## move the files to the right place
				sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), self.utils.get_path_to_reference_file(user.id, reference.id), reference.reference_fasta_name)
				self.utils.copy_file(file, sz_file_to)
				reference.reference_fasta.name = os.path.join(self.utils.get_path_to_reference_file(user.id, reference.id), reference.reference_fasta_name)
				
				### genbank files
				sz_file_to = os.path.join(getattr(settings, "MEDIA_ROOT", None), self.utils.get_path_to_reference_file(user.id, reference.id), reference.reference_genbank_name)
				if (gbk_file_name is None):
					temp_dir = software.run_prokka(file, os.path.basename(file))
					self.utils.move_file(os.path.join(temp_dir, reference.reference_genbank_name), sz_file_to)
					self.utils.remove_dir(temp_dir)
				else:
					self.utils.copy_file(gbk_file_name, sz_file_to)

				reference.reference_genbank.name = os.path.join(self.utils.get_path_to_reference_file(user.id, reference.id), reference.reference_genbank_name)
				reference.hash_reference_genbank = self.utils.md5sum(sz_file_to)
				reference.save()
				
				### create bed and index for genbank
				self.utils.from_genbank_to_bed(sz_file_to, reference.get_reference_bed(TypePath.MEDIA_ROOT))
				software.create_index_files_from_igv_tools(reference.get_reference_bed(TypePath.MEDIA_ROOT))
				software.run_genbank2gff3(sz_file_to, reference.get_gff3(TypePath.MEDIA_ROOT))
				with_gene_annotation = True
				software.run_genbank2gff3(sz_file_to, reference.get_gff3_with_gene_annotation(TypePath.MEDIA_ROOT), with_gene_annotation)
				software.run_genbank2gff3_positions_comulative(sz_file_to, reference.get_gff3_comulative_positions(TypePath.MEDIA_ROOT))
				### save in database the elements and coordinates
				self.utils.get_elements_from_db(reference, user)
				self.utils.get_elements_and_cds_from_db(reference, user)

				## create the index before commit in database, throw exception if something goes wrong
				software.create_fai_fasta(os.path.join(getattr(settings, "MEDIA_ROOT", None), reference.reference_fasta.name))
		
				print("Upload reference: {}".format(name))
				n_upload += 1
			
			if (b_test and n_upload > 2): break
		print("#Upload references: {}".format(n_upload))
