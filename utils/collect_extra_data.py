'''
Created on Nov 27, 2017

@author: mmp
'''
import os, csv, time, json, logging
import plotly.graph_objs as go
from utils.utils import Utils 
from managing_files.manage_database import ManageDatabase
from managing_files.models import Project, TagNames, ProcessControler
from managing_files.models import Software as SoftwareModel
from constants.meta_key_and_values import MetaKeyAndValue
from utils.result import DecodeObjects
from constants.constants import TypePath, Constants, FileType, FileExtensions
from constants.software_names import SoftwareNames
from utils.tree import CreateTree
from plotly.offline import plot
from settings.default_software_project_sample import DefaultProjectSoftware
from settings.default_parameters import DefaultParameters
from utils.result import Coverage, Result, SoftwareDesc
from utils.software_pangolin import SoftwarePangolin
from utils.parse_out_files import ParseOutFiles
from utils.process_SGE import ProcessSGE
from django.conf import settings
from settings.constants_settings import ConstantsSettings
from utils.result import MaskingConsensus
from Bio import SeqIO

class CollectExtraData(object):
	'''
	classdocs
	'''

	HEADER_SAMPLE_OUT_CSV = "id,fastq1,fastq2,data set,vaccine status,week,onset date,collection date,lab reception date,latitude,longitude,type-subtype,putative mixed-infection"
	HEADER_SAMPLE_OUT_CSV_simple = "id,data set,vaccine status,week,onset date,collection date,lab reception date,latitude,longitude,type-subtype,putative mixed-infection"

	utils = Utils()
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	def collect_extra_data_for_project(self, project, user):
		"""
		"""
		### make it running 
		process_controler = ProcessControler()
		process_SGE = ProcessSGE()
		process_SGE.set_process_controler(user, process_controler.get_name_project(project), ProcessControler.FLAG_RUNNING)
		
		## need to add a delay for the test in command line
		if (settings.RUN_TEST_IN_COMMAND_LINE): time.sleep(4)
		
		## run collect data
		self.__collect_extra_data_for_project(project, user)

	def collect_update_extra_metadata_for_project(self, project, user):
		"""
		Only for update metadata
		"""
		### make it running 
		process_controler = ProcessControler()
		process_SGE = ProcessSGE()
		process_SGE.set_process_controler(user, process_controler.get_name_project(project), ProcessControler.FLAG_RUNNING)
		
		## need to add a delay for the test in command line
		if (settings.RUN_TEST_IN_COMMAND_LINE): time.sleep(4)
		
		## run collect data
		self.__collect_update_extra_metadata_for_project(project, user)
		
	def collect_update_pangolin_lineage(self, project, user):
		"""
		Only for update metadata
		"""
		### make it running 
		process_controler = ProcessControler()
		process_SGE = ProcessSGE()
		process_SGE.set_process_controler(user, process_controler.get_name_project(project), ProcessControler.FLAG_RUNNING)
		
		## need to add a delay for the test in command line
		if (settings.RUN_TEST_IN_COMMAND_LINE): time.sleep(4)
		
		## need to do this because old projects doesn't have defined the "project_sample.seq_name_all_consensus" field
		## collect all consensus files for a project_sample
		self.calculate_global_files(Project.PROJECT_FILE_NAME_SAMPLE_RESULT_all_consensus, project, user)
			
		## run pangolin collect data
		self.__collect_update_pangolin_lineage(project, user, False)
		
		## run collect data
		self.__collect_update_extra_metadata_for_project(project, user)


	def __collect_update_pangolin_lineage(self, project, user, b_mark_sge_success = True):
		"""
		Only for update metadata
		"""
		### get the taskID and seal it
		software_pangolin = SoftwarePangolin()
		process_controler = ProcessControler()
		process_SGE = ProcessSGE()
		try:
			## collect sample table with plus type and subtype, mixed infection, equal to upload table
			file_pangolin_output = project.get_global_file_by_project(TypePath.MEDIA_ROOT,
										Project.PROJECT_FILE_NAME_Pangolin_lineage)
			file_consensus = project.get_global_file_by_project(TypePath.MEDIA_ROOT,
										Project.PROJECT_FILE_NAME_SAMPLE_RESULT_all_consensus)
			## test if is necesserary to run pangolin lineage
			if (os.path.exists(file_pangolin_output) or software_pangolin.is_ref_sars_cov_2(
					project.reference.get_reference_fasta(TypePath.MEDIA_ROOT))):
				## process pangolin
				software_pangolin.run_pangolin(file_consensus, file_pangolin_output)
				
				try:
					software = SoftwareModel.objects.get(name=SoftwareNames.SOFTWARE_Pangolin_name)
					
					### set last version of the Pangolin run in this project
					dt_result_version = software.get_versions()
					result_all = Result()
					manage_database = ManageDatabase()
					software_names = SoftwareNames()
					result_all.add_software(SoftwareDesc(software_names.get_pangolin_name(),
							dt_result_version.get(software_names.get_pangolin_name(), ""), ""))
					result_all.add_software(SoftwareDesc(software_names.get_pangolin_learn_name(),
							dt_result_version.get(software_names.get_pangolin_learn_name(), ""), ""))
					result_all.add_software(SoftwareDesc(software_names.get_pangolin_designation_name(),
							dt_result_version.get(software_names.get_pangolin_designation_name(), ""), ""))
					manage_database.set_project_metakey(project, user,
							MetaKeyAndValue.META_KEY_Identify_pangolin,\
							MetaKeyAndValue.META_VALUE_Success,\
							result_all.to_json() )
				except SoftwareModel.DoesNotExist:	## need to create with last version
					self.logger_production.error('ProjectID: {} Fail to detect software model'.format(project.id))
					self.logger_debug.error('ProjectID: {} Fail to detect software model'.format(project.id))
				except:
					self.logger_production.error('ProjectID: {}  Fail to run pangolin'.format(project.id))
					self.logger_debug.error('ProjectID: {}  '.format(project.id))
		except:
			## finished with error
			process_SGE.set_process_controler(user, process_controler.get_name_project(project), ProcessControler.FLAG_ERROR)
			return
		
		if (b_mark_sge_success):
			process_SGE.set_process_controler(user, process_controler.get_name_project(project), ProcessControler.FLAG_FINISHED)


	def __collect_update_extra_metadata_for_project(self, project, user):
		"""
		Only for update metadata
		"""
		### get the taskID and seal it
		process_controler = ProcessControler()
		process_SGE = ProcessSGE()
		manage_database = ManageDatabase()
		metaKeyAndValue = MetaKeyAndValue()
		
		try:
			## collect sample table with plus type and subtype, mixed infection, equal to upload table
			self.calculate_global_files(Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV, project, user)
			self.calculate_global_files(Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV_simple, project, user)
			self.calculate_global_files(Project.PROJECT_FILE_NAME_SAMPLE_RESULT_TSV, project, user)
			## IMPORTANT -> this need to be after of Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV
			self.calculate_global_files(Project.PROJECT_FILE_NAME_SAMPLE_RESULT_json, project, user)
		except:
			## finished with error
			process_SGE.set_process_controler(user, process_controler.get_name_project(project), ProcessControler.FLAG_ERROR)
			return
		
		## seal the tag		
		meta_project = manage_database.get_project_metakey_last(project, metaKeyAndValue.get_meta_key(\
						MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id), MetaKeyAndValue.META_VALUE_Queue)
		if (meta_project != None):
			manage_database.set_project_metakey(project, user, metaKeyAndValue.get_meta_key(\
					MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id),
					MetaKeyAndValue.META_VALUE_Success, meta_project.description)
				
		### finished
		process_SGE.set_process_controler(user, process_controler.get_name_project(project), ProcessControler.FLAG_FINISHED)
		
##	@transaction.atomic
	def __collect_extra_data_for_project(self, project, user):
		"""
		Everything that is necessary to do in the project
		Collect all extra data after all samples are finished
		only run after all the vect_taskID are finished
		"""
		### get the taskID and seal it
		metaKeyAndValue = MetaKeyAndValue()
		manage_database = ManageDatabase()
		process_controler = ProcessControler()
		process_SGE = ProcessSGE()
		
		#### create variation graph, png and html
		## Obsolete, is to make a html graph, now it is with chart.js
# 		(out_file_html, out_file_png) = self.create_graph_minor_variants(project, user)
# 		file_destination = project.get_global_file_by_project(TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_GRAPH_MINO_VAR_HTML)
# 		if (out_file_html != None):
# 			self.utils.copy_file(out_file_html, file_destination)
# 			os.unlink(out_file_html)
# 		elif (os.path.exists(file_destination)): os.unlink(file_destination)
# 		file_destination = project.get_global_file_by_project(TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_GRAPH_MINO_VAR_PNG)
# 		if (out_file_png != None):
# 			self.utils.copy_file(out_file_png, file_destination)
# 			os.unlink(out_file_png)
# 		elif (os.path.exists(file_destination)): os.unlink(file_destination)

		try:
			## calculate the max sample label size of the samples that belong to this project
			## used in MSA viewer 
			b_calculate_again = True
			manage_database.get_max_length_label(project, user, b_calculate_again)
			
			## masks all consensus first, must be before of merge AllConsensus (Project.PROJECT_FILE_NAME_SAMPLE_RESULT_all_consensus)
			self.calculate_global_files(Project.PROJECT_FILE_NAME_SAMPLE_mask_all_consensus, project, user)
			
			## collect all consensus files for a project_sample
			self.calculate_global_files(Project.PROJECT_FILE_NAME_SAMPLE_RESULT_all_consensus, project, user)
			
			## calculate the lineage, if necessary
			## This need to be after the PROJECT_FILE_NAME_SAMPLE_RESULT_all_consensus
			self.__collect_update_pangolin_lineage(project, user, False)
			
			### calculate global file
			self.calculate_global_files(Project.PROJECT_FILE_NAME_COVERAGE, project, user)
			## collect tab variations snippy
			self.calculate_global_files(Project.PROJECT_FILE_NAME_TAB_VARIATIONS_SNIPPY, project, user)
			## collect tab variations freebayes, <50%
			## remove del or ins
			self.calculate_global_files(Project.PROJECT_FILE_NAME_TAB_VARIATIONS_FREEBAYES, project, user)
			## with snp, del and ins
			self.calculate_global_files(Project.PROJECT_FILE_NAME_TAB_VARIATIONS_FREEBAYES_with_snps_indels, project, user)
			## collect sample table with plus type and subtype, mixed infection, equal to upload table
			self.calculate_global_files(Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV, project, user)
			self.calculate_global_files(Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV_simple, project, user)
			self.calculate_global_files(Project.PROJECT_FILE_NAME_SAMPLE_RESULT_TSV, project, user)
			## IMPORTANT -> this need to be after of Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV
			self.calculate_global_files(Project.PROJECT_FILE_NAME_SAMPLE_RESULT_json, project, user)
		
			## calculate global variations for a project
			self.calculate_count_variations(project)
			
			### create trees
			createTree = CreateTree()
			createTree.create_tree_and_alignments(project, user)
			
			meta_project = manage_database.get_project_metakey_last(project, metaKeyAndValue.get_meta_key(\
						MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id), MetaKeyAndValue.META_VALUE_Queue)
			if (meta_project != None):
				manage_database.set_project_metakey(project, user, metaKeyAndValue.get_meta_key(\
						MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id),
						MetaKeyAndValue.META_VALUE_Success, meta_project.description)
		except:
			## finished with error
			process_SGE.set_process_controler(user, process_controler.get_name_project(project), ProcessControler.FLAG_ERROR)
			return
		
		### finished
		process_SGE.set_process_controler(user, process_controler.get_name_project(project), ProcessControler.FLAG_FINISHED)
	

	def create_graph_minor_variants(self, project, user):
		"""
		OBSOLETE
		Create a html graph with minor variants data
		return: html and png file
		"""
		manageDatabase = ManageDatabase()
		vect_data_sample = []
		
		for project_sample in project.project_samples.all():
			if (not project_sample.get_is_ready_to_proccess()): continue
			
			meta_data = manageDatabase.get_project_sample_metakey_last(project_sample, MetaKeyAndValue.META_KEY_Count_Hits, MetaKeyAndValue.META_VALUE_Success)
			if (meta_data == None): continue
			decodeCoverage = DecodeObjects()
			count_hits = decodeCoverage.decode_result(meta_data.description)
			vect_data_sample.append([project_sample.sample.name, count_hits.get_total(), count_hits.get_hits_less_50(), count_hits.get_hits_50_90()])
		
		### there's no data to create the graph
		if (len(vect_data_sample) == 0):
			meta_data = manageDatabase.set_project_metakey(project, user, MetaKeyAndValue.META_KEY_Count_Samples_Var_Graph, MetaKeyAndValue.META_VALUE_Success, '0') 
			return None
		
		### sort the data by bigger
		##	in: [[1, 2, 4], [2, 34, 5], [23, 54, 65]]
		##  in: sorted(a, key=lambda tup: (tup[0]), reverse=True)
		##  out: [[23, 54, 65], [2, 34, 5], [1, 2, 4]]
		vect_data_sample = sorted(vect_data_sample, key=lambda tup: (tup[1]))
		
		### create the graph
		(temp_file_html, temp_file_png) = self.create_graph_plotly(vect_data_sample)
		
		## set the meta_Key
		meta_data = manageDatabase.set_project_metakey(project, user, MetaKeyAndValue.META_KEY_Count_Samples_Var_Graph, MetaKeyAndValue.META_VALUE_Success,	
									'{}'.format(len(vect_data_sample)))
		return (temp_file_html, temp_file_png)
	
	
	def create_graph_plotly(self, vect_data_sample):
		"""
		create a file html with data from variation information
		# https://plot.ly/python/horizontal-bar-charts/
		"""
		trace1 = go.Bar(
		    y=[vect_[0] for vect_ in vect_data_sample],
		    x=[vect_[2] for vect_ in vect_data_sample],
		    name='Var. <50',
		    orientation = 'h',
		    marker = dict(
		        color = 'rgba(50, 171, 96, 0.6)',
		        line = dict(
		            color = 'rgba(50, 171, 96, 1.0)',
		            width = 3)
		    )
		)
		trace2 = go.Bar(
		    y=[vect_[0] for vect_ in vect_data_sample],
		    x=[vect_[3] for vect_ in vect_data_sample],
		    name='50<Var.<90',
		    orientation = 'h',
		    marker = dict(
		        color = 'rgba(50, 96, 171, 0.6)',
		        line = dict(
		            color = 'rgba(50, 96, 171, 1.0)',
		            width = 3)
		    )
		)
		
		data = [trace1, trace2]
		layout = go.Layout(
		    barmode='stack',
		)
		
		fig = go.Figure(data=data, layout=layout)
		## save static image
		temp_file_png = None
		#temp_file_png = self.utils.get_temp_file("minor_var_graph", '.png')
		#py.image.save_as(fig, temp_file_png)	## online method, need to pay
		
		## save in a html file
		temp_file_html = self.utils.get_temp_file("minor_var_graph", '.html')
		plot(fig, filename=temp_file_html, auto_open=False, image_width=400, image_height=300, show_link=False)
		return (temp_file_html, temp_file_png)

	def calculate_global_files(self, type_file, project, user):
		"""
		Collect extra files
 		"""
		out_file = None
		out_file_file_system = None
		if (type_file == Project.PROJECT_FILE_NAME_COVERAGE):
			## collect coverage file for all samples
			out_file = self.create_coverage_file(project, user)
			out_file_file_system = project.get_global_file_by_project(TypePath.MEDIA_ROOT, type_file)
			
		elif (type_file == Project.PROJECT_FILE_NAME_TAB_VARIATIONS_SNIPPY):
			## snippy variations
			out_file = self.collect_variations_snippy(project)
			out_file_file_system = project.get_global_file_by_project(TypePath.MEDIA_ROOT, type_file)
		
		elif (type_file == Project.PROJECT_FILE_NAME_TAB_VARIATIONS_FREEBAYES):
			## freebayes remove >50%
			vect_type_remove = ['ins', 'del']
			out_file = self.collect_variations_freebayes(project, vect_type_remove)
			out_file_file_system = project.get_global_file_by_project(TypePath.MEDIA_ROOT, type_file)
		elif (type_file == Project.PROJECT_FILE_NAME_TAB_VARIATIONS_FREEBAYES_with_snps_indels):
			## freebayes remove >50%, keep snp and indel
			vect_type_remove = []
			out_file = self.collect_variations_freebayes(project, vect_type_remove)
			out_file_file_system = project.get_global_file_by_project(TypePath.MEDIA_ROOT, type_file)
		elif (type_file == Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV):
			## samples csv
			out_file = self.collect_sample_table(project, Constants.SEPARATOR_COMMA)
			out_file_file_system = project.get_global_file_by_project(TypePath.MEDIA_ROOT, type_file)
			
		elif (type_file == Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV_simple):
			## samples csv, simple version
			b_simple = True		## doesn't have software information and fastq1, fastq2
			out_file = self.collect_sample_table(project, Constants.SEPARATOR_COMMA, b_simple)
			out_file_file_system = project.get_global_file_by_project(TypePath.MEDIA_ROOT, type_file)
		
		elif (type_file == Project.PROJECT_FILE_NAME_SAMPLE_RESULT_TSV):
			## samples tsv
			out_file = self.collect_sample_table(project, Constants.SEPARATOR_TAB)
			out_file_file_system = project.get_global_file_by_project(TypePath.MEDIA_ROOT, type_file)
			
		elif (type_file == Project.PROJECT_FILE_NAME_SAMPLE_RESULT_json):
			## tree json
			out_file = self.create_json_file_from_sample_csv(project)
			out_file_file_system = project.get_global_file_by_project(TypePath.MEDIA_ROOT, type_file)
		
		elif (type_file == Project.PROJECT_FILE_NAME_SAMPLE_RESULT_all_consensus):
			out_file = self.merge_all_consensus_files(project)
			out_file_file_system = project.get_global_file_by_project(TypePath.MEDIA_ROOT, type_file)
		elif (type_file == Project.PROJECT_FILE_NAME_SAMPLE_mask_all_consensus):
			self.mask_all_consensus_files(project)	## mask all consensus for all projects, defined by user
			
		if (not out_file is None):
			self.utils.copy_file(out_file, out_file_file_system)
			self.utils.remove_file(out_file)
		elif (not out_file_file_system is None and os.path.exists(out_file_file_system)): self.utils.remove_file(out_file_file_system)


	def create_json_file_from_sample_csv(self, project):
		"""
		Create JSON file to insaPhylo
		"""
		
		file_name_root_sample = project.get_global_file_by_project(TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV_simple)
		if (os.path.exists(file_name_root_sample)):
			out_file = self.utils.get_temp_file('json_sample_file', FileExtensions.FILE_JSON)
			with open(out_file, 'w', encoding='utf-8') as handle_write, open(file_name_root_sample) as handle_in_csv:
				reader = csv.DictReader(handle_in_csv)
				all_data = json.loads(json.dumps(list(reader)))
				dt_result = {}
				for dict_data in all_data:
					if ('id' in dict_data):
						dt_out = dict_data.copy()
						del dt_out['id']
						dt_result[dict_data['id']] = dt_out
				if len(dt_result) == len(all_data):
					handle_write.write(json.dumps(dt_result))
				else:
					os.unlink(out_file)
					self.logger_production.error('ProjectID: {}  different number of lines processing Sample {} -> JSON {}'.format(project.id, len(dt_result), len(all_data)))
					self.logger_debug.error('ProjectID: {}  different number of lines processing Sample {} -> JSON {}'.format(project.id, len(dt_result), len(all_data)))
					return None
			return out_file
		else:
			self.logger_production.error('Sample csv file does not exist: {}'.format(file_name_root_sample))
			self.logger_debug.error('Sample csv file does not exist {}'.format(file_name_root_sample))
		return None


	def create_coverage_file(self, project, user):
		"""
		collect all coverage and make a file
		"""
		decode_coverage = DecodeObjects()
		manageDatabase = ManageDatabase()
		geneticElement = self.utils.get_elements_and_cds_from_db(project.reference, user)
		if (geneticElement == None): return None
		
		vect_ratios = ["% of size covered by at least 1-fold", "% of size covered by at least X-fold"]
		vect_reference = geneticElement.get_sorted_elements()
		out_file = self.utils.get_temp_file('coverage_file', FileExtensions.FILE_TSV)
		n_count = 0
		with open(out_file, "w", newline='') as output_file_handle_csv:
			csv_writer = csv.writer(output_file_handle_csv, delimiter=Constants.SEPARATOR_TAB, quotechar='"', quoting=csv.QUOTE_ALL)
			
			### write headers
			csv_writer.writerow([])
			csv_writer.writerow(['Chromosome'])
			vect_out = ['Name']
			vect_out.extend(vect_reference)
			csv_writer.writerow(vect_out)
			
			## length chromosome
			vect_out = ['Length']
			for element_name in vect_reference:
				vect_out.append(geneticElement.get_size_element(element_name))
			csv_writer.writerow(vect_out)
			
			csv_writer.writerow([])
			csv_writer.writerow([])
			vect_out = ['', 'Mean depth of coverage']
			vect_out.extend([''] * len(vect_reference))
			for ratio in vect_ratios: 
				vect_out.append(ratio)
				vect_out.extend([''] * len(vect_reference))
			csv_writer.writerow(vect_out)
			
			## chr names
			vect_out = ['Samples']
			vect_out.extend(vect_reference)
			for _, ratio in enumerate(vect_ratios): 
				vect_out.append('')
				if _ == len(vect_ratios) - 1: vect_out.append("X-Fold")
				vect_out.extend(vect_reference)
			csv_writer.writerow(vect_out)
			
			for project_sample in project.project_samples.all():
				if (not project_sample.get_is_ready_to_proccess()): continue
				vect_out = [project_sample.sample.name]
				
				meta_data = manageDatabase.get_project_sample_metakey_last(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
				if (meta_data == None): continue
				coverage = decode_coverage.decode_result(meta_data.description)
				
				for element_name in vect_reference:
					vect_out.append(coverage.get_coverage(element_name, Coverage.COVERAGE_ALL))
				
				vect_out.append('')
				for element_name in vect_reference:
					vect_out.append(coverage.get_coverage(element_name, Coverage.COVERAGE_MORE_0))
				
				### last coverage X-Fold
				vect_out.append('')
				vect_out.append(coverage.get_middle_limit())
				for element_name in vect_reference:
					vect_out.append(coverage.get_coverage_by_middle_tag(element_name))
				csv_writer.writerow(vect_out)
				n_count += 1
		if (n_count == 0):
			os.unlink(out_file)
			return None
		return out_file


	def collect_variations_snippy(self, project):
		"""
		collect snippy variations
		
		add flag if to snippy when coverage is yellow or red add column 'VARIANTS IN INCOMPLETE LOCUS' Yes or empty
		"""
		
		out_file = self.utils.get_temp_file('variations_snippy', FileExtensions.FILE_TSV)
		parse_out_files = ParseOutFiles()
		n_count = 0
		vect_type_out = []
		with open(out_file, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=Constants.SEPARATOR_TAB, quotechar='"', quoting=csv.QUOTE_ALL)
			for project_sample in project.project_samples.all():
				if (not project_sample.get_is_ready_to_proccess()): continue
				
				### different type of SNVs 
				if (project_sample.is_sample_illumina()):
					tab_file_to_process = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_SNIPPY_name)
				elif (project_sample.is_sample_ont()):
					tab_file_to_process = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_Medaka_name)
				else:
					continue
				if (not os.path.exists(tab_file_to_process)): continue
				parse_out_files.parse_tab_files_snippy(project_sample.sample.name, tab_file_to_process, csv_writer, vect_type_out, True if n_count == 0 else False)
				n_count += 1
		if (n_count == 0):
			os.unlink(out_file)
			return None
		return out_file
	
	def mask_all_consensus_files(self, project):
		"""
		merge all consensus files
		"""
		### read masking consensus
		manageDatabase = ManageDatabase()
		decode_masking_consensus = DecodeObjects()
		meta_value = manageDatabase.get_project_metakey_last(project, MetaKeyAndValue.META_KEY_Masking_consensus, MetaKeyAndValue.META_VALUE_Success)
		masking_consensus_original = None
		if not meta_value is None:
			masking_consensus_original = decode_masking_consensus.decode_result(meta_value.description)
		
		### if exist masking data
		if not masking_consensus_original is None and masking_consensus_original.has_masking_data():
			for project_sample in project.project_samples.all():
				if (not project_sample.get_is_ready_to_proccess() or project_sample.is_deleted): continue
				if not os.path.exists(project_sample.get_consensus_file(TypePath.MEDIA_ROOT)): continue

				""" check if exist backup """
				if (not os.path.exists(project_sample.get_backup_consensus_file())):
					self.utils.copy_file(project_sample.get_consensus_file(TypePath.MEDIA_ROOT),
										project_sample.get_backup_consensus_file())

				vect_record_out = []
				## always work with the backup	
				with open(project_sample.get_backup_consensus_file(), "rU") as handle_fasta:
					for record in SeqIO.parse(handle_fasta, "fasta"):
						masking_consensus = masking_consensus_original.dt_elements_mask.get(record.id, MaskingConsensus())
						if masking_consensus.has_data():
							vect_record_out.append(self.utils.mask_sequence(record,
								masking_consensus.mask_sites, masking_consensus.mask_from_beginning, 
								masking_consensus.mask_from_ends, masking_consensus.mask_regions))
						else: vect_record_out.append(record)
				if (len(vect_record_out) > 0):
					temp_file = self.utils.get_temp_file("masked_seq_", ".fasta")
					with open(temp_file, "w") as handle_fasta_out:
						SeqIO.write(vect_record_out, handle_fasta_out, "fasta")

					### move temp consensus to original position, if has info
					if os.stat(temp_file).st_size > 0:
						self.utils.move_file(temp_file, project_sample.get_consensus_file(TypePath.MEDIA_ROOT))
					else: os.unlink(temp_file)
	
	def merge_all_consensus_files(self, project):
		"""
		merge all consensus files
		"""
		
		out_file = self.utils.get_temp_file('all_consensus', FileExtensions.FILE_FASTA)
		vect_to_process = []
		for project_sample in project.project_samples.all():
			if (not project_sample.get_is_ready_to_proccess()): continue
			vect_to_process.append([
				project_sample.get_consensus_file(TypePath.MEDIA_ROOT),\
				project_sample.sample.name, project_sample.id])
		
		### set the number sequences that passed 	
		project.number_passed_sequences = len(vect_to_process)
		project.save()
		
		self.utils.merge_fasta_files(vect_to_process, out_file)
		return out_file

	def collect_variations_freebayes(self, project, vect_type_remove):
		"""
		collect freebayes variations
		"""
		out_file = self.utils.get_temp_file('variations_freebayes', FileExtensions.FILE_TSV)
		parse_out_files = ParseOutFiles()
		n_count = 0
		vect_type_out = []
		#vect_type_remove = ['ins', 'del']	### evan if it's in complex,ins or del,snp
		with open(out_file, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=Constants.SEPARATOR_TAB, quotechar='"', quoting=csv.QUOTE_ALL)
			for project_sample in project.project_samples.all():
				if (not project_sample.get_is_ready_to_proccess()): continue
				tab_file_to_process = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)
				if (not os.path.exists(tab_file_to_process)): continue
				parse_out_files.parse_tab_files(project_sample.sample.name, tab_file_to_process, csv_writer, vect_type_out,\
								vect_type_remove, Project.PERCENTAGE_validated_minor_variants, True if n_count == 0 else False)
				n_count += 1
		if (n_count == 0):
			os.unlink(out_file)
			return None
		return out_file
	
	def collect_variations_freebayes_only_one_file(self, file_in, file_out, vect_type_remove):
		"""
		collect freebayes variations
		"""
		parse_out_files = ParseOutFiles()
		vect_type_out = []
		##vect_type_remove = ['ins', 'del']

		with open(file_out, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=Constants.SEPARATOR_TAB, quotechar='"', quoting=csv.QUOTE_ALL)
			count_print = parse_out_files.parse_tab_files(None, file_in, csv_writer, vect_type_out, vect_type_remove,
					Project.PERCENTAGE_validated_minor_variants, True)

		### remove file if there's no out sequences
		if (count_print == 0):
			if (os.path.exists(file_out)): os.unlink(file_out)
			return None
		return file_out


	def collect_sample_table(self, project, column_separator, b_simple = False):
		"""
		collect sample table
		column_separator : COMMA or TAB
		id,fastq1,fastq2,data set,vaccine status,week,onset date,collection date,lab reception date,latitude,longitude
		:param b_simple == True, doesn't have fastq1,fastq2 and all software information, created for the tree view
		
		"""
		manage_database = ManageDatabase()
		default_software = DefaultProjectSoftware()
		out_file = self.utils.get_temp_file('sample_out', FileExtensions.FILE_CSV if\
					column_separator == Constants.SEPARATOR_COMMA else FileExtensions.FILE_TSV)
		
		### test if exist pangolin lineage, if exist parse the file.
		parse_pangolin = ParsePangolinResult(project.get_global_file_by_project(TypePath.MEDIA_ROOT,
								Project.PROJECT_FILE_NAME_Pangolin_lineage))
		
		with open(out_file, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=column_separator, quotechar='"',
						quoting=csv.QUOTE_MINIMAL if column_separator == Constants.SEPARATOR_COMMA else csv.QUOTE_ALL)
			if (b_simple): vect_out_header = CollectExtraData.HEADER_SAMPLE_OUT_CSV_simple.split(',')
			else: vect_out_header = CollectExtraData.HEADER_SAMPLE_OUT_CSV.split(',')
			
			### extra tags
			vect_tags = self.get_tags_for_samples_in_projects(project)
			vect_out_header.extend(vect_tags)
			
			if parse_pangolin.has_data():
				vect_out_header.append("Lineage (Pangolin)")
				vect_out_header.append("Scorpio (Pangolin)") 
			if (not b_simple):
				vect_out_header.append("Technology")
				vect_out_header.append("Sample Downsized")
			
			### all information about the softwares
			if (not b_simple):
				### set software names and versions
				### array with software names
				(vect_tags_out_sample, vect_tags_out_project_sample, b_trimmomatic_stats, 
					b_nanostat_sats, dict_all_results) = self.get_list_of_softwares_used(project)
				
				star_stats_tag = 0
				if (b_trimmomatic_stats):
					star_stats_tag = len(vect_out_header)
					vect_out_header.extend([tag_name.replace(":", "")\
						for tag_name in SoftwareNames.SOFTWARE_TRIMMOMATIC_vect_info_to_collect]) 
				
				if (b_nanostat_sats):
					if (star_stats_tag == 0): star_stats_tag = len(vect_out_header)
					vect_out_header.extend([tag_name + " - ONT"\
						for tag_name in SoftwareNames.SOFTWARE_NANOSTAT_vect_info_to_collect])
				
				vect_out_header.extend(vect_tags_out_sample)
				vect_out_header.extend(vect_tags_out_project_sample)
				
				### set default parameters
				vect_out_header.extend(SoftwareNames.VECT_INSAFLU_PARAMETER)
				
				### write reference name for this project name on a single line
				if (not b_trimmomatic_stats and not b_nanostat_sats):
					vect_out = ["Reference name", project.reference.name] + [''] * (len(vect_out_header) - 2)
				else: 
					vect_out = ["Reference name:", project.reference.name] +\
						[''] * (star_stats_tag - 2)
					if (b_trimmomatic_stats): vect_out.extend(["Trimmo. stats."] * len(SoftwareNames.SOFTWARE_TRIMMOMATIC_vect_info_to_collect))
					if (b_nanostat_sats): vect_out.extend(["NanoStat stats."] * len(SoftwareNames.SOFTWARE_NANOSTAT_vect_info_to_collect))
					vect_out.extend([''] * (len(vect_tags_out_sample) + len(vect_tags_out_project_sample) +\
									len(SoftwareNames.VECT_INSAFLU_PARAMETER)))
				csv_writer.writerow(vect_out)
			
			### write header
			csv_writer.writerow(vect_out_header)
			
			### if file to need to have the reference under the header
			if (b_simple): csv_writer.writerow([project.reference.name] + [''] * (len(vect_out_header) - 1))
			
			n_count = 0
			for project_sample in project.project_samples.all():
				if (not project_sample.get_is_ready_to_proccess()): continue
				vect_out = [project_sample.sample.name]
				
				### don't include file names, don't go to show in the tree
				if (not b_simple):
					### fastq1
					vect_out.append(project_sample.sample.file_name_1)
					### fastq2
					vect_out.append(project_sample.sample.file_name_2 if project_sample.sample.file_name_2 != None and len(project_sample.sample.file_name_2) > 0 else '')
				### dataset
				vect_out.append(project_sample.sample.data_set.name if project_sample.sample.data_set != None else '')
				### vaccine status
				vect_out.append(project_sample.sample.vaccine_status.name if project_sample.sample.vaccine_status != None else '')
				### week
				vect_out.append(str(project_sample.sample.week) if project_sample.sample.week != None else '')
				### onset date
				vect_out.append(project_sample.sample.date_of_onset.strftime('%Y-%m-%d') if project_sample.sample.date_of_onset != None else '')
				### collection date
				vect_out.append(project_sample.sample.date_of_collection.strftime('%Y-%m-%d') if project_sample.sample.date_of_collection != None else '')
				### lab reception date
				vect_out.append(project_sample.sample.date_of_receipt_lab.strftime('%Y-%m-%d') if project_sample.sample.date_of_receipt_lab != None else '')
				### latitude
				vect_out.append(str(project_sample.sample.geo_local.coords[0]) if project_sample.sample.geo_local != None else '')
				### longitude
				vect_out.append(str(project_sample.sample.geo_local.coords[1]) if project_sample.sample.geo_local != None else '')
				### type_subtype
				vect_out.append(project_sample.sample.type_subtype if project_sample.sample.type_subtype != None else '')
				### mixedinfection
				vect_out.append(project_sample.mixed_infections.tag.name if project_sample.mixed_infections != None else '')

				### print extra informations
				query_set = TagNames.objects.filter(sample=project_sample.sample)
				for tag_name_to_test in vect_tags:
					b_print = False
					for tag_names in query_set:
						if (tag_names.tag_name.name == tag_name_to_test):
							vect_out.append(tag_names.value)
							b_print = True
							break
					if (not b_print): vect_out.append('')

				### information about the software
				if (not b_simple):
					### pangolin if exists, not simple
					if parse_pangolin.has_data():
						vect_out.append(parse_pangolin.get_value(project_sample.seq_name_all_consensus, ParsePangolinResult.KEY_LINEAGE))
						vect_out.append(parse_pangolin.get_value(project_sample.seq_name_all_consensus, ParsePangolinResult.KEY_SCORPIO))

					### print info about technology	
					vect_out.append(project_sample.get_type_technology())
					### downsized
					vect_out.append("True" if manage_database.is_sample_downsized(project_sample.sample) else "False")
					
					###  BEGIN info about software versions  ####
					if project_sample.id in dict_all_results:
					
						### trimmomatic stats	
						if (b_trimmomatic_stats):
							self._get_info_from_trimmomatic_stats(vect_out,
								project_sample.id, dict_all_results)
						### nanostat stats	
						if (b_nanostat_sats):
							self._get_info_from_nanostats_stats(vect_out,
								project_sample.id, dict_all_results)
						
						### samples NAME_sample
						self._get_info_for_software(MetaKeyAndValue.NAME_sample, vect_out, dict_all_results,
							project_sample.id, vect_tags_out_sample)
	
						### samples NAME_project_sample
						self._get_info_for_software(MetaKeyAndValue.NAME_project_sample, vect_out, dict_all_results,
							project_sample.id, vect_tags_out_project_sample)
						
					else: vect_out.extend([''] * (len(vect_tags_out_sample) +\
						len(vect_tags_out_project_sample) +\
						len(SoftwareNames.SOFTWARE_TRIMMOMATIC_vect_info_to_collect if b_trimmomatic_stats else 0) +\
						len(SoftwareNames.SOFTWARE_NANOSTAT_vect_info_to_collect if b_nanostat_sats else 0)))
					###		second project_sample
					###  END info about software versions  ####
	
					### BEGIN save global parameters
					for parameter_name in SoftwareNames.VECT_INSAFLU_PARAMETER:
						### 
						if (parameter_name == SoftwareNames.INSAFLU_PARAMETER_MASK_CONSENSUS_name):
							if (project_sample.is_mask_consensus_sequences):
								vect_out.append(default_software.get_mask_consensus_single_parameter(project_sample,\
									DefaultParameters.MASK_CONSENSUS_threshold,
									ConstantsSettings.TECHNOLOGY_illumina \
									if project_sample.is_sample_illumina() else ConstantsSettings.TECHNOLOGY_minion))
							else:
								vect_out.append("Not applied")
						else:
							vect_out.append("")
				else: ## simple, for the tree
					### pangolin if exists, simple
					if parse_pangolin.has_data():
						vect_out.append(parse_pangolin.get_value(project_sample.seq_name_all_consensus, ParsePangolinResult.KEY_LINEAGE))
						vect_out.append(parse_pangolin.get_value(project_sample.seq_name_all_consensus, ParsePangolinResult.KEY_SCORPIO))

				### END save global parameters
				csv_writer.writerow(vect_out)
				n_count += 1
			
			### if there are previous lines add masking consensus regions
			if (n_count > 0 and not b_simple):
				### save masking limits
				manageDatabase = ManageDatabase()
				decode_masking_consensus = DecodeObjects()
				meta_value = manageDatabase.get_project_metakey_last(project, MetaKeyAndValue.META_KEY_Masking_consensus, MetaKeyAndValue.META_VALUE_Success)
				masking_consensus_original = None
				if not meta_value is None:
					masking_consensus_original = decode_masking_consensus.decode_result(meta_value.description)
			
					masking_consensus = MaskingConsensus()
					csv_writer.writerow([])
					csv_writer.writerow(["Masking consensus regions"])
					csv_writer.writerow(["Elements"] + masking_consensus.get_vect_header())
					for element in masking_consensus_original.get_sorted_elements(): 
						vect_out = [element]
						vect_out.append(masking_consensus_original.dt_elements_mask[element].mask_sites)
						vect_out.append(masking_consensus_original.dt_elements_mask[element].mask_from_beginning)
						vect_out.append(masking_consensus_original.dt_elements_mask[element].mask_from_ends)
						vect_out.append(masking_consensus_original.dt_elements_mask[element].mask_regions)
						csv_writer.writerow(vect_out)

		if (n_count == 0):
			os.unlink(out_file)
			return None
		return out_file

	
	def _get_info_for_software(self, tag_to_search, vect_out, dict_all_results, project_sample_id, vect_tags_info):
		"""
		"""
		if tag_to_search in dict_all_results[project_sample_id]:
			for software_name in vect_tags_info:
				software_result = ""
				for result in dict_all_results[project_sample_id][tag_to_search]:
					software_result = result.get_software(software_name)
					if (len(software_result) > 0): break
				vect_out.append(software_result)
		else: vect_out.extend([''] * len(vect_tags_info))
	
	def _get_info_from_trimmomatic_stats(self, vect_out, project_sample_id, dict_all_results):
		"""
		return the stats of trimmomatic
		"""
		if MetaKeyAndValue.NAME_sample in dict_all_results[project_sample_id]:
			for result in dict_all_results[project_sample_id][MetaKeyAndValue.NAME_sample]:
				software = result.get_software_instance(SoftwareNames.SOFTWARE_TRIMMOMATIC_name)
				if (not software is None and not software.get_vect_key_values() is None):
					for stat_names in SoftwareNames.SOFTWARE_TRIMMOMATIC_vect_info_to_collect:
						b_found = False
						for key_value in software.get_vect_key_values():
							if (key_value.key == stat_names):
								vect_out.append(key_value.value)
								b_found = True
								break
						if (not b_found): vect_out.append("")
					return
				
		vect_out.extend([''] * len(SoftwareNames.SOFTWARE_TRIMMOMATIC_vect_info_to_collect))
		
	def _get_info_from_nanostats_stats(self, vect_out, project_sample_id, dict_all_results):
		"""
		return the stats of trimmomatic
		"""
		if MetaKeyAndValue.NAME_sample in dict_all_results[project_sample_id]:
			for result in dict_all_results[project_sample_id][MetaKeyAndValue.NAME_sample]:
				vect_soft = result.get_list_software_instance(SoftwareNames.SOFTWARE_NanoStat_name)
				if (not vect_soft is None and len(vect_soft) == 2):
					for stat_names in SoftwareNames.SOFTWARE_NANOSTAT_vect_info_to_collect:
						b_found = False
						for key_value in vect_soft[1].get_vect_key_values():
							if (key_value.key == stat_names):
								vect_out.append(key_value.value)
								b_found = True
								break
						if (not b_found): vect_out.append("")
					return
		vect_out.extend([''] * len(SoftwareNames.SOFTWARE_NANOSTAT_vect_info_to_collect))
					
		
	def get_tags_for_samples_in_projects(self, project):
		"""
		tags for samples in projects
		"""
		dict_tags_out = {}
		vect_tags_out = []
		for project_sample in project.project_samples.all():
			if (not project_sample.get_is_ready_to_proccess()): continue
			for tag_name in project_sample.sample.tag_names.all():
				if (tag_name.name in dict_tags_out): continue
				dict_tags_out[tag_name.name] = 1
				vect_tags_out.append(tag_name.name)
		return vect_tags_out

	def get_list_of_softwares_used(self, project):
		"""
		get list of software, in sample and sample_project
		out: [Software_name, Software_name_1, Software_name_2, ...]
				all results by id project_sample and key software
		"""
		manage_database = ManageDatabase()
		decode_result = DecodeObjects()
		
		dict_all_results = {}		### all results in for project_sample
		vect_tags_out_sample = []
		vect_tags_out_project_sample = []
		b_trimmomatic_tags = False	### test if has trimmomatic statistics
		b_nanostat_sats = False		### test if has NanoFilt statistics
		
		### START First, test pangolin version, if exist
		result_pangolin = None
		meta_sample = manage_database.get_project_metakey_last(project,
							MetaKeyAndValue.META_KEY_Identify_pangolin,\
							MetaKeyAndValue.META_VALUE_Success)
		if (not meta_sample is None):
			result_pangolin = decode_result.decode_result(meta_sample.description)
			for software_name in result_pangolin.get_all_software_names():
				if not software_name in vect_tags_out_project_sample:
					vect_tags_out_project_sample.append(software_name)
		## END pangolin
	
		for project_sample in project.project_samples.all():
			if (not project_sample.get_is_ready_to_proccess()): continue
			### return software name, illumina and ONT
			dt_out = {MetaKeyAndValue.NAME_sample: []}
			for technology in MetaKeyAndValue.VECT_TECHNOLOGIES_OUT_REPORT:
				if (MetaKeyAndValue.NAME_sample in MetaKeyAndValue.DICT_SOFTWARE_SHOW_IN_RESULTS[technology]):
					for keys_for_sample in MetaKeyAndValue.DICT_SOFTWARE_SHOW_IN_RESULTS[technology].get(MetaKeyAndValue.NAME_sample):
						meta_sample = manage_database.get_sample_metakey_last(project_sample.sample,
								keys_for_sample, MetaKeyAndValue.META_VALUE_Success)
						if (meta_sample is None): continue
						result = decode_result.decode_result(meta_sample.description)
						dt_out[MetaKeyAndValue.NAME_sample].append(result)
						
						if (not b_trimmomatic_tags and technology == MetaKeyAndValue.TECHNOLOGY_illumina):
							soft_desc = result.get_software_instance(SoftwareNames.SOFTWARE_TRIMMOMATIC_name)
							if (not soft_desc is None and not soft_desc.get_vect_key_values() is None and\
									len(soft_desc.get_vect_key_values()) > 0):
								b_trimmomatic_tags = True
						if (not b_nanostat_sats and technology == MetaKeyAndValue.TECHNOLOGY_ONT):
							soft_desc = result.get_software_instance(SoftwareNames.SOFTWARE_NanoStat_name)
							if (not soft_desc is None and not soft_desc.get_vect_key_values() is None and\
									len(soft_desc.get_vect_key_values()) > 0):
								b_nanostat_sats = True
							
						for software_name in result.get_all_software_names():
							if not software_name in vect_tags_out_sample: vect_tags_out_sample.append(software_name)
			
			### key of project sample
			dt_out[MetaKeyAndValue.NAME_project_sample] = []
			for technology in MetaKeyAndValue.VECT_TECHNOLOGIES_OUT_REPORT:
				if (MetaKeyAndValue.NAME_project_sample in MetaKeyAndValue.DICT_SOFTWARE_SHOW_IN_RESULTS[technology]):
					for keys_for_proj_sample in MetaKeyAndValue.DICT_SOFTWARE_SHOW_IN_RESULTS[technology].get(MetaKeyAndValue.NAME_project_sample):
						meta_sample = manage_database.get_project_sample_metakey_last(project_sample,
								keys_for_proj_sample, MetaKeyAndValue.META_VALUE_Success)
						if (meta_sample is None): continue
						result = decode_result.decode_result(meta_sample.description)
						dt_out[MetaKeyAndValue.NAME_project_sample].append(result)
						for software_name in result.get_all_software_names():
							if not software_name in vect_tags_out_project_sample: vect_tags_out_project_sample.append(software_name)

			## add pangolin if exists...
			if (not result_pangolin is None): dt_out[MetaKeyAndValue.NAME_project_sample].append(result_pangolin)
			## add everything to all	
			dict_all_results[project_sample.id] = dt_out
			
		return (vect_tags_out_sample, vect_tags_out_project_sample, b_trimmomatic_tags, b_nanostat_sats, dict_all_results)


	def calculate_count_variations(self, project):
		"""
		calculate global variations for a project
		"""
		data_out = []		## [[#<50, 50<var<90, sample name], [#<50, 50<var<90, sample name], ....]
		for project_sample in project.project_samples.all():
			if (project_sample.is_deleted): continue
			if (project_sample.is_error): continue
			if (not project_sample.is_finished): continue
			if (not project_sample.is_sample_illumina()): continue
			data_out.append([project_sample.count_variations.var_less_50, project_sample.count_variations.var_bigger_50_90, project_sample.sample.name])

		destination_file = project.get_global_file_by_project(TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_TOTAL_VARIATIONS)
		if (len(data_out) > 0):
			data_out = sorted(data_out, key=lambda data_temp: data_temp[0] + data_temp[1], reverse=True)
			
			out_file = self.utils.get_temp_file('count_variations', FileExtensions.FILE_TSV)
			with open(out_file, 'w', newline='') as handle_out:
				csv_writer = csv.writer(handle_out, delimiter=Constants.SEPARATOR_TAB, quotechar='"', quoting=csv.QUOTE_MINIMAL)
				vect_out = ['Sample', 'Less 50', 'Between 90 and 50']
				csv_writer.writerow(vect_out)
				for data_temp in data_out:
					csv_writer.writerow([data_temp[2], data_temp[0], data_temp[1]])
			### move file
			self.utils.copy_file(out_file, destination_file)
			self.utils.remove_file(out_file);
			return destination_file
		else: self.utils.remove_file(destination_file)
		return None


class ParsePangolinResult(object):
	
	utils = Utils()
	HEADER_PANGOLIN_file = "taxon,lineage"
	KEY_LINEAGE = 'lineage'
	KEY_SCORPIO = 'scorpio_call'
	VECT_CALL = [KEY_LINEAGE, KEY_SCORPIO]
	
	def __init__(self, file_pangolin_output):
		self.file_pangolin_output = file_pangolin_output
		self.dt_data = {}	### SAmple : [lineage, probability, status]
		self.dt_header = {}
		self.process_file()
		
	def has_data(self):
		return len(self.dt_data) > 0

	def get_value(self, sample_name_starts_with, value_to_call):
		if (sample_name_starts_with is None): return ""
		vect_match = [key for key in self.dt_data.keys() if key.startswith(sample_name_starts_with)]
		return ";".join([self.dt_data.get(key, [""]).get(value_to_call, "") for key in vect_match ] )
	
	def process_file(self):
		"""
		### pangolin ouput
		taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,version,pangolin_version,pangoLEARN_version,pango_version,status,note
		ULSLA_1183718_2021_S16__MN908947,AY.4,0.0,1.0,Delta (AY.4-like),0.957400,0.000000,PLEARN-v1.2.86,3.1.15,2021-10-13,v1.2.86,passed_qc,scorpio call: Altalleles 45; Ref alleles 0; Amb alleles 2; Oth alleles 0
		ULSLA_1183440_2021_S15__MN908947,B.1.617.2,0.0,1.0,Delta (B.1.617.2-like),1.000000,0.000000,PLEARN-v1.2.86,3.1.15,2021-10-13,v1.2.86,passed_qc,scorpiocall: Alt alleles 13; Ref alleles 0; Amb alleles 0; Oth alleles 0
		"""

		## start parsing
		dt_positions = {}
		if (os.path.exists(self.file_pangolin_output)):
			vect_data = self.utils.read_text_file(self.file_pangolin_output)
			b_header = False
			for line in vect_data:
				if not b_header:
					if line.startswith(ParsePangolinResult.HEADER_PANGOLIN_file): b_header = True
					for i, value in enumerate(line.split(',')):
						if value.lower() in ParsePangolinResult.VECT_CALL: self.dt_header[value] = i
				else:	## start processing pangolin data
					lst_data = line.split(',')
					if len(lst_data) > 4:
						dt_value = {}
						for value in self.dt_header:
							dt_value[value] = lst_data[self.dt_header[value]]
						self.dt_data[lst_data[0]] = dt_value
		
		
					