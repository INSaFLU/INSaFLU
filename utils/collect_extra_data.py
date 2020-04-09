'''
Created on Nov 27, 2017

@author: mmp
'''

from utils.utils import Utils 
from managing_files.manage_database import ManageDatabase
from managing_files.models import Project, TagNames, ProcessControler
from constants.meta_key_and_values import MetaKeyAndValue
from utils.result import DecodeObjects
from constants.constants import TypePath, Constants, FileType, FileExtensions
from constants.software_names import SoftwareNames
from utils.tree import CreateTree
import os, csv, time, json, logging
import plotly.graph_objs as go
from plotly.offline import plot
from django.db import transaction
from utils.result import Coverage
from utils.parse_out_files import ParseOutFiles
from utils.process_SGE import ProcessSGE
from django.conf import settings

class CollectExtraData(object):
	'''
	classdocs
	'''

	HEADER_SAMPLE_OUT_CSV = "id,fastq1,fastq2,data set,vaccine status,week,onset date,collection date,lab reception date,latitude,longitude,type-subtype,putative mixed-infection"

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
		
	def __collect_update_extra_metadata_for_project(self, project, user):
		"""
		Only for update metadata
		"""
		### get the taskID and seal it
		process_controler = ProcessControler()
		process_SGE = ProcessSGE()
		try:
			
			## collect sample table with plus type and subtype, mixed infection, equal to upload table
			self.calculate_global_files(Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV, project, user)
			self.calculate_global_files(Project.PROJECT_FILE_NAME_SAMPLE_RESULT_TSV, project, user)
			## IMPORTANT -> this need to be after of Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV
			self.calculate_global_files(Project.PROJECT_FILE_NAME_SAMPLE_RESULT_json, project, user)
		except:
			## finished with error
			process_SGE.set_process_controler(user, process_controler.get_name_project(project), ProcessControler.FLAG_ERROR)
			return
		
		### finished
		process_SGE.set_process_controler(user, process_controler.get_name_project(project), ProcessControler.FLAG_FINISHED)
		
#	@transaction.atomic
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
			
			### calculate global file
			self.calculate_global_files(Project.PROJECT_FILE_NAME_COVERAGE, project, user)
			## collect tab variations snippy
			self.calculate_global_files(Project.PROJECT_FILE_NAME_TAB_VARIATIONS_SNIPPY, project, user)
			## collect tab variations freebayes, <50%
			## remove del or ins
			self.calculate_global_files(Project.PROJECT_FILE_NAME_TAB_VARIATIONS_FREEBAYES, project, user)
			## collect sample table with plus type and subtype, mixed infection, equal to upload table
			self.calculate_global_files(Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV, project, user)
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
						MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id), MetaKeyAndValue.META_VALUE_Success, meta_project.description)
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
			
			meta_data = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Count_Hits, MetaKeyAndValue.META_VALUE_Success)
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
			## freebayes <50
			out_file = self.collect_variations_freebayes(project)
			out_file_file_system = project.get_global_file_by_project(TypePath.MEDIA_ROOT, type_file)
		
		elif (type_file == Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV):
			## samples csv
			out_file = self.collect_sample_table(project, Constants.SEPARATOR_COMMA)
			out_file_file_system = project.get_global_file_by_project(TypePath.MEDIA_ROOT, type_file)
		
		elif (type_file == Project.PROJECT_FILE_NAME_SAMPLE_RESULT_TSV):
			## samples tsv
			out_file = self.collect_sample_table(project, Constants.SEPARATOR_TAB)
			out_file_file_system = project.get_global_file_by_project(TypePath.MEDIA_ROOT, type_file)
			
		elif (type_file == Project.PROJECT_FILE_NAME_SAMPLE_RESULT_json):
			## tree json
			out_file = self.create_json_file_from_sample_csv(project)
			out_file_file_system = project.get_global_file_by_project(TypePath.MEDIA_ROOT, type_file)
		
		if (out_file != None):
			self.utils.copy_file(out_file, out_file_file_system)
			os.unlink(out_file)
		elif (out_file_file_system != None and os.path.exists(out_file_file_system)): os.unlink(out_file_file_system)


	def create_json_file_from_sample_csv(self, project):
		"""
		Create JSON file to insaPhylo
		"""
		
		file_name_root_sample = project.get_global_file_by_project(TypePath.MEDIA_ROOT, Project.PROJECT_FILE_NAME_SAMPLE_RESULT_CSV)
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
		
		vect_ratios = ["% of size covered by at least 1-fold", "% of size covered by at least 10-fold"]
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
			vect_out = ['Mean depth of coverage']
			vect_out.extend([''] * len(vect_reference))
			vect_out.append('')
			for ratio in vect_ratios: 
				vect_out.append(ratio)
				vect_out.extend([''] * len(vect_reference))
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
					
				vect_out.append('')
				for element_name in vect_reference:
					vect_out.append(coverage.get_coverage(element_name, Coverage.COVERAGE_MORE_9))
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
				tab_file_to_process = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_SNIPPY_name)
				if (not os.path.exists(tab_file_to_process)): continue
				parse_out_files.parse_tab_files_snippy(project_sample.sample.name, tab_file_to_process, csv_writer, vect_type_out, True if n_count == 0 else False)
				n_count += 1
		if (n_count == 0):
			os.unlink(out_file)
			return None
		return out_file


	def collect_variations_freebayes(self, project):
		"""
		collect freebayes variations
		"""
		out_file = self.utils.get_temp_file('variations_freebayes', FileExtensions.FILE_TSV)
		parse_out_files = ParseOutFiles()
		n_count = 0
		vect_type_out = []
		vect_type_remove = ['ins', 'del']	### evan if it's in complex,ins or del,snp
		with open(out_file, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=Constants.SEPARATOR_TAB, quotechar='"', quoting=csv.QUOTE_ALL)
			for project_sample in project.project_samples.all():
				if (not project_sample.get_is_ready_to_proccess()): continue
				tab_file_to_process = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)
				if (not os.path.exists(tab_file_to_process)): continue
				parse_out_files.parse_tab_files(project_sample.sample.name, tab_file_to_process, csv_writer, vect_type_out,\
								vect_type_remove, 50, True if n_count == 0 else False)
				n_count += 1
		if (n_count == 0):
			os.unlink(out_file)
			return None
		return out_file
	
	def collect_variations_freebayes_only_one_file(self, file_in, file_out):
		"""
		collect freebayes variations
		"""
		parse_out_files = ParseOutFiles()
		vect_type_out = []
		vect_type_remove = ['ins', 'del']

		with open(file_out, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=Constants.SEPARATOR_TAB, quotechar='"', quoting=csv.QUOTE_ALL)
			count_print = parse_out_files.parse_tab_files(None, file_in, csv_writer, vect_type_out, vect_type_remove, 50, True)

		### remove file if there's no out sequences
		if (count_print == 0):
			if (os.path.exists(file_out)): os.unlink(file_out)
			return None
		return file_out


	def collect_sample_table(self, project, column_separator):
		"""
		collect sample table
		column_separator : COMMA or TAB
		id,fastq1,fastq2,data set,vaccine status,week,onset date,collection date,lab reception date,latitude,longitude
		"""
		
		out_file = self.utils.get_temp_file('sample_out', FileExtensions.FILE_CSV if\
					column_separator == Constants.SEPARATOR_COMMA else FileExtensions.FILE_TSV)
		
		with open(out_file, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=column_separator, quotechar='"',
						quoting=csv.QUOTE_MINIMAL if column_separator == Constants.SEPARATOR_COMMA else csv.QUOTE_ALL)
			vect_out = CollectExtraData.HEADER_SAMPLE_OUT_CSV.split(',')
			
			### extra tags
			vect_tags = self.get_tags_for_samples_in_projects(project)
			vect_out.extend(vect_tags)
			csv_writer.writerow(vect_out)

			vect_out = [project.reference.name] + [''] * (len(vect_out) - 1)
			csv_writer.writerow(vect_out)
			n_count = 0
			for project_sample in project.project_samples.all():
				if (not project_sample.get_is_ready_to_proccess()): continue
				vect_out = [project_sample.sample.name]
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

				csv_writer.writerow(vect_out)
				n_count += 1
		if (n_count == 0):
			os.unlink(out_file)
			return None
		return out_file
	
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
			
	def calculate_count_variations(self, project):
		"""
		calculate global variations for a project
		"""
		data_out = []		## [[#<50, 50<var<90, sample name], [#<50, 50<var<90, sample name], ....]
		for project_sample in project.project_samples.all():
			if (project_sample.is_deleted): continue
			if (project_sample.is_error): continue
			if (not project_sample.is_finished): continue
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
			return destination_file
		else: os.unlink(destination_file)
		return None

