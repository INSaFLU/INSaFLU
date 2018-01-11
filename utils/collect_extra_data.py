'''
Created on Nov 27, 2017

@author: mmp
'''

from utils.utils import Utils 
from managing_files.manage_database import ManageDatabase
from managing_files.models import Project, ProjectSample, TagNames
from constants.meta_key_and_values import MetaKeyAndValue
from utils.result import DecodeObjects
from constants.constants import TypePath, Constants, FileType, FileExtensions
from constants.software_names import SoftwareNames
from utils.tree import CreateTree
import os, time, csv
import plotly.graph_objs as go
from plotly.offline import plot
from django.db import transaction
from utils.result import Coverage
from utils.parse_out_files import ParseOutFiles

class CollectExtraData(object):
	'''
	classdocs
	'''

	HEADER_SAMPLE_OUT_CSV = "id,fastq1,fastq2,data set,vaccine status,week,onset date,collection date,lab reception date,latitude,longitude,type-subtype,putative mixed-infection"
	SEPARATOR_COMMA = ','
	SEPARATOR_TAB = '\t'
	
	utils = Utils()
	
	def __init__(self):
		'''
		Constructor
		'''
		pass
	
	@transaction.atomic
	def collect_extra_data_for_project(self, project, user, vect_taskID):
		"""
		Everything that is necessary to do in the project
		Collect all extra data after all samples are finished
		only run after all the vect_taskID are finished
		"""
		### get the taskID and seal it
		metaKeyAndValue = MetaKeyAndValue()
		manage_database = ManageDatabase()
		
		if (vect_taskID != None and len(vect_taskID) > 0):
			while not self.utils.is_all_tasks_finished(vect_taskID):
				time.sleep(Constants.WAIT_TIME_TASKS_FINISHED)
	
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
		
		### create trees
		createTree = CreateTree()
		createTree.create_tree_and_alignments(project, user)
		
		meta_project = manage_database.get_project_metakey_last(project, metaKeyAndValue.get_meta_key(\
					MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id), MetaKeyAndValue.META_VALUE_Queue)
		if (meta_project != None):
			manage_database.set_project_metakey(project, user, metaKeyAndValue.get_meta_key(\
					MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id), MetaKeyAndValue.META_VALUE_Success, meta_project.description)

	
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
			out_file = self.collect_sample_table(project, CollectExtraData.SEPARATOR_COMMA)
			out_file_file_system = project.get_global_file_by_project(TypePath.MEDIA_ROOT, type_file)
		
		elif (type_file == Project.PROJECT_FILE_NAME_SAMPLE_RESULT_TSV):
			## samples tsv
			out_file = self.collect_sample_table(project, CollectExtraData.SEPARATOR_TAB)
			out_file_file_system = project.get_global_file_by_project(TypePath.MEDIA_ROOT, type_file)
		
		if (out_file != None):
			self.utils.copy_file(out_file, out_file_file_system)
			os.unlink(out_file)
		elif (out_file_file_system != None and os.path.exists(out_file_file_system)): os.unlink(out_file_file_system)


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
		with open(out_file, "w") as output_file_handle:

			### write headers
			output_file_handle.write("\nChromosome\nName\t" + "\t".join(vect_reference) + "\nLength")
			for element_name in vect_reference:
				output_file_handle.write("\t{}".format(geneticElement.get_size_element(element_name)))
			
			output_file_handle.write("\n\nMean depth of coverage\t" + "\t" * len(vect_reference))
			print("\t" * len(vect_reference))
			for ratio in vect_ratios: output_file_handle.write("\t%s" % (ratio) + "\t" * len(vect_reference))
			output_file_handle.write("\n")
			
			for project_sample in project.project_samples.all():
				if (not project_sample.get_is_ready_to_proccess()): continue
				sz_out = project_sample.sample.name
				
				meta_data = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
				if (meta_data == None): continue
				coverage = decode_coverage.decode_result(meta_data.description)
				
				for element_name in vect_reference:
					sz_out += "\t{}".format(coverage.get_coverage(element_name, Coverage.COVERAGE_ALL))
				
				sz_out += "\t"
				for element_name in vect_reference:
					sz_out += "\t{}".format(coverage.get_coverage(element_name, Coverage.COVERAGE_MORE_0))
					
				sz_out += "\t"
				for element_name in vect_reference:
					sz_out += "\t{}".format(coverage.get_coverage(element_name, Coverage.COVERAGE_MORE_9))
				output_file_handle.write(sz_out + "\n")
				n_count += 1
		if (n_count == 0):
			os.unlink(out_file)
			return None
		return out_file


	def collect_variations_snippy(self, project):
		"""
		collect snippy variations
		"""
		
		out_file = self.utils.get_temp_file('variations_snippy', FileExtensions.FILE_TSV)
		parse_out_files = ParseOutFiles()
		n_count = 0
		vect_type_out = ['snp', 'del', 'ins']
		with open(out_file, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=CollectExtraData.SEPARATOR_TAB, quotechar='"', quoting=csv.QUOTE_ALL)
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
		vect_type_out = ['snp']
		with open(out_file, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=CollectExtraData.SEPARATOR_TAB, quotechar='"', quoting=csv.QUOTE_ALL)
			for project_sample in project.project_samples.all():
				if (not project_sample.get_is_ready_to_proccess()): continue
				tab_file_to_process = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)
				if (not os.path.exists(tab_file_to_process)): continue
				parse_out_files.parse_tab_files(project_sample.sample.name, tab_file_to_process, csv_writer, vect_type_out, 50, True if n_count == 0 else False)
				n_count += 1
		if (n_count == 0):
			os.unlink(out_file)
			return None
		return out_file


	def collect_sample_table(self, project, column_separator):
		"""
		collect sample table
		column_separator : COMMA or TAB
		id,fastq1,fastq2,data set,vaccine status,week,onset date,collection date,lab reception date,latitude,longitude
		"""
		
		out_file = self.utils.get_temp_file('sample_out', FileExtensions.FILE_CSV if\
					column_separator == CollectExtraData.SEPARATOR_COMMA else FileExtensions.FILE_TSV)
		
		with open(out_file, 'w', newline='') as handle_out:
			csv_writer = csv.writer(handle_out, delimiter=column_separator, quotechar='"',
						quoting=csv.QUOTE_MINIMAL if column_separator == CollectExtraData.SEPARATOR_COMMA else csv.QUOTE_ALL)
			vect_out = CollectExtraData.HEADER_SAMPLE_OUT_CSV.split(',')
			
			### extra tags
			vect_tags = self.get_tags_for_samples_in_projects(project)
			vect_out.extend(vect_tags)
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

				### print extra informatios
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
			
			
			