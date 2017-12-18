'''
Created on Nov 27, 2017

@author: mmp
'''

from utils.utils import Utils 
from managing_files.manage_database import ManageDatabase
from managing_files.models import Project
from constants.meta_key_and_values import MetaKeyAndValue
from utils.result import DecodeObjects
from constants.constants import TypePath, Constants
from utils.tree import CreateTree
import os, time
import plotly.graph_objs as go
from plotly.offline import plot
from django.db import transaction

class CollectExtraData(object):
	'''
	classdocs
	'''

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
		manageDatabase = ManageDatabase()
		
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

		### create trees
		createTree = CreateTree()
		createTree.create_tree_and_alignments(project, user)
		
		meta_project = manageDatabase.get_project_metakey_last(project, metaKeyAndValue.get_meta_key(\
					MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id), MetaKeyAndValue.META_VALUE_Queue)
		if (meta_project != None):
			manageDatabase.set_project_metakey(project, user, metaKeyAndValue.get_meta_key(\
					MetaKeyAndValue.META_KEY_Queue_TaskID_Project, project.id), MetaKeyAndValue.META_VALUE_Success, meta_project.description)

	
	def create_graph_minor_variants(self, project, user):
		"""
		OBSOLETE
		Create a html graph with minor variants data
		return: html and png file
		"""
		manageDatabase = ManageDatabase()
		vect_data_sample = []
		
		for project_sample in project.project_sample.all():
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

