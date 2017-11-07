'''
Created on Nov 5, 2017

@author: mmp
'''
import unittest, math, time 
from django_q.tasks import async, result, fetch
from django_q.humanhash import humanize
from django_q.cluster import Cluster

class Test(unittest.TestCase):

	def setUp(self):
		self.c = Cluster()
		self.c.start()

	def tearDown(self):
		self.c.stop()
	
	def testDjangoQ(self):
		task_id = async(math.copysign, 2, -2)
		## print(humanize(task_id))
		while True:
			task_result = result(task_id)
			if (task_result is None):
				time.sleep(7)
			else: 
				self.assertTrue('-2.0', str(task_result))
				break
		
			
	def testDjangoSyncQ(self):
		# create a synchronous task
		task_id = async(math.copysign, 2, -2, sync=True)
		# the task will then be available immediately
		task = fetch(task_id)
		# and can be examined
		if not task.success: self.fail('An error occurred: {}'.format(task.result))
		
		task_result = result(task_id)
		self.assertTrue('-2.0', str(task_result))
		
		
