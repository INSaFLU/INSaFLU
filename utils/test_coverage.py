'''
Created on Nov 21, 2017

@author: mmp
'''
import unittest
from utils.coverage import Coverage
import os, filecmp
from django.conf import settings 
from utils.constantsTestsCase import ConstantsTestsCase

class Test(unittest.TestCase):


	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass
	
	def test_coverage(self):
		
		coverage = Coverage()
		vect_coverage = [1,2,4,5,6,7]
		vect_coverage.extend([2] * 100)
		vect_coverage.extend([20] * 100)
		vect_coverage.extend([200] * 100)
		vect_coverage.extend([500] * 100)
		vect_coverage.extend([300] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([330] * 100)
		vect_coverage.extend([3] * 100)
		vect_coverage.extend([33] * 100)
		vect_genes = [[10, 300, 'xpto', 1], [500, 800, 'xpto2', -1]]
		var_more_50 = [2, 100, 500, 800]
		var_less_50 = [6, 200, 350]
		output_image = "/tmp/xpto_zpt.png"
		average_coverage = "23432"
		ratio_more_zero = "100"
		rati_more_nine = "98"
		sample_name = "xpto"
		sequence_name = "1"
		coverage.create_coverage(vect_coverage, vect_genes, var_more_50, var_less_50, output_image,
					average_coverage, ratio_more_zero, rati_more_nine, sample_name, sequence_name)

		self.assertTrue(os.path.exists(output_image))
		self.assertTrue(filecmp.cmp(output_image, os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_IMAGES, 'first.png' )))
		os.unlink(output_image)
