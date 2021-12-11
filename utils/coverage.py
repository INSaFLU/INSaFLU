'''
Created on Nov 21, 2017

@author: mmp
'''
from PIL import Image, ImageDraw, ImageFont
from constants.constants import FileType, TypePath, FileExtensions
from utils.result import DecodeObjects
from utils.parse_coverage_file import GetCoverage
from managing_files.manage_database import ManageDatabase
from managing_files.models import ProjectSample
from constants.meta_key_and_values import MetaKeyAndValue
from utils.utils import Utils
from utils.result import Coverage
import numpy, os, logging
from constants.software_names import SoftwareNames

class DrawAllCoverage(object):
	
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	def __init__(self):
		'''
		Constructor
		'''
		pass


	def draw_all_coverages(self, project_sample, software_name = SoftwareNames.SOFTWARE_SNIPPY_name):
		"""
		draw all coverage images for a coverage sample
		
		.vect_coverage, vect_genes, var_more_50, var_less_50, output_image,
					average_coverage, ratio_more_zero, ratio_more_nine, sample_name, sequence_name):
		"""
		utils = Utils()
		
		### get coverage vectors from deep file
		get_coverage = GetCoverage()
		coverage_file = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ, software_name)
		if (not os.path.exists(coverage_file)):
			self.logger_production.error("File doesn't exist: " + coverage_file)
			self.logger_debug.error("File doesn't exist: " + coverage_file)
			raise IOError("File doesn't exist: " + coverage_file)
		dict_coverage = get_coverage.get_dict_with_coverage(coverage_file)
		
		### get all elements and gene names
		geneticElement = utils.get_elements_and_genes(project_sample.project.reference.get_reference_gbk(TypePath.MEDIA_ROOT))
		
		### get positions of variations
		(dict_less_50, dict_more_50, dict_more_90) = ({}, {}, {})
		vect_count_type = ['snp']
		tab_file_from_freebayes = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)
		if (os.path.exists(tab_file_from_freebayes)):
			(dict_less_50, dict_more_50, dict_more_90) = utils.get_variations_by_freq_from_tab(tab_file_from_freebayes, vect_count_type)
		else:
			tab_file_from_medaka = project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_Medaka_name)
			(dict_less_50, dict_more_50, dict_more_90) = utils.get_variations_by_freq_from_tab(tab_file_from_medaka, vect_count_type)
		
		### get coverage
		manageDatabase = ManageDatabase()
		meta_value = manageDatabase.get_project_sample_metakey_last(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
		decode_coverage = DecodeObjects()
		coverage = decode_coverage.decode_result(meta_value.description)
		
		draw_coverage = DrawCoverage(coverage.limit_defined_by_user)
		for sequence_name in geneticElement.get_sorted_elements():
			draw_coverage.create_coverage(dict_coverage[sequence_name], geneticElement.get_genes(sequence_name),\
					dict_more_90.get(sequence_name, []),\
					dict_more_50.get(sequence_name, []),\
					dict_less_50.get(sequence_name, []),\
					project_sample.get_global_file_by_element(TypePath.MEDIA_ROOT, ProjectSample.PREFIX_FILE_COVERAGE,
					sequence_name.replace('/', '_'), FileExtensions.FILE_PNG),\
					coverage.get_coverage(sequence_name, Coverage.COVERAGE_ALL),\
					coverage.get_coverage(sequence_name, Coverage.COVERAGE_MORE_0),\
					coverage.get_coverage(sequence_name, Coverage.COVERAGE_MORE_DEFINED_BY_USER) if\
					coverage.is_exist_limit_defined_by_user() else\
					coverage.get_coverage(sequence_name, Coverage.COVERAGE_MORE_9),\
					project_sample.sample.name, sequence_name)

class DrawCoverage(object):
	'''
	classdocs
	'''
	PATH_FONT = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
	PATH_FONT_BOLD = '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'
	FONT_SIZE_BIG = 20
	FONT_SIZE_MEDIUM = 14
	FONT_SIZE_SMALL = 9
	
	WITHE_COLOR = (255, 255, 255, 0)
	
	GAP_BETWEEN_GENES = 17
	GAP_START_X = 20
	GAP_START_Y = 25
	START_DRAW_HEADER = 5
	DRAW_HEADER_Y = GAP_START_Y
	DRAW_COVERAGE_Y = DRAW_HEADER_Y + GAP_START_Y
	SIZE_COVERAGE_Y = 200
	SHRINK_GENE_BAR_UTR = 2
	SIZE_TOTAL_COVERAGE_Y = SIZE_COVERAGE_Y + GAP_START_Y	## plus the ruler...
	GAP_START_GENES_Y = SIZE_TOTAL_COVERAGE_Y + DRAW_COVERAGE_Y	+ 10 ## position to draw the exons, green part
	GAP_START_GENES_X = GAP_START_X + 20
	DRAW_HEIGHT_GENE_SQUARE_Y = 20
	DRAW_HEIGHT_GENE_STRAND_Y = 17
	DRAW_STRAND_Y = GAP_START_GENES_Y + DRAW_HEIGHT_GENE_SQUARE_Y + 5
	DRAW_BOTTOM_Y = GAP_START_GENES_Y + DRAW_HEIGHT_GENE_SQUARE_Y + DRAW_HEIGHT_GENE_STRAND_Y + (GAP_START_Y >> 1) + (START_DRAW_HEADER << 2) + 40
	DRAW_HEADER_BOTTOM = DRAW_STRAND_Y + (START_DRAW_HEADER << 2)
	LIMIT_MINUNUM_COVERAGE = 10
	DRAW_MINIMUN_LENGTH_X = 600
	GAP_MARK_VARIATIONS = 3
	GAP_MARK_VARIATIONS_BIG = 5
	
	COLOR_SEQUENCE = (0, 102, 0)
	COLOR_GENE = (0, 0, 153)
	COLOR_RGBGrey_64_64_64 = (64, 64, 64)
	COLOR_RGBGrey_32_32_32 = (32, 32, 32)
	COLOR_RGBGrey_225_225_225 = (225, 225, 225)
	COLOR_RGBGrey_225_225_225_255 = (225, 225, 225, 255)
	COLOR_RGBGrey_169_169_169_255 = (169, 169, 169, 255)
	COLOR_RGBGreen_0_153_0 = (0, 153, 0)
	COLOR_RGBGreen_0_102_0 = (0, 102, 0)
	COLOR_RGBRed_153_0_0 = (153, 0, 0)
	COLOR_RGBBlue_0_0_153 = (0, 0, 153)
	COLOR_RGBBlack = (0, 0, 0)

	rateImage = 2			## rate of the image
	rateDrawCoverage = 2	## rate to draw the coverage
	
	def __init__(self, limit_defined_by_user):
		'''
		Constructor
		'''
		self.limit_defined_by_user = self.LIMIT_MINUNUM_COVERAGE if limit_defined_by_user is None else limit_defined_by_user
		self.start_image_x = DrawCoverage.GAP_START_GENES_X
	
	def create_coverage(self, vect_coverage, vect_genes, var_more_90, var_more_50, var_less_50, output_image,
					average_coverage, ratio_more_zero, ratio_more_nine, sample_name, sequence_name):
		"""
		vect_genes = [Gene(name, pos_start, pos_end, strand), Gene(...), ...]
		var_more_90 = [pos, pos1, pos2, ...]
		var_more_50 = [pos, pos1, pos2, ...]
		var_less_50 = [pos, pos1, pos2, ...]
		"""
		## size of the image
		self.maxSizeImage = ((self.GAP_START_GENES_X << 1) + int((len(vect_coverage) / self.rateImage)), self.DRAW_BOTTOM_Y)
		
		self.start_image_x = (DrawCoverage.GAP_START_GENES_X + ((DrawCoverage.DRAW_MINIMUN_LENGTH_X - self.maxSizeImage[0]) >> 1))\
			if self.maxSizeImage[0] < DrawCoverage.DRAW_MINIMUN_LENGTH_X else DrawCoverage.GAP_START_GENES_X
		size_image_x = DrawCoverage.DRAW_MINIMUN_LENGTH_X if self.maxSizeImage[0] < DrawCoverage.DRAW_MINIMUN_LENGTH_X else self.maxSizeImage[0]
		im = Image.new("RGB", (size_image_x, self.maxSizeImage[1]), self.WITHE_COLOR)
		
		draw = ImageDraw.Draw(im) 
		self.draw_genes_and_coverage(draw, vect_coverage, vect_genes)
		self.draw_legend_coverage(draw, self.get_start_x(), self.get_start_x() + int(len(vect_coverage) / self.rateImage), int(len(vect_coverage) / self.rateImage))
		self.draw_variants(draw, var_more_50, var_less_50)
		(max_coverage, min_coverage) = self.get_max_min(vect_coverage)
		(position_first_header, position_second_header, position_third_header, position_x_fourth_header) =\
				self.draw_header(draw, 0, 0, 0, 0, average_coverage, ratio_more_zero, ratio_more_nine,\
				sample_name, sequence_name, max_coverage, min_coverage, True)		### get only the size
		self.draw_header(draw, 
				(size_image_x >> 1) - (position_first_header >> 1),\
				(size_image_x >> 1) - (position_second_header >> 1),\
				(size_image_x >> 1) - (position_third_header >> 1),\
				(size_image_x >> 1) - (position_x_fourth_header >> 1), 
				average_coverage, ratio_more_zero, ratio_more_nine, sample_name,\
				sequence_name, max_coverage, min_coverage, False)
		
		### test path 
		utils = Utils()
		utils.make_path(os.path.dirname(output_image))
		im.save(output_image)
		
	def get_start_x(self): return self.start_image_x
	
	def get_max_min(self, vect_coverage):
		(max_coverage, min_coverage) = (0, 99999999)
		for value in vect_coverage:
			if (value > max_coverage): max_coverage = value
			if (value < min_coverage): min_coverage = value
		return (max_coverage, min_coverage)


	def draw_genes_and_coverage(self, draw, vect_coverage, vect_genes):
		"""
		draw gene and coverage
		"""
		self.draw_coverage(draw, vect_coverage)
		
		color_squares = ColorSquares(draw)
		
		### paint main bar
		color_squares.set_color(self.COLOR_RGBBlue_0_0_153)
		color_squares.paint_bar((self.get_start_x(), self.GAP_START_GENES_Y + (self.SHRINK_GENE_BAR_UTR << 1) + 1),\
									(self.get_start_x() + int(len(vect_coverage) / self.rateImage),\
									(self.GAP_START_GENES_Y + self.DRAW_HEIGHT_GENE_SQUARE_Y) - (self.SHRINK_GENE_BAR_UTR << 1) - 1), False)
		for gene in vect_genes:
# 			y_start_top = self.GAP_START_GENES_Y - (self.SHRINK_GENE_BAR_UTR << 2) if gene.is_forward() else self.GAP_START_GENES_Y + (self.SHRINK_GENE_BAR_UTR << 3)
# 			y_start_bottom = (self.GAP_START_GENES_Y + self.DRAW_HEIGHT_GENE_SQUARE_Y) - (self.SHRINK_GENE_BAR_UTR << 2) if gene.is_forward() else (self.GAP_START_GENES_Y + self.DRAW_HEIGHT_GENE_SQUARE_Y) + (self.SHRINK_GENE_BAR_UTR << 2)
# 			color_squares.set_color(self.COLOR_RGBGreen_0_153_0)
# 			color_squares.paint_bar((int(gene.start / self.rateImage) + self.get_start_x(), y_start_top),\
# 									(self.get_start_x() + int(gene.end / self.rateImage),\
# 									y_start_bottom), False)
			
			self.draw_strand(draw, 1 if gene.is_forward() else -1, int(gene.start / self.rateImage) + self.get_start_x(),
							(int(gene.end / self.rateImage) + self.get_start_x()) - (int(gene.start / self.rateImage) + self.get_start_x()),
							(self.GAP_START_GENES_Y + self.DRAW_HEIGHT_GENE_SQUARE_Y) if gene.is_forward() else\
							self.GAP_START_GENES_Y )
			color_squares.draw_text(int(gene.start / self.rateImage) + self.get_start_x(),\
								self.get_start_x() + int(gene.end / self.rateImage),
								(self.GAP_START_GENES_Y + self.DRAW_HEIGHT_GENE_SQUARE_Y) if gene.is_forward() else\
								self.GAP_START_GENES_Y - self.DRAW_HEIGHT_GENE_SQUARE_Y, gene.name)
			
			
	def draw_strand(self, draw, direction, start_x, nLength_X, start_y):
		"""
		direction 1|-1
		"""
		draw.line((start_x, start_y, start_x + nLength_X, start_y), fill=self.COLOR_RGBBlack, width=1)
		if ((nLength_X >> 4) > 1):
			for i in range(0, nLength_X >> 4):
				slice_pos = (i * (nLength_X / (nLength_X >> 4)))
				draw.point((start_x + slice_pos, start_y - 1), fill=self.COLOR_RGBBlack)
				draw.point((start_x + slice_pos + direction, start_y - 2), fill=self.COLOR_RGBBlack)
				draw.point((start_x + slice_pos, start_y + 1), fill=self.COLOR_RGBBlack)
				draw.point((start_x + slice_pos + direction, start_y + 2), fill=self.COLOR_RGBBlack)
		else: ## draw in the middle
			draw.point((start_x + (nLength_X >> 1), start_y - 1), fill=self.COLOR_RGBBlack)
			draw.point((start_x + (nLength_X >> 1) + direction, start_y - 2), fill=self.COLOR_RGBBlack)
			draw.point((start_x + (nLength_X >> 1), start_y + 1), fill=self.COLOR_RGBBlack)
			draw.point((start_x + (nLength_X >> 1) + direction, start_y + 2), fill=self.COLOR_RGBBlack)


	def draw_coverage(self, draw, vect_coverage):
		"""
		draw all the coverage
		"""
		start_x = self.get_start_x()
		draw_y = 0
		draw_y_last = 0
		posYBottom = (self.DRAW_COVERAGE_Y + self.SIZE_COVERAGE_Y) - 1;
		bMoreThanLimit = False
		bMoreThanLimitLast = False
		bLessThanLimitLast = False
		vect_coverage_average = []
		for i in range(0, len(vect_coverage), self.rateImage):
			vect_coverage_average.append(int(numpy.mean(vect_coverage[i:i+self.rateImage])))
		
		draw.line((start_x - 1, self.DRAW_COVERAGE_Y, start_x - 1, self.DRAW_COVERAGE_Y + self.SIZE_COVERAGE_Y), fill=self.COLOR_RGBGrey_64_64_64, width=1)
		draw.line((start_x + len(vect_coverage_average), self.DRAW_COVERAGE_Y, start_x + len(vect_coverage_average), self.DRAW_COVERAGE_Y + self.SIZE_COVERAGE_Y),\
				fill=self.COLOR_RGBGrey_64_64_64, width=1)
		draw.line((start_x - 1, self.DRAW_COVERAGE_Y + self.SIZE_COVERAGE_Y + 1, start_x + len(vect_coverage_average), self.DRAW_COVERAGE_Y + self.SIZE_COVERAGE_Y + 1),\
				fill=self.COLOR_RGBGrey_64_64_64, width=1)
		for coverage in vect_coverage_average:
			bLessThanLimitLast = False
			bMoreThanLimit = False
			if (coverage >= (self.rateDrawCoverage * self.SIZE_COVERAGE_Y)):
				draw_y = self.DRAW_COVERAGE_Y
				bMoreThanLimit = True
			else:
				draw_y = self.DRAW_COVERAGE_Y + self.SIZE_COVERAGE_Y - (coverage / self.rateDrawCoverage)
				if (draw_y == (self.DRAW_COVERAGE_Y + self.SIZE_COVERAGE_Y)): draw_y = posYBottom;
				if (coverage <= self.limit_defined_by_user): bLessThanLimitLast = True;

			draw.line((start_x, posYBottom, start_x, draw_y), fill=self.COLOR_RGBGrey_225_225_225, width=1)
			if (draw_y_last == 0): draw.point((start_x, draw_y), fill=self.COLOR_RGBGrey_64_64_64)
			else: draw.line((start_x - 1, draw_y_last, start_x, draw_y), fill=self.COLOR_RGBGrey_64_64_64, width=1)
			
			### limits of coverage
			if (bMoreThanLimitLast and draw_y_last > 0):
				draw.line((start_x - 1, draw_y_last, start_x - 1, draw_y_last - 2), fill=self.COLOR_RGBGreen_0_153_0, width=1)
			elif (bLessThanLimitLast and draw_y_last >  0 and draw_y_last >= (self.DRAW_COVERAGE_Y + self.SIZE_COVERAGE_Y - (self.limit_defined_by_user / self.rateDrawCoverage))):
				draw.line((start_x - 1, draw_y_last, start_x - 1, draw_y_last - 2), fill=self.COLOR_RGBRed_153_0_0, width=1)
			bMoreThanLimitLast = bMoreThanLimit;
			draw_y_last = draw_y;
			start_x += 1;
		
		## exception for last point
		if (bMoreThanLimitLast and draw_y_last > 0):
			draw.line((start_x - 1, draw_y_last, start_x - 1, draw_y_last - 2), fill=self.COLOR_RGBGreen_0_153_0, width=1)
		elif (bLessThanLimitLast):
			draw.line((start_x - 1, draw_y_last, start_x - 1, draw_y_last - 2), fill=self.COLOR_RGBRed_153_0_0, width=1)
		
		
	def draw_header(self, draw, position_x_first_header, position_x_second_header, position_x_third_header, position_x_fourth_header, average_coverage,\
				ratio_more_zero, rati_more_nine, sample_name, sequence_name, max_coverage, min_coverage, b_only_calculate_size):
		"""
		draw header
		"""
		lines_size = 18
		font_12 = ImageFont.truetype(DrawCoverage.PATH_FONT_BOLD, 12)
		font_16 = ImageFont.truetype(DrawCoverage.PATH_FONT_BOLD, 16)
		
		## draw coverage text
		if (not b_only_calculate_size):
			start_coverage = (font_12.getsize("Coverage")[0] >> 1)
			start_coverage = 5 if (self.get_start_x() - start_coverage) < 0 else (self.get_start_x() - start_coverage) 
			self.draw_text_header(draw, font_12, "Coverage".format(sample_name, sequence_name), start_coverage, self.DRAW_COVERAGE_Y - 17, DrawCoverage.COLOR_RGBGrey_32_32_32, b_only_calculate_size)
			
		position_x = position_x_first_header
		position_x += self.draw_text_header(draw, font_16, "Sample: {}".format(sample_name), position_x, self.START_DRAW_HEADER + 5, DrawCoverage.COLOR_RGBGrey_32_32_32, b_only_calculate_size)
		position_x += self.draw_text_header(draw, font_16, "Sequence: {}".format(sequence_name), position_x, self.START_DRAW_HEADER + 5, DrawCoverage.COLOR_RGBGrey_32_32_32, b_only_calculate_size)
		position_x_first_header = position_x
		
		point_y = self.DRAW_HEADER_BOTTOM
		position_x = position_x_fourth_header;
		position_x += self.draw_text_header(draw, font_12, "Mean depth of coverage: {}".format(average_coverage), position_x, point_y, DrawCoverage.COLOR_RGBGrey_32_32_32, b_only_calculate_size)
		position_x += self.draw_text_header(draw, font_12, "% of size covered by at least 1-fold: {}%".format(ratio_more_zero), position_x, point_y, DrawCoverage.COLOR_RGBGrey_32_32_32, b_only_calculate_size)
		position_x_fourth_header = position_x
		position_x = position_x_third_header
		point_y += (self.START_DRAW_HEADER << 1) + font_12.getsize("text")[1]
		position_x += self.draw_text_header(draw, font_12, "% of size covered by at least {}-fold: {}%".format(\
				self.limit_defined_by_user, rati_more_nine),\
				position_x, point_y, DrawCoverage.COLOR_RGBGrey_32_32_32, b_only_calculate_size)
		position_x_third_header = position_x
		
		### second line
		position_x = position_x_second_header;
		point_y += (self.START_DRAW_HEADER << 1) + font_12.getsize("text")[1]
#		position_x += self.draw_text_header(draw, font_12, "Max. Coverage: {}".format(max_coverage), position_x, point_y, DrawCoverage.COLOR_RGBGrey_32_32_32, b_only_calculate_size)
#		position_x += self.draw_text_header(draw, font_12, "Min. Coverage: {}".format(min_coverage), position_x, point_y, DrawCoverage.COLOR_RGBGrey_32_32_32, b_only_calculate_size)
		if (not b_only_calculate_size):
			draw.line((position_x, point_y + (font_12.getsize("text")[1] >> 1), position_x + lines_size, point_y + (font_12.getsize("text")[1] >> 1)), fill=self.COLOR_RGBRed_153_0_0, width=3)
		position_x += 7
		position_x += self.draw_text_header(draw, font_12, "Cov. <{}".format(self.limit_defined_by_user), position_x  + lines_size, point_y, self.COLOR_RGBRed_153_0_0, b_only_calculate_size)
		position_x += lines_size
		if (not b_only_calculate_size):
			draw.line((position_x, point_y + (font_12.getsize("text")[1] >> 1), position_x + lines_size, point_y + (font_12.getsize("text")[1] >> 1)), fill=self.COLOR_RGBGreen_0_153_0, width=3)
		position_x += 7
		position_x += self.draw_text_header(draw, font_12, "Cov. >{}".format((self.rateDrawCoverage * self.SIZE_COVERAGE_Y)), position_x  + lines_size, point_y, self.COLOR_RGBGreen_0_153_0, b_only_calculate_size)
#		position_x_second_header = position_x
		
		### third line line
#		position_x = position_x_third_header;
#		point_y += (self.START_DRAW_HEADER << 1) + font_12.getsize("text")[1]
		position_x += 8
		if (not b_only_calculate_size):
			draw.line((position_x, point_y + (font_12.getsize("text")[1] >> 1), position_x + lines_size, point_y + (font_12.getsize("text")[1] >> 1)), fill=self.COLOR_RGBRed_153_0_0, width=3)
			draw.ellipse((position_x + lines_size - 4, point_y + (font_12.getsize("text")[1] >> 1) - self.GAP_MARK_VARIATIONS - 1,\
							position_x + lines_size + 4, point_y + (font_12.getsize("text")[1] >> 1) + self.GAP_MARK_VARIATIONS + 1),
							fill = self.COLOR_RGBRed_153_0_0, outline = self.COLOR_RGBRed_153_0_0)
		position_x += lines_size
		position_x += self.draw_text_header(draw, font_12, "Variants AF <50%", position_x  + 10, point_y, DrawCoverage.COLOR_RGBGrey_32_32_32, b_only_calculate_size)
		
		position_x += 8
		if (not b_only_calculate_size):
			draw.line((position_x, point_y + (font_12.getsize("text")[1] >> 1), position_x + lines_size, point_y + (font_12.getsize("text")[1] >> 1)), fill=self.COLOR_RGBBlack, width=3)
			draw.ellipse((position_x + lines_size - 4, point_y + (font_12.getsize("text")[1] >> 1) - self.GAP_MARK_VARIATIONS - 1 ,\
							position_x + lines_size + 4, point_y + (font_12.getsize("text")[1] >> 1) + self.GAP_MARK_VARIATIONS + 1),
							fill = self.COLOR_RGBBlack, outline = self.COLOR_RGBBlack)
		position_x += lines_size
		position_x += self.draw_text_header(draw, font_12, "Variants 50%< AF <90%", position_x  + 10, point_y, DrawCoverage.COLOR_RGBGrey_32_32_32, b_only_calculate_size)
		position_x_second_header = position_x
#		position_x_third_header = position_x
		return (position_x_first_header, position_x_second_header, position_x_third_header, position_x_fourth_header)
	
	
	def draw_text_header(self, draw, font_, text, pointX, pointY, color, bOnlyCalculatesize):
		gap_between_text = 20
		if (not bOnlyCalculatesize): draw.text((pointX, pointY), text, fill=color, font = font_)
		return font_.getsize(text)[0] + gap_between_text
	
	
	def draw_legend_coverage(self, draw, startDraw, endDraw, length):
		"""
		draw lagend
		start_x is to center the text
		"""
		tick_length = 5
		number_space = 100	## the rateImage division is implicit 
		fontsize = 13
		font_ = ImageFont.truetype(DrawCoverage.PATH_FONT_BOLD, fontsize)
		for i in range(0, int(length / number_space)):
			draw.line((self.get_start_x() + i * number_space, self.DRAW_COVERAGE_Y + self.SIZE_COVERAGE_Y,\
					self.get_start_x() + i * number_space, self.DRAW_COVERAGE_Y + self.SIZE_COVERAGE_Y + tick_length), fill=self.COLOR_RGBGrey_64_64_64, width=1)
			middle_size = font_.getsize("{}".format(i * number_space * self.rateImage))[0] >> 1
			draw.text((self.get_start_x() + i * number_space - middle_size, self.DRAW_COVERAGE_Y + self.SIZE_COVERAGE_Y + tick_length),\
					"{}".format(i * number_space * self.rateImage), fill=DrawCoverage.COLOR_RGBGrey_32_32_32, font = font_)
				
		smallOffset = 4
		self.draw_legend_coverage_and_text(draw, font_, startDraw, endDraw, self.DRAW_COVERAGE_Y + smallOffset, "{}".format(self.SIZE_COVERAGE_Y * self.rateDrawCoverage))
		self.draw_legend_coverage_and_text(draw, font_, startDraw, endDraw, self.DRAW_COVERAGE_Y + (self.SIZE_COVERAGE_Y >> 2) + smallOffset, "{}".format(int(((self.SIZE_COVERAGE_Y * self.rateDrawCoverage) / 4) * 3)))
		self.draw_legend_coverage_and_text(draw, font_, startDraw, endDraw, self.DRAW_COVERAGE_Y + (self.SIZE_COVERAGE_Y >> 1) + smallOffset, "{}".format((self.SIZE_COVERAGE_Y * self.rateDrawCoverage) >> 1))
		self.draw_legend_coverage_and_text(draw, font_, startDraw, endDraw, self.DRAW_COVERAGE_Y + ((self.SIZE_COVERAGE_Y >> 2) * 3) + smallOffset, "{}".format(int((self.SIZE_COVERAGE_Y * self.rateDrawCoverage) / 4)))

		
	def draw_legend_coverage_and_text(self, draw, font_, startDraw, endDraw, pointY, value):
		
		smallOffset = 3
		smallOffset_y = (-1 * (font_.getsize("123")[1] >> 1)) - smallOffset
		nLength_X = endDraw - startDraw
		draw.line((startDraw - smallOffset, pointY, endDraw + smallOffset, pointY), fill=self.COLOR_RGBGrey_64_64_64, width=1)
		step = int(nLength_X / 170)
		if (step > 1):
			slice_ = nLength_X / step;
			for i in range(0, step - 1):
				draw.text((self.get_start_x() + (i + 1) * slice_, pointY + smallOffset_y), value, fill=DrawCoverage.COLOR_RGBGrey_32_32_32, font = font_)
		else:
			draw.text((self.get_start_x() + (nLength_X >> 1), pointY + smallOffset_y), value, fill=DrawCoverage.COLOR_RGBGrey_32_32_32, font = font_)

	def draw_variants(self, draw, var_more_50, var_less_50):
		
		for pos in var_more_50:
			self.draw_variant(draw, int(pos / self.rateImage), False)
		for pos in var_less_50:
			self.draw_variant(draw, int(pos / self.rateImage), True)
		
	def draw_variant(self, draw, position, b_up):
		if (b_up):
			draw.line((self.get_start_x() + position, self.GAP_START_GENES_Y - self.GAP_MARK_VARIATIONS, 
					self.get_start_x() + position, self.GAP_START_GENES_Y + self.DRAW_HEIGHT_GENE_SQUARE_Y - (self.GAP_MARK_VARIATIONS << 1)),\
					fill=self.COLOR_RGBRed_153_0_0, width=1)
			draw.ellipse((self.get_start_x() + position - 3, self.GAP_START_GENES_Y - (self.GAP_MARK_VARIATIONS * 3),\
						self.get_start_x() + position + 3, self.GAP_START_GENES_Y - self.GAP_MARK_VARIATIONS),
						fill = self.COLOR_RGBRed_153_0_0, outline = self.COLOR_RGBRed_153_0_0)
		else:
			draw.line((self.get_start_x() + position, self.GAP_START_GENES_Y + (self.GAP_MARK_VARIATIONS << 1), 
					self.get_start_x() + position, self.GAP_START_GENES_Y + self.DRAW_HEIGHT_GENE_SQUARE_Y + self.GAP_MARK_VARIATIONS),\
					fill=self.COLOR_RGBBlack, width=1)
			draw.ellipse((self.get_start_x() + position - 3, self.GAP_START_GENES_Y + self.DRAW_HEIGHT_GENE_SQUARE_Y + self.GAP_MARK_VARIATIONS,\
						self.get_start_x() + position + 3, self.GAP_START_GENES_Y + self.DRAW_HEIGHT_GENE_SQUARE_Y + (self.GAP_MARK_VARIATIONS * 3)),
						fill = self.COLOR_RGBBlack, outline = self.COLOR_RGBBlack)
		
		
class ColorSquares(object):

	def __init__(self, draw):
		self.draw = draw
		self.set_color((0, 0, 204))

	def set_color(self, color):
		self.m_crColor = color
		self.set_colors();
		
	def set_colors(self):
		self.m_crColorLight = self.lighten_color(self.m_crColor, 51)
		self.m_crColorLighter = self.lighten_color(self.m_crColorLight, 51)
		self.m_crColorLightest = self.lighten_color(self.m_crColorLighter, 51)
		self.m_crColorDark = self.darken_color(self.m_crColor, 51)
		self.m_crColorDarker = self.darken_color(self.m_crColorDark, 51)
		self.m_crDkShadow = (59, 238, 244)
		self.m_crLiteShadow = (160, 160, 160)

		## Get a color halfway between COLOR_3DDKSHADOW and COLOR_3DSHADOW
		byRed3DDkShadow = self.m_crDkShadow[0]
		byRed3DLiteShadow = self.m_crLiteShadow[0]
		byGreen3DDkShadow = self.m_crDkShadow[1]
		byGreen3DLiteShadow = self.m_crLiteShadow[1]
		byBlue3DDkShadow = self.m_crDkShadow[2]
		byBlue3DLiteShadow = self.m_crLiteShadow[2]

		self.m_crShadow = (byRed3DLiteShadow + ((byRed3DDkShadow - byRed3DLiteShadow) >> 1),
							  byGreen3DLiteShadow + ((byGreen3DDkShadow - byGreen3DLiteShadow) >> 1),
							  byBlue3DLiteShadow + ((byBlue3DDkShadow - byBlue3DLiteShadow) >> 1))
	
	
	def lighten_color(self, cr_color, by_increase_val):
		byRed = cr_color[0]
		byGreen = cr_color[1]
		byBlue = cr_color[2]
	
		if ((byRed + by_increase_val) <= 255): byRed = byRed + by_increase_val
		if ((byGreen + by_increase_val)	<= 255): byGreen = byGreen + by_increase_val
		if ((byBlue + by_increase_val) <= 255): byBlue = byBlue + by_increase_val
	
		return (byRed, byGreen, byBlue);

	
	def darken_color(self, crColor, byReduceVal):
		byRed = crColor[0]
		byGreen = crColor[1]
		byBlue = crColor[2]
	
		if (byRed >= byReduceVal): byRed = byRed - byReduceVal
		if (byGreen >= byReduceVal): byGreen = byGreen - byReduceVal
		if (byBlue >= byReduceVal): byBlue = byBlue - byReduceVal
	
		return (byRed, byGreen, byBlue)

	def paint_bar(self, rect_point_top, rect_point_bottom, b_vertical):
		if (b_vertical): self.draw.rectangle((rect_point_top, rect_point_bottom), fill=self.m_crColorLightest)
		else: self.draw.rectangle((rect_point_top[0], rect_point_top[1], rect_point_bottom[0], rect_point_bottom[1] - 1), fill=self.m_crColorLightest)

		if (b_vertical): self.draw_vertical_bar(rect_point_top, rect_point_bottom)
		else: self.draw_horizontal_bar(rect_point_top, rect_point_bottom)

	def draw_horizontal_bar(self, rect_point_top, rect_point_bottom):
		height = rect_point_bottom[1] - rect_point_top[1]
		if (height == 0): return;

		nLeft = rect_point_top[0]
		nTop = rect_point_top[1]
		nRight = rect_point_bottom[0]
		nBottom = rect_point_bottom[1]
		
		if (height > 7):
			self.draw.line((nLeft + 2, nBottom - 4, nRight - 1, nBottom - 4), fill=self.m_crColorLight, width=1)
			self.draw.line((nLeft + 1, nTop + 2, nRight - 1, nTop + 2), fill=self.m_crColorLight, width=1)
			self.draw.point((nLeft + 1, nBottom - 3), fill=self.m_crColorLight)
			self.draw.point((nLeft + 1, nTop + 1), fill=self.m_crColorLight)
	
			self.draw.line((nLeft + 2, nBottom - 5, nRight - 2, nBottom - 5), fill=self.m_crColorLighter, width=1)
			self.draw.line((nRight - 2, nBottom - 5, nRight - 2, nTop + 3), fill=self.m_crColorLighter, width=1)
			self.draw.line((nRight - 2, nTop + 3, nLeft + 1, nTop + 3), fill=self.m_crColorLighter, width=1)
			self.draw.point((nLeft + 1, nBottom - 4), fill=self.m_crColorLighter)
			self.draw.point((nLeft + 1, nTop + 2), fill=self.m_crColorLighter)
	
			self.draw.line((nLeft, nBottom - 1, nLeft, nTop), fill=self.m_crColor, width=1)
			self.draw.line((nLeft, nTop, nLeft + 2, nTop), fill=self.m_crColor, width=1)
			self.draw.line((nLeft + 1, nBottom - 3, nRight - 1, nBottom - 3), fill=self.m_crColor, width=1)
			self.draw.line((nLeft + 1, nTop + 1, nRight, nTop + 1), fill=self.m_crColor)
			self.draw.point((nLeft + 1, nBottom - 2), fill=self.m_crColor)
		
			self.draw.line((nLeft + 1, nBottom - 2, nRight - 1, nBottom - 2), fill=self.m_crColorDark, width=1)
			self.draw.line((nRight - 1, nBottom - 2, nRight - 1, nTop + 1), fill=self.m_crColorDark, width=1)
			self.draw.line((nLeft + 1, nTop, nRight, nTop), fill=self.m_crColorDark, width=1)
			self.draw.point((nLeft + 1, nBottom - 1), fill=self.m_crColorDark)
	
			self.draw.line((nLeft + 1, nBottom - 1, nRight, nBottom - 1), fill=self.m_crColorDarker, width=1)
			self.draw.line((nRight, nBottom - 1, nRight, nTop), fill=self.m_crColorDarker, width=1)
		else:
			self.draw.rectangle((nLeft, nTop, nRight, nBottom), fill=self.m_crColorDark)


	def draw_vertical_bar(self, rect_point_top, rect_point_bottom):
		height = rect_point_bottom[1] - rect_point_top[1]
		if (height == 0): return;
	
		nLeft = rect_point_top[0]
		nTop = rect_point_top[1]
		nRight = rect_point_bottom[0]
		nBottom = rect_point_bottom[1]
	
		if (height > 7):
			self.draw.line((nLeft, nTop + 1, nLeft, nTop), fill=self.m_crColor, width=1)
			self.draw.line((nLeft, nTop, nRight, nTop), fill=self.m_crColor, width=1)
			self.draw.line((nLeft + 1, nBottom - 2, nLeft + 1, nTop + 1), fill=self.m_crColor, width=1)
			self.draw.line((nRight - 2, nBottom - 3, nRight - 2, nTop + 1), fill=self.m_crColor, width=1)
			self.draw.point((nRight - 1, nTop + 1), fill=self.m_crColor)
	
			self.draw.line((nLeft + 2, nBottom - 3, nLeft + 2, nTop + 1), fill=self.m_crColorLight, width=1)
			self.draw.line((nRight - 3, nBottom - 3, nRight - 3, nTop + 1), fill=self.m_crColorLight, width=1)
			self.draw.point((nLeft + 1, nTop + 1), fill=self.m_crColorLight)
			self.draw.point((nRight - 2, nTop + 1), fill=self.m_crColorLight)
			
			self.draw.line((nLeft + 3, nBottom - 3, nLeft + 3, nTop + 1), fill=self.m_crColorLighter, width=1)
			self.draw.line((nRight - 4, nBottom - 3, nRight - 4, nTop + 1), fill=self.m_crColorLighter, width=1)
			self.draw.point((nLeft + 2, nTop + 1), fill=self.m_crColorLighter)
			self.draw.point((nRight - 3, nTop + 1), fill=self.m_crColorLighter)
	
			self.draw.line((nLeft, nBottom - 1, nLeft, nTop + 1), fill=self.m_crColorDark, width=1)
			self.draw.line((nLeft + 2, nBottom - 2, nRight - 1, nBottom - 2), fill=self.m_crColorDark, width=1)
			self.draw.line((nRight - 1, nBottom - 2, nRight - 1, nTop + 1), fill=self.m_crColorDark, width=1)
			self.draw.point((nRight, nTop + 1), fill=self.m_crColorDark)

			self.draw.line((nLeft + 1, nBottom - 1, nRight, nBottom - 1), fill=self.m_crColorDarker, width=1)
			self.draw.line((nRight, nBottom - 1, nRight, nTop + 1), fill=self.m_crColorDarker, width=1)
		else:
			self.draw.rectangle((nLeft, nTop, nRight, nBottom), fill=self.m_crColorDark)


	def paintBarsEx(self, rect_point_top, rect_point_bottom, nPos, rgb1, rgb2, shrinkFirst, shrinkLast, text, bVertical):

		if (bVertical):
			rect_point_top_temp = rect_point_top
			rect_point_bottom_temp = (rect_point_bottom[0], nPos)
		else:
			rect_point_top_temp = (rect_point_top[0], rect_point_top[1] + (shrinkFirst << 1))
			rect_point_bottom_temp = (nPos, rect_point_bottom[1] - (shrinkFirst << 1))
		self.set_color(rgb1)
		self.paint_bar(rect_point_top_temp, rect_point_bottom_temp, bVertical);

		if (bVertical):
			rect_point_top_temp = (rect_point_top[0], rect_point_top[1] + nPos)
			rect_point_bottom_temp = (rect_point_bottom[0], rect_point_bottom[1] - nPos)
		else:
			rect_point_top_temp = (nPos + 1, rect_point_top[1] + (shrinkLast << 1))
			rect_point_bottom_temp = (rect_point_bottom[0], rect_point_bottom[1] - (shrinkLast << 1))
		self.set_color(rgb2)
		self.paint_bar(rect_point_top_temp, rect_point_bottom_temp, bVertical)
	
		if (text is not None and len(text)):
			fontsize = 14
			font_ = ImageFont.truetype(DrawCoverage.PATH_FONT_BOLD, fontsize)
			size_x = font_.getsize(text)[0]
			if (bVertical): self.draw.text((rect_point_top[0] - 2, rect_point_top[1] + 5), text, fill=DrawCoverage.COLOR_RGBGrey_32_32_32, font = font_)
			else: 
				middle = nPos + ((rect_point_bottom[0] - nPos) >> 1) - (size_x >> 1) 
				self.draw.text((middle, rect_point_top[1] + 2), text, fill=DrawCoverage.COLOR_RGBGrey_32_32_32, font = font_)

	def draw_text(self, start_x, end_x, y, text):
		fontsize = 14
		font_ = ImageFont.truetype(DrawCoverage.PATH_FONT, fontsize)
		size_x = font_.getsize(text)[0]
		middle = start_x + ((end_x - start_x) >> 1) - (size_x >> 1)
		self.draw.text((middle, y + 2), text, fill=DrawCoverage.COLOR_RGBGrey_32_32_32, font = font_)
				
