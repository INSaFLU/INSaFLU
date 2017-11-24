'''
Created on Nov 21, 2017

@author: mmp
'''
from PIL import Image, ImageDraw, ImageFont
import numpy

class Coverage(object):
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
	GAP_START_Y = 40
	START_DRAW_HEADER = 5
	DRAW_HEADER_Y = GAP_START_Y
	DRAW_COVERAGE_Y = DRAW_HEADER_Y + GAP_START_Y
	SIZE_COVERAGE_Y = 200
	SHRINK_GENE_BAR_UTR = 2
	SIZE_TOTAL_COVERAGE_Y = SIZE_COVERAGE_Y + GAP_START_Y	## plus the ruler...
	GAP_START_GENES_Y = SIZE_TOTAL_COVERAGE_Y + DRAW_COVERAGE_Y - 10	## position to draw the exons, green part
	GAP_START_GENES_X = GAP_START_X + 20
	DRAW_HEIGHT_GENE_SQUARE_Y = 20
	DRAW_HEIGHT_GENE_STRAND_Y = 17
	DRAW_STRAND_Y = GAP_START_GENES_Y + DRAW_HEIGHT_GENE_SQUARE_Y + 2
	DRAW_LEGENDE_Y = GAP_START_GENES_Y + DRAW_HEIGHT_GENE_SQUARE_Y + DRAW_HEIGHT_GENE_STRAND_Y + (GAP_START_Y >> 1)
	DRAW_BOTTOM_Y = DRAW_LEGENDE_Y
	LIMIT_MINUNUM_COVERAGE = 30
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
	
	def __init__(self):
		'''
		Constructor
		'''
		
	def create_coverage(self, vect_coverage, vect_genes, var_more_50, var_less_50, output_image,
					average_coverage, ratio_more_zero, ratio_more_nine, sample_name, sequence_name):
		"""
		vect_genes = [[pos_start, pos_end, name, strand 1|-1], [...], ...]
		var_more_50 = [pos, pos1, pos2, ...]
		var_less_50 = [pos, pos1, pos2, ...]
		"""
		## size of the image
		self.maxSizeImage = ((self.GAP_START_GENES_X << 1) + int((len(vect_coverage) / self.rateImage)), self.DRAW_BOTTOM_Y)
		im = Image.new("RGB", (Coverage.DRAW_MINIMUN_LENGTH_X if self.maxSizeImage[0] < Coverage.DRAW_MINIMUN_LENGTH_X\
					else self.maxSizeImage[0], self.maxSizeImage[1]), self.WITHE_COLOR)
		
		draw = ImageDraw.Draw(im) 
		
		self.draw_genes_and_coverage(draw, vect_coverage, vect_genes)
		self.draw_legend_coverage(draw, self.get_start_x(), self.get_start_x() + int(len(vect_coverage) / self.rateImage), int(len(vect_coverage) / self.rateImage))
		self.draw_variants(draw, var_more_50, var_less_50)
		(max_coverage, min_coverage) = self.get_max_min(vect_coverage)
		(position_first_header, position_second_header, position_third_header) = self.draw_header(draw, 0, 0, 0, average_coverage, ratio_more_zero, ratio_more_nine,\
					sample_name, sequence_name, max_coverage, min_coverage, True)		### get only the size
		self.draw_header(draw, (self.maxSizeImage[0] - position_first_header) >> 1, (self.maxSizeImage[0] - position_second_header) >> 1,\
				(self.maxSizeImage[0] - position_third_header) >> 1, average_coverage, ratio_more_zero, ratio_more_nine, sample_name,\
				sequence_name, max_coverage, min_coverage, False)
		im.save(output_image)
		
	def get_start_x(self): return Coverage.GAP_START_GENES_X;
	
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
		if (len(vect_genes) == 0):	## all without genes
			color_squares.set_color(self.COLOR_RGBBlue_0_0_153)
			color_squares.paint_bar((self.get_start_x(), self.GAP_START_GENES_Y), (int(len(vect_coverage) / self.rateImage), self.DRAW_HEIGHT_GENE_SQUARE_Y), False)
		else:
			stop_position = 0
			for vect_info_genes in vect_genes:
				color_squares.paintBarsEx((stop_position + self.get_start_x(), self.GAP_START_GENES_Y),\
								(self.get_start_x() + int(vect_info_genes[1] / self.rateImage),\
								self.GAP_START_GENES_Y + self.DRAW_HEIGHT_GENE_SQUARE_Y),
								self.get_start_x() + int(vect_info_genes[0] / self.rateImage), self.COLOR_RGBBlue_0_0_153, self.COLOR_RGBGreen_0_153_0,\
								self.SHRINK_GENE_BAR_UTR, 0, vect_info_genes[2], False)
				stop_position = int(vect_info_genes[1] / self.rateImage) + 1
				self.draw_strand(draw, vect_info_genes[3], self.get_start_x() + int(vect_info_genes[0] / self.rateImage),\
							int((vect_info_genes[1] - vect_info_genes[0]) / self.rateImage))
			
			if ( int(len(vect_coverage) / self.rateImage) - stop_position > 10):
				color_squares.set_color(self.COLOR_RGBBlue_0_0_153)
				color_squares.paint_bar((stop_position + self.get_start_x(), self.GAP_START_GENES_Y + (self.SHRINK_GENE_BAR_UTR << 1)),\
									(self.get_start_x() + int(len(vect_coverage) / self.rateImage),\
									(self.GAP_START_GENES_Y + self.DRAW_HEIGHT_GENE_SQUARE_Y) - (self.SHRINK_GENE_BAR_UTR << 1)), False)
				
	def draw_strand(self, draw, direction, start_x, nLength_X):
		"""
		"""
		draw.line((start_x, self.DRAW_STRAND_Y + (self.DRAW_HEIGHT_GENE_STRAND_Y >> 1), 
				start_x + nLength_X, self.DRAW_STRAND_Y + (self.DRAW_HEIGHT_GENE_STRAND_Y >> 1)), fill=self.COLOR_RGBGreen_0_102_0, width=1)
		if ((nLength_X >> 4) > 1):
			for i in range(0, nLength_X >> 4):
				slice_pos = (i * (nLength_X / (nLength_X >> 4)))
				draw.point((start_x + slice_pos, self.DRAW_STRAND_Y + (self.DRAW_HEIGHT_GENE_STRAND_Y >> 1) - 1), fill=self.COLOR_RGBGreen_0_102_0)
				draw.point((start_x + slice_pos + direction, self.DRAW_STRAND_Y + (self.DRAW_HEIGHT_GENE_STRAND_Y >> 1) - 2), fill=self.COLOR_RGBGreen_0_102_0)
				draw.point((start_x + slice_pos, self.DRAW_STRAND_Y + (self.DRAW_HEIGHT_GENE_STRAND_Y >> 1) + 1), fill=self.COLOR_RGBGreen_0_102_0)
				draw.point((start_x + slice_pos + direction, self.DRAW_STRAND_Y + (self.DRAW_HEIGHT_GENE_STRAND_Y >> 1) + 2), fill=self.COLOR_RGBGreen_0_102_0)
		else: ## draw in the middle
			draw.point((start_x + (nLength_X >> 1), self.DRAW_STRAND_Y + (self.DRAW_HEIGHT_GENE_STRAND_Y >> 1) - 1), fill=self.COLOR_RGBGreen_0_102_0)
			draw.point((start_x + (nLength_X >> 1) + direction, self.DRAW_STRAND_Y + (self.DRAW_HEIGHT_GENE_STRAND_Y >> 1) - 2), fill=self.COLOR_RGBGreen_0_102_0)
			draw.point((start_x + (nLength_X >> 1), self.DRAW_STRAND_Y + (self.DRAW_HEIGHT_GENE_STRAND_Y >> 1) + 1), fill=self.COLOR_RGBGreen_0_102_0)
			draw.point((start_x + (nLength_X >> 1) + direction, self.DRAW_STRAND_Y + (self.DRAW_HEIGHT_GENE_STRAND_Y >> 1) + 2), fill=self.COLOR_RGBGreen_0_102_0)


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
				if (coverage <= self.LIMIT_MINUNUM_COVERAGE): bLessThanLimitLast = True;

			draw.line((start_x, posYBottom, start_x, draw_y), fill=self.COLOR_RGBGrey_225_225_225, width=1)
			if (draw_y_last == 0): draw.point((start_x, draw_y), fill=self.COLOR_RGBGrey_64_64_64)
			else: draw.line((start_x - 1, draw_y_last, start_x, draw_y), fill=self.COLOR_RGBGrey_64_64_64, width=1)
			
			### limits of coverage
			if (bMoreThanLimitLast and draw_y_last > 0):
				draw.line((start_x - 1, draw_y_last, start_x - 1, draw_y_last - 2), fill=self.COLOR_RGBGreen_0_153_0, width=1)
			elif (bLessThanLimitLast and draw_y_last >  0 and draw_y_last >= (self.DRAW_COVERAGE_Y + self.SIZE_COVERAGE_Y - (self.LIMIT_MINUNUM_COVERAGE / self.rateDrawCoverage))):
				draw.line((start_x - 1, draw_y_last, start_x - 1, draw_y_last - 2), fill=self.COLOR_RGBRed_153_0_0, width=1)
			bMoreThanLimitLast = bMoreThanLimit;
			draw_y_last = draw_y;
			start_x += 1;
		
		## exception for last point
		if (bMoreThanLimitLast and draw_y_last > 0):
			draw.line((start_x - 1, draw_y_last, start_x - 1, draw_y_last - 2), fill=self.COLOR_RGBGreen_0_153_0, width=1)
		elif (bLessThanLimitLast):
			draw.line((start_x - 1, draw_y_last, start_x - 1, draw_y_last - 2), fill=self.COLOR_RGBRed_153_0_0, width=1)
		
		
	def draw_header(self, draw, position_x_first_header, position_x_second_header, position_x_third_header, average_coverage,\
				ratio_more_zero, rati_more_nine, sample_name, sequence_name, max_coverage, min_coverage, b_only_calculate_size):
		"""
		draw header
		"""
		lines_size = 14 
		fontsize = 12
		font_ = ImageFont.truetype(Coverage.PATH_FONT_BOLD, fontsize)
		
		## draw coverage text
		if (not b_only_calculate_size):
			self.draw_text_header(draw, font_, "Coverage".format(sample_name, sequence_name), 5, self.DRAW_COVERAGE_Y - 14, b_only_calculate_size)
			
		position_x = position_x_first_header
		position_x += self.draw_text_header(draw, font_, "Sample/Seq.: {}/{}".format(sample_name, sequence_name), position_x, self.START_DRAW_HEADER, b_only_calculate_size)
		position_x += self.draw_text_header(draw, font_, "Coverage: {}".format(average_coverage), position_x, self.START_DRAW_HEADER, b_only_calculate_size)
		position_x += self.draw_text_header(draw, font_, "Ratio >0: {}%".format(ratio_more_zero), position_x, self.START_DRAW_HEADER, b_only_calculate_size)
		position_x += self.draw_text_header(draw, font_, "Ratio >9: {}%".format(rati_more_nine), position_x, self.START_DRAW_HEADER, b_only_calculate_size)
		position_x_first_header = position_x
		
		### second line
		position_x = position_x_second_header;
		point_y = (self.START_DRAW_HEADER << 1) + font_.getsize("text")[1]
		position_x += self.draw_text_header(draw, font_, "Max. Coverage: {}".format(max_coverage), position_x, point_y, b_only_calculate_size)
		position_x += self.draw_text_header(draw, font_, "Min. Coverage: {}".format(min_coverage), position_x, point_y, b_only_calculate_size)
		if (not b_only_calculate_size):
			draw.line((position_x, point_y + (font_.getsize("text")[1] >> 1), position_x + lines_size, point_y + (font_.getsize("text")[1] >> 1)), fill=self.COLOR_RGBRed_153_0_0, width=3)
		position_x += 7
		position_x += self.draw_text_header(draw, font_, "Cov. <{}".format(self.LIMIT_MINUNUM_COVERAGE), position_x  + lines_size, point_y, b_only_calculate_size)
		position_x += lines_size
		if (not b_only_calculate_size):
			draw.line((position_x, point_y + (font_.getsize("text")[1] >> 1), position_x + lines_size, point_y + (font_.getsize("text")[1] >> 1)), fill=self.COLOR_RGBGreen_0_153_0, width=3)
		position_x += 7
		position_x += self.draw_text_header(draw, font_, "Cov. >{}".format((self.rateDrawCoverage * self.SIZE_COVERAGE_Y)), position_x  + lines_size, point_y, b_only_calculate_size)
		position_x_second_header = position_x
		
		### third line line
		position_x = position_x_third_header;
		point_y += (self.START_DRAW_HEADER << 1) + font_.getsize("text")[1]
		if (not b_only_calculate_size):
			draw.line((position_x, point_y + (font_.getsize("text")[1] >> 1), position_x + lines_size, point_y + (font_.getsize("text")[1] >> 1)), fill=self.COLOR_RGBRed_153_0_0, width=1)
			draw.ellipse((position_x + lines_size - 2, point_y + (font_.getsize("text")[1] >> 1) - self.GAP_MARK_VARIATIONS,\
							position_x + lines_size + 2, point_y + (font_.getsize("text")[1] >> 1) + self.GAP_MARK_VARIATIONS),
							fill = self.COLOR_RGBRed_153_0_0, outline = self.COLOR_RGBRed_153_0_0)
		position_x += 10
		position_x += self.draw_text_header(draw, font_, "Variants AF <50%", position_x  + 10, point_y, b_only_calculate_size)
		
		position_x += 10
		if (not b_only_calculate_size):
			draw.line((position_x, point_y + (font_.getsize("text")[1] >> 1), position_x + lines_size, point_y + (font_.getsize("text")[1] >> 1)), fill=self.COLOR_RGBBlack, width=1)
			draw.ellipse((position_x + lines_size - 2, point_y + (font_.getsize("text")[1] >> 1) - self.GAP_MARK_VARIATIONS,\
							position_x + lines_size + 2, point_y + (font_.getsize("text")[1] >> 1) + self.GAP_MARK_VARIATIONS),
							fill = self.COLOR_RGBBlack, outline = self.COLOR_RGBBlack)
		position_x += 10
		position_x += self.draw_text_header(draw, font_, "Variants AF >50%", position_x  + 10, point_y, b_only_calculate_size)
		position_x_third_header = position_x
		return (position_x_first_header, position_x_second_header, position_x_third_header)
	
	
	def draw_text_header(self, draw, font_, text, pointX, pointY, bOnlyCalculatesize):
		gap_between_text = 20
		if (not bOnlyCalculatesize): draw.text((pointX, pointY), text, fill=Coverage.COLOR_RGBGrey_32_32_32, font = font_)
		return font_.getsize(text)[0] + gap_between_text
	
	
	def draw_legend_coverage(self, draw, startDraw, endDraw, length):
		"""
		draw lagend
		start_x is to center the text
		"""
		tick_length = 5
		number_space = 100	## the rateImage division is implicit 
		fontsize = 13
		font_ = ImageFont.truetype(Coverage.PATH_FONT_BOLD, fontsize)
		for i in range(0, int(length / number_space)):
			draw.line((self.get_start_x() + i * number_space, self.DRAW_COVERAGE_Y + self.SIZE_COVERAGE_Y,\
					self.get_start_x() + i * number_space, self.DRAW_COVERAGE_Y + self.SIZE_COVERAGE_Y + tick_length), fill=self.COLOR_RGBGrey_64_64_64, width=1)
			middle_size = font_.getsize("{}".format(i * number_space * self.rateImage))[0] >> 1
			draw.text((self.get_start_x() + i * number_space - middle_size, self.DRAW_COVERAGE_Y + self.SIZE_COVERAGE_Y + tick_length),\
					"{}".format(i * number_space * self.rateImage), fill=Coverage.COLOR_RGBGrey_32_32_32, font = font_)
				
		smallOffset = 4
		self.draw_legend_coverage_and_text(draw, font_, startDraw, endDraw, self.DRAW_COVERAGE_Y + smallOffset, "{}".format(self.SIZE_COVERAGE_Y * self.rateDrawCoverage))
		self.draw_legend_coverage_and_text(draw, font_, startDraw, endDraw, self.DRAW_COVERAGE_Y + (self.SIZE_COVERAGE_Y >> 2) + smallOffset, "{}".format(int(((self.SIZE_COVERAGE_Y * self.rateDrawCoverage) / 4) * 3)))
		self.draw_legend_coverage_and_text(draw, font_, startDraw, endDraw, self.DRAW_COVERAGE_Y + (self.SIZE_COVERAGE_Y >> 1) + smallOffset, "{}".format((self.SIZE_COVERAGE_Y * self.rateDrawCoverage) >> 1))
		self.draw_legend_coverage_and_text(draw, font_, startDraw, endDraw, self.DRAW_COVERAGE_Y + ((self.SIZE_COVERAGE_Y >> 2) * 3) + smallOffset, "{}".format(int((self.SIZE_COVERAGE_Y * self.rateDrawCoverage) / 4)))

		
	def draw_legend_coverage_and_text(self, draw, font_, startDraw, endDraw, pointY, value):
		
		smallOffset = 3
		smallOffset_y = -1 * (font_.getsize("123")[1] >> 1)
		nLength_X = endDraw - startDraw
		draw.line((startDraw - smallOffset, pointY, endDraw + smallOffset, pointY), fill=self.COLOR_RGBGrey_64_64_64, width=1)
		step = (nLength_X >> 9)
		if (step <= 0): step = int(nLength_X / 150)
		slice_ = nLength_X / step;
		if (step > 1):
			for i in range(0, step - 1):
				draw.text((self.get_start_x() + (i + 1) * slice_, pointY + smallOffset_y), value, fill=Coverage.COLOR_RGBGrey_32_32_32, font = font_)
		else:
			draw.text((self.get_start_x() + (nLength_X >> 1), pointY + smallOffset_y), value, fill=Coverage.COLOR_RGBGrey_32_32_32, font = font_)

	def draw_variants(self, draw, var_more_50, var_less_50):
		
		for pos in var_more_50:
			self.draw_variant(draw, int(pos / self.rateImage), True)
		for pos in var_less_50:
			self.draw_variant(draw, int(pos / self.rateImage), False)
		
	def draw_variant(self, draw, position, b_up):
		if (b_up):
			draw.line((self.get_start_x() + position, self.GAP_START_GENES_Y - self.GAP_MARK_VARIATIONS, 
					self.get_start_x() + position, self.GAP_START_GENES_Y + self.DRAW_HEIGHT_GENE_SQUARE_Y - (self.GAP_MARK_VARIATIONS << 1)),\
					fill=self.COLOR_RGBRed_153_0_0, width=1)
			draw.ellipse((self.get_start_x() + position - 2, self.GAP_START_GENES_Y - self.GAP_MARK_VARIATIONS - self.GAP_MARK_VARIATIONS,\
						self.get_start_x() + position + 2, self.GAP_START_GENES_Y - self.GAP_MARK_VARIATIONS),
						fill = self.COLOR_RGBRed_153_0_0, outline = self.COLOR_RGBRed_153_0_0)
		else:
			draw.line((self.get_start_x() + position, self.GAP_START_GENES_Y + (self.GAP_MARK_VARIATIONS << 1), 
					self.get_start_x() + position, self.GAP_START_GENES_Y + self.DRAW_HEIGHT_GENE_SQUARE_Y + self.GAP_MARK_VARIATIONS),\
					fill=self.COLOR_RGBBlack, width=1)
			draw.ellipse((self.get_start_x() + position - 2, self.GAP_START_GENES_Y + self.DRAW_HEIGHT_GENE_SQUARE_Y + self.GAP_MARK_VARIATIONS,\
						self.get_start_x() + position + 2, self.GAP_START_GENES_Y + self.DRAW_HEIGHT_GENE_SQUARE_Y + self.GAP_MARK_VARIATIONS_BIG),
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
			font_ = ImageFont.truetype(Coverage.PATH_FONT_BOLD, fontsize)
			size_x = font_.getsize(text)[0]
			if (bVertical): self.draw.text((rect_point_top[0] - 2, rect_point_top[1] + 5), text, fill=Coverage.COLOR_RGBGrey_32_32_32, font = font_)
			else: 
				middle = nPos + ((rect_point_bottom[0] - nPos) >> 1) - (size_x >> 1) 
				self.draw.text((middle, rect_point_top[1] + 2), text, fill=Coverage.COLOR_RGBGrey_32_32_32, font = font_)

