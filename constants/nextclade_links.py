'''
Created on 11/12/2022

@author: mmp
'''
from managing_files.models import Reference
from constants.constants import Constants
from django.conf import settings

def get_constext_nextclade(media_url_path, context, current_site, specie_identification):
	""" 
	:param specie_identification """
	
	## sarscov 2
	if (specie_identification == Reference.SPECIES_SARS_COV_2):
		context['nextclade_link_covid'] = "{}{}://{}{}".format(
			Constants.NEXTCLADE_LINK_sars_cov_2,
			settings.WEB_SITE_HTTP_NAME,
			current_site,
			media_url_path)
	elif (specie_identification == Reference.SPECIES_MPXV):
		context['nextclade_link_mpxv_hmpxv_b1'] = "{}{}://{}{}".format(
				Constants.NEXTCLADE_LINK_hMPXV_B1,
				settings.WEB_SITE_HTTP_NAME,
				current_site,
				media_url_path)
		context['nextclade_link_mpxv_hmpxv'] = "{}{}://{}{}".format(
				Constants.NEXTCLADE_LINK_hMPXV,
				settings.WEB_SITE_HTTP_NAME,
				current_site,
				media_url_path)
		context['nextclade_link_mpxv_all_clades'] = "{}{}://{}{}".format(
				Constants.NEXTCLADE_LINK_MPXV_All_clades,
				settings.WEB_SITE_HTTP_NAME,
				current_site,
				media_url_path)
	elif (specie_identification == Reference.SPECIES_INFLUENZA):
		context['nextclade_link_a_h3n2'] = "{}{}://{}{}".format(
				Constants.NEXTCLADE_LINK_A_H3N2,
				settings.WEB_SITE_HTTP_NAME,
				current_site,
				media_url_path)
		context['nextclade_link_a_h1n1'] = "{}{}://{}{}".format(
				Constants.NEXTCLADE_LINK_A_H1N1,
				settings.WEB_SITE_HTTP_NAME,
				current_site,
				media_url_path)
		context['nextclade_link_b_yamagata'] = "{}{}://{}{}".format(
				Constants.NEXTCLADE_LINK_B_Yamagata,
				settings.WEB_SITE_HTTP_NAME,
				current_site,
				media_url_path)
		context['nextclade_link_b_victoria'] = "{}{}://{}{}".format(
				Constants.NEXTCLADE_LINK_B_Victoria,
				settings.WEB_SITE_HTTP_NAME,
				current_site,
				media_url_path)
	return context

