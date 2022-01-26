'''
Created on 23/01/2022

@author: mmp
Has several methods to support Django template language 
'''
import os
from django.utils.safestring import mark_safe

def get_link_for_dropdown_item(link_file, download_file_name = None):
    """ returm HTML to dropdown item """
    return mark_safe('<a rel="nofollow" href="' + \
            link_file + '" download="' + \
            (os.path.basename(link_file) if download_file_name is None else download_file_name) + \
            '" class="dropdown-item"> Download - ' +\
            os.path.basename(link_file) + '</a>')
