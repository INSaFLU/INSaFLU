'''
Created on 31/05/2022

@author: mmp
'''
from utils.utils import Utils
from constants.constants import Constants

def is_all_check_box_in_session(vect_check_to_test, request):
    """
    test if all check boxes are in session 
    If not remove the ones that are in and create the new ones all False
    """
    utils = Utils()
    dt_data = {}
    
    ## get the dictonary
    for key in request.session.keys():
        if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and utils.is_integer(key.split('_')[2])):
            dt_data[key] = True
    
    b_different = False
    if (len(vect_check_to_test) != len(dt_data)): 
        b_different = True
    
    ## test the vector
    if (not b_different):
        for key in vect_check_to_test:
            if (key not in dt_data):
                b_different = True
                break
        
    if (b_different):
        ## remove all
        for key in dt_data:
            del request.session[key] 
    
        ## create new
        for key in vect_check_to_test:
            request.session[key] = False
        return False
    return True

    
    
def clean_check_box_in_session(request):
    """
    check all check boxes on samples/references to add samples
    """
    utils = Utils()
    ## clean all check unique
    if (Constants.CHECK_BOX_ALL in request.session): del request.session[Constants.CHECK_BOX_ALL]
    vect_keys_to_remove = []
    for key in request.session.keys():
        if (key.startswith(Constants.CHECK_BOX) and len(key.split('_')) == 3 and utils.is_integer(key.split('_')[2])):
            vect_keys_to_remove.append(key)
    for key in vect_keys_to_remove:
        del request.session[key] 