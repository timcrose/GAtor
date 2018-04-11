"""                                                                            
If any part of this module is used for a publication please cite:              
                                                                               
F. Curtis, X. Li, T. Rose, A. Vazquez-Mayagoitia, S. Bhattacharya,             
L. M. Ghiringhelli, and N. Marom "GAtor: A First-Principles Genetic            
Algorithm for Molecular Crystal Structure Prediction",                         
J. Chem. Theory Comput., DOI: 10.1021/acs.jctc.7b01152;                        
arXiv 1802.08602 (2018)                                                        
"""

from configparser import ConfigParser
from io import StringIO
import time
import os

from core.file_handler import read_data, ui_conf
from structures.structure import Structure, StoicDict
from core.file_handler import structure_dir

__author__ = "Farren Curtis, Xiayue Li, and Timothy Rose"
__copyright__ = "Copyright 2018, Carnegie Mellon University and "+\
                "Fritz-Haber-Institut der Max-Planck-Gessellschaft"
__credits__ = ["Farren Curtis", "Xiayue Li", "Timothy Rose",
               "Alvaro Vazquez-Mayagoita", "Saswata Bhattacharya",
               "Luca M. Ghiringhelli", "Noa Marom"]
__license__ = "BSD-3"
__version__ = "1.0"
__maintainer__ = "Timothy Rose"
__email__ = "trose@andrew.cmu.edu"
__url__ = "http://www.noamarom.com"

def determine_stoic():
    list_of_stoics = [name for name in os.listdir(structure_dir) \
                    if os.path.isdir(os.path.join(structure_dir, name))]
    if len(list_of_stoics) == 1:
        return string_to_stoic(list_of_stoics[0])
    else:
        raise Exception 

def string_to_stoic(stoic_string):
    '''
    Takes a string in the form Mg:2_O:5 and returns the StoicDict representation
    '''
    stoic_dict = StoicDict(int)
    for item in stoic_string.split('_'):
        [key, val] = item.split(':')
        stoic_dict[key] = int(val)
    return stoic_dict

def is_stoic_ini(stoic_string):
    ''' checks if stoic configuration is in normal format'''
    stoic_string_stream = StringIO(stoic_string)
    config = ConfigParser()
    try:  config.readfp(stoic_string_stream)
    except: return False  # not an ini file
    if not config.has_section('stoichiometry'): return False
    options = config.items('stoichiometry')
    for item in options: 
        if item[1].isdigit() == False: return False
    return True

def ini_to_stoic(stoic_string):
    '''gets stoic from config in StoicDict format'''
    stoic_string_stream = StringIO(stoic_string)
    config = ConfigParser()
    config.readfp(stoic_string_stream)
    s_dict = dict(config.items('stoichiometry'))
    stoic_dict = StoicDict(int)
    for element, number in list(s_dict.items()):
        stoic_dict[element.capitalize()] = int(number)
    return stoic_dict    


def atom_to_stoic(stoic_string):
    struct = Structure()
    success = struct.build_geo_whole_atom_format(stoic_string)
    return struct.get_stoic()

def smile_to_stoic(stoic_string):
    
    stoic_string_stream = StringIO(stoic_string)
    config = ConfigParser()
    config.readfp(stoic_string_stream)
    smile_to_take = config.get('molecule_data', 'smile')
    
    from pybel import readstring      # @UnresolvedImport
    mol = readstring("smi", smile_to_take)
    mol.addh()
    aims_string = mol.write(format='fhiaims')
    
    return atom_to_stoic(aims_string)
