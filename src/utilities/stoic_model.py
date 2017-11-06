'''
Created on Nov 12, 2013

@author: newhouse
'''
from ConfigParser import ConfigParser
from StringIO import StringIO
import time
import os

from core.file_handler import read_data, ui_conf
from structures.structure import Structure, StoicDict
from core.file_handler import structure_dir

#def determine_stoic(stoic_path=None):
#    if stoic_path == None: stoic_path = ui_conf
#    try: 
#        stoic_string = read_data(stoic_path)
#        try: return ini_to_stoic(stoic_string)
#        except: pass
#        try: return atom_to_stoic(stoic_string)
#        except: pass
#        try: return smile_to_stoic(stoic_string)
#        except:pass
#        # nothing works
#        raise Exception
#    except: print 'could not determine stoichiometry from input file' 
    # if all else fails

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
    for element, number in s_dict.items():
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
