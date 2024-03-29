"""                                                                            
If any part of this module is used for a publication please cite:              

F. Curtis, X. Li, T. Rose, A. Vazquez-Mayagoitia, S. Bhattacharya,             
L. M. Ghiringhelli, and N. Marom "GAtor: A First-Principles Genetic            
Algorithm for Molecular Crystal Structure Prediction",                         
J. Chem. Theory Comput., DOI: 10.1021/acs.jctc.7b01152;                        
arXiv 1802.08602 (2018)                                                        
"""
import sys
import ast, os
from core.file_handler import default_config, ui_conf
from select import select
from configparser import SafeConfigParser
from io import StringIO

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

DEFAULT_CONFIG_REPLICA = -1

class ListSafeConfigParser(SafeConfigParser):
    '''Inherits SafeConfigParser and provides list parsing with json'''
    def get_list(self, section, option,eval=False):
        '''provides list parsing with json'''
        if not eval:
            return self.get(section, option).split()
        else:
            return [ast.literal_eval(x) for x in self.get(section,option).split()]
    
    def get_eval(self, section, option):
        return ast.literal_eval(self.get(section, option))

    def get_boolean(self,section,option):
        '''
        Check and see if the section and option is specified
        if specified, has to set to "TRUE", else an error will be raised
        '''
        if self.has_option(section,option):
            if self.get(section,option) == "TRUE":
                return True
        return False

    def get_list_of_booleans(self,section,option):
        '''
		Allows only TRUE and FALSE in the list
		'''
        l = self.get_list(section,option)
        result = []
        for k in l:
            if k == "TRUE":
                result.append(True)
            elif k == "FALSE":
                result.append(False)
        return result

    def get_section_as_dict(self,section,eval=False):
        '''
		Return all the option under a section as a dictionary
		'''
        dicti = {}
        for key in self.options(section):
            if eval:
                try:
                    dicti[key] = self.get_eval(section,key)
                except:
                    dicti[key] = self.get(section,key)
            else:
                dicti[key] = self.get(section,key)
        return dicti

    def get_atom_pair_list(self, section, option):
        return self.get(section, option).split()

    def get_bundled_run_info(self):
        sname = "bundled_run_settings"
        run_names = self.get_list(sname,"run_names")
        runs = []
        for sname in run_names:
            if not self.has_section(sname):
                message = "Bundled run missing section for "
                message += "run: " + sname
                raise KeyError(message)
            d = self.get_section_as_dict(sname,eval=True)
            d["run_name"] = sname
            runs.append(d)
        return runs

    def get_replica_name(self):
        return self.get("parallel_settings","replica_name")
    
    def get_property_to_optimize(self):
        return self.get("run_settings","property_to_optimize")

    def get_nmpc(self):
        return self.get_eval("run_settings","num_molecules")

    def verbose(self):
        return self.get_boolean("run_settings","verbose")

    def all_geo(self):
        return self.get_boolean("run_settings","output_all_geometries")

    def ortho(self):
        return self.get_boolean("run_settings","orthogonalize_unit_cell")

    def is_master_process(self):
        return not self.get_boolean("parallel_settings","im_not_master_process")

    def set_false(self,section,option):
        if self.has_option(section,option):
            self.remove_option(section,option)

    def set_true(self,section,option):
        self.set(section,option,"TRUE")

    def set_working_dir(self,wdir):
        self.set("GAtor_master","working_directory",wdir)

    def __deepcopy__(self,memo):
        '''
        Due to the inability to deepcopy a configuration file
        Will generate a temporary config file and read it back in
        '''
        config_string = StringIO()
        self.write(config_string)
        config_string.seek(0)
        copied = ListSafeConfigParser()
        copied.read_file(config_string)
        return copied

def get_config():
	'''
	Reads in default and user defined UI from the filesystem
	'''
	config = ListSafeConfigParser()

	default_config_file = open(default_config, 'r')
	config.read_file(default_config_file)
	default_config_file.close()

	local_config_file = open(ui_conf, 'r')
	config.read_file(local_config_file)
	local_config_file.close()
	return config

def get_config_with_path(path):
	'''
	Reads in default and user defined UI from the filesystem
	'''
	config = ListSafeConfigParser()

	default_config_file = open(default_config, 'r')
	config.read_file(default_config_file)
	default_config_file.close()

	local_config_file = open(path, 'r')
	config.read_file(local_config_file)
	local_config_file.close()
	return config

def keyboard_input(prompt,allowed_answers=None,time_out=86400, attempts=10):
	'''
	Allows interactive user input
	'''
	user_answer = None
	while attempts>0:
		sys.stdout.write(prompt+" ")
		sys.stdout.flush()
		rlist, _, _ = select([sys.stdin], [], [], time_out)
		if rlist:
			user_answer = sys.stdin.readline()
		else:
			return None
		if allowed_answers==None or user_answer[0:-1] in list(allowed_answers):
			return user_answer[0:-1]
		attempts -= 1
		
