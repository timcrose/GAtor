"""                                                                            
If any part of this module is used for a publication please cite:              
                                                                               
F. Curtis, X. Li, T. Rose, A. Vazquez-Mayagoitia, S. Bhattacharya,             
L. M. Ghiringhelli, and N. Marom "GAtor: A First-Principles Genetic            
Algorithm for Molecular Crystal Structure Prediction",                         
J. Chem. Theory Comput., DOI: 10.1021/acs.jctc.7b01152;                        
arXiv 1802.08602 (2018)                                                        
"""   
import os
import shutil
import user_input

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

class ConfigurationChecker():
    def __init__(self):
        self.ui = user_input.get_config()
        self.control_dir = self.ui.get("FHI-aims","control_in_directory")
        self.control_files = self.ui.get_list("FHI-aims","control_in_filelist")
        self.aims_x = self.ui.get("FHI-aims","path_to_aims_executable")
        self.initial_pool_dir = self.ui.get("initial_pool","user_structures_dir")
        self.selection_module = self.ui.get("modules","selection_module")
        self.skip_energies = self.ui.get_boolean("run_settings", "skip_energy_evaluations")

    def run_checks(self):
        self.check_init_pool_paths()
        self.check_control_paths()
        self.check_aims_executable()
        self.check_selection_parameters()

    def check_init_pool_paths(self):
        if os.path.isdir(self.initial_pool_dir):
            if os.listdir(self.initial_pool_dir) == []:
                msg = ('Initial pool directory %s is empty.' 
                                  % (self.initial_pool_dir))
                raise Exception(msg)
        if not os.path.isdir(self.initial_pool_dir):
            msg = ('Initial pool directory %s doesnt exist.' 
                                  % (self.initial_pool_dir))
            raise IOError(msg)

    def check_control_paths(self):
        if not os.path.isdir(self.control_dir):
            msg = 'Control directory %s doesnt exist.' % (self.control_dir)
            raise IOError(msg)
        for f in self.control_files:
            path = os.path.join(self.control_dir, f)
            if not os.path.isfile(path):
                msg = 'Control file %s doesnt exist.' % (path)
                raise IOError(msg)

    def check_aims_executable(self):
        if not self.skip_energies:
            if not os.path.isfile(self.aims_x):            
                msg = "Aims executable %s does not exist" % (self.aims_x)
                raise IOError(msg)

    def check_selection_parameters(self):
        if self.selection_module == "tournament_selection":
            if not self.ui.has_option("selection","tournament_size"):
                msg = "Tournament selection called but [selection]/"
                msg +="tournament_size not set in .conf file"
                raise ValueError(msg)
