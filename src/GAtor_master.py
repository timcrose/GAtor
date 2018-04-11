"""
Master script of the GAtor genetic algorithm

Runs main procedures designated in conf file

If any part of this module is used for a publication please cite:              
                                                                               
F. Curtis, X. Li, T. Rose, A. Vazquez-Mayagoitia, S. Bhattacharya,             
L. M. Ghiringhelli, and N. Marom "GAtor: A First-Principles Genetic            
Algorithm for Molecular Crystal Structure Prediction",                         
J. Chem. Theory Comput., DOI: 10.1021/acs.jctc.7b01152;                        
arXiv 1802.08602 (2018)                                                        
"""

import os, subprocess, shutil
import sys,socket
import time
import core.file_handler as fh
from core import user_input, output, check_conf
from utilities import parallel_run, stoic_model, misc
import imp

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

def main():
    src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 
                                                          os.pardir))
    sys.path.append(src_dir)
    GAtor()

class GAtor():
    """
    This is the master class of GAtor
    Inputs the path to the configuration file
    """
    def __init__(self):
        sname = "GAtor_master"
        self.ui = user_input.get_config()
        if self.ui.has_option(sname,"working_directory"):
            os.chdir(self.ui.get(sname,"working_directory"))
            reload_modules()
        if self.ui.get_boolean(sname,"testing_mode"):
            self.testing_mode()
            return 
        if self.ui.get_boolean(sname,"check_conf_file") and\
            self.ui.is_master_process():
                self.check_conf_file()
        if self.ui.get_boolean(sname,"fill_initial_pool") and\
            self.ui.is_master_process():
                self.check_conf_file()
                self.fill_initial_pool()
        if self.ui.get_boolean(sname,"run_ga"):
            self.check_initial_pool_filled()
            self.run_ga()

    def check_conf_file(self):
        CC = check_conf.ConfigurationChecker() 
        CC.run_checks()

    def fill_initial_pool(self):
        IP_module = fh.my_import(self.ui.get("modules",
                                 "initial_pool_module"),
                                 package="initial_pool")
        fh.mkdir_p(fh.out_tmp_dir)
        IP_module.main()

    def check_initial_pool_filled(self):
        s = "           "
        IP_dat = os.path.join(fh.tmp_dir,"num_IP_structs.dat")   
        msg = "Initial pool has not been filled properly.\n" + s
        msg += "Make sure [GAtor_master]/fill_initial_pool=TRUE has been run.\n" +s
        msg += "Check [initial_pool]/stored_energy_name is correct if not 'energy'."
        if not os.path.isfile(IP_dat):
             raise Exception(msg)
        if os.path.isfile(IP_dat):
            if os.stat(IP_dat).st_size == 0:
                raise Exception(msg)

    def run_ga(self):
        from core import run_GA
        sname = "parallel_settings"
        time.sleep(5)
        if self.ui.get(sname,"parallelization_method") != "serial":
            #Launch parallelism
            parallel_run.launch_parallel()
            return
        if self.ui.get_replica_name() == "master":
            #Need to assign new name
            self.ui.set(sname,"replica_name",misc.get_random_index())
            conf_path = os.path.join(fh.conf_tmp_dir,
                        self.ui.get_replica_name()+".conf")
            f = open(conf_path,"w")
            self.ui.write(f)
            f.close()
            sys.argv.append(conf_path)
            reload_modules()
        message = "GAtor instance reporting from " + \
                                socket.gethostname()
        if self.ui.has_option(sname,"allocated_nodes"):
            nodes = self.ui.get_eval(sname,"allocated_nodes")
            message += "; controlling node(s): " + \
					", ".join(map(str,nodes))
        if self.ui.has_option(sname,"processes_per_replica"):
            ppr = self.ui.get_eval(sname, "processes_per_replica")
            message += "; controlling %i process(es)" % ppr
        if self.ui.has_option(sname,"runjob_block"):
            message += "; block: " + self.ui.get(sname,"runjob_block")
        if self.ui.has_option(sname,"runjob_corner"):
            message += "; corner: " + self.ui.get(sname,"runjob_corner")

        output.time_log(message)
        stoic = stoic_model.determine_stoic()
        ga = run_GA.RunGA(self.ui.get_replica_name(),stoic)
        ga.run()
        output.move_to_shared_output(self.ui.get_replica_name())
	
    def testing_mode(self):
        from utilities import test_and_debug
        test_procedure = self.ui.get("test_and_debug","testing_procedure")
        if self.ui.is_master_process():
            output.time_log("Testing and debugging mode enabled")
            output.time_log("Testing procedure used: "+test_procedure)
        getattr(test_and_debug,test_procedure)()
        return
def reload_modules():
	'''
	These modules need to be reloaded after the configuration file is changed
	'''
	imp.reload(fh)
	imp.reload(user_input)
	imp.reload(output)
	imp.reload(stoic_model)
	imp.reload(misc)
	imp.reload(kill)

if __name__ == "__main__":
    main()
