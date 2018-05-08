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
import imp
import core.file_handler as fh
from core import user_input, output, check_conf
from utilities import stoic_model, misc
from mpi4py import MPI

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
    """
    Initializes the initial pool
    Checks the user configuration file for errors
    Runs parallel GA replicas using MPI4py
    """
    gator = GAtor()
    gator.initialize_MPI_communicator()
    gator.check_conf_file()
    gator.fill_initial_pool()
    gator.run_ga_replicas_mpi4py()

class GAtor():
    """
    This is the master class of GAtor
    Inputs the path to the user configuration file
    and stores parameters in self.ui
    """
    def __init__(self):
        self.ui = user_input.get_config()
        self.src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),          
                                                          os.pardir))          
        sys.path.append(self.src_dir)
        if self.ui.has_option("GAtor_master","working_directory"):
            os.chdir(self.ui.get("GAtor_master","working_directory"))
            self.reload_modules()
        if self.ui.get_boolean("GAtor_master","testing_mode"):
            self.testing_mode()
            return 

    def initialize_MPI_communicator(self):
        """
        Initializes global MPI communicator 
        """
        self.comm = MPI.COMM_WORLD                                                      
        size = self.comm.Get_size()
        if self.comm.Get_rank() == 0:
            output.time_log("Size of MPI communicator %s" % (size)) 
        self.comm.Barrier()

    def check_conf_file(self):
        """                                                                        
        Runs ui.conf consistency checks
        (e.g. make sure paths exist, etc)

        Will be replaced by full test suite                               
        """ 
        if self.comm.Get_rank() == 0:
            CC = check_conf.ConfigurationChecker() 
            CC.run_checks()
            output.time_log("Done with checking conf file")
        self.comm.Barrier()


    def fill_initial_pool(self):
        """                                                                        
        Fills the user-definied initial pool
        of structures into the GAtor database
        stored in ./structures                                
        """ 
        if self.comm.Get_rank() == 0:
            IP_module = fh.my_import(self.ui.get("modules",
                                 "initial_pool_module"),
                                 package="initial_pool")
            fh.mkdir_p(fh.out_tmp_dir)
            IP_module.main()
        self.comm.Barrier()

    def check_initial_pool_filled(self):
        """
        Checks the user-defined initial pool
        has been properly filled into the GAtor 
        database
        """
        if self.comm.Get_rank() == 0:   
            s = "           "
            IP_dat = os.path.join(fh.tmp_dir,"num_IP_structs.dat")   
            msg = "Initial pool has not been filled properly.\n" + s
            msg += "Check [GAtor_master]/fill_initial_pool=TRUE has been run.\n" +s
            msg += "Check [initial_pool]/stored_energy_name is correct if not 'energy'."
            if not os.path.isfile(IP_dat):
                 raise Exception(msg)
            if os.path.isfile(IP_dat):
                if os.stat(IP_dat).st_size == 0:
                    raise Exception(msg)
 

    def run_ga_replicas_mpi4py(self):                                        
        """                                                                        
        Runs parallel GA replicas using mpi4py                                 
                                                                               
        The core modules such as selection, mutation, and crossover            
        are run using the RunGA class from /src/core/run_GA_mpi4py             
                                                                                                                  
        """                                                                    
        from core import run_GA_mpi

        # Get user-defined number of GA replicas
        num_replicas = self.ui.get_eval("run_settings","number_of_replicas")

        # Assign MPI color for every replica
        color = self.comm.Get_rank() % num_replicas

        # Split MPI communicator for each replica
        replica_comm  = self.comm.Split(color)
        replica_comm.Barrier()

        # Check initial pool, broadcast replica name to all ranks
        if replica_comm.Get_rank() == 0:                                       
            self.check_initial_pool_filled()                               
            replica_name = misc.get_random_index()                         
            self.make_replica_conf_file(replica_name)                      
        else:                                                              
            replica_name = None      
        replica_comm.Barrier()                                      
        replica_name = replica_comm.bcast(replica_name, root=0) 

        # Output parallelization info to GAtor.log
        output.time_log("Replica %s_%i reporting from %s" %(replica_name, 
                                            replica_comm.Get_rank(), 
                                            MPI.Get_processor_name()))                             
        replica_comm.Barrier()          

        # Run Run_GA module for every replica
        stoic = stoic_model.determine_stoic()
        ga = run_GA_mpi.RunGA(replica_name, stoic, replica_comm)         
        replica_comm.Barrier()
        ga.run()         
        replica_comm.Barrier()
                      
    def make_replica_conf_file(self, replica_name):
        """
        Generates replica-specific conf file for every replica
        """
        self.ui.set("parallel_settings","replica_name",replica_name)
        try:
            if not os.path.isdir(fh.conf_tmp_dir):
                os.mkdir(fh.conf_tmp_dir)
        except: pass
        conf_path = os.path.join(fh.conf_tmp_dir, 
                                 self.ui.get_replica_name()+".conf")
        f = open(conf_path,"w")
        self.ui.write(f)
        f.close()
        sys.argv.append(conf_path)
        self.reload_modules()            

    def testing_mode(self):
        """
        TODO: Can be used for developing test and debug features
        """
        from utilities import test_and_debug
        test_procedure = self.ui.get("test_and_debug","testing_procedure")
        output.time_log("Testing and debugging mode enabled")
        output.time_log("Testing procedure used: "+test_procedure)
        getattr(test_and_debug,test_procedure)()
        return

    def reload_modules(self):
        '''
        These modules are reloaded if the configuration file is modified
        '''
        imp.reload(fh)
        imp.reload(user_input)
        imp.reload(output)
        imp.reload(stoic_model)
        imp.reload(misc)

if __name__ == "__main__":
    main()
