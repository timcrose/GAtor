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
    # Check conf file, fill initial pool, then spawn
    # parallel GA replicas
    gator = GAtor()
    gator.initialize_MPI_communicator()
    gator.check_conf_file()
    gator.fill_initial_pool()
    gator.spawn_ga_replicas_mpi4py()
#    return

class GAtor():
    """
    This is the master class of GAtor
    Inputs the path to the configuration file
    """
    def __init__(self):
        sname = "GAtor_master"
        self.ui = user_input.get_config()
        self.src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),          
                                                          os.pardir))          
        sys.path.append(self.src_dir)
        if self.ui.has_option(sname,"working_directory"):
            os.chdir(self.ui.get(sname,"working_directory"))
            self.reload_modules()
        if self.ui.get_boolean(sname,"testing_mode"):
            self.testing_mode()
            return 

    def initialize_MPI_communicator(self):
        #from mpi4py import MPI
        self.comm = MPI.COMM_WORLD                                                      
        size = self.comm.Get_size()
        if self.comm.Get_rank() == 0:
            print ("Size of MPI communicator %s" % (size)) 

    def check_conf_file(self):
        """                                                                        
        Runs ui.conf consistency checks
        (e.g. make sure paths exist, etc)

        Will be replaced by full test suite                               
        """ 
        if self.comm.Get_rank() == 0:
            CC = check_conf.ConfigurationChecker() 
            CC.run_checks()
            print ("Done with checking conf file")
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
            print ("Done with filling IP")
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
            print ("Initial pool already filled sucessfully")

    def spawn_ga_replicas_mpi4py(self):
        """                                                                        
        Runs parallel GA replicas using mpi4py
        
        The core modules such as selection, mutation, and crossover    
        are run using the RunGA class from /src/core/run_GA_mpi4py  

        The RunGA class inputs the random                                      
        """  
        from core import run_GA_dummy

        master_group = self.comm.Get_group()

        # First Group
        new_group = master_group.Range_incl([(0, 1, 1)])
        new_comm = self.comm.Create(new_group)
        try:
            print ("first %s" % (new_comm.Get_size()))
            self.newcomm = new_comm
        except: pass

        # Second Group
        #new_group2 = master_group.Range_incl([(2, 3, 1)])    
        new_group2 = master_group.Range_excl([(0, 1, 1)])                   
        new_comm2 = self.comm.Create(new_group)                                 
        try:                                                                   
            print ("second %s" % (new_comm2.Get_size()))
            self.newcomm = new_comm2                         
        except: pass 


        # List of communicators
        replica_comms = [new_comm, new_comm2]
        for comm in replica_comms:
            try:
                if comm.Get_rank() == 0:
                    self.check_initial_pool_filled()
                    stoic = stoic_model.determine_stoic()
                    replica_name = misc.get_random_index()
                    self.make_replica_conf_file(replica_name)
                else:
                    replica_name = None
                replica_name = comm.bcast(replica_name, root=0) 
            except: pass

        self.comm.Barrier()
        try:
            print (replica_name)
        except: pass
        # NEED TO FIGURE OUT HOW TO RUN BOTH GROUPS OF COMMS AT SAME TIME

#        if self.newcomm.Get_rank() == 0:                                      
#            self.check_initial_pool_filled()                          
#            stoic = stoic_model.determine_stoic()                     
#            replica_name = misc.get_random_index()                    
#            self.make_replica_conf_file(replica_name)                 
#        else:                                                         
#            replica_name = None                                       
#        replica_name = self.newcomm.bcast(replica_name, root=0)               
#        print (replica_name)  

# NEED TO FIGURE OUT HOW TO RUNGA BOTH GROUPS OF COMMS AT SAME TIME ######
        stoic = stoic_model.determine_stoic()
        ga = run_GA_dummy.RunGA(replica_name, stoic, new_comm)
        ga.run()
        output.move_to_shared_output(self.ui.get_replica_name())

                  

    def make_replica_conf_file(self, replica_name):
        self.ui.set("parallel_settings","replica_name",replica_name)
        if not os.path.isdir(fh.conf_tmp_dir):
            os.mkdir(fh.conf_tmp_dir)
        conf_path = os.path.join(fh.conf_tmp_dir, self.ui.get_replica_name()+".conf")
        f = open(conf_path,"w")
        self.ui.write(f)
        f.close()
        sys.argv.append(conf_path)
        self.reload_modules()            
        print ("Replica conf file created") 

    def testing_mode(self):
        """
        Can be used for developing test and debug features
        """
        from utilities import test_and_debug
        test_procedure = self.ui.get("test_and_debug","testing_procedure")
        output.time_log("Testing and debugging mode enabled")
        output.time_log("Testing procedure used: "+test_procedure)
        getattr(test_and_debug,test_procedure)()
        return

    def reload_modules(self):
        '''
        These modules need to be reloaded if the configuration file is changed
        '''
        imp.reload(fh)
        imp.reload(user_input)
        imp.reload(output)
        imp.reload(stoic_model)
        imp.reload(misc)

if __name__ == "__main__":
    main()
