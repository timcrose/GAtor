"""                                                                            
If any part of this module is used for a publication please cite:              
                                                                               
F. Curtis, X. Li, T. Rose, A. Vazquez-Mayagoitia, S. Bhattacharya,             
L. M. Ghiringhelli, and N. Marom "GAtor: A First-Principles Genetic            
Algorithm for Molecular Crystal Structure Prediction",                         
J. Chem. Theory Comput., DOI: 10.1021/acs.jctc.7b01152;                        
arXiv 1802.08602 (2018)                                                        
"""  
import os
import time
from core.file_handler import cwd, output_file, log_file, tmp_dir, out_tmp_dir
from core import user_input
from external_libs.filelock import FileLock

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

ui = user_input.get_config()

def restart_message(message):
    out_file = os.path.join(cwd, 'restart_relaxations.dat')
    data_file = open(out_file, 'a')
    data_file.write(str(message) + '\n')
    data_file.close()

def local_message(message, replica=ui.get_replica_name()):
    out_file = os.path.join(out_tmp_dir, str(replica) + '.out')
    data_file = open(out_file, 'a')
    data_file.write(str(message) + '\n')
    data_file.close()

def time_log(message,replica=ui.get_replica_name(),logfile=log_file):
	message=time.strftime("%Y-%m-%d %H:%M:%S")+' '+replica+" : "+message+"\n"
	with FileLock(logfile,cwd,3600):
		if not os.path.exists(os.path.join(cwd,logfile)):
			f=open(os.path.join(cwd,logfile),"w")
			f.close()
		f=open(os.path.join(cwd,logfile),"a")
		f.write(message)
		f.close()

def error(message, replica=ui.get_replica_name()):
    if replica == None: r = ''
    else: r = str(replica) + ' '
    out_file = os.path.join(cwd, 'error.out')
    data_file = open(out_file, 'a')
    data_file.write(r + str(message) + '\n')
    data_file.close()

def reset_local(replica=ui.get_replica_name()):
    out_file = os.path.join(out_tmp_dir, str(replica) + '.out')
    data_file = open(out_file, 'w')
    data_file.write(str('') + '\n')
    data_file.close()
    
def move_to_shared_output(replica=ui.get_replica_name(),output_file=output_file):
    local_out_file = os.path.join(out_tmp_dir, str(replica) + '.out')
    if not os.path.exists(local_out_file): pass
    else: 
        d_file = open(local_out_file, 'r')
        contents_string = d_file.read()
        d_file.close()
        
        data_file = open(output_file, 'a')
        #data_file.write('Replica: ' + str(replica)+" --------------------------------\n" + str(contents_string) + '\n')
        data_file.write(contents_string)
        data_file.close()
    reset_local(replica)
