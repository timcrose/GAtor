"""                                                                            
If any part of this module is used for a publication please cite:              
                                                                               
F. Curtis, X. Li, T. Rose, A. Vazquez-Mayagoitia, S. Bhattacharya,             
L. M. Ghiringhelli, and N. Marom "GAtor: A First-Principles Genetic            
Algorithm for Molecular Crystal Structure Prediction",                         
J. Chem. Theory Comput., DOI: 10.1021/acs.jctc.7b01152;                        
arXiv 1802.08602 (2018)                                                        
"""

from core import user_input,output
import socket

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

def test_parallelism():
	'''
	Tests the parallelization method on different platforms
	'''
	sname = "parallel_settings"
	if ui.get_boolean(sname,"im_not_master_process"):
		message = "Replica %s reporting from %s" % (ui.get(sname,"replica_name"),socket.gethostname())
		if ui.has_option(sname,"allocated_nodes"):
			message += "; controlling nodes: "+", ".join(map(str,ui.get_eval(sname,"allocated_nodes")))
		if ui.has_option(sname,"processes_per_replica"):
			message += "; controlling %i processes" % (ui.get_eval(sname,"processes_per_replica"))
		output.time_log(message)
		return
	message = "test_parallelism launched to test parallelization method: "+ui.get(sname,"parallelization_method")
	output.time_log(message)


