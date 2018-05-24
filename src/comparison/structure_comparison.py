"""                                                                            
If any part of this module is used for a publication please cite:              
                                                                               
F. Curtis, X. Li, T. Rose, A. Vazquez-Mayagoitia, S. Bhattacharya,             
L. M. Ghiringhelli, and N. Marom "GAtor: A First-Principles Genetic            
Algorithm for Molecular Crystal Structure Prediction",                         
J. Chem. Theory Comput., DOI: 10.1021/acs.jctc.7b01152;                        
arXiv 1802.08602 (2018)                                                        
"""


import os
import numpy as np
import multiprocessing
from copy import deepcopy
from core import user_input, output
from core.file_handler import cwd, tmp_dir
from time import time
from pymatgen.analysis.structure_matcher import StructureMatcher,ElementComparator
from pymatgen.analysis.structure_matcher import SpeciesComparator,FrameworkComparator
from core import user_input, output                                            
from core.file_handler import cwd, tmp_dir                                     
from time import time                                                          
from structures import structure_collection                                    
from structures.structure_collection import StructureCollection

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

def main(struct, structure_coll, replica, comm, comparison_type="pre_relaxation_comparison"):
    '''
    Inputs a structure and and structure collection for comparison.
    comparison_type 
    Returns: True/False if the structure passes/fails test for uniqueness and for energy
    '''
    if comparison_type == "pre_relaxation_comparison":
        comp = Comparison(struct, structure_coll, replica, comm, comparison_type)
        structs_to_compare = comp.get_all_structures()
        comm.Barrier()
        dup_result = comp.check_if_duplicate(struct, structs_to_compare, comparison_type)
    elif comparison_type == "post_relaxation_comparison":
        comp = Comparison(struct, structure_coll, replica, comm, comparison_type)
        structs_to_compare = comp.get_similar_energy_structures()
        dup_result = comp.check_if_duplicate(struct, structs_to_compare, comparison_type)
    if not dup_result:
        if comm.Get_rank() == 0:
            output.local_message("-- Structure %s compared is unique. " % (struct.struct_id), replica)
    return dup_result # Boolean

class Comparison:

    def __init__(self, struct, structure_coll, replica, comm, comparison_type):
        if struct is False or struct is None: raise Exception
        self.replica = replica
        self.struct = deepcopy(struct)
        self.structure_coll = structure_coll
        self.ui = user_input.get_config()
        self.comparison_type = comparison_type
        self.stoic = struct.get_stoic()
        self.comm = comm
        self.comparison_type = comparison_type
    
    def output(self, message): output.local_message(message, self.replica)

    def get_all_structures(self):
        '''
        returns full list of structures 
        '''
        struct_list = []
        for index, struct in self.structure_coll:
            struct_list.append((index, struct))
        return struct_list
        
    def get_similar_energy_structures(self):
        '''
        reduces the list of structures that are checked for duplicates within
        a certain window of energy defined by the user
        '''
        e_tol = float(self.ui.get_eval(self.comparison_type, 'energy_comp_window'))
        if e_tol == None: raise Exception

        sim_list = []
        en = float(self.struct.get_property(self.ui.get_property_to_optimize()))
        for index, comp_struct in self.structure_coll:
            comp_en = float(comp_struct.get_property(self.ui.get_property_to_optimize()))
            if en - e_tol <= comp_en <= en + e_tol:
                sim_list.append((index, comp_struct)) 
        if self.comm.Get_rank() == 0:
            self.output("-- Number of structures w/in duplicate energy tolerance: %i" % len(sim_list))
        return sim_list


    def check_if_duplicate(self, struct, comp_list, comparison_type):          
        processes = self.comm.Get_size()                                       
        if self.comm.Get_rank() == 0:                                          
            self.output("-- Comparison done with %i parallel processes" % processes)
                                                                              
        i = 0                
        is_dup = False                                                  
        while is_dup is False and i < len(comp_list):                                                            
            compare = comp_list[i:self.comm.Get_size()+i]                
            if len(compare) < self.comm.Get_size():                            
                while len(compare) < self.comm.Get_size():                     
                    compare.append(comp_list[0])                               
                                                                               
            ID, structc = self.comm.scatter(compare, root=0)                   
            fit = self.compute_pymatgen_fit(struct, structc)                   
            fit_list = self.comm.gather(fit, root=0)         
            if self.comm.Get_rank() == 0:                  
                for fit in fit_list:
                    if fit is False:
                        continue
                    else:
                        ID = fit
                        self.output("-- Structure is a duplicate of another in common pool")   
                        self.output("-- Structure ID in Common pool is: %s" % ID)              
                        index = structure_collection.add_structure(struct, struct.get_stoic() , 'duplicates')
                        self.output("-- Duplicate Structure ID in duplicates pool is: %s" % index)
                        pair = ("0/"+ str(ID),"duplicates/"+str(index))                    
                        dup_output = open(os.path.join(tmp_dir, "GA_duplicates.dat"),'a')      
                        dup_output.write('\t'.join(str(s) for s in pair) + '\n') 
                        is_dup = True   
                        break
            self.comm.Barrier()
            is_dup = self.comm.bcast(is_dup, root=0)                                                                                      
            i += self.comm.Get_size()                                                                                                                 
                                                                               
        return is_dup 

    def set_comp_structure_matcher(self, comparison_type):                     
        '''                                                                    
        Args: self                                                             
        Returns: Pymatgen StructureMatcher object                              
        '''                                                                    
        ui = self.ui                                                           
        L_tol = ui.get_eval(self.comparison_type, 'ltol')                           
        S_tol = ui.get_eval(self.comparison_type, 'stol')                           
        Angle_tol = ui.get_eval(self.comparison_type, 'angle_tol')                  
        Scale = ui.get_boolean(self.comparison_type, 'scale_vol')                   
        sm = (StructureMatcher(ltol=L_tol, 
                               stol=S_tol, 
                               angle_tol=Angle_tol, 
                               primitive_cell=True,
                               scale=Scale, 
                               attempt_supercell=False, 
                               comparator=SpeciesComparator()))
        return sm                                                              
                                                                               
    def compute_pymatgen_fit(self, s1, s2):                             
        ui = user_input.get_config()                                               
        L_tol =ui.get_eval(self.comparison_type, 'ltol')                                
        S_tol = ui.get_eval(self.comparison_type, 'stol')                               
        Angle_tol = ui.get_eval(self.comparison_type, 'angle_tol')                      
        Scale = ui.get_boolean(self.comparison_type, 'scale_vol')                       
        sm =  (StructureMatcher(ltol=L_tol,                                        
                    stol=S_tol,                                                    
                    angle_tol=Angle_tol,                                           
                    primitive_cell=True,                                           
                    scale=Scale,                                                   
                    attempt_supercell=False,                                       
                    comparator=SpeciesComparator()))                               
                                                                                   
        sp1 = s1.get_pymatgen_structure()                                          
        sp2 = s2.get_pymatgen_structure()                                          
        fit = sm.fit(sp1, sp2)                                                     
                                                                                   
        if fit:                                                                    
            return s2.struct_id                                                    
        return fit         
