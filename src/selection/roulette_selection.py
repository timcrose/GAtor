"""                                                                            
If any part of this module is used for a publication please cite:              
                                                                               
F. Curtis, X. Li, T. Rose, A. Vazquez-Mayagoitia, S. Bhattacharya,             
L. M. Ghiringhelli, and N. Marom "GAtor: A First-Principles Genetic            
Algorithm for Molecular Crystal Structure Prediction",                         
J. Chem. Theory Comput., DOI: 10.1021/acs.jctc.7b01152;                        
arXiv 1802.08602 (2018)                                                        
"""    
from __future__ import division
from collections import defaultdict
import numpy as np
from random import choice
import random
import time
import os
from core import user_input, output
from core.file_handler import INITIAL_POOL_REFID, tmp_dir
from structures import structure_collection
from structures.structure import StoicDict
from structures import structure_handling

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

def main(structure_supercoll, replica_stoic, replica=user_input.get_config().get_replica_name()):
    '''
    For a structure selection module, the method main() must be implemented.
    It must take as arguments a 'Structure supercollection'
    and a 'replica_stoic' which defines the main stoichiometry of the replica calling this module.
    It must return a list of 2 or more Structure objects intended for crossing.
    '''
    struct_sel = StructureSelection(structure_supercoll, replica_stoic, replica)
    structure_collection.update_supercollection(structure_supercoll)
    [struct_a, struct_b] = struct_sel.get_structures()
    output.local_message("---- Structure selection ----", replica)
    output.local_message("-- Parent A's ID:      "+ struct_a.struct_id, replica)
    output.local_message("-- Parent A's energy: "+ str(struct_a.get_property('energy')), replica)
    output.local_message("-- Parent B's ID:      "+ struct_b.struct_id, replica) 
    output.local_message("-- Parent B's energy: "+ str(struct_b.get_property('energy')), replica)
    return [struct_a.struct_id, struct_b.struct_id]

class StructureSelection():
    def __init__(self, structure_supercoll, replica_stoic,replica):
        self.ui = user_input.get_config()
        self.structure_supercoll = structure_supercoll
        self.replica = replica
        self.replica_stoic = replica_stoic
        self.index = 0
        self.percent = self.ui.get_eval('selection','percent_best_structs_to_select')
        self.prop = self.ui.get("run_settings","property_to_optimize")
        self.op_style = self.ui.get("run_settings","optimization_style")
        if self.op_style!="maximize" and self.op_style!="minimize":
            message = "Unknown type of optimization style in run_settings; supporing maximize and minimize"
            raise ValueError(message)

    def output(self, message): output.local_message(message, self.replica)

    def get_structures(self):
        banned = []
        structure_coll = self.structure_supercoll.get((self.replica_stoic, self.index))
        struct_a, fit_a = self.select_structures(structure_coll)
        while True: 
            struct_b, fit_b = self.select_structures(structure_coll)
            if not struct_a == struct_b: break  
        struct_a.struct_id, struct_b.struct_id
        return [struct_a, struct_b]

    def select_structures(self, structure_coll):
        '''
        Will calculate fitness for each structure in a collection 
        and select parents
        '''
        if len(structure_coll.structures) == 1:
            return (structure_coll.get_struct(0), 0)
        else:
            if self.ui.get('selection', 'fitness_function') == 'standard_cluster':
                self.output("-- Using shared fitness scheme")
                fitness_dict = self.get_shared_fitness(structure_coll)
                return self.select_best_from_fitness(fitness_dict)
            else:
                fitness_dict = self.get_fitness(structure_coll)
                self.output("-- Using normal fitness scheme")
                return self.select_best_from_fitness(fitness_dict)

    def get_fitness(self, structure_coll):
        '''
        Take a structure collection of 2 or more structures and 
        returns a fitness dictionary
        '''
        reverse = np.random.random() < self.ui.get_eval('selection', 'fitness_reversal_probability')
        prop_list = np.array([])
        for index, structure in structure_coll:
            try:
                prop = structure.get_property(self.prop)
                if self.op_style=="maximize":
                    prop = -prop
                prop_list = np.append(prop,prop_list)
            except KeyError:
                ID = structure.struct_id
                prop = self.prop
                self.output("Structure %s missing the property: %s" % (ID, prop),self.replica)

        prop_list= np.sort(prop_list.reshape(len(prop_list),1),axis=0)
        min_prop = prop_list[0][0]
        max_prop = prop_list[-1][0] 
        fitness = {}
        for index, struct in structure_coll:
            try: prop = float(struct.get_property(self.prop))
            except: pass
            rho = (max_prop - prop) / (max_prop - min_prop)
            if reverse: rho = 1 - rho
            if self.ui.get('selection', 'fitness_function') == 'standard_energy':
                fitness[struct] = rho
        return fitness

    def get_shared_fitness(self, structure_coll):
        '''
        Take a structure collection of 2 or more structures and returns 
        a dictionary of fitness 
        '''
        # First get max and min value of cluster-fitness
        reverse = np.random.random() < (self.ui.get_eval('selection',
                                      'fitness_reversal_probability'))
        prop_list = np.array([])
        for index, structure in structure_coll:
            try:
                prop = structure.get_property(self.prop)
                if self.op_style=="maximize":
                    prop = -prop
                clust_mem = structure.get_property('cluster_members')
                if clust_mem == None:
                    clust_mem = 1
                prop_clus = prop/clust_mem
                prop_list = np.append(prop_clus,prop_list)
            except KeyError:
                ID = structure.struct_id
                prop = self.prop
                message = "Structure %s missing the property: %s" % (ID, prop)
                output.local_message(message, self.replica)
        # Outputs
        prop_list= np.sort(prop_list.reshape(len(prop_list),1),axis=0)
        min_prop = prop_list[0][0]
        max_prop = prop_list[-1][0]
        # Compute relative fitness for all structs
        fitness = {}
        for index, struct in structure_coll:
            prop = struct.get_property(self.prop)
            if self.op_style=="maximize":
                prop = -prop
            clust_mem = struct.get_property('cluster_members')
            if clust_mem == None:
                clust_mem = 1
            prop_clus = prop/clust_mem
            rho = (max_prop - prop_clus) / (max_prop - min_prop)
            if reverse: rho = 1 - rho
            fitness[struct] = rho
        return fitness

 
    def sorted_fitness(self, fitness_dict):
        '''
        returns fitness as a sorted list of tuples.
        '''
        sorted_fitness = sorted(fitness_dict.iteritems(), key=lambda x:x[1])
        return sorted_fitness

    def normalized_fitness(self, sorted_fit):
        '''
        takes a sorted list of fitness tuples and returns the list of tuples normalized
        '''
        f_sum = 0
        tmp = []
        total = 0
        normalized_fit = []
    
        # sum and reduce all fitnesses
        for index, fitness in sorted_fit: f_sum += fitness
        for index, fitness in sorted_fit: tmp.append((index, fitness / f_sum))
        for index, t_fitness in tmp: 
            normalized_fit.append((index, t_fitness + total))
            total = t_fitness + total
        return normalized_fit
    
    def select_best_from_fitness(self, fitness_dict):
        fitness = self.sorted_fitness(fitness_dict)
        fitness = self.normalized_fitness(fitness)
        dec = float(self.percent*0.01)
        dec_comp = 1-dec
        random_num = np.random.uniform(dec_comp,1.0)
        # selects first element to be greater than random number
        return list(filter((lambda x : x[1] > random_num), fitness))[0]
    
