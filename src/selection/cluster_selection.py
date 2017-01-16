'''
@authors newhouse, farren
'''
from __future__ import division
from collections import defaultdict
from compiler.ast import Print
import math
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
#Returns ID of two parents instead of structure object.  This is to be used to pass to multiprocessing.

def main(structure_supercoll, replica_stoic,replica=user_input.get_config().get_replica_name()):
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

    def get_structures(self):
        structure_coll = self.structure_supercoll.get((self.replica_stoic, self.index))
        struct_a, fit_a = self.select_structure(structure_coll)
        while True: 
            struct_b, fit_b = self.select_structure(structure_coll)
            if not struct_a == struct_b: break  # remove this to allow double-selection
        return [struct_a, struct_b]

    def compute_RDFs(self, structure_coll):

        for index, struct in structure_coll:
            struct = structure_handling.compute_RDF_vector(struct)
            RDF = struct.get_property('RDF_vector')
            print RDF
    def select_structure(self, structure_coll):
        '''
        Will calculate fitness for each structure in a collection and select parents
        '''
        # take care of single-structure case. return only structure
        if len(structure_coll.structures) == 1: 
            return (structure_coll.get_struct(0), 0)
        else: 
            # this can be altered to select for different fitness
            fitness_dict = self.get_fitness(structure_coll)
            return self.select_best_from_fitness(fitness_dict) 

    def get_fitness(self, structure_coll):
        '''
        Take a structure collection of 2 or more structures and returns 
        a dictionary of fitness 
        '''
        # First get max and min value of cluster-fitness
        rand_num = np.random.random()
        reverse = rand_num< self.ui.get_eval('selection', 'fitness_reversal_probability')
        prop_list = np.array([])
        for index, structure in structure_coll:
            try:
                prop = structure.get_property(self.prop)
                if self.op_style=="maximize":
                    prop = -prop

                clust_mem = structure.get_property('cluster_members')
                prop_clus = prop/clust_mem
                prop_list = np.append(prop_clus,prop_list)
            except KeyError:
                ID = structure.struct_id
                prop = self.prop
                output.local_message("Structure %s missing the property: %s" % (ID, prop),self.replica)

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
            prop_clus = prop/clust_mem
            rho = (max_prop - prop_clus) / (max_prop - min_prop)
            if reverse: rho = 1 - rho
            if self.ui.get('selection', 'fitness_function') == 'standard':
                fitness[struct] = rho
            if self.ui.get('selection', 'fitness_function') == 'exponential':
                fitness[struct] = math.exp(-self.ui.get('selection', 'alpha') * rho)

        return fitness

  
    def select_best_from_fitness(self, fitness_dict):
        fitness = sorted(fitness_dict.iteritems(), key=lambda x:x[1])
        print fitness[0][0].struct_id, fitness[0][1]
        print fitness[-1][0].struct_id, fitness[-1][1]
        dec = float(self.percent*0.01)
        dec_comp = 1-dec
        random_num = np.random.uniform(dec_comp,1.0)

        # selects first element to be greater than random number
        return list(filter((lambda x : x[1] > random_num), fitness))[0]
    
