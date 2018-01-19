'''
@author Farren Curtis
'''
from __future__ import division
import numpy as np
import random
import time
import os
from core import user_input, output
from core.file_handler import INITIAL_POOL_REFID, tmp_dir
from structures import structure_collection
from structures.structure import StoicDict
from structures import structure_handling

def main(structure_supercoll, replica_stoic,
         replica=user_input.get_config().get_replica_name()):
    '''
    For a structure selection module, 
    the method main() must be implemented.
    It must return a list of 2 or more Structure ID's 
    '''
    struct_sel = StructureSelection(structure_supercoll, replica_stoic, replica)
    structure_collection.update_supercollection(structure_supercoll)
    [struct_a, struct_b] = struct_sel.get_structures()
    en_a = struct_a.get_property('energy')
    en_b =  struct_b.get_property('energy')
    output.local_message("---- Structure selection ----", replica)
    output.local_message("-- Parent A's ID:      %s" % struct_a.struct_id, replica)
    output.local_message("-- Parent A's energy: %s" % en_a, replica)
    output.local_message("-- Parent B's ID:      %s" % struct_b.struct_id, replica) 
    output.local_message("-- Parent B's energy: %s" % en_b, replica)
    return [struct_a.struct_id, struct_b.struct_id]

class StructureSelection():
    def __init__(self, structure_supercoll, replica_stoic,replica):
        self.index = 0
        self.replica = replica
        self.replica_stoic = replica_stoic
        self.ui = user_input.get_config()
        self.structure_supercoll = structure_supercoll
        self.percent = (self.ui.get_eval('selection',
                        'percent_best_structs_to_select'))
        self.prop = self.ui.get("run_settings","property_to_optimize")
        self.op_style = self.ui.get("run_settings","optimization_style")
        self.tournament_size = self.ui.get("selection", "tournament_size")
        if self.op_style!="maximize" and self.op_style!="minimize":
            message = "Unknown type of optimization style in run_settings; "
            message += "supporing maximize and minimize"
            raise ValueError(message)

    def output(self, message): output.local_message(message, self.replica)

    def get_structures(self):
        structure_coll = self.structure_supercoll.get((self.replica_stoic, self.index))
        a, b = self.select_structures(structure_coll)
        struct_a = a[0]
        struct_b = b[0]
        #print "selected a", struct_a.struct_id, a[1]
        #print "selected b", struct_b.struct_id, b[1]
        return [struct_a, struct_b]

    def select_structures(self, structure_coll):
        '''
        Will calculate fitness for each structure in a collection 
        and select parents
        '''
        if len(structure_coll.structures) == 1: 
            return (structure_coll.get_struct(0), 0)
        else:
            if self.ui.get_boolean("clustering","cluster_pool"):
                #print "shared fitness"
                self.output("-- Using shared fitness scheme")
                fitness_dict = self.get_shared_fitness(structure_coll)
                return self.select_best_from_fitness(fitness_dict)
            else:
                fitness_dict = self.get_fitness(structure_coll)
                #print "normal"
                self.output("-- Using normal fitness scheme")
                return self.select_best_from_fitness(fitness_dict) 

    def get_fitness(self, structure_coll):
        '''
        Take a structure collection of 2 or more structures and 
        returns a dictionary of fitness 
        '''
        reverse = np.random.random() < (self.ui.get_eval('selection', 
                                      'fitness_reversal_probability'))
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
                message = "Structure %s missing the property: %s" % (ID, prop)
                self.output(message)
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
        # sorted_fitness = sorted_fitness[::-1] # sort fitness 1 to 0
        return sorted_fitness

    def normalized_fitness(self, sorted_fit):
        '''
        inputs a sorted list of fitness tuples and 
        returns the list of tuples normalized
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

        size = int(self.tournament_size)
        tournament_participants = random.sample(list(fitness), size)
        tournament_participants.sort(key=lambda x: x[1])
        all_participants = list(fitness)
        all_participants.sort(key=lambda x: x[1])

        winner = tournament_participants[-1]
        runner_up = tournament_participants[-2]
        if self.ui.get_boolean('selection','select_in_cluster'):
            win_group = winner[0].get_property('cluster_label')
            run_group = runner_up[0].get_property('cluster_label')
            if win_group == run_group:
                self.output("-- Winner group: %s" % (win_group))
                self.output("-- Runner up group equal to winner group")
                pass    
            elif win_group != run_group:
                self.output("-- Runner up group NOT equal to winner group")
                self.output("-- Selecting another structure with same group")
                for structfit in all_participants:
                    if structfit[0] != winner[0]:
                        group = structfit[0].get_property('cluster_label')
                        if group == win_group:
                            runner_up = structfit
                            self.output("-- Selected ID: %s" %(runner_up[0].struct_id))
                            break    
        return winner, runner_up
    
