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
from core import user_input, output
from core.file_handler import INITIAL_POOL_REFID
from structures import structure_collection
from structures.structure import StoicDict


def main(structure_supercoll, replica_stoic, replica):
    '''
    For a structure selection module, the method main() must be implemented.
    It must take as arguments a 'Structure supercollection'
    and a 'replica_stoic' which defines the main stoichiometry of the replica calling this module.
    It must return a list of 2 or more Structure objects intended for crossing.
    '''
    struct_sel = StructureSelection(replica, structure_supercoll, replica_stoic)
    structure_collection.update_supercollection(structure_supercoll)
    [struct_a, struct_b] = struct_sel.get_structures()  
    return [struct_a, struct_b]

class StructureSelection():
    def __init__(self, replica, structure_supercoll, replica_stoic):
        self.ui = user_input.get_config()
        self.replica = replica
	self.structure_supercoll = structure_supercoll
	self.replica_stoic = replica_stoic
	self.control_list = self.ui.get_list('control', 'control_in_filelist')
        self.max_cascade = len(self.control_list) - 1
	self.index = 0
	self.percent = self.ui.get_eval('selection','percent_best_structs_to_select')

    def get_structures(self):
	control_list = self.ui.get_list('control', 'control_in_filelist')
	structure_coll_a = self.structure_supercoll.get((self.replica_stoic, self.index))
	structure_coll_b = self.structure_supercoll.get((self.replica_stoic, self.index))
	struct_a, fit_a = self.select_structure(structure_coll_a)
        while True: 
            struct_b, fit_b = self.select_structure(structure_coll_b)
            if not struct_a == struct_b: break  # remove this to allow double-selection
        return [struct_a, struct_b]

    def select_structure(self, structure_coll):
        '''
        Will calculate fitness for each structure in a collection and select one structure based on a distribution
        '''
        # take care of single-structure case. return only structure
        if len(structure_coll.structures) == 1: 
            return (structure_coll.get_struct(0), 0)
        else: 
            fitness_dict = self.get_energy_fitness(structure_coll)  # this can be altered to select for different fitness
            return self.select_best_from_fitness(fitness_dict) 
            
    def get_energy_fitness(self, structure_coll):
        '''
        Takes a structure collection of 2 or more structures and returns a dictionary of fitness based on energy
        '''
	reverse = np.random.random() < self.ui.get_eval('selection', 'fitness_reversal_probability')
	e_list = np.array([])
        for index, structure in structure_coll:
		try:
        		energy = structure.get_property('energy')
                	e_list = np.append(energy,e_list)
		except:

                        energy = structure.get_property('energy_tier_1')
                        e_list = np.append(energy,e_list)
#		except ValueError:
#			output.local_message("Structure has no 'energy' property",self.replica)
        e_list= np.sort(e_list.reshape(len(e_list),1),axis=0)
	#output.local_message("e_list" +str(e_list), self.replica)
        min_e = e_list[0][0]
  	max_e = e_list[-1][0] 
#	output.local_message("min e" +str(min_e),self.replica)
#	output.local_message("max " +str(max_e), self.replica) 
        fitness = {}
        for index, struct in structure_coll:
            try: energy = float(struct.get_property('energy'))
            except: pass
            rho = (max_e - energy) / (max_e - min_e)
            if reverse: rho = 1 - rho
            if self.ui.get('selection', 'fitness_function') == 'standard':
                fitness[struct] = rho
            if self.ui.get('selection', 'fitness_function') == 'exponential':
                fitness[struct] = math.exp(-self.ui.get('selection', 'alpha') * rho)  
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
        #random_num = np.random.random()
	random_num = np.random.uniform(dec_comp,1.0)
	output.local_message("selection random num: "+str(random_num),self.replica)
        # selects first element to be greater than random number
        return list(filter((lambda x : x[1] > random_num), fitness))[0]
    
