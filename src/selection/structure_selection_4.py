'''
Created on Oct 16, 2013

@author: newhouse
'''
from __future__ import division

from collections import defaultdict
from compiler.ast import Print
import math
import numpy
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
    
    It must take as arguments a 'Structure supercollection' (a StoicDict of structure_collection)
    and a 'replica_stoic' which defines the main stoichiometry of the replica calling this module.
    
    It must return a list of 2 or more Structure objects intended for crossing.
    
    This module has the responsibility of selecting stoichiometries in cross-stoic crossings,
    calculating fitness for those stoichiometries, and selecting a structure from each.
    
    The following is an example of a possible implementation
    '''
#    print 'selection mod'
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

    def get_structures(self):

	control_list = self.ui.get_list('control', 'control_in_filelist')
        self.MAX_CASCADE = len(control_list) -1

	#select from collection -1 (relaxed IP structures)
#	structure_coll_a = self.structure_supercoll.get((self.replica_stoic, INITIAL_POOL_REFID))
#        structure_coll_b = self.structure_supercoll.get((self.replica_stoic, INITIAL_POOL_REFID))

	#select from collection 0 which contains all relaxed IP structures as well as children
		#need to change this if I increase levels of the cascade
	structure_coll_a = self.structure_supercoll.get((self.replica_stoic, INITIAL_POOL_REFID+1))
	structure_coll_b = self.structure_supercoll.get((self.replica_stoic, INITIAL_POOL_REFID+1))
	struct_a, fit_a = self.select_structure(structure_coll_a)
        while True: 
            struct_b, fit_b = self.select_structure(structure_coll_b)
            if not struct_a == struct_b: break  # remove this to allow double-selection
        return [struct_a, struct_b]

    
    def select_structure(self, structure_coll):
        '''
        Will calculate fitness for each structure in a collection and select one structure based on a distribution
        
        An optional feature can select from a lower fitness structure. Yet to be implemented.
        '''
        # take care of single-structure case. return only structure
        if len(structure_coll.structures) == 1: 
            return (structure_coll.get_struct(0), 0)
        else: 
            fitness_dict = self.get_energy_fitness(structure_coll)  # this can be altered to select for different fitness
            return self.select_best_from_fitness(fitness_dict)  # TODO: change the selection method
            
    def get_energy_fitness(self, structure_coll):
        '''
        Takes a structure collection of 2 or more structures and returns a dictionary of fitness based on energy
        '''
        # decides whether to reverse fitness values
        #numpy.random.seed(int(time.time() * 83 + len(structure_coll.structures) * 19))
 #    	print numpy.random.seed(int(time.time() * 83 + len(structure_coll.structures) * 19))

	reverse = numpy.random.random() < self.ui.get_eval('selection', 'fitness_reversal_probability')
        # find min and max energy
        # random energy in collection just a base for comparison
        min_e = max_e = structure_coll.itervalues().next().get_property('energy')
        for index, struct in structure_coll:
            try: e = float(struct.get_property('energy'))
            except: 
                message = 'no energy: structure ' + str(struct.get_struct_id()) + '\n'
                message += str(struct.properties)
                output.error(message, self.replica)
            if max_e < e: max_e = e
            if min_e > e: min_e = e
    
        fitness = {}
        for index, struct in structure_coll:
            try: energy = float(struct.get_property('energy'))
            except: pass
            rho = (max_e - energy) / (max_e - min_e)
            if reverse: rho = 1 - rho
            if self.ui.get('selection', 'fitness_function') == 'standard':
                fitness[struct] = rho
            if self.ui.get('selection', 'fitness_function') == 'exponential':
                fitness[struct] = math.exp(-self.ui.get('selection', 'alpha') * rho)  # TODO: implement fitness
        return fitness
        
    def sorted_fitness(self, fitness_dict):
        '''
        returns fitness as a sorted list of tuples.
        '''
        sorted_fitness = sorted(fitness_dict.iteritems(), key=lambda x:x[1])
        # sorted_fitness = sorted_fitness[::-1] # sort fitness 1 to 0
#	print "sorted fitness complete"
#	print sorted_fitness
        return sorted_fitness
    
    def normalized_fitness(self, sorted_fit):
        '''
        takes a sorted list of fitness tuples and returns the list of tuples normalized
        '''
        # sum all fitness
        f_sum = 0
        # reduce all fitnesses
        for index, fitness in sorted_fit: f_sum += fitness
        tmp = []
        for index, fitness in sorted_fit: tmp.append((index, fitness / f_sum))
        
        total = 0
        normalized_fit = []
        for index, t_fitness in tmp: 
            normalized_fit.append((index, t_fitness + total))
            total = t_fitness + total
        return normalized_fit
    
    def select_best_from_fitness(self, fitness_dict):
        fitness = self.sorted_fitness(fitness_dict)
        fitness = self.normalized_fitness(fitness)
       
	#TODO later why is this always done below to seed the generator? why doesnt it work? 
      #  numpy.random.seed(int(time.time() * 87 + len(fitness) * 11))
        random_num = numpy.random.random()
        # selects first element to be greater than random number
	print "best selected from fitness"
        return list(filter((lambda x : x[1] > random_num), fitness))[0]
        
                      
    
    
    
    
