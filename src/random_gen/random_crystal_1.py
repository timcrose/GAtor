'''
Created on Jul 30, 2013

@author: newhouse
'''

from __future__ import division

import time
import math

from core import user_input
import numpy
from structures.structure import Structure
from utilities.element_radii import radii


def main(stoic, seed):
    # decides what type of mutation to execute with a weighted random
    rand_struct = RandomStructure(stoic,seed)
    return rand_struct.create_random_structure()
    
class RandomStructure(object):
    '''
    classdocs
    TODO Refactor me!
    '''
    
    def __init__(self, stoic, seed):
        '''
        Constructor
        '''
        self.ui = user_input.get_config()
       # MAX_PARAM = 1.1  # this seems to work for all dimensions. could be larger #TODO: make ui
       # self.dimension = self.ui.getint('random_gen', 'rand_struct_dimension')
        #self.min_dist = self.ui.getfloat('random_gen', 'minimum_bond_length')
       # self.closeness_ratio = self.ui.get_eval('random_gen', 'min_distance_ratio')
        self.atom_list = self.stoic_to_list(stoic)

#	self.atom_list =  ['C','C','H','H','N','H','H','H','O','O','C','C','H','H','N','H','H','H','O','O']

        self.seed = abs(int(seed.__hash__()))
#        self.max_dist = len(self.atom_list) * self.min_dist * MAX_PARAM / self.dimension
#	if self.dimension == 1:
#		self.max_dist = 1.2 * self.closeness_ratio * self.min_dist * len(self.atom_list) / 2
#	elif self.dimension == 2:
#		self.max_dist = MAX_PARAM * self.closeness_ratio * self.min_dist * (len(self.atom_list) / math.pi) ** (1.0 / self.dimension)
##		self.max_dist = self.closeness_ratio * self.min_dist * len(self.atom_list) / (2 * 2)
#	else:
 #               self.max_dist = MAX_PARAM * self.closeness_ratio * self.min_dist * (len(self.atom_list) * 3 / (4 * math.pi)) ** (1.0 / self.dimension)
#		self.max_dist = self.closeness_ratio * self.min_dist * len(self.atom_list) / (3 * 2)
        # make array of random numbers to select from in the process
       # self.set_random_array()
        
    def create_random_structure(self):
        while True:
   #         # I'm unsure why the structures sometimes are not filled properly
    #        # But this loop should account for any anomalies 
            struct = self.generate_structure()
      #      acceptable = True
       #     if struct == False: continue
        #    for atom in struct.geometry:
         #       if atom['element'] == '':
          #          acceptable = False
           # if acceptable == True:
            return struct   
                
   # def set_random_array(self):
    #    self.array_index = 0
     #   # multiply by some number to avoid clashes
      #  numpy.random.seed(self.seed)
       # self.random_array = numpy.random.random(size=len(self.atom_list) * self.dimension * 3)
        # times 3 to account for unacceptable coordinates
        
#    def rand_coordinate(self, max_dist):
 #       coordinate = self.random_array[self.array_index] * max_dist
#	print coordinate
 #       self.array_index += 1
  #      # positive or negative?
   #     sign = 1 * numpy.random.randint(0, 2)
    #    if sign < 1: sign = -1
     #   return coordinate * sign
        
    def generate_structure(self):
        struct = Structure()
#	if self.dimension == 1:
#	        dimension_modified = self.dimension
#		max_dist_modified = self.max_dist
#	elif self.dimension == 2:
#	        random_number1 = numpy.random.random()
#	        if random_number1 < 0.2:
#	                dimension_modified = 1
#			max_dist_modified  = 1.2 * self.closeness_ratio * self.min_dist * len(self.atom_list) / 2
#	        else:
#	                dimension_modified = 2
#			max_dist_modified  = self.max_dist
#	else:
#		if self.dimension == 3:
#		        random_number2 = numpy.random.random()
#		        if random_number2 < 0.1:
#		                dimension_modified = 1
#				max_dist_modified  = 1.2 * self.closeness_ratio * self.min_dist * len(self.atom_list) / 2
#		        elif random_number2 >= 0.1 and random_number2 < 0.3:
#		                dimension_modified = 2
##				max_dist_modified = self.closeness_ratio * self.min_dist * len(self.atom_list) / (2 * 1.6)
#		                max_dist_modified = 1.1 * self.closeness_ratio * self.min_dist * (len(self.atom_list) / math.pi) ** (1.0 / dimension_modified)
#		        else:
#		                dimension_modified = 3
#				max_dist_modified  = self.max_dist
#		else:
#			dimension_modified = 3
#			max_dist_modified  = self.max_dist
        counter = 0
#        while True:
        for element in self.atom_list:
            # x, y, z, element
            atom = [0, 0, 0]
            condition = True

	for j in range(2):
                        # if structure is in 2 dimensions, geometry will lie on plane
                 #       atom[j] = self.rand_coordinate(max_dist_modified)
        	atom[j] = numpy.random.random_sample()
	struct.build_geo_by_atom(atom[0], atom[1], atom[2], element)

	return struct
    
    def is_too_close(self, x, y, z, struct):
        is_too_close = False
#        for item in struct.get_geometry():
#            if item['element'] == '': continue  # may remove, never called
#            bx = item['x']
#            by = item['y']
#            bz = item['z']
#            dist = numpy.sqrt((x - bx) ** 2 + (y - by) ** 2 + (z - bz) ** 2)
#            if dist < self.min_dist * self.closeness_ratio: is_too_close = True
        return is_too_close

    def stoic_to_list(self, stoic):
        '''
        retrieves user input structure in the form Mg:2 O:5 
        and returns elements in usable list form
        '''
        # form: StoicDict(<type 'int'>, {'Mg': 2, 'O': 5})
#        atom_list = ['C','C','H','H','N','H','H','H','O','O','C','C','H','H','N','H','H','H','O','O']
	atom_list = []
        for element in stoic:
            to_add = [element] * int(stoic.get(element))
            atom_list.extend(to_add)
        atom_list.sort
        return atom_list  # form: list:[Mg, Mg, O, O, O, O, O]

    def acceptable_distance(self, atom, element, struct):
        rad_a = radii.get(element)
        for atom_b in struct.get_geometry():
            rad_b = radii.get(atom_b['element'])
            dist = math.sqrt((atom[0] - atom_b['x']) ** 2 +
                             (atom[1] - atom_b['y']) ** 2 +
                             (atom[2] - atom_b['z']) ** 2)
            threshold = self.closeness_ratio * (rad_a + rad_b)
            if dist < threshold: return False
        return True
