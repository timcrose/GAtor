'''
@author: farren
'''
from __future__ import division
from copy import deepcopy
from math import cos, sin
import math
import numpy as np
import random
from core import user_input, output
from structures.structure import StoicDict, Structure
#from structures import structure_handling

def main(list_of_structures, replica):
    '''
    Args: list of 2 Structures() to crossover, the replica name running the crossover instance.
    Returns: A single Structure() if crossover is successful or False if crossover fails 
    '''

    num_mols = user_input.get_config().get_eval('unit_cell_settings', 'num_molecules')
#    max_attempts = int(user_input.get_config().get_eval('crossover', 'max_attempts'))
    parent_a = list_of_structures[0]
    parent_b = list_of_structures[1]
 
#   cross_attempts = 0
#   while cross_attempts != max_attempts:
#    output.local_message("Attempting crossover: "+str(cross_attempts+1), replica)
    if num_mols == 2:
        geo_opts = [1, 2, 3, 4]
    elif num_mols == 4 or num_mols == 8:
        geo_opts = [1, 2, 3, 4, 5, 6, 7, 8]
    lat_opts = [1, 2, 3, 4, 5, 6, 7]
    cross_method = [random.choice(geo_opts), random.choice(lat_opts)]
    output.local_message("Crossover type:  " + str(cross_method), replica)
    output_parent_properties(parent_a, parent_b, replica)
    cross_obj = Crossover(cross_method, parent_a, parent_b, num_mols, replica)
    child_struct = cross_obj.cross()
#    cell_check = structure_handling.cell_check(child_struct, replica)
#    if cell_check == False:
#        cross_attempts = cross_attempts + 1
#    if cell_check == True:
#    output.local_message("Crossover:  " + str(cross_attempts+1) +" successful" , replica)
#                cross_attempts = max_attempts
    output.local_message("-- Crossover completed" , replica)
    return child_struct

class Crossover(object):
    ''' Takes 2 parent structures and combines them via different crossover options.'''

    def __init__(self, cross_method, parent_a, parent_b, num_mols, replica):
        '''__init__ will always run when a class is initialized. '''
        self.ui = user_input.get_config()
        self.replica = replica
        self.geometry_a = deepcopy(parent_a.get_geometry())
        self.geometry_b = deepcopy(parent_b.get_geometry())
        self.latA_a = deepcopy(parent_a.get_property('lattice_vector_a'))
        self.latB_a = deepcopy(parent_a.get_property('lattice_vector_b'))
        self.latC_a = deepcopy(parent_a.get_property('lattice_vector_c'))
        self.latA_b = deepcopy(parent_b.get_property('lattice_vector_a'))
        self.latB_b = deepcopy(parent_b.get_property('lattice_vector_b'))
        self.latC_b = deepcopy(parent_b.get_property('lattice_vector_c'))
	self.cross_method = cross_method
        self.num_mols = num_mols

    def cross(self):
        child_geometry = (self.combine_mols(self.cross_method, self.num_mols, 
                                self.geometry_a, self.geometry_b))
	child_lattice = (self.combine_lattice(self.cross_method, self.num_mols,
                                self.latA_a, self.latB_a, self.latC_a,
                                self.latA_b, self.latB_b, self.latC_b))
	child_struct = self.create_child_struct(self.cross_method, child_geometry, child_lattice)
	self.output_child_properties(child_struct)
	return child_struct

    def combine_mols(self, cross_method, num_mols, geo_a, geo_b):
        '''combines molecules of two parents based on crossmethod type'''

        num_atom_per_mol = int(len(geo_a)/num_mols)
        mol_listA = [geo_a[x:x+num_atom_per_mol] for x in range(0, len(geo_a), num_atom_per_mol)]
        mol_listB = [geo_b[x:x+num_atom_per_mol] for x in range(0, len(geo_b), num_atom_per_mol)]

        if cross_method[0] == 1: 
            child_geometry = geo_a
        elif cross_method[0] == 2:
            child_geometry = geo_b
        elif cross_method[0] == 3:
            if num_mols == 2:
                child_geometry = np.concatenate((mol_listA[0], mol_listB[1]))
            elif num_mols == 4:
                child_geometry = (np.concatenate((mol_listA[0], mol_listA[1],
                                                 mol_listB[2], mol_listB[3])))
            elif num_mols == 8:
                child_geometry = (np.concatenate((mol_listA[0], mol_listA[1],
                                                 mol_listA[2], mol_listA[3],
                                                 mol_listB[4], mol_listB[5],
                                                 mol_listB[6], mol_listB[7])))
        elif cross_method[0] == 4:
            if num_mols == 2:
                child_geometry = np.concatenate((mol_listB[0], mol_listA[1]))
            elif num_mols == 4:
                child_geometry = (np.concatenate((mol_listB[0], mol_listB[1],
                                                 mol_listA[2], mol_listA[3])))           
            elif num_mols == 8:
                child_geometry = (np.concatenate((mol_listB[0], mol_listB[1],
                                                 mol_listB[2], mol_listB[3],
                                                 mol_listA[4], mol_listA[5],
                                                 mol_listA[6], mol_listA[7])))
        elif cross_method[0] == 5: # every other
            if num_mols == 4:
                child_geometry = (np.concatenate((mol_listA[0], mol_listB[1],
                                                 mol_listA[2], mol_listB[3])))    
            elif num_mols == 8:
                child_geometry = (np.concatenate((mol_listA[0], mol_listB[1],
                                                 mol_listA[2], mol_listB[3],
                                                 mol_listA[4], mol_listB[5],
                                                 mol_listA[6], mol_listB[7])))
        elif cross_method[0] == 6: # every other reversed 
            if num_mols == 4:
                child_geometry = (np.concatenate((mol_listB[0], mol_listA[1],
                                                 mol_listB[2], mol_listA[3])))
            elif num_mols == 8:
                child_geometry = (np.concatenate((mol_listB[0], mol_listA[1],
                                                 mol_listB[2], mol_listA[3],
                                                 mol_listB[4], mol_listA[5],
                                                 mol_listB[6], mol_listA[7])))
        elif cross_method[0] == 7: # every other two
            if num_mols == 4:
                child_geometry = (np.concatenate((mol_listA[0], mol_listB[1],
                                                 mol_listB[2], mol_listA[3])))
	    if num_mols == 8:
                child_geometry = (np.concatenate((mol_listA[0], mol_listA[1],
                                                 mol_listB[2], mol_listB[3],
                                                 mol_listA[4], mol_listA[5],
                                                 mol_listB[6], mol_listB[7])))
        elif cross_method[0] == 8: # every other two reversed
            if num_mols == 4:
                child_geometry = (np.concatenate((mol_listB[0], mol_listA[1],
                                                 mol_listA[2], mol_listB[3])))
            if num_mols == 8:
                child_geometry = (np.concatenate((mol_listB[0], mol_listB[1],
                                                 mol_listA[2], mol_listA[3],
                                                 mol_listB[4], mol_listB[5],
                                                 mol_listA[6], mol_listA[7])))
	return child_geometry


    def combine_lattice(self, cross_method, num_mols, latA_a, latB_a, latC_a,
                                                latA_b, latB_b, latC_b): 
        if cross_method[1] == 1:
		lattice = [latA_a, latB_a, latC_a]
	elif cross_method[1] == 2:
		lattice = [latA_b, latB_b, latC_b]
	elif cross_method[1] == 3:
		lattice = [latA_a, latB_a, latC_b]
	elif cross_method[1] == 4:
		lattice = [latA_a, latB_b, latC_b]
	elif cross_method[1] == 5:
                lattice = [latA_b, latB_b, latC_a]
        elif cross_method[1] == 6:
                lattice = [latA_b, latB_a, latC_a]
	elif cross_method[1] == 7:
                rand_vec = [random.uniform(0.25,0.75) for i in range(3)]
                self.output("Random Frac to Combine Lattices: " +str(rand_vec))
                latA = rand_vec[0]*np.array(latA_a) + (1 - rand_vec[0])*np.array(latA_b)
                latB = rand_vec[1]*np.array(latB_a) + (1 - rand_vec[1])*np.array(latB_b)
                latC = rand_vec[2]*np.array(latC_a) + (1 - rand_vec[2])*np.array(latC_b)
		lattice = [latA, latB, latC]
	return lattice

    def create_child_struct(self,cross_type, child_geometry, child_lattice):
	child_latA = child_lattice[0]
	child_latB = child_lattice[1]
	child_latC = child_lattice[2]
        temp_vol = np.dot(np.cross(child_latA, child_latB), child_latC)
        alpha = self.angle(child_latB, child_latC)
        beta = self.angle(child_latA, child_latC)
        gamma = self.angle(child_latA, child_latB)
        new_struct = Structure() #Create new Structure() for child
        new_struct.build_geo_whole(child_geometry)
        new_struct.set_property('lattice_vector_a', child_latA)
        new_struct.set_property('lattice_vector_b', child_latB)
        new_struct.set_property('lattice_vector_c', child_latC)
        new_struct.set_property('a', np.linalg.norm(child_latA))
        new_struct.set_property('b', np.linalg.norm(child_latB))
        new_struct.set_property('c', np.linalg.norm(child_latC))
        new_struct.set_property('cell_vol', temp_vol)
        new_struct.set_property('crossover_type', cross_type)
        new_struct.set_property('alpha',alpha)
        new_struct.set_property('beta', beta)
        new_struct.set_property('gamma', gamma)
        return new_struct

    def angle(self,v1,v2):
        numdot = np.dot(v1,v2)
        anglerad = np.arccos(numdot/(self.leng(v1)*self.leng(v2)))
        angledeg = anglerad*180/np.pi
        return (angledeg)

    def leng(self,v):
        length = np.linalg.norm(v)
        return length

    def output_child_properties(self, child):
        child_vecs = child.get_lattice_magnitudes()
	child_angs = child.get_lattice_angles()
        message = ('Child lattice vectors:  ' + str(child_vecs)+
                   '\nChild lattice angles:  ' + str(child_angs))
	self.output(message)

    def output(self, message): output.local_message(message, self.replica)

def output_parent_properties(parent_a, parent_b, replica):
    parent_a_vecs = parent_a.get_lattice_magnitudes()
    parent_b_vecs = parent_b.get_lattice_magnitudes()
    parent_a_angs = parent_a.get_lattice_angles()
    parent_b_angs = parent_b.get_lattice_angles()
    message = ('Parent A lattice vectors:  ' + str(parent_a_vecs)+
               '\nParent A lattice angles:  ' + str(parent_a_angs)+
               '\nParent B lattice vectors:  ' + str(parent_b_vecs)+
               '\nParent B lattice angles:  ' + str(parent_b_angs))
    output.local_message(message, replica)
