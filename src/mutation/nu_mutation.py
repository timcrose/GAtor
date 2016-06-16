'''
@author: farren
'''
from __future__ import division

from copy import deepcopy
import numpy as np
import random
import time

from core import user_input,output
from structures.structure import Structure

def main(struct, replica):
    '''
    Args: Single Structure() to mutate, the replica name running the crossover instance.
    Returns: A single Structure() if mutation is successful or False if mutation fails 
    '''
    input_struct = deepcopy(struct)
    num_mols = user_input.get_config().get_eval('unit_cell_settings', 'num_molecules')

    mutate_obj = select_mutator(input_struct, num_mols, replica)
    mutated_struct = mutate_obj.mutate()
    return mutated_struct
   
def select_mutator(input_struct, num_mols, replica):
    '''
    In this mutation implementation, there are several classes, each performing a 
    different mutation. This method is responsible for reading the preferences set
    by the user and selecting which mutation to empoly, or no mutation at all.
    Expects: Structure, number of molecules per cell, replica name
    Returns: Mutation Class
    '''
    mutation_list = ["Trans_mol","Rot_mol","Strain_rand_mols","Strain_rand","Strain_sym_mols","Strain_sym"]
    mutation_list = ["Trans_mol"]

    try:
        mut_choice = np.random.choice(mutation_list)
    except:
        mut_choice = mutation_list[int(np.random.random()*len(mutation_list))]
    message = "Mutation Choice:    " +str(mut_choice)
    output.local_message(message, replica)

    if mut_choice == "Trans_mol":
        mutator = RandomTranslationMutation(input_struct, num_mols, replica)
    elif mut_choice == "Rot_mol":
        mutator = RandomRotationMolMutation(input_struct, num_mols, replica)
    elif mut_choice == "Strain_rand_mols":
        mutator = RandomStrainMutationMoveMols(input_struct, num_mols, replica)
    elif mut_choice == "Strain_rand":
        mutator = RandomStrainMutation(input_struct, num_mols, replica)
    elif mut_choice == "Strain_sym_mols":
        mutator = RandomSymmetryStrainMutationMoveMols(input_struct, num_mols, replica)
    elif mut_choice == "Strain_sym":
        mutator = RandomSymmetryStrainMutation(input_struct, num_mols, replica)
    elif mut_choice == "Comp_cell":
        mutator = CompressCellMutationMoveMols(input_struct, num_mols, replica)
    return mutator

class RandomTranslationMutation(object):
    '''
    This mutation gives a random translation to the COM of one or more of the molecules 
    in the unit cell.
    '''
    def __init__(self, input_struct, num_mols, replica):
        self.st_dev = self.ui.get_eval('mutation', 'stand_dev_trans')
        self.geometry = deepcopy(input_struct.get_geometry())
        self.ui = user_input.get_config()
        self.input_struct = input_struct
        self.num_mols = num_mols
        self.replica = replica

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.random_translation()

    def random_translation(self):
        '''Performs random translation to each molecule's COM and returns a Structure'''

        temp_geo = center_geometry(self.geometry) 
        num_atom_per_mol = int(len(temp_geo)/self.num_mols)
        mol_list = [temp_geo[x:x+num_atom_per_mol] for x in range(0, len(temp_geo), num_atom_per_mol)]

        translated_geometry = self.translate_molecules(self, mol_list, st_dev)
        mutated_struct = self.create_mutated_struct(translated_geometry)
        return mutated_struct

    def translate_molecules(self, mol_list, st_dev):
        ''' Randomly displaces the COM of each molecule within gaussian dist'''
        for mol in mol_list:
            rand_disp = np.random.standard_normal(3) * st_dev            
            self.output(rand_disp)
            for atom in mol:
                for i in range(2):
                    atom[i] = atom[i] - rand_disp[i]
        translated_geometry = np.concatenate(np.array(mol_list))
        return translated_geometry

    def create_mutated_struct(translated_geometry):
        ''' Creates Structure from mutated geometry'''
        struct = Structure()
        struct.build_geo_whole(translated_geometry)
        struct.set_property('lattice_vector_a', self.input_struct.get_property('lattice_vector_a'))
        struct.set_property('lattice_vector_b', self.input_struct.get_property('lattice_vector_b'))
        struct.set_property('lattice_vector_c', self.input_struct.get_property('lattice_vector_c'))
        struct.set_property('a', np.linalg.norm(self.input_struct.get_property('lattice_vector_a')))
        struct.set_property('b', np.linalg.norm(self.input_struct.get_property('lattice_vector_b')))
        struct.set_property('c', np.linalg.norm(self.input_struct.get_property('lattice_vector_c')))
        struct.set_property('cell_vol', self.input_struct.get_property('cell_vol'))
        struct.set_property('crossover_type', self.input_struct.get_property('crossover_type'))
        struct.set_property('alpha',self.input_struct.get_property('alpha'))
        struct.set_property('beta', self.input_struct.get_property('beta'))
        struct.set_property('gamma', self.input_struct.get_property('gamma'))
        struct.set_property('mutation_type', 'Trans_mol')
        return struct

#---- Functions Common to All Classes ----#
def center_geometry(self, geometry):
    ''' Centers the origin in relation to the max and min of each axis '''

    for i in range(2):  # x, y, and z
        coordinate_sum = 0
        counter = 0
        for atom in geometry:
            coordinate_sum += atom[i]
            counter += 1
        average = coordinate_sum / counter
        for atom in geometry:
            atom[i] = atom[i] - average
    return geometry



























