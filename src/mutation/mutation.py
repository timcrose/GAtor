'''
@author: farren
'''
from __future__ import division

from copy import deepcopy
import numpy as np
import random
import time
import itertools

from core import user_input,output
from structures.structure import Structure
from utilities import space_group_utils as sgu

def main(struct, replica):
    '''
    Args: Single Structure() to mutate, the replica name running the 
    crossover instance.
    Returns: A single Structure() if mutation is successful or False if 
    mutation fails 
    '''
    input_struct = deepcopy(struct)
    ui = user_input.get_config()
    num_mols = ui.get_eval('unit_cell_settings', 'num_molecules')
    napm = int(len(input_struct.geometry)/num_mols)
    tapc = napm*num_mols
    if ui.get_boolean("mutation","enable_symmetry"):
	mutated_struct = symmetric_mutation(input_struct, ui, num_mols, napm, tapc, replica)
    else:  
        mutate_obj = select_mutator(input_struct, num_mols, replica, 
                                    reduced_cell=False, pure_strain=False)
        mutated_struct = mutate_obj.mutate()
    return mutated_struct

def symmetric_mutation(input_struct, ui, num_mols, napm, tapc, replica):
    '''
    Uses space group utilities to perform mutation on the input 
    structure's asymmetric unit
    Returns: Mutated Structure
    '''
    #Will reduce the cell first before mutation
    output.local_message("Symmetry is enabled for mutation",replica)
    if random.random() < 0.99: #50 percent chance of symmetric strain
	output.local_message("Performing strain on whole unit cell")
	mutated_obj = select_mutator(input_struct, num_mols, replica, 
                                     reduced_cell=False, pure_strain=True)
	return mutated_obj.mutate()
    else:
        red_mutated_struct = sgu.reduce_by_symmetry(input_struct)
        symmetry_operations = red_mutated_struct.properties["symmetry_operations"]
        if len(red_mutated_struct.geometry) % napm != 0:
            output.local_message("Structure reduction by symmetry failed", replica)
            return False

        num_mols = int(len(red_mutated_struct.geometry)/napm)
        mutate_obj = select_mutator(red_mutated_struct, num_mols, replica, pure_strain=False)
        red_mutated_struct = mutate_obj.mutate()
        red_mutated_struct.properties["symmetry_operations"] = symmetry_operations
        if red_mutated_struct == False:
            return False

        mutated_struct = sgu.rebuild_by_symmetry(red_mutated_struct,napm=napm,create_duplicate=True)
        if mutated_struct == False:
            output.local_message('Structure reconstruction by symmetry failed',replica)
            return False
        if len(mutated_struct.geometry)!=tapc:
            output.local_message('Structure reconstruction by symmetry failed',replica)
            return False
        mutated_struct.set_property('mutation_type', red_mutated_struct.get_property('mutation_type'))
        mutated_struct.set_property('crossover_type', red_mutated_struct.get_property('crossover_type'))
        return mutated_struct

def select_mutator(input_struct, num_mols, replica, reduced_cell=True, pure_strain=True):
    '''
    In this mutation implementation, there are several classes, each performing a 
    different mutation. This method selecting which mutation to empoly
    Expects: Structure, number of molecules per cell, replica name
    Returns: Mutation Class
    '''
    if reduced_cell and num_mols == 1 and not pure_strain:
        mutation_list = ["Trans_mol", "Rot_mol"]
    elif reduced_cell and num_mols > 1 and not pure_strain:
	mutation_list = ["Trans_mol", "Rot_mol", "Sym_rot_mol", "Swap_mol"]
    elif not reduced_cell and pure_strain:
	mutation_list = ["Strain_sym_mols"]
    elif not reduced_cell and not pure_strain:
        mutation_list = (["Trans_mol","Rot_mol", "Sym_rot_mol",
                          "Strain_rand_mols","Strain_sym_mols", "Swap_mol"])
    try:
        mut_choice = np.random.choice(mutation_list)
    except:
        mut_choice = mutation_list[int(np.random.random()*len(mutation_list))]
    message = "-- Mutation Choice: %s" % (mut_choice)
    message += "\n-- Reduced Cell: %s" % (reduced_cell)
    output.local_message(message, replica)

    if mut_choice == "Trans_mol":
        mutator = RandomTranslationMutation(input_struct, num_mols, replica)
    elif mut_choice == "Rot_mol":
        mutator = RandomRotationMolMutation(input_struct, num_mols, replica)
    elif mut_choice == "Sym_rot_mol":
	mutator = SymRotationMolMutation(input_struct, num_mols, replica)
    elif mut_choice == "Swap_mol":
	mutator = SwapMolMutation(input_struct, num_mols, replica)
    elif mut_choice == "Strain_rand_mols":
        mutator = RandomStrainMutationMoveMols(input_struct, num_mols, replica)
    elif mut_choice == "Strain_rand":
        mutator = RandomStrainMutation(input_struct, num_mols, replica)
    elif mut_choice == "Strain_sym_mols":
        mutator = RandomSymmetryStrainMutationMoveMols(input_struct, num_mols, replica)
    elif mut_choice == "Strain_sym":
        mutator = RandomSymmetryStrainMutation(input_struct, num_mols, replica)
    return mutator

class RandomTranslationMutation(object):
    '''
    This mutation gives a random translation to the COM of one or more of the molecules 
    in the unit cell.
    '''
    def __init__(self, input_struct, num_mols, replica):
        self.input_struct = input_struct
        self.num_mols = num_mols
        self.replica = replica
        self.geometry = deepcopy(input_struct.get_geometry())
        self.ui = user_input.get_config()
        self.st_dev = self.ui.get_eval('mutation', 'stand_dev_trans')
        self.A = self.input_struct.get_property('lattice_vector_a')
        self.B = self.input_struct.get_property('lattice_vector_b')      
        self.C = self.input_struct.get_property('lattice_vector_c')
        self.alpha = self.input_struct.get_property('alpha')
        self.beta = self.input_struct.get_property('beta')
        self.gamma = self.input_struct.get_property('gamma')
        self.cell_vol = self.input_struct.get_unit_cell_volume()
        self.cross_type = self.input_struct.get_property('crossover_type')

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.random_translation()

    def random_translation(self):
        '''Calls for a random translation to each molecule's COM and returns a Structure'''

        temp_geo = self.geometry 
        num_atom_per_mol = int(len(temp_geo)/self.num_mols)
        mol_list = [temp_geo[x:x+num_atom_per_mol] for x in range(0, len(temp_geo), num_atom_per_mol)]

        translated_geometry = self.translate_molecules(mol_list, self.st_dev)
        mutated_struct = self.create_mutated_struct(translated_geometry)
        return mutated_struct

    def translate_molecules(self, mol_list, st_dev):
        ''' Randomly displaces the COM of each molecule within gaussian dist'''
        for mol in mol_list:
            rand_disp = np.random.standard_normal(3) * st_dev            
            #self.output(rand_disp)
            for atom in mol:
                for i in range(2):
                    atom[i] = atom[i] - rand_disp[i]
        translated_geometry = np.concatenate(np.array(mol_list))
        return translated_geometry

    def create_mutated_struct(self, translated_geometry):
        ''' Creates Structure from mutated geometry'''
        struct = Structure()
        struct.build_geo_whole(translated_geometry)
        struct.set_property('lattice_vector_a', self.A)
        struct.set_property('lattice_vector_b', self.B)
        struct.set_property('lattice_vector_c', self.C)
        struct.set_property('a', leng(self.A))
        struct.set_property('b', leng(self.B))
        struct.set_property('c', leng(self.C))
        struct.set_property('cell_vol', self.cell_vol)
        struct.set_property('crossover_type', self.cross_type)
        struct.set_property('alpha',self.alpha)
        struct.set_property('beta', self.beta)
        struct.set_property('gamma', self.gamma)
        struct.set_property('mutation_type', 'trans_mol')
        return struct

class RandomRotationMolMutation(object):
    ''' Gives a random rotation to the COM of the molecules in the unit cell.'''

    def __init__(self, input_struct, num_mols, replica):
        self.input_struct = input_struct
        self.num_mols = num_mols
        self.replica = replica
        self.geometry = deepcopy(input_struct.get_geometry())
        self.ui = user_input.get_config()
        self.st_dev = self.ui.get_eval('mutation', 'stand_dev_rot')
        self.A = self.input_struct.get_property('lattice_vector_a')
        self.B = self.input_struct.get_property('lattice_vector_b')      
        self.C = self.input_struct.get_property('lattice_vector_c')
        self.alpha = self.input_struct.get_property('alpha')
        self.beta = self.input_struct.get_property('beta')
        self.gamma = self.input_struct.get_property('gamma')
        self.cell_vol = self.input_struct.get_unit_cell_volume()
        self.cross_type = self.input_struct.get_property('crossover_type')

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.random_rotation()

    def random_rotation(self):
        '''Calls for a random rotation to each molecule's COM and returns a Structure'''
        temp_geo = self.geometry 
        num_atom_per_mol = int(len(temp_geo)/self.num_mols)
        mol_list = [temp_geo[x:x+num_atom_per_mol] for x in range(0, len(temp_geo), num_atom_per_mol)]
        atom_type_list = [temp_geo[i][3] for i in range(len(temp_geo))]
	mol_list_COM = get_COM_mol_list(mol_list)
        rotated_geometry = self.rotate_molecules(mol_list, mol_list_COM, self.st_dev)
        rotated_struct = self.create_rotated_struct(rotated_geometry, atom_type_list)
        return rotated_struct
	
    def rotate_molecules(self, mol_list, mol_list_COM, st_dev):
        ''' Randomly rotates each molecule within gaussian dist'''
        rot_geometry = []
	i = 0 
        for mol in mol_list:
            rand_vec = np.random.normal(scale=st_dev,size=3)
	    #rand_vec = [0, 0, 90]
	    self.output("-- Random Rotation %s" % (rand_vec))
            theta= (np.pi/180)*rand_vec[0]
            psi = (np.pi/180)*rand_vec[1]
            phi= (np.pi/180)*rand_vec[2]
	    mol_COM = np.array([mol_list_COM[i][0], mol_list_COM[i][1],mol_list_COM[i][2]])
            for atom in mol:
                atom_vec = np.array([atom[0], atom[1], atom[2]]) - mol_COM
                rot_geometry.append(np.dot(np.asarray(rotation_matrix(theta, psi, phi)),atom_vec) + mol_COM)
	    i+=1
        return rot_geometry

    def create_rotated_struct(self, rotated_geo, atom_types):
        ''' Creates Structure from mutated geometry'''
        struct = Structure()
        for i in range(len(rotated_geo)):
            struct.build_geo_by_atom(float(rotated_geo[i][0]), float(rotated_geo[i][1]),
                                     float(rotated_geo[i][2]), atom_types[i])
        struct.set_property('lattice_vector_a', self.A)
        struct.set_property('lattice_vector_b', self.B)
        struct.set_property('lattice_vector_c', self.C)
        struct.set_property('a', leng(self.A))
        struct.set_property('b', leng(self.B))
        struct.set_property('c', leng(self.C))
        struct.set_property('cell_vol', self.cell_vol)
        struct.set_property('crossover_type', self.cross_type)
        struct.set_property('alpha',self.alpha)
        struct.set_property('beta', self.beta)
        struct.set_property('gamma', self.gamma)
        struct.set_property('mutation_type', 'rot_mol')
        return struct

class SymRotationMolMutation(object):
    ''' 
    Gives a symmetric rotation to the COM of the molecules in the unit cell.
    For example, for Z=4 it may rotate molecules 1 & 2 the same amount, leaving
    3 & 4 alone
    '''

    def __init__(self, input_struct, num_mols, replica):
        self.input_struct = input_struct
        self.num_mols = num_mols
        self.replica = replica
        self.geometry = deepcopy(input_struct.get_geometry())
        self.ui = user_input.get_config()
        self.st_dev = self.ui.get_eval('mutation', 'stand_dev_rot')
        self.A = self.input_struct.get_property('lattice_vector_a')
        self.B = self.input_struct.get_property('lattice_vector_b')
        self.C = self.input_struct.get_property('lattice_vector_c')
        self.alpha = self.input_struct.get_property('alpha')
        self.beta = self.input_struct.get_property('beta')
        self.gamma = self.input_struct.get_property('gamma')
        self.cell_vol = self.input_struct.get_unit_cell_volume()
        self.cross_type = self.input_struct.get_property('crossover_type')

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.sym_rotation()

    def sym_rotation(self):
        '''Calls for a random rotation to each molecule's COM and returns a Structure'''
        temp_geo = self.geometry
        num_atom_per_mol = int(len(temp_geo)/self.num_mols)
        mol_list = [temp_geo[x:x+num_atom_per_mol] for x in range(0, len(temp_geo), num_atom_per_mol)]
        atom_type_list = [temp_geo[i][3] for i in range(len(temp_geo))]
        mol_list_COM = get_COM_mol_list(mol_list)
        rotated_geometry = self.rotate_molecules_sym(mol_list, mol_list_COM, self.st_dev)
        rotated_struct = self.create_rotated_struct(rotated_geometry, atom_type_list)
        return rotated_struct

    def rotate_molecules_sym(self, mol_list, mol_list_COM, st_dev):
        ''' Randomly rotates each molecule within gaussian dist'''
        rot_geometry = []
	ang = np.random.normal(scale=30,size=1)
	rots = ([ang[0], 0, 0], [0, ang[0], 0], [0, 0, ang[0]],
                     [ang[0], ang[0], 0], [0, ang[0], ang[0]], [ang[0], 0 , ang[0]])
	rot = random.choice(rots)
	self.output("-- Rotation %s" % (rot))
	mol_in = 0
        mol_info = []
        mol_combos = list(itertools.combinations(range(self.num_mols), 2))
        for mol in mol_list:
	    mol_info.append([mol, mol_list_COM[mol_in]])
	    mol_in +=1
	mol_combo = random.choice(mol_combos)

	self.output("-- On molecules: ")
	self.output(mol_combo)
	i = 0
        for mol in mol_list:
	    if mol_combo[0] == i:
                mol_COM = np.array([
                                   mol_list_COM[mol_combo[0]][0], 
                                   mol_list_COM[mol_combo[0]][1],
                                   mol_list_COM[mol_combo[0]][2]])
                for atom in mol:
                    atom_vec = np.array([atom[0], atom[1], atom[2]]) - mol_COM
                    rot_geometry.append(np.dot(np.asarray(rotation_matrix(rot[0], rot[1], rot[2])),atom_vec) + mol_COM)
	    elif mol_combo[1] == i:
                mol_COM = np.array([
                                   mol_list_COM[mol_combo[1]][0], 
                                   mol_list_COM[mol_combo[1]][1],
                                   mol_list_COM[mol_combo[1]][2]])
                for atom in mol:
                    atom_vec = np.array([atom[0], atom[1], atom[2]]) - mol_COM
                    rot_geometry.append(np.dot(np.asarray(rotation_matrix(rot[0], rot[1], rot[2])),atom_vec) + mol_COM)
	    else:
                for atom in mol:
                    rot_geometry.append(np.array([atom[0], atom[1], atom[2]]))
	    i+=1
        return rot_geometry

    def create_rotated_struct(self, rotated_geo, atom_types):
        ''' Creates Structure from mutated geometry'''
        struct = Structure()
        for i in range(len(rotated_geo)):
            struct.build_geo_by_atom(float(rotated_geo[i][0]), float(rotated_geo[i][1]),
                                     float(rotated_geo[i][2]), atom_types[i])
        struct.set_property('lattice_vector_a', self.A)
        struct.set_property('lattice_vector_b', self.B)
        struct.set_property('lattice_vector_c', self.C)
        struct.set_property('a', leng(self.A))
        struct.set_property('b', leng(self.B))
        struct.set_property('c', leng(self.C))
        struct.set_property('cell_vol', self.cell_vol)
        struct.set_property('crossover_type', self.cross_type)
        struct.set_property('alpha',self.alpha)
        struct.set_property('beta', self.beta)
        struct.set_property('gamma', self.gamma)
        struct.set_property('mutation_type', 'sym_rot_mol')
        return struct

class SwapMolMutation(object):
    ''' Swaps any two molecules in the unit cell'''

    def __init__(self, input_struct, num_mols, replica):
        self.input_struct = input_struct
        self.num_mols = num_mols
        self.replica = replica
        self.geometry = deepcopy(input_struct.get_geometry())
        self.ui = user_input.get_config()
        self.st_dev = self.ui.get_eval('mutation', 'stand_dev_rot')
        self.A = self.input_struct.get_property('lattice_vector_a')
        self.B = self.input_struct.get_property('lattice_vector_b')
        self.C = self.input_struct.get_property('lattice_vector_c')
        self.alpha = self.input_struct.get_property('alpha')
        self.beta = self.input_struct.get_property('beta')
        self.gamma = self.input_struct.get_property('gamma')
        self.cell_vol = self.input_struct.get_unit_cell_volume()
        self.cross_type = self.input_struct.get_property('crossover_type')

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.swapped_mutation()

    def swapped_mutation(self):
        '''Calls for a random rotation to each molecule's COM and returns a Structure'''
        temp_geo = self.geometry
        num_atom_per_mol = int(len(temp_geo)/self.num_mols)
        mol_list = [temp_geo[x:x+num_atom_per_mol] for x in range(0, len(temp_geo), num_atom_per_mol)]
        mol_list_COM = get_COM_mol_list(mol_list)
        atom_type_list = [temp_geo[i][3] for i in range(len(temp_geo))]
        swapped_geometry = self.swap_molecules(mol_list, mol_list_COM)
        swapped_struct = self.create_swapped_struct(swapped_geometry, atom_type_list)
        return swapped_struct

    def swap_molecules(self, mol_list, mol_list_COM):
        ''' Randomly rotates each molecule within gaussian dist'''
        mol_in = 0
        swap_geometry = []
        mol_info = []
        for mol in mol_list:
            mol_info.append([mol, mol_list_COM[mol_in]])
            mol_in +=1
	mol_choice1 = random.choice(mol_info)
	while True:
            mol_choice2 = random.choice(mol_info)
            if mol_choice2[0][0][0] != mol_choice1[0][0][0]:	
                break
        self.output("-- Swapping 2 molecules")
        for mol in mol_list:
            if mol[0][0] == mol_choice1[0][0][0]:
	        for atom in mol:
                    atom_vec = np.array([atom[0], atom[1], atom[2]]) - mol_choice1[1] + mol_choice2[1]
                    swap_geometry.append(atom_vec)
            elif mol[0][0] == mol_choice2[0][0][0]:
                for atom in mol:
                    atom_vec = np.array([atom[0], atom[1], atom[2]]) - mol_choice2[1] + mol_choice1[1]
                    swap_geometry.append(atom_vec)
            else:
                for atom in mol:
                    swap_geometry.append(np.array([atom[0], atom[1], atom[2]]))
        return swap_geometry

    def create_swapped_struct(self, swapped_geo, atom_types):
        ''' Creates Structure from mutated geometry'''
        struct = Structure()
        for i in range(len(swapped_geo)):
            struct.build_geo_by_atom(float(swapped_geo[i][0]), float(swapped_geo[i][1]),
                                     float(swapped_geo[i][2]), atom_types[i])
        struct.set_property('lattice_vector_a', self.A)
        struct.set_property('lattice_vector_b', self.B)
        struct.set_property('lattice_vector_c', self.C)
        struct.set_property('a', leng(self.A))
        struct.set_property('b', leng(self.B))
        struct.set_property('c', leng(self.C))
        struct.set_property('cell_vol', self.cell_vol)
        struct.set_property('crossover_type', self.cross_type)
        struct.set_property('alpha',self.alpha)
        struct.set_property('beta', self.beta)
        struct.set_property('gamma', self.gamma)
        struct.set_property('mutation_type', 'swap_mol')
        return struct

class RandomStrainMutation(object):
    '''Gives a random strain to the lattice and doesn't move the COM of the molecules'''
    def __init__(self, input_struct, num_mols, replica):
        self.ui = user_input.get_config()
        self.input_struct = input_struct
        self.replica = replica
        self.geometry = deepcopy(input_struct.get_geometry())
        self.A = np.asarray(deepcopy(input_struct.get_property('lattice_vector_a')))
        self.B = np.asarray(deepcopy(input_struct.get_property('lattice_vector_b')))
        self.C = np.asarray(deepcopy(input_struct.get_property('lattice_vector_c')))
        self.st_dev = self.ui.get_eval('mutation', 'stand_dev_strain')
        self.cross_type = self.input_struct.get_property('crossover_type')

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.strain_lat()

    def strain_lat(self):
        lat_mat = np.zeros((3,3))
        lat_mat[0] = self.A
        lat_mat[1] = self.B
        lat_mat[2] = self.C
		
	if self.ui.verbose():
		self.output("Original structure's lattices: ")
		self.output(str(self.input_struct.get_lattice_vectors()))
		self.output("This is lat_mat:")
		self.output(str(lat_mat))	

        strain_A, strain_B, strain_C = self.rand_strain(lat_mat)
	if self.ui.verbose():
		self.output("strain_A: " +str(strain_A))
		self.output("strain_B: " +str(strain_B))
		self.output("strain_C: " +str(strain_C))

        strained_struct = self.create_strained_struct(strain_A, strain_B, strain_C)
        return strained_struct

    def rand_strain(self, lat_mat):

#	rand_sleep = random.random()
#	self.output("Randomly sleeping: " + str(rand_sleep))
#	time.sleep(rand_sleep)

        strain_list = np.random.normal(scale=self.st_dev, size=6)
	self.output("Strain_list" +str(strain_list))
        strain_mat = get_strain_mat(strain_list)
        self.output("strain_mat" + str(strain_mat))
        strain = np.asarray(np.mat(lat_mat) + np.mat(strain_mat)*np.mat(lat_mat))
        self.output("strained_lattice" + str(strain))
        return strain[0], strain[1], strain[2]

    def create_strained_struct(self, lat_A, lat_B, lat_C):
        struct = Structure()
        struct.build_geo_whole(self.geometry)
        struct.set_property('lattice_vector_a', lat_A)
        struct.set_property('lattice_vector_b', lat_B)
        struct.set_property('lattice_vector_c', lat_C)
        struct.set_property('a', leng(lat_A))
        struct.set_property('b', leng(lat_B))
        struct.set_property('c', leng(lat_C))
        struct.set_property('cell_vol', np.dot(lat_A, np.cross(lat_B, lat_C)))
        struct.set_property('crossover_type', self.cross_type)
        struct.set_property('alpha', angle(lat_B, lat_C))
        struct.set_property('beta', angle(lat_A, lat_C))
        struct.set_property('gamma', angle(lat_A, lat_B))
        struct.set_property('mutation_type', 'rand_strain')
        return struct

class RandomSymmetryStrainMutation(object):
    '''
    Gives a random symmetric strain to the lattice and doesn't move the COM 
    of the molecules
    '''
    def __init__(self, input_struct, num_mols, replica):
        self.ui = user_input.get_config()
        self.input_struct = input_struct
        self.replica = replica
        self.geometry = deepcopy(input_struct.get_geometry())
        self.A = np.asarray(deepcopy(input_struct.get_property('lattice_vector_a')))
        self.B = np.asarray(deepcopy(input_struct.get_property('lattice_vector_b')))
        self.C = np.asarray(deepcopy(input_struct.get_property('lattice_vector_c')))
        self.st_dev = self.ui.get_eval('mutation', 'stand_dev_strain')
        self.cross_type = self.input_struct.get_property('crossover_type')

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.strain_lat()

    def strain_lat(self):
        lat_mat = np.zeros((3,3))
        lat_mat[0] = self.A
        lat_mat[1] = self.B
        lat_mat[2] = self.C
        if self.ui.verbose():
                self.output("Original structure's lattices: ")
                self.output(str(self.input_struct.get_lattice_vectors()))
                self.output("This is lat_mat:")
                self.output(str(lat_mat))

        strain_A, strain_B, strain_C = self.rand_sym_strain(lat_mat)
        if self.ui.verbose():
                self.output("strain_A: " +str(strain_A))
                self.output("strain_B: " +str(strain_B))
                self.output("strain_C: " +str(strain_C))

        strained_struct = self.create_strained_struct(strain_A, strain_B, strain_C)
        return strained_struct

    def rand_sym_strain(self, lat_mat):
        strain_param = np.random.normal(scale=self.st_dev, size=1)
        strain_list = strain_param*get_rand_sym_strain(lat_mat)
        strain_mat = get_strain_mat(strain_list)
        self.output("Strain parameter: " + str(strain_param).strip('[').strip(']'))
        self.output("Strain_matrix: \n" + str(strain_mat).strip('[').strip(']'))
	strain = np.asarray(np.mat(lat_mat) + np.mat(strain_mat)*np.mat(lat_mat))
        self.output("Strained_lattice:" + str(strain))
        return strain[0], strain[1], strain[2]

    def create_strained_struct(self, lat_A, lat_B, lat_C):
        struct = Structure()
        struct.build_geo_whole(self.geometry)
        struct.set_property('lattice_vector_a', lat_A)
        struct.set_property('lattice_vector_b', lat_B)
        struct.set_property('lattice_vector_c', lat_C)
        struct.set_property('a', leng(lat_A))
        struct.set_property('b', leng(lat_B))
        struct.set_property('c', leng(lat_C))
        struct.set_property('cell_vol', np.dot(lat_A, np.cross(lat_B, lat_C)))
        struct.set_property('crossover_type', self.cross_type)
        struct.set_property('alpha', angle(lat_B, lat_C))
        struct.set_property('beta', angle(lat_A, lat_C))
        struct.set_property('gamma', angle(lat_A, lat_B))
        struct.set_property('mutation_type', 'sym_strain')
        return struct

class RandomStrainMutationMoveMols(object):
    '''Gives a random strain to the lattice and moves the COM of the molecules'''
    def __init__(self, input_struct, num_mols, replica):
        self.ui = user_input.get_config()
        self.input_struct = input_struct
        self.replica = replica
        self.geometry = deepcopy(input_struct.get_geometry())
        self.num_mols = num_mols
        self.A = np.asarray(deepcopy(input_struct.get_property('lattice_vector_a')))
        self.B = np.asarray(deepcopy(input_struct.get_property('lattice_vector_b')))
        self.C = np.asarray(deepcopy(input_struct.get_property('lattice_vector_c')))
        self.st_dev = self.ui.get_eval('mutation', 'stand_dev_strain')
        self.cross_type = self.input_struct.get_property('crossover_type')

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.strain_lat_and_mols()

    def strain_lat_and_mols(self):
        temp_geo = self.geometry
        num_atom_per_mol = int(len(temp_geo)/self.num_mols)
        atom_type_list = [temp_geo[i][3] for i in range(len(temp_geo))] 
        lat_mat = set_lat_mat(self.A, self.B, self.C)
        lat_mat_f = np.linalg.inv(lat_mat)
        strain_A, strain_B, strain_C = self.rand_strain(lat_mat)
        temp_geo = self.geometry
        mol_list = [temp_geo[x:x+num_atom_per_mol] for x in range(0, len(temp_geo), num_atom_per_mol)]
        mol_list_COM = get_COM_mol_list(mol_list)
        mol_list_COM_f = get_COM_mol_list_f(lat_mat_f, np.array(mol_list_COM))
        strain_lat_mat = set_lat_mat(strain_A, strain_B, strain_C)
        strained_geometry = strain_geometry(strain_lat_mat, mol_list, mol_list_COM, mol_list_COM_f)
        strained_struct = (self.create_strained_struct(strain_A, strain_B, strain_C, 
                                                          strained_geometry, atom_type_list))
        return strained_struct

    def rand_strain(self, lat_mat):
        strain_list = np.random.normal(scale=self.st_dev, size=6)
        strain_mat = get_strain_mat(strain_list)
        self.output("strain_mat" + str(strain_mat))
        strain = np.asarray(np.mat(lat_mat) + np.mat(strain_mat)*np.mat(lat_mat))
        return strain[0], strain[1], strain[2]

    def create_strained_struct(self, lat_A, lat_B, lat_C, strained_geo, atom_types):
        ''' Creates Structure from mutated geometry'''
        struct = Structure()
        for i in range(len(strained_geo)):
            struct.build_geo_by_atom(float(strained_geo[i][0]), float(strained_geo[i][1]),
                                     float(strained_geo[i][2]), atom_types[i])
        struct.set_property('lattice_vector_a', lat_A)
        struct.set_property('lattice_vector_b', lat_B)
        struct.set_property('lattice_vector_c', lat_C)
        struct.set_property('a', leng(lat_A))
        struct.set_property('b', leng(lat_B))
        struct.set_property('c', leng(lat_C))
        struct.set_property('cell_vol', np.dot(lat_A, np.cross(lat_B, lat_C)))
        struct.set_property('crossover_type', self.cross_type)
        struct.set_property('alpha', angle(lat_B, lat_C))
        struct.set_property('beta', angle(lat_A, lat_C))
        struct.set_property('gamma', angle(lat_A, lat_B))
        struct.set_property('mutation_type', 'rand_strain_mol')
        return struct

class RandomSymmetryStrainMutationMoveMols(object):
    '''Gives a random strain to the lattice and moves the COM of the molecules'''
    def __init__(self, input_struct, num_mols, replica):
        self.ui = user_input.get_config()
        self.input_struct = input_struct
        self.replica = replica
        self.geometry = deepcopy(input_struct.get_geometry())
        self.num_mols = num_mols
        self.A = np.asarray(deepcopy(input_struct.get_property('lattice_vector_a')))
        self.B = np.asarray(deepcopy(input_struct.get_property('lattice_vector_b')))
        self.C = np.asarray(deepcopy(input_struct.get_property('lattice_vector_c')))
        self.st_dev = self.ui.get_eval('mutation', 'stand_dev_strain')
        self.cross_type = self.input_struct.get_property('crossover_type')

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.strain_lat_and_mols()

    def strain_lat_and_mols(self):
        temp_geo = self.geometry
        num_atom_per_mol = int(len(temp_geo)/self.num_mols)
        atom_type_list = [temp_geo[i][3] for i in range(len(temp_geo))] 
        lat_mat = set_lat_mat(self.A, self.B, self.C)
        lat_mat_f = np.linalg.inv(lat_mat)
        strain_A, strain_B, strain_C = self.rand_sym_strain(lat_mat)
        temp_geo = self.geometry
        mol_list = [temp_geo[x:x+num_atom_per_mol] for x in range(0, len(temp_geo), num_atom_per_mol)]
        mol_list_COM = get_COM_mol_list(mol_list)
        mol_list_COM_f = get_COM_mol_list_f(lat_mat_f, np.array(mol_list_COM))
        strain_lat_mat = set_lat_mat(strain_A, strain_B, strain_C)
        strained_geometry = strain_geometry(strain_lat_mat, mol_list, mol_list_COM, mol_list_COM_f)
        strained_struct = (self.create_strained_struct(strain_A, strain_B, strain_C, 
                                                          strained_geometry, atom_type_list))
        return strained_struct

    def rand_sym_strain(self, lat_mat):
        strain_param = np.random.normal(scale=self.st_dev, size=1)
        strain_list = strain_param*get_rand_sym_strain(lat_mat)
        strain_mat = get_strain_mat(strain_list)
        self.output("Strain parameter: " + str(strain_param))
        self.output("Strain_matrix: \n" + str(strain_mat))
        strain = np.asarray(np.mat(lat_mat) + np.mat(strain_mat)*np.mat(lat_mat))
        self.output("strained_lattice" + str(strain))
        return strain[0], strain[1], strain[2]

    def create_strained_struct(self, lat_A, lat_B, lat_C, strained_geo, atom_types):
        ''' Creates Structure from mutated geometry'''
        struct = Structure()
        for i in range(len(strained_geo)):
            struct.build_geo_by_atom(float(strained_geo[i][0]), float(strained_geo[i][1]),
                                     float(strained_geo[i][2]), atom_types[i])
        struct.set_property('lattice_vector_a', lat_A)
        struct.set_property('lattice_vector_b', lat_B)
        struct.set_property('lattice_vector_c', lat_C)
        struct.set_property('a', leng(lat_A))
        struct.set_property('b', leng(lat_B))
        struct.set_property('c', leng(lat_C))
        struct.set_property('cell_vol', np.dot(lat_A, np.cross(lat_B, lat_C)))
        struct.set_property('crossover_type', self.cross_type)
        struct.set_property('alpha', angle(lat_B, lat_C))
        struct.set_property('beta', angle(lat_A, lat_C))
        struct.set_property('gamma', angle(lat_A, lat_B))
        struct.set_property('mutation_type', 'sym_strain_mol')
        return struct

#---- Functions Shared Between Several Mutation Classes ----#
def rotation_matrix(theta, psi, phi):
    Rxyz = np.matrix([ ((np.cos(theta) * np.cos(psi)),
                (-np.cos(phi) * np.sin(psi)) + (np.sin(phi) * np.sin(theta) * np.cos(psi)),
                (np.sin(phi) * np.sin(psi)) + (np.cos(phi) * np.sin(theta) * np.cos(psi))),

                ((np.cos(theta) * np.sin(psi)),
                (np.cos(phi) * np.cos(psi)) + (np.sin(phi) * np.sin(theta) * np.sin(psi)),
                (-np.sin(phi) * np.cos(psi)) + (np.cos(phi) * np.sin(theta) * np.sin(psi))),

                ((-np.sin(theta)),
                (np.sin(phi) * np.cos(theta)),
                (np.cos(phi) * np.cos(theta)))])
    return Rxyz

def get_strain_mat(strain_list):
    sl = strain_list
    s_mat = np.zeros((3,3))
    s_mat[0][0] = sl[0]
    s_mat[0][1] = sl[5]/2.
    s_mat[0][2] = sl[4]/2.
    s_mat[1][0] = sl[5]/2.
    s_mat[1][1] = sl[1]
    s_mat[1][2] = sl[3]/2.
    s_mat[2][0] = sl[4]/2.
    s_mat[2][1] = sl[3]/2.
    s_mat[2][2] = sl[2]
    return s_mat

def get_rand_sym_strain(lat_mat):
    strain_dict={                       \
    '1':[ 1., 1., 1., 0., 0., 0.],\
    '2':[ 1., 0., 0., 0., 0., 0.],\
    '3':[ 0., 1., 0., 0., 0., 0.],\
    '4':[ 0., 0., 1., 0., 0., 0.],\
    '5':[ 0., 0., 0., 2., 0., 0.],\
    '6':[ 0., 0., 0., 0., 2., 0.],\
    '7':[ 0., 0., 0., 0., 0., 2.],\
    '8':[ 1., 1., 0., 0., 0., 0.],\
    '9':[ 1., 0., 1., 0., 0., 0.],\
    '10':[ 1., 0., 0., 2., 0., 0.],\
    '11':[ 1., 0., 0., 0., 2., 0.],\
    '12':[ 1., 0., 0., 0., 0., 2.],\
    '13':[ 0., 1., 1., 0., 0., 0.],\
    '14':[ 0., 1., 0., 2., 0., 0.],\
    '15':[ 0., 1., 0., 0., 2., 0.],\
    '16':[ 0., 1., 0., 0., 0., 2.],\
    '17':[ 0., 0., 1., 2., 0., 0.],\
    '18':[ 0., 0., 1., 0., 2., 0.],\
    '19':[ 0., 0., 1., 0., 0., 2.],\
    '20':[ 0., 0., 0., 2., 2., 0.],\
    '21':[ 0., 0., 0., 2., 0., 2.],\
    '22':[ 0., 0., 0., 0., 2., 2.],\
    '23':[ 0., 0., 0., 2., 2., 2.],\
    '24':[-1., .5, .5, 0., 0., 0.],\
    '25':[ .5,-1., .5, 0., 0., 0.],\
    '26':[ .5, .5,-1., 0., 0., 0.],\
    '27':[ 1.,-1., 0., 0., 0., 0.],\
    '28':[ 1.,-1., 0., 0., 0., 2.],\
    '29':[ 0., 1.,-1., 0., 0., 2.],}
    rand_choice = str(np.random.randint(1,29))
    return np.array(strain_dict[rand_choice])

def get_COM_mol_list(mol_list):
    ''' Returns COM as np array '''
    COM_mol_list = []
    for mol in mol_list:
        N = len(mol)
        rsum = [0, 0, 0]
        for atom in mol:
            rsum[0] +=atom[0]
            rsum[1] +=atom[1]
            rsum[2] +=atom[2]
        COM = np.array(rsum)/N
        COM_mol_list.append(COM)
    return np.asarray(COM_mol_list)

def get_COM_mol_list_f(lat_mat_f, COM_mol_list):
    ''' Returns COM as np array '''
    COM_mol_list_f = []
    for COM in COM_mol_list:
        COM_f = np.dot(lat_mat_f, COM)
        COM_mol_list_f.append(COM_f)
    return COM_mol_list_f

def strain_geometry(strain_lat_mat, mol_list, mol_list_COM, mol_list_COM_f):
    A = strain_lat_mat[0]
    B = strain_lat_mat[1]
    C = strain_lat_mat[2]
    strain_list_COM = []
    strained_geometry = []
    for mol_f in mol_list_COM_f:
        strain_COM = mol_f[0]* A + mol_f[1] * B + mol_f[2] * C
        strain_list_COM.append(strain_COM)
    COM_diff = mol_list_COM - strain_list_COM
    i = 0
    for mol in mol_list:
	xyz = []
	for j in range(len(mol)):
            xyz.append((mol[j][0], mol[j][1], mol[j][2]))
        for atom in xyz:
            strained_geometry.append((atom - COM_diff[i])) 
        i += 1
    return np.asarray(strained_geometry)

def set_lat_mat(lat_A, lat_B, lat_C):
    lat_mat = np.zeros((3,3))
    lat_mat[0] = lat_A
    lat_mat[1] = lat_B
    lat_mat[2] = lat_C
    return lat_mat
    
def angle(v1,v2):
    numdot = np.dot(v1,v2)
    anglerad = np.arccos(numdot/(leng(v1)*leng(v2)))
    angledeg = anglerad*180/np.pi
    return (angledeg)

def leng(v):
        length = np.linalg.norm(v)
        return length
