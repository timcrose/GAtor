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
    output.local_message(input_struct.get_geometry_atom_format(), replica)
    mutate_obj = select_mutator(input_struct, num_mols, replica)
    mutated_struct = mutate_obj.mutate()
    output.local_message(mutated_struct.get_geometry_atom_format(), replica)
    return mutated_struct
   
def select_mutator(input_struct, num_mols, replica):
    '''
    In this mutation implementation, there are several classes, each performing a 
    different mutation. This method is responsible for reading the preferences set
    by the user and selecting which mutation to empoly, or no mutation at all.
    Expects: Structure, number of molecules per cell, replica name
    Returns: Mutation Class
    '''
    mutation_list = ["Trans_mol","Rot_mol","Strain_rand","Strain_rand_mols","Strain_sym_mols","Strain_sym"]
    mutation_list = ["Strain_rand"]

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
        self.cell_vol = self.input_struct.get_property('cell_vol')
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
            self.output(rand_disp)
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
        self.cell_vol = self.input_struct.get_property('cell_vol')
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
        rotated_geometry = self.rotate_molecules(mol_list, self.st_dev)
        rotated_struct = self.create_rotated_struct(rotated_geometry, atom_type_list)
        return rotated_struct
	
    def rotate_molecules(self, mol_list, st_dev):
        ''' Randomly rotates each molecule within gaussian dist'''
        rot_geometry = []
        for mol in mol_list:
            rand_vec = np.random.normal(scale=st_dev,size=3)
            for atom in mol:
                atom_vec = np.array([atom[0], atom[1], atom[2]])
                rot_geometry.append(np.dot(np.asarray(self.rotation_matrix(rand_vec)),atom_vec))
        return rot_geometry

    def rotation_matrix(self, rand_vec):
        theta= (np.pi/180)*rand_vec[0]
        psi = (np.pi/180)*rand_vec[1]
        phi= (np.pi/180)*rand_vec[2]
        self.output("Random Rotation Angles:  ")
        self.output(str(theta*180/np.pi)+' '+str(psi*180/np.pi)+' '+str(phi*180/np.pi))
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

class RandomStrainMutation(object):
    '''Gives a random strain to the lattice and doesn't move the COM of the molecules'''
    def __init__(self, input_struct, target_stoic, replica):
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
        strain_A, strain_B, strain_C = self.rand_strain(lat_mat)
        strained_struct = self.create_strained_struct(strain_A, strain_B, strain_C)
        return strained_struct

    def rand_strain(self, lat_mat):
        strain_list = np.random.normal(scale=self.st_dev, size=6)
        strain_mat = get_strain_mat(strain_list)
        self.output("strain_mat" + str(strain_mat))
        strain_A = np.dot(lat_mat.transpose()[0], strain_mat)
        strain_B = np.dot(lat_mat.transpose()[1], strain_mat)
        strain_C = np.dot(lat_mat.transpose()[2], strain_mat)
        return strain_A, strain_B, strain_C

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
        struct.set_property('mutation_type', 'Strain_rand')
        return struct

class RandomSymmetryStrainMutation(object):
    '''
    Gives a random symmetric strain to the lattice and doesn't move the COM 
    of the molecules
    '''
    def __init__(self, input_struct, target_stoic, replica):
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
        strain_A, strain_B, strain_C = self.rand_sym_strain(lat_mat)
        strained_struct = self.create_strained_struct(strain_A, strain_B, strain_C)
        return strained_struct

    def rand_sym_strain(self, lat_mat):
        strain_param = np.random.normal(scale=self.st_dev, size=1)
        strain_list = strain_param*get_rand_sym_strain(lat_mat)
        strain_mat = get_strain_mat(strain_list)
        seld.output("Strain parameter: " + str(strain_param))
        self.output("Strain_matrix: " + str(strain_mat))
        strain_A = np.dot(lat_mat.transpose()[0], strain_mat)
        strain_B = np.dot(lat_mat.transpose()[1], strain_mat)
        strain_C = np.dot(lat_mat.transpose()[2], strain_mat)
        return strain_A, strain_B, strain_C

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
        struct.set_property('mutation_type', 'Strain_sym')
        return struct


#---- Functions Shared Between Several Mutation Classes ----#
def get_strain_mat(strain_list):
    sl = strain_list
    s_mat = np.zeros((3,3))
    s_mat[0][0] = 1.0 + sl[0]
    s_mat[0][1] = sl[5]/2.
    s_mat[0][2] = sl[4]/2.
    s_mat[1][0] = sl[5]/2.
    s_mat[1][1] = 1.0 + sl[1]
    s_mat[1][2] = sl[3]/2.
    s_mat[2][0] = sl[4]/2.
    s_mat[2][1] = sl[3]/2.
    s_mat[2][2] = 1.0 + sl[2]
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
    '29':[ 0., 1.,-1., 0., 0., 2.],\
    '30':[ 1., 2., 3., 4., 5., 6.],\
    '31':[-2., 1., 4.,-3., 6.,-5.],\
    '32':[ 3.,-5.,-1., 6., 2.,-4.],\
    '33':[-4.,-6., 5., 1.,-3., 2.],\
    '34':[ 5., 4., 6.,-2.,-1.,-3.],\
    '35':[-6., 3.,-2., 5.,-4., 1.]}
    rand_choice = str(np.random.randint(1,35))
    return np.array(strain_dict[rand_choice])

def angle(v1,v2):
    numdot = np.dot(v1,v2)
    anglerad = np.arccos(numdot/(leng(v1)*leng(v2)))
    angledeg = anglerad*180/np.pi
    return (angledeg)

def leng(v):
        length = np.linalg.norm(v)
        return length

def center_geometry(geometry):
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

