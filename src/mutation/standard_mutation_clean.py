'''
@author: farren
'''
from __future__ import division

from copy import deepcopy
import numpy as np
import random
import time
import itertools
import math
import sys

sys.path.append("/lustre/project/nmarom/fcurtis/NEW_gator/latest_git/src")
from core import user_input,output
from structures.structure import Structure
from structures import structure_handling
from utilities import space_group_utils as sgu
from pymatgen import Molecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer as pga


def main(struct, replica, mut_choice):
    '''
    Args: Single Structure() to mutate, the replica name running the 
    crossover instance.
    Returns: A single Structure() if mutation is successful or False if 
    mutation fails 
    '''
    input_struct = deepcopy(struct)
    vol = struct.get_unit_cell_volume()
    output.local_message("Input Structure's Volume: %s" % (vol), replica)
    ui = user_input.get_config()
    num_mols = ui.get_eval('run_settings', 'num_molecules')

    Mutate = select_mutator(input_struct, num_mols, replica, mut_choice)
    mutated_struct = Mutate.mutate()
    return mutated_struct


def select_mutator(input_struct, num_mols, replica, mut_choice):
    '''
    In this mutation implementation, there are several classes, each performing a 
    different mutation. This method selects which mutation to employ
    Expects: Structure, number of molecules per cell, replica name
    Returns: Mutation Object
    '''
    message = "-- Mutation Choice: %s" % (mut_choice)
    output.local_message(message, replica)

    # Translation Mutation Classes
    if mut_choice == "Rand_trans":
        mutator = RandomTranslationMutation(input_struct, num_mols, replica)
    elif mut_choice == "Frame_trans":
        mutator = RandomTranslationFrameMutation(input_struct, num_mols, replica)

    # Rotation Mutation Classes
    elif mut_choice == "Rand_rot":
        mutator = RandomRotationMutation(input_struct, num_mols, replica)
    elif mut_choice == "Frame_rot":
        mutator = RandomRotationFrameMutation(input_struct, num_mols, replica)

    # Permutation Mutation Classes
    elif mut_choice =="Permute_mol": 
        mutator = PermutationMutation(input_struct, num_mols, replica)
    elif mut_choice == "Permute_rot":
        mutator = PermutationRotationFrameMutation(input_struct, num_mols, replica)
    elif mut_choice == "Permute_ref":
        mutator = PermutationReflectionMutation(input_struct, num_mols, replica)

    # Strain Mutations
    elif mut_choice == "Rand_strain":
        mutator = RandomStrainMutation(input_struct, num_mols, replica)
    elif mut_choice == "Sym_strain":
        mutator = RandomSymmetryStrainMutation(input_struct, num_mols, replica)
    elif mut_choice =="Vol_strain":
	    mutator = RandomVolumeStrainMutation(input_struct, num_mols, replica)
    elif mut_choice =="Angle_strain":
        mutator = RandomCellAngleMutation(input_struct, num_mols, replica)
    return mutator

class RandomTranslationMutation(object):
    ''' 
    Gives a random rotation to the COM of randomly-selected 
    molecules in the unit cell.
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
        self.mol_list = self.input_struct.get_molecules(self.num_mols)
        self.atom_types = self.input_struct.get_atom_types()

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.random_translation()

    def random_translation(self):
        '''
        Calls for random to rotation to random molecule's COM 
        and returns a Structure
        '''
        translated_geometry = self.translate_molecules()
        translated_struct = self.build_child_struct(translated_geometry)
        return translated_struct

    def translate_molecules(self):
        ''' Randomly rotates all molecules COM by some random angle'''
        # Select which molecules are receiving a rotation
        mol_indices = [i for i in range(self.num_mols)]
        while True:
            select_mols = np.random.choice(mol_indices) + 1
            if select_mols % 2 == 0:
                break
        trans_indices = sorted(np.random.choice(mol_indices, size=select_mols, replace=False))

        self.output("-- Performing translation on molecules %s " % (trans_indices))
        rot_geometry = []
        for i in range(len(self.mol_list)):
            mol = self.mol_list[i]
            if i in trans_indices:
                rand_disp = np.random.standard_normal(3) * self.st_dev
                self.output("-- Random_displacement: %s to molecule %i" %(rand_disp, i))
                for atom in mol:
                    atom_vec = np.array([atom[0], atom[1], atom[2]]) + rand_disp
                    rot_geometry.append(atom_vec)
            else:
                for atom in mol:
                    atom_vec = np.array([atom[0], atom[1], atom[2]])
                    rot_geometry.append(atom_vec)
        return rot_geometry

    def build_child_struct(self, rotated_geo):
        ''' Creates Structure from mutated geometry'''
        struct = Structure()
        for i in range(len(rotated_geo)):
            struct.build_geo_by_atom(float(rotated_geo[i][0]), float(rotated_geo[i][1]),
                                     float(rotated_geo[i][2]), self.atom_types[i])
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
        struct.set_property('mutation_type', 'Rand_trans')
        return struct

class RandomTranslationFrameMutation(object):
    '''
    This mutation gives a random translation to the COM of one or more of the molecules 
    in the unit cell in their inertial frame.
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
        self.orientation_info = []
        self.mol_list = self.input_struct.get_molecules(self.num_mols)
        self.atom_types = self.input_struct.get_atom_types()

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.random_translation()

    def random_translation(self):
        '''Calls for a random translation to each molecule's COM and returns a Structure'''
        mol_index = 0
        rand_disp = np.random.standard_normal(3) * self.st_dev
        self.output("-- Random_displacement: %s" %(rand_disp))
        for mol in self.mol_list:
            COM = self.get_COM_frac(mol, [self.A, self.B, self.C])
            centered_mol, types = self.get_centered_molecule(mol)
            rot, aligned_mol = self.translate_get_orientation_info(centered_mol, mol_index, types, rand_disp)
            self.orientation_info.append([rot, COM, aligned_mol])
        lattice, child_coords = self.reconstruct_child_lat_coords(self.orientation_info)
        mutated_struct = self.build_child_struct(child_coords, lattice, self.atom_types)
        return mutated_struct
    
    def reconstruct_child_lat_coords(self, orientation_info):
        '''
        Reconstructs the child's atomic positions and lattice
        from its inherited genes 
        '''
        lattice = [self.A, self.B, self.C]
        child_coordinates = []
        COM_xyz = []
        for mol_info in orientation_info:
            mol_coords = []
            rot, COM, centered_mol = mol_info
            COM_xyz = np.dot(lattice, COM)
            for atom in centered_mol:
                mol_coords.append(np.dot(rot, 
                        np.array(atom).reshape((3,1))).tolist())
            for coord in mol_coords:
                new_coords = [coord[0][0] + COM_xyz[0], 
                              coord[1][0] + COM_xyz[1], 
                              coord[2][0] + COM_xyz[2]]
                child_coordinates.append(new_coords)
        return lattice, np.array(child_coordinates)

    def get_centered_molecule(self, mol):
        ''' Returns array of COM and centered geometry for each molecule '''
        atoms = []
        types = []
        centered_mol = []
        for atom in mol:
            atoms.append([float(atom[0]), float(atom[1]), float(atom[2])])
            types.append(atom[3])
        molp = Molecule(types, atoms)
        centered_molp = molp.get_centered_molecule()
        for site in centered_molp.sites:
            centered_mol.append(list(site._coords))
        return centered_mol, types

    def get_COM_frac(self, mol, lattice_vectors):
        ''' Returns array of COM for a given molecule '''
        atoms = []
        types = []
        for atom in mol:
            atoms.append([float(atom[0]), float(atom[1]), float(atom[2])])
            types.append(atom[3])
        molp = Molecule(types, atoms)
        COM = molp.center_of_mass #cartesian
        latinv = np.linalg.inv(lattice_vectors)
        frac_COM = np.dot(latinv, COM)
        return frac_COM

    def translate_get_orientation_info(self, mol, mol_index, types, rand_disp):
        ''' Computes the principal axes for each molecule and the corresponding
            rotation matrices. Translates each molecule in inertial reference frame
            Returns: rotations about z, y, x axis and molecule with principal axes 
            aligned. '''
        atoms = []
        for atom in mol:
            atoms.append([float(atom[0]), float(atom[1]), float(atom[2])])
           # types.append(atom[3])
        molp = Molecule(types, atoms)
        centered_molp = molp.get_centered_molecule()

        #compute principal axes and eignvalues
        PGA = pga(centered_molp)
        ax1, ax2, ax3 = PGA.principal_axes
        eig1, eig2, eig3 = PGA.eigvals
        eigen_sys = [[eig1, ax1],[eig2, ax2],[eig3, ax3]]
        sort_eig = sorted(eigen_sys)

        #Construct rotation matrix and its inverse
        rot = np.column_stack((sort_eig[0][1], sort_eig[1][1], sort_eig[2][1]))
        rot_trans = np.linalg.inv(np.array(rot))

        #Aligned geometry
        centered_sites = []
        aligned_mol = []
        for site in centered_molp.sites:
            centered_sites.append(list(site._coords))
        for atom in centered_sites:
            new_atom = np.dot(rot_trans, np.array(atom))
            new_atom = new_atom + rand_disp #random translation
            aligned_mol.append(new_atom)
        return rot, aligned_mol

    def angle(self,v1,v2):
        numdot = np.dot(v1,v2)
        anglerad = np.arccos(numdot/(leng(v1)*leng(v2)))
        angledeg = anglerad*180/np.pi
        return (anglerad)

    def build_child_struct(self, child_geometry, child_lattice, atom_types):
        ''' Creates a Structure() for the child to be added to the collection'''
        child_A, child_B, child_C = child_lattice
        temp_vol = np.dot(np.cross(child_A, child_B), child_C)
        alpha = self.angle(child_B, child_C)*180./np.pi
        beta = self.angle(child_C, child_A)*180./np.pi
        gamma = self.angle(child_A, child_B)*180./np.pi
        struct = Structure()
        for i in range(len(child_geometry)):
            struct.build_geo_by_atom(float(child_geometry[i][0]), float(child_geometry[i][1]),
                                     float(child_geometry[i][2]), atom_types[i])
        struct.set_property('lattice_vector_a', child_A)
        struct.set_property('lattice_vector_b', child_B)
        struct.set_property('lattice_vector_c', child_C)
        struct.set_property('a', np.linalg.norm(child_A))
        struct.set_property('b', np.linalg.norm(child_B))
        struct.set_property('c', np.linalg.norm(child_C))
        struct.set_property('cell_vol', temp_vol)
        struct.set_property('mutation_type', "Frame_trans")
        struct.set_property('alpha',alpha)
        struct.set_property('beta', beta)
        struct.set_property('gamma', gamma)
        return struct

class RandomRotationMutation(object):
    ''' 
    Gives a random rotation to the COM of randomly-selected 
    molecules in the unit cell.
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
        self.mol_list = self.input_struct.get_molecules(self.num_mols)
        self.atom_types = self.input_struct.get_atom_types()

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.random_rotation()

    def random_rotation(self):
        '''
        Calls for random to rotation to random molecule's COM 
        and returns a Structure
        '''
        rotated_geometry = self.rotate_molecules()
        rotated_struct = self.build_child_struct(rotated_geometry)
        return rotated_struct

    def rotate_molecules(self):
        ''' Randomly rotates all molecules COM by some random angle'''
        mol_list_COM = get_COM_mol_list(self.mol_list)

        # Select which molecules are receiving a rotation
        mol_indices = [i for i in range(self.num_mols)]
        while True:
            select_mols = np.random.choice(mol_indices) + 1
            if select_mols % 2 == 0:
                break
        rot_indices = sorted(np.random.choice(mol_indices, size=select_mols, replace=False))

        self.output("-- Performing rotation on molecules %s " % (rot_indices))
        rot_geometry = []
        for i in range(len(self.mol_list)):
            mol = self.mol_list[i]
            if i in rot_indices:
                rand_rot = np.random.normal(scale=self.st_dev,size=3)
                self.output("--Random Rotation %s to molecule %i" % (rand_rot, i))
                theta, psi, phi = [(np.pi/180)*rand_rot[i] for i in range(3)]
                rot_mat = rotation_matrix(theta, psi, phi)
                mol_COM = np.array(mol_list_COM[i])
                for atom in mol:
                    atom_vec = np.array([atom[0], atom[1], atom[2]]) - mol_COM
                    rot_geometry.append(np.dot(np.asarray(rot_mat),atom_vec) + mol_COM)
            else:
                for atom in mol:
                    atom_vec = np.array([atom[0], atom[1], atom[2]])
                    rot_geometry.append(atom_vec)
        return rot_geometry

    def build_child_struct(self, rotated_geo):
        ''' Creates Structure from mutated geometry'''
        struct = Structure()
        for i in range(len(rotated_geo)):
            struct.build_geo_by_atom(float(rotated_geo[i][0]), float(rotated_geo[i][1]),
                                     float(rotated_geo[i][2]), self.atom_types[i])
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
        struct.set_property('mutation_type', 'Rand_rot')
        return struct

class RandomRotationFrameMutation(object):
    '''
    This mutation gives a random rotation to the COM of one or more of the molecules 
    in the unit cell in their inertial frame.
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
        self.orientation_info = []
        self.mol_list = self.input_struct.get_molecules(self.num_mols)
        self.atom_types = self.input_struct.get_atom_types()

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.random_rotation()

    def random_rotation(self):
        '''Calls for a random translation to each molecule's COM and returns a Structure'''
        mol_index = 0
        for mol in self.mol_list:
            COM = self.get_COM_frac(mol, [self.A, self.B, self.C])
            centered_mol, types = self.get_centered_molecule(mol)
            z, y, x, aligned_mol = self.get_orientation_info(centered_mol, mol_index, types)
            self.orientation_info.append([z, y, x, COM, aligned_mol])
        lattice, child_coords = self.random_rotation_frame(self.orientation_info)
        mutated_struct = self.build_child_struct(child_coords, lattice, self.atom_types)
        return mutated_struct

    def get_COM_frac(self, mol, lattice_vectors):
        ''' Returns array of COM for a given molecule '''
        atoms = []
        types = []
        for atom in mol:
            atoms.append([float(atom[0]), float(atom[1]), float(atom[2])])
            types.append(atom[3])
        molp = Molecule(types, atoms)
        COM = molp.center_of_mass #cartesian
        latinv = np.linalg.inv(lattice_vectors)
        frac_COM = np.dot(latinv, COM)
        return frac_COM

    def get_centered_molecule(self, mol):
        ''' Returns array of COM and centered geometry for each molecule '''
        atoms = []
        types = []
        centered_mol = []
        for atom in mol:
            atoms.append([float(atom[0]), float(atom[1]), float(atom[2])])
            types.append(atom[3])
        molp = Molecule(types, atoms)
        centered_molp = molp.get_centered_molecule()
        for site in centered_molp.sites:
            centered_mol.append(list(site._coords))
        return centered_mol, types

    def random_rotation_frame(self, orientation_info):
        '''
        Reconstructs the child's atomic positions and lattice
        from its inherited genes 
        '''
        rand_vec = np.random.normal(scale=self.st_dev,size=3)
        self.output("-- Random Rotation %s" % (rand_vec))
        dx, dy, dz = [rand_vec[i] for i in range(3)]
        self.output("-- Mutation rotation angles: %s %s %s" %(dx, dy, dz))

        lattice = [self.A, self.B, self.C]
        child_coordinates = []
        COM_xyz = []
        for mol_info in orientation_info:
            mol_coords = []
            z, y, x, COM, centered_mol = mol_info
            rot_from_euler = euler2mat(z + dz, y + dy, x + dx)
            COM_xyz = np.dot(lattice, COM)
            for atom in centered_mol:
                mol_coords.append(np.dot(rot_from_euler, 
                                  np.array(atom).reshape((3,1))).tolist())
            for coord in mol_coords:
                new_coords = [coord[0][0] + COM_xyz[0], 
                              coord[1][0] + COM_xyz[1], 
                              coord[2][0] + COM_xyz[2]]
                child_coordinates.append(new_coords)
        return lattice, np.array(child_coordinates)


    def get_orientation_info(self, mol, mol_index, types):
        ''' 
        Computes the principal axes for each molecule and the corresponding
        rotation matrices.
        Returns: rotations about z, y, x axis and molecule with principal axes 
        aligned. 
        '''
        # Center molecule at origin about center of mass
        atoms = []
        for atom in mol:
            atoms.append([float(atom[0]), float(atom[1]), float(atom[2])])
        molp = Molecule(types, atoms)
        centered_molp = molp.get_centered_molecule()

        #Xompute principal axes and eignvalues
        PGA = pga(centered_molp)
        ax1, ax2, ax3 = PGA.principal_axes
        eig1, eig2, eig3 = PGA.eigvals
        eigen_sys = [[eig1, ax1],[eig2, ax2],[eig3, ax3]]
        sort_eig = sorted(eigen_sys)

        #Construct rotation matrix and its inverse
        rot = np.column_stack((sort_eig[0][1], sort_eig[1][1], sort_eig[2][1]))
        rot_trans = np.linalg.inv(np.array(rot))

        #Aligned geometry
        centered_sites = []
        aligned_mol = []
        for site in centered_molp.sites:
            centered_sites.append(list(site._coords))
        for atom in centered_sites: #chiral
            new_atom = np.dot(rot_trans, np.array(atom))
            aligned_mol.append(new_atom)
        
        #Compute euler angles of orientation
        z, y, x =  mat2euler(rot)
        return z, y, x, aligned_mol

    def angle(self,v1,v2):
        numdot = np.dot(v1,v2)
        anglerad = np.arccos(numdot/(leng(v1) * leng(v2)))
        angledeg = anglerad*180/np.pi
        return (anglerad)

    def build_child_struct(self, child_geometry, child_lattice, atom_types):
        ''' Creates a Structure() for the child to be added to the collection'''
        child_A, child_B, child_C = child_lattice
        temp_vol = np.dot(np.cross(child_A, child_B), child_C)
        alpha = self.angle(child_B, child_C)*180./np.pi
        beta = self.angle(child_C, child_A)*180./np.pi
        gamma = self.angle(child_A, child_B)*180./np.pi
        struct = Structure()
        for i in range(len(child_geometry)):
            struct.build_geo_by_atom(float(child_geometry[i][0]), float(child_geometry[i][1]),
                                     float(child_geometry[i][2]), atom_types[i])
        struct.set_property('lattice_vector_a', child_A)
        struct.set_property('lattice_vector_b', child_B)
        struct.set_property('lattice_vector_c', child_C)
        struct.set_property('a', np.linalg.norm(child_A))
        struct.set_property('b', np.linalg.norm(child_B))
        struct.set_property('c', np.linalg.norm(child_C))
        struct.set_property('cell_vol', temp_vol)
        struct.set_property('mutation_type', "Frame_rot")
        struct.set_property('alpha',alpha)
        struct.set_property('beta', beta)
        struct.set_property('gamma', gamma)
        return struct

class PermutationMutation(object):
    ''' Chnages the COM of molecules in the structure'''

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
        self.mol_list = self.input_struct.get_molecules(self.num_mols)
        self.atom_types = self.input_struct.get_atom_types()

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.permutation_mutation()

    def permutation_mutation(self):
        '''Permutes the location of each molecule eg. mol 1234 -> 4123'''
        mol_list_COM = get_COM_mol_list(self.mol_list)
        permuted_geometry = self.permute_molecules(self.mol_list, mol_list_COM)
        permutated_struct = self.create_permutated_struct(permuted_geometry, self.atom_types)
        return permutated_struct

    def permute_molecules(self, mol_list, mol_list_COM):
        ''' Randomly permutes the location of each molecule'''

        self.output("--- Permuting COM of molecules")
        permuted_geometry = []
        perm = [i for i in range(len(mol_list))]
        np.random.shuffle(perm)
        self.output("-- New ordering %s" %(perm))
 
        for i in range(len(mol_list)):
            mol = mol_list[i]
            for atom in mol:
                xyz = [atom[0], atom[1], atom[2]]
                geo = xyz - mol_list_COM[i] + mol_list_COM[perm[i]]
                permuted_geometry.append(np.array(geo))
        return permuted_geometry

    def create_permutated_struct(self, swapped_geo, atom_types):
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
        struct.set_property('mutation_type', 'Permute_mol')
        return struct

class PermutationRotationFrameMutation(object):
    '''
    This mutation gives a random rotation to the COM of one or more of the molecules 
    in the unit cell in their inertial frame.
    '''
    def __init__(self, input_struct, num_mols, replica):
        self.input_struct = input_struct
        self.num_mols = num_mols
        self.replica = replica
        self.geometry = deepcopy(input_struct.get_geometry())
        self.ui = user_input.get_config()
        self.st_dev = self.ui.get_eval('mutation', 'stand_dev_rot')
        self.cross_type = self.input_struct.get_property('crossover_type')
        self.A = self.input_struct.get_property('lattice_vector_a')
        self.B = self.input_struct.get_property('lattice_vector_b')
        self.C = self.input_struct.get_property('lattice_vector_c')
        self.orientation_info = []
        self.mol_list = self.input_struct.get_molecules(self.num_mols)
        self.atom_types = self.input_struct.get_atom_types()

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.random_permutation_rotation_mutation()

    def random_permutation_rotation_mutation(self):
        '''
        Calls for a random permutation of the molecules followed by
        a rotation about each molecules principal axes
        '''
        mol_index = 0
        for mol in self.mol_list:
            COM = self.get_COM_frac(mol, [self.A, self.B, self.C])
            centered_mol, mol_types = self.get_centered_molecule(mol)
            z, y, x, aligned_mol = self.get_orientation_info(centered_mol, mol_index, mol_types)
            self.orientation_info.append([z, y, x, COM, aligned_mol])
        lattice, child_coords = self.random_permutation_rotation(self.orientation_info)
        mutated_struct = self.build_child_struct(child_coords, lattice, self.atom_types)
        return mutated_struct

    def get_COM_frac(self, mol, lattice_vectors):
        ''' 
        Returns array of COMs in fractional coordinates 
        for a given molecule and lattice
        '''
        atoms = []
        types = []
        for atom in mol:
            atoms.append([float(atom[0]), float(atom[1]), float(atom[2])])
            types.append(atom[3])
        molp = Molecule(types, atoms)
        COM = molp.center_of_mass #cartesian
        latinv = np.linalg.inv(lattice_vectors)
        frac_COM = np.dot(latinv, COM)
        return frac_COM

    def get_centered_molecule(self, mol):
        ''' 
        Returns array of COM and centered at the origin 
        geometry for each molecule 
        '''
        atoms = []
        types = []
        centered_mol = []
        for atom in mol:
            atoms.append([float(atom[0]), float(atom[1]), float(atom[2])])
            types.append(atom[3])
        molp = Molecule(types, atoms)
        centered_molp = molp.get_centered_molecule()
        for site in centered_molp.sites:
            centered_mol.append(list(site._coords))
        return centered_mol, types

    def random_permutation_rotation(self, orientation_info):
        '''
        Performs a permuation of the molecule's COMs and
        a random rotation in each molecules reference frame
        '''
        # Select Permutation
        self.output("-- Permuting COM of molecules")
        perm = [i for i in range(len(self.mol_list))]
        np.random.shuffle(perm)
        self.output("-- New ordering %s" %(perm))

        # Random Rotation and Permutation
        rand_vec = np.random.normal(scale=self.st_dev,size=3)
        dx, dy, dz  = [rand_vec[i] for i in range(len(rand_vec))]
        self.output("-- Mutation rotation angles: %s %s %s" %(dx, dy, dz))

        lattice = [self.A, self.B, self.C]
        child_coordinates = []
        COMs = []
        for mol_info in orientation_info:
            mol_coords = []
            z, y, x, COM, centered_mol = mol_info
            rot_from_euler = euler2mat(z + dz, y + dy, x + dx)
            COM_xyz = np.dot(lattice, COM)
            COMs.append(COM_xyz)

        i = 0
        for mol_info in orientation_info:
            mol_coords = []
            z, y, x, COM, centered_mol = mol_info
            rot_from_euler = euler2mat(z + dz, y + dy, x + dx)

            for atom in centered_mol:
                mol_coords.append(np.dot(rot_from_euler,
                                  np.array(atom).reshape((3,1))).tolist())
            for coord in mol_coords:
                new_coords = [coord[0][0] + COMs[perm[i]][0],
                              coord[1][0] + COMs[perm[i]][1],
                              coord[2][0] + COMs[perm[i]][2]]
                child_coordinates.append(new_coords)
            i+=1
        return lattice, np.array(child_coordinates)

    def get_orientation_info(self, mol, mol_index, types):
        ''' 
        Computes the principal axes for each molecule and the corresponding
        rotation matrices.
        Returns: rotations about z, y, x axis and molecule with principal axes 
        aligned. 
        '''
        # Center molecule at origin about center of mass
        atoms = []
        for atom in mol:
            atoms.append([float(atom[0]), float(atom[1]), float(atom[2])])
        molp = Molecule(types, atoms)
        centered_molp = molp.get_centered_molecule()

        #Xompute principal axes and eignvalues
        PGA = pga(centered_molp)
        ax1, ax2, ax3 = PGA.principal_axes
        eig1, eig2, eig3 = PGA.eigvals
        eigen_sys = [[eig1, ax1],[eig2, ax2],[eig3, ax3]]
        sort_eig = sorted(eigen_sys)

        #Construct rotation matrix and its inverse
        rot = np.column_stack((sort_eig[0][1], sort_eig[1][1], sort_eig[2][1]))
        rot_trans = np.linalg.inv(np.array(rot))

        #Aligned geometry
        centered_sites = []
        aligned_mol = []
        for site in centered_molp.sites:
            centered_sites.append(list(site._coords))
        for atom in centered_sites: #chiral
            new_atom = np.dot(rot_trans, np.array(atom))
            aligned_mol.append(new_atom)

        #Compute euler angles of orientation
        z, y, x =  mat2euler(rot)
        #self.output("-- Molecule rotation angles: %s %s %s" %  (z, y, x))
        return z, y, x, aligned_mol

    def angle(self,v1,v2):
        numdot = np.dot(v1,v2)
        anglerad = np.arccos(numdot/(leng(v1) * leng(v2)))
        angledeg = anglerad*180/np.pi
        return (anglerad)

    def build_child_struct(self, child_geometry, child_lattice, atom_types):
        ''' Creates a Structure() for the child to be added to the collection'''
        child_A, child_B, child_C = child_lattice
        temp_vol = np.dot(np.cross(child_A, child_B), child_C)
        alpha = self.angle(child_B, child_C)*180./np.pi
        beta = self.angle(child_C, child_A)*180./np.pi
        gamma = self.angle(child_A, child_B)*180./np.pi
        struct = Structure()
        for i in range(len(child_geometry)):
            struct.build_geo_by_atom(float(child_geometry[i][0]), float(child_geometry[i][1]),
                                     float(child_geometry[i][2]), atom_types[i])
        struct.set_property('lattice_vector_a', child_A)
        struct.set_property('lattice_vector_b', child_B)
        struct.set_property('lattice_vector_c', child_C)
        struct.set_property('a', np.linalg.norm(child_A))
        struct.set_property('b', np.linalg.norm(child_B))
        struct.set_property('c', np.linalg.norm(child_C))
        struct.set_property('cell_vol', temp_vol)
        struct.set_property('mutation_type', "Permute_rot")
        struct.set_property('alpha',alpha)
        struct.set_property('beta', beta)
        struct.set_property('gamma', gamma)
        return struct


class PermutationReflectionMutation(object):
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
        self.mol_list = self.input_struct.get_molecules(self.num_mols)
        self.atom_types = self.input_struct.get_atom_types()

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.permutation_mutation()

    def permutation_mutation(self):
        '''Permutes the location of each molecule eg. mol 1234 -> 4123'''
        mol_list_COM = get_COM_mol_list(self.mol_list)
        permuted_geometry = self.permute_reflect_cartesian(self.mol_list, mol_list_COM)
        permutated_struct = self.build_child_struct(permuted_geometry, self.atom_types)
        return permutated_struct

    def permute_reflect_cartesian(self, mol_list, mol_list_COM):
        ''' 
        Randomly permutes the location of each molecule
        and the performsa reflection
        '''

        self.output("-- Permuting COM of molecules")
        permuted_geometry = []
        perm = [i for i in range(len(mol_list))]
        np.random.shuffle(perm)
        self.output("-- New ordering %s" %(perm))
        reflections = [[-1, 1, 1],
                       [1, -1, 1],
                       [1, 1, -1]]
        choice = np.random.choice(range(len(reflections)))
        self.output("-- Reflection %s" %(reflections[choice]))
        for i in range(len(mol_list)):
            mol = mol_list[i]
            for atom in mol:
                xyz = [atom[0], atom[1], atom[2]]
                geo_COM = xyz - mol_list_COM[i] 
                geo_ref = np.multiply(reflections[choice], np.array(geo_COM))
                geo_final = geo_ref + np.array(mol_list_COM[perm[i]])
                permuted_geometry.append(np.array(geo_final))

        return permuted_geometry

    def build_child_struct(self, swapped_geo, atom_types):
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
        struct.set_property('mutation_type', 'Permute_ref')
        return struct

class RandomStrainMutation(object):
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
        self.mol_list = self.input_struct.get_molecules(self.num_mols)
        self.atom_types = self.input_struct.get_atom_types()

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.strain_lat_and_mols()

    def strain_lat_and_mols(self):
        # Strain lattice
        lat_mat = set_lat_mat(self.A, self.B, self.C)
        lat_mat_f = np.linalg.inv(lat_mat)
        strain_A, strain_B, strain_C = self.rand_strain(lat_mat)
        # Move COM of molecules accordingly
        mol_list_COM = get_COM_mol_list(self.mol_list)
        mol_list_COM_f = get_COM_mol_list_f(lat_mat_f, np.array(mol_list_COM))
        strain_lat_mat = set_lat_mat(strain_A, strain_B, strain_C)
        strained_geometry = strain_geometry(strain_lat_mat, self.mol_list, mol_list_COM, mol_list_COM_f)
        strained_struct = (self.create_strained_struct(strain_A, strain_B, strain_C, 
                                                          strained_geometry, self.atom_types))
        return strained_struct

    def rand_strain(self, lat_mat):
        strain_list = np.random.normal(scale=self.st_dev, size=6)
        strain_mat = get_strain_mat(strain_list)
        self.output("-- Strain_mat" + str(strain_mat))
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
        struct.set_property('beta', angle(lat_C, lat_A))
        struct.set_property('gamma', angle(lat_A, lat_B))
        struct.set_property('mutation_type', 'Rand_strain')
        return struct

class RandomSymmetryStrainMutation(object):
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
        self.mol_list = self.input_struct.get_molecules(self.num_mols)
        self.atom_types = self.input_struct.get_atom_types()

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.strain_lat_and_mols()

    def strain_lat_and_mols(self):
        # Strain lattice
        lat_mat = set_lat_mat(self.A, self.B, self.C)
        lat_mat_f = np.linalg.inv(lat_mat)
        strain_A, strain_B, strain_C = self.rand_sym_strain(lat_mat)

        # Move COM of molecules accordingly
        mol_list_COM = get_COM_mol_list(self.mol_list)
        mol_list_COM_f = get_COM_mol_list_f(lat_mat_f, np.array(mol_list_COM))
        strain_lat_mat = set_lat_mat(strain_A, strain_B, strain_C)
        strained_geometry = strain_geometry(strain_lat_mat, self.mol_list, mol_list_COM, mol_list_COM_f)
        strained_struct = (self.create_strained_struct(strain_A, strain_B, strain_C, 
                                                          strained_geometry, self.atom_types))
        return strained_struct

    def rand_sym_strain(self, lat_mat):
        strain_param = np.random.normal(scale=self.st_dev, size=1)
        strain_list = strain_param*get_rand_sym_strain(lat_mat)
        strain_mat = get_strain_mat(strain_list)
        self.output("-- Strain parameter: " + str(strain_param))
        self.output("-- Strain_matrix: \n" + str(strain_mat))
        strain = np.asarray(np.mat(lat_mat) + np.mat(strain_mat)*np.mat(lat_mat))
        self.output("-- Strained_lattice: \n" + str(strain))
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
        struct.set_property('beta', angle(lat_C, lat_A))
        struct.set_property('gamma', angle(lat_A, lat_B))
        struct.set_property('mutation_type', 'Sym_strain')
        return struct

class RandomVolumeStrainMutation(object):
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
        self.a = input_struct.get_property('a')
        self.b = input_struct.get_property('b')
        self.c = input_struct.get_property('c')
        self.alpha = np.deg2rad(input_struct.get_property('alpha'))
        self.beta = np.deg2rad(input_struct.get_property('beta'))
        self.gamma = np.deg2rad(input_struct.get_property('gamma'))
        self.st_dev = self.ui.get_eval('mutation', 'stand_dev_strain')
        self.cross_type = self.input_struct.get_property('crossover_type')
        self.cell_vol = self.input_struct.get_unit_cell_volume()
        self.mol_list = self.input_struct.get_molecules(self.num_mols)
        self.atom_types = self.input_struct.get_atom_types()

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.strain_lat_and_mols()

    def strain_lat_and_mols(self):
        # Strain lattice and preserve original unit cell volume
        lat_mat = set_lat_mat(self.A, self.B, self.C)
        lat_mat_f = np.linalg.inv(lat_mat)
        selection = np.random.choice(['a','b','c'])
        change = np.random.normal(scale=self.st_dev)
        scale_fac = (np.sqrt(1 + 2*np.cos(self.alpha)*np.cos(self.beta)*np.cos(self.gamma)
                               -np.cos(self.alpha)**2-np.cos(self.beta)**2-np.cos(self.gamma)**2))
        rand = np.random.random()
        if selection == 'a':
            new_a = self.a * (1 - change)
            diff_a = new_a - self.a
            diff_b = rand * diff_a
            diff_c = (1 - rand) * diff_a
            new_b = self.b - diff_b
            new_c = self.cell_vol/((new_a * new_b) * scale_fac)
        if selection == 'b':
            new_b = self.b * (1 - change)
            diff_b = new_b - self.b
            diff_c = rand * diff_b
            diff_a = (1 - rand) * diff_b
            new_c = self.c - diff_c
            new_a = self.cell_vol/((new_b * new_c) * scale_fac)
        if selection == 'c':
            new_c = self.c * (1 - change)
            diff_c = new_c - self.c
            diff_a = rand * diff_c
            diff_b = (1 - rand) * diff_c
            new_a = self.a - diff_a
            new_b = self.cell_vol/((new_c * new_a) * scale_fac)
        self.output("-- Changing lattice vector %s by %.2f percent" % (selection, change*100))
        self.output("-- Updating other lattice vectors to preserve %.2f A^3 volume" % (self.cell_vol))
        strain_A = [new_a, 0, 0]
        strain_B = [np.cos(self.gamma) * new_b, np.sin(self.gamma) * new_b, 0]
        cx = np.cos(self.beta)* new_c
        cy = (new_b * new_c * np.cos(self.alpha) - strain_B[0] * cx)/strain_B[1]
        cz = (new_c**2 - cx**2 - cy**2)**0.5
        strain_C = [cx, cy, cz]
        
        # Move COM of molecules accordingly
        mol_list_COM = get_COM_mol_list(self.mol_list)
        mol_list_COM_f = get_COM_mol_list_f(lat_mat_f, np.array(mol_list_COM))
        strain_lat_mat = set_lat_mat(strain_A, strain_B, strain_C)
        strained_geometry = strain_geometry(strain_lat_mat, self.mol_list, mol_list_COM, mol_list_COM_f)
        strained_struct = (self.build_child_struct(strain_A, strain_B, strain_C,
                                                          strained_geometry, self.atom_types))
        return strained_struct

    def rand_sym_strain(self, lat_mat):
        strain_param = np.random.normal(scale=self.st_dev, size=1)
        strain_list = strain_param*get_rand_sym_strain(lat_mat)
        strain_mat = get_strain_mat(strain_list)
        self.output("-- Strain parameter: " + str(strain_param))
        self.output("-- Strain_matrix: \n" + str(strain_mat))
        strain = np.asarray(np.mat(lat_mat) + np.mat(strain_mat)*np.mat(lat_mat))
        self.output("-- Strained_lattice: \n" + str(strain))
        return strain[0], strain[1], strain[2]

    def build_child_struct(self, lat_A, lat_B, lat_C, strained_geo, atom_types):
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
        struct.set_property('beta', angle(lat_C, lat_A))
        struct.set_property('gamma', angle(lat_A, lat_B))
        struct.set_property('mutation_type', 'Vol_strain')
        return struct

class RandomCellAngleMutation(object):
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
        self.a = input_struct.get_property('a')
        self.b = input_struct.get_property('b')
        self.c = input_struct.get_property('c')
        self.alpha = input_struct.get_property('alpha')
        self.beta = input_struct.get_property('beta')
        self.gamma = input_struct.get_property('gamma')
        self.st_dev = self.ui.get_eval('mutation', 'stand_dev_cell_angle')
        self.cross_type = self.input_struct.get_property('crossover_type')
        self.cell_vol = self.input_struct.get_unit_cell_volume()
        self.mol_list = self.input_struct.get_molecules(self.num_mols)
        self.atom_types = self.input_struct.get_atom_types()

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.strain_lat_and_mols()

    def strain_lat_and_mols(self):
        # Strain lattice by changing a single unit cell angle
        lat_mat = set_lat_mat(self.A, self.B, self.C)
        lat_mat_f = np.linalg.inv(lat_mat)
        count = 0
        while True:
            selection = np.random.choice(['alpha','beta','gamma'])
            degrees = np.random.normal(scale=self.st_dev, size=1)
            if selection == 'alpha':
                self.alpha = self.alpha + degrees
            elif selection == 'beta':
                self.beta = self.beta + degrees
            elif selection == 'gamma':
                self.gamma = self.gamma + degrees

            self.output("-- Changing unit cell angle %s by %.2f degrees" % (selection, degrees))
            self.output("-- alpha=%.2f beta=%.2f gamma=%.2f " % (self.alpha, self.beta, self.gamma))

            alpha = np.deg2rad(self.alpha)
            beta = np.deg2rad(self.beta)
            gamma = np.deg2rad(self.gamma)
            strain_A = [self.a, 0, 0]
            strain_B = [float(np.cos(gamma) * self.b), float(np.sin(gamma)) * self.b, 0]
            cx = float(np.cos(beta))* self.c
            cy = (self.b * self.c * float(np.cos(alpha)) - strain_B[0] * cx)/strain_B[1]
            cz = self.c**2 - cx**2 - cy**2
            if cz > 0:
                cz = (self.c**2 - cx**2 - cy**2)**0.5
                strain_C = [cx, cy, cz]
                break
        # Move COM of molecules accordingly
        mol_list_COM = get_COM_mol_list(self.mol_list)
        mol_list_COM_f = get_COM_mol_list_f(lat_mat_f, np.array(mol_list_COM))
        strain_lat_mat = set_lat_mat(strain_A, strain_B, strain_C)
        strained_geometry = strain_geometry(strain_lat_mat, self.mol_list, mol_list_COM, mol_list_COM_f)
        strained_struct = (self.build_child_struct(strain_A, strain_B, strain_C,
                                                          strained_geometry, self.atom_types))
        return strained_struct

    def rand_sym_strain(self, lat_mat):
        strain_param = np.random.normal(scale=self.st_dev, size=1)
        strain_list = strain_param*get_rand_sym_strain(lat_mat)
        strain_mat = get_strain_mat(strain_list)
        self.output("-- Strain parameter: " + str(strain_param))
        self.output("-- Strain_matrix: \n" + str(strain_mat))
        strain = np.asarray(np.mat(lat_mat) + np.mat(strain_mat)*np.mat(lat_mat))
        self.output("-- Strained_lattice: \n" + str(strain))
        return strain[0], strain[1], strain[2]

    def build_child_struct(self, lat_A, lat_B, lat_C, strained_geo, atom_types):
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
        struct.set_property('beta', angle(lat_C, lat_A))
        struct.set_property('gamma', angle(lat_A, lat_B))
        struct.set_property('mutation_type', 'Angle_strain')
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
    strain_dict={                 \
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
    rand_choice = str(np.random.randint(5,29))
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

def mat2euler(M, cy_thresh=None):
        M = np.asarray(M)
        if cy_thresh is None:
            try:
                cy_thresh = np.finfo(M.dtype).eps * 4
            except ValueError:
                cy_thresh = _FLOAT_EPS_4
        r11, r12, r13, r21, r22, r23, r31, r32, r33 = M.flat
        cy = math.sqrt(r33 * r33 + r23 * r23)
        if cy > cy_thresh: # cos(y) not close to zero, standard form
            z = math.atan2(-r12,  r11) # atan2(cos(y)*sin(z), cos(y)*cos(z))
            y = math.atan2(r13,  cy) # atan2(sin(y), cy)
            x = math.atan2(-r23, r33) # atan2(cos(y)*sin(x), cos(x)*cos(y))
        else: # cos(y) (close to) zero, so x -> 0.0 (see above)
            # so r21 -> sin(z), r22 -> cos(z) and
            z = math.atan2(r21,  r22)
            y = math.atan2(r13,  cy) # atan2(sin(y), cy)
            x = 0.0
        return z, y, x

def euler2mat(z=0, y=0, x=0):
        Ms = []
        if z:
            cosz = math.cos(z)
            sinz = math.sin(z)
            Ms.append(np.array(
                 [[cosz, -sinz, 0],
                 [sinz, cosz, 0],
                 [0, 0, 1]]))
        if y:
            cosy = math.cos(y)
            siny = math.sin(y)
            Ms.append(np.array(
                [[cosy, 0, siny],
                 [0, 1, 0],
                 [-siny, 0, cosy]]))
        if x:
            cosx = math.cos(x)
            sinx = math.sin(x)
            Ms.append(np.array(
                [[1, 0, 0],
                 [0, cosx, -sinx],
                 [0, sinx, cosx]]))
        if Ms:
            return reduce(np.dot, Ms[::-1])
        return np.eye(3)







if __name__ == "__main__":
    main()
