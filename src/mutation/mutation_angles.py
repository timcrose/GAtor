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
from pymatgen import Molecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer as pga

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
    vol = struct.get_unit_cell_volume()
    output.local_message("Input Structure's Volume: %s" % (vol), replica)
    mutate_obj = select_mutator(input_struct, num_mols, replica)
    mutated_struct = mutate_obj.mutate()
    return mutated_struct


def select_mutator(input_struct, num_mols, replica):
    '''
    In this mutation implementation, there are several classes, each performing a 
    different mutation. This method selecting which mutation to empoly
    Expects: Structure, number of molecules per cell, replica name
    Returns: Mutation Class
    '''
    mutation_list = (["Frame_trans_mol", "Pair_trans_mol", 
                      "Frame_rot_mol", "Pair_rot_mol", "Rot_mol",
                      "Strain_rand_mols","Strain_sym_mols", 
                      "Swap_mol", "Permute_mol","Strain_vol"])
    mutation_list = ["Permute_mol"]
    try:
        mut_choice = np.random.choice(mutation_list)
    except:
        mut_choice = mutation_list[int(np.random.random()*len(mutation_list))]
    message = "-- Mutation Choice: %s" % (mut_choice)
    output.local_message(message, replica)

    if mut_choice == "Frame_trans_mol":
        mutator = RandomTranslationFrameMutation(input_struct, num_mols, replica)
    elif mut_choice == "Frame_rot_mol":
        mutator = RandomRotationFrameMutation(input_struct, num_mols, replica)
    elif mut_choice == "Rot_mol":
        mutator = RandomRotationMolMutation(input_struct, num_mols, replica)
    elif mut_choice == "Pair_trans_mol":
        mutator = PairTranslationMutation(input_struct, num_mols, replica)
    elif mut_choice == "Pair_rot_mol":
        mutator = PairRotationMolMutation(input_struct, num_mols, replica)
    elif mut_choice == "Swap_mol":
	    mutator = SwapMolMutation(input_struct, num_mols, replica)
    elif mut_choice =="Permute_mol":
        mutator = PermutationMutation(input_struct, num_mols, replica)
    elif mut_choice == "Strain_rand_mols":
        mutator = RandomStrainMutationMoveMols(input_struct, num_mols, replica)
    elif mut_choice == "Strain_sym_mols":
        mutator = RandomSymmetryStrainMutationMoveMols(input_struct, num_mols, replica)
    elif mut_choice =="Strain_vol":
	    mutator = RandomVolumeStrainMutationMoveMols(input_struct, num_mols, replica)
    return mutator

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
    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.random_translation()

    def random_translation(self):
        '''Calls for a random translation to each molecule's COM and returns a Structure'''
        temp_geo = self.geometry
        num_atom_per_mol = int(len(temp_geo)/self.num_mols)
        atom_types = [temp_geo[i][3] for i in range(len(temp_geo))]
        mol_list = [temp_geo[x:x+num_atom_per_mol] for x in range(0, len(temp_geo), num_atom_per_mol)]
        mol_index = 0
        for mol in mol_list:
            COM = self.get_COM_frac(mol, [self.A, self.B, self.C])
            centered_mol, types = self.get_centered_molecule(mol)
            rot, aligned_mol = self.get_orientation_info(centered_mol, mol_index, types)
            self.orientation_info.append([rot, COM, aligned_mol])
        lattice, child_coords = self.translate_in_frame(self.orientation_info)
        mutated_struct = self.create_child_struct(child_coords, lattice, atom_types)
        return mutated_struct
    
    def translate_in_frame(self, orientation_info):
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
                mol_coords.append(np.dot(rot, np.array(atom).reshape((3,1))).tolist())
            for coord in mol_coords:
                new_coords = [coord[0][0] + COM_xyz[0], coord[1][0] + COM_xyz[1], coord[2][0]+COM_xyz[2]]
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

    def get_orientation_info(self, mol, mol_index, types):
        ''' Computes the principal axes for each molecule and the corresponding
            rotation matrices.
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
        rot = np.column_stack((sort_eig[2][1], -sort_eig[1][1], sort_eig[0][1]))
        rot_trans = np.transpose(np.array(rot))

        #Aligned geometry
        centered_sites = []
        aligned_mol = []
        rand_disp = np.random.standard_normal(3) * self.st_dev
        if ( (mol_index % 2) == 0 ): # racemic 
            for site in centered_molp.sites:
                centered_sites.append(list(site._coords))
            for atom in centered_sites:
                new_atom = np.dot(rot_trans, np.array(atom))
                new_atom = new_atom + rand_disp
                aligned_mol.append(new_atom)
        else: 
            for site in centered_molp.sites:
                centered_sites.append(list(site._coords))
            for atom in centered_sites: #chiral
                new_atom = np.dot(rot_trans, np.array(atom))
                atom = new_atoms - rand_disp
                aligned_mol.append(new_atom)
        return rot, aligned_mol

    def angle(self,v1,v2):
        numdot = np.dot(v1,v2)
        anglerad = np.arccos(numdot/(leng(v1)*leng(v2)))
        angledeg = anglerad*180/np.pi
        return (anglerad)

    def create_child_struct(self, child_geometry, child_lattice, atom_types):
        ''' Creates a Structure() for the child to be added to the collection'''
        child_A, child_B, child_C = child_lattice
        temp_vol = np.dot(np.cross(child_A, child_B), child_C)
        alpha = self.angle(child_A, child_B)*180./np.pi
        beta = self.angle(child_B, child_C)*180./np.pi
        gamma = self.angle(child_C, child_A)*180./np.pi
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
        struct.set_property('mutation_type', "trans_mol_frame")
        struct.set_property('alpha',alpha)
        struct.set_property('beta', beta)
        struct.set_property('gamma', gamma)
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

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.random_rotation()

    def random_rotation(self):
        '''Calls for a random translation to each molecule's COM and returns a Structure'''
        temp_geo = self.geometry
        atom_types = [self.geometry[i][3] for i in range(len(self.geometry))]
        num_atom_per_mol = int(len(temp_geo)/self.num_mols)
        mol_list = [temp_geo[x:x+num_atom_per_mol] for x in range(0, len(temp_geo), num_atom_per_mol)]
        mol_index = 0
        for mol in mol_list:
            COM = self.get_COM_frac(mol, [self.A, self.B, self.C])
            centered_mol, types = self.get_centered_molecule(mol)
            rot, aligned_mol = self.get_orientation_info(centered_mol, mol_index, types)
            self.orientation_info.append([rot, COM, aligned_mol])
        lattice, child_coords = self.rotate_in_frame(self.orientation_info)
        mutated_struct = self.create_child_struct(child_coords, lattice, atom_types)
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

    def rotate_in_frame(self, orientation_info):
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
                mol_coords.append(np.dot(rot, np.array(atom).reshape((3,1))).tolist())
            for coord in mol_coords:
                new_coords = [coord[0][0] + COM_xyz[0], coord[1][0] + COM_xyz[1], coord[2][0]+COM_xyz[2]]
                child_coordinates.append(new_coords)
        return lattice, np.array(child_coordinates)


    def get_orientation_info(self, mol, mol_index, types):
        ''' 
        Computes the principal axes for each molecule and the corresponding
        rotation matrices.
        Returns: rotations about z, y, x axis and molecule with principal axes 
        aligned. 
        '''
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
        rot = np.column_stack((sort_eig[2][1], -sort_eig[1][1], sort_eig[0][1]))
        rot_trans = np.transpose(np.array(rot))

        #Apply random rotation  
        centered_sites = []
        aligned_mol = []
        rand_angs = [0, 5, 10, 
                    15, 30, 45, 
                    60, 75, 90, 
                    105, 130, 145, 
                    160, 170, 175, 
                    180]

        u = random.choice([-1,1])
        t = random.choice([-1,1])
        v = random.choice([-1,1])
        theta = u*random.choice(rand_angs)
        psi = t*random.choice(rand_angs)
        phi = v*random.choice(rand_angs)
        rot_mut = rotation_matrix(theta, psi, phi)
     #   if ( (mol_index % 2) == 0 ): # racemic 
     #       for site in centered_molp.sites:
     #           centered_sites.append(list(site._coords))
     #       for atom in centered_sites:
     #           new_atom = np.dot(rot_trans, np.array(atom))
     #           new_atom = np.dot(rot_mut, new_atom)
     #           aligned_mol.append(new_atom)
     #   else:
     #       for site in centered_molp.sites:
     #           centered_sites.append(list(site._coords))
     #       for atom in centered_sites: #chiral
     #           new_atom = np.dot(rot_trans, np.array(atom))
     #           new_atom = np.dot(np.linalg.inv(rot_mut), new_atom)
     #           aligned_mol.append(new_atom)

        for site in centered_molp.sites:
            centered_sites.append(list(site._coords))
        for atom in centered_sites: #chiral
            new_atom = np.dot(rot_trans, np.array(atom))
            new_atom = np.dot(rot_mut, new_atom)
            aligned_mol.append(new_atom)
        return rot, aligned_mol

    def angle(self,v1,v2):
        numdot = np.dot(v1,v2)
        anglerad = np.arccos(numdot/(leng(v1) * leng(v2)))
        angledeg = anglerad*180/np.pi
        return (anglerad)

    def create_child_struct(self, child_geometry, child_lattice, atom_types):
        ''' Creates a Structure() for the child to be added to the collection'''
        child_A, child_B, child_C = child_lattice
        temp_vol = np.dot(np.cross(child_A, child_B), child_C)
        alpha = self.angle(child_A, child_B)*180./np.pi
        beta = self.angle(child_B, child_C)*180./np.pi
        gamma = self.angle(child_C, child_A)*180./np.pi
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
        struct.set_property('mutation_type', "rot_mol_frame")
        struct.set_property('alpha',alpha)
        struct.set_property('beta', beta)
        struct.set_property('gamma', gamma)
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
        ''' Randomly rotates all molecules COM by some random angle'''
        rot_geometry = []
        i = 0
        rand_vec = np.random.normal(scale=st_dev,size=3)
        rand_vec = random.choice([rand_vec, -rand_vec])
        theta= (np.pi/180)*rand_vec[0]
        psi = (np.pi/180)*rand_vec[1]
        phi= (np.pi/180)*rand_vec[2] 
        for mol in mol_list:
            rot_mat = rotation_matrix(theta, psi, phi)
            mol_COM = np.array([mol_list_COM[i][0], mol_list_COM[i][1],mol_list_COM[i][2]])
            for atom in mol:
                atom_vec = np.array([atom[0], atom[1], atom[2]]) - mol_COM
                rot_geometry.append(np.dot(np.asarray(rot_mat),atom_vec) + mol_COM)
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

class PairTranslationMutation(object):
    ''' 
    Gives a translation to the COM of pairs of molecules in the unit cell.
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
        return self.pair_translation()

    def pair_translation(self):
        '''Calls for a random rotation to each molecule's COM and returns a Structure'''
        temp_geo = self.geometry
        num_atom_per_mol = int(len(temp_geo)/self.num_mols)
        mol_list = [temp_geo[x:x+num_atom_per_mol] for x in range(0, len(temp_geo), num_atom_per_mol)]
        atom_type_list = [temp_geo[i][3] for i in range(len(temp_geo))]
        mol_list_COM = get_COM_mol_list(mol_list)
        translated_geometry = self.translate_molecules_pair(mol_list, mol_list_COM, self.st_dev)
        translated_struct = self.create_translated_struct(translated_geometry, atom_type_list)
        return translated_struct

    def translate_molecules_pair(self, mol_list, mol_list_COM, st_dev):
        ''' Randomly rotates each molecule within gaussian dist'''
        #Generate Translation 
        trans_geometry = []
        disp = np.random.normal(scale=st_dev)
        transs = ([disp, 0, 0],
                [0, disp, 0],
                [0, 0, disp],
                [disp/np.sqrt(2), disp/np.sqrt(2), 0],
                [0, disp/np.sqrt(2), disp/np.sqrt(2)],
                [disp/np.sqrt(2), 0 , disp/np.sqrt(2)],
                [disp/np.sqrt(3), disp/np.sqrt(3), disp/np.sqrt(3)])
        trans = np.array(random.choice(transs))
        if np.random.random() < 0.5:
            direction_lat = random.choice([self.A, self.B, self.C])
            self.output("-- Translation along lattice vector %s " %(direction_lat))
            direction = direction_lat/np.linalg.norm(direction_lat)  
            trans = np.multiply(direction, trans)
        self.output("-- Translation is %s" % (trans))

        # Select molecule pairs to translate
        mol_in = 0
        mol_info = []
        mol_combos = list(itertools.combinations(range(self.num_mols), 2))
        for mol in mol_list:
            mol_info.append([mol, mol_list_COM[mol_in]])
            mol_in +=1
        mol_combo = random.choice(mol_combos)
        self.output("-- On molecules: ")
        self.output(mol_combo)

        # Translate those pairs
        i = 0
        for mol in mol_list:
            if mol_combo[0] == i:
                for atom in mol:
                    coord = np.array([atom[0], atom[1], atom[2]])
                    trans_geometry.append(coord + trans)
            elif mol_combo[1] == i:
                if np.random.random() < 0.5:
                    trans = -trans
                    self.output("Opposite translation b/t molecules")
                for atom in mol:
                    coord = np.array([atom[0], atom[1], atom[2]]) 
                    trans_geometry.append(coord + trans)
            else:
                for atom in mol:
                    coord = np.array([atom[0], atom[1], atom[2]])
                    trans_geometry.append(np.array(coord))
            i+=1
        return trans_geometry

    def create_translated_struct(self, trans_geo, atom_types):
        ''' Creates Structure from mutated geometry'''
        struct = Structure()
        for i in range(len(trans_geo)):
            struct.build_geo_by_atom(float(trans_geo[i][0]), float(trans_geo[i][1]),
                                     float(trans_geo[i][2]), atom_types[i])
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
        struct.set_property('mutation_type', 'pair_trans_mol')
        return struct
class PairRotationMolMutation(object):
    ''' 
    Gives a rotation to the COM of the pairs of molecules in the unit cell.
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
        return self.pair_rotation()

    def pair_rotation(self):
        '''Calls for a random rotation to each molecule's COM and returns a Structure'''
        temp_geo = self.geometry
        num_atom_per_mol = int(len(temp_geo)/self.num_mols)
        mol_list = [temp_geo[x:x+num_atom_per_mol] for x in range(0, len(temp_geo), num_atom_per_mol)]
        atom_type_list = [temp_geo[i][3] for i in range(len(temp_geo))]
        mol_list_COM = get_COM_mol_list(mol_list)
        rotated_geometry = self.rotate_molecules_pair(mol_list, mol_list_COM, self.st_dev)
        rotated_struct = self.create_rotated_struct(rotated_geometry, atom_type_list)
        return rotated_struct

    def rotate_molecules_pair(self, mol_list, mol_list_COM, st_dev):
        ''' Randomly rotates each molecule within gaussian dist'''
        rot_geometry = []
        ang = random.choice([5, 10, 15, 
                            30, 45, 60, 
                            90, 120, 150, 
                            170, 175, 180])
        ang = random.choice([ang, -ang])
        rots = ([ang, 0, 0], 
                [0, ang, 0], 
                [0, 0, ang],
                [ang, ang, 0], 
                [0, ang, ang], 
                [ang, 0 , ang])
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
                    rot_mat = np.asarray(rotation_matrix(rot[0], rot[1], rot[2]))
                    rot_geometry.append(np.dot(rot_mat, atom_vec) + mol_COM)
            elif mol_combo[1] == i:
                mol_COM = np.array([
                                   mol_list_COM[mol_combo[1]][0], 
                                   mol_list_COM[mol_combo[1]][1],
                                   mol_list_COM[mol_combo[1]][2]])
                for atom in mol:
                    atom_vec = np.array([atom[0], atom[1], atom[2]]) - mol_COM
                    rot_mat = np.asarray(rotation_matrix(rot[0], rot[1], rot[2]))
                    rot_geometry.append(np.dot(rot_mat,atom_vec) + mol_COM)
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
        struct.set_property('mutation_type', 'pair_rot_mol')
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
        ''' Randomly swaps pairs of molecules'''
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
                    xyz = [atom[0], atom[1], atom[2]]
                    atom_vec = np.array(xyz) - mol_choice1[1] + mol_choice2[1]
                    swap_geometry.append(atom_vec)
            elif mol[0][0] == mol_choice2[0][0][0]:
                for atom in mol:
                    xyz = [atom[0], atom[1], atom[2]]
                    atom_vec = np.array(xyz) - mol_choice2[1] + mol_choice1[1]
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

class PermutationMutation(object):
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
        return self.permutation_mutation()

    def permutation_mutation(self):
        '''Permutes the location of each molecule eg. mol 1234 -> 4123'''
        temp_geo = self.geometry
        num_atom_per_mol = int(len(temp_geo)/self.num_mols)
        mol_list = [temp_geo[x:x+num_atom_per_mol] for x in range(0, len(temp_geo), num_atom_per_mol)]
        mol_list_COM = get_COM_mol_list(mol_list)
        atom_type_list = [temp_geo[i][3] for i in range(len(temp_geo))]
        permuted_geometry = self.permute_molecules(mol_list, mol_list_COM)
        permutated_struct = self.create_permutated_struct(permuted_geometry, atom_type_list)
        return permutated_struct

    def permute_molecules(self, mol_list, mol_list_COM):
        ''' Randomly permutes the location of each molecule'''

        self.output("Permuting COM of molecules")
        permuted_geometry = []
        permute = [i for i in range(1, len(mol_list)-1)]
        choice = random.choice(permute)    

        for i in range(len(mol_list)):
            mol = mol_list[i]
            for atom in mol:
                xyz = [atom[0], atom[1], atom[2]]
                geo = xyz - mol_list_COM[i] + mol_list_COM[choice]
                permuted_geometry.append(np.array(geo))
            choice +=1
            if choice % len(mol_list) == 0:
                choice = 0

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
        struct.set_property('mutation_type', 'permute_mol')
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
        struct.set_property('alpha', angle(lat_A, lat_B))
        struct.set_property('beta', angle(lat_B, lat_C))
        struct.set_property('gamma', angle(lat_C, lat_A))
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
        struct.set_property('alpha', angle(lat_A, lat_B))
        struct.set_property('beta', angle(lat_B, lat_C))
        struct.set_property('gamma', angle(lat_C, lat_A))
        struct.set_property('mutation_type', 'sym_strain_mol')
        return struct

class RandomVolumeStrainMutationMoveMols(object):
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

    def output(self, message): output.local_message(message, self.replica)

    def mutate(self):
        return self.strain_lat_and_mols()

    def strain_lat_and_mols(self):
        temp_geo = self.geometry
        num_atom_per_mol = int(len(temp_geo)/self.num_mols)
        atom_type_list = [temp_geo[i][3] for i in range(len(temp_geo))]
        lat_mat = set_lat_mat(self.A, self.B, self.C)
        lat_mat_f = np.linalg.inv(lat_mat)

    	selection = np.random.choice(['a','b','c'])
        change = np.random.normal(scale=0.15)
        self.output("Decimal change %s" % (change))
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
        self.output("Biggest change on lattice vector %s" % (selection))
        strain_A = [new_a, 0, 0]
        strain_B = [np.cos(self.gamma) * new_b, np.sin(self.gamma) * new_b, 0]
        cx = np.cos(self.beta)* new_c
        cy = (new_b * new_c * np.cos(self.alpha) - strain_B[0] * cx)/strain_B[1]
        cz = (new_c**2 - cx**2 - cy**2)**0.5
        strain_C = [cx, cy, cz]
        
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
        struct.set_property('alpha', angle(lat_A, lat_B))
        struct.set_property('beta', angle(lat_B, lat_C))
        struct.set_property('gamma', angle(lat_C, lat_A))
        struct.set_property('mutation_type', 'vol_strain')
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
