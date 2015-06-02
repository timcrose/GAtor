'''
Created on May 1, 2015

@authors: farren
'''
from __future__ import division

from copy import deepcopy
import numpy as np
import random
import time

from core import user_input
from structures.structure import Structure

def main(struct, r_stoic, replica):
    '''
    Every mutation module must have a method named "main" that will take
    as arguments one Structure and a target stoichiometry (in the form of
    a StoicDict). It must return a single Structure with defined geometry
    or return False if the mutation fails and a new Structure is needed.
    '''
    # deepcopy makes copy of data instead of altering the obect being referenced
    input_structure = deepcopy(struct)
    replica_stoic = deepcopy(r_stoic)
    
    # decides what type of mutation to execute
    mutate_object = select_mutator(input_structure, replica_stoic, replica)
    return mutate_object.mutate()

def select_mutator(input_structure, replica_stoic, replica):
    '''
    In this mutation implementation, there are several classes, each performing a 
    different mutation. This method is responsible for reading the preferences set
    by the user and selecting which mutation to empoly, or no mutation at all.
    Expects: Structure, target_stoic
    Returns: Class
    '''
    #Make sure if single parent is selected at crossover that it is forced to get a mutation
    if input_structure.get_property('cross_type') == [1,1]:
	mutation_list = ["Trans_mol", "Rot_mol", "Strain_rand_mols", "Strain_rand", "Strain_sym"]
    elif input_structure.get_property('cross_type') == [2,2]:
	mutation_list = ["Trans_mol", "Rot_mol", "Strain_rand_mols", "Strain_rand", "Strain_sym"]
    else:
	mutation_list = ["None", "Trans_mol", "Rot_mol", "Strain_rand_mols", "Strain_rand", "Strain_sym"]

    mut_choice = np.random.choice(mutation_list)
    mutator = object
    print "Mutation_Choice: ", mut_choice
   
    if mut_choice == "None":
    	mutator = NoMutation(input_structure, replica_stoic, replica)
    elif mut_choice == "Trans_mol":
	mutator = RandomTranslationMutation(input_structure, replica_stoic, replica)
    elif mut_choice == "Rot_mol":
	mutator = RandomRotationMolMutation(input_structure, replica_stoic, replica)
    elif mut_choice == "Strain_rand_mols":
	mutator = RandomStrainMutationMoveMols(input_structure, replica_stoic, replica)
    elif mut_choice == "Strain_rand":
	mutator = RandomStrainMutation(input_structure, replica_stoic, replica)
    elif mut_choice == "Strain_sym":
	mutator = RandomSymmetryStrainMutationMoveMols(input_structure, replica_stoic, replica)
    return mutator	

###########################################################################################    
class NoMutation(object):
    '''
     This class leaves the geometry un-mutilated
    '''
    def __init__(self, input_structure, target_stoic, replica):
        self.geometry = deepcopy(input_structure.get_geometry())
	self.input_structure = input_structure	
    def mutate(self):
        new_struct = Structure()
        new_struct.build_geo_whole(self.geometry)
        new_struct.set_property('lattice_vector_a', self.input_structure.get_property('lattice_vector_a'))
        new_struct.set_property('lattice_vector_b', self.input_structure.get_property('lattice_vector_b'))
        new_struct.set_property('lattice_vector_c', self.input_structure.get_property('lattice_vector_c'))
        new_struct.set_property('cell_vol', self.input_structure.get_property('cell_vol'))
        new_struct.set_property('crossover_type', self.input_structure.get_property('crossover_type'))
        new_struct.set_property('alpha',self.input_structure.get_property('alpha'))
        new_struct.set_property('beta', self.input_structure.get_property('beta'))
        new_struct.set_property('gamma', self.input_structure.get_property('gamma'))
        return new_struct

####################################################################################################
class RandomTranslationMutation(object):
    '''
     This mutation gives a random translation to the COM of one or more of the molecules in the unit cell.
    '''
    def __init__(self, input_structure, target_stoic, replica):
        self.geometry = deepcopy(input_structure.get_geometry())
        self.ui = user_input.get_config()
        self.input_structure = input_structure
	self.num_mols = self.ui.get_eval('unit_cell_settings', 'num_molecules')
    	self.st_dev = self.ui.get_eval('mutation', 'stand_dev_trans')
    
    def mutate(self):
        return self.random_translation()  

    def random_translation(self):
        ''' 
        make a random translation within reasonable bounds and returns Structure
        '''
        temp_geo = self.center_geometry(deepcopy(self.geometry))  # center the geometry
	atom_num_per_mol = int(len(temp_geo)/self.num_mols)
	mol_list = [temp_geo[x:x+atom_num_per_mol] for x in range(0, len(temp_geo), atom_num_per_mol)]
 
	#Displace molecules (CHANGE FOR >2 mols)
	displace_mol1 = self.displace_center_geometry(mol_list[0], self.st_dev)
	displace_mol2 = self.displace_center_geometry(mol_list[1], self.st_dev)
	return_geo =  np.concatenate((displace_mol1, displace_mol2), axis=0)

	#Set new structure
        new_struct = Structure()
        new_struct.build_geo_whole(return_geo)
        new_struct.set_property('lattice_vector_a', self.input_structure.get_property('lattice_vector_a'))
        new_struct.set_property('lattice_vector_b', self.input_structure.get_property('lattice_vector_b'))
        new_struct.set_property('lattice_vector_c', self.input_structure.get_property('lattice_vector_c'))
        new_struct.set_property('cell_vol', self.input_structure.get_property('cell_vol'))
        new_struct.set_property('crossover_type', self.input_structure.get_property('crossover_type'))
        new_struct.set_property('alpha',self.input_structure.get_property('alpha'))
        new_struct.set_property('beta', self.input_structure.get_property('beta'))
        new_struct.set_property('gamma', self.input_structure.get_property('gamma'))

	return new_struct

    def displace_center_geometry(self, geometry, st_dev):
        '''
        randomly displaces the COM of a molecule within gaussian dist
        '''
   	rand_disp = np.random.standard_normal(3) * st_dev
        print rand_disp	
        for atom in geometry:
	    atom[0] = atom[0] - rand_disp[0]
	    atom[1] = atom[1] - rand_disp[1]
	    atom[2] = atom[2] - rand_disp[2]		
        return geometry
    
    def center_geometry(self, geometry):
        '''
        the centers the origin in relation to the max and min of each axis
        '''
        for i in range(2):  # x, y, and z
            coordinate_sum = 0
            counter = 0
            for atom in geometry:
                coordinate_sum += atom[i]
                counter += 1
            # average axis value
            average = coordinate_sum / counter 
            for atom in geometry:
                # shift all towards center
                atom[i] = atom[i] - average 
        return geometry

############################################################################################
class RandomRotationMolMutation(object):
    '''
     This mutation gives a random rotation to the COM of one or more of the molecules in the unit cell.
    '''
    def __init__(self, input_structure, target_stoic, replica):
        self.geometry = deepcopy(input_structure.get_geometry())
        self.ui = user_input.get_config()
        self.input_structure = input_structure
        self.num_mols = self.ui.get_eval('unit_cell_settings', 'num_molecules')
        self.st_dev = self.ui.get_eval('mutation', 'stand_dev_rot')

    def mutate(self):
        return self.random_rotation()

    def random_rotation(self):
        '''
        make a random translation within reasonable bounds and returns Structure
        '''
	#Center the geometry and makes lists of the molecules and atom types
        temp_geo = self.center_geometry(deepcopy(self.geometry))  # center the geometry
        atom_num_per_mol = int(len(temp_geo)/self.num_mols)
        mol_list = [temp_geo[x:x+atom_num_per_mol] for x in range(0, len(temp_geo), atom_num_per_mol)]
	atoms = [i for i in range(len(temp_geo))]
	for i in range(len(temp_geo)):
        	atoms[i] = temp_geo[i][3] 

        #Rotate molecules (CHANGE FOR >2 mols)
        rot_mol1 = self.rot_center_geometry(mol_list[0], self.st_dev)
        rot_mol2 = self.rot_center_geometry(mol_list[1], self.st_dev)
        return_geo =  self.center_geometry(np.concatenate((rot_mol1, rot_mol2), axis=0))

        #Set new structure
        new_struct = Structure()
        for i in range(len(return_geo)):
                 new_struct.build_geo_by_atom(float(return_geo[i][0]),
                                                 float(return_geo[i][1]),
                                                 float(return_geo[i][2]),
                                                 atoms[i])
        new_struct.set_property('lattice_vector_a', self.input_structure.get_property('lattice_vector_a'))
        new_struct.set_property('lattice_vector_b', self.input_structure.get_property('lattice_vector_b'))
        new_struct.set_property('lattice_vector_c', self.input_structure.get_property('lattice_vector_c'))
        new_struct.set_property('cell_vol', self.input_structure.get_property('cell_vol'))
        new_struct.set_property('crossover_type', self.input_structure.get_property('crossover_type'))
        new_struct.set_property('alpha',self.input_structure.get_property('alpha'))
        new_struct.set_property('beta', self.input_structure.get_property('beta'))
        new_struct.set_property('gamma', self.input_structure.get_property('gamma'))

        return new_struct

    def rot_center_geometry(self, geometry, st_dev):
        '''
        randomly rotates each molecule
        '''	
	#Initialize and switch to COM frame of mol
	xyz = [i for i in range(len(geometry))]
	rot = [i for i in range(len(geometry))] 
	
	#Rotate the mol by randome angles with standard dev
	rand_angles = np.random.standard_normal(3) * st_dev
	for i in range(len(geometry)):
		xyz[i]= [geometry[i][0],geometry[i][1], geometry[i][2]]
		rot[i] = np.dot(np.asarray(self.rotation_matrix(rand_angles)), np.asarray(xyz[i]))

	result_struct = np.asarray(rot)	
	return result_struct

    def rotation_matrix(self, angles):
	theta= (np.pi/180)*angles[0]
        psi = (np.pi/180)*angles[1]     
        phi= (np.pi/180)*angles[2]
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

    def center_geometry(self, geometry):
        '''
        the centers the origin in relation to the max and min of each axis
        '''
        for i in range(2):  # x, y, and z
            coordinate_sum = 0
            counter = 0
            for atom in geometry:
                coordinate_sum += atom[i]
                counter += 1
            # average axis value
            average = coordinate_sum / counter
            for atom in geometry:
                # shift all towards center
                atom[i] = atom[i] - average
        return geometry

##############################################################################################
	

##############################################################################################
class RandomStrainMutationMoveMols(object):
    '''
     This mutation gives a random strain to the lattice and moves the COM of the molecules
    '''
    def __init__(self, input_structure, target_stoic, replica):
        self.geometry = deepcopy(input_structure.get_geometry())
        self.A = np.asarray(deepcopy(input_structure.get_property('lattice_vector_a')))
        self.B = np.asarray(deepcopy(input_structure.get_property('lattice_vector_b')))
        self.C = np.asarray(deepcopy(input_structure.get_property('lattice_vector_c')))
	self.ui = user_input.get_config()
        self.input_structure = input_structure
	self.st_dev = self.ui.get_eval('mutation', 'stand_dev_strain')
	self.num_mols = self.ui.get_eval('unit_cell_settings', 'num_molecules')

    def mutate(self):
        return self.rstrain()

    def rstrain(self):
	A = self.A
	B = self.B
	C = self.C 	
        lat_mat = np.zeros((3,3))
        lat_mat[0][0] = A[0]; lat_mat[0][1] = B[0]; lat_mat[0][2] = C[0]
        lat_mat[1][0] = A[1]; lat_mat[1][1] = B[1]; lat_mat[1][2] = C[1]
        lat_mat[2][0] = A[2]; lat_mat[2][1] = B[2]; lat_mat[2][2] = C[2]

	#Fractional Transformation Matrix
        lat_mat_f = np.linalg.inv(lat_mat)

	#Strained Lattice Vecctors
	strain_A, strain_B, strain_C = self.rand_vstrain(lat_mat)

	#Lengths and angles of strained vectors
        a = self.leng(strain_A)
        b = self.leng(strain_B)
        c = self.leng(strain_C)
        alpha = self.angle(strain_B,strain_C)
        beta = self.angle(strain_A,strain_C)
        gamma = self.angle(strain_A,strain_B)

	#Realign axis to have a along x
        rad=float(np.pi/180)
        ax=a;   ay=0.0;   az=0.0
        bx=np.cos(gamma*rad)*b; by=np.sin(gamma*rad)*a; bz=0.0
        cx=c*np.cos(beta*rad);  cy=(b*c*np.cos(alpha*rad)-bx*cx)/by;    cz=np.sqrt(np.absolute(c**2-cx**2-cy**2))
        lata_out = np.zeros(3); latb_out = np.zeros(3); latc_out = np.zeros(3)
        lata_out[0] = ax;       lata_out[1] = ay;       lata_out[2] = az
        latb_out[0] = bx;       latb_out[1] = by;       latb_out[2] = bz
        latc_out[0] = cx;       latc_out[1] = cy;       latc_out[2] = cz

	#Move COM of molecules to new lattice but DONT stretch/rotate geometry
        xyz = [i for i in range(len(self.geometry))]
        strain_xyz = [i for i in range(len(self.geometry))]
	atoms = [i for i in range(len(self.geometry))]	

        for i in range(len(self.geometry)):
       		xyz[i]= [self.geometry[i][0], self.geometry[i][1], self.geometry[i][2]]
	for i in range(len(self.geometry)):
                atoms[i] = self.geometry[i][3]
	atom_num_per_mol = int(len(xyz)/self.num_mols)
        mol_list = [xyz[x:x+atom_num_per_mol] for x in range(0, len(xyz), atom_num_per_mol)]		
	mol1_COM_f = np.dot(lat_mat_f, self.return_COM(mol_list[0]))
	mol2_COM_f = np.dot(lat_mat_f, self.return_COM(mol_list[1]))
	strain_COM_1 = mol1_COM_f[0]*lata_out + mol1_COM_f[1]*latb_out + mol1_COM_f[2]*latc_out
	strain_COM_2 = mol2_COM_f[0]*lata_out + mol2_COM_f[1]*latb_out + mol2_COM_f[2]*latc_out

	for i in range(len(strain_xyz)):
		if i<atom_num_per_mol:
			strain_xyz[i] = xyz[i]-strain_COM_1
		if i>=atom_num_per_mol:
			strain_xyz[i] = xyz[i]-strain_COM_2
        #Set new structure
        new_struct = Structure()
        for i in range(len(strain_xyz)):
                 new_struct.build_geo_by_atom(float(strain_xyz[i][0]),
                                                 float(strain_xyz[i][1]),
                                                 float(strain_xyz[i][2]),
                                                 atoms[i])
        new_struct.set_property('lattice_vector_a', lata_out)
        new_struct.set_property('lattice_vector_b', latb_out)
        new_struct.set_property('lattice_vector_c', latc_out)
        new_struct.set_property('cell_vol', np.dot(lata_out, np.cross(latb_out, latc_out)))
        new_struct.set_property('crossover_type', self.input_structure.get_property('crossover_type'))
        new_struct.set_property('alpha',self.angle(latb_out, latc_out))
        new_struct.set_property('beta', self.angle(lata_out, latc_out))
        new_struct.set_property('gamma', self.angle(lata_out, latb_out))
        return new_struct

    def rand_vstrain(self, lat_mat):
	strain_list = np.random.standard_normal(6) * self.st_dev
	strain_mat = self.get_strain_mat(strain_list)
	print "strain_mat", strain_mat

	strain_A = np.dot(lat_mat.transpose()[0], strain_mat)
	strain_B = np.dot(lat_mat.transpose()[1], strain_mat)
	strain_C = np.dot(lat_mat.transpose()[2], strain_mat)
	return strain_A, strain_B, strain_C
	
    def get_strain_mat(self, strain_list):
	e = strain_list
	s_mat = np.zeros((3,3))
	s_mat[0][0] = 1.0 + e[0]
        s_mat[0][1] = e[5]/2.
        s_mat[0][2] = e[4]/2.

    	s_mat[1][0] = e[5]/2.
        s_mat[1][1] = 1.0 + e[1]
        s_mat[1][2] = e[3]/2.

        s_mat[2][0] = e[4]/2.
        s_mat[2][1] = e[3]/2.
        s_mat[2][2] = 1.0 + e[2]
	return s_mat

    def leng(self,v):
        length = np.linalg.norm(v)
        return length

    def angle(self,v1,v2):
        numdot = np.dot(v1,v2)
        anglerad = np.arccos(numdot/(self.leng(v1)*self.leng(v2)))
        angledeg = anglerad*180/np.pi
        return (angledeg)

    def return_COM(self, geometry):
        '''
        returns COM as np array 
        '''
	N = len(geometry)
	xsum = 0 ; ysum = 0 ; zsum = 0;
	for atom in geometry:
	    xsum +=atom[0]
	    ysum +=atom[1]
	    zsum +=atom[2]
	com = [xsum%N, ysum%N, zsum%N]
        return np.asarray(com)

##############################################################################################
class RandomStrainMutation(object):
    '''
     This mutation gives a random strain to the lattice and DOESNT move the COM of the molecules
    '''
    def __init__(self, input_structure, target_stoic, replica):
        self.geometry = deepcopy(input_structure.get_geometry())
        self.A = np.asarray(deepcopy(input_structure.get_property('lattice_vector_a')))
        self.B = np.asarray(deepcopy(input_structure.get_property('lattice_vector_b')))
        self.C = np.asarray(deepcopy(input_structure.get_property('lattice_vector_c')))
	self.ui = user_input.get_config()
        self.input_structure = input_structure
	self.st_dev = self.ui.get_eval('mutation', 'stand_dev_strain')

    def mutate(self):
        return self.rstrain()

    def rstrain(self):
	A = self.A
	B = self.B
	C = self.C 	
        lat_mat = np.zeros((3,3))
        lat_mat[0][0] = A[0]; lat_mat[0][1] = B[0]; lat_mat[0][2] = C[0]
        lat_mat[1][0] = A[1]; lat_mat[1][1] = B[1]; lat_mat[1][2] = C[1]
        lat_mat[2][0] = A[2]; lat_mat[2][1] = B[2]; lat_mat[2][2] = C[2]

	#Fractional Transformation Matrix
        lat_mat_f = np.linalg.inv(lat_mat)

	#Strained Lattice Vecctors
	strain_A, strain_B, strain_C = self.rand_vstrain(lat_mat)

	#Lengths and angles of strained vectors
        a = self.leng(strain_A)
        b = self.leng(strain_B)
        c = self.leng(strain_C)
        alpha = self.angle(strain_B,strain_C)
        beta = self.angle(strain_A,strain_C)
        gamma = self.angle(strain_A,strain_B)

	#Realign axis to have a along x
        rad=float(np.pi/180)
        ax=a;   ay=0.0;   az=0.0
        bx=np.cos(gamma*rad)*b; by=np.sin(gamma*rad)*a; bz=0.0
        cx=c*np.cos(beta*rad);  cy=(b*c*np.cos(alpha*rad)-bx*cx)/by;    cz=np.sqrt(np.absolute(c**2-cx**2-cy**2))
        lata_out = np.zeros(3); latb_out = np.zeros(3); latc_out = np.zeros(3)
        lata_out[0] = ax;       lata_out[1] = ay;       lata_out[2] = az
        latb_out[0] = bx;       latb_out[1] = by;       latb_out[2] = bz
        latc_out[0] = cx;       latc_out[1] = cy;       latc_out[2] = cz

        #Set new structure
        new_struct = Structure()
        new_struct.build_geo_whole(self.geometry)
	new_struct.set_property('lattice_vector_a', lata_out)
        new_struct.set_property('lattice_vector_b', latb_out)
        new_struct.set_property('lattice_vector_c', latc_out)
        new_struct.set_property('cell_vol', np.dot(lata_out, np.cross(latb_out, latc_out)))
        new_struct.set_property('crossover_type', self.input_structure.get_property('crossover_type'))
        new_struct.set_property('alpha',self.angle(latb_out, latc_out))
        new_struct.set_property('beta', self.angle(lata_out, latc_out))
        new_struct.set_property('gamma', self.angle(lata_out, latb_out))
        return new_struct

    def rand_vstrain(self, lat_mat):
	strain_list = np.random.standard_normal(6) * self.st_dev
	strain_mat = self.get_strain_mat(strain_list)
	print "strain_mat", strain_mat

	strain_A = np.dot(lat_mat.transpose()[0], strain_mat)
	strain_B = np.dot(lat_mat.transpose()[1], strain_mat)
	strain_C = np.dot(lat_mat.transpose()[2], strain_mat)
	return strain_A, strain_B, strain_C
	
    def get_strain_mat(self, strain_list):
	e = strain_list
	s_mat = np.zeros((3,3))
	s_mat[0][0] = 1.0 + e[0]
        s_mat[0][1] = e[5]/2.
        s_mat[0][2] = e[4]/2.

    	s_mat[1][0] = e[5]/2.
        s_mat[1][1] = 1.0 + e[1]
        s_mat[1][2] = e[3]/2.

        s_mat[2][0] = e[4]/2.
        s_mat[2][1] = e[3]/2.
        s_mat[2][2] = 1.0 + e[2]
	return s_mat

    def leng(self,v):
        length = np.linalg.norm(v)
        return length

    def angle(self,v1,v2):
        numdot = np.dot(v1,v2)
        anglerad = np.arccos(numdot/(self.leng(v1)*self.leng(v2)))
        angledeg = anglerad*180/np.pi
        return (angledeg)

    def return_COM(self, geometry):
        '''
        returns COM as np array 
        '''
	N = len(geometry)
	xsum = 0 ; ysum = 0 ; zsum = 0;
	for atom in geometry:
	    xsum +=atom[0]
	    ysum +=atom[1]
	    zsum +=atom[2]
	com = [xsum%N, ysum%N, zsum%N]
        return np.asarray(com)

##################################################################################################
class RandomSymmetryStrainMutationMoveMols(object):
    '''
     This mutation gives a random strain to the lattice and moves the COM of the molecules
    '''
    def __init__(self, input_structure, target_stoic, replica):
        self.geometry = deepcopy(input_structure.get_geometry())
        self.A = np.asarray(deepcopy(input_structure.get_property('lattice_vector_a')))
        self.B = np.asarray(deepcopy(input_structure.get_property('lattice_vector_b')))
        self.C = np.asarray(deepcopy(input_structure.get_property('lattice_vector_c')))
	self.ui = user_input.get_config()
        self.input_structure = input_structure
	self.st_dev = self.ui.get_eval('mutation', 'stand_dev_strain')
	self.num_mols = self.ui.get_eval('unit_cell_settings', 'num_molecules')

    def mutate(self):
        return self.rstrain()

    def rstrain(self):
	A = self.A
	B = self.B
	C = self.C 	
        lat_mat = np.zeros((3,3))
        lat_mat[0][0] = A[0]; lat_mat[0][1] = B[0]; lat_mat[0][2] = C[0]
        lat_mat[1][0] = A[1]; lat_mat[1][1] = B[1]; lat_mat[1][2] = C[1]
        lat_mat[2][0] = A[2]; lat_mat[2][1] = B[2]; lat_mat[2][2] = C[2]

	#Fractional Transformation Matrix
        lat_mat_f = np.linalg.inv(lat_mat)

	#Strained Lattice Vecctors
	strain_A, strain_B, strain_C = self.rand_vstrain(lat_mat)

	#Lengths and angles of strained vectors
        a = self.leng(strain_A)
        b = self.leng(strain_B)
        c = self.leng(strain_C)
        alpha = self.angle(strain_B,strain_C)
        beta = self.angle(strain_A,strain_C)
        gamma = self.angle(strain_A,strain_B)

	#Realign axis to have a along x
        rad=float(np.pi/180)
        ax=a;   ay=0.0;   az=0.0
        bx=np.cos(gamma*rad)*b; by=np.sin(gamma*rad)*a; bz=0.0
        cx=c*np.cos(beta*rad);  cy=(b*c*np.cos(alpha*rad)-bx*cx)/by;    cz=np.sqrt(np.absolute(c**2-cx**2-cy**2))
        lata_out = np.zeros(3); latb_out = np.zeros(3); latc_out = np.zeros(3)
        lata_out[0] = ax;       lata_out[1] = ay;       lata_out[2] = az
        latb_out[0] = bx;       latb_out[1] = by;       latb_out[2] = bz
        latc_out[0] = cx;       latc_out[1] = cy;       latc_out[2] = cz

        #Set new structure
        new_struct = Structure()
        new_struct.build_geo_whole(self.geometry)
        new_struct.set_property('lattice_vector_a', lata_out)
        new_struct.set_property('lattice_vector_b', latb_out)
        new_struct.set_property('lattice_vector_c', latc_out)
        new_struct.set_property('cell_vol', np.dot(lata_out, np.cross(latb_out, latc_out)))
        new_struct.set_property('crossover_type', self.input_structure.get_property('crossover_type'))
        new_struct.set_property('alpha',self.angle(latb_out, latc_out))
        new_struct.set_property('beta', self.angle(lata_out, latc_out))
        new_struct.set_property('gamma', self.angle(lata_out, latb_out))
        return new_struct

    def rand_vstrain(self, lat_mat):
	eta = np.random.standard_normal(1) * self.st_dev 
        strain_list = eta*self.choose_rand_sym_strain(lat_mat)
        strain_mat = self.get_strain_mat(strain_list)
#  	print "eta", eta
        print "strain_mat", strain_mat

        strain_A = np.dot(lat_mat.transpose()[0], strain_mat)
        strain_B = np.dot(lat_mat.transpose()[1], strain_mat)
        strain_C = np.dot(lat_mat.transpose()[2], strain_mat)
        return strain_A, strain_B, strain_C

    def choose_rand_sym_strain(self, lat_mat):
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
#	print "rand_choice", rand_choice
#	print strain_dict[rand_choice] 
	return np.array(strain_dict[rand_choice])
	
    def get_strain_mat(self, strain_list):
	e = strain_list
	s_mat = np.zeros((3,3))
	s_mat[0][0] = 1.0 + e[0]
        s_mat[0][1] = e[5]/2.
        s_mat[0][2] = e[4]/2.

    	s_mat[1][0] = e[5]/2.
        s_mat[1][1] = 1.0 + e[1]
        s_mat[1][2] = e[3]/2.

        s_mat[2][0] = e[4]/2.
        s_mat[2][1] = e[3]/2.
        s_mat[2][2] = 1.0 + e[2]
	return s_mat

    def leng(self,v):
        length = np.linalg.norm(v)
        return length

    def angle(self,v1,v2):
        numdot = np.dot(v1,v2)
        anglerad = np.arccos(numdot/(self.leng(v1)*self.leng(v2)))
        angledeg = anglerad*180/np.pi
        return (angledeg)

    def return_COM(self, geometry):
        '''
        returns COM as np array 
        '''
	N = len(geometry)
	xsum = 0 ; ysum = 0 ; zsum = 0;
	for atom in geometry:
	    xsum +=atom[0]
	    ysum +=atom[1]
	    zsum +=atom[2]
	com = [xsum%N, ysum%N, zsum%N]
        return np.asarray(com)

##############################################################################################