'''
Created on Aug 5, 2013

@author: newhouse
'''
from __future__ import division

from copy import deepcopy
from math import cos, sin
import math
import numpy 
import random

from core import user_input, output
from structures.structure import StoicDict, Structure


def main(list_of_structures, targit_stoic, replica):
    '''
    Every crossover module must have a method named "main" that will take
    as arguments a list of 2 or more Structures and a target stoichiometry (in the form of
    a StoicDict). It must return a single Structure with defined geometry
    or return False if the crossover fails and new structures are needed.
    
    This example uses a class to handle the crossover process.
    
    See the Structure class documentation for instructions on easily
    building a Structure's geometry.
    '''
    cross_object = Crossover(list_of_structures[0], list_of_structures[1], targit_stoic, replica)
    return_struct = cross_object.cross()
    if cross_object.verbose and not return_struct is False: 
         cross_object.output('geometries to cross: \n' + list_of_structures[0].get_geometry_atom_format() +\
                        '\n' + list_of_structures[1].get_geometry_atom_format())
         cross_object.output('crossed_geometry: \n' + return_struct.get_geometry_atom_format())
   
  #  print 'crossover mod resultant structure'
  #  print return_struct.get_geometry()
  #  print return_struct.get_lattice_vectors()
    return return_struct
    
class Crossover(object):
    '''
    Crossover takes two structures, randomly rotates, and splices 
    two halves together to create a new un-relaxed geometry
    '''
    
    def __init__(self, struct_a, struct_b, stoic, replica):
        '''
        __init__ will always run when a class is initialized. 
        This is also how the class can take arguments.
        '''
        self.ui = user_input.get_config()
        self.replica = replica
        # deepcopy makes copy of data instead of altering the obect being referenced
        self.geometry_a = deepcopy(struct_a.get_geometry())
        self.geometry_b = deepcopy(struct_b.get_geometry())
        self.target_stoic = deepcopy(stoic)
        self.verbose = self.ui.get_eval('run_settings', 'verbose')
	# Added by farren 12/16/14
	self.latA_a = deepcopy(struct_a.get_property('lattice_vector_a'))
	self.latB_a = deepcopy(struct_a.get_property('lattice_vector_b'))
        self.latC_a = deepcopy(struct_a.get_property('lattice_vector_c'))
        self.latA_b = deepcopy(struct_b.get_property('lattice_vector_a'))
        self.latB_b = deepcopy(struct_b.get_property('lattice_vector_b'))
        self.latC_b = deepcopy(struct_b.get_property('lattice_vector_c'))


    def output(self, message): output.local_message(message, self.replica)

    def cross(self):
        '''
        Rotates and combines the structures so that the resulting structure
        matches the target stoichiometry. 
        '''   
	# Create copies of each parent structure and center their geometries (if not already so)
	temp_a = deepcopy(self.geometry_a)
        temp_b = deepcopy(self.geometry_b)
	self.center_geometry(temp_a)
        self.center_geometry(temp_b)

	#Enter number of molecules in cell- this assumes ordering in geometry file by whole molecules
	num_mol = 2
 	atom_num_per_mol = int(len(temp_a)/num_mol)

	#Separate Molecules from each cell for crossing (need to modify for >2 mol)
	mol1_a = temp_a[0:atom_num_per_mol]
	mol2_a = temp_a[atom_num_per_mol:len(temp_a)]
	mol1_b = temp_b[0:atom_num_per_mol]
	mol2_b = temp_b[atom_num_per_mol:len(temp_b)]
	mol1_a_COM = self.center_of_geo(mol1_a)
	mol2_a_COM = self.center_of_geo(mol2_a)
	mol1_b_COM = self.center_of_geo(mol1_b)
	mol2_b_COM = self.center_of_geo(mol2_b)

	#Take only xyz comps of geometry for manipulation

        mol1a = zip(mol1_a['x'], mol1_a['y'], mol1_a['z'])
	mol2a = zip(mol2_a['x'], mol2_a['y'], mol2_a['z'])
	mol1b = zip(mol1_b['x'], mol1_b['y'], mol1_b['z'])
	mol2b = zip(mol2_b['x'], mol2_b['y'], mol2_b['z'])
	
	print mol1_a
	print mol1a
	print mol2_a
	print mol2a
	print mol1_b
	print mol1b
	print mol2_b
	print mol2b

	#Create copies of each parents lattice vectors and components
        A_a = numpy.asarray(deepcopy(self.latA_a))
        B_a = numpy.asarray(deepcopy(self.latB_a))
        C_a = numpy.asarray(deepcopy(self.latC_a))
        A_b = numpy.asarray(deepcopy(self.latA_b))
        B_b = numpy.asarray(deepcopy(self.latB_b))
        C_b = numpy.asarray(deepcopy(self.latC_b))

	#Create Matrices of Lats to find fractional rep.
        lat_matA = numpy.zeros((3,3))
	lat_matB = numpy.zeros((3,3))
        lat_matA[0][0] = A_a[0]; lat_matA[0][1] = B_a[0]; lat_matA[0][2] = C_a[0]
        lat_matA[1][0] = A_a[1]; lat_matA[1][1] = B_a[1]; lat_matA[1][2] = C_a[1]
        lat_matA[2][0] = A_a[2]; lat_matA[2][1] = B_a[2]; lat_matA[2][2] = C_a[2]
	lat_matB[0][0] = A_b[0]; lat_matB[0][1] = B_b[0]; lat_matB[0][2] = C_b[0]
        lat_matB[1][0] = A_b[1]; lat_matB[1][1] = B_b[1]; lat_matB[1][2] = C_b[1]
        lat_matB[2][0] = A_b[2]; lat_matB[2][1] = B_b[2]; lat_matB[2][2] = C_b[2]
	
	lat_matA_f = numpy.linalg.inv(lat_matA)
	lat_matB_f = numpy.linalg.inv(lat_matB)
	
	#Fractional representation of each COM of molcule
	mol1_a_COM_f = numpy.dot(lat_matA_f, mol1_a_COM)
	mol2_a_COM_f = numpy.dot(lat_matA_f, mol2_a_COM)
	mol1_b_COM_f = numpy.dot(lat_matB_f, mol1_b_COM)
	mol2_b_COM_f = numpy.dot(lat_matB_f, mol2_b_COM)


	#Choose Random Crossover type
	cross_type  = self.choose_crossover_type()
	print "cross_type="
        print cross_type

	#Setup lattice for child 
	child_lats = self.lat_options(cross_type, A_a, B_a, C_a, A_b, B_b, C_b)
	child_latA = child_lats[0]
	child_latB = child_lats[1]
	child_latC = child_lats[2]
	lat_mat = numpy.zeros((3,3))
	lat_mat[0][0] = child_latA[0]; lat_mat[0][1] = child_latB[0]; lat_mat[0][2] = child_latC[0]
	lat_mat[1][0] = child_latA[1]; lat_mat[1][1] = child_latB[1]; lat_mat[1][2] = child_latC[1]
	lat_mat[2][0] = child_latA[2]; lat_mat[2][1] = child_latB[2]; lat_mat[2][2] = child_latC[2]

	#Setup geometry for child 
#	child_geo = self.geo_options(cross_type, mol1_a, mol2_a, mol1_b, mol2_b)
	child_geo = self.geo_options(cross_type, mol1a, mol2a, mol1b, mol2b, mol1_a_COM_f, mol2_a_COM_f, mol1_b_COM_f, mol2_b_COM_f)
	
	#Print parents + child geom  for checking
	print "geo parent A"
	print temp_a
	print "geo parent B"
	print temp_b
	print "child geo"
	print child_geo
	#Printout parent +child lattices for checking
	print "lattice vectors parent A"
        print A_a
	print B_a      
        print C_a
        print "lattice vectors parent B"
        print A_b
        print B_b
        print C_b 
	print "child lattice vectors"
	print child_latA
	print child_latB
	print child_latC


	#Set new structure
	new_struct = Structure()
	new_struct.build_geo_whole(child_geo)
	new_struct.set_property('lattice_vector_a', child_latA)
	new_struct.set_property('lattice_vector_b', child_latB)
	new_struct.set_property('lattice_vector_c', child_latC)
		
	return new_struct



    def choose_crossover_type(self):
	geom_opts = [1, 2, 3, 4]
        lat_opts = [1]
	cross_types = [[1,1],[2,1],[3,1],[4,1]]
	method = random.choice(cross_types)
	return method

    def geo_options(self, crosstype, mol1A, mol2A, mol1B, mol2B, COM_1A, COM_2A, COM_1B, COM_2B):
	'''
        geo opt1 = take both mols from mol A
        geo opt2 = take both mols from mol B
        geo opt3 = take first from A, second from B
        geo opt4 = take the second from A, first from B
	 '''
	geo_choice = crosstype[0]
		
	if geo_choice == 1:
		mol1_geo = numpy.multiply(mol1A, COM_1A)
		mol2_geo = numpy.multiply(mol2A, COM_2A)
	elif geo_choice == 2:
		mol1_geo = numpy.multiply(mol1B, COM_1B)
		mol2_geo = numpy.multiply(mol2B, COM_2B)
	elif geo_choice == 3:
		mol1_geo = numpy.multiply(mol1A, COM_1A)
		mol2_geo = numpy.multiply(mol2B, COM_2B)
	elif geo_choice == 4:
		mol1_geo = numpy.multiply(mol1B, COM_1B)
		mol2_geo = numpy.multiply(mol2A, COM_2A)

	geometry = numpy.concatenate((mol1_geo, mol2_geo), axis=0)
	return geometry
		
    def lat_options(self, crosstype, lat1A, lat1B, lat1C, lat2A, lat2B, lat2C):
        '''
	lat opt1 = take lattice of A
        lat opt2 = take lattice of B
        lat opt3 = choose random fractions to combine A and B
        '''
	lat_choice = crosstype[1]

	if lat_choice == 1:
	#take lattice of random combo of parent A + B
		rand_vec = [random.random() for i in range(3)]

		latA = rand_vec[0]*lat1A + (1 - rand_vec[0])*lat2A
		latB = rand_vec[1]*lat1B + (1 - rand_vec[1])*lat2B
		latC = rand_vec[2]*lat1C + (1 - rand_vec[2])*lat2C


		print "BUTT"
		lats =[latA, latB, latC]
		print lats
	else:
		print "Haven't defined this lattice crossover type yet"
	return lats
#TODO: maybe make other crossover options? Swapping vectors?

	
    def center_of_geo(self, geometry):
        '''
        finds center of geometry but does not shift to origin
        '''
	x_sum = 0
	y_sum = 0
	z_sum = 0
	for i in range(len(geometry)):
		x_sum += geometry[i][0]
		y_sum += geometry[i][1]
		z_sum += geometry[i][2]

	x_sum_av = x_sum/len(geometry)
	y_sum_av = y_sum/len(geometry)
	z_sum_av = z_sum/len(geometry)

        COM = [x_sum_av, y_sum_av, z_sum_av]
	return COM


#MAY NOT NEED ANYMORE

	
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


