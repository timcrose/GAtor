'''
Created on Aug 5, 2013

@author: farren
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
    '''
    cross_object = Crossover(list_of_structures[0], list_of_structures[1], targit_stoic, replica)
    return_struct = cross_object.cross()
    if cross_object.verbose and not return_struct is False: 
         cross_object.output('geometries to cross: \n' + list_of_structures[0].get_geometry_atom_format() +\
                        '\n' + list_of_structures[1].get_geometry_atom_format())
         cross_object.output('crossed_geometry: \n' + return_struct.get_geometry_atom_format())
  
  #  print return_struct.get_property('crossover_type') 
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
	self.num_mols = self.ui.get_eval('unit_cell_settings', 'num_molecules')

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
	num_mol = self.num_mols
 	atom_num_per_mol = int(len(temp_a)/num_mol)
	mol_listA= [temp_a[x:x+atom_num_per_mol] for x in range(0, len(temp_a), atom_num_per_mol)]
	mol_listB = [temp_b[x:x+atom_num_per_mol] for x in range(0, len(temp_b), atom_num_per_mol)]	

	#Separate Molecules from each cell for crossing (need to modify for >2 mol)
	mol1_a = mol_listA[0]
	mol2_a = mol_listA[1]
	mol1_b = mol_listB[0]
	mol2_b = mol_listB[1]

	#Create copies of each parents lattice vectors 
        A_a = numpy.asarray(self.latA_a)
        B_a = numpy.asarray(self.latB_a)
        C_a = numpy.asarray(self.latC_a)
        A_b = numpy.asarray(self.latA_b)
        B_b = numpy.asarray(self.latB_b)
        C_b = numpy.asarray(self.latC_b)

	#Choose Random Crossover type
	cross_type  = self.choose_crossover_type()
	print "cross_type:     ", cross_type

	#Setup geo for child 
	child_geo = self.geo_options(cross_type, mol1_a, mol2_a, mol1_b, mol2_b)
        mol1_c = child_geo[0:atom_num_per_mol]
        mol2_c = child_geo[atom_num_per_mol:len(temp_a)]
	tooclose = self.is_too_close(mol1_c, mol2_c)
	if tooclose is True:
		print "!child mols too close together!"
		return False
	#Setup lattice for child
	child_lats = self.lat_options(cross_type, A_a, B_a, C_a, A_b, B_b, C_b)
	
	child_latA = child_lats[0]
	child_latB = child_lats[1]
	child_latC = child_lats[2]
	
	temp_vol = numpy.dot(numpy.cross(child_latA, child_latB), child_latC)
	alpha = self.angle(child_latB, child_latC)
	beta = self.angle(child_latA, child_latC)
	gamma = self.angle(child_latA, child_latB)

	#Printout parents + child geom  for checking
#	print "geo parent A", temp_a
#	print "geo parent B", temp_b
#	print "child geo", child_geo

	#Printout parent +child lattices for checking
#	print "lattice vectors parent A"
#       print A_a
#	print B_a      
#       print C_a
#       print "lattice vectors parent B"
#       print A_b
#       print B_b
#       print C_b 
#	print "child lattice vectors"
#	print child_latA
#	print child_latB
#	print child_latC

	#Printout parent + child lengths and angles
	print "Parent A lats:  ", numpy.sqrt(numpy.dot(A_a, A_a)), numpy.sqrt(numpy.dot(B_a, B_a)), numpy.sqrt(numpy.dot(C_a, C_a))
	print "Parent B lats:  ", numpy.sqrt(numpy.dot(A_b, A_b)), numpy.sqrt(numpy.dot(B_b, B_b)), numpy.sqrt(numpy.dot(C_b, C_b))
	print "Child lats:     ", numpy.sqrt(numpy.dot(child_latA, child_latA)), numpy.sqrt(numpy.dot(child_latB, child_latB)), numpy.sqrt(numpy.dot(child_latC, child_latC))
	print "Child angles:   ", alpha, beta, gamma


	

	#Set new structure
	new_struct = Structure()
	new_struct.build_geo_whole(child_geo)
	new_struct.set_property('lattice_vector_a', child_latA)
	new_struct.set_property('lattice_vector_b', child_latB)
	new_struct.set_property('lattice_vector_c', child_latC)
	new_struct.set_property('cell_vol', temp_vol)
	new_struct.set_property('crossover_type', cross_type)
	new_struct.set_property('alpha',alpha)
	new_struct.set_property('beta', beta)
	new_struct.set_property('gamma', gamma)	
	return new_struct



    def choose_crossover_type(self):
	geom_opts = [1, 2, 3, 4]
        lat_opts = [1, 2, 3]
	cross_types = [[1,1],[1,2],[1,3],[2,1],[2,2],[2,3],[3,1],[3,2],[3,3],[4,1],[4,2],[4,3]]
	method = random.choice(cross_types)
	return method

    def geo_options(self, crosstype, mol1A, mol2A, mol1B,  mol2B):
	'''
        geo opt1 = take both mols from mol A
        geo opt2 = take both mols from mol B
        geo opt3 = take first from A, second from B
        geo opt4 = take the second from A, first from B
	 '''
	geo_choice = crosstype[0]
		
	if geo_choice == 1:
		geometry = numpy.concatenate((mol1A, mol2A), axis=0)
	elif geo_choice == 2:
		geometry = numpy.concatenate((mol1B, mol2B), axis=0)
	elif geo_choice == 3:
		geometry = numpy.concatenate((mol1A, mol2B), axis=0)
	else:
		geometry = numpy.concatenate((mol1B, mol2A), axis=0)

	return geometry
		
    def lat_options(self, crosstype, lat1A, lat1B, lat1C, lat2A, lat2B, lat2C):
        '''
	lat opt3 = choose random fractions to combine A and B
        lat opt1 = take whole lattice of A
        lat opt2 = take whole lattice of B
        '''
	lat_choice = crosstype[1]

	if lat_choice == 3:
	#take lattice of random combo of parent A + B
		rand_vec = [random.random() for i in range(3)]
		latA = rand_vec[0]*lat1A + (1 - rand_vec[0])*lat2A
		latB = rand_vec[1]*lat1B + (1 - rand_vec[1])*lat2B
		latC = rand_vec[2]*lat1C + (1 - rand_vec[2])*lat2C
		lats =[latA, latB, latC]	
	elif lat_choice == 1:
		lats = [lat1A, lat1B, lat1C]
	elif lat_choice == 2:
		lats = [lat2A, lat2B, lat2C]
	else:
		print "Haven't defined this lattice crossover type yet"
	return lats
#TODO: maybe make other crossover options? Swapping vectors?

    def is_too_close(self, mol1, mol2):
        '''
        only compares the distances between structures and not within each half.
        relaxed geometries may have distances smaller than the minimum distance.
        '''
        is_too_close = False
        min_dist = float(self.ui.get('crossover', 'crossover_minimum_interface_distance'))
        for atom_a in mol1:
            ax = atom_a[0]
            ay = atom_a[1]
            az = atom_a[2]
            for atom_b in mol2:
                bx = atom_b[0]
                by = atom_b[1]
                bz = atom_b[2]
                dist = numpy.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)
                if dist < min_dist:
                    is_too_close = True
        return is_too_close	

    def angle(self,v1,v2):
    	numdot = numpy.dot(v1,v2)
        anglerad = numpy.arccos(numdot/(self.leng(v1)*self.leng(v2)))
        angledeg = anglerad*180/numpy.pi
        return (angledeg)
    def leng(self,v):
        return((v[0]**2+v[1]**2+v[2]**2)**0.5)	


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
