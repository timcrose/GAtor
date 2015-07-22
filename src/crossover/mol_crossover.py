'''
Created 2015

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
from structures import structure_handling

def main(list_of_structures, targit_stoic, replica):
    '''
    Every crossover module must have a method named "main" that will take
    as arguments a list of 2 or more Structures and a target stoichiometry (in the form of
    a StoicDict). It must return a single Structure with defined geometry
    or return False if the crossover fails and new structures are needed.
    '''
	
    num_mols = user_input.get_config().get_eval('unit_cell_settings', 'num_molecules')
    max_attempts = int(user_input.get_config().get_eval('crossover', 'max_attempts'))	
    cross_attempts = 0

    while cross_attempts != max_attempts:
	output.local_message("Attempting crossover: "+str(cross_attempts+1), replica)
	if num_mols == 2:
        	cross_object = Crossover_2mol(list_of_structures[0], list_of_structures[1], targit_stoic, replica)
        	return_struct = cross_object.cross()
    	elif num_mols == 4:
        	cross_object = Crossover_4mol(list_of_structures[0], list_of_structures[1], targit_stoic, replica)
        	return_struct = cross_object.cross()
    	elif num_mols == 8:
        	cross_object = Crossover_8mol(list_of_structures[0], list_of_structures[1], targit_stoic, replica)
        	return_struct = cross_object.cross()
	cell_check = structure_handling.cell_check(return_struct, replica) #unit cell considered not acceptable
	if cell_check == False:
		cross_attempts = cross_attempts + 1
	if cell_check == True:
		output.local_message("Crossover: "+str(cross_attempts+1)+" successful" , replica)
		cross_attempts = max_attempts

    return return_struct

######## Two Molecule per cell Crossover ########
class Crossover_2mol(object):
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
	message = 'cross_type:     ' + str(cross_type)     
	self.output(message)

	#Setup geo for child 
	child_geo = self.geo_options(cross_type, mol1_a, mol2_a, mol1_b, mol2_b)
        mol1_c = child_geo[0:atom_num_per_mol]
        mol2_c = child_geo[atom_num_per_mol:len(temp_a)]
	
	#Closeness Check
#	tooclose = self.is_too_close(mol1_c, mol2_c)
#	if tooclose is True:
#		self.output("!child mols too close together!")
#		return False

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
        message = 'Parent A lattice vectors:  ' + str(numpy.linalg.norm(A_a))+' '+str(numpy.linalg.norm(B_a))+' '+str(numpy.linalg.norm(C_a))+\
		  '\nParent B lattice vectors:  ' + str(numpy.linalg.norm(A_b))+' '+str(numpy.linalg.norm(B_b))+' '+str(numpy.linalg.norm(C_b))+\
		  '\nChild lattice vectors:  ' + str(numpy.linalg.norm(child_latA))+' '+str(numpy.linalg.norm(child_latB))+' '+str(numpy.linalg.norm(child_latC))+\
		  '\nChild angles:    ' +str(alpha)+' '+str(beta)+' '+str(gamma) 		
        self.output(message)

	#Set new structure
	new_struct = Structure()
	new_struct.build_geo_whole(child_geo)
	new_struct.set_property('lattice_vector_a', child_latA)
	new_struct.set_property('lattice_vector_b', child_latB)
	new_struct.set_property('lattice_vector_c', child_latC)
	new_struct.set_property('a', numpy.linalg.norm(child_latA))
        new_struct.set_property('b', numpy.linalg.norm(child_latB))
        new_struct.set_property('c', numpy.linalg.norm(child_latC))
	new_struct.set_property('cell_vol', temp_vol)
	new_struct.set_property('crossover_type', cross_type)
	new_struct.set_property('alpha',alpha)
	new_struct.set_property('beta', beta)
	new_struct.set_property('gamma', gamma)	
	return new_struct

    def choose_crossover_type(self):
	geom_opts = [1, 2, 3, 4]
        lat_opts = [1, 2, 3]
	cross_types = [[1,2],[1,3],[2,1],[2,3],[3,1],[3,2],[3,3],[4,1],[4,2],[4,3]]
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
		rand_vec = [random.uniform(0.25,0.75) for i in range(3)]
		self.output("Random Frac: " +str(rand_vec))	
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
	length = numpy.linalg.norm(v)
        return length	
	
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


######## Four Molecule per cell Crossover ########
class Crossover_4mol(object):
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

	#Separate Molecules from each cell for crossing (need to modify for number of  mol)
	mol1_a = mol_listA[0]
	mol2_a = mol_listA[1]
	mol3_a = mol_listA[2]
	mol4_a = mol_listA[3]
	# ADD 3rd and 4th molecules here ####
	mol1_b = mol_listB[0]
	mol2_b = mol_listB[1]
	mol3_b = mol_listB[2]
	mol4_b = mol_listB[3]

	#Create copies of each parents lattice vectors 
        A_a = numpy.asarray(self.latA_a)
        B_a = numpy.asarray(self.latB_a)
        C_a = numpy.asarray(self.latC_a)
        A_b = numpy.asarray(self.latA_b)
        B_b = numpy.asarray(self.latB_b)
        C_b = numpy.asarray(self.latC_b)

	#Choose Random Crossover type
	cross_type  = self.choose_crossover_type()
	self.output( "cross_type:     "+str(cross_type))

	#Setup geo for child 
	child_geo = self.geo_options(cross_type, mol1_a, mol2_a, mol3_a, mol4_a, mol1_b, mol2_b, mol3_b, mol4_b)
        mol1_c = child_geo[0:atom_num_per_mol]
	mol2_c = child_geo[atom_num_per_mol: 2*atom_num_per_mol]
	mol3_c = child_geo[2*atom_num_per_mol: 3*atom_num_per_mol]
        mol4_c = child_geo[3*atom_num_per_mol:len(temp_a)]

	#Closeness Checks
#        tooclose = self.is_too_close(mol1_c, mol2_c)
#	if tooclose is True:
#		self.output("!child mols 1&2 too close together!")
#		return False
#	elif self.is_too_close(mol2_c, mol3_c) is True:
#		self.output("!child mols 2&3 too close together!")
#                return False
#        elif self.is_too_close(mol3_c, mol4_c) is True:
#                self.output("!child mols 3&4 too close together!")
#                return False
#        elif self.is_too_close(mol4_c, mol1_c) is True:
#                self.output("!child mols 4&1 too close together!")
#                return False

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

	#Printout parent + child lengths and angles
        message = 'Parent A lattice vectors:  ' + str(numpy.linalg.norm(A_a))+' '+str(numpy.linalg.norm(B_a))+' '+str(numpy.linalg.norm(C_a))+\
                  '\nParent B lattice vectors:  ' + str(numpy.linalg.norm(A_b))+' '+str(numpy.linalg.norm(B_b))+' '+str(numpy.linalg.norm(C_b))+\
                  '\nChild lattice vectors:  ' + str(numpy.linalg.norm(child_latA))+' '+str(numpy.linalg.norm(child_latB))+' '+str(numpy.linalg.norm(child_latC))+\
                  '\nChild angles:    ' +str(alpha)+' '+str(beta)+' '+str(gamma)
        self.output(message)

	#Set new structure
	new_struct = Structure()
	new_struct.build_geo_whole(child_geo)
	new_struct.set_property('lattice_vector_a', child_latA)
	new_struct.set_property('lattice_vector_b', child_latB)
	new_struct.set_property('lattice_vector_c', child_latC)
	new_struct.set_property('a', numpy.linalg.norm(child_latA))
        new_struct.set_property('b', numpy.linalg.norm(child_latB))
        new_struct.set_property('c', numpy.linalg.norm(child_latC))
	new_struct.set_property('cell_vol', temp_vol)
	new_struct.set_property('crossover_type', cross_type)
	new_struct.set_property('alpha',alpha)
	new_struct.set_property('beta', beta)
	new_struct.set_property('gamma', gamma)	
	return new_struct

    def choose_crossover_type(self):
	geom_opts = [1, 2, 3, 4, 5, 6, 7, 8]
        lat_opts = [1, 2, 3]
	cross_types = [[1,2],[1,3],[2,1],[2,3],[3,1],[3,2],[3,3],[4,1],[4,2],[4,3],[5,1],[5,2],[5,3],[6,1],[6,2],[6,3],[7,1],[7,2],[7,3],[8,1],[8,2],[8,3]]
	method = random.choice(cross_types)
	return method

    def geo_options(self, crosstype, mol1A, mol2A, mol3A, mol4A, mol1B,  mol2B, mol3B, mol4B):
	'''
        geo opt1 = take 4 mols from mol A
        geo opt2 = take 2 mols from mol B
        geo opt3 = take top 2 from A, bot 2 from B
        geo opt4 = take bottom 2 from A, top 2 from B
	geo opt5 = take mol 1 and 4 from A, 2 and 3 from B
	geo opt6 = take mol 1 and 4 from B, 2 and 3 from A
	geo opt7 = take mol 1 and 3 from A, 2 and 4 from B
	geo opt8 = take mol 1 and 3 from B, 2 and 4 from B
	 '''
	geo_choice = crosstype[0]
		
	if geo_choice == 1:
		geometry = numpy.concatenate((mol1A, mol2A, mol3A, mol4A), axis=0)
	elif geo_choice == 2:
		geometry = numpy.concatenate((mol1B, mol2B, mol3B, mol4B), axis=0)
	elif geo_choice == 3:
		geometry = numpy.concatenate((mol1A, mol2A, mol3B, mol4B), axis=0)
	elif geo_choice == 4:
		geometry = numpy.concatenate((mol1B, mol2B, mol3A, mol4A), axis=0)
        elif geo_choice == 5:
                geometry = numpy.concatenate((mol1A, mol2B, mol3B, mol4A), axis=0)
        elif geo_choice == 6:
                geometry = numpy.concatenate((mol1B, mol2A, mol3A, mol4B), axis=0)
        elif geo_choice == 7:
                geometry = numpy.concatenate((mol1A, mol2B, mol3A, mol4B), axis=0)
        elif geo_choice == 8:
                geometry = numpy.concatenate((mol1B, mol2A, mol3B, mol4A), axis=0)
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
		rand_vec = [random.uniform(0.25,0.75) for i in range(3)]
                self.output("Random Frac: " +str(rand_vec))
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
	length = numpy.linalg.norm(v)
        return length	
	
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

####################################################################################
######## Eight Molecule per cell Crossover ########
class Crossover_8mol(object):
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

	#Separate Molecules from each cell for crossing (need to modify for number of  mol)
	mol1_a = mol_listA[0]
	mol2_a = mol_listA[1]
	mol3_a = mol_listA[2]
	mol4_a = mol_listA[3]
	mol5_a = mol_listA[4]
	mol6_a = mol_listA[5]	
	mol7_a = mol_listA[6]
	mol8_a = mol_listA[7]

	# ADD Extra molecules here ####
	mol1_b = mol_listB[0]
	mol2_b = mol_listB[1]
	mol3_b = mol_listB[2]
	mol4_b = mol_listB[3]
        mol5_b = mol_listB[4]
        mol6_b = mol_listB[5]
        mol7_b = mol_listB[6]
        mol8_b = mol_listB[7]

	#Create copies of each parents lattice vectors 
        A_a = numpy.asarray(self.latA_a)
        B_a = numpy.asarray(self.latB_a)
        C_a = numpy.asarray(self.latC_a)
        A_b = numpy.asarray(self.latA_b)
        B_b = numpy.asarray(self.latB_b)
        C_b = numpy.asarray(self.latC_b)

	#Choose Random Crossover type
	cross_type  = self.choose_crossover_type()
	self.output( "cross_type:     "+str(cross_type))

	#Setup geo for child 
	child_geo = self.geo_options(cross_type, mol1_a, mol2_a, mol3_a, mol4_a, mol1_b, mol2_b, mol3_b, mol4_b)
        mol1_c = child_geo[0: 1*atom_num_per_mol]
	mol2_c = child_geo[1*atom_num_per_mol: 2*atom_num_per_mol]
	mol3_c = child_geo[2*atom_num_per_mol: 3*atom_num_per_mol]
        mol4_c = child_geo[3*atom_num_per_mol: 4*atom_num_per_mol]
	mol5_c = child_geo[4*atom_num_per_mol: 5*atom_num_per_mol]
	mol6_c = child_geo[5*atom_num_per_mol: 6*atom_num_per_mol]
	mol7_c = child_geo[6*atom_num_per_mol: 7*atom_num_per_mol]
	mol8_c = child_geo[7*atom_num_per_mol: len(temp_a)]


	# Check if Child Molecules are too close together
#	tooclose = self.is_too_close(mol1_c, mol2_c)
#	if tooclose is True:
#		self.output("!child mols 1&2 too close together!")
#		return False
#	elif self.is_too_close(mol2_c, mol3_c) is True:
#                self.output("!child mols 2&3 too close together!")
#                return False
#        elif self.is_too_close(mol3_c, mol4_c) is True:
#                self.output("!child mols 3&4 too close together!")
#                return False
#	elif self.is_too_close(mol4_c, mol5_c) is True:
#                self.output("!child mols 4&5 too close together!")
#                return False
#        elif self.is_too_close(mol5_c, mol6_c) is True:
#	        self.output("!child mols 5&6 too close together!")
#        	return False
#        elif self.is_too_close(mol6_c, mol7_c) is True:
#                self.output("!child mols 6&7 too close together!")
#                return False
#        elif self.is_too_close(mol7_c, mol8_c) is True:
#                self.output("!child mols 7&8 too close together!")
#                return False
#        elif self.is_too_close(mol8_c, mol1_c) is True:
#                self.output("!child mols 8&1 too close together!")
#                return False

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

	#Printout parent + child lengths and angles
        message = 'Parent A lattice vectors:  ' + str(numpy.linalg.norm(A_a))+' '+str(numpy.linalg.norm(B_a))+' '+str(numpy.linalg.norm(C_a))+\
                  '\nParent B lattice vectors:  ' + str(numpy.linalg.norm(A_b))+' '+str(numpy.linalg.norm(B_b))+' '+str(numpy.linalg.norm(C_b))+\
                  '\nChild lattice vectors:  ' + str(numpy.linalg.norm(child_latA))+' '+str(numpy.linalg.norm(child_latB))+' '+str(numpy.linalg.norm(child_latC))+\
                  '\nChild angles:    ' +str(alpha)+' '+str(beta)+' '+str(gamma)
        self.output(message)

	#Set new structure
	new_struct = Structure()
	new_struct.build_geo_whole(child_geo)
	new_struct.set_property('lattice_vector_a', child_latA)
	new_struct.set_property('lattice_vector_b', child_latB)
	new_struct.set_property('lattice_vector_c', child_latC)
	new_struct.set_property('a', numpy.linalg.norm(child_latA))
        new_struct.set_property('b', numpy.linalg.norm(child_latB))
        new_struct.set_property('c', numpy.linalg.norm(child_latC))
	new_struct.set_property('cell_vol', temp_vol)
	new_struct.set_property('crossover_type', cross_type)
	new_struct.set_property('alpha',alpha)
	new_struct.set_property('beta', beta)
	new_struct.set_property('gamma', gamma)	
	return new_struct

    def choose_crossover_type(self):
	geom_opts = [1, 2, 3, 4, 5, 6, 7, 8]
        lat_opts = [1, 2, 3]
	cross_types = [[1,2],[1,3],[2,1],[2,3],[3,1],[3,2],[3,3],[4,1],[4,2],[4,3],[5,1],[5,2],[5,3],[6,1],[6,2],[6,3],[7,1],[7,2],[7,3],[8,1],[8,2],[8,3]]
	method = random.choice(cross_types)
	return method

    def geo_options(self, crosstype, mol1A, mol2A, mol3A, mol4A, mol5A, mol6A, mol7A, mol8A, mol1B, mol2B, mol3B, mol4B, mol5B, mol6B, mol7B, mol8B):
	'''
        geo opt1 = take 8 mols from mol A
        geo opt2 = take 8 mols from mol B
        geo opt3 = take top 4 from A, bot 4 from B
        geo opt4 = take bottom 4 from A, top 4 from B
	geo opt5 = take every other from A, other from B
	geo opt6 = take every other from B, other from A
	geo opt7 = take every other two from A, other from B
	geo opt8 = take every other two from B, other from A
	 '''
	geo_choice = crosstype[0]
		
	if geo_choice == 1:
		geometry = numpy.concatenate((mol1A, mol2A, mol3A, mol4A, mol5A, mol6A, mol7A, mol8A), axis=0)
	elif geo_choice == 2:
		geometry = numpy.concatenate((mol1B, mol2B, mol3B, mol4B, mol5B, mol6B, mol7B, mol8B), axis=0)
	elif geo_choice == 3:
		geometry = numpy.concatenate((mol1A, mol2A, mol3A, mol4A, mol5B, mol6B, mol7B, mol8B), axis=0)
	elif geo_choice == 4:
		geometry = numpy.concatenate((mol1B, mol2B, mol3B, mol4B, mol5A, mol6A, mol7A, mol8A), axis=0)
        elif geo_choice == 5:
                geometry = numpy.concatenate((mol1A, mol2B, mol3A, mol4B, mol5A, mol6B, mol7A, mol8B), axis=0)
        elif geo_choice == 6:
                geometry = numpy.concatenate((mol1B, mol2A, mol3B, mol4A, mol5B, mol6A, mol7B, mol8A), axis=0)
        elif geo_choice == 7:
                geometry = numpy.concatenate((mol1A, mol2A, mol3B, mol4B, mol5A, mol6A, mol7B, mol8B), axis=0)
        elif geo_choice == 8:
                geometry = numpy.concatenate((mol1B, mol2B, mol3A, mol4A, mol5B, mol6B, mol7A, mol8A), axis=0)
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
		rand_vec = [random.uniform(0.25,0.75) for i in range(3)]
                self.output("Random Frac: " +str(rand_vec))
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
	length = numpy.linalg.norm(v)
        return length	
	
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

