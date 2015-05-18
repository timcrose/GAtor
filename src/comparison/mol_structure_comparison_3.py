'''
Created on Oct 16, 2013

@author: newhouse modified by farren- THIS version as opposed to 1 tracks the same atom in each molecule (not too largest)
'''
from __future__ import division

import math
import numpy

from core import user_input, output


def main(struct, structure_coll, replica):
    '''
    For a structure comparison module, the method main() must be implemented.
    
    The module will take a new structure and compare it to the given structure collection.
    It will return True if the structure passes a test for uniqueness and for energy
    '''

   # print struct.get_property('crossover_type')		
   # a = datetime.datetime.now()
    comp = Comparison(struct, structure_coll, replica)
    result = comp.is_acceptable()
   # b = datetime.datetime.now()
   # if comp.verbose: self.output("time taken to compare: " + str(struct.struct_id) + ' -- ' + str(b - a))
    print result	
    return result

class Comparison:
    def __init__(self, struct, structure_coll, replica):
        if struct is False or struct is None: raise Exception
        self.replica = replica
        self.struct_a = struct
        self.structure_coll = structure_coll
        self.ui = user_input.get_config()
        self.verbose = self.ui.get_eval('run_settings', 'verbose')
	self.angle_up_bound = float(self.ui.get_eval('comparison', 'angle_up_bound'))
	self.angle_low_bound =float(self.ui.get_eval('comparison', 'angle_low_bound'))

    def output(self, message): output.local_message(message, self.replica)

    def is_acceptable(self):
	#Structure goes through series of checks and is accepted or rejected (T/F)
	if self.acceptable_angles() is True:
		if self.acceptable_vol() is False:
			return False
		else:
			return True
	else:
		return False



#	if self.acceptable_fine_vol() is False: return False
#        return True
############# CHECK alpha, beta, gamma 60<*<120 degrees ##########

    def acceptable_angles(self):
	alpha = self.struct_a.get_property('alpha')
	beta = self.struct_a.get_property('beta')
        gamma = self.struct_a.get_property('gamma')
	
	if alpha < self.angle_low_bound or alpha > self.angle_up_bound:
		print "Alpha NOT in range: ", alpha
		return False
	if beta < self.angle_low_bound or beta > self.angle_up_bound:
		print "Beta NOT in range: ", beta
		return False
	if gamma < self.angle_low_bound or gamma > self.angle_up_bound:
		print "Gamma NOT in range: ", gamma
		return False
	else:
		print "Angles in range"
		return True
	

############# VOLUME UNIQUENESS ###############

    def acceptable_vol(self):
        '''
        Modified by Farren: Make sure volume isn't too big or small
        '''
        #User-Specified Initial Volume and percentage tolerance
        vol0 = float(self.ui.get_eval('comparison', 'initial_vol'))
        tol = float(self.ui.get_eval('comparison', 'vol_decimal_tol'))
        tolf = float(self.ui.get_eval('comparison', 'fine_vol_tol'))
        #Get volume of structure in question
        vol = self.struct_a.get_property('cell_vol')

        #If volume is not within bounds reject it
        if vol > vol0 + tol*vol0:
                print "Volume is too big: ", vol
                return False
        elif vol < vol0 - tol*vol0:
                print "Volume is too small: ", vol
                return False
        else:
                print "Volume is within large range: ", vol
                vol_list = []
                for struct in self.structure_coll.itervalues():
			cell_vol = struct.get_property(('cell_vol'))
			diff = numpy.absolute(cell_vol- vol)
                        vol_list.append(diff)
		delta_list = [d for d in vol_list if d <= tolf]
#		print vol_list
		if len(delta_list) == 0:
	        	print "Fine Volume unique!"
                        return True
		else:
                        print "Fine Volume not unique"
                        return False
	


    def acceptable_vol_3(self):
        '''
        Modified by Farren: Make sure volume isn't too big or small
        '''
	#User-Specified Initial Volume and percentage tolerance
	vol0 = float(self.ui.get_eval('comparison', 'initial_vol'))
	tol = float(self.ui.get_eval('comparison', 'vol_decimal_tol'))
	tolf = float(self.ui.get_eval('comparison', 'fine_vol_tol'))    
        #Get volume of structure in question
        vol = self.struct_a.get_property('cell_vol')
	
	#If volume is not within bounds reject it
	if vol > vol0 + tol*vol0:
		print "Volume is too big: ", vol
		return False
	elif vol < vol0 - tol*vol0:
		print "Volume is too small: ", vol
		return False
	else:
		print "Volume is within large range: ", vol
        	vol_list = []
        	for struct in self.structure_coll.itervalues():
            		vol_list.append(struct.get_property(('cell_vol')))
		print vol_list	
		for v in vol_list:
			if v <= vol + tolf and v >= vol - tolf:
				print "Fine Volume not unique"
				return False
			else:
				print "Fine Volume unique!"
				return True

    def acceptable_fine_vol(self, vol):
        '''
        Modified by Farren: Make sure volume isn't equivalent to other structures
        '''

	#Make list of volumes in the collection
        vol_list = []
        for struct in self.structure_coll.itervalues():
            vol_list.append(struct.get_property(('cell_vol'))) 

        #If volume is not within bounds reject it
	vol_list_equal = [v for v in vol_list if v == vol + tolf]
        if len(vol_list_equal)== 0:
                print "Fine volume is unique"
                return True
        else:
                print "Fine volume is equal to another pool"
                return False







############# PREVIOUSLY USED FUNCTIONS ###############
    def compare_distance_array(self, struct_a, list_to_compare):
        '''
        sums the difference in the sorted distances and compares to the user-definied minimum
        '''
        dist_array_tolerance = self.ui.get_eval('comparison', 'dist_array_tolerance')
        da_a = self.calc_distance_array(struct_a)
        for struct_b in list_to_compare:
	    if struct_a == struct_b: raise Exception # debugging feature added the wrong structure to comparison list
            da_b = self.calc_distance_array(struct_b)
            dist_sum = 0.0				
            if len(da_a) != len(da_b): 
                continue  # different distance array lengths. obviously not the same structure
            for i in range(len(da_a)):
                dist_sum += abs(numpy.sqrt((da_a[i] - da_b[i])**2))
            if self.verbose: self.output('distance sum = ' + str(dist_sum) + ' dist_array_tolerance = ' + str(dist_array_tolerance))
	    print "~DIST SUM~"
	    print dist_sum		
            if float(dist_sum) < float(dist_array_tolerance):
                print 'signature not unique'
                return False
       	    else:
		print 'signature unique'
        	return True
        
    def energy_window(self, struct_a, list_to_compare):
        '''
        reduces the list of structures that warrant comparison by filtering those not within
        a certain window of energy defined by the user
        '''
        e_window = self.ui.get_eval('comparison', 'energy_window')
        if e_window == None: return list_to_compare
        return_list = []
        e_a = float(struct_a.get_property('energy'))
        for struct_b in list_to_compare:
            e_b = float(struct_b.get_property('energy'))
            if e_b <= e_a + e_window and e_b >= e_a + e_window: return_list.append(struct_b) 
        return return_list
     
    def acceptable_energy(self):
        '''
        For each structure in the collection, check if the new energy is greater than the maximum existing energy
        '''
        e_a = self.struct_a.get_property('energy')
        if e_a == None: return True  # structure has not yet been relaxed. energy is unknown
        energy_list = []
        for struct in self.structure_coll.itervalues():
            energy_list.append(struct.get_property(('energy')))
        all_energies = int(len(energy_list))
        if all_energies == 0: return True
        # filter for energies lower than e_a
        energy_list_2 = [e_b for e_b in energy_list if e_a > e_b]
        lower_energies = int(len(energy_list_2))
        ratio = lower_energies / all_energies 
#         if ratio < 1.0: # accept energy as long as it's not greater than max
#         if ratio < 0.5: # accept energy as long as it's not greater than half of the existing energies
#         if ratio <= 0.0: # accept energy only if it's lower than all the energies
        if ratio <= self.ui.get_eval('comparison', 'energy_comparison'):  
            return True
        else: 
            if self.verbose: self.output('Energy too high')
            return False


    def acceptable_energy_2(self):
        '''
        Modified by Farren: Check to make sure new structures energy is not identical to existing energies
        '''
	#Get energy of structure in question
        e_a = self.struct_a.get_property('energy')
        if e_a == None: return True  # structure has not yet been relaxed. energy is unknown
       
	#Make list of current energies in collection
	energy_list = []
        for struct in self.structure_coll.itervalues():
            energy_list.append(struct.get_property(('energy')))

	#Make list of structures for which e_a = e_b
	energy_list_2 = [e_b for e_b in energy_list if e_a == e_b]

	#Check to make sure e_a is not equal to any of these energies
	if len(energy_list_2)== 0:
		print "Energy is unique"
                return True
	else:
		print "Energy is equal to another energy in pool"
		return False

		 
    def calc_distance_array(self, struct):  # TODO: this a particular instance of uniqueness checking. move it to a module.
        '''
        creates an UNSORTED list of distances between atoms in the structure
        used for existence checks
        '''
        geo = struct.get_geometry() 
	lats = struct.get_lattice_vectors()
	distance_array = []

	#lattice vector signatures ie a, b, a+b etc
	lata = lats[0]
        latb = lats[1]
        latc = lats[2]
	ax = lata[0]
	ay = lata[1]
	az = lata[2]
	bx = latb[0]
	by = latb[1]
	bz = latb[2]
	cx = latc[0]
	cy = latc[1]
	cz = latc[2]

	a_dist = numpy.sqrt(ax**2+ay**2+az**2)
        b_dist = numpy.sqrt(bx**2+by**2+bz**2)
        c_dist = numpy.sqrt(cx**2+cy**2+cz**2)
        ab_dist = numpy.sqrt((ax+bx)**2+(ay+by)**2+(az+bz)**2)
        ac_dist = numpy.sqrt((ax+cx)**2+(ay+cy)**2+(az+cz)**2)
        bc_dist = numpy.sqrt((bx+cx)**2+(by+cy)**2+(bz+cz)**2)
        abc_dist = numpy.sqrt((ax+bx+cx)**2+(ay+by+cy)**2+(az+bz+cz)**2)

	#new version for molecular crystals r is length of each atom from origin
	r = [float(i) for i in range(len(geo))]
	x = [float(i) for i in range(len(geo))]
	y = [float(i) for i in range(len(geo))]
	z = [float(i) for i in range(len(geo))]

	r = numpy.array(r)
	for i in range(len(geo)):
		r[i] = numpy.sqrt(geo[i]['x']**2+geo[i]['y']**2+geo[i]['z']**2)
		x[i] = geo[i]['x']
		y[i] = geo[i]['y']
		z[i] = geo[i]['z']
	
	#sort list and find indices of two largest distances
#	sorted_r = numpy.sort(r)
   #     max1_r = sorted_r[-1]
  #      max2_r = sorted_r[-2]
 #       max1_ind = numpy.where(r == max1_r)
#	max2_ind = numpy.where(r == max2_r)
#	max1_ind = float(max1_ind[0])
#	max2_ind = float(max2_ind[0])


	#use 8th and 18th molecule (one oxy from each mol) for tracking
	#NOTE: CHANGE FOR NEW SYSTEM

	max1_ind = 8
	max2_ind = 18	

	#Make new distance array from distances of from two max coords to next periodic image
	xyz = numpy.asarray(zip(x,y,z))
	
#	print xyz
#	print xyz[max1_ind]
#	print xyz[max2_ind]

        a_d1 = numpy.sqrt((ax-xyz[max1_ind][0])**2 + (ay-xyz[max1_ind][1])**2 +
                         (az-xyz[max1_ind][2])**2)
	b_d1 = numpy.sqrt((bx-xyz[max1_ind][0])**2 + (by-xyz[max1_ind][1])**2 +
                         (bz-xyz[max1_ind][2])**2)
        c_d1 = numpy.sqrt((cx-xyz[max1_ind][0])**2 + (cy-xyz[max1_ind][1])**2 +
                         (cz-xyz[max1_ind][2])**2)
        ab_d1 = numpy.sqrt((ax+bx-2*xyz[max1_ind][0])**2+(ay+by-2*xyz[max1_ind][1])**2+
                        (az+bz-2*xyz[max1_ind][2])**2)
        ac_d1 = numpy.sqrt((ax+cx-2*xyz[max1_ind][0])**2+(ay+cy-2*xyz[max1_ind][1])**2+
                        (az+cz-2*xyz[max1_ind][2])**2)
        bc_d1 = numpy.sqrt((bx+cx-2*xyz[max1_ind][0])**2+(by+cy-2*xyz[max1_ind][1])**2+
                        (bz+cz-2*xyz[max1_ind][2])**2)
        abc_d1 = numpy.sqrt((ax+bx+cx-3*xyz[max1_ind][0])**2+(ay+by+cy-3*xyz[max1_ind][1])**2+
                        (az+bz+cz-3*xyz[max1_ind][2])**2)

        a_d2 = numpy.sqrt((ax-xyz[max2_ind][0])**2 + (ay-xyz[max2_ind][1])**2 +
                         (az-xyz[max2_ind][2])**2)
        b_d2 = numpy.sqrt((bx-xyz[max2_ind][0])**2 + (by-xyz[max2_ind][1])**2 +
                         (bz-xyz[max2_ind][2])**2)
        c_d2 = numpy.sqrt((cx-xyz[max2_ind][0])**2 + (cy-xyz[max2_ind][1])**2 +
                         (cz-xyz[max2_ind][2])**2)
        ab_d2 = numpy.sqrt((ax+bx-2*xyz[max2_ind][0])**2+(ay+by-2*xyz[max2_ind][1])**2+
                        (az+bz-2*xyz[max2_ind][2])**2)
        ac_d2 = numpy.sqrt((ax+cx-2*xyz[max2_ind][0])**2+(ay+cy-2*xyz[max2_ind][1])**2+
                        (az+cz-2*xyz[max2_ind][2])**2)
        bc_d2 = numpy.sqrt((bx+cx-2*xyz[max2_ind][0])**2+(by+cy-2*xyz[max2_ind][1])**2+
                        (bz+cz-2*xyz[max2_ind][2])**2)
        abc_d2 = numpy.sqrt((ax+bx+cx-3*xyz[max2_ind][0])**2+(ay+by+cy-3*xyz[max2_ind][1])**2+
                        (az+bz+cz-3*xyz[max2_ind][2])**2)

        sig_dists=[a_dist, b_dist, c_dist, ab_dist, ac_dist, bc_dist, abc_dist,
                        a_d1, b_d1, c_d1, ab_d1, ac_d1, bc_d1, abc_d1,
                        a_d2, b_d2, c_d2, ab_d2, ac_d2, bc_d2, abc_d2]
        sort_dists = numpy.sort(sig_dists)


       #old version for just atomic xyz positions
	#could modify to include lat vectors somehow in future
#        for i in range(len(geo)):
#           for j in range(len(geo)):
#              if i > j:
#                 distance = (geo[i]['x'] - geo[j]['x']) ** 2 
#                 + (geo[i]['y'] - geo[j]['y']) ** 2 
#                + (geo[i]['z'] - geo[j]['z']) ** 2
#                distance = math.sqrt(distance)
#              # distance_array[i][j] = distance
#                distance_array.append(distance)
	
#	print sig_dists
	
	distance_array = list(sort_dists)
	biggest_distance = distance_array[0]
	struct.set_property('dist_array', distance_array)
	struct.set_property('biggest_distance', biggest_distance)
 	return distance_array
    

#   def center_geometry(self, geometry):
#
 #       for i in range(2):  # x, y, and z
 #           coordinate_sum = 0
 #           counter = 0
 #           for atom in geometry:
 #               coordinate_sum += atom[i]
 #               counter += 1
 #           # average axis value
 #           average = coordinate_sum / counter
 #        for atom in geometry:
 #             # shift all towards center
 #              atom[i] = atom[i] - average
 #       return geometry
