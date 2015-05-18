'''
Created on Oct 16, 2013

@author: newhouse modified by farren
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
   # a = datetime.datetime.now()
    comp = Comparison(struct, structure_coll, replica)
    result = comp.is_acceptable()
   # b = datetime.datetime.now()
   # if comp.verbose: self.output("time taken to compare: " + str(struct.struct_id) + ' -- ' + str(b - a))
    print 'comparison mod'
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

    def output(self, message): output.local_message(message, self.replica)

    def is_acceptable(self):
        # structure has not yet been relaxed. energy is unknown
        # unwise to spend time comparing to existing structures, likely random crossover is not contained here
        # comment out to remove this behavior
        if self.struct_a.get_property('energy') == None: return True 

        if self.acceptable_energy() is False: 
            return False
        if self.is_unique() is False: return False
        return True
    
    def is_unique(self):
        '''
        determines whether a structure is unique to the current pool
        '''
        list_to_compare_1 = self.structure_coll.structures.values()
        
	# filter by energy criteria
        list_to_compare_2 = self.energy_window(self.struct_a, list_to_compare_1)

        # filter those which passed by distance array
        if self.verbose: 
            try: self.output('comparing ' + str(int(len(list_to_compare_2) / len(list_to_compare_1) * 100)) + '% of structures at fine level')
            except: pass
        result = self.compare_distance_array(self.struct_a, list_to_compare_2)
	
	

	print list_to_compare_1
	print list_to_compare_2
        return result
    
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
                dist_sum += abs(da_a[i] - da_b[i])
            if self.verbose: self.output('distance sum = ' + str(dist_sum) + ' dist_array_tolerance = ' + str(dist_array_tolerance))
            if float(dist_sum) < float(dist_array_tolerance):
                if self.verbose: self.output('not unique')
                return False
        if self.verbose: self.output('structure unique')
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
        
   # def get_sorted_distance_array(self, struct):
   #     distance_array = self.get_dist_array(struct)
   #     distance_array.sort()
   #     return distance_array
    
#    def get_dist_array(self, struct):
#         '''
#         Returns the distance array from the structures properties.
#         If not found, calculates and adds to properties
#         '''
#         struct.set_property('dist_array', self.calc_distance_array(struct))
#         return struct.get_property('dist_array')
            
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
	sorted_r = numpy.sort(r)
        max1_r = sorted_r[-1]
        max2_r = sorted_r[-2]
        max1_ind = numpy.where(r == max1_r)
	max2_ind = numpy.where(r == max2_r)
	max1_ind = float(max1_ind[0])
	max2_ind = float(max2_ind[0])
	

	#Make new distance array from distances of from two max coords to next periodic image
	xyz = numpy.asarray(zip(x,y,z))
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
	struct.set_property('dist_array', distance_array)
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
