"""
Created on Tues May 10 11:03:12 2016

@author: Farren Curtis
"""
from __future__ import division

import os
import math
import numpy as np
from copy import deepcopy
from core import user_input, output
from core.file_handler import cwd, tmp_dir
from time import time
from structures import structure_collection
from structures.structure_collection import StructureCollection
from pymatgen import Lattice as LatticeP
from pymatgen import Structure as StructureP
from pymatgen.analysis.structure_matcher import StructureMatcher,ElementComparator
from pymatgen.analysis.structure_matcher import SpeciesComparator,FrameworkComparator
from core import output

import multiprocessing

def main(struct, structure_coll, replica, comparison_type):
    '''
    The module takes a structure and compare it to the given structure collection.
    Returns: True/False if the structure passes/fails test for uniqueness and for energy
    '''
    if comparison_type == "pre_relaxation_comparison":
	comp = Comparison(struct, structure_coll, replica, comparison_type)
	structs_to_compare = comp.get_all_structures()
        dup_result = comp.check_if_duplicate_multiprocessing(struct, structs_to_compare, comparison_type)
    elif comparison_type == "post_relaxation_comparison":
        comp = Comparison(struct, structure_coll, replica, comparison_type)
        # make sure energy is higher than the worst in the collection
#        en_result = comp.acceptable_energy()
#        if en_result is False:
#            return False
	structs_to_compare = comp.get_similar_energy_structures(comparison_type)
	dup_result = comp.check_if_duplicate_multiprocessing(struct, structs_to_compare, comparison_type)
    if dup_result:
        output.local_message("-- The structure compared is unique. ", replica)
    return dup_result # Boolean


class Comparison:

    def __init__(self, struct, structure_coll, replica, comparison_type):
        if struct is False or struct is None: raise Exception
        self.replica = replica
        self.struct = deepcopy(struct)
        self.structure_coll = structure_coll
        self.ui = user_input.get_config()
	self.comparison_type = comparison_type

    def output(self, message): output.local_message(message, self.replica)

    def acceptable_energy(self):
        energy = self.struct.get_property('energy')
        if energy == None: raise Exception

        energies = []
        for index, comp_struct in self.structure_coll:
            energies.append(comp_struct.get_property('energy'))
        sorted_ens = np.sort(np.array(energies))
        worst_energy =sorted_ens[-1] 

	self.output("Energy of current structure: " + str(energy))
        self.output("Highest energy in collection: " + str(worst_energy))

        if energy > worst_energy:
            self.output("Structure has unacceptable energy that is higher than entire collection")
            return False
        elif energy <= worst_energy:
	    self.output("Structure has acceptable energy.")
            return True

    def get_all_structures(self):
        '''
	returns full list of structures w/o index
        '''
        struct_list = []
        for index, struct in self.structure_coll:
                struct_list.append((index, struct))
        return struct_list
        
    def get_similar_energy_structures(self, comparison_type):
        '''
        reduces the list of structures that are checked for duplicates within
        a certain window of energy defined by the user
        '''
        e_tol = float(self.ui.get_eval(comparison_type, 'energy_comp_window'))
        if e_tol == None: raise Exception

        sim_list = []
        en = float(self.struct.get_property(self.ui.get_property_to_optimize()))
	for index, comp_struct in self.structure_coll:
            comp_en = float(comp_struct.get_property(self.ui.get_property_to_optimize()))
            if en - e_tol <= comp_en <= en + e_tol:
#		self.output("comp en: " +str(comp_en)) 
                sim_list.append((index, comp_struct)) 
        self.output("-- Number of structures w/in duplicate energy tolerance: "+str(len(sim_list)))
        return sim_list

    def check_if_duplicate(self, struct, comp_list, comparison_type):
	dup_pair = []
	dup_output = open(os.path.join(tmp_dir, "GA_duplicates.dat"),'a')
	for indexc, structc in comp_list:
            fit = self.compute_pymatgen_fit(struct, structc, comparison_type)
            if fit:
                self.output("-- Structure is a duplicate of another in common pool")
		self.output("-- Structure ID in Common pool is: %s" % indexc)
		index = structure_collection.add_structure(struct, struct.get_stoic(), 'duplicates')
		self.output("-- Duplicate Structure ID in duplicates pool is: %s" % index)
		dup_pair.append(("0/"+ str(indexc),"duplicates/"+str(index)))
	        for pair in dup_pair:
	            dup_output.write('\t'.join(str(s) for s in pair) + '\n')
	        return False
        self.ui.grant_permission(os.path.join(tmp_dir,"GA_duplicates.dat"))
	return True

    def check_if_duplicate_multiprocessing(self, struct, comp_list, comparison_type):
        global pool
	processes = self.ui.get_multiprocessing_processes()
        pool = multiprocessing.Pool(processes)
	self.output("-- Comparison done with %i parallel processes" % processes)

        global is_dup
	is_dup = None

	runs = []
        for indexc, structc in comp_list:
            runs.append(pool.apply_async(compute_pymatgen_fit, 
                                         args = [struct, structc, comparison_type], 
                                         callback = compute_pymatgen_fit_callback))
        pool.close()
        pool.join()

        if is_dup == None:
	#Not a duplicate
            return True

	dup_output = open(os.path.join(tmp_dir, "GA_duplicates.dat"),'a')
        self.output("-- Structure is a duplicate of another in common pool")
        self.output("-- Structure ID in Common pool is: %s" % is_dup)
        index = structure_collection.add_structure(struct, struct.get_stoic(), 'duplicates')
	self.output("-- Duplicate Structure ID in duplicates pool is: %s" % index)
        pair = ("0/"+ str(is_dup),"duplicates/"+str(index))
	dup_output.write('\t'.join(str(s) for s in pair) + '\n')
	return False
        

#    def compute_pymatgen_fit(self, struct, structc, comparison_type):
#            sm = self.set_comp_structure_matcher(comparison_type)
#            structp = self.get_pymatgen_structure(struct.get_frac_data())
#            structpc = self.get_pymatgen_structure(structc.get_frac_data())
#            fit = sm.fit(structp, structpc)
#            return fit

    def check_if_duplicate_2(self, comp_list, comparison_type):
        '''
        Args: list of Structures() to compare
        Returns: T/F is structure is duplicate
        '''
        sm = self.set_comp_structure_matcher(comparison_type)
 
	TF_list = []
        structp = self.get_pymatgen_structure(self.struct)
        for comp_struct in comp_list: 
            comp_frac_data = comp_struct.get_frac_data()
            comp_structp = self.get_pymatgen_structure(comp_struct)
            fitTF = sm.fit(structp,comp_structp)
            TF_list.append(fitTF)
                    #print TF_list
        try:            
            if True not in TF_list:
                print "Structure is non-duplicate!"
                print "Total Checked: "+str(len(comp_list)) 
                return True
        except:
            self.output("Structure compared found to be a duplicate")
            self.output("Total Checked"+ str(len(comp_list))) 
            return False


    def set_comp_structure_matcher(self, comparison_type):
        '''
        Args: self
        Returns: Pymatgen StructureMatcher object
        '''
        ui= self.ui
        L_tol =ui.get_eval(comparison_type, 'ltol')
        S_tol = ui.get_eval(comparison_type, 'stol')
        Angle_tol = ui.get_eval(comparison_type, 'angle_tol')
        Scale = ui.get_boolean(comparison_type, 'scale_vol')
        sm = (StructureMatcher(ltol=L_tol, stol=S_tol, angle_tol=Angle_tol, primitive_cell=True, 
                          scale=Scale, attempt_supercell=False, comparator=SpeciesComparator()))
        return sm

    def get_pymatgen_structure(self, frac_data):
        '''
        Args: self, Geometric data from GAtor's Structure() object
        Returns: A pymatgen StructureP() object with the same geometric properties
        '''
        coords = frac_data[0] # frac coordinates
        atoms = frac_data[1] # site labels
        lattice = (LatticeP.from_parameters(a=frac_data[2], b=frac_data[3], c=frac_data[4], 
                                  alpha=frac_data[5],beta=frac_data[6], gamma=frac_data[7]))
        structp = StructureP(lattice, atoms, coords)
        return structp
    

        
def compute_pymatgen_fit(s1, s2, comparison_type):
	ui = user_input.get_config()
	L_tol =ui.get_eval(comparison_type, 'ltol')
	S_tol = ui.get_eval(comparison_type, 'stol')
	Angle_tol = ui.get_eval(comparison_type, 'angle_tol')
	Scale = ui.get_boolean(comparison_type, 'scale_vol')
	sm =  (StructureMatcher(ltol=L_tol, 
				stol=S_tol, 
				angle_tol=Angle_tol, 
				primitive_cell=True,
				scale=Scale, 
				attempt_supercell=False, 
				comparator=SpeciesComparator()))

	sp1 = s1.get_pymatgen_structure()
	sp2 = s2.get_pymatgen_structure()
	fit = sm.fit(sp1, sp2)

	if fit:
		return s2.struct_id
	return fit

def compute_pymatgen_fit_callback(result):
	global pool
	global is_dup
	if result != False:
		pool.terminate()
		#Stop further comparison from being done
		is_dup = result
		#Record struct_id of duplicate structure
