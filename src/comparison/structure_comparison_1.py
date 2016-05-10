"""
Created on Tues May 10 11:03:12 2016

@author: Farren Curtis
"""
from __future__ import division

import math
import numpy as np

from core import user_input, output
from datetime import datetime

from pymatgen import Lattice as LatticeP
from pymatgen import Structure as StructureP
from pymatgen.analysis.structure_matcher import StructureMatcher,ElementComparator,SpeciesComparator,FrameworkComparator

def main(struct, structure_coll, replica):
    '''
    The module takes a structure and compare it to the given structure collection.
    Returns: True if the structure passes a test for uniqueness and for energy
    '''

    comp = Comparison(struct, structure_coll, replica)

    t1 = datetime.now()
  
    en_result = comp.acceptable_energy() # make sure energy is higher than the worst in the collection
    if en_result is False:
        return False

    structs_to_compare = comp.get_similar_structures() # return list of structures within a difference tolerance of comparison (.5 eV)
    dup_result = comp.check_if_duplicate(structs_to_compare)

    t2 = datetime.now()
    output.local_message("The structure compared is unique. ",replica)
    output.local_message("Time taken to compare structure to collection: " + str(struct.struct_id) + ' -- ' + str(t2 - t1),replica)

    return dup_result # Boolean


class Comparison:

    def __init__(self, struct, structure_coll, replica):
        if struct is False or struct is None: raise Exception
        self.replica = replica
        self.struct = struct
        self.structure_coll = structure_coll
        self.ui = user_input.get_config()

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
    
        
    def get_similar_structures(self):
        '''
        reduces the list of structures that are checked for duplicates within
        a certain window of energy defined by the user
        '''
        e_tol = float(self.ui.get_eval('comparison_settings', 'dup_energy_tol'))
        if e_tol == None: raise Exception

        sim_list = []
        en = float(self.struct.get_property('energy'))
        for comp_struct in self.structure_coll.itervalues():
            comp_en = float(comp_struct.get_property('energy'))
            if comp_en <= en + e_tol and comp_en >= en + e_tol: 
                sim_list.append(comp_struct) 
        self.output("Number of Structures w/in duplicate energy window: "+str(len(sim_list)))
        self.output("Structures w/in energy window: "+str(sim_list))
        return sim_list

    def check_if_duplicate(self, comp_list):
        '''
        Args: list of Structures() to compare
        Returns: T/F is structure is duplicate
        '''
        sm = self.set_comp_structure_matcher()
 
        structp = self.get_pymatgen_structure(struct)
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


    def set_comp_structure_matcher(self):
        '''
        Args: self
        Returns: Pymatgen StructureMatcher object
        '''
        ui= self.ui
        L_tol =ui.get_eval('comparison_settings', 'ltol')
        S_tol = ui.get_eval('comparison_settings', 'stol')
        Angle_tol = ui.get_eval('comparison_settings', 'angle_tol')
        Scale = ui.get_eval('comparison_settings', 'scale_vol')
        sm = StructureMatcher(ltol=L_tol, stol=S_tol, angle_tol=Angle_tol, primitive_cell=True, scale=Scale, attempt_supercell=False, comparator=SpeciesComparator())
        return sm

    def get_pymatgen_structure(self, struct):
        '''
        Args: self, Geometric data from GAtor's Structure() object
        Returns: A pymatgen StructureP() object with the same geometric properties
        '''
        frac_data = struct.get_frac_data()
        coords = frac_data[0] # frac coordinates
        atoms = frac_data[1] # site labels
        lattice = LatticeP.from_parameters(a=frac_data[2], b=frac_data[3], c=frac_data[4], alpha=frac_data[5],beta=frac_data[6], gamma=frac_data[7])
        structp = StructureP(lattice, atoms, coords)
        return structp
    

        

