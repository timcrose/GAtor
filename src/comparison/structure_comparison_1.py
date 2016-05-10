"""
Created on Tues May 10 11:03:12 2016

@author: Farren Curtis
"""
from __future__ import division

import math
import numpy as np

from core import user_input, output


def main(struct, structure_coll, replica):
    '''
    The module takes a structure and compare it to the given structure collection.
    Returns: True if the structure passes a test for uniqueness and for energy
    '''

    comp = Comparison(struct, structure_coll, replica)

    t1 = datetime.datetime.now()
  
    en_result = comp.acceptable_energy() # make sure energy is higher than the worst in the collection
    if en_result is False:
        return False


    structs_to_compare = comp.get_similar_structures() # return list of structures within a difference tolerance of comparison (.5 eV)
    output.local_message(("Number of Structures w/in duplicate energy window: "+len(structs_to_compare),replica)
    output.local_message(("Structures w/in energy window: "+str(structs_to_compare),replica)

    #is_duplicate = comp.check_if_duplicate(structs_to_compare) #Boolean
   # if is_duplicate is False:
    #    return False
  
    t2 = datetime.datetime.now()
   
    output.local_message(("The structure compared is unique. ",replica)
    output.local_message(("Time taken to compare structure to collection: " + str(struct.struct_id) + ' -- ' + str(t2 - t1),replica)
    return True


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
        sorted_ens = np.sort(np.array(worst_energy))
        worst_energy =sorted_ens[-1] 

        self.output("worst energy: " +str(worst_energy))

        if energy < worst_energy:
            self.output("Structure has unacceptable energy, higher than entire collection")
            return False
        elif energy >= worst_energy:
            return True
    
        
    def get_similar_structures(self,list_to_compare):
        '''
        reduces the list of structures that are checked for duplicates within
        a certain window of energy defined by the user
        '''
        e_tol = float(self.ui.get_eval('comparison', 'dup_energy_tol'))
        if e_tol == None: raise Exception

        sim_list = []
        en = float(self.struct.get_property('energy'))
        for comp_struct in list_to_compare:
            comp_en = float(comp_struct.get_property('energy'))
            if comp_en <= en + e_tol and comp_en >= en + e_tol: 
                sim_list.append(comp_struct) 
        return sim_list
    

        

