#!/usr/bin/python

'''
Created on Jul 19, 2013

@author: newhouse
'''
from __future__ import division

import datetime
import os
import sys
import time
import numpy as np
# add source directory to python path
src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.append(src_dir)

from core import user_input, data_tools, output
from file_handler import *
from kill import *
from structures import structure_collection
from structures.structure import Structure
from structures.structure_collection import StructureCollection, string_to_stoic
from utilities.stoic_model import determine_stoic



# from external_libs import asizeof


# from profilestats import profile  # for time profiling
# @profile
def main(stoic):
    # setup
    mkdir_p(tmp_dir)
    mkdir_p(structure_dir)
    set_unkill()
    # check if going to be multiprocess
    ui = user_input.get_config()
    n_processors = ui.get_eval('run_settings', 'parallel_on_core')
    if not isinstance(n_processors, (int, long)) : begin(stoic)
    else: begin_multiprocess(stoic, int(n_processors))


def begin(stoic):
    replica = get_random_index()
    # run genetic algorithm
    ga = RunGA(replica, stoic)
    # catch crashes
    if ga.ui.get_eval('run_settings', 'recover_from_crashes') is not True: 
        ga.start()  # no exception catching
    else:
        while True:
            try:
                ga.start() 
                break
            except Exception, e: 
                print output.error(e, replica)
    output.move_to_shared_output(replica)
    
class RunGA():
    '''
    DOCUMENT
    '''
    def __init__(self, replica, stoic):
        '''
        Document
        '''

        # Define class fields
        self.replica = replica
        self.ui = user_input.get_config()
        self.replica_stoic = stoic
        self.structure_supercoll = {}
        self.number_of_structures = int(self.ui.get('run_settings', 'number_of_structures'))
        self.working_dir = os.path.join(tmp_dir, str(self.replica))
        self.verbose = self.ui.get_eval('run_settings', 'verbose')        
        # read in all other structures to local memory from storage
        self.control_list = self.ui.get_list('control', 'control_in_filelist')
        self.max_cascade = len(self.control_list) - 1
        for i in range(INITIAL_POOL_REFID, len(self.control_list)):  # initializes the structure collections necessary for cascade
            self.structure_supercoll[(self.replica_stoic, i)] = StructureCollection(self.replica_stoic, i)
        structure_collection.update_supercollection(self.structure_supercoll)
        self.ip_coll = self.structure_supercoll.get((self.replica_stoic, INITIAL_POOL_REFID))
        self.structure_coll = self.structure_supercoll.get((self.replica_stoic, self.max_cascade))
	self.child_counter = 0
	self.success_counter = len(self.structure_coll)
	self.min_energies_0 = []
	self.delta_convg = float(self.ui.get('run_settings', 'delta_convergence'))   
	self.doublemutate = self.ui.get_eval('mutation', 'double_mutate_prob')
 
    def start(self):
        '''
        Performs main genetic algorithm operations
        Loads necessary modules based on UI at runtime
        '''  # TODO: improve documentation
        # Dynamically load modules outside of loop

        initial_pool_module = my_import(self.ui.get('modules', 'initial_pool_module'), package='initial_pool')
        selection_module = my_import(self.ui.get('modules', 'selection_module'), package='selection')
        crossover_module = my_import(self.ui.get('modules', 'crossover_module'), package='crossover')
        mutation_module = my_import(self.ui.get('modules', 'mutation_module'), package='mutation')
        relaxation_module = my_import(self.ui.get('modules', 'relaxation_module'), package='relaxation')
        comparison_module = my_import(self.ui.get('modules', 'comparison_module'), package='comparison')

        ########## Fill initial pool ##########
	print "--Filling IP--"
	initial_pool_module.main(self.replica, self.replica_stoic)
        # update all
        structure_collection.update_supercollection(self.structure_supercoll)
        print "***************** USER INITIAL POOL RELAXED *****************"
        #rm this break after initial pool tests
        #    break            


        while True:
	 ########## Beginning of Iteration Tasks ##########
            output.move_to_shared_output(self.replica)
            if self.verbose: self.output('Beginning iteration')
            self.ui = user_input.get_config()
            begin_time = datetime.datetime.now()
            
            ########## Check if finished ##########
            try: 
                if len(self.structure_coll.structures) >= self.number_of_structures:
 			print "Length of collection"
			print len(self.structure_coll.structures)
			print self.structure_coll.structures
			return
            except: pass
            if get_kill() == "kill": 
                print 'replica ' + str(self.replica) + ' killed'
                self.output('replica ' + str(self.replica) + ' killed')
                return
            	  
	 ########## Perform Genetic Algorithm Tasks ##########
            
            ########## Structure Selection ##########
            # Expects: dictionary_of_on_or_more<Stoic, StructureCollection> #Returns: list_of_2_or_more<Structure>
            if self.verbose: self.output('Beginning structure selection')
	    print "--Structure selection--"	
            structures_to_cross = selection_module.main(self.structure_supercoll, self.replica_stoic, self.replica)
            if structures_to_cross is False: self.output('Selection failure'); continue
            if self.verbose: 
                stc_message = ''
                for stc in structures_to_cross: stc_message += stc.get_path() + ' + '
                self.output('Crossing structures: ' + stc_message[:-3])
             
	    ########## Mutation Decision ##########
            # Expects: None #Returns: stoichiometry or None
            # target_stoic in StoicDict format
            if not hasattr(mutation_module, 'get_targ_stoic'): targ_stoic = self.replica_stoic
            else: targ_stoic = mutation_module.get_targ_stoic()

            ########## Crossover ##########
	 # Expects: list_of_2_or_more<Structure>, target stochiometry #Returns Structure or False
	    print "--Crossover--"
            if self.verbose: self.output('Beginning crossover')
            new_struct = crossover_module.main(structures_to_cross, targ_stoic, self.replica)
	    if new_struct is False: 
           	if self.verbose: self.output('Crossover failure')
             	continue  # crossover failed, start with new selection
#	    print "post crossover geo: "
#	    print new_struct.get_geometry_atom_format()       	    
			 
            ########## Mutation Execution ##########
            # Expects: Structure, target_stoichiometry [decision] #Returns: Structure
       
	    print "--Mutation--"  	
            if self.verbose: self.output('Beginning mutation')
            new_struct = mutation_module.main(new_struct, targ_stoic, self.replica)
            if new_struct is False: 
                 if self.verbose: self.output('Mutation failure')
                 continue  # mutation failed, start with new selection
#	    print "post mutation geo: "	
#	    print new_struct.get_geometry_atom_format()	
          
	    ##### Choose possibility of second Mutation ##########
	    randnum = np.random.random()
	    if randnum < self.doublemutate:
		print "--Second Mutation--"
		new_struct = mutation_module.main(new_struct, targ_stoic, self.replica)
#		print "post second mutation geo: "
#            	print new_struct.get_geometry_atom_format()
	        if new_struct is False:
                	if self.verbose: self.output('Mutation failure')
                 	continue  # mutation failed, start with new selection

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
            ########### Begin Cascade #############
            cascade_counter = 0
 	    if self.max_cascade > 0: self.output('Beginning cascade')
            ########## Comparison ##########
            # Expects: Structure, Structure:collection #Returns: Boolean
            # Checks if structure is unique or does not meet constraints (energy, etc.)
 	    print "--Comparison--"  	  
            structure_collection.update_supercollection(self.structure_supercoll)
            is_acceptable = comparison_module.main(new_struct, self.structure_supercoll.get((self.replica_stoic, cascade_counter)), self.replica)
            if is_acceptable is False: 
                if self.verbose: self.output('Structure to be relaxed is not acceptable')
                continue  # structure not acceptable start with new selection
    
            structures_to_add = {}
#	    child_number = 0
            while True:
                input_string = read_data(os.path.join(cwd, self.ui.get('control', 'control_in_directory')),
                                         self.control_list[cascade_counter])
                ########## Relaxation ##########
                # Expects: Structure, working_dir, input_path #Returns: Structure
		print "--Relaxation--"
                self.output('Beginning relaxation')
		struct = relaxation_module.main(new_struct, self.working_dir, input_string, self.replica)
                if struct is False: 
                    if self.verbose: self.output('Relaxation failure')
                    break  # optimization failed, start with new selection
                    
                ########## Comparison ##########
		print "--Comparison--"
                structure_collection.update_supercollection(self.structure_supercoll)
                is_acceptable = comparison_module.main(struct, self.structure_supercoll.get((struct.get_stoic(), cascade_counter)), self.replica)
                if is_acceptable is False:
                    if self.verbose: self.output('Newly relaxed structure is not acceptable') 
                    break  # structure not acceptable start with new selection
	
		 # Add structure to a list of structures to be put in to collection
                struct.input_ref = cascade_counter
                self.set_parents(structures_to_cross, struct)
                structures_to_add[(struct.get_stoic(), cascade_counter)] = struct

		# Set current minimum energy of collection (0 for now)
		coll = self.structure_supercoll.get((self.replica_stoic, cascade_counter))
 		temp_e_list = []
		for index, structure in coll:
			energy = structure.get_property('energy')
			temp_e_list += [index, energy]

		temp_min_e = min(temp_e_list)

                # End of cascade tasks
                cascade_counter = cascade_counter + 1  # move to next input
                new_struct = struct  # set new structure and return to beginning of cascade
                if cascade_counter > self.max_cascade: break
            	    
            if is_acceptable is False or struct is False: continue  # start with new selection
            else: self.end_of_iteration_tasks(structures_to_add, temp_min_e, self.success_counter)	
            end_time = datetime.datetime.now()
            if self.verbose: self.output("Iteration time: -- " + str(end_time - begin_time))

    def end_of_iteration_tasks(self, structures_to_add, min_e, num_success_IP):
        prev_struct_index = None
	self.child_counter = self.child_counter + 1
	self.struct_counter = self.child_counter + num_success_IP 
	energy_list = []

        for key, struct in structures_to_add.items():
            # Set IDs and count numbers
            struct.set_property('prev_struct_id', prev_struct_index)  # tracks cascade sequence
	    struct.set_property('child_count', self.child_counter)
	    struct.set_property('ID', self.struct_counter)	
		
	    #Check for new Minimum Energy		
            e_new = struct.get_property('energy')
	    if e_new < min_e:
		diff = min_e - e_new
		print "*************************** NEW MINIMUM FOUND **********************************"
		print "old minima:  ", min_e
		print "new minima:  ", e_new
		print "difference:  ", diff
		struct.set_property('new_local_minima', "True")
		struct.set_property('new_local_minima_diff', min_e - e_new)

		if diff < .001:
			print "~CONVERGED ENERGY!!!~"
		else:
			print "~New minimum not converged yet~"

	    #Output success message to screen and write energy hierarchy	
	    index = structure_collection.add_structure(struct, key[0], key[1])
            prev_struct_index = str(key) + str(index)
            message = 'Success: \n  stoichiometry-- ' + key[0].get_string() + \
                      '\n  cascade-- ' + str(key[1]) + \
                      '\n  structure index-- ' + str(index) + \
		      '\n  child count-- ' + str(self.child_counter) + \
		      '\n  collection count -- ' + str(self.struct_counter) + \
		      '\n  replica-- ' + str(self.replica)
            print message
            data_tools.write_energy_hierarchy(self.structure_coll)
            self.output(message)
#             print str(asizeof.asizeof(self.structure_supercoll) / 1048576) + ' MB'
            
    def set_parents(self, structures_to_cross, struct):
        # set the the parent structures as a property for family tree tracing
        for i in range(len(structures_to_cross)): 
            par_st = structures_to_cross[i]
            struct.set_property('parent_' + str(i), par_st.get_stoic_str() + '/' \
                                + str(par_st.get_input_ref()) + '/' + str(par_st.get_struct_id()))

    def output(self, message): output.local_message(message, self.replica)


if __name__ == '__main__':
    '''
    This command is important. If the module is run directly instead of imported,
    it will execute the main() method. This allows for a single replica to join 
    a current genetic algorithm search.
    '''
    try: 
        stoic_filename = sys.argv[1]
        stoic = determine_stoic(os.path.join(cwd, stoic_filename))
    except: 
        stoic = determine_stoic()
    if stoic == None: raise Exception
        
    main(stoic)

