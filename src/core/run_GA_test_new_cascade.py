#!/usr/bin/python

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
from structures import structure_collection, structure_handling
from structures.structure import Structure
from structures.structure_collection import StructureCollection, string_to_stoic
from utilities.stoic_model import determine_stoic
# from external_libs import asizeof
# from profilestats import profile  # for time profiling

def main(replica,stoic):
    # setup
    mkdir_p(tmp_dir)	
    mkdir_p(structure_dir)
    set_unkill()
    # check if going to be multiprocess
    ui = user_input.get_config()
    n_processors = ui.get_eval('run_settings', 'parallel_on_core')
    begin(replica,stoic)

def begin(replica,stoic):
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
                output.error(e, replica)
    output.move_to_shared_output(replica)
    
class RunGA():
    '''
    This is the core algorithm which runs the genetic algorithm tasks
    '''
    def __init__(self, replica, stoic):
        '''
        Initialization
        '''
        # Define class fields
        self.replica = replica
        self.ui = user_input.get_config()
        self.replica_stoic = stoic
        self.structure_supercoll = {}
        self.working_dir = os.path.join(tmp_dir, str(self.replica))
        self.verbose = self.ui.get_eval('run_settings', 'verbose')        
        # read in all other structures to local memory from storage
        self.control_list = self.ui.get_list('control', 'control_in_filelist')
        self.max_cascade = len(self.control_list) - 1
        for i in range(INITIAL_POOL_REFID, len(self.control_list)):  # initializes the structure collections necessary for cascade
            self.structure_supercoll[(self.replica_stoic, i)] = StructureCollection(self.replica_stoic, i)
        structure_collection.update_supercollection(self.structure_supercoll)
        self.ip_coll = self.structure_supercoll.get((self.replica_stoic, INITIAL_POOL_REFID))
        self.structure_coll = self.structure_supercoll.get((self.replica_stoic, 0))
	self.child_counter = 0
	self.success_counter = len(self.structure_coll)
	self.min_energies_0 = []
	self.doublemutate = self.ui.get_eval('mutation', 'double_mutate_prob')
	#Convergence settings from ui.conf  
	self.full_relax_tol = float(self.ui.get('run_settings', 'full_relax_tol'))   
	self.top_en_count = int(self.ui.get('run_settings', 'number_of_top_energies')) 
	self.max_en_it = int(self.ui.get('run_settings', 'max_iterations_energy'))
	self.number_of_structures = int(self.ui.get('run_settings', 'number_of_structures'))
	self.mod_iteration_counter = 0

    def start(self):
        '''
        Performs main genetic algorithm operations
        Loads necessary modules based on UI at runtime
        ''' 
	initial_pool_module = my_import(self.ui.get('modules', 'initial_pool_module'), package='initial_pool')
        selection_module = my_import(self.ui.get('modules', 'selection_module'), package='selection')
        crossover_module = my_import(self.ui.get('modules', 'crossover_module'), package='crossover')
        mutation_module = my_import(self.ui.get('modules', 'mutation_module'), package='mutation')
        relaxation_module = my_import(self.ui.get('modules', 'relaxation_module'), package='relaxation')
        comparison_module = my_import(self.ui.get('modules', 'comparison_module'), package='comparison')

        ########## Fill initial pool ##########
	self.output("--Replica %s updating local pool--" %(self.replica))
	#initial_pool_module.main(self.replica, self.replica_stoic) #The stoichiometry passed in here is collected from ui.conf
        structure_collection.update_supercollection(self.structure_supercoll)
        self.output("***************** USER INITIAL POOL FILLED *****************") 

        while True:
            ########## Beginning of Iteration Tasks ##########
            output.move_to_shared_output(self.replica)
            self.output('Beginning iteration')
            self.ui = user_input.get_config()
            begin_time = datetime.datetime.now()
            
            ########## Check if finished/converged ##########
            try: 
                if len(self.structure_coll.structures) >= self.number_of_structures:
 			self.output("Length of collection has reached user-specified number:")
			self.output(str(len(self.structure_coll.structures)))
			return
		if convergeTF == True:
			self.output("~*~*~*~*~*~*~*~*~*~*~*~*~*~*~* GA CONVERGED *~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*")
			self.output("Top energies haven't changed in user-specfied number of iterations")
			self.output("Length of Collection")
                        self.output((len(self.structure_coll.structures)))
                        return
            except: pass
            if get_kill() == "kill": 
                self.output('replica ' + str(self.replica) + ' killed')
                return
            	  
	     ########## Structure Selection #######
            # Expects: dictionary_of_on_or_more<Stoic, StructureCollection> #Returns: list_of_2_or_more<Structure>
	    self.output("--Structure selection--")	
            structures_to_cross = selection_module.main(self.structure_supercoll, self.replica_stoic, self.replica)
            if structures_to_cross is False: self.output('Selection failure'); continue
	    #stc_message = ''
            #for stc in structures_to_cross: stc_message += stc.get_path() + ' + '
            #self.output('Crossing structures: ' + stc_message[:-3])
             
	    ########## Mutation Decision ##########
            # Expects: None #Returns: stoichiometry or None. target_stoic in StoicDict format
            if not hasattr(mutation_module, 'get_targ_stoic'): targ_stoic = self.replica_stoic
            else: targ_stoic = mutation_module.get_targ_stoic()

            ########## Crossover ##########
            # Expects: list_of_2_or_more<Structure>, target stochiometry #Returns Structure or False
	    self.output("--Crossover--")
            new_struct = crossover_module.main(structures_to_cross, targ_stoic, self.replica)
	    if new_struct is False: 
           	self.output("Crossover failure")
             	continue  # crossover failed, start with new selection
	    #self.output("post crossover geo: ")
	    #self.output(str(new_struct.get_geometry_atom_format()))       	    
			 
            ########## Mutation Execution ##########
            # Expects: Structure, target_stoichiometry [decision] #Returns: Structure
       	    self.output("--Mutation--")  	
	    new_struct = mutation_module.main(new_struct, targ_stoic, self.replica)
            if new_struct is False: 
                 self.output('Mutation failure')
                 continue  # mutation failed, start with new selection
	    #self.output("post mutation geo: ")	
	    #self.output(str(new_struct.get_geometry_atom_format()))	
          
	    ##### Choose possibility of second Mutation ##########
	    randnum = np.random.random()
	    if randnum < self.doublemutate:
		self.output("--Second Mutation--")
		new_struct = mutation_module.main(new_struct, targ_stoic, self.replica)
		#self.output("post second mutation geo: ")
            	#self.output(str(new_struct.get_geometry_atom_format()))
	        if new_struct is False:
                	self.output('Mutation failure')
                 	continue  # mutation failed, start with new selection

	    ####Structure modification of angles. Checks reasonable structure is fed into relaxation####
	    #self.output("within GA, this is new_struct.properties"+str(new_struct.properties))
	    self.output("--Cell Checks--")	
	    structure_handling.cell_modification(new_struct, self.replica,create_duplicate=False)
	    if not structure_handling.cell_check(new_struct,self.replica):
		#unit cell considered not acceptable
		continue
 
	     ########### Begin Cascade #############
            cascade_counter = 0
 	    if self.max_cascade > 0: self.output('Beginning cascade') 
            structures_to_add = {}
            control_check_SPE_string = read_data(os.path.join(cwd, self.ui.get('control', 'control_in_directory')),
                                         self.control_list[0])
            control_relax_full_string = read_data(os.path.join(cwd, self.ui.get('control', 'control_in_directory')),
                                         self.control_list[1])
		
            ########## Check SPE and perform Full Relaxation ##########
            # Expects: Structure, working_dir, input_path #Returns: Structure
	    self.output("--Relaxation--")
	    struct = relaxation_module.main(new_struct, self.working_dir, control_check_SPE_string, control_relax_full_string, self.replica)
            if struct is False: 
	    	self.output('Relaxation failure for replica'+ str(self.replica))
		is_acceptable = False
                break  # optimization failed, start with new selection
                    
            ########## Comparison #############
	    self.output("--Comparison--")
            structure_collection.update_supercollection(self.structure_supercoll)
            is_acceptable = comparison_module.main(struct, self.structure_supercoll.get((struct.get_stoic(), cascade_counter)), self.replica)
            if is_acceptable is False:
            	self.output('Newly relaxed structure is not acceptable') 
                break  # structure not acceptable start with new selection
		
	    #Sort energies of collection for convergence and global minima checks
	    coll = self.structure_supercoll.get((self.replica_stoic, cascade_counter))
            temp_e_list = np.array([])
	    tot_e_list = np.array([])
	    for index, structure in coll:
		energy = structure.get_property('energy')
		tot_e_list = np.append(energy,tot_e_list)
		tot_e_list= np.sort(tot_e_list.reshape(len(tot_e_list),1),axis=0)
		temp_e_list = tot_e_list[:self.top_en_count]
		temp_min_e = temp_e_list[0][0]
		
	    #######Check if needs to be fully relaxed and compare######
	    e_struct = struct.get_property('energy')
	    if np.absolute(temp_min_e)-np.absolute(e_struct) <= self.full_relax_tol:
	    	self.output("Structure close to global minima")
		self.output("E_stuct: "+str(e_struct)+"  GM: "+str(temp_min_e))
		########## Full Relaxation #############
                self.output("--Full Relaxation--")
#                struct = relaxation_module.main(new_struct, self.working_dir, control_relax_full_string, self.replica)
 #               if struct is False:
  #              	self.output('Relaxation failure for replica'+ str(self.replica))
   #                 	is_acceptable = False
    #                	break  # optimization failed, start with new selection
	        ########## Comparison #############
            	self.output("--Full Comparison--")
            	structure_collection.update_supercollection(self.structure_supercoll)
            	is_acceptable = comparison_module.main(struct, self.structure_supercoll.get((struct.get_stoic(), 0)), self.replica)
            	if is_acceptable is False:
                	self.output('Newly relaxed structure is not acceptable')
                	break  # structure not acceptable start with new selection
#		struct.set_property('Relax_type','Full_Relax')
#	    else: struct.set_property('Relax_type','Geom_Relax')
		
            #Add structure to a list of structures to be put in to collection
            self.set_parents(structures_to_cross, struct)
            structures_to_add[(struct.get_stoic(), 0)] = struct

	    ######End of iteration tasks and convergence checks########	
            if is_acceptable is False or struct is False: 
		convergeTF = False
		continue  # start with new selection
            else: 
		convergeTF = self.end_of_iteration_tasks(structures_to_add, temp_min_e, self.success_counter, temp_e_list, tot_e_list)	
       	    end_time = datetime.datetime.now()
            self.output("Iteration time: -- " + str(end_time - begin_time))

    def end_of_iteration_tasks(self, structures_to_add, min_e, num_success_common, old_en_list, tot_en_list):
        prev_struct_index = None #none for now because only using GA for FF or aims not both
	ID = len(self.structure_coll.structures) + 1
	self.child_counter = self.child_counter + 1 
	energy_list = []
	self.top_en_count 
        for key, struct in structures_to_add.items():
            # Set IDs and count numbers
            struct.set_property('prev_struct_id', prev_struct_index)  # tracks cascade sequence
	    struct.set_property('child_count', self.child_counter)
	    struct.set_property('ID', ID)	
	    struct.set_property('replica', self.replica)		

	    #Check if new structure's energy is a new global minima
            e_new = struct.get_property('energy')
	    if e_new < min_e:
		diff = min_e - e_new
	        message = '*********** NEW GLOBAL MINIMUM FOUND ************' + \
                      '\n  old minima:  ' + str(min_e) + \
                      '\n  new minima:  ' + str(e_new) + \
                      '\n  difference:  ' + str(diff)  
		self.output(message)
		struct.set_property('new_local_minima', "True")

	    #Output success message to screen and write energy hierarchy	
	    index = structure_collection.add_structure(struct, key[0], key[1])
            prev_struct_index = str(key) + str(index)
            message = 'Success!: \n  stoichiometry-- ' + key[0].get_string() + \
                      '\n  cascade-- ' + str(key[1]) + \
                      '\n  structure index-- ' + str(index) + \
		      '\n  replica child count-- ' + str(self.child_counter) + \
		      '\n  collection count -- ' + str(ID) + \
		      '\n  replica-- ' + str(self.replica)
            self.output(message)
	    self.output("writing hierachy")	
            data_tools.write_energy_hierarchy(self.structure_coll)

	    #Check for Energy Convergence of GA 
	    old_list_top_en = old_en_list
	    new_list_top_en = np.append(tot_en_list, e_new)
	    new_list_top_en= np.sort(new_list_top_en.reshape(len(new_list_top_en),1),axis=0)
            new_list_top_en = new_list_top_en[:self.top_en_count]	
	    self.output("old top energies:    "+ str(old_list_top_en))
	    self.output("new top energies:    "+ str(new_list_top_en))
	    converged = self.check_convergence(old_list_top_en,new_list_top_en,len(self.structure_coll.structures)) 
	    if converged is "not_converged":
                self.output("GA not converged yet.")
                pass	
	    elif converged is "converged":
		return True	

    def check_convergence(self, old_list_top_en, new_list_top_en, length_coll):
	#Check if top N energies havent changed in X iterations
	self.mod_iteration_counter = self.mod_iteration_counter + 1
	for en_new in new_list_top_en:
 		if en_new in old_list_top_en:
			continue
		else:
                	self.output("Top "+str(self.top_en_count)+" energies have changed.")
			self.mod_iteration_counter = 0
			self.output("Convergence counter reset.")
			self.output("Convergence iteration:  "+ str(self.mod_iteration_counter))
                        return "not_converged"
	self.output("Convergence iteration:  "+ str(self.mod_iteration_counter))
	if self.mod_iteration_counter == self.max_en_it:
		return "converged" 

    def set_parents(self, structures_to_cross, struct):
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
	(options,argv)=argument_opt()
	replica=options.replica
	if replica==None:
		replica=get_random_index()
	print "This is the replica received: "+ str(replica)
	stoic = determine_stoic()
	if stoic == None: raise Exception
	main(replica,stoic)
