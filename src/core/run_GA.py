#!/usr/bin/python

from __future__ import division
import datetime
import os
import sys
import time
import numpy as np
src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)) # add source directory to python path
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
        self.control_list = self.ui.get_list('control', 'control_in_filelist')
        self.max_cascade = len(self.control_list) - 1
        for i in range(INITIAL_POOL_REFID, len(self.control_list)):  # initializes the structure collections necessary for cascade
            self.structure_supercoll[(self.replica_stoic, i)] = StructureCollection(self.replica_stoic, i)
        structure_collection.update_supercollection(self.structure_supercoll)
        self.ip_coll = self.structure_supercoll.get((self.replica_stoic, INITIAL_POOL_REFID))
        self.structure_coll = self.structure_supercoll.get((self.replica_stoic, 0))
	self.child_counter = 0
	
	self.success_counter = len(self.structure_coll)
	#self.min_energies_0 = []
	self.singlemutate = self.ui.get_eval('mutation', 'mutation_probability')
	self.doublemutate = self.ui.get_eval('mutation', 'double_mutate_prob')
	#Convergence settings from ui.conf  
	self.full_relax_tol = float(self.ui.get('run_settings', 'full_relax_tol'))   
	self.top_en_count = int(self.ui.get('run_settings', 'number_of_top_energies')) 
	self.max_en_it = int(self.ui.get('run_settings', 'max_iterations_energy'))
	self.number_of_structures = int(self.ui.get('run_settings', 'number_of_structures'))
	self.number_of_IP = int(self.ui.get('run_settings', 'number_of_IP'))
	self.number_of_replicas =int(self.ui.get('parallel_settings', 'number_of_multiprocesses'))
	self.mod_iteration_counter = 0
	self.global_iteration_counter = 0
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
	#initial_pool_module.main(self.replica, self.replica_stoic) #Old spot for IP filling
        structure_collection.update_supercollection(self.structure_supercoll)
        self.output("***************** USER INITIAL POOL FILLED *****************") 

	restart_counter = 0
        while True:
            ########## Beginning of Iteration Tasks ##########
            output.move_to_shared_output(self.replica)
            self.output('Beginning iteration')
            self.ui = user_input.get_config()
            begin_time = datetime.datetime.now()
#            self.restart(str(self.replica)+' '+str(restart_counter)+' started_iteration:  ' +str(datetime.datetime.now()))            
            ########## Check if finished/converged ##########
            try: 
                if len(self.structure_coll.structures)-self.number_of_IP >= self.number_of_structures:
			self.output("~*~*~*~*~*~*~*~*~*~*~*~ GA CONVERGED *~*~*~*~*~*~*~*~*~*~*~*~*~*~*")
 			self.output("Number of additions to pool has reached user-specified number:")
			self.output(str(len(self.structure_coll.structures)-self.number_of_IP))
			self.output("Total size of collection:")
			self.output(str(len(self.structure_coll.structures)))
			self.output("GA ended at: "+str(datetime.datetime.now()))
			return
		if convergeTF == True:
			self.output("~*~*~*~*~*~*~*~*~*~*~*~ GA CONVERGED *~*~*~*~*~*~*~*~*~*~*~*~*~*~*")
			self.output("Top energies haven't changed in user-specfied number of iterations")
			self.output("Total size of collection:")
                        self.output((len(self.structure_coll.structures)))
			self.output("GA ended at: "+str(datetime.datetime.now()))
                        return
            except: pass
            if get_kill() == "kill": 
                self.output('replica ' + str(self.replica) + ' killed')
                return
            	  
	    ########## Structure Selection ########
            # Expects: dictionary_of_on_or_more<Stoic, StructureCollection> #Returns: list_of_2_or_more<Structure>
	    self.output("--Structure selection--")	
            structures_to_cross = selection_module.main(self.structure_supercoll, self.replica_stoic, self.replica)
            if structures_to_cross is False: self.output('Selection failure'); continue
 
	    ########## Mutation Decision ##########
            # Expects: None #Returns: stoichiometry or None. target_stoic in StoicDict format
            if not hasattr(mutation_module, 'get_targ_stoic'): targ_stoic = self.replica_stoic
            else: targ_stoic = mutation_module.get_targ_stoic()

            ############# Crossover ###############
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
	    randnum = np.random.random()	
	    randnum2 = np.random.random()
	    #Single Parents have to be mutated	
	    if new_struct.get_property('cross_type') == [1,1] or new_struct.get_property('cross_type') == [2,2]:
	    	new_struct = mutation_module.main(new_struct, targ_stoic, self.replica)
		self.output("--Second Mutation--")
		if randnum2 < self.doublemutate:
			new_struct = mutation_module.main(new_struct, targ_stoic, self.replica) 
	    else: #Otherwise stick to probabilty set in conf file
		if randnum < self.singlemutate:
			new_struct = mutation_module.main(new_struct, targ_stoic, self.replica)
			if randnum2 < self.doublemutate:
				self.output("--Second Mutation--")
                        	new_struct = mutation_module.main(new_struct, targ_stoic, self.replica)
		else:
			self.output('No mutation applied.')
			new_struct.set_property('mutation_type', 'No_mutation')
            if new_struct is False: 
                 self.output('Mutation failure')
                 continue  # mutation failed, start with new selection

	    ####Structure modification of angles. Checks reasonable structure is fed into relaxation####
	    #self.output("within GA, this is new_struct.properties"+str(new_struct.properties))
	    self.output("--Cell Checks--")	
	    structure_handling.cell_modification(new_struct, self.replica,create_duplicate=False)
	    if not structure_handling.cell_check(new_struct,self.replica): #unit cell considered not acceptable
		continue
 
	    ########### Begin 'Cascade' ##############################
            cascade_counter = 0 
            structures_to_add = {}
            control_check_SPE_string = read_data(os.path.join(cwd, self.ui.get('control', 'control_in_directory')),
                                         self.control_list[0])
            control_relax_full_string = read_data(os.path.join(cwd, self.ui.get('control', 'control_in_directory')),
                                         self.control_list[1])
		
            ########## Check SPE and perform Full Relaxation ##########
            # Expects: Structure, working_dir, input_path #Returns: Structure
	    self.output("--SPE Check and Full Relaxation--")
            self.restart(str(self.replica)+' '+str(restart_counter)+' started_relaxing:    ' +str(datetime.datetime.now()))

	    struct = relaxation_module.main(new_struct, self.working_dir, control_check_SPE_string, control_relax_full_string, self.replica)
            if struct is False: 
	    	self.output('SPE check not passed or full relaxation failure for '+ str(self.replica))
		is_acceptable = False
                continue  # optimization failed, start with new selection

	    self.output("FHI-aims relaxation wall time (s):   " +str(struct.get_property('relax_time')))		
	    self.restart(str(self.replica)+' '+str(restart_counter)+' finished_relaxing:   ' +str(datetime.datetime.now()))
	
	    ########### Comparison Post Relaxation ####################
	    self.output("--Comparison--")
            structure_collection.update_supercollection(self.structure_supercoll)
            is_acceptable = comparison_module.main(struct, self.structure_supercoll.get((struct.get_stoic(), cascade_counter)), self.replica)
            if is_acceptable is False:
            	self.output('Newly relaxed structure is not acceptable') 
                continue  # structure not acceptable start with new selection
		
            #Add structure to a list of structures to be put in to collection
            self.set_parents(structures_to_cross, struct)
            structures_to_add[(struct.get_stoic(), 0)] = struct

	    ########## End of iteration tasks/convergence checks#######
            if is_acceptable is False or struct is False: 
		convergeTF = False
		continue  # start with new selection
            else:
		restart_counter = restart_counter +1 
		convergeTF = self.end_of_iteration_tasks(structures_to_add, self.success_counter, restart_counter)
       	    end_time = datetime.datetime.now()
            self.output("GA Iteration time: -- " + str(end_time - begin_time))
	    self.output("Current Wallclock: -- " + str(end_time))	


############################# Helper Functions used in start function above #######################################################  

    def end_of_iteration_tasks(self, structures_to_add, num_success_common, restart_counter):
        prev_struct_index = None #none for now because only using GA for FF or aims not both
	ID = len(self.structure_coll.structures) + 1
	self.child_counter = self.child_counter + 1 
	energy_list = []

        for key, struct in structures_to_add.items():
            # Set IDs and count numbers
            struct.set_property('prev_struct_id', prev_struct_index)  # tracks cascade sequence
	    struct.set_property('ID', ID)	
	    struct.set_property('replica', self.replica)		

	   #Sort energies of collection before adding new replica 
	    coll = self.structure_supercoll.get((self.replica_stoic,0))
            old_list_top_en = np.array([])
            old_e_list = np.array([])
            for index, structure in coll:
                energy = structure.get_property('energy')
                old_e_list = np.append(energy,old_e_list)
                old_e_list= np.sort(old_e_list.reshape(len(old_e_list),1),axis=0)
                old_list_top_en = old_e_list[:self.top_en_count]
                old_min_e = old_list_top_en[0][0]		

	    #Check if new structure's energy is a new global minima
            e_new = struct.get_property('energy') 		
            self.output("E_stuct: "+str(e_new)+"  GM: "+str(old_min_e))
	    self.check_if_global_minima(old_min_e, e_new) 		

	    #Add Structure to Collection
	    struct.set_property('child_counter', self.child_counter)
	    index = structure_collection.add_structure(struct, key[0], key[1])			
	    structure_collection.update_supercollection(self.structure_supercoll)

	    #Check for Energy Convergence of GA
	    new_list_top_en = np.array([])
	    new_e_list = np.array([])
            for index, structure in self.structure_supercoll.get((self.replica_stoic,0)):
                energy = structure.get_property('energy')
                new_e_list = np.append(energy,new_e_list)
                new_e_list= np.sort(new_e_list.reshape(len(new_e_list),1),axis=0)
                new_list_top_en = new_e_list[:self.top_en_count]
	    self.output("old top energies:    "+ str(old_list_top_en))
	    self.output("new top energies:    "+ str(new_list_top_en))
	    converged = self.check_convergence(old_list_top_en,new_list_top_en,len(self.structure_coll.structures)) 

	    #Output success message to screen and write energy hierarchy        
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

	    #End of Iteration Outputs
	    additions = len(self.structure_coll.structures)-self.number_of_IP	
	    avg_add = float(additions)/self.number_of_replicas	
	    self.restart(str(self.replica)+' '+str(restart_counter-1)+' finished_iteration:  ' +str(datetime.datetime.now()))
	    self.output(str(self.replica)+' finished iteration')
	    self.output('Cumulative additions to common pool: '+str(additions))
	    self.output('Avg addition per replica:            '+str(avg_add))
	    self.output('Total size of common pool:           '+str(len(self.structure_coll.structures)))		 	
	    if converged is "not_converged":
                self.output("GA not converged yet.")
                pass	
	    elif converged is "converged":
		return True	

    def check_if_global_minima(self,min_e, e_new):	
        if e_new < min_e:
        	diff = min_e - e_new
        	message = '*********** NEW GLOBAL MINIMUM FOUND ************' + \
           	'\n  old minima:  ' + str(min_e) + \
           	'\n  new minima:  ' + str(e_new) + \
           	'\n  difference:  ' + str(diff)
        	self.output(message)

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
    def restart(self, message): output.restart_message(message)	

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

