#!/usr/bin/python

from __future__ import division
import datetime
import os
import sys
import time
import random
import numpy as np
src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)) # add source directory to python path
sys.path.append(src_dir)

from core import user_input, data_tools, output, activity
from file_handler import *
from kill import *
from structures import structure_collection, structure_handling
from structures.structure import Structure
from structures.structure_collection import StructureCollection, string_to_stoic
from utilities.stoic_model import determine_stoic
from utilities.parse_aims_output import ParseOutput
from utilities import misc
from utilities import space_group_utils as sgu
from selection import structure_selection
import copy, shutil, multiprocessing


def main(replica,stoic):
	mkdir_p(tmp_dir)	
	mkdir_p(structure_dir)
	ga = RunGA(replica, stoic)	
	if not ga.ui.get_boolean('run_settings', 'recover_from_crashes'): 
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
    '''This class controls the main the genetic algorithm tasks each replica runs'''
    def __init__(self, replica, stoic):
        self.replica = replica
        self.ui = user_input.get_config()
        self.replica_stoic = stoic
        self.working_dir = os.path.join(tmp_dir, str(self.replica))
        self.GA_module_init()
        self.clustering_module = my_import(self.ui.get('modules','clustering_module'), package='clustering')
        self.verbose = self.ui.verbose()
        self.prop = self.ui.get("run_settings","property_to_optimize")
        self.op_style = self.ui.get("run_settings","optimization_style")
        self.number_of_GA_structures = int(self.ui.get('run_settings', 'end_GA_structures_added'))
        self.number_of_tot_structures = int(self.ui.get('run_settings', 'end_GA_structures_total'))
        self.top_count = self.ui.get_eval('run_settings', 'followed_top_structures')
        self.max_it = self.ui.get_eval('run_settings', 'max_iterations_no_change')
        self.replica_child_count = 0
        self.convergence_count= 0
        self.structure_supercoll = {}
        self.structure_supercoll[(self.replica_stoic, 0)] = \
        structure_collection.get_collection(self.replica_stoic, 0)
        self.structure_supercoll[(self.replica_stoic, 'duplicates')] = \
        structure_collection.get_collection(self.replica_stoic, 'duplicates')
        data_tools.write_energy_hierarchy(self.structure_supercoll[(self.replica_stoic,0)])
        structure_collection.update_supercollection(self.structure_supercoll)
        self.structure_coll = structure_collection.stored_collections[(self.replica_stoic, 0)]
        #sname = "parallel_settings"
        self.processes = self.ui.get_multiprocessing_processes()
        if self.processes > 1:
            self.worker_pool = multiprocessing.Pool(processes=self.processes)

    def output(self, message):
        ''' Output messages to replica output in replica_out folder'''
        output.local_message(message, self.replica)

    def start(self):
        '''
        Performs main genetic algorithm operations
        Loads necessary modules based on UI at runtime
        ''' 
        # Report Replica to Common Output 
        self.output("----Replica %s running GA on common pool----" %(self.replica))

        # Intialiaze restarts
        restart_replica = self.ui.get_boolean("run_settings","restart_replicas")
        restart_count = 0
        convergeTF = None

        # Optionally Cluster Input Collection
        if self.ui.get_boolean("clustering","cluster_pool"):
            struct_coll = self.structure_supercoll.get((self.replica_stoic, 0))
            self.clustering_module.main(struct_coll, self.replica)
            data_tools.write_energy_hierarchy(struct_coll)
        # Update Supercollection of Structures
        structure_collection.update_supercollection(self.structure_supercoll)

        # Main loop of GA
        while True:
            #----- Beginning of Iteration Tasks -----#
            begin_time = self.beginning_tasks(restart_count)
            
            #----- Check if GA finished/converged -----#
            end = self.check_finished(convergeTF)
            if end: return
            
            #----- Restart Scheme -----#
            struct = self.restart_scheme(restart_replica)
            if struct == False:
                #-----Generate Trial Structure -----#
                struct = self.generate_trial_structure()
                if struct == False:
                    continue

            #----- Compare Pre-relaxed Structure to Collection -----#
            if self.structure_comparison(struct, "pre_relaxation_comparison") == False:
                convergeTF = False
                rmdir_silence(self.working_dir)
                continue

            #----- Relax Trial Structure -----#
            if self.ui.get_boolean("run_settings", "skip_energy_evaluations"):
                struct.set_property('energy',0.0)
            else:
                mkdir_p(self.working_dir) # make relaxation directory in tmp
                struct = self.structure_relax(struct)
                if struct == False: #relaxation failed, start with new selection
                    rmdir_silence(self.working_dir)
                    continue
                #----- Compare Post-relaxed Structure to Collection -----#
                if self.structure_comparison(struct, "post_relaxation_comparison") == False:
                    convergeTF = False
                    rmdir_silence(self.working_dir)
                    continue

            #---- Compute Spacegroup of Relaxed Structure ----#
            struct = sgu.identify_space_group(struct)
            self.output("-- Structure's space group: %s" % (struct.get_property('space_group')))

            #---- Check If Energy is Global Minimum -----#
            ref_label = 0
            coll = self.structure_supercoll.get((self.replica_stoic,ref_label))
            top_prop_list = self.return_property_array(coll)
            self.check_global_optimization(struct,top_prop_list)

            #----- Add Structure to collection -----#
            self.struct_index = self.add_to_collection(struct, ref_label)

            #----- Optionally Save Output Data -----#
            try:
                self.save_unrelaxed_geometry()
                self.save_relaxation_data()
                self.save_aims_output()
                self.save_replica_output()
            except: pass
            #----- Success Message -----#
            self.success_message()

            #----- Optional Clustering ----#
            if self.ui.get_boolean("clustering","cluster_pool"):
                struct_coll = self.structure_supercoll.get((self.replica_stoic, 0))
                self.clustering_module.main(struct_coll, self.replica)
                structure_collection.update_supercollection(self.structure_supercoll)

            #----- End of Iteration Data Tasks -----#
            restart_count += 1
            convergeTF = self.end_of_iteration_tasks(restart_count, top_prop_list, begin_time, coll)


    def save_replica_output(self):
        input_ref = "0"
        out_path = os.path.join(structure_dir, self.replica_stoic.get_string(), \
                       input_ref, self.struct_index, "creation.out")
        replica_out_path = os.path.join(out_tmp_dir, self.replica + ".out")
        shutil.copyfile(replica_out_path, out_path)
        self.output("-- Saving replica output")

    def save_unrelaxed_geometry(self):
        input_ref = "0"
        geo_orig = "geometry.unrelaxed"
        out_path = os.path.join(structure_dir, self.replica_stoic.get_string(), \
                       input_ref, self.struct_index)
        geo_orig_path = os.path.join(self.working_dir, "geometry.in")
        shutil.copyfile(geo_orig_path, os.path.join(out_path, geo_orig))
        self.output("-- Saving original/unrelaxed child geometry") 

    def save_relaxation_data(self):
        input_ref = "0" 
        list_of_Properties_to_get = self.ui.get_list("FHI-aims", "save_relaxation_data")
        aims_out = os.path.join(self.working_dir, "aims.out")
        po = ParseOutput(list_of_Properties_to_get, self.working_dir)
        po.parseFile(aims_out)
        save_path = os.path.join(tmp_dir, self.replica, "relaxation_data")
        out_path = os.path.join(structure_dir, self.replica_stoic.get_string(), \
                   input_ref, self.struct_index, "relaxation_data")
        shutil.copytree(save_path, out_path)
        self.output("-- Saving relaxation data")

    def save_aims_output(self):
        input_ref = "0"
        if self.ui.get_boolean("FHI-aims","save_aims_output"):
            path = os.path.join(tmp_dir, self.replica)
            files = [i for i in os.listdir(path) if 'aims' in i]
            out_path = os.path.join(structure_dir, self.replica_stoic.get_string(), \
                       input_ref, self.struct_index)
            for file in files:
                save_path = os.path.join(tmp_dir, self.replica, file)
                shutil.copyfile(save_path, os.path.join(out_path,file))
            self.output("-- Saving full aims relaxation ouput")

    def GA_module_init(self):
        '''This routine reads in the modules defined in ui.conf'''
        self.relaxation_module = my_import(self.ui.get('modules', 'relaxation_module'), package='relaxation')
        self.comparison_module = my_import(self.ui.get('modules', 'comparison_module'), package='comparison')
        try: self.clustering_module = my_import(self.ui.get('modules','clustering_module'), package='clustering')
        except: pass

    def beginning_tasks(self, restart_count):
        st = ' -------------------------------------------------------------------------'
        output.move_to_shared_output(self.replica)
        begin_time = datetime.datetime.now()
        self.output(st)
        self.output('|  Replica %s Beginning New Iteration: %s |' % (self.replica, begin_time))
        self.output(st)
        return begin_time

    def return_property_array(self, coll):
        prop_list = np.array([])
        for index, structure in coll:
            if self.op_style == "minimize":
                prop = structure.get_property(self.prop)
            else:
                prop = -structure.get_property(self.prop)
                #Reverse for sorting
            prop_list = np.append(prop,prop_list)
            prop_list = np.sort(prop_list.reshape(len(prop_list),1),axis=0)
        if self.op_style == "maximize":
            prop_list = [-x for x in prop_list] #Reverse it back
        return prop_list

    def check_global_optimization(self,struct,prop_list):
		prop = struct.get_property(self.prop)
		glob = prop_list[0][0] #Best current property
		self.output("-- Structure's %s: %f eV \n-- Previous global minimum: %f eV" % 
		(self.prop, prop, glob))
		diff = abs(prop-glob)
		message = ""
		if self.op_style=="minimize" and prop<glob:
			message = '*********** NEW GLOBAL MINIMUM FOUND ************' + \
			'\n  Old minimum:  ' + str(glob) + \
			'\n  New minimum:  ' + str(prop) + \
			'\n  Difference:  ' + str(diff)
		if self.op_style=="maximize" and prop>glob:
			message = '*********** NEW GLOBAL MAXIMUM FOUND ************' + \
			'\n  Old maximum:  ' + str(glob) + \
			'\n  New maximum:  ' + str(prop) + \
			'\n  Difference:  ' + str(diff)
		self.output(message)
    
    def add_to_collection(self, struct, ref_label):
        '''
        Return index of structure sucessfully added 
        to the common pool
        '''	
        prev_struct_index = None #None for now
        ID = len(StructureCollection(self.replica_stoic, 0).structures) + 1
        self.replica_child_count = self.replica_child_count + 1
        struct.set_property('ID', ID) 
        try:
            struct.set_property('prev_struct_id', prev_struct_index)
            struct.set_property('ID', ID)
            struct.set_property('replica', self.replica)
            struct.set_property('child_counter', self.replica_child_count)
        except: pass

        #UpdateSupercollection/Database
        struct_index = structure_collection.add_structure(struct, self.replica_stoic, ref_label)
        struct_coll = structure_collection.get_collection(self.replica_stoic, ref_label)
        struct_coll.update_local()
        structure_collection.update_supercollection(self.structure_supercoll)
        data_tools.write_energy_vs_addition(struct_coll)
        data_tools.write_energy_hierarchy(struct_coll)
        return struct_index

    def check_finished(self, convergeTF):
        end = False
        IP_dat = os.path.join(tmp_dir,"num_IP_structs.dat")
        try: number_of_IP = open(IP_dat).read()
        except: raise ValueError("Initial pool has not been filled")
        struct_coll = structure_collection.get_collection(self.replica_stoic, 0)
        struct_coll.update_local()
        size_of_common = len(struct_coll.structures)
        size_of_added = size_of_common - int(number_of_IP)

        self.output(' Total size of common pool: %i   Total number of GA-added structures: %i'
                    % (size_of_common, size_of_added))
        cst = '|------------------------------GA CONVERGED-------------------------------|'
        st =  '|-------------------------------------------------------------------------|'
        header = '\n'+st+'\n'+cst+'\n'
        if size_of_added >= self.number_of_GA_structures:
            message = ''
            message += header
            message +='Number of GA additions to pool has reached user-specified number: '
            message += str(size_of_added) + '\n'
            message +='Total size of pool: '
            message += str(size_of_common) +'\n'
            message += 'GAtor ended at: '+str(datetime.datetime.now())  + '\n' + st
            self.output(message)
            end = True
        if size_of_common >= self.number_of_tot_structures:
            message = ''
            message += header
            message +='Total size of pool has at least reached user-specified number: '
            message += str(size_of_common) +'\n'
            message +='Number of GA additions to pool: '
            message += str(size_of_added) +'\n'
            message += 'GAtor ended at: '+str(datetime.datetime.now()) + '\n' + st
            self.output(message)
            end = True
        if convergeTF == True:
            message = ''
            message += header + 'Top '+str(self.top_count)
            message +=' energies havent changed in user-specfied number of iterations: '
            message += str(self.max_it) +'\n'
            message +='Number of GA additions to pool: '
            message += str(size_of_added) +'\n'
            message +='Total size of collection: '
            message += str(size_of_common) +'\n'
            message += 'GAtor ended at: '+str(datetime.datetime.now()) + '\n' + st
            self.output(message)
            end = True
        if os.path.isfile(os.path.join(tmp_dir, "kill.dat")):
            self.output("User has killed GAtor")
            os.remove(os.path.join(tmp_dir, "kill.dat"))
            end = True
        return end

    def check_local_convergence(self, old_prop_list):
        '''
        Check if top N structures havent 
        changed in X iterations
        '''
        new_e_list = np.array([])
        coll_new = self.structure_supercoll.get((self.replica_stoic,0))
        new_prop_list = self.return_property_array(coll_new)
        new_list_top_prop = new_prop_list[:self.top_count]
        old_list_top_prop = old_prop_list[:self.top_count]
        self.convergence_count= self.convergence_count + 1
        for prop_new in new_list_top_prop:
            if prop_new in old_list_top_prop:
                continue
            else:
                self.output("-- Top "+str(self.top_count)+" structures have changed.")
                self.convergence_count= 0
                self.output("-- Local convergence counter reset.")
                self.output("-- Replica's local iteration:  "+ str(self.convergence_count))
                return "not_converged"
                self.output("-- Replica's local iteration:  "+ str(self.convergence_count))
                if self.convergence_count== self.max_it:
                    return "converged"

    def restart_scheme(self,restart_replica):
        self.output("-- Restart of replicas has been requested")
        self.output(str(restart_replica))
        structures_to_add = {}
        struct = False
        if struct == False and restart_replica:
            folder_to_scavenge = activity.request_folder_to_check()
            while struct == False and folder_to_scavenge!= False:
                struct = self.structure_scavenge(os.path.join(tmp_dir,folder_to_scavenge))
                if struct == False:
                    folder_to_scavenge = activity.request_folder_to_check()
        return struct

    def structure_scavenge(self,folder,next_step=False,cleanup=True):
        '''
		This routine takes a folder (directory) that should be an fhi-aims job
        directory and salvages a structure from it if next_step=True, 
        geometry.in.next_step has to be present. Will attempt to read
        struct.json to update other properties that might be lost in restart
        if scavenge failure, returns False
		if cleanup=True and scavenging is successful, removes the folder afterwards
		WARNING: make sure the folder is inactive before calling this function 
        '''
        geometry = os.path.join(folder,"geometry.in")
        geometry_step = os.path.join(folder,"geometry.in.next_step")
        json_dat = os.path.join(folder,"struct.json")

        #Check for geometry.in.next_step
        self.output("--Scavenging folder %s--" % folder)
        if next_step and not os.path.isfile(geometry_step):
            self.output("next_step=True, but no geometry.in.next_step found")
            return False
        if os.path.isfile(geometry_step):
            geostring = read_data(folder,"geometry.in.next_step")
        elif os.path.isfile(geometry):
            geostring = read_data(folder,"geometry.in")
        else:
            message = "Neither geometry.in.next_step nor "
            message +="geometry.in is found in the folder"
            self.output(message)
            return False

        #Build Structure() from found geometry
        struct = Structure()
        try:
            struct.build_geo_whole_atom_format(geostring)
        except:
            message = "Attempt to build_geo_from_atom_file failed."
            self.output(message)
            return False

        #Check for json property file
        if os.path.isfile(json_dat):
            infostring = read_data(folder,"struct.json")
            struct_info = Structure()
            success = False
            try:
                struct_info.loads(infostring)
                success = True
            except:
                message = "struct.json found but can't be read"
                message += "Continuing scavenging..."
                self.output(message)
            if success:
                if len(struct_info.geometry)!=len(struct_info.geometry):
                    message = "struct.json found but the number of "
                    message +="atoms is different from geometry.in\n"
                    self.output(message)
                    self.output("Recommend manually removing folder %s" % folder)
                    return False
                else:
                    for key in struct_info.properties:
                        if (not key in struct.properties) or (struct.properties[key]==None):
                            struct.properties[key]=struct_info.properties[key]
                    self.output("struct.json found and information extracted")
                self.output("Scavenge folder success!")
                fdir = os.path.abspath(os.path.join(folder,os.pardir))

        #Backup scavenged folder
        if activity.bk_folder(fdir,folder[len(fdir)+1:],scavenge_dir,"random"):
            self.output("Successfully backed up scavenged folder")
        else:
            self.output("Failed to back up scavenged folder")

        #Clean scavenged folder
        if cleanup:
            try:
                shutil.rmtree(folder)
                self.output("Folder %s removed" % folder)
            except:
                self.output("Scavenged folder %s clean-up failure" % folder)
                #Probably due to file permission issue
                if os.path.exists(os.path.join(folder,"active.info")):
                    try:
                        os.remove(os.path.join(folder,"active.info"))
                        self.output("Removing active.info instead")
                    except:
                        self.output("active.info remains there!")
                        output.time_log("Folder %s unable to clean-up" % folder,self.replica)
                else:
                    self.output("active.info already removed")
        return struct

    def generate_trial_structure(self):
        sname = "run_settings"
        begin_time = time.time()
        structure_supercoll = {}
        structure_supercoll[(self.replica_stoic, 0)] = structure_collection.get_collection(self.replica_stoic, 0)
        total_attempts = self.ui.get_eval(sname,"failed_generation_attempts")
        self.output("Generating trial structure with %i processes" % self.processes)
        selection_module = my_import(self.ui.get('modules', 'selection_module'), package='selection')
        
        # Select ID's of two parents chosen for selection
        parent_a_ID, parent_b_ID = selection_module.main(structure_supercoll, 
                                                         self.replica_stoic, 
                                                         self.replica)
        
        # Choose probability of crossover, mutation, and symmetric crossover
        rand_cross = np.random.random()
        cross_prob = self.ui.get_eval('crossover', 'crossover_probability')
        try: sym_cross_prob = self.ui.get_eval('crossover', 'symmetric_crossover_probability')
        except: sym_cross_prob = 0.0
        mut_prob = 1.0 - float(cross_prob) - float(sym_cross_prob)
        mutation_list = self.ui.get_list('mutation','specific_mutations')
        mutation_choice = np.random.choice(mutation_list) 

        # Setup structure either in serial or mutliprocessing
        count = 0
        struct = False
        while count < total_attempts and struct == False:
            if self.processes == 1: #Serial
                struct = structure_create_for_multiprocessing((self.replica, self.replica_stoic, 
                                                               parent_a_ID, parent_b_ID, 
                                                               rand_cross, cross_prob, 
                                                               sym_cross_prob, mutation_choice))
            else:
                arglist = [(self.replica+"_"+str(x), self.replica_stoic, parent_a_ID, parent_b_ID,\
                rand_cross, cross_prob, sym_cross_prob, mutation_choice) for x in range (self.processes)]
                results = self.worker_pool.map(structure_create_for_multiprocessing, arglist)

                for i in range (self.processes): #Find a success
                    if results[i]!=False:
                        struct = results[i]
                        break
                if struct != False: #Found a success
                    output.move_to_shared_output(self.replica+"_"+str(i), os.path.join(out_tmp_dir,self.replica+".out"))
                #else:
                    #output.local_message("-- "+str(self.processes)+" attempts have failed to create a new structure")
                for i in range (self.processes):
                    try:
                        os.remove(os.path.join(out_tmp_dir, self.replica+"_"+str(i)+".out"))
                    except OSError:
                        pass
            count += self.processes

        end_time = time.time()
        if count == total_attempts and struct==False:
            self.output("-- Generating structure maxed out on generation attempts: %s" % (total_attempts))
            self.output("-- Selecting new parents --")
            return struct
        self.output("-- New trial structure generated:")
        self.output("-- Number of attempts for structure generation: %s" % (count))
        self.output("-- Time for structure generation: %s seconds" % (end_time-begin_time))
        return struct

	def return_energy_array(self, coll):
		e_list = np.array([])
		for index, structure in coll:
			energy = structure.get_property('energy')
			e_list = np.append(energy,e_list)
		e_list = np.sort(e_list.reshape(len(e_list),1),axis=0)
		return e_list

    def structure_comparison(self, struct, comparison_type):
        '''
        This routine takes a structure, updates self.structure_supercoll, and 
                does comparison on the structure
        '''
        t1 = time.time()
        if comparison_type == "pre_relaxation_comparison":
            self.output("\n---- Pre-relaxation Comparison  ----")
        elif comparison_type == "post_relaxation_comparison":
            self.output("\n---- Post-relaxation Comparison ----")
        structure_collection.update_supercollection(self.structure_supercoll)
        is_acceptable = (self.comparison_module.main(struct, self.structure_supercoll.get((self.replica_stoic, 0)),
                                                                                            self.replica, comparison_type))
        t2 = time.time()
        self.output("-- Time taken to compare structure to collection: %0.3f seconds" % (t2-t1))
        if is_acceptable is False:
            self.output('-- Structure is not acceptable')
            return False  # structure not acceptable start with new selection

    def structure_relax(self,struct):
        '''
        This routines takes a structure and relaxes it
        Returns False if doesn't meet SPE criteria or relaxation fails
        '''
        mkdir_p_clean(self.working_dir) # make relaxation directory in tmp

        sname = "save_structures"
        if self.ui.get_boolean(sname,"pre-evaluation"):
            self.add_structure(struct,"pre-evaluation")

        struct, stat = self.relaxation_module.main(struct)
        if stat=="failed" and self.ui.get_boolean(sname,"evaluation_failed"):
            self.add_structure(struct,"evaluation_failed")
        if stat=="rejected" and self.ui.get_boolean(sname,"evaluation_rejected"):
            self.add_structure(struct,"evaluation_rejected")
        if stat=="failed" or stat=="rejected":
            return False

        #----- Make sure cell is lower triangular -----#
        self.output("-- Ensuring cell is lower triangular and orthogonal")
        nmpc = self.ui.get_eval('unit_cell_settings','num_molecules')
        napm = int(struct.get_n_atoms()/nmpc)
        structure_handling.cell_modification(struct, napm,create_duplicate=False)
        struct=structure_handling.cell_lower_triangular(struct,False)	
        a=struct.get_property('lattice_vector_a')
        b=struct.get_property('lattice_vector_b')
        c=struct.get_property('lattice_vector_c')
        struct.set_property('lattice_vector_a',list(a))
        struct.set_property('lattice_vector_b',list(b))
        struct.set_property('lattice_vector_c',list(c))
        if self.ui.all_geo():
            self.output("Final Structure's geometry:\n" +
            struct.get_geometry_atom_format())
        return struct
	
    def success_message(self):
        time = datetime.datetime.now()
        message = ""
        message += "---- Structure Successfully Added to Common Pool! ----"
        message += "\n-- Time: %s " % (time)
        message += "\n-- Structure's index: %s " % (self.struct_index)
        message += "\n-- Added from replica: %s " % (self.replica)
        self.output(message)

    def end_of_iteration_tasks(self, begin_time, old_list_top_en, restart_count, coll):
        # Remove working directory in /tmp
        rmdir_silence(self.working_dir)

        #Check convergence
        converged = self.check_local_convergence(old_list_top_en)

        #End of Iteration Outputs
        IP_dat = os.path.join(tmp_dir,"num_IP_structs.dat")
        number_of_IP = open(IP_dat).read()
        size_of_common = len(self.structure_supercoll.get((self.replica_stoic, 0)).structures)
        size_of_added = size_of_common - int(number_of_IP)
        self.output(' Total size of common pool: %i   Total number of GA-added structures: %i'
                    % (size_of_common, size_of_added))
        if converged is "not_converged":
            self.output("GA not converged yet.")
            pass
        elif converged is "converged":
            return True	

	def add_structure(self, struct, input_ref):
		structure_collection.add_structure(struct,struct.self.replica_stoic(),"pre-evaluation")
		self.output("--Added-- structure %s to input_ref %s"
		% (struct.struct_id,str(input_ref)))

    def restart(self, message): output.restart_message(message)

def structure_create_for_multiprocessing(args):
    '''
    This is a function for structure creation reading straight from the ui.conf file
    '''
    ui = user_input.get_config()
    nmpc = ui.get_eval('unit_cell_settings','num_molecules')
    replica, stoic, parent_a_id, parent_b_id, rand_cross, cross_prob, sym_cross_prob, mut_choice= args
    output.local_message("\n|----------------------- Structure creation process ----------------------|",replica)
    crossover_module = my_import(ui.get('modules', 'crossover_module'), package='crossover')
    try: alt_crossover_module = my_import(ui.get('modules', 'alt_crossover_module'), package='crossover')
    except: pass
    mutation_module = my_import(ui.get('modules', 'mutation_module'), package='mutation')
    struct_coll = structure_collection.get_collection(stoic, 0)
 
    #---- Get structures of selected parent ID's ----#
    struct_coll.update_local()
    parent_a = struct_coll.get_struct(parent_a_id)
    parent_b = struct_coll.get_struct(parent_b_id)
    structures_to_cross = [parent_a, parent_b]
    if structures_to_cross is False:
        output.local_message('Selection failure',replica)
        return False
    if parent_a == None or parent_b == None:
        output.local_message('Parent ID returned NoneType structure',replica)
        return False

    #----- Crossover -----#
    if rand_cross <= cross_prob:
        output.local_message("\n---- Crossover ----", replica)
        new_struct = crossover_module.main(structures_to_cross, replica)
        if new_struct is False:
            output.local_message("Crossover failure", replica)
            return False
        if ui.all_geo():
            output.local_message("\n-- Parent A's geometry --\n" +
            structures_to_cross[0].get_geometry_atom_format(), replica)
            output.local_message("-- Parent B's geometry --\n" +
            structures_to_cross[1].get_geometry_atom_format(), replica)
            output.local_message("-- Child's geometry --\n" +
            new_struct.get_geometry_atom_format(),replica)
        new_struct.set_property('mutation_type', 'No_mutation')
        for i in range(len(structures_to_cross)):
            par_st = structures_to_cross[i]
            new_struct.set_property('parent_' + str(i), par_st.get_stoic_str()+'/'
            + str(par_st.get_input_ref()) + '/' + str(par_st.get_struct_id()))

    #----- Alternative Crossover -----#
    elif rand_cross > cross_prob and rand_cross <= cross_prob + sym_cross_prob:
        output.local_message("\n---- Alternative Crossover ----", replica)
        new_struct = alt_crossover_module.main(structures_to_cross, replica)
        if new_struct is False:
            output.local_message("Crossover failure", replica)
            return False
        if ui.all_geo():
            output.local_message("\n-- Parent A's geometry --\n" +
            structures_to_cross[0].get_geometry_atom_format(), replica)
            output.local_message("-- Parent B's geometry --\n" +
            structures_to_cross[1].get_geometry_atom_format(), replica)
            output.local_message("-- Child's geometry --\n" +
            new_struct.get_geometry_atom_format(),replica)
        new_struct.set_property('mutation_type', 'No_mutation')
        for i in range(len(structures_to_cross)):
            par_st = structures_to_cross[i]
            new_struct.set_property('parent_' + str(i), par_st.get_stoic_str()+'/'
                + str(par_st.get_input_ref()) + '/' + str(par_st.get_struct_id()))

    #----- Mutation  -----#
    elif rand_cross > cross_prob + sym_cross_prob:
        output.local_message("\n---- Mutation ----",replica)
        choice_struct = random.choice([structures_to_cross[0], structures_to_cross[1]])
        if ui.all_geo():
            output.local_message("Single Parent geometry:\n"
            + choice_struct.get_geometry_atom_format(),replica)
        new_struct = mutation_module.main(choice_struct, replica, mut_choice)
        if new_struct!=False and ui.all_geo():
            output.local_message("Mutated geometry:\n"
            + new_struct.get_geometry_atom_format(),replica)
        if new_struct is False:
            output.local_message('-- Mutation failure',replica)
            return False
        new_struct.set_property('parent_1', choice_struct.get_stoic_str() + '/'
            + str(choice_struct.get_input_ref()) + '/' + str(choice_struct.get_struct_id()))
        if random.random() < 0.0:
            output.local_message("\n---- Second Mutation ----",replica)
            choice_struct = new_struct
            new_struct = mutation_module.main(choice_struct, replica)
            if new_struct!=False and ui.all_geo():
                output.local_message("Mutated geometry:\n"
                + new_struct.get_geometry_atom_format(),replica)
            if new_struct is False:
                output.local_message('-- Mutation failure',replica)
                return False

    #---- Orthogonalization ----#
    if ui.ortho():
        output.local_message("\n---- Checking Cell Orthogonalization ----",replica)
        napm = int(new_struct.get_n_atoms()/nmpc)
        success = \
        structure_handling.cell_modification(new_struct, 
						     napm,
						     create_duplicate=False)
        if success == False and ui.verbose():
            message = "--Niggli reduction of lattice failed"
            message += "\n--Lattice vectors: \n"
            message += "\n".join(map(str,
                new_struct.get_lattice_vectors()))
            message += "\n--Only setting to lower triangular"
            output.local_message(message+"\n",replica)
        if success == False:
            structure_handling.cell_lower_triangular(new_struct,
                False)
        if ui.all_geo():
            output.local_message(new_struct.get_geometry_atom_format(),
					     replica)

	#----- Cell Check -----#
	output.local_message("\n---- Cell Checks ----",replica)
	if not structure_handling.cell_check(new_struct, replica): #unit cell considered not acceptable
		return False

    #----- Assign ID -----#
	output.local_message("---- Assign structure ID ----",replica)
	new_struct.struct_id = misc.get_random_index()
	output.local_message("-- ID assigned: "+new_struct.struct_id,replica)
	return new_struct

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
