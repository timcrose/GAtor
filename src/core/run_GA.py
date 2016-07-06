#!/usr/bin/python

from __future__ import division
import datetime
import os
import sys
import time
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
from utilities import misc, compute_spacegroup_pymatgen
from selection import structure_selection
import copy, shutil, multiprocessing


def main(replica,stoic):
	mkdir_p(tmp_dir)	
	mkdir_p(structure_dir)
	ga = RunGA(replica, stoic)	# run genetic algorithm
	if ga.ui.get_eval('run_settings', 'recover_from_crashes') is not True: # catch crashes
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
		'''Initialization of replica class fields'''
		self.replica = replica
		self.ui = user_input.get_config()
		self.replica_stoic = stoic
		self.working_dir = os.path.join(tmp_dir, str(self.replica))
		self.GA_module_init() # initializes GA modules specified in .conf file
		self.verbose = self.ui.get_eval('run_settings', 'verbose')        
#		self.control_list = self.ui.get_list('FHI-aims', 'control_in_filelist')
		self.singlemutate = self.ui.get_eval('mutation', 'mutation_probability')
		self.doublemutate = self.ui.get_eval('mutation', 'double_mutate_prob')

		self.prop = self.ui.get("run_settings","property_to_optimize")
		self.op_style = self.ui.get("run_settings","optimization_style")

		self.top_count = self.ui.get_eval('run_settings', 'followed_top_structures')
		self.max_it = self.ui.get_eval('run_settings', 'max_iterations')
		self.number_of_structures = int(self.ui.get('run_settings', 'number_of_structures'))
		# Initialize Supercollection
		self.replica_child_count = 0
		self.convergence_count= 0
		self.structure_supercoll = {}
		self.structure_supercoll[(self.replica_stoic, 0)] = structure_collection.get_collection(self.replica_stoic, 0)
		self.structure_supercoll[(self.replica_stoic, 'duplicates')] = structure_collection.get_collection(self.replica_stoic, 'duplicates')

		structure_collection.update_supercollection(self.structure_supercoll)
		self.structure_coll = structure_collection.stored_collections[(self.replica_stoic, 0)]
		
#------------------------- MAIN TASKS PERFORMED BY EVERY REPLICA -------------------------#
	def start(self):
		'''
		Performs main genetic algorithm operations
		Loads necessary modules based on UI at runtime
		''' 

		# Report Replica to Common Output and Update Supercollection
		self.output("\n--Replica %s updating local pool--" %(self.replica))
		structure_collection.update_supercollection(self.structure_supercoll)

		# Intialiaze restarts
		restart_replica = self.ui.get_boolean("parallel_settings","restart_replica")
		restart_count = 0
		convergeTF = None
		while True:
			#----- Beginning of Iteration Tasks -----#
			begin_time = self.beginning_tasks(restart_count)

			#----- Check if GA finished/converged -----#
			end = self.check_finished(convergeTF)
			if end: return

			#----- Restart Scheme -----#
			self.restart_scheme(restart_replica)

			#-----Generate Trial Structure -----#
			struct = self.generate_trial_structure()
			
			#----- Compare Pre-relaxed Structure to Collection -----#
			if self.structure_comparison(struct, "pre_relaxation_comparison") == False:
				convergeTF = False
				rmdir_silence(self.working_dir)
				continue

			#----- Relax Trial Structure -----#
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
			struct = compute_spacegroup_pymatgen.main(struct)
			self.output("Space group %s" % (struct.get_property('space_group')))

			#---- Compute Spacegroup of Relaxed Structure ----#
			struct = compute_spacegroup_pymatgen.main(struct)
			self.output("Space group %s" % (struct.get_property('space_group')))

			#---- Check If Energy is Global Minimum -----#
			ref_label = 0
			coll = self.structure_supercoll.get((self.replica_stoic,ref_label))
			top_prop_list = self.return_property_array(coll)
			self.check_global_optimization(struct,top_prop_list)
			
			#----- Add Structure to collection -----#
			self.add_to_collection(struct, ref_label)

			#----- End of Iteration Data Tasks -----#
			restart_count += 1 
			convergeTF = self.end_of_iteration_tasks(restart_count, top_prop_list, begin_time, coll)

#----------------------- END OF MAIN TASKS PERFORMED BY EVERY REPLICA ----------------------#


#----------------------- FUNCTIONS USED WITHIN MAIN REPLICA TASKS --------------------------#
	def GA_module_init(self):
		'''
		This routine reads in the modules defined in ui.conf
		'''
#		print(self.ui.get('modules','initial_pool_module'))
#		self.initial_pool_module = my_import(self.ui.get('modules', 'initial_pool_module'), package='initial_pool')
		self.selection_module = my_import(self.ui.get('modules', 'selection_module'), package='selection')
		self.crossover_module = my_import(self.ui.get('modules', 'crossover_module'), package='crossover')
		self.mutation_module = my_import(self.ui.get('modules', 'mutation_module'), package='mutation')
		self.relaxation_module = my_import(self.ui.get('modules', 'relaxation_module'), package='relaxation')
		self.comparison_module = my_import(self.ui.get('modules', 'comparison_module'), package='comparison')

#	def initialize_restart(self):
#		restart_count = 0
#		restart_replica = self.ui.get_eval("parallel_settings","restart_replicas")
#		return restart_replica

	def beginning_tasks(self, restart_count):
		output.move_to_shared_output(self.replica)
		self.output('Beginning iteration')
		begin_time = datetime.datetime.now()
		return begin_time

	def check_finished(self, convergeTF):
		end = False
		IP_dat = os.path.join(tmp_dir,"num_IP_structs.dat")
		number_of_IP = open(IP_dat).read()
		
                struct_coll = structure_collection.get_collection(self.replica_stoic, 0)
                struct_coll.update_local()
		size_of_common = len(struct_coll.structures)
                size_of_added = size_of_common - int(number_of_IP)

		self.output("GA added structures: "+ str(size_of_added))
		self.output("Total structures: "+ str(size_of_common))

		if size_of_added >= self.number_of_structures:
			message = ''
			message +=' ~*~*~*~*~* GA CONVERGED *~*~*~*~*~\n'
			message +='Number of additions to pool has reached user-specified number: '
			message += str(size_of_added) +'\n'
			message +='Total size of collection: '
			message += str(size_of_common) +'\n'
			message += 'GA ended at: '+str(datetime.datetime.now())
			self.output(message)
			end = True
		if convergeTF == True:
			message = ''
			message +=' ~*~*~*~*~* GA CONVERGED *~*~*~*~*~\n'
			message +='Top energies havent changed in user-specfied number of iterations: '
			message += str(self.max_en_it) +'\n'
			message +='Total size of collection: '
			message += str(total_structs) 
			self.output(message)
			end = True
		return end

	def restart_scheme(self,restart_replica):	
			structures_to_add = {}	
			struct = False
			if os.path.isdir(self.working_dir): 
				#First check if there is a risk of overwriting a name_sake replica's leftover folder
				struct = self.structure_scavenge_old(self.working_dir)
			if struct == False and restart_replica:
				folder_to_scavenge = activity.request_folder_to_check()
				while struct == False and folder_to_scavenge!= False:
					struct = self.structure_scavenge_old(os.path.join(tmp_dir,folder_to_scavenge))
					if struct == False:
						folder_to_scavenge = activity.request_folder_to_check()

			return struct

	def generate_trial_structure(self):
		struct = False 
		sname = "parallel_settings"
		processes = 1
		if self.ui.has_option(sname,"processes_per_replica"):
			processes = self.ui.get_eval(sname,"processes_per_replica")
			if self.ui.has_option(sname,"processes_per_node"):
				processes = min(self.ui.get_eval(sname,"processes_per_node"), processes)

		
		sname = "run_settings"
		total_attempts = self.ui.get_eval(sname,"failed_generation_attempts")
		count = 0
		begin_time = time.time()
		while count<total_attempts and struct == False:
			if processes == 1: #Serial
				struct = self.structure_create_new()
			else:
				p = multiprocessing.Pool(processes=processes)
				arglist=[(self.replica+"_"+str(x),self.replica_stoic)\
				for x in range (processes)]
				results = p.map(structure_create_for_multiprocessing,\
						arglist)
				for i in range (processes): #Find a success
					if results[i]!=False:
						struct = results[i]
						break
				if struct!=False: #Found a success
					output.move_to_shared_output(self.replica+"_"+str(i),self.replica+".out")
				else:
					output.local_message(str(processes)+" attempts have failed to create a new structure")
				for i in range (processes):
					try:
						os.remove(os.path.join(cwd,self.replica+"_"+str(x)+".out"))
					except OSError:
						pass

			count += processes

		end_time = time.time()
		if count == total_attempts and struct==False:
			raise RuntimeError("Generating structure maxed out on generation attempts.")
		output.local_message("-- Number of attempts for structure generation: "+str(counter))
		output.local_message("-- Time for structure generation: "+str(end_time-begin_time)+" s")
		
		return struct

	def return_energy_array(self, coll):
		e_list = np.array([])
		for index, structure in coll:
			energy = structure.get_property('energy')
			e_list = np.append(energy,e_list)
		e_list = np.sort(e_list.reshape(len(e_list),1),axis=0)
		return e_list

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
		self.output("structure's %s: %f; previous global best: %f" % 
		(self.prop, prop, glob))
		diff = abs(prop-glob)
		if self.op_style=="minimize" and prop<glob:
			message = '*********** NEW GLOBAL MINIMUM FOUND ************' + \
			'\n  old minimum:  ' + str(glob) + \
			'\n  new minimum:  ' + str(prop) + \
			'\n  difference:  ' + str(diff)
		if self.op_style=="maximize" and prop>glob:
			message = '*********** NEW GLOBAL MAXIMUM FOUND ************' + \
			'\n  old maximum:  ' + str(glob) + \
			'\n  new maximum:  ' + str(prop) + \
			'\n  difference:  ' + str(diff)
		self.output(message)

	def check_local_convergence(self, old_prop_list):
		#Check if top N structures havent changed in X iterations

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
				self.output("Top "+str(self.top_count)+" best structures have changed.")
				self.convergence_count= 0
				self.output("Convergence counter reset.")
				self.output("Convergence iteration:  "+ str(self.convergence_count))
				return "not_converged"
		self.output("Convergence iteration:  "+ str(self.convergence_count))
		if self.convergence_count== self.max_it:
			return "converged" 


	def add_to_collection(self, struct, ref_label):	
		prev_struct_index = None #none for now because only using GA for FF not both
		ID = len(StructureCollection(self.replica_stoic, 0).structures) + 1
		#ID = len(self.structure_supercoll.get((self.replica_stoic,0)).structures) + 1
		self.replica_child_count = self.replica_child_count + 1
		struct.set_property('ID', ID) 
		try:
			struct.set_property('prev_struct_id', prev_struct_index)  # tracks cascade sequence
			struct.set_property('ID', ID)	
			struct.set_property('replica', self.replica)	
			struct.set_property('child_counter', self.replica_child_count)		
		except: pass

		struct_index = structure_collection.add_structure(struct, self.replica_stoic, ref_label)
		struct_coll = StructureCollection(self.replica_stoic, ref_label)
		struct_coll.update_local()
		structure_collection.update_supercollection(self.structure_supercoll) #UpdateSupercollection/Database		
                data_tools.write_energy_vs_addition(struct_coll)
                data_tools.write_energy_hierarchy(struct_coll)
		data_tools.write_spe_vs_addition(struct_coll)
		self.output("Structure index: "+str(struct_index))
		return

	def structure_create_new(self):        	  
		'''
		This is the normal process to create a new structure through crossover and mutation
		'''
		#----- Structure Selection -----#
		self.output("-------- Structure Creation --------")
		self.output("--Structure selection--")	
		structures_to_cross = self.selection_module.main(self.structure_supercoll, self.replica_stoic, self.replica)
		if structures_to_cross is False: 
			self.output('Selection failure')
			return False

		#----- Crossover -----#
		self.output("\n--Crossover--")
		new_struct = self.crossover_module.main(structures_to_cross, self.replica)
		if new_struct is False: 
			self.output("Crossover failure")
			return False  
	
		#----- Mutation Execution -----#
		self.output("\n--Mutation--")  
		randnum = np.random.random()	
		randnum2 = np.random.random()
		#Single Parents have to be mutated	
		if new_struct.get_property('crossover_type') == [1,1] or new_struct.get_property('crossover_type') == [2,2] or randnum<self.singlemutate:
			new_struct = self.mutation_module.main(new_struct, self.replica)
			if new_struct!=False and randnum2 < self.doublemutate:
				self.output("--Second Mutation--")
				new_struct = self.mutation_module.main(new_struct, self.replica) 
		else:
			self.output('No mutation applied.')
			new_struct.set_property('mutation_type', 'no_mutation')
		if new_struct is False: 
			self.output('Mutation failure')
			return False

		#-----Structure modification of angles. Checks reasonable structure is created -----#
		self.output("\n--Cell Checks--")	
		structure_handling.cell_modification(new_struct, self.replica,create_duplicate=False)
		if not structure_handling.cell_check(new_struct,self.replica): #unit cell considered not acceptable
			return False
		self.set_parents(structures_to_cross, new_struct)

		self.output("\n--Assign structure ID--")
		new_struct.struct_id = misc.get_random_index()
		self.output("ID assigned: "+new_struct.struct_id)
		return new_struct
	
	def structure_scavenge_old(self,folder,next_step=False,cleanup=True):
		'''
		This routine takes a folder (directory) that should be an fhi-aims job directory and salvages a structure from it
		if next_step=True, geometry.in.next_step has to be present
		Will attempt to read struct.json to update other properties that might be lost in restart
		if scavenge failure, returns False
		if cleanup=True and scavenging is successful, removes the folder afterwards
		WARNING: make sure the folder is inactive before calling this function 
		'''
		self.output("--Scavenging folder %s--" % folder)
		if next_step and not os.path.isfile(os.path.join(folder,"geometry.in.next_step")):
			self.output("next_step=True, but no geometry.in.next_step found")
			return False

		if os.path.isfile(os.path.join(folder,"geometry.in.next_step")):
			geostring = read_data(folder,"geometry.in.next_step")
		elif os.path.isfile(os.path.join(folder,"geometry.in")):
			geostring = read_data(folder,"geometry.in")
		else:
			self.output("Neither geometry.in.next_step nor geometry.in is found in the folder")
			return False
		struct = Structure()
		try:
			struct.build_geo_whole_atom_format(geostring)
		except:
			self.output("Attempt to build_geo_from_atom_file failed. File possibly corrupted")
			return False

		if os.path.isfile(os.path.join(folder,"struct.json")):
			infostring = read_data(folder,"struct.json")
			struct_info = Structure()
			success = False
			try:
				struct_info.loads(infostring)
				success = True
			except:
				self.output("struct.json found but corrupted ; moving on with scavenging")
			if success:
				if len(struct_info.geometry)!=len(struct_info.geometry):
					self.output("struct.json found but the length of its geometry is not the same as from geometry.in")
					self.output("File possibly corrupted")
					self.output("Recommend manually removing folder %s" % folder)
					return False
				else:
					for key in struct_info.properties:
						if (not key in struct.properties) or (struct.properties[key]==None):
							struct.properties[key]=struct_info.properties[key]
					self.output("struct.json found and information extracted")
		
				self.output("Scavenge folder success!")
				fdir = os.path.abspath(os.path.join(folder,os.pardir))
		if activity.bk_folder(fdir,folder[len(fdir)+1:],scavenge_dir,"random"):
			self.output("Successfully backed up scavenged folder")
		else:
			self.output("Failed to back up scavenged folder")
		
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
						output.time_log("WARNING! Folder %s unable to clean-up" % folder,self.replica)
				else:
					self.output("active.info already gone")		
		
		return struct	
			
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
		self.output("Ensuring cell is lower triangular...")
		struct=structure_handling.cell_lower_triangular(struct,False)	
		a=struct.get_property('lattice_vector_a')
		b=struct.get_property('lattice_vector_b')
		c=struct.get_property('lattice_vector_c')		
		struct.set_property('lattice_vector_a',list(a))
		struct.set_property('lattice_vector_b',list(b))
		struct.set_property('lattice_vector_c',list(c))
		self.output("post second check geo: ")
		#self.output(str(struct.get_geometry_atom_format()))
		return struct

	def structure_comparison(self, struct, comparison_type):
		'''
		This routine takes a structure, updates self.structure_supercoll, and does comparison on the structure
		'''
		t1 = time.time()
		if comparison_type == "pre_relaxation_comparison":
			self.output("\n--Pre-relaxation Comparison--")
		elif comparison_type == "post_relaxation_comparison":
			self.output("\n--Post-relaxation Comparison--")
		structure_collection.update_supercollection(self.structure_supercoll)
		is_acceptable = (self.comparison_module.main(struct, self.structure_supercoll.get((self.replica_stoic, 0)), self.replica, comparison_type))
                t2 = time.time()
                self.output("Time taken to compare structure to collection: %0.3f seconds" % (t2-t1))
                if is_acceptable is False:
                        self.output('Newly relaxed structure is not acceptable')
			
                        return False  # structure not acceptable start with new selection

	def end_of_iteration_tasks(self, begin_time, old_list_top_en, restart_count, coll):
		# Remove working directory in /tmp
		rmdir_silence(self.working_dir)

		#Check convergence
		converged = self.check_convergence(old_list_top_en)	
		message = 'Success!: \n  stoichiometry-- '+str(self.replica_stoic)+\
			'\n  cascade-- ' + str("key[1]") + \
			'\n  structure index-- ' + str("struct_index") + \
			'\n  replica child count-- ' + str(self.replica_child_count) + \
			'\n  collection count -- ' + str("ID") + \
			'\n  replica-- ' + str(self.replica)

		self.output(message)

		#write energy hierarchy 
		self.output("Writing hierachy and data files")
		#data_tools.write_energy_hierarchy(self.structure_coll)		   	
	        #data_tools.write_energy_vs_iteration(self.structure_coll)
		#End of Iteration Outputs
		IP_dat = os.path.join(tmp_dir,"num_IP_structs.dat")
		number_of_IP = open(IP_dat).read()
		size_of_common = len(self.structure_supercoll.get((self.replica_stoic, 0)).structures)
		size_of_added = size_of_common - int(number_of_IP)
		self.output('Total size of common pool: '+str(size_of_common))
		self.output('Total number of GA-added structures: '+str(size_of_added))	
		if converged is "not_converged":
			self.output("GA not converged yet.")
			pass	
		elif converged is "converged":
			return True	


#	def check_if_global_minima(self, e_new, min_e):	
#		if e_new < min_e:
#			diff = min_e - e_new
#			message = '*********** NEW GLOBAL MINIMUM FOUND ************' + \
#			'\n  old minima:  ' + str(min_e) + \
#			'\n  new minima:  ' + str(e_new) + \
#			'\n  difference:  ' + str(diff)
#			self.output(message)

#	def get_top_energies(self, coll):
#		list_top_ens = np.array([])
#		list_ens = np.array([])
#		coll_new = self.structure_supercoll.get((self.replica_stoic,0)) 
#		for index, structure in coll_new:
#			energy = structure.get_property('energy')
#			list_ens = np.append(energy,new_e_list)
#		list_top_ens= np.sort(new_e_list.reshape(len(new_e_list),1),axis=0)
#		list_top_ens = new_e_list[:self.top_en_count]
#		return list_top_ens

#	def check_convergence(self, old_list_top_en):
#		#Check if top N energies havent changed in X iterations
#		new_e_list = np.array([])
#		coll_new = self.structure_supercoll.get((self.replica_stoic,0)) 
#		for index, structure in coll_new:
#			energy = structure.get_property('energy')
#			new_e_list = np.append(energy,new_e_list)
#		new_e_list= np.sort(new_e_list.reshape(len(new_e_list),1),axis=0)
#		new_list_top_en = new_e_list[:self.top_en_count]
#		min_e = new_e_list[0][0]
#		max_e = new_e_list[-1][0]
		#self.output("old top energies:    "+ str(old_list_top_en))
		#self.output("new top energies:    "+ str(new_list_top_en))


#		self.convergence_count= self.convergence_count+ 1
#		for en_new in new_list_top_en:
#			if en_new in old_list_top_en:
#				continue
#			else:
#				self.output("Top "+str(self.top_en_count)+" energies have changed.")
#				self.convergence_count= 0
#				self.output("Convergence counter reset.")
#				self.output("Convergence iteration:  "+ str(self.convergence_count))
#				return "not_converged"
#		self.output("Convergence iteration:  "+ str(self.convergence_count))
#		if self.convergence_count== self.max_en_it:
#			return "converged" 

	def avg_fitness(self, min_e, max_e, structure_coll):
		fitness = {}
		for index, struct in structure_coll:
			try: energy = float(struct.get_property('energy'))
			except: pass
			rho = (max_e - energy) / (max_e - min_e)
			if self.ui.get('selection', 'fitness_function') == 'standard':
				fitness[struct] = rho
				sorted_fit = sorted(fitness.iteritems(), key=lambda x:x[1])
				f_sum = 0
				tmp = []
				total = 0
				normalized_fit = []

		# sum and reduce all fitnesses
		for index, fitness in sorted_fit: f_sum += fitness
		for index, fitness in sorted_fit: tmp.append((index, fitness / f_sum))
		for index, t_fitness in tmp:
			normalized_fit.append((index, t_fitness + total))
			total = t_fitness + total
		return normalized_fit


	def set_parents(self, structures_to_cross, struct):
		for i in range(len(structures_to_cross)): 
			par_st = structures_to_cross[i]
			struct.set_property('parent_' + str(i), par_st.get_stoic_str() + '/' \
			+ str(par_st.get_input_ref()) + '/' + str(par_st.get_struct_id()))

	def output(self, message): output.local_message(message, self.replica)

	def add_structure(self,struct,input_ref):
		structure_collection.add_struture(struct,struct.self.replica_stoic(),"pre-evaluation")
		self.output("--Added-- structure %s to input_ref %s"
		% (struct.struct_id,str(input_ref)))
	

	def restart(self, message): output.restart_message(message)

	#--------------------------------------- END OF FUNCTIONS USED WITHIN MAIN REPLICA TASKS ----------------------------------------#




def structure_create_for_multiprocessing(args):        	  
	'''
	This is a function for structure creation reading straight from the ui.conf file
	'''
	ui = user_input.get_config()
	replica, replica_stoic = args
	#----- Structure Selection -----#
	output.local_message("--Beginning normal structure creation process--",replica)
	output.local_message("--Structure selection--", replica)
	selection_module = my_import(ui.get('modules', 'selection_module'), package='selection')
	crossover_module = my_import(ui.get('modules', 'crossover_module'), package='crossover')
	mutation_module = my_import(ui.get('modules', 'mutation_module'), package='mutation')
	structure_supercoll = {}
	structure_supercoll[(stoic, 0)] = structure_collection.get_collection(stoic, 0)
	structures_to_cross = selection_module.main(structure_supercoll, stoic)
	replica = ui.get_replica_name()

	if structures_to_cross is False: 
		output.local_message('Selection failure',replica)
		return False

	#----- Crossover -----#
	output.local_message("\n--Crossover--", replica)
	new_struct = crossover_module.main(structures_to_cross, replica)
	if new_struct is False:
		output.local_message("Crossover failure", replica)
		return False
	
	#----- Mutation Execution -----#
	output.local_message("\n--Mutation--",replica)
	randnum = np.random.random()	
	randnum2 = np.random.random()
	#Single Parents have to be mutated	
	if new_struct.get_property('crossover_type') == [1,1] or new_struct.get_property('crossover_type') == [2,2] or randnum<self.singlemutate:
		new_struct = mutation_module.main(new_struct, replica)
		if new_struct!=False and randnum2 < ui.get_eval('mutation', 'double_mutate_prob'):
			output.local_message("--Second Mutation--",replica)
			new_struct = mutation_module.main(new_struct, replica) 
	else:
		output.local_message('No mutation applied.',replica)
		new_struct.set_property('mutation_type', 'No_mutation')
	if new_struct is False: 
		output.local_message('Mutation failure'.replica)
		return False

	#----- Cell Check -----#
	output.local_message("--Cell Checks--",replica)
	structure_handling.cell_modification(new_struct, replica,create_duplicate=False)
	if not structure_handling.cell_check(new_struct, replica): #unit cell considered not acceptable
		return False

	#----- Set parents ------#
	for i in range(len(structures_to_cross)):  
		par_st = structures_to_cross[i]
		struct.set_property('parent_' + str(i), par_st.get_stoic_str() + '/' \
		+ str(par_st.get_input_ref()) + '/' + str(par_st.get_struct_id()))

	output.local_message("\n--Assign structure ID--",replica)
	new_struct.struct_id = misc.get_random_index()
	output.local_message("ID assigned: "+new_struct.struct_id,replica)
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
