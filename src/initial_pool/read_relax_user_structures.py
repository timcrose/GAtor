'''
@author: farren
'''
import os
import shutil
import time
import numpy
from core import user_input, data_tools
from core.file_handler import cwd, set_progress, my_import, tmp_dir, read_data
from external_libs.filelock import FileLock
from structures import structure_collection
from structures.structure import get_geo_from_file, Structure
from structures.structure_collection import StructureCollection


def main(replica, stoic):
    # Fills Initial Pool (level 0) and returns the number of successful relaxations
    fip = FillInitialPool(replica, stoic)    
    num_success = fip.copy_user_structures()
    return num_success	

class FillInitialPool():
    
    def __init__(self, replica, stoic):
        self.replica = replica
	self.stoic = stoic
        self.ui = user_input.get_config()
	self.FFcommonpool_index = 0
	self.structure_coll0 = StructureCollection(stoic, self.FFcommonpool_index)
	self.ip_count = 0
	self.reject_ip_count = 0
        self.initial_pool_relaxation_module = my_import(self.ui.get('modules', 'initial_pool_relaxation_module'), package='relaxation')
        self.comparison_module = my_import(self.ui.get('modules', 'initial_pool_comparison_module'), package='comparison')
        self.num_initial_pool_structures = int(self.ui.get('initial_pool', 'num_initial_pool_structures'))
        self.control_in = read_data(os.path.join(cwd, self.ui.get('control', 'control_in_directory')),
                                    self.ui.get('control', 'initial_pool'))
        self.verbose = self.ui.get_eval('run_settings', 'verbose')
    	self.working_dir = os.path.join(tmp_dir, str(self.replica))   
	self.control_list = self.ui.get_list('control', 'control_in_filelist') 
    
    def copy_user_structures(self):
        '''
        Adds structures defined by the users to the shared storage.
        '''
        # setup
        user_structures_dir = os.path.join(cwd, self.ui.get('initial_pool', 'user_structures_dir'))
        added_user_structures = os.path.join(tmp_dir, 'added_user_structures.dat')
        open(added_user_structures, 'a').close()  # creates if non-existant
        # read files
        try:files = os.listdir(user_structures_dir)
        except Exception: return  # no directory
        if len(files) == 0:
            return
        # copy files and add name to list of copied files to avoid duplicates
        files_to_add = []
        with FileLock(added_user_structures):
            ffile = open(added_user_structures, 'r')
            file_list = ffile.read()
            ffile.close()
            ffile = open(added_user_structures, 'a')
            for item in files:
                filename = os.path.join(user_structures_dir, item)
                if not filename in file_list:
                    ffile.write(filename + '\n')
                    files_to_add.append(filename)
            ffile.close()
            
        for file in files_to_add:	
            # Tests passed. add structure to collection.
            try:
                struct = Structure()
                struct.build_geo_from_atom_file(file)
            except: continue  # could not build structure from geometry
            # keeps track of which user defined structures have been added
            struct.set_property('to_be_relaxed', True)
            struct.set_property('u_def_filename', filename)

		#####RELAX structures before adding to level 0 ***NEW****** (add to own def later)	
            print "--Considering new IP structure--"	
	    print "--Relaxation--"
	    input_string = read_data(os.path.join(cwd, self.ui.get('control', 'control_in_directory')),self.control_list[0]) 
	    relax_struct = self.initial_pool_relaxation_module.main(struct, self.working_dir, input_string, self.replica)
            if relax_struct is False:
            	if self.verbose: self.output('Relaxation failure')	
	    else:
		print "--Comparison--"
	 	if self.comparison_module.main(relax_struct, self.structure_coll0, self.replica) is True:   	
			self.ip_count = self.ip_count +1
			print "Success! IP Number:  ", self.ip_count
			relax_struct.set_property('ID', self.ip_count) 	
			structure_collection.add_structure(relax_struct, relax_struct.get_stoic(), 0) #add struct to common FF pool
			self.structure_coll0.update_local() 
			data_tools.write_energy_hierarchy(self.structure_coll0) #starts writing energy hierachy file
		else:
			self.reject_ip_count = self.reject_ip_count +1
#	return self.ip_count  #returns total number of relaxed structures added to common pool- should redo this a better way later
	return

 
if __name__ == '__main__':
    main()
