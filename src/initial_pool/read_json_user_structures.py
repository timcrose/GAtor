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
    success = fip.copy_user_structures()
    return success	

class FillInitialPool():
    
    def __init__(self, replica, stoic):
        self.replica = replica
	self.stoic = stoic
        self.ui = user_input.get_config()
	self.FFcommonpool_index = 0
	self.structure_coll0 = StructureCollection(stoic, self.FFcommonpool_index)
	self.ip_count = 0
	self.reject_ip_count = 0
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
            
#	for file in files_to_add:
#	    print "TEST!!!!!"	
#	    print read_data(file)
#	    struct = Structure()
#           print struct.build_geo_from_json_file(file)
	



        for file in files_to_add:	
            # add structure to collection.
            struct = Structure()
            struct.build_geo_from_json_file(file)
#		print struct
	    self.ip_count = self.ip_count +1
	    struct.set_property('u_def_filename', filename)
	    struct.set_property('ID', self.ip_count)	
	    struct.set_property('replica', 'init__pool')
	    structure_collection.add_structure(struct, struct.get_stoic(), 0) #add struct to common FF pool
            self.structure_coll0.update_local()
	    print "IP structure added. ID: ", self.ip_count
#            self.structure_coll0.update_local()
           # data_tools.write_energy_hierarchy(self.structure_coll0) #starts writing energy hierachy file


        return	

if __name__ == '__main__':
    main()
