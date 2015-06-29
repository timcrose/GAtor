'''
@author: farren & Patrick
Fills the user defined pool into common storage
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


def main():
    # Fills Initial Pool (level 0) and returns the number of successful relaxations
	FFcommonpool_index=0
	ip_count=0
	ui=user_input.get_config()
        # setup
        user_structures_dir = os.path.join(cwd, ui.get('initial_pool', 'user_structures_dir'))
        added_user_structures = os.path.join(tmp_dir, 'added_user_structures.dat')
        open(added_user_structures, 'a').close()  # creates if non-existant
        # read files
        try:files = os.listdir(user_structures_dir)
        except:
		raise RuntimeError('IP directory defined but not found; omit -i if IP filling unnecessary')
	if len(files) == 0:
		print "Empty initial pool directory; omit -i if IP filling unnecessary"
		return 0
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
	os.system("chmod 775 "+added_user_structures)

	for file in files_to_add:	
		# add structure to collection.
		struct = Structure()
		struct.build_geo_from_json_file(file)
		ip_count += 1 #WHy is this necessary?
		struct.set_property('u_def_filename', filename)
		struct.set_property('ID', ip_count)	
		struct.set_property('replica', 'init__pool')
		structure_collection.add_structure(struct, struct.get_stoic(), 0) #add struct to common FF pool
		print "IP structure added. ID: ", ip_count
        return	ip_count

if __name__ == '__main__':
    main()
