'''
@author: farren 
'''

import os
from core import user_input, data_tools, output
from core.file_handler import cwd, tmp_dir
from external_libs.filelock import FileLock
from structures import structure_collection, structure_handling
from structures.structure import get_geo_from_file, Structure
from structures.structure_collection import StructureCollection
from pymatgen import Lattice as LatticeP
from pymatgen import Structure as StructureP
from pymatgen.analysis.structure_matcher import StructureMatcher,ElementComparator,SpeciesComparator,FrameworkComparator


def main():
	'''
	Called by core/master.py when -i option is present, requires the user-defined initial pool (or additional user structures)
	to be added into the common 0 storage. An optional duplicate check is performed if requested in ui.conf
	Returns: Total number of initial structures added to the collection.
    	'''
	ui = user_input.get_config()
    	user_structures_dir = os.path.join(cwd, ui.get('initial_pool', 'user_structures_dir'))
    	added_user_structures = os.path.join(tmp_dir, 'added_user_structures.dat')
	num_IP_structures = os.path.join(tmp_dir, 'num_IP_structs.dat')
	files_to_add = file_lock_structures(user_structures_dir, added_user_structures)
	initial_list = convert_to_structures(files_to_add)
    	if ui.get_eval('initial_pool', 'duplicate_check'):
    		print "Checking Initial pool for duplicates"
    		ip_count = return_IP_non_duplicates(initial_list, ui)
		print "Final Initial Pool Count: "+ str(ip_count)
                with open(num_IP_structures,'a') as f:
                        f.write(str(ip_count))
                        f.close()
		return	ip_count
	else:
		ip_count = return_all_user_structures(initial_list)
                print "Final Initial Pool Count: "+ str(ip_count)
                with open(num_IP_structures,'a') as f:
                        f.write(str(ip_count))
                        f.close()
		return ip_count

def file_lock_structures(user_structures_dir, added_user_structures):
	'''
	Args: Path to user defined structures directory, 
        path name to save filepaths of user added structures 
	Returns: List of files to add to collection
    	'''
	open(added_user_structures, 'a').close()  
    	# read files in user-defined structures directory
	try:files = os.listdir(user_structures_dir)
	except:
		raise RuntimeError('IP directory defined but not found; omit -i in run script if IP filling unnecessary')
	if len(files) == 0:
		print "Empty initial pool directory; omit -i in run script if IP filling unnecessary"
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
	return files_to_add

def convert_to_structures(files_to_add):
	'''
	Args: List of files to add to collection
	Returns: List of Structures(), 
	'''
	initial_list = []
	print "Converting all user-input geometries to Structure() instances"
	for file in files_to_add:	
		struct = Structure()
		struct.build_geo_from_json_file(file)
		struct.set_property('file_path', file)
		struct.set_property('replica', 'init_pool')
		message= "Stoic of IP struct: "  +str(struct.get_stoic())
		mod_struct = structure_handling.cell_modification(struct, "init_pool",create_duplicate=False)
		initial_list.append(mod_struct)
	print initial_list
	return initial_list

def set_IP_structure_matcher(ui):
	'''
	Args: The UI for setting tolerances for pymatgen's structure matcher
	Returns: Pymatgen StructureMatcher object
	'''
	L_tol =ui.get_eval('initial_pool', 'ltol')
	S_tol = ui.get_eval('initial_pool', 'stol')
	Angle_tol = ui.get_eval('initial_pool', 'angle_tol')
	Scale = ui.get_eval('initial_pool', 'scale_vol')
	sm = (StructureMatcher(ltol=L_tol, stol=S_tol, angle_tol=Angle_tol, primitive_cell=True, 
                          scale=Scale, attempt_supercell=False, comparator=SpeciesComparator()))
	return sm

def return_all_user_structures(initial_list):
	'''
	Called when no duplicate check is required. Adds all user-defined structures into common (0) storage
	Args: The initial list of Structures() from the user-defined folder
	Returns: The total number of structures added
	'''
	ip_count = 0
	for struct in initial_list:
		structure_collection.add_structure(struct, struct.get_stoic(), 0)
		ip_count += 1 
	return ip_count

    	
def return_IP_non_duplicates(initial_list, ui):
	'''
	Called when duplicate check is required before adding user-defined structures into common (0) storage
	Args: The initial list of Structures() from the user-defined folder
	Returns: The total number of non-duplicate structures added
	'''
	check_count = 0
	ip_count = 0
	struct_list = []
	sm = set_IP_structure_matcher(ui)
        structure_supercoll = {}
	for struct in initial_list:
		stoic = struct.get_stoic()
		structure_supercoll[(stoic, 0)] = StructureCollection(stoic, 0)
		frac_data = struct.get_frac_data()
		check_count += 1
		structp = get_pymatgen_structure(frac_data)
		if len(struct_list) == 0:
			struct_list.append(struct)
                        struct.set_property('ID',0)	
                        structure_collection.add_structure(struct, stoic, 0)
			print "Added First Structure"
		elif len(struct_list) > 0:
			TF_list = []
			fitTF= False		
			print "Checking Next Structure"
			for comp_struct in struct_list:	
				comp_frac_data = comp_struct.get_frac_data()
        		       	comp_structp = get_pymatgen_structure(comp_frac_data)
				fitTF = sm.fit(structp,comp_structp)
				TF_list.append(fitTF)
				#print TF_list
			if True not in TF_list:
				struct_list.append(struct)
				struct.set_property('ID',0)		
				structure_collection.add_structure(struct, stoic, 0)	
				#print "Found Non-Duplicate!"
				print "Total Non-Duplicates: "+str(len(struct_list))
				print "Total Checked: "+str(check_count) 
				print "Total Duplicates Found: " + str(check_count - len(struct_list))
                #check_count +=1
	ip_count = str(len(struct_list))
        structure_collection.update_supercollection(structure_supercoll)
	data_tools.write_energy_hierarchy(StructureCollection(stoic, 0))
	return ip_count


def get_pymatgen_structure(frac_data):
	'''
	Args: Geometric data from GAtor's Structure() object
	Returns: A pymatgen StructureP() object with the same geometric properties
	'''
	coords = frac_data[0] # frac coordinates
	atoms = frac_data[1] # site labels
	lattice = LatticeP.from_parameters(a=frac_data[2], b=frac_data[3], c=frac_data[4], alpha=frac_data[5],beta=frac_data[6], gamma=frac_data[7])
	structp = StructureP(lattice, atoms, coords)
	return structp

if __name__ == '__main__':
    main()
