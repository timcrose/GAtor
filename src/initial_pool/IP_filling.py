'''
@author: farren 
'''
from __future__ import division
import os
import itertools
import numpy as np
from math import ceil
from core import user_input, data_tools, output
from core.file_handler import cwd, tmp_dir, my_import
from external_libs.filelock import FileLock
from structures import structure_collection, structure_handling
from structures.structure import get_geo_from_file, Structure
from structures.structure_collection import StructureCollection
from utilities import compute_spacegroup_pymatgen
from pymatgen import Lattice as LatticeP
from pymatgen import Structure as StructureP
from pymatgen.analysis.structure_matcher import StructureMatcher,ElementComparator,SpeciesComparator,FrameworkComparator


def main():
    '''
    Called by GAtor_master.py when fill_initial_pool = TRUE
    Adds the user-defined initial pool (or additional user structures)
    into the common 0 storage. 
    An optional duplicate check is performed if requested in ui.conf
    Duplicate check uses RDF differences and the Pymatgen Structure Comparer
    Returns: Total number of initial structures added to the collection.
    '''
    ui = user_input.get_config()
    replica = "init_pool"
    user_structures_dir = os.path.join(cwd, ui.get('initial_pool', 'user_structures_dir'))
    added_user_structures = os.path.join(tmp_dir, 'added_user_structures.dat')
    num_IP_structures = os.path.join(tmp_dir, 'num_IP_structs.dat')
    energy_name = ui.get('initial_pool', 'stored_energy_name')
    files_to_add = file_lock_structures(user_structures_dir, added_user_structures)
    st = ' ------------------------------------------------------------------------'
    gst = '|                       GAtor filling initial pool                       |'
    output.local_message(st, replica)
    output.local_message(gst, replica)
    output.local_message(st, replica)
    output.local_message("Initial pool: %s" % (user_structures_dir), replica)
    initial_list = convert_to_structures(files_to_add, energy_name)
    if len(initial_list) == 0:
        output.time_log("Initial pool already filled")
        output.local_message("Initial pool already filled",replica)
        return 0
    if ui.get_boolean('initial_pool', 'duplicate_check'):
        output.time_log("Checking initial pool for duplicates")
        output.local_message("---- Checking initial pool for duplicates ----",replica)
        ip_count = return_non_duplicates(initial_list, replica, ui)
        output.time_log("Final initial pool count: %i" % ip_count)
        output.local_message("Final initial pool count: %i" % ip_count,replica)
    else:
		ip_count = return_all_user_structures(initial_list, replica, ui)
		output.time_log("Final initial pool count: %i" % ip_count)
		output.local_message("Final initial pool count: %i" % ip_count,replica)
    output.move_to_shared_output(replica="init_pool")
    if ip_count!=0 and ip_count!=None:
        with open(num_IP_structures,'w') as f:
            f.write(str(ip_count))
            f.close()
        return ip_count
    else:
        return 

def file_lock_structures(user_structures_dir, added_user_structures):
	'''
    Args: Path to user defined structures directory, 
    path name to save filepaths of user added structures 
    Returns: List of files to add to collection
    '''
	open(added_user_structures, 'a').close()  
	try:files = os.listdir(user_structures_dir)
	except:
		raise RuntimeError('Initial pool directory defined but not found')
	if len(files) == 0:
		print "Empty initial pool directory"
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
	return files_to_add

def convert_to_structures(files_to_add, energy_name="energy"):
    '''
    Args: List of files to add to collection
    Returns: List of Structures(), 
    '''
    initial_list = []
    ui = user_input.get_config()
    for file in files_to_add:	
        struct = Structure()
        struct.build_geo_from_json_file(file)
        struct.set_property('file_path', file)
        struct.set_property('replica', 'init_pool') 
        en = struct.get_property(energy_name)
        if en is not None:
            struct.set_property('energy', en)
            if ui.ortho():
                napm = int(struct.get_n_atoms()/ui.get_nmpc())
                struct = structure_handling.cell_modification(struct,
							     napm,
							     create_duplicate=True)
            initial_list.append(struct)
        else:
            raise ValueError('Energy name %s not found for  %s' % \
                             (en, file))
    return initial_list

def set_IP_structure_matcher(ui):
	'''
	Args: The UI for setting tolerances for pymatgen's structure matcher
	Returns: Pymatgen StructureMatcher object
	'''
	L_tol =ui.get_eval('initial_pool', 'ltol')
	S_tol = ui.get_eval('initial_pool', 'stol')
	Angle_tol = ui.get_eval('initial_pool', 'angle_tol')
	Scale = ui.get_boolean('initial_pool', 'scale_vol')
	sm = (StructureMatcher(ltol=L_tol, stol=S_tol, angle_tol=Angle_tol, primitive_cell=True, 
                          scale=Scale, attempt_supercell=False, comparator=SpeciesComparator()))
	return sm

def return_all_user_structures(initial_list, replica, ui):
    '''
	Called when no duplicate check is required. Adds all user-defined structures into common (0) storage
	Args: The initial list of Structures() from the user-defined folder
	Returns: The total number of structures added
	'''
    ip_count = 0
    structure_supercoll = {}

    if ui.get_boolean("clustering","cluster_pool"):
        message = "Clustering is requested."  
        output.local_message(message,replica)
    for struct in initial_list:
        stoic = struct.get_stoic()
        struct.set_property('ID',0)
        struct = compute_spacegroup_pymatgen.main(struct)
        struct.set_property('crossover_type', '')
        struct.set_property('mutation_type', '')
        if ui.get_boolean("clustering","cluster_pool"):
            clustering_mod = my_import(ui.get('modules','clustering_module'), package='clustering')
            AFV = clustering_mod.AssignFeatureVector(struct)
            struct = AFV.compute_feature_vector() 
        structure_collection.add_structure(struct, stoic, 0)
        ip_count += 1
    return ip_count


#def compute_feature_vector(struct, ui):
#    feature_vector = ui.get("clustering", "feature_vector")
#    if struct.get_property('feature_vector') is not None:
#        return struct
#    if feature_vector == "RDF_vector":
#        struct = structure_handling.compute_RDF_vector(struct)
#        return struct
#    elif feature_vector == "PCA_RDF_vector":
#        struct = structure_handling.compute_RDF_vector(struct)
#        return struct
#    elif feature_vector == "Lat_vol_vector":
#        a = np.linalg.norm(struct.get_property('lattice_vector_a'))
#        b = np.linalg.norm(struct.get_property('lattice_vector_b'))
#        c = np.linalg.norm(struct.get_property('lattice_vector_c'))
#        vol = struct.get_unit_cell_volume()
#        vol = np.cbrt(vol)
#        lat_vol = [a/vol, b/vol, c/vol]
#        struct.set_property(feature_vector, lat_vol)
#        return struct
#    elif feature_vector == "RCD_vector":  
#        RCD_struct = rcd_vector_calculation(struct)
#        RCD_vector = RCD_struct.get_property(self.feature_type)
#        struct.set_property(feature_vector, RCD_vector)
#    else:
#        message = "Clustering for %s is not availble" % (feature_vector)
#        raise RuntimeError(message)
        


def return_non_duplicates(initial_list, replica, ui):
    '''
    Called when duplicate check is required. 
    Uses RDF and pymatgen similarity comparison
    Args: The initial list of Structures() from the user-defined folder
    Returns: The total number of structures added
    '''
    ui = user_input.get_config()
    remove_list = []
    structure_supercoll = {}
    dup_pairs = return_duplicate_pairs(initial_list, ui, replica)
    if ui.get_boolean("clustering","cluster_pool"):
        message = "Clustering is requested."
        output.local_message(message,replica)
    for path, path_dup in dup_pairs:
        remove_list.append(path_dup)
    for struct in initial_list:
        stoic = struct.get_stoic()
        struct_fp = struct.get_property('file_path')
        if struct_fp in remove_list:
            continue
        struct.set_property('ID',0)
        struct = compute_spacegroup_pymatgen.main(struct)
        struct.set_property('crossover_type', '')
        struct.set_property('mutation_type', '')
        if ui.get_boolean("clustering","cluster_pool"):
           struct = compute_feature_vector(struct, ui)        
        structure_collection.add_structure(struct, stoic, 0)
    if len(initial_list)!=0:
        struct_coll = StructureCollection(stoic, 0)
        structure_collection.update_supercollection(structure_supercoll)
        if ui.verbose(): 
            output.local_message("Total duplicate pairs found: %d" % (len(remove_list)),replica)
            output.local_message("Total checked: %d" % (len(initial_list)),replica)
            output.local_message("Total unique structures added: %d" % (len(struct_coll.structures)),replica)
        return len(struct_coll.structures)

def return_duplicate_pairs(initial_list, ui, replica):
    dup_pairs = []
    ip_dup_output = open(os.path.join(tmp_dir, "IP_duplicates.dat"),'w')
    vector_name = ui.get_eval('initial_pool', 'vector_for_comparison') 
    check_if_vec = initial_list[0].get_property(vector_name)
    if check_if_vec is not None:
        for struct, structc in itertools.combinations(initial_list, 2):
            rdf_tol = compute_rdf_diff(struct, structc, ui)
            if rdf_tol < ui.get_eval('initial_pool', 'vector_cosdiff_threshold'):
                fit = compute_pymatgen_fit(struct, structc, ui)
                if fit:
                    struct_fp = struct.get_property('file_path')
                    structc_fp = structc.get_property('file_path')
                    if ui.verbose():
                        output.local_message("Found duplicate pair", replica)
                        output.local_message("-- %s" % struct_fp, replica)
                        output.local_message("-- %s" % structc_fp, replica)
                    dup_pairs.append((struct_fp, structc_fp))
    else: dup_pairs = return_duplicate_pymat_pairs(initial_list, ui)
    for pair in dup_pairs:
        ip_dup_output.write('\t'.join(str(s) for s in pair) + '\n')
    return dup_pairs

def compute_rdf_diff(struct, structc, ui):
    vector_name = ui.get_eval('initial_pool', 'vector_for_comparison') 
    rd = np.asarray(struct.get_property(vector_name))
    rdc = np.asarray(structc.get_property(vector_name))
    rdf_tol = 1-(np.dot(rd,rdc))/(np.linalg.norm(rd)*np.linalg.norm(rdc))
    return rdf_tol

def compute_pymatgen_fit(struct, structc, ui):
    sm = set_IP_structure_matcher(ui)
    structp = get_pymatgen_structure(struct.get_frac_data())
    structpc = get_pymatgen_structure(structc.get_frac_data())
    fit = sm.fit(structp, structpc)
    return fit
	
def return_duplicate_pymat_pairs(initial_list, ui):
    dup_pairs = []
    sm = set_IP_structure_matcher(ui)
    replica = ui.get_replica_name() 
    for struct, structc in itertools.combinations(initial_list, 2):
        structp = get_pymatgen_structure(struct.get_frac_data())
        structpc = get_pymatgen_structure(structc.get_frac_data())
        fit = sm.fit(structp, structpc)
        if fit:
            struct_fp = struct.get_property('file_path')
            structc_fp = structc.get_property('file_path')
            output.local_message("Found duplicate pair", replica)
            output.local_message("-- %s" % struct_fp, replica)
            output.local_message("-- %s" % structc_fp, replica)
            dup_pairs.append((struct_fp, structc_fp))
    return dup_pairs

def get_pymatgen_structure(frac_data):
	'''
	Args: Geometric data from GAtor's Structure() object
	Returns: A pymatgen StructureP() object with the same geometric properties
	'''
	coords = frac_data[0] # frac coordinates
	atoms = frac_data[1] # site labels
	lattice = (LatticeP.from_parameters(a=frac_data[2], b=frac_data[3], c=frac_data[4], 
                                alpha=frac_data[5],beta=frac_data[6], gamma=frac_data[7]))
	structp = StructureP(lattice, atoms, coords)
	return structp


if __name__ == '__main__':
    main()
