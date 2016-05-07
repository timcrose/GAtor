'''
@author: farren & Patrick
Fills the user defined pool into common storage
'''
import os
import shutil
import time
import numpy
from core import user_input, data_tools, output
from core.file_handler import cwd, set_progress, my_import, tmp_dir, read_data
from external_libs.filelock import FileLock
from structures import structure_collection,structure_handling
from structures.structure import get_geo_from_file, Structure
from structures.structure_collection import StructureCollection
#from utilities.duplicate_check import *
from pymatgen import Lattice as LatticeP
from pymatgen import Structure as StructureP
from pymatgen.analysis.structure_matcher import StructureMatcher,ElementComparator,SpeciesComparator,FrameworkComparator
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
	initial_list = []
	sm = StructureMatcher(ltol=.1, stol=.2, angle_tol=5, primitive_cell=True, scale=True, attempt_supercell=False, comparator=SpeciesComparator())
	for file in files_to_add:	
		# add structure to collection.
		struct = Structure()
		struct.build_geo_from_json_file(file)
		ip_count += 1 #WHy is this necessary?
		struct.set_property('file_path', file)
		struct.set_property('ID', ip_count)	
		struct.set_property('replica', 'init_pool')
	        message= "Stoic of IP struct: "  +str(struct.get_stoic())
		mod_struct = structure_handling.cell_modification(struct, "init_pool",create_duplicate=False)
		initial_list.append(mod_struct)
	print "Done with converting files to Structure() format"

        print "Checking for duplicates"
        struct_list = []
	check_count = 0
	for struct in initial_list:
		check_count = check_count + 1
		frac_data = struct.get_frac_data()
		structp = get_pymatgen_structure(frac_data)
                if len(struct_list) == 0:
                        struct_list.append(struct)	
			print "Added first structure"
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
				structure_collection.add_structure(struct, struct.get_stoic(), 0)	
				print "Found Non-Duplicate!"
				print "Length of Non-Duplicates: "+str(len(struct_list))
				print "Total Checked: "+str(check_count) 
				print "Duplicate Count: " + str(check_count - len(struct_list))
	print "Final size of added collection: "+str(len(struct_list))
        return	ip_count


def get_pymatgen_structure(frac_data):
        '''
        Inputs: A np.ndarry structure with standard "geometry" format
        Outputs: A pymatgen core structure object with basic geometric properties
        '''
#        frac_data = self.get_frac_data()
        coords = frac_data[0] # frac coordinates
        atoms = frac_data[1] # site labels
        lattice = LatticeP.from_parameters(a=frac_data[2], b=frac_data[3], c=frac_data[4], alpha=frac_data[5],beta=frac_data[6], gamma=frac_data[7])
        structp = StructureP(lattice, atoms, coords)
        return structp

if __name__ == '__main__':
    main()
