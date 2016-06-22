from core import output
from structures.structure import Structure
from structures.structure_collection import StructureCollection

import os
import numpy as np
#--PYMATGEN MODULES--#
from pymatgen.analysis.structure_matcher import StructureMatcher,ElementComparator,SpeciesComparator,FrameworkComparator
from pymatgen import Lattice as LatticeP
from pymatgen import Structure as StructureP
from pymatgen.core.structure import SiteCollection as sc
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga




def duplicate_check(structp, structp_coll, replica, ip_count):
	print structp_coll
	if len(structp_coll)!=0:
                sm = (StructureMatcher(ltol=.3, stol=.4, angle_tol=5, primitive_cell=True, 
                      scale=False, attempt_supercell=False, comparator=SpeciesComparator())
                for current_structp in structp_coll:
			print "current", current_structp
                	rms = sm.get_rms_dist(structp,structp_coll[current_structp])
                        fitTF = sm.fit(structp,current_structp)
			if fitTF == False:
				output.local_message("no dup",replica)
			if fitTF == True:
				output.local_message("dupfoun!",replica)
				

			return fitTF 

if __name__ == '__main__':
        duplicate_check(struct, struct_coll)

