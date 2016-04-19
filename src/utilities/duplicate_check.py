from structures.structure import Structure
from structures.structure_collection import StructureCollection

import os
import numpy as np
#--PYMATGEN MODULES--#
from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder as vcf
from pymatgen.analysis.structure_matcher import StructureMatcher,ElementComparator,SpeciesComparator,FrameworkComparator
from pymatgen import Lattice as LatticeP
from pymatgen import Structure as StructureP
from pymatgen.core.structure import SiteCollection as sc
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga




def duplicate_check(struct):	
	structp = struct.get_pymatgen_structure()
	print structp

if __name__ == '__main__':
        duplicate_check(struct, struct_coll)

