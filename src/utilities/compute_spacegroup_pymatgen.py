'''
Computes Spacegroup of Structure()

Created by Farren Curtis on July5 5th, 2016
'''
from pymatgen import Lattice as LatticeP
from pymatgen import Structure as StructureP
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as SGA

def main(struct):
	''' 
        Args: Structure() from the GA
        Returns: Structure with pymatgen spacegroup attached a property
	'''
	structp = get_pymatgen_structure(struct.get_frac_data())
	SG = str(return_spacegroup(structp))[:-11]
	#print SG
	struct.set_property('space_group', SG)
	return struct

def get_pymatgen_structure(frac_data):
        '''
        Args: self, Geometric data from GAtor's Structure() object
        Returns: A pymatgen StructureP() object with the same geometric properties
        '''
        coords = frac_data[0] # frac coordinates
        atoms = frac_data[1] # site labels
        lattice = (LatticeP.from_parameters(a=frac_data[2], b=frac_data[3], c=frac_data[4],
                                  alpha=frac_data[5],beta=frac_data[6], gamma=frac_data[7]))
        structp = StructureP(lattice, atoms, coords)
        return structp

def return_spacegroup(structp):
 	return SGA(structp, symprec= 0.7).get_spacegroup()

if __name__ == "__main__":
        main()