'''
Computes the spacegroup for the given structure
using the python pacakge pymatgen

Created by Farren Curtis on July 5th, 2016
'''
from pymatgen import Lattice as LatticeP
from pymatgen import Structure as StructureP
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as SGA

def main(struct):
    ''' 
    Args: Structure() from the GA
    Returns: Structure with approx spacegroup attached a property
    '''
    structp = get_pymatgen_structure(struct.get_frac_data())
    SG = return_spacegroup(structp)
    struct.set_property("space_group",SG)
    return struct

def get_pymatgen_structure(frac_data):
    '''
    Args: self, Geometric data from GAtor's Structure() object
    Returns: A pymatgen StructureP() object with the same 
    geometric properties
    '''
    coords = frac_data[0] # frac coordinates
    atoms = frac_data[1] # site labels
    lattice = (LatticeP.from_parameters(a=frac_data[2], 
                                        b=frac_data[3], 
                                        c=frac_data[4], 
                                        alpha=frac_data[5],
                                        beta=frac_data[6], 
                                        gamma=frac_data[7]))
    structp = StructureP(lattice, atoms, coords)
    return structp

def return_spacegroup(structp):
    try:
        SG = SGA(structp, symprec= 1.0).get_spacegroup_number()
    except: #different for newer versions of pymatgen
        SG = SGA(structp, symprec= 1.0).get_space_group_number()
    return SG

if __name__ == "__main__":
    main()
