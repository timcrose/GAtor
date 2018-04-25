"""                                                                            
If any part of this module is used for a publication please cite:              
                                                                               
F. Curtis, X. Li, T. Rose, A. Vazquez-Mayagoitia, S. Bhattacharya,             
L. M. Ghiringhelli, and N. Marom "GAtor: A First-Principles Genetic            
Algorithm for Molecular Crystal Structure Prediction",                         
J. Chem. Theory Comput., DOI: 10.1021/acs.jctc.7b01152;                        
arXiv 1802.08602 (2018)                                                        
"""

from pymatgen import Lattice as LatticeP
from pymatgen import Structure as StructureP
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as SGA
from pymatgen.symmetry.groups import in_array_list
import pymatgen
from structures import structure, structure_handling
import numpy as np
import copy
import spglib

__author__ = "Farren Curtis, Xiayue Li, and Timothy Rose"
__copyright__ = "Copyright 2018, Carnegie Mellon University and "+\
                "Fritz-Haber-Institut der Max-Planck-Gessellschaft"
__credits__ = ["Farren Curtis", "Xiayue Li", "Timothy Rose",
               "Alvaro Vazquez-Mayagoita", "Saswata Bhattacharya",
               "Luca M. Ghiringhelli", "Noa Marom"]
__license__ = "BSD-3"
__version__ = "1.0"
__maintainer__ = "Timothy Rose"
__email__ = "trose@andrew.cmu.edu"
__url__ = "http://www.noamarom.com"

def identify_space_group(struct,key="space_group"):
	''' 
	Args: Structure() from the GA
	Returns: Structure with pymatgen spacegroup attached a property
	'''
	structp = struct.get_pymatgen_structure()
	SG = return_spacegroup(structp)
	struct.set_property(key,SG)
	return struct

def return_spacegroup(structp):
    #try: #newer pymatgen
    spacegroup = SGA(structp, symprec= 1.0).get_space_group_number()
    #except: #older pymatgen
    #    spacegroup = SGA(structp, symprec= 1.0).get_spacegroup_number()
    return spacegroup


def reduce_by_symmetry(struct,create_duplicate=True):
    '''
    Reduce a strucutre with multiple molecules to single molecules
    '''
    if create_duplicate:
        struct = copy.deepcopy(struct)
    structp = struct.get_pymatgen_structure()
    sga = SGA(structp, symprec=1.0)

	#Generate unqiue index for each element
    elements = set([geo[3] for geo in struct.geometry])
    edict = {} 
    i = 0
    for element in elements:
        i += 1
        edict[element] = i
    cell = structp.lattice.matrix,structp.frac_coords,\
        [0 for x in range(len(structp))]
    info = spglib.get_symmetry_dataset(cell,symprec=sga._symprec)
    if info is None:
        return False
    sg = [] #Store symmetry operation for reconstruction
    for rot, trans in zip(info["rotations"],info["translations"]):
        sg.append([np.ndarray.tolist(rot),np.ndarray.tolist(trans)])

    new_geo = None
    symst = sga.get_symmetrized_structure()
    nstruc = structure.Structure()
    for sites in symst.equivalent_sites:
        site = sites[0] #Select only 1 site among the equivalent
        k = structp.index(site)
        nstruc.build_geo_by_atom(site.coords[0],
        site.coords[1],
        site.coords[2],
        site.specie,
        struct.geometry[k]["spin"],
        struct.geometry[k]["charge"],
        struct.geometry[k]["fixed"])

    #Update the geometry of the old structure
    struct.geometry = copy.deepcopy(nstruc.geometry)
    new_lat = symst.equivalent_sites[0][0].lattice.matrix
    struct.properties["lattice_vector_a"] = new_lat[0]
    struct.properties["lattice_vector_b"] = new_lat[1]
    struct.properties["lattice_vector_c"] = new_lat[2]
    structure_handling.cell_lower_triangular(struct,False)
    struct.properties["symmetry_operations"] = sg
    return struct

def get_orbit(self, p, tol=1e-3):
	"""
	Used to overwrite the get_orbit function in the SpaceGroup class of Pymatgen
	In order to not chop up the molecules
	Returns the orbit for a point.

	Args:
		p: Point as a 3x1 array.
		tol: Tolerance for determining if sites are the same. 1e-5 should
		be sufficient for most purposes. Set to 0 for exact matching
		(and also needed for symbolic orbits).

	Returns:
		([array]) Orbit for point.
	"""
	orbit = []
	for o in self.symmetry_ops:
		pp = o.operate(p)
#		pp = np.mod(np.round(pp, decimals=10), 1)
		#Commented out to not tear up the molecules
		if not in_array_list(orbit, pp, tol=tol):
			orbit.append(pp)
	return orbit

def rebuild_by_symmetry(struct,symmops=None,napm=None,create_duplicate=True):
	'''
	Args:
		struct: Input structure
		symmops: A list of symmetry operations
		napm: Number of atoms per molecule
	Returns:
		Reconstructed structure
	'''
	
	if create_duplicate:
		struct = copy.deepcopy(struct)
	structp = struct.get_pymatgen_structure()
	nstruct = structure.Structure()

	if napm == None and not "NAPM" in struct.properties:
		raise ValueError("Missing napm in both argument "+
		"(napm) and properties (NAPM)")
	if symmops == None and not "symmetry_operations" in struct.properties:
		raise ValueError("Missing symmetry operation "+
		"in both argument (symmops) and properties (symmetry_operations)")

	lats = struct.get_lattice_vectors()
	if symmops == None:
		symmops = []
		for rot, trans in struct.properties["symmetry_operations"]:
			op = pymatgen.core.operations.SymmOp.\
			     from_rotation_and_translation(rot, trans)
			symmops.append(op)
		del struct.properties["symmetry_operations"]

#	for op in symmops:
#		if not is_compatible(lats,op.rotation_matrix):
#			print "Rebuild symmetry operation incompatible"
#			return False
		
	for op in symmops:
		for old_site in structp:
			k = structp.index(old_site)
			#Build new site from operation
			site = (pymatgen.symmetry.analyzer.
				PeriodicSite(old_site.species_and_occu,
				op.operate(old_site.frac_coords),old_site.lattice))

			#Build new geo from pymatgen site
			nstruct.build_geo_by_atom(site.coords[0],
						 site.coords[1],
						 site.coords[2],
						 site.specie,
						 struct.geometry[k]["spin"],
						 struct.geometry[k]["charge"],
						 struct.geometry[k]["fixed"])

	#Update geometry
	#struct.geometry = copy.deepcopy(nstruc.geometry)
	new_lat = structp[0].lattice.matrix
	nstruct.properties["lattice_vector_a"] = new_lat[0]
	nstruct.properties["lattice_vector_b"] = new_lat[1]
	nstruct.properties["lattice_vector_c"] = new_lat[2]
	structure_handling.cell_lower_triangular(nstruct,False)
	structure_handling.move_molecule_in(nstruct,
					    len(nstruct.geometry)/napm,
					    False)
	return nstruct

def is_transformation_matrix(mat,tol=0.01):
	'''
	Determine whether or not a matrix is a rotation matrix
	'''
	det = np.linalg.det(mat)
	if abs(det - 1) > tol and abs(det + 1) > tol:
		return False

	n = np.dot(np.transpose(mat),mat)
	for i in range(len(mat)):
		for j in range(len(mat)):
			if (i==j and abs(n[i][i]-1) > tol)\
			or (i!=j and abs(n[i][j]) > tol):
				return False
	return True

def is_compatible(lat,rot,tol=0.05):
	'''
	Determine whether a rotation operation is compatible with a given lattice
	Lattice should be given as row vectors
	'''
	latt = np.transpose(lat)
	return is_transformation_matrix(np.dot(latt,
					np.dot(rot,np.linalg.inv(latt))))

def are_symmops_compatible(lat,symmops,tol=0.01):
	'''
	Determine whether a set of symmetry operations (with rot and trans)
	is compatible with a given lattice
	'''
	for op in symmops:
		if not is_compatible(lat,op[0]):
			return False
	return True

