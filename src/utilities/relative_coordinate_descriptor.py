"""                                                                            
If any part of this module is used for a publication please cite:              
                                                                               
F. Curtis, X. Li, T. Rose, A. Vazquez-Mayagoitia, S. Bhattacharya,             
L. M. Ghiringhelli, and N. Marom "GAtor: A First-Principles Genetic            
Algorithm for Molecular Crystal Structure Prediction",                         
J. Chem. Theory Comput., DOI: 10.1021/acs.jctc.7b01152;                        
arXiv 1802.08602 (2018)                                                        
"""

import os
import sys
sys.path.append(os.path.join(os.getcwd(),".."))
import multiprocessing
from structures.structure_handling import *
from structures.structure import Structure
import numpy as np
from copy import deepcopy
from utilities import misc

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

def _flatten_rcd_vector(v):
    return [entry for mol in v for rel in mol for entry in rel]

def generate_relative_coordinate_descriptor(struct,nmpc,napm,axes,close_picks=8,key="RCD_vector"):
    '''
    Generates and stores the rcd for the given structure
    '''
    ref_struct = create_ref_struct(struct, nmpc, napm, close_picks)
    axes_list = [_calculate_molecule_axes(ref_struct.geometry[x * napm:(x + 1) * napm],axes) 
                 for x in range(close_picks + 1)]
    COM = [cm_calculation(ref_struct,range(x*napm,(x+1)*napm)) 
                 for x in range(close_picks + 1)]
    diff = [np.subtract(COM[x], COM[0]) 
                 for x in range(1, close_picks+1)]
    axes_i = np.linalg.inv(axes_list[0].T)
    diff_t = [list(np.dot(axes_i, x)) 
                 for x in diff] #Relative coordinate in axes basis
    orien = [[float(np.dot(axes_list[x][y], axes_list[0][y])) 
                 for y in range(3)] 
                 for x in range(1, close_picks + 1)]
    vector = [[diff_t[x], orien[x]] for x in range(close_picks)]
    struct.properties[key] = vector
    return struct

def _generate_relative_coordinate_descriptor_multiprocess_wrapper(arglist):
    return generate_relative_coordinate_descriptor(*arglist)

def create_ref_struct(struct, nmpc, napm, close_picks=8):
    '''
    Using the 1st molecule in the structure as reference
    Select and return a structure with the molecules that are closest to it

    struct: input structure
    nmpc: number of molecupe per cell
    napm: number of atom per cell
    close_picks: number of closest molecules to select

    '''
    ref_struct = cell_modification(struct, napm)

    COM = [cm_calculation(ref_struct,
                          range(x * napm, (x + 1) * napm))
           for x in range(nmpc)]

    lat_mat = np.transpose(ref_struct.get_lattice_vectors())
    COMt = [(k, x, y, z, 
            np.linalg.norm(np.subtract(np.add(COM[k], np.dot(lat_mat, [x,y,z])),
                                       COM[0]))) 
            for z in range(-2,3)
            for y in range(-2,3)
            for x in range(-2,3)
            for k in range(nmpc)] #Calculate the distance between the COM's
    COMt.sort(key=lambda com: com[4])

    all_geo = deepcopy(ref_struct.geometry)
    ref_struct.geometry = np.delete(ref_struct.geometry,
                                    range(napm,nmpc*napm))

    for i in range(1, close_picks+1): #Skip the 1st original molecule
        k, x, y, z, dist = COMt[i]
        ref_struct.geometry = np.concatenate((ref_struct.geometry,
                                              all_geo[k * napm:(k + 1) * napm]))

        mole_translation(ref_struct, i, napm, frac=[x,y,z], create_duplicate=False)

#    print ref_struct.get_geometry_atom_format()i
    return ref_struct

def _calculate_molecule_axes(geo,axes):
    '''
    Calculates the orientation axes of the molecule
    geo: geometry slice of the molecule
    axes: a list of 2 tuples, each tuple include the indeces of two atoms
    The first 2 axes will be established as the 2 vectors formed by the pairs of atoms
    The third will be their cross product
    '''
    a1 = [geo[axes[0][0]][x]-geo[axes[0][1]][x] for x in range(3)]
    a2 = [geo[axes[1][0]][x]-geo[axes[1][1]][x] for x in range(3)]
    a3 = list(np.cross(a1,a2))
    X = np.array([a1,a2,a3])
    X = _gram_schmidt(X)
    return X


def _gram_schmidt(X, row_vecs=True, norm=True):
    if not row_vecs:
        X = X.T
    Y = X[0:1,:].copy()
    for i in range(1, X.shape[0]):
        proj = np.diag((X[i,:].dot(Y.T)/np.linalg.norm(Y,axis=1)**2).flat).dot(Y)
        Y = np.vstack((Y, X[i,:] - proj.sum(0)))
    if norm:
        Y = np.diag(1/np.linalg.norm(Y,axis=1)).dot(Y)
    if row_vecs:
        return Y
    else:
        return Y.T

def rcd_difference_calculation(inst):
    '''
    Main module of rcd difference
    '''
    sname = "rcd_difference_calculation"
    key = inst.get_with_default(sname,"stored_property_key","RCD_vector")
    ratio = inst.get_with_default(sname,"contribution_ratio",1,eval=True)
    pairs = inst.get_with_default(sname,"select_pairs",4,eval=True)
    dis_en = inst.get_boolean(sname,"disable_enantiomer")
    diff_list_output = inst.get_with_default(sname,"diff_list_output",
                                             "./rcd_difference_list.info")
    diff_matrix_output = inst.get_with_default(sname,"diff_matrix_output","")
    processes = inst.get_processes_limit(sname)

    coll = misc.load_collection_with_inst(inst,sname)
    coll.sort(key=lambda struct: struct.struct_id)

    result, diff_mat = calculate_collection_rcd_difference(coll, key=key,
                                                           ratio=ratio, 
                                                           select_pairs=pairs,
                                                           allow_enantiomer=not dis_en, 
                                                           processes=processes)
  
    f = open(diff_list_output,"a") 

    for k in result:
        f.write("%s %s %f \n" % (coll[k[0]].struct_id, coll[k[1]].struct_id, k[2]))
    f.close()

    if diff_matrix_output!="":
        f = open(diff_matrix_output,"a")
        for k in diff_mat:
            f.write(" ".join(map(str,k))+"\n")
        f.close()

    return result, diff_mat

def calculate_collection_rcd_difference(coll, key="RCD_vector",
                                        ratio=1, select_pairs=4, allow_enantiomer=True,
                                        processes=1):
    '''
    Calculates the rcd difference between each pairs of structures
    Outputs: a tuple of difference list and difference matrix
    '''
    arglist = [(i, j, coll[i].properties[key], coll[j].properties[key],
                ratio, select_pairs, allow_enantiomer)
                for i in range(len(coll)) for j in range(i,len(coll))]

    if processes > 1:
        p = multiprocessing.Pool(processes)
        result = p.map(_calculate_rcd_difference_multiprocess_wrapper,arglist)
    else:
        result = [_calculate_rcd_difference_multiprocess_wrapper(args) 
                                                    for args in arglist]

    diff_mat = [[0]*len(coll) for x in range(len(coll))]
    for k in result: 
        diff_mat[k[0]][k[1]] = k[2]
        diff_mat[k[1]][k[0]] = k[2]
    return result, diff_mat

   

def calculate_rcd_difference(v1, v2, ratio=1, select_pairs=4, 
                             allow_enantiomer=True):
    '''
    Given two rcd vectors, calculate their RCD_vector difference
    ratio: The ratio between the contributions to the final difference
           by relative orientation difference and relative coordinate difference
    allow_enantiomer: if True, then s2 will be mirror reflected 
                      by flipping the sign of its third relative coordinate.
                      This is necessary b/c the third reference axes is
                      generated with cross product, and thus will acquire a 
                      sign flip under mirror reflection  
    '''

    v1 = copy.deepcopy(v1)
    v2 = copy.deepcopy(v2)
    result = _calculate_rcd_vector_difference(v1, v2, ratio=ratio,
                                              select_pairs=select_pairs)
    if not allow_enantiomer:
        return result
    
    for x in v2:
        x[0][2] = - x[0][2]
    resulte = _calculate_rcd_vector_difference(v1, v2, ratio=ratio,
                                               select_pairs=select_pairs)

    return min(result, resulte)

def _calculate_rcd_difference_multiprocess_wrapper(arglist):
    return (arglist[0], arglist[1], calculate_rcd_difference(*arglist[2:]))

    

def _calculate_rcd_vector_difference(v1, v2, ratio=1, select_pairs=4):
    '''
    Given 2 RCD vectors, return their difference
    v1: RCD vector 1
    v2: RCD vector 2
    ratio: The difference will be the normalized Euclidean distance of relative coordinates 
           + Euclidean distance of relative orientation
    '''
    
    dist = [(x,y,_calculate_diff(v1[x],v2[y],ratio)) for y in range(len(v2))
                                                     for x in range(len(v1))]
    dist.sort(key=lambda x: x[2])

    result = 0
    s1 = [False] * len(v1)
    s2 = [False] * len(v2)
    dist_iter = iter(dist)
    for l in range(select_pairs):
        try:
            while True:
                k = dist_iter.next()
                if not s1[k[0]] and not s2[k[1]]:
                #Avoiding selecting same molecule from cell
                    result += k[2]
                    s1[k[0]] = True
                    s2[k[1]] = True
                   #print "Here is the selection", k
                    break
        except:
            raise RuntimeError("Not enough close picks to select the amount of pairs")

    return result


def _calculate_diff(v1, v2, ratio=1):
    relative_coordinate_diff = (np.linalg.norm(np.subtract(v1[0],v2[0]))**2
                                /np.linalg.norm(v1[0])/np.linalg.norm(v2[0]))
    relative_orientation_diff = np.linalg.norm(np.subtract(v1[1],v2[1]))**2/3
    #Divided by 6 to normalize (b/c the cosine value range from -1 to 1)
#    print "These are relative coordinates", v1[0], v2[0]
#    print "This is relative_coordinate_diff", relative_coordinate_diff
#    print "These are relative_orientations", v1[1], v2[1]
#    print "This is relative_orientation_diff", relative_orientation_diff
    return relative_coordinate_diff + ratio*relative_orientation_diff
    
