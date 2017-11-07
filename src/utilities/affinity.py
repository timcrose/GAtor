'''
Created by Patrick Kilecdi on 12/31/2016

Conducts Affinity Propagation for a given collection
'''

from sklearn.cluster import AffinityPropagation
import numpy as np
from utilities import misc

    

def AP_distance_matrix(coll, dist_mat, affinity_type=["exponential",1],
                       damping=0.5, convergence_iter=75, max_iter=200,
                       preference=None, stored_property_key="AP_cluster"):
    '''
    Given a collection of structures and their distance matrix
    Conduct the Affinity Propagation on the collection
    '''
    affinity_mat = _convert_distance_to_affinity_mat(dist_mat, affinity_type)
    return AP_affinity_matrix(coll, affinity_mat, damping=damping,
                              convergence_iter=convergence_iter, 
                              max_iter=max_iter, preference=preference,
                              stored_property_key=stored_property_key)


def AP_affinity_matrix(coll, affinity_mat, damping=0.5, convergence_iter=75,
                       max_iter=200, preference=None, stored_property_key="AP_cluster"):
    '''
    Given a collection of structure and their affinity matrix
    Conduct the Affinity Propagation on the collection
    Returns a full list of collection with cluster assigned,
    as well as a list of examplars
    '''
    ap = AffinityPropagation(damping=damping, max_iter=max_iter, 
                             convergence_iter=convergence_iter,
                             copy=True, preference=preference,
                             affinity="precomputed",verbose=False)

    
    result = ap.fit(affinity_mat)
    for x in range(len(coll)):
        coll[x].properties[stored_property_key] = result.labels_[x]

    centers = [coll[x] for x in result.cluster_centers_indices_]
    return coll, centers

def _convert_distance_to_affinity_mat(dist_mat, affinity_type=["exponential",1]):
    '''
    Converts a given distance matrix to affinity matrix
    '''
    m = len(dist_mat); n = len(dist_mat[0])
    affinity_mat = [[0]*n for x in range(m)]

    print(affinity_type)
    #print m, n
    for i in range(m):
        for j in range(n):
            #print i, j
            if affinity_type[0]=="exponential":
                affinity_mat[i][j] = - np.exp(dist_mat[i][j]*affinity_type[1])
            elif affinity_type[0]=="power":
                affinity_mat[i][j] = - dist_mat[i][j]**affinity_type[1]

    return affinity_mat

      
