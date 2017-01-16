# take in a structure collection and perform clustering on it,

# assign cluster label, total members in cluster  to each structure as property

# then in selection normalize fitness by dividing by this number of members
# could also look into only mating structures from same clusters (in selection)





"""
Created on Tues Jan 10 12:59:12 2017

@author: Farren Curtis
"""
from __future__ import division


import numpy as np
from core import user_input, output
from time import time
from structures import structure_collection
from structures.structure_collection import StructureCollection
from sklearn.cluster import KMeans

import multiprocessing

def main(structure_coll, replica):
    ''' 
    Inputs:  a StructureCollection and performs clustering
    Returns: a StructureCollection labeled with cluster indices
    '''
    start = time()

    ui = user_input.get_config()
    cluster_type = ui.get("clustering","clustering_algorithm")
    if cluster_type == "KMeansRDF":
        ClusterColl = KMeansClusteringRDF(structure_coll, replica)
        clustered_coll = ClusterColl.return_clusters()

    end = time()
    message = "-- Time for clustering: %s seconds" % (end-start)
    output.local_message(message, replica)
    print "Time %s" %(end-start)
    return clustered_coll


class KMeansClusteringRDF():
    '''
    This class performs K means clustering for a fixed number of 
    user-defined clusters on previously defined RDF vectors in each
    structure in a structure collection
    '''
    def __init__(self, structure_coll, replica):
        self.ui = user_input.get_config()
        self.num_clusters = self.ui.get_eval("clustering", "num_clusters")
        self.replica = replica
        self.struct_coll = structure_coll

    def output(self, message):
        output.local_message(message, self.replica)

    def return_clusters(self):
        feature_list = self.return_RDF_list()
        KM = KMeans(n_clusters=self.num_clusters, init='k-means++')
        clustered_data = KM.fit_predict(feature_list)
        clustered_coll = self.cluster_coll(clustered_data)
        return self.struct_coll

    def return_RDF_list(self):
        RDFs = []
        for index, struct in self.struct_coll:
            RDF = struct.get_property("RDF_smooth")
            RDFs.append(RDF)
        return RDFs

    def cluster_coll(self, clustered_data):
        ''' 
        Takes cluster labels and assigns as properties 
        of each structure in collection
        '''    
        # Find clusters and number of members in each
        clusters = clustered_data.tolist()
        info = [[x, clusters.count(x)] for x in set(clusters)]
        infox = [x[0] for x in info]
        sorted_info = sorted(info, key=lambda x: x[1])

        # Assign cluster label as property of each structure
        i = 0
        for index, struct in self.struct_coll:
            label = clusters[i]
            struct.set_property("cluster_label", label)
            i += 1
            for j in info:
                if label == j[0]:
                    struct.set_property('cluster_members', j[1])

        # Output info
        self.output("Number of clusters %s" % (len(info)))
        self.output("Distribution of clusters %s" % (sorted_info))






        
