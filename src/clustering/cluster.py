
"""
Created on Tues Jan 10 12:59:12 2017

@author: Farren Curtis
"""
from __future__ import division


import numpy as np
from core import user_input, output
from time import time
from structures import structure_collection, structure_handling
from structures.structure_collection import StructureCollection
from utilities.relative_coordinate_descriptor import *
from utilities.affinity import *
from utilities.affinity import _convert_distance_to_affinity_mat
from sklearn.cluster import KMeans, AffinityPropagation
from sklearn import metrics
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
    elif cluster_type == "AffinityPropagationRDF":
        ClusterColl = AffinityPropagationClusteringRDF(structure_coll, replica)
        clustered_coll = ClusterColl.return_clusters()
    elif cluster_type == "AffinityPropagationRCD":
        ClusterColl = AffinityPropagationClusteringRCD(structure_coll, replica)
        clustered_coll = ClusterColl.return_clusters()
    elif cluster_type == "AffinityPropagationBendAngle":
        ClusterColl = AffinityPropagationClusteringBendAngle(structure_coll, replica)
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
        self.feature_type = self.ui.get("clustering", "feature_vector")


    def return_clusters(self):
        feature_list = np.array(self.return_RDF_list())
        KM = KMeans(n_clusters=self.num_clusters, init='k-means++')
        clustered_data = KM.fit_predict(feature_list)
        clustered_coll = self.cluster_coll(clustered_data)

        print("Silhouette Coefficient: %0.3f"
            % metrics.silhouette_score(feature_list, clustered_data, metric='sqeuclidean'))
        return self.struct_coll

    def output(self, message):
        output.local_message(message, self.replica)

    def return_RDF_list(self):
        RDFs = []
        for index, struct in self.struct_coll:
            RDF = struct.get_property(self.feature_type)
            if RDF is not None:
                RDFs.append(RDF)
            elif RDF is None:
                struct = structure_handling.compute_RDF_vector(struct)
                RDF = struct.get_property(self.feature_type)
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
        print info
        # Output info
        self.output("-- Clustering Feature vector %s" % (self.feature_type))
        self.output("-- Number of clusters %s" % (len(info)))
        self.output("-- Distribution of clusters %s" % (sorted_info))
        print "-- Number of clusters %s" % (len(info))

class AffinityPropagationClusteringRDF():
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
        self.feature_type = self.ui.get("clustering", "feature_vector")


    def return_clusters(self):
        feature_list = np.array(self.return_RDF_list())

        af = AffinityPropagation().fit(feature_list)
        cluster_centers_indices = af.cluster_centers_indices_
        clustered_data = af.labels_
          
        n_clusters_ = len(cluster_centers_indices)
        print "here"
        print('Estimated number of clusters: %d' % n_clusters_)
        print("Silhouette Coefficient: %0.3f"
            % metrics.silhouette_score(feature_list, clustered_data, metric='sqeuclidean'))



        clustered_coll = self.cluster_coll(clustered_data)
        return self.struct_coll

    def output(self, message):
        output.local_message(message, self.replica)

    def return_RDF_list(self):
        RDFs = []
        for index, struct in self.struct_coll:
            RDF = struct.get_property(self.feature_type)
            if RDF is not None:
                RDFs.append(RDF)
            elif RDF is None:
                struct = structure_handling.compute_RDF_vector(struct)
                RDF = struct.get_property(self.feature_type)
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
        print info

        # Output info
        self.output("-- Clustering Feature vector %s" % (self.feature_type))
        self.output("-- Number of clusters %s" % (len(info)))
        self.output("-- Distribution of clusters %s" % (sorted_info))

class AffinityPropagationClusteringBendAngle():
    '''
    This class performs K means clustering for a fixed number of 
    user-defined clusters on previously defined RDF vectors in each
    structure in a structure collection
    '''
    def __init__(self, structure_coll, replica):
        self.ui = user_input.get_config()
        self.replica = replica
        self.struct_coll = structure_coll
        self.feature_type = self.ui.get("clustering", "feature_vector")


    def return_clusters(self):
        feature_list = np.array(self.return_RDF_list())

        af = AffinityPropagation().fit(feature_list)
        cluster_centers_indices = af.cluster_centers_indices_
        clustered_data = af.labels_

        n_clusters_ = len(cluster_centers_indices)
        print "here"
        print('Estimated number of clusters: %d' % n_clusters_)
        print("Silhouette Coefficient: %0.3f"
            % metrics.silhouette_score(feature_list, clustered_data, metric='sqeuclidean'))



        clustered_coll = self.cluster_coll(clustered_data)
        return self.struct_coll

    def output(self, message):
        output.local_message(message, self.replica)

    def return_descriptor_list(self):
        RDFs = []
        for index, struct in self.struct_coll:
            RDF = struct.get_property(self.feature_type)
            if RDF is not None:
                RDFs.append(RDF)
            elif RDF is None:
                struct = structure_handling.compute_RDF_vector(struct)
                RDF = struct.get_property(self.feature_type)
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
        print info

        # Output info
        self.output("-- Clustering Feature vector %s" % (self.feature_type))
        self.output("-- Number of clusters %s" % (len(info)))
        self.output("-- Distribution of clusters %s" % (sorted_info))


class AffinityPropagationClusteringRCD():
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
        self.feature_type = self.ui.get("clustering", "feature_vector")
        self.affinity_type = "exponential"
        self.damping = 0.5
        self.convergence_iter = 15
        self.max_iter = 200
        self.preference = None
        self.stored_property_key = "cluster_label"
        self.dist_from_center = "AP_cluster_distance"
        self.cluster_mem = "cluster_members"

    def return_clusters(self):
        # Find RCD clusters
        self.return_RCD_collection()
        coll = []
        for index, struct in self.struct_coll:
            coll.append(struct)
        coll.sort(key=lambda struct: struct.struct_id)
        result, diff_mat = calculate_collection_rcd_difference(coll)

        coll, centers = AP_distance_matrix(coll, diff_mat, 
                                           affinity_type=self.affinity_type,
                                           damping=self.damping, 
                                           convergence_iter=self.convergence_iter,
                                           max_iter=self.max_iter, 
                                           preference=self.preference,
                                           stored_property_key=\
                                           self.stored_property_key)

        # Assign Cluster Labels and Distances to Centers to each Structure
        clusters = []
        for index, struct in self.struct_coll:
            coll_ind = coll.index(struct)
            coll_struct = coll[coll_ind]
            AP_cluster_dist = (diff_mat[coll.index(struct)]
                              [coll.index(centers[struct.properties\
                              [self.stored_property_key]])])
            cluster_label = coll_struct.properties[self.stored_property_key]
            clusters.append(cluster_label)
            struct.set_property(self.stored_property_key, cluster_label)  
            struct.set_property(self.dist_from_center, AP_cluster_dist)

        # Assign Number of Cluster members to each Structure
        info = [[x, clusters.count(x)] for x in set(clusters)]
        infox = [x[0] for x in info]
        sorted_info = sorted(info, key=lambda x: x[1])

        for index, struct in self.struct_coll:
            label = struct.get_property(self.stored_property_key)
            for j in info:
                if label == j[0]:
                    struct.set_property(self.cluster_mem, j[1])

        # Outputs
        print ("Number of clusters generated: %i\n" % len(centers))
        print ("List of cluster centers:\n" +
            "\n".join(map(str,[struct.struct_id for struct in centers])) + "\n")
        print("Assigned cluster labels:\n")
        print("\n".join(map(str,[struct.struct_id + " " +
                            str(struct.properties[self.stored_property_key])+" "+
                            str(diff_mat[coll.index(struct)]
                            [coll.index(centers[struct.properties\
                            [self.stored_property_key]])])
                            for struct in coll])))

        self.output("-- Clustering Feature vector %s" % (self.feature_type))
        self.output("-- Number of clusters %s" % (len(info)))
        self.output("-- Distribution of clusters %s" % (sorted_info))
        return self.struct_coll

    def return_RCD_collection(self):
        start = time()
        for index, struct in self.struct_coll:
            RCD_vector = struct.get_property(self.feature_type)
            if RCD_vector is None:
                RCD_struct = rcd_vector_calculation(struct)
                RCD_vector = RCD_struct.get_property(self.feature_type)
                struct.set_property(self.feature_type, RCD_vector)
        end = time()
        message = "-- Time for clustering: %s seconds" % (end-start)
        print message
        return 
       
    def output(self, message):
        output.local_message(message, self.replica)


        
