
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
from sklearn.decomposition import PCA
from sklearn import metrics
import multiprocessing

def main(structure_coll, replica):
    ''' 
    Inputs:  a StructureCollection and performs clustering
    Returns: a StructureCollection labeled with cluster indices
    '''
    start = time()
    
    # Get and check clustering algorithm and feature vector
    ui = user_input.get_config()
    feature_vector = ui.get("clustering", "feature_vector")
    cluster_type = ui.get("clustering","clustering_algorithm")
    if feature_vector is None:
        message = "Clustering algorithm called but feature vector is not"+\
                  "defined in ui.conf"
        raise ValueError(message)
    if cluster_type is None:
        message = "Clustering algorithm called but feature vector is not"+\
                  "defined in ui.conf"
        raise ValueError(message)

    # Create cluster instance 
    if cluster_type == "AffinityPropagation":
        if feature_vector == "RDF_vector":
            ClusterColl = AffinityPropagationClusteringRDF(structure_coll, replica)
        elif feature_vector == "PCA_RDF_vector":
            ClusterColl = AffinityPropagationPCAClusteringRDF(structure_coll, replica)
        elif feature_vector == "RCD_vector":
            ClusterColl = AffinityPropagationClusteringRCD(structure_coll, replica)
        elif feature_vector == "Lat_vol_vector":
            ClusterColl = AffinityPropagationClusteringLatVol(structure_coll, replica)
        else:
            message = "Affinity Propagation not implemented for %s" % (feature_vector)
            raise RuntimeError(message)
    elif cluster_type == "Kmeans":
        if feature_vector == "RDF_vector":
            ClusterColl = KMeansClusteringRDF(structure_coll, replica)
        elif feature_vector == "Lat_vol_vector":
            ClusterColl = KMeansClusteringLatVol(structure_coll, replica)
        else:
            message = "K-means clustering not implemented for %s" % (feature_vector)
            raise RuntimeError(message)
    else:
        message = "Clustering type %s not implemented. Check spelling" % (cluster_type)
        raise RuntimeError(message)

    # Return Clustered Collection
    clustered_coll = ClusterColl.return_clusters()

    #Outputs
    end = time()
    message = "-- Time for clustering: %s seconds" % (end-start)
    output.local_message(message, replica)
    
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
                AFV = AssignFeatureVector(struct)
                struct = AFV.compute_feature_vector()
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
            struct.set_property("total_clusters", len(info))
            i += 1
            for j in info:
                if label == j[0]:
                    struct.set_property('cluster_members', j[1])
        # Output info
        self.output("-- Clustering Feature vector %s" % (self.feature_type))
        self.output("-- Number of clusters %s" % (len(info)))
        self.output("-- Distribution of clusters %s" % (sorted_info))

class KMeansClusteringLatVol():
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
        feature_list = np.array(self.return_descriptor_list())
        KM = KMeans(n_clusters=self.num_clusters, init='k-means++')
        clustered_data = KM.fit_predict(feature_list)
        clustered_coll = self.cluster_coll(clustered_data)
        return self.struct_coll

    def output(self, message):
        output.local_message(message, self.replica)

    def return_descriptor_list(self):
        lat_vols = []
        for index, struct in self.struct_coll:
            lat_vol = struct.get_property(self.feature_type)
            if lat_vol is None:
                AFV = AssignFeatureVector(struct)
                struct = AFV.compute_feature_vector()
                lat_vol = struct.get_property(self.feature_type)
            lat_vols.append(lat_vol)
        return lat_vols

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
            struct.set_property("total_clusters", len(info))
            i += 1
            for j in info:
                if label == j[0]:
                    struct.set_property('cluster_members', j[1])
        # Output info
        self.output("-- Clustering Feature vector %s" % (self.feature_type))
        self.output("-- Number of clusters %s" % (len(info)))
        self.output("-- Distribution of clusters %s" % (sorted_info))

class AffinityPropagationClusteringRDF():
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

        af = AffinityPropagation(convergence_iter=75).fit(feature_list)
        cluster_centers_indices = af.cluster_centers_indices_
        clustered_data = af.labels_
          
        n_clusters_ = len(cluster_centers_indices)
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
                AFV = AssignFeatureVector(struct)
                struct = AFV.compute_feature_vector()
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
            struct.set_property("total_clusters", len(info))
            i += 1
            for j in info:
                if label == j[0]:
                    struct.set_property('cluster_members', j[1])

        # Output info
        self.output("-- Clustering Feature vector %s" % (self.feature_type))
        self.output("-- Number of clusters %s" % (len(info)))
        self.output("-- Distribution of clusters %s" % (sorted_info))

class AffinityPropagationPCAClusteringRDF():
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
        self.pca_components = self.ui.get_eval("clustering", "pca_components")

    def return_clusters(self):
        feature_list = np.array(self.return_RDF_list())

        pca = PCA(n_components=self.pca_components)
        pca_reduced_features = pca.fit_transform(feature_list)
        af = AffinityPropagation(convergence_iter=75).fit(pca_reduced_features)
        cluster_centers_indices = af.cluster_centers_indices_
        clustered_data = af.labels_

        n_clusters_ = len(cluster_centers_indices)
        clustered_coll = self.cluster_coll(clustered_data)
        return self.struct_coll

    def output(self, message):
        output.local_message(message, self.replica)

    def return_RDF_list(self):
        RDFs = []
        for index, struct in self.struct_coll:
            RDF = struct.get_property("RDF_vector")
            if RDF is not None:
                RDFs.append(RDF)
            elif RDF is None:
                AFV = AssignFeatureVector(struct)
                struct = AFV.compute_feature_vector()
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
            struct.set_property("total_clusters", len(info))
            i += 1
            for j in info:
                if label == j[0]:
                    struct.set_property('cluster_members', j[1])
        
        # Output info
        self.output("-- Clustering Feature vector %s" % (self.feature_type))
        self.output("-- Number of clusters %s" % (len(info)))
        self.output("-- Distribution of clusters %s" % (sorted_info))


class AffinityPropagationClusteringLatVol():
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
        feature_list = np.array(self.return_descriptor_list())
        af = AffinityPropagation(convergence_iter=75).fit(feature_list)
        cluster_centers_indices = af.cluster_centers_indices_
        clustered_data = af.labels_

        n_clusters_ = len(cluster_centers_indices)
        self.output('Estimated number of clusters: %d' % n_clusters_)
        self.output("Silhouette Coefficient: %0.3f"
            % metrics.silhouette_score(feature_list, clustered_data, metric='sqeuclidean'))

        clustered_coll = self.cluster_coll(clustered_data)
        return self.struct_coll

    def output(self, message):
        output.local_message(message, self.replica)

    def return_descriptor_list(self):
        lat_vols = []
        for index, struct in self.struct_coll:
            lat_vol = struct.get_property(self.feature_type)
            if lat_vol is None:
                AFV = AssignFeatureVector(struct)
                struct = AFV.compute_feature_vector()
                lat_vol = struct.get_property(self.feature_type)
            lat_vols.append(lat_vol)
        return lat_vols

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
            struct.set_property("total_clusters", len(info))
            i += 1
            for j in info:
                if label == j[0]:
                    struct.set_property('cluster_members', j[1])

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
        self.replica = replica
        self.struct_coll = structure_coll
        self.feature_type = self.ui.get("clustering", "feature_vector")
        self.affinity_type = "exponential"
        self.damping = 0.5
        self.convergence_iter = 75
        self.max_iter = 200
        self.preference = None
        self.stored_property_key = "cluster_label"
        self.dist_from_center = "AP_cluster_distance"
        self.cluster_mem = "cluster_members"
        self.tot_clusters = "total_clusters"

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
            struct.set_property(self.tot_clusters, len(info))
            for j in info:
                if label == j[0]:
                    struct.set_property(self.cluster_mem, j[1])

        # Outputs
        self.output("-- Clustering Feature vector %s" % (self.feature_type))
        self.output("-- Number of clusters %s" % (len(info)))
        self.output("-- Distribution of clusters %s" % (sorted_info))
        return self.struct_coll

    def return_RCD_collection(self):
        start = time()
        for index, struct in self.struct_coll:
            RCD_vector = struct.get_property(self.feature_type)
            if RCD_vector is None:
                AFV = AssignFeatureVector(struct)
                struct = AFV.compute_feature_vector()
                struct.set_property(self.feature_type, RCD_vector)
        end = time()
        message = "-- Time for clustering: %s seconds" % (end-start)
        self.output(message)
        return 
       
    def output(self, message):
        output.local_message(message, self.replica)


class AssignFeatureVector():
    '''
    Computes the feature vector for a given structure
    Returns Structure with feature assigned as a property
    '''
    def __init__(self, struct):
        self.ui = user_input.get_config()
        self.feature_type = self.ui.get("clustering", "feature_vector")
        self.num_mols = self.ui.get_eval('run_settings', 'num_molecules') 
        self.struct = struct
        self.num_atoms_per_mol = int(len(self.struct.geometry)/self.num_mols)

    def compute_feature_vector(self):
        if self.feature_type == "RDF_vector":
            struct = structure_handling.compute_RDF_vector(self.struct)
            return struct
        elif self.feature_type == "PCA_RDF_vector":
            struct = structure_handling.compute_RDF_vector(self.struct)
            return struct
        elif self.feature_type == "Lat_vol_vector":
            a = np.linalg.norm(self.struct.get_property('lattice_vector_a'))
            b = np.linalg.norm(self.struct.get_property('lattice_vector_b'))
            c = np.linalg.norm(self.struct.get_property('lattice_vector_c'))
            vol = self.struct.get_unit_cell_volume()
            #vol = np.cbrt(vol)
            vol = vol**(1./3)
            lat_vol = [a/vol, b/vol, c/vol]
            self.struct.set_property(self.feature_type, lat_vol)
            return self.struct
        elif self.feature_type == "RCD_vector": 
            axes_indices = list(self.ui.get_list('clustering','rcd_axes_indices'))
            axes_indices = [int(i) for i in axes_indices]
            axes_indices = ([(axes_indices[0], axes_indices[1]), 
                             (axes_indices[2], axes_indices[3])])   
            close_pics = int(self.ui.get_eval('clustering','rcd_close_picks'))
            RCD_struct = generate_relative_coordinate_descriptor(self.struct, 
                                              self.num_mols, 
                                              self.num_atoms_per_mol, 
                                              axes_indices, 
                                              close_picks=8)
            RCD_vector = RCD_struct.get_property(self.feature_type)
            self.struct.set_property(self.feature_type, RCD_vector)
            return self.struct
        else:
            message = "Clustering for %s is not availble" % (self.feature_type)
            raise RuntimeError(message)
