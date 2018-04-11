"""                                                                            
If any part of this module is used for a publication please cite:              
                                                                               
F. Curtis, X. Li, T. Rose, A. Vazquez-Mayagoitia, S. Bhattacharya,             
L. M. Ghiringhelli, and N. Marom "GAtor: A First-Principles Genetic            
Algorithm for Molecular Crystal Structure Prediction",                         
J. Chem. Theory Comput., DOI: 10.1021/acs.jctc.7b01152;                        
arXiv 1802.08602 (2018)                                                        
"""   
import os
import subprocess
from core import user_input
from core.file_handler import db_file, mkdir_p, cwd, structure_dir, read_data, \
    tmp_dir, write_data
from structures.structure import StoicDict
from structures.structure_collection import StructureCollection
from utilities import stoic_model
from structures import structure_collection
from utilities.stoic_model import determine_stoic
from external_libs.filelock import FileLock

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

ui = user_input.get_config()
INPUT_REF = -1

def main(argv):
    supercoll = {} 
    structure_collection.update_supercollection(supercoll)
    stoic = determine_stoic(os.path.join(cwd, argv[1]))
    if stoic == False: stoic = determine_stoic()
    if stoic == False: raise Exception
    structure_coll = StructureCollection(stoic, INPUT_REF)
    structure_coll.update_local()
    write_energy_hierarchy(structure_coll)
 
def get_energy_list(structure_coll):
    energy_list = []
    for index, structure in structure_coll:
        ID = structure.get_property('ID')
        replica = structure.get_property('replica')
        try: energy = '{:.3f}'.format(structure.get_property(ui.get_property_to_optimize()))
        except: energy = 0
        vol = '{:.1f}'.format(structure.get_unit_cell_volume())
        A, B, C = structure.get_lattice_magnitudes()
        a = '{:.2f}'.format(A)
        b = '{:.2f}'.format(B)
        c = '{:.2f}'.format(C)
        Alpha, Beta, Gamma = structure.get_lattice_angles()
        alpha = '{:.1f}'.format(Alpha)
        beta = '{:.1f}'.format(Beta)
        gamma = '{:.1f}'.format(Gamma)
        spacegroup = structure.get_property("space_group")
        mut = structure.get_property('mutation_type')
        parent0 = structure.get_property('parent_0')
        parent1 = structure.get_property('parent_1')
        cluster = structure.get_property('cluster_label')
        mem = structure.get_property('cluster_members')
        tot_clusters = structure.get_property('total_clusters')
        if parent1 is None:
            parent1 = ""
        if parent0 is None:
            parent0 = ""
        if mem is None:
            mem = ""
        if cluster is None:
            cluster = ""
        if tot_clusters is None:
            tot_clusters = ""
        if energy is not None:
            energy_list.append([ID, replica, index, energy, \
                                vol, a, b, c, alpha, beta, \
                                gamma, spacegroup, mut, \
                                parent0, parent1,\
                                cluster, mem, tot_clusters])
    return energy_list

def write_energy_hierarchy(structure_coll):
    to_write = ''
    ranked_energy_list = []
    energy_list = get_energy_list(structure_coll)
    energy_list.sort(key=lambda x: float(x[3]))
    for count in range(len(energy_list)):
        ranked_energy_list.append([count + 1] + energy_list[count])
    header = (['#Rank', 'Added', 'Replica', 'Index', 'Energy (eV)', 
                'Volume', 'A', 'B', 'C', 'Alpha', 'Beta', 'Gamma', 
                'SG', 'Mutation','ParentA', 'ParentB',
                'Cluster', 'Cluster_Members', 'Total_Clusters'])
    form = '{:<5} {:<5} {:<12} {:12} {:<12} '+\
           '{:<7} {:<6} {:<6} {:<6} {:<6} {:<6} {:<6} '+\
           '{:<3} {:<12} {:<11} {:<11} '+\
           '{:<7} {:<7} {:9}'
    to_write += form.format(*header) + '\n'
    for line in ranked_energy_list:
        to_write += form.format(*line) + '\n'
        stoic = structure_coll.get_stoic()
        input_ref = structure_coll.get_input_ref()
        filename = "energy_hierarchy_%s_%s.dat" % (stoic.get_string(),str(input_ref))
        f = open(os.path.join(tmp_dir,filename),"w")
        f.write(to_write)
        f.close()
    
