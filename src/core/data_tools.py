'''
Created on Oct 30, 2013

@author: newhouse

This module's purpose is to provide tools to the user wishing to access the data stored in the database

Implemented:

To Implement:
Print data to filesystem
Get lowest energy structure

'''
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


INPUT_REF = -1

def main(argv):
    '''
    enter here the functions you wish to execute when this module is run
    '''
    ui = user_input.get_config()
    supercoll = {} 
    structure_collection.update_supercollection(supercoll)
#     write_to_filesystem()
    stoic = determine_stoic(os.path.join(cwd, argv[1]))
    if stoic == False: stoic = determine_stoic()
    if stoic == False: raise Exception
    structure_coll = StructureCollection(stoic, INPUT_REF)
    structure_coll.update_local()
    energy_tuples = get_energy_tuples(structure_coll)
    energy_tuples.sort(key=lambda x: x[1])
    loe = [item[1] for item in energy_tuples]
    index = energy_tuples[0][0]
#     index, ignore, energy = structure_index(STOIC, 633)
    best_structure = structure_coll.get_struct(index)
    try:
        jmol(best_structure, str(index))
    except: pass
    try:
        plot_energies(loe)  # ,len(loe) - 50)
    except: pass
    print structure_coll.get_struct(energy_tuples[0][0]).get_geometry_atom_format()
    write_energy_hierarchy(structure_coll)


#modified from original to include max distances for checking    
def get_energy_tuples(structure_coll):
    energy_tuples = []
    for index, structure in structure_coll:
        energy = structure.get_property('energy') 
	vol = structure.get_property('cell_vol')
	a = structure.get_property('a')
	b = structure.get_property('b')
	c = structure.get_property('c')
	alpha = structure.get_property('alpha')
	beta = structure.get_property('beta')
	gamma = structure.get_property('gamma')
	parent0 = structure.get_property('parent_0')
	parent1 = structure.get_property('parent_1')
	crosstype = structure.get_property('crossover_type')
	localmindiff = structure.get_property('new_local_minima_diff')
	childnum = structure.get_property('child_count')
	ID = structure.get_property('ID')
	if energy is not None: energy_tuples.append((index, energy, ID, childnum, localmindiff, vol, a, b, c, alpha, beta, gamma, crosstype, str(parent0)[16:], str(parent1)[16:]))
#	print energy_tuples
    return energy_tuples
 

def make_user_structure_directory(structure_coll, energy_tuples, n_structures=10):
    path = os.path.join(tmp_dir, 'user_structures')
    mkdir_p(path)
    for i in n_structures:
        try:
            struct = structure_coll.get_struct(energy_tuples[i][0])
            write_data(path, struct.get_stoic_str() + '_' + energy_tuples[i][0], struct.get_geo_atom_format())
        except: raise Exception
     

def write_energy_hierarchy(structure_coll):
    energy_tuples = get_energy_tuples(structure_coll)
 #   biggest_dists = get_biggest_dist(structure_coll)
    to_write = ''
    energy_tuples.sort(key=lambda x: x[1])
    for index, Id, energy, vol, a, b, c, al, be, ga, p0, p1, ct, md, cn in energy_tuples: 
#	to_write += structure_coll.get_stoic().get_string() + '/'
        to_write += str(structure_coll.get_input_ref()) + '/'
        to_write += str(index) + '/'
	to_write +='	' + str(Id)
        to_write +='    ' + str(energy)
	to_write +='    ' + str(vol)
	to_write +='    ' + str(a)
	to_write +='    ' + str(b)
	to_write +='    ' + str(c)
	to_write +='    ' + str(al)
	to_write +='    ' + str(be)
	to_write +='    ' + str(ga)
	to_write +='    ' + str(p0)
	to_write +='    ' + str(p1)	
	to_write +='    ' + str(ct)
	to_write +='    ' + str(md)
	to_write +='    ' + str(cn)
        to_write += '\n'
    with open(os.path.join(tmp_dir, 'energy_hierarchy.' + str(structure_coll.get_input_ref()) + '.dat'), 'w') as f: f.write(to_write)
        
def write_geometries_to_file(structure_coll):
    pass
    
def plot_energies(energy_list, skip_first_n=0):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    energy_list.sort(reverse=True)
    x = []
    y = []
    if skip_first_n < len(energy_list):
        for i in range(skip_first_n):
            energy_list.pop(0)        
    for i in range(len(energy_list)):
        x.append(i)
        y.append(energy_list[i])
    ax1 = plt.gcf().add_subplot(111).scatter(x, y, c=y, s=50)
    plt.title("Best Energy: " + str(energy_list[-1]))
    plt.savefig(os.path.join(cwd,'energy_plot.eps'), format='eps', dpi=1000)
    subprocess.Popen(['xdg-open', os.path.join(cwd,'energy_plot.eps')])

def jmol(struct, name):
    path = os.path.join(cwd, 'tmp', 'jmol')
    mkdir_p(os.path.join(cwd, 'tmp', 'jmol'))
    atom_file = open(os.path.join(path , name + '.dat'), 'w')
    atom_file.write(struct.get_geometry_atom_format())
    atom_file.close()
    input_file = open(os.path.join(path, 'exe.sh'), 'w')
    input_file.write('jmol ' + os.path.join(path, name + '.dat') + '  >/dev/null' + '&')  # supress out
    input_file.close()
    # creates shell script with proper syntax and runs it as a subprocess
    p = subprocess.Popen(['sh', 'exe.sh'], cwd=os.path.join(cwd, 'tmp', 'jmol'))
    p.wait()
#     if os.path.exists(path): rmtree(path)
    
if __name__ == '__main__':
    main('stoic.conf')
