'''
@authors: newhouse, farren

This module's purpose is to provide tools to the user wishing to access the data stored in the database

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
from external_libs.filelock import FileLock


ui = user_input.get_config()
INPUT_REF = -1

def main(argv):
    '''
    enter here the functions you wish to execute when this module is run
    '''
    supercoll = {} 
    structure_collection.update_supercollection(supercoll)
    stoic = determine_stoic(os.path.join(cwd, argv[1]))
    if stoic == False: stoic = determine_stoic()
    if stoic == False: raise Exception
    structure_coll = StructureCollection(stoic, INPUT_REF)
    structure_coll.update_local()
    energy_list = get_energy_list(structure_coll)
    energy_list.sort(key=lambda x: x[1])
    loe = [item[1] for item in energy_list]
    index = energy_list[0][0]
    best_structure = structure_coll.get_struct(index)
    try:
        jmol(best_structure, str(index))
    except: pass
    try:
        plot_energies(loe)  # ,len(loe) - 50)
    except: pass
    print structure_coll.get_struct(energy_list[0][0]).get_geometry_atom_format()
    write_energy_hierarchy(structure_coll)

def write_avg_fitness(Ig, fit_avg, structure_coll):
    to_write = ''
    to_write += str(ID)+ ' '+ str(fit_avg)
    to_write += '\n'	
    with open(os.path.join(tmp_dir, 'fitness.' + str(structure_coll.get_input_ref()) + '.dat'), 'a+') as f: f.write(to_write)
    os.system("chmod g=u "+os.path.join(tmp_dir,'fitness.'+str(structure_coll.get_input_ref())+'.dat'))
 
def get_energy_list(structure_coll):
    energy_list = []
    for index, structure in structure_coll:
        ID = structure.get_property('ID')
        replica = structure.get_property('replica')
        energy = '{:.3f}'.format(structure.get_property(ui.get_property_to_optimize()))
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
        crosstype = structure.get_property('crossover_type')
        parent0 = structure.get_property('parent_0')
        parent1 = structure.get_property('parent_1')
        cluster = structure.get_property('cluster_label')
        #apcluster = structure.get_property('af_cluster_label')
        #apmem = structure.get_property('af_cluster_members')
        mem = structure.get_property('cluster_members')
        rev = structure.get_property('revisits')
        if energy is not None:
            energy_list.append([ID, replica, index, energy, \
                                vol, a, b, c, alpha, beta, \
                                gamma, spacegroup, mut, \
                                str(parent0)[15:], str(parent1)[15:],\
                                crosstype, cluster, mem, rev])
    return energy_list

def write_energy_vs_addition(structure_coll):
    energy_list = get_energy_list(structure_coll)
    to_write = ''
    energy_list.sort(key=lambda x: x[0])
    for  Id, rep, index, energy, vol, a, b, c, al, be, ga, sg, mut, par0, par1, crosst, clus, mem, st in energy_list:
        if rep == 'init_pool':
            continue
        to_write += str(index)+"    "+str(Id)+"    "+str(energy)+'\n' 
        file_path = os.path.join(tmp_dir, 'energy_vs_addition.' + str(structure_coll.get_input_ref()) + '.dat')
        with open(file_path, 'w') as f:
            f.write(to_write)
            ui.grant_permission(file_path)

def write_spe_vs_addition(structure_coll):
    energy_list = get_energy_list(structure_coll)
    to_write = ''
    energy_list.sort(key=lambda x: x[0])
    for  Id, rep, index, energy, spe, vol, a, b, c, al, be, ga, sg, mut, par0, par1, crosst, ap in energy_list:
        if rep == 'init_pool':
            continue
        to_write += str(Id)+'    '+str(spe)+'\n'
        with open(os.path.join(tmp_dir, 'sp_energy_vs_addition.' + str(structure_coll.get_input_ref()) + '.dat'), 'w') as f:
            f.write(to_write)

def write_energy_hierarchy(structure_coll):
	to_write = ''
	ranked_energy_list = []
	energy_list = get_energy_list(structure_coll)
	energy_list.sort(key=lambda x: float(x[3]))
	for count in range(len(energy_list)):
		ranked_energy_list.append([count + 1] + energy_list[count])
	header = (['#Rank', 'Added', 'Replica', 'Index', 'Energy (eV)', 
                'Volume', 'A', 'B', 'C', 'Alpha', 'Beta', 'Gamma', 
                'Spacegroup', 'Mutation','ParentA', 'ParentB',
                'Crossover', "Cluster", "Cluster_Members", "Revisits"])
	form = '{:<5} {:<5} {:<12} {:20} {:<12} {:<7} {:<6} {:<6} {:<6} {:<6} {:<6} {:<6} {:<10} {:<16} {:<8} {:<8} {:<30} {:<7} {:<7} {:<7}'
	to_write += form.format(*header) + '\n'
	for line in ranked_energy_list:
		to_write += form.format(*line) + '\n'
    	stoic = structure_coll.get_stoic()
    	input_ref = structure_coll.get_input_ref()
    	filename = "energy_hierarchy_%s_%s.dat" % (stoic.get_string(),str(input_ref))
        f = open(os.path.join(tmp_dir,filename),"w")
        f.write(to_write)
        f.close()
        ui.grant_permission(os.path.join(tmp_dir,filename))

def make_user_structure_directory(structure_coll, energy_list, n_structures=10):
    path = os.path.join(tmp_dir, 'user_structures')
    mkdir_p(path)
    for i in n_structures:
        try:
            struct = structure_coll.get_struct(energy_list[i][0])
            write_data(path, struct.get_stoic_str() + '_' + energy_list[i][0], struct.get_geo_atom_format())
        except: raise Exception
             
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
