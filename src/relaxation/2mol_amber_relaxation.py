'''
Created on Jul 29, 2013

@author: newhouse
'''
import os
import re
import numpy as np
from shutil import copyfile
import shutil
import subprocess
import time
from copy import deepcopy

from core import user_input, output
from core.file_handler import mkdir_p_clean, res_dir
from structures.structure import Structure
from utilities.element_masses import masses


def main(input_structure, working_dir, control_in_string, replica):
    '''
    For a relaxation module, the method main() must be implemented.
    
    It must take as arguments a 'Structure' and an 'int' which defines 
    in which working directory the relaxation will take place.
    
    It must return a Structure object. (Please see structure documentation
    for information on how to easily setup a Structure object) 
    '''
    try: 
        r = AmberRelaxation(input_structure, working_dir, control_in_string, replica)
        r.setup()
        r.execute()
        return r.extract()
    except Exception, e: print str(e); return False # if there is an error in lammps this will prevent crashing
    pass

class AmberRelaxation():
    '''
    Creates an relaxation object using lammps to execute a geometry.
    
    This is an example of using a class to define a state, in this case the
    state of the relaxation process throughout setup, execution, and extraction
    '''
    def __init__(self, input_structure, working_dir, control_in_string, replica):
        '''
        __init__ will always run when a class is initialized. 
        This is also how the class can take arguments.
        '''
        # defining 'self' fields gives access to these values throughout the class 
        self.replica = replica
        self.working_dir = working_dir
        self.input_structure = input_structure
        self.control_in_string = control_in_string
        # this gives easy access to retrieve preferences from the user_input file
        self.ui = user_input.get_config()



    def output(self, message): output.local_message(message, self.replica)

    def setup(self):
        '''
        Creates the environment for the optimizer and constructs the input files.
        In this case, TINKER will take a geometry.dat file and an geometry.key
        '''
        mkdir_p_clean(self.working_dir)
        self.set_element_list()
        self.create_geometry_in()
        self.create_input()
#        shutil.copy(os.path.join(res_dir, 'ffield', 'reax.ffield'), os.path.join(self.working_dir, 'ffield'))
    
    def execute(self):
        '''
        Creates an executable shell script with proper LAMMPS syntax 
        then executes in a subprocess
        
        '''
        # creates shell script with proper syntax and runs it as a subprocess
        exe_file = open(os.path.join(self.working_dir, 'exe.sh'), 'w')
#	exe_file.write('module purge\n')
#        exe_file.write('module load intel/13.1\n')
#        exe_file.write('module load impi\n')
#        exe_file.write('module load lammps\n')
        exe_file.write(self.ui.get('TINKER', 'path_to_tinker_executable') + ' geometry_amber.xyz -k geometry_amber.key .01' + ' >' + os.path.join(self.working_dir, 'log.tinker'))  # write log        
        exe_file.close()

        p = subprocess.Popen(['sh', 'exe.sh'], cwd=self.working_dir)
        p.wait()    
        
    def extract(self):
        '''
        Reads the output files from relaxation and specifies properties of the structure (e.g. energy)
        
        Returns: Structure or False
        '''
        # creates empty structure
        self.result_struct = Structure()
        
        # reads output file to extract energy
        energy = self.extract_energy()
        if energy == False: return False  # did not converge
        self.result_struct.set_property('energy', energy)
        
        # sets the result structure's geometry based on output file
        extract_success = self.extract_geometry()
        if not extract_success: return False
        if not self.input_structure.get_stoic() == self.result_struct.get_stoic(): return False
 # 	self.result_struct.set_property('lattice_vector_a', alat)
#	self.result_struct.set_property('lattice_vector_b', blat)
#	self.result_struct.set_property('lattice_vector_c', clat)      
        return self.result_struct
        
    '''
    The following are helper methods that do the actual work of setup and extraction.
    Helper methods are not required by the interface but may be implemented by the developer
    to make the module more organized and readable.
    '''
            
    def set_element_list(self):
        '''
        Keeps track of the elements in the structure.
        This is because LAMMPS uses indices instead of strings to specify elements.
        '''
        element_list = []
        for i in range(self.input_structure.get_n_atoms()):
            if self.input_structure.geometry[i]['element'] not in element_list:
                element_list.append(self.input_structure.geometry[i]['element'])
#        element_list.sort()
        self.element_list = element_list
        
	#SOME OF THESE FUNCTIONS PULLED FROM PATRICKS SCRIPTS

    def leng(self,v):
        return((v[0]**2+v[1]**2+v[2]**2)**0.5)	

    def angle(self,v1,v2):
    	numdot = np.dot(v1,v2)
	anglerad = np.arccos(numdot/(self.leng(v1)*self.leng(v2))) 
	angledeg = anglerad*180/np.pi
	return (angledeg)

    def create_geometry_in(self):
        '''
        Writes the data (geometry.xyz) file called by TINKER. Builds with no template.
        Format is very particular. atom types and number change for change in parameter set.
        Returns: None
        '''
	geo = self.input_structure.geometry
	atomnum = [i+1 for i in range(len(geo))]
	molnum = 2
	
	repeat = len(geo)/molnum
	ambertypes = ['C', 'CT','H1', 'H1', 'N3', 'H', 'H','H', 'O2', 'O2','C', 'CT','H1', 'H1', 'N3', 'H', 'H','H', 'O2', 'O2']
	ambernum = [503, 351, 355, 355, 350, 353, 353, 353, 505, 505, 503, 351, 355, 355, 350, 353, 353, 353, 505, 505 ]
	connect1 = [2, 1, 2, 2, 2, 5, 5, 5, 1, 1,  12, 11, 12, 12, 12, 15, 15, 15, 11, 11]
	connect2 = [' ']*len(geo)
	connect3 = [' ']*len(geo)
	connect4 = [' ']*len(geo)
	connect5 = [' ']*len(geo)

	connect2[0] = 9
	connect3[0] = 10
 
	connect2[1] = 3
	connect3[1] = 4
	connect4[1] = 5

	connect2[4] = 6
	connect3[4] = 7
	connect4[4] = 8

	connect2[0+repeat] = connect2[0]+repeat
	connect3[0+repeat] = connect3[0]+repeat
	
	connect2[1+repeat] = connect2[1]+repeat 
        connect3[1+repeat] = connect3[1]+repeat 
        connect4[1+repeat] = connect4[1]+repeat

	connect2[4+repeat] = connect2[4]+repeat 
        connect3[4+repeat] = connect3[4]+repeat
        connect4[4+repeat] = connect4[4]+repeat

	lata = np.asarray(self.input_structure.get_property('lattice_vector_a'))
	latb = np.asarray(self.input_structure.get_property('lattice_vector_b'))
        latc = np.asarray(self.input_structure.get_property('lattice_vector_c'))

	a = self.leng(lata)
	b = self.leng(latb)
	c = self.leng(latc)	
	alpha = self.angle(latb,latc)
	beta = self.angle(lata,latc)
	gamma = self.angle(lata,latb)

        data_string = ''
	data_string += str(geo.size) + ' atoms (B glycine crystal 2 mols per cell)\n'
	data_string += str(a) + ' ' + str(b) + ' ' + str(c) + ' ' +str(alpha) + ' ' + str(beta) + ' ' +str(gamma)
	data_string += '\n'
	for j in range(len(ambertypes)):
		data_string += str(j+1) + ' ' + str(ambertypes[j])  + ' ' + str(geo[j]['x']) + ' ' + str(geo[j]['y']) + ' ' +  str(geo[j]['z']) + ' ' +str(ambernum[j])+ ' ' +str(connect1[j])+ ' ' + str(connect2[j]) + ' ' +str(connect3[j])+ ' '+ str(connect4[j]) 
                data_string += '\n'

        data_file = open(os.path.join(self.working_dir, 'geometry_amber.xyz'), 'w')
	data_file.write(data_string)
        data_file.close()
       
    def create_input(self):
        '''
        Writes the input .key file required by TINKEr
        Returns: None
        '''
        input_file = open(os.path.join(self.working_dir, 'geometry_amber.key'), 'w')
        input_file.write(self.control_in_string)
        input_file.close()

    def extract_geometry(self):
        '''
        Reads the dump file in xyz format and builds the structure's geometry
        Returns: True if successful
        '''
        # reads dump file from tinker output
        dump = open(os.path.join(self.working_dir, 'geometry_amber.xyz_2'), 'r')
        dump_list = []
	dump_array = []
	for line in dump:
            dump_list.append(line)

	# Reads xyz coordinates from dump file
	element_lst = ['C','C','H','H','N','H','H','H','O','O','C','C','H','H','N','H','H','H','O','O']	
	for i in range(len(self.input_structure.geometry)):
            xyz = dump_list[i+2].split()
            dump_array +=[[xyz[2],xyz[3],xyz[4], element_lst[i]],]
        for i in range(len(self.input_structure.geometry)):
                 self.result_struct.build_geo_by_atom(float(dump_array[i][0]),
                                                 float(dump_array[i][1]),
                                                 float(dump_array[i][2]),
                                                 dump_array[i][3])

	#Reads lattice vector components and angles
	[a, b, c, alpha, beta, gamma] = dump_list[1].split()
	
	a=float(a); b=float(b); c=float(c); 
	rad=float(np.pi/180)
	alpha=float(alpha); beta=float(beta); gamma=float(gamma)

	ax=a;	ay=0;	az=0
	bx=np.cos(gamma*rad)*b;	by=np.sin(gamma*rad)*b;	bz=0
	cx=c*np.cos(beta*rad);	cy=(b*c*np.cos(alpha*rad)-bx*cx)/by;	cz=np.sqrt(np.absolute(c**2-cx**2-cy**2))
	lata_out = np.zeros(3);	latb_out = np.zeros(3);	latc_out = np.zeros(3)

	lata_out[0] = ax;	lata_out[1] = ay;	lata_out[2] = az
	latb_out[0] = bx;	latb_out[1] = by;	latb_out[2] = bz
	latc_out[0] = cx;	latc_out[1] = cy;	latc_out[2] = cz

	cell_vol = np.dot(lata_out, np.cross(latb_out, latc_out))
        self.result_struct.set_property('lattice_vector_a',lata_out.tolist())
        self.result_struct.set_property('lattice_vector_b',latb_out.tolist())
        self.result_struct.set_property('lattice_vector_c',latc_out.tolist())	

	self.result_struct.set_property('a', a)
        self.result_struct.set_property('b', b)
        self.result_struct.set_property('c', c)
	self.result_struct.set_property('alpha', alpha)
        self.result_struct.set_property('beta', beta)
        self.result_struct.set_property('gamma', gamma)
	self.result_struct.set_property('cell_vol', cell_vol)

	ct = self.input_structure.get_property('crossover_type')
	mt = self.input_structure.get_property('mutation_type')
	self.result_struct.set_property('crossover_type',ct)
	self.result_struct.set_property('mutation_type',mt)

#	print self.result_struct.get_lattice_vectors()
#	print self.result_struct.geometry()
        return True
    
    def extract_energy(self):
        '''
        reads the resulting energy from the tinker log in the tmp dir
        Returns: float if successful, False if unsuccessful
        '''
        logtinker = open(os.path.join(self.working_dir, 'log.tinker'), 'r')
		
	while True:
		line = logtinker.readline().split()
		if 'Potential' in line:
			print "Energy:  ", line[5]
			self.result_struct.set_property('energy',float(line[5]))
			return float(line[5])  # final energy


