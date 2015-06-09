'''
Created on Jul 29, 2013

@author: newhouse
'''
import os
import subprocess
import time
import numpy as np
from core import user_input, output
from core.file_handler import mkdir_p_clean
from structures.structure import Structure
from copy import deepcopy


def main(input_structure, working_dir, control_in_string, replica):
    '''
    For a relaxation module, the method main() must be implemented.
    
    It must take as arguments a 'Structure', a string which defines 
    in which working directory the relaxation will take place, and a string
    representation of the full control.in file
    
    It must return a Structure object. (Please see structure documentation
    for information on how to easily setup a Structure object) 
    '''
    r = FHIAimsRelaxation(input_structure, working_dir, control_in_string, replica)
    print "setting up FHI-aims input files"	
    r.setup()
    print "executing FHI-aims"	
    r.execute()
    if r.is_successful(): return r.extract()
    else: 
	print "Relaxation failed!"
#	raise RuntimeError('relaxation failed!')
	return False

class FHIAimsRelaxation():
    '''
    Creates an optimizer object using FHI-Aims to execute a geometry.
    
    This is an example of using a class to define a state, in this case the
    state of the relaxation process throughout setup, execution, and extraction
    '''
    def __init__(self , input_structure, working_dir, control_in_string, replica):
        '''
        __init__ runs when class is initialized. 
        '''
        # this gives easy access to retrieve preferences from the user_input file
        self.ui = user_input.get_config()
        self.replica = replica
        # defining 'self' fields gives access to these values throughout the class 
        self.input_structure = input_structure
        self.working_dir = working_dir
        self.control_in_string = control_in_string
    
    def output(self, message): output.local_message(message, self.replica)    

    def setup(self):
        '''
        Creates the environment for the optimizer and constructs the input files.
        In this case, FHI-aims will take a geometry.in file and an comtrol.in file.
        '''
        mkdir_p_clean(self.working_dir)
        self.create_geometry_in()
        self.create_control_in()
        
    def create_geometry_in(self):
        '''
        Writes the geometry.in file called by FHI-aims. 
        Calls from FHI-aims control settings found in res directory.
        Returns: None
        '''
        geometry_file = open(os.path.join(self.working_dir, 'geometry.in'), 'w')
        geometry_file.write(self.input_structure.get_geometry_atom_format())
        geometry_file.close()
        
    def create_control_in(self):
        '''
        Writes the control.in input file required by FHI-aims.
        Returns: None
        '''
        control_file = open(os.path.join(self.working_dir, 'control.in'), 'w')
        control_file.write(self.control_in_string)
        control_file.close()

    def execute(self):
        '''
        Creates an executable shell script with proper FHI-aims syntax 
        then executes in a subprocess
        '''
	out_location = str(self.working_dir)
	exe_string = """#!/bin/bash\n
#SBATCH --qos=normal\n
#SBATCH --job-name=fhi_aims\n

#SBATCH --time=0-2:00:00\n
#SBATCH -o aims.log
#SBATCH -e aims.err
#SBATCH --nodes=12\n
#SBATCH --ntasks-per-node=20\n
#SBATCH --workdir="""
	exe_string += out_location 

	exe_string +="""\nexport OMP_NUM_THREADS=1\n
export MKL_NUM_THREADS=1\n
export MKL_DYNAMIC=FALSE\n
module load intel/mic/runtime/3.3\n
module load intel/mic/sdk/3.3\n
module load tulane/intel-psxe/2015\n"""
        exe_string += 'mpirun /home/xli20/fhi-aims/fhi-aims.071914/bin/aims.071914_1.scalapack.mpi.x > ' + out_location+ '/aims.out'
	exe_script = open(os.path.join(self.working_dir, 'exe.sh'), 'w')
        exe_script.write(exe_string)
	exe_script.close()
#	os.system('sbatch ' + os.path.join(self.working_dir, 'exe.sh'))
        while True:
    	    p = subprocess.Popen(['sbatch','exe.sh'],shell=False,cwd=self.working_dir,stdout=subprocess.PIPE)
            p.wait()
	    out,err=p.communicate()
	    print out
	#print "in fhi_relaxation_cypress, this is what is captured by screen", out	
  	    if out.split()[0]!='Submitted' or out.split()[1]!='batch':
		raise RuntimeError("Is the job run on Cypress?")
	    self.jobid=out.split()[3]
	    time.sleep(1)
	    p = subprocess.Popen(['squeue --job '+self.jobid],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            p.wait()
            out,err=p.communicate()
            try:
                    l=out.split()[12]
            except: #When a job is done, squeue will no longer return the information about the job
                    print "Job submission failed! ", self.jobid
                    continue
	    break
	self.wait_on_job()
	time.sleep(2)

    def extract(self):
        '''
        Reads the output files from relaxation and specifies properties of the structure (e.g. energy)
        Returns: Structure or False
       	'''
        # creates empty structure
        self.result_struct = Structure()
        #try: self.result_struct.set_lattice_vectors(self.input_structure.get_lattice_vectors())

	# reads output file to extract energy
        energy = self.extract_energy()
	print "energy", energy
        if energy == False: return False  # did not converge
        self.result_struct.set_property('energy', energy)
       
	# sets the result structure's lattice vectors based on output file
	lats= self.extract_lats()

	latA = [float(lats[0][0]), float(lats[0][1]), float(lats[0][2])]
        latB = [float(lats[1][0]), float(lats[1][1]), float(lats[1][2])]
        latC = [float(lats[2][0]), float(lats[2][1]), float(lats[2][2])]
	self.result_struct.set_property('lattice_vector_a', latA)
        self.result_struct.set_property('lattice_vector_b', latB)
        self.result_struct.set_property('lattice_vector_c', latC)
	
	latA = np.asarray(latA)
	latB = np.asarray(latB)
	latC = np.asarray(latC)

	temp_vol = np.dot(np.cross(latA, latB), latC)
        alpha = self.angle(latB, latC)
        beta = self.angle(latA, latC)
        gamma = self.angle(latA, latB)
	self.result_struct.set_property('cell_vol', temp_vol)
	self.result_struct.set_property('alpha',alpha)
        self.result_struct.set_property('beta', beta)
        self.result_struct.set_property('gamma', gamma)

	print "latA",latA
	print "latB", latB
	print "latc", latC
	print "alpha", alpha
	print "beta", beta
	print "gamma", gamma
 
        # sets the result structure's geometry based on output file
        extract_success = self.extract_geometry()

        # extract the spin moment
        if self.ui.get('FHI-aims','initial_moment') == 'custom': self.extract_spin_moment()
        
        return self.result_struct

    def extract_energy(self):
        '''
        reads the resulting energy from the FHI-aims log
        Returns: float if successful, False if unsuccessful
        '''
        aims_out = open(os.path.join(self.working_dir, 'aims.out'))
        while True:
            line = aims_out.readline()
            if not line: return False  # energy not converged
            if '  | Total energy corrected        :' in line:
                tokens = line.split()
                energy = float(tokens[5])  # converts from SI string to float
                return energy

    def extract_lats(self):
	print "extract lats..."
	aims_out = open(os.path.join(self.working_dir, 'aims.out'))
	while True:
        	line = aims_out.readline()
        	if 'Final atomic structure' in line: break
        aims_out.readline()  # skip one line
        atom_string = ''
	while True:
        	line = aims_out.readline()
        	if 'atom' in line: break
        	else: atom_string += line

	latout = atom_string.split()

#	print latout
#	lat_A = np.asarray(latout[0]),latout[1], float(latout[2]))
#	lat_B = np.asarray(float(latout[5]),float(latout[6]), float(latout[7]))
#	lat_C = np.asarray(float(latout[9]),float(latout[10]), float(latout[11]))
	lat_A = [latout[1], latout[2], latout[3]]
	lat_B = [latout[5], latout[6], latout[7]]
	lat_C = [latout[9], latout[10], latout[11]]
	lats = [lat_A, lat_B, lat_C]

	return lats

    def angle(self,v1,v2):
        numdot = np.dot(v1,v2)
        anglerad = np.arccos(numdot/(self.leng(v1)*self.leng(v2)))
        angledeg = anglerad*180/np.pi
        return (angledeg)
   
    def leng(self,v):
        length = np.linalg.norm(v)
        return length
  
    def extract_geometry(self):
        '''
        Reads the FHI-aims output file and builds the structure's geometry
        Returns: True if successful
        '''
        aims_out = open(os.path.join(self.working_dir, 'aims.out'))
        while True:
            line = aims_out.readline()
            if not line: return False  # not converged
            if 'Final atomic structure' in line: break
        aims_out.readline()  # skip one line
        aims_out.readline()  # skip one line
        aims_out.readline()  # skip one line
        aims_out.readline()  # skip one line
        atom_string = ''
        while True:
            line = aims_out.readline()
            if 'Fractional coordinates:' in line: break
            else: atom_string += line
        self.result_struct.build_geo_whole_atom_format(atom_string)
        return True
    
    def extract_spin_moment(self):
        try:
            os.system("grep '|   Hirshfeld spin moment' aims.out > spin_moment.out")
            spin_file = open('spin_moment.out', 'r')
            total_spin_list = spin_file.read()
            spin_file.close()
        except: return False
        
        reduced_spin_list = total_spin_list[-len(self.result_struct.get_geometry()):]
        for i in range(len(reduced_spin_list)):
            self.result_struct.geometry[i]['spin'] = float(reduced_spin_list[i].split()[-1])
        return True
        
    def extract_charge(self):
        try:
            os.system("grep '|   Hirshfeld charge' aims.out > charge.out")
            charge_file = open('charge.out', 'r')
            total_charge_list = charge_file.read()
            charge_file.close()
        except: return False
        
        reduced_charge_list = total_charge_list[-len(self.result_struct.get_geometry()):]
        for i in range(len(reduced_charge_list)):
            self.result_struct.geometry[i]['charge'] = float(reduced_charge_list[i].split()[-1])
        return True
    
    def is_successful(self):
        '''
        checks if relaxation/optimization was successful
        '''
        # TODO: need proper check
#	aims_out = open(os.path.join(self.working_dir, 'aims.out'))
        aims_out = os.path.join(self.working_dir, 'aims.out')
#	self.wait_on_file(aims_out)
	aims_out = open(aims_out)
        while True:
            line = aims_out.readline()
	    if line == '':
#		print "In is_successful, still waiting for aims to finish output."
#	    else:
#		print "Read: ", line
		break
        #    if not line: return False  # energy not converged
            if 'Have a nice day' in line:
                return True
        return False

    def wait_on_file(self,wait_file,sleeptime=1):
	"""
   	Simple function that sleeps until a file path exists and any writes
   	to it finish.
   	"""
	lasttime=0
#	print "Waiting for the file!"
   # waiting for the file to exist
   	while not os.path.exists(wait_file):
#		print "Still waiting!"
       		time.sleep(sleeptime)
   # wait for the file to stop updating
#	print "Waiting for the file to be done ", os.path.getmtime(wait_file)
   	while (os.path.getmtime(wait_file)-lasttime!=0):
#		print "Still waiting! Last edited", os.path.getmtime(wait_file)
       		lasttime=os.path.getmtime(wait_file)
       		time.sleep(sleeptime)

    def wait_on_job(self,sleeptime=1):
	'''
	Waits until the submitted job is finished
	'''
	time.sleep(sleeptime*5)
	while True:
		p = subprocess.Popen(['squeue --job '+self.jobid],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
       		p.wait()
		time.sleep(sleeptime)
        	out,err=p.communicate()
		try:
			l=out.split()[12]
		except: #When a job is done, squeue will no longer return the information about the job
			print "FHI-aims Job finished! Jobid:", self.jobid
			break
		time.sleep(sleeptime)
	
