'''
Created on Jul 29, 2013
@authors: newhouse/patrick/farren
'''
import os
import subprocess
import time
import datetime
import numpy as np
import shutil
from core import user_input, output
from core.file_handler import mkdir_p_clean, fail_dir
from core.activity import *

from structures.structure import Structure
from copy import deepcopy


def main(input_structure, working_dir, control_check_SPE_string, control_relax_full_string, replica):
	'''	
	It must take as arguments a 'Structure', a string which defines 
	in which working directory the relaxation will take place, and a string
	representation of the full control.in file
	
	It must return a Structure object. (Please see structure documentation
	for information on how to easily setup a Structure object) 
	'''
	ui = user_input.get_config()
	SPE_threshold = ui.get_eval('run_settings','SPE_threshold')

	SPE = FHIAimsRelaxation(input_structure, working_dir, control_check_SPE_string, replica)
	output.local_message("Replica " +str(replica)+" executing FHI-aims and checking the SPE.", replica)
	checkSPE = SPE.execute()
	if SPE.is_successful_spe():
		SP_energy = SPE.extract_energy()
		input_structure.set_property('spe_energy',SP_energy)
		output.local_message("SPE extracted"+str(SP_energy), replica)
		if SP_energy >= SPE_threshold:
			return False
		elif SP_energy <= SPE_threshold:
			message = "Structure passed SPE check"
			output.local_message(message, replica)
			fullrelax = FHIAimsRelaxation(input_structure, working_dir, control_relax_full_string, replica)
			attempts = 3
			for i in range (attempts):
				output.local_message("Replica " +str(replica)+" executing FHI-aims for a full unit cell relax ; attempt "+str(i),replica)
				output.time_log("Attempt %i for full unit cell relax" % i,replica)
				fullrelax.execute()
				if fullrelax.is_successful(): 
					return fullrelax.extract()
		       		output.local_message("Relaxation failed!",replica)
				output.time_log("aims relaxation determined as failed",replica)
				fdir = os.path.abspath(os.path.join(working_dir,os.pardir))
		                if bk_folder(fdir,working_dir[len(fdir)+1:],fail_dir,"random"):
		                        output.local_message("Successfully backed up failed folder",replica)
		                else:
                		        output.local_message("Failed to back up failed folder",replica)

                		if i==attempts-1:
					output.time_log("aims job maxed out on attempts",replica)
					return False
				if os.path.isfile(os.path.join(working_dir,"geometry.in.next_step")):
					output.time_log("geometry.in.next_step found",replica)
					try:
						os.remove(os.path.join(working_dir,"geometry.in"))
					except:
						output.time_log("WARNING! geometry.in file missing for some reason",replica)
					os.rename(os.path.join(working_dir,"geometry.in.next_step"),os.path.join(working_dir,"geometry.in"))

				
		else:
			print "Error with energy threshold"
	return False

class FHIAimsRelaxation():
    '''
    Creates an optimizer object using FHI-Aims to execute a geometry.
    
    This is an example of using a class to define a state, in this case the
    state of the relaxation process throughout setup, execution, and extraction
    '''
    def __init__(self , input_structure, working_dir, control_in_string, replica):
        '''
	Creates the environment for the optimizer and constructs the input files.
	In this case, FHI-aims will take a geometry.in file and an comtrol.in file.
        '''
        # this gives easy access to retrieve preferences from the user_input file
        self.ui = user_input.get_config()
        self.replica = replica
        # defining 'self' fields gives access to these values throughout the class 
        self.input_structure = input_structure
        self.working_dir = working_dir
        self.control_in_string = control_in_string
	mkdir_p_clean(self.working_dir)
	self.create_geometry_in()
	self.create_control_in()
	self.set_permission()
    
    def output(self, message): output.local_message(message, self.replica)   
        
    def create_geometry_in(self):
        '''
        Writes the geometry.in file called by FHI-aims. 
        Calls from FHI-aims control settings found in res directory.
        Returns: None
        '''
        geometry_file = open(os.path.join(self.working_dir, 'geometry.in'), 'w')
        geometry_file.write(self.input_structure.get_geometry_atom_format())
        geometry_file.close()
	more_file = open(os.path.join(self.working_dir, "struct.json"),"w")
	more_file.write(self.input_structure.dumps())
	more_file.close()
        
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
	Directly calls mpirun to run the aims executable        
	'''
	out_location = str(self.working_dir)
        ui=user_input.get_config()
        bin=ui.get('FHI-aims','path_to_aims_executable')
	self.bin=bin
	environment=ui.get('parallel_settings','system')
	output.time_log("Beginning aims execution",self.replica)
	output.local_message("Aims relaxation being called. out_location=%s" % (out_location),self.replica)
	output.local_message("Binary location is"+bin,self.replica)

	if environment=='Cypress_login' or environment=="cypress_login" or environment=="Cypress-login" or environment=="cypress-login" or environment=="Cypress" or environment=="cypress":
		arglist=["mpirun","-wdir",self.working_dir,self.bin]

	elif environment=='Edison_login':
		nodes = ui.get_eval('parallel_settings','nodes_per_replica')
                ppn_edison = 24
                width = ppn_edison*nodes
		arglist=["aprun","-n",str(width),self.bin]
		current_dir=os.getcwd()
		os.chdir(out_location) #aprun does not allow -wdir option

	elif environment=="cetus" or environment=="Cetus" or environment=="mira" or environment=="Mira":
		block_size=ui.get_eval('parallel_settings','nodes_per_replica')
		#Will run it with modes=4 and thre=4
		modes=4; thres=4
		try: #Replica name is defined as block%core%shape%random_index
			l=self.replica.index("%")
			block=self.replica[0:l]
			rest=self.replica[l+1:]
			l=rest.index("%")
			corner=rest[:l]
			rest=rest[l+1:]
			l=rest.index("%")
			shape=rest[:l]
			arglist=["runjob","--np",str(modes*block_size),"-p",str(modes),"--envs","OMP_NUM_THREADS="+str(thres),"--verbose","INFO","--block",block,"--corner",corner,"--shape",shape,"--cwd",self.working_dir,"--exe",self.bin]
		except: #Only has a block name
			l=self.replica.index("%")
			block=self.replica[0:l]
			arglist=["runjob","--np",str(modes*block_size),"-p",str(modes),"--envs","OMP_NUM_THREADS="+str(thres),"--verbose","INFO","--block",block,"--cwd",self.working_dir,"--exe",self.bin]

	aimsout=os.path.join(self.working_dir,"aims.out")
	for i in range (10):
		outfile=open(aimsout,"w")
		get_execute_clearance(request_folder=self.working_dir)
		output.time_log("aims job execute clearance acquired",self.replica)
		p=subprocess.Popen(arglist,stdout=outfile)
		time.sleep(1)
		try:
			status=p.poll()
		except: #OSError Errno 3 Process does not exist
			output.time_log("Nodes failure ; replica will pass out from now on")
			time.sleep(86400)
			return False
			
		time_limit=60
		for j in range (time_limit): #Allow 60 seconds for aims to start outputting
			if (p.poll()!=None) or (os.stat(aimsout).st_size>512):
				break
			write_active(self.working_dir)
			self.set_permission()
			time.sleep(1)
		if (os.stat(aimsout).st_size>512):
			output.time_log("aims.out begins output", self.replica)
			break
		outfile.close()
		output.time_log("aims job launch failure",self.replica)
		try:
			p.send_signal(2)
		except:
			output.time_log("Unable to kill process ; possible node failures", self.replica)
			time.sleep(86400)
		active_sleep(60,self.working_dir)
		try:
			self.set_permission()
		except:
			pass


		if i==9:
			output.time_log("WARNING: Repeated launch failure ; exiting",self.replica)
			return False
	counter=0; last=os.stat(aimsout).st_size
	while counter<60 and p.poll()==None: #The output file needs to update at least once in every 5 minutes
		write_active(self.working_dir)
		self.set_permission()
		time.sleep(10)
		if os.stat(aimsout).st_size>last:
			last=os.stat(aimsout).st_size
			counter=0
		else:
			counter+=1
	if counter==60:
		output.time_log("aims job hung",self.replica)
		try:
			p.send_signal(2)
		except:
			output.time_log("Unable to kill process ; possible node failures", self.replica)
			time.sleep(86400)
		active_sleep(60,self.working_dir)
		try:
			self.set_permission()
		except:
			pass

			
	outfile.close()
	output.local_message("aims job exited with status "+str(p.poll()),self.replica)
	output.time_log("aims job exited with status "+str(p.poll()),self.replica)
	
	if environment=="Edison_login": #change it back to the current directory
		os.chdir(current_dir)

    def output(self, message): output.local_message(message, self.replica)

    def is_successful(self):
        '''
        checks if relaxation/optimization was successful
        '''
        aims_path = os.path.join(self.working_dir, 'aims.out')
        aims_out = open(aims_path,"r")
        counter = 0
        while True:
                line = aims_out.readline()
                if "Have a nice day" in line:
			return True
                elif line == '':
                        counter += 1
		else:
			counter = 0
                if counter > 10: #Allowing 10 empty lines in a row before determining eof
                        break
	return False



    def is_successful_2(self):
	fileBytePos = 0
	aims_out = os.path.join(self.working_dir, 'aims.out')
	file = open(aims_out,"r")
	while 1:
		#file = open(aims_out,"r")
    		where = file.tell()
    		line = file.readline()
    		if not line:
        		time.sleep(30)
        		file.seek(where)
			line = file.readline()
			if not line:
				self.output("Output is stuck")
    		else:
        		self.output(line) # already has newline



    def is_successful_spe(self):
        '''
        checks if relaxation/optimization was successful
        '''
        aims_out = os.path.join(self.working_dir, 'aims.out')
        aims_out = open(aims_out)
        while True:
            line = aims_out.readline()
            if line == '':
	    	break
            if 'Leaving FHI-aims.' in line:
                return True
        return False

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
	self.output("extracted energy" +str(energy))
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
	a = np.linalg.norm(latA)
	b = np.linalg.norm(latB)
	c = np.linalg.norm(latC)
	self.result_struct.set_property('cell_vol', temp_vol)
	self.result_struct.set_property('alpha',alpha)
        self.result_struct.set_property('beta', beta)
        self.result_struct.set_property('gamma', gamma)
	self.result_struct.set_property('a',a)
        self.result_struct.set_property('b', b)
        self.result_struct.set_property('c', c)

        ct = self.input_structure.get_property('crossover_type')
	mt = self.input_structure.get_property('mutation_type')
	spe = self.input_structure.get_property('spe_energy')
	self.result_struct.set_property('mutation_type',mt)
        self.result_struct.set_property('crossover_type',ct)
	self.result_struct.set_property('spe_energy',spe)
	
        # sets the result structure's geometry based on output file
        extract_success = self.extract_geometry()

	# extracts time from aims output
	wall_time = self.extract_time()
	self.result_struct.set_property('relax_time', wall_time)
        # extract the spin moment
        if self.ui.get('FHI-aims','initial_moment') == 'custom': self.extract_spin_moment()
        
        return self.result_struct

    def extract_time(self):
        '''
        reads the resulting time from the FHI-aims log
        Returns: float if successful, False if unsuccessful
        '''
        aims_out = open(os.path.join(self.working_dir, 'aims.out'))
        while True:
            line = aims_out.readline()
            if not line: return False  # energy not converged
            if '| Total time                                 :' in line:
                tokens = line.split()
                time = tokens[6]  # converts from SI string to float
                return time

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

    def set_permission(self):
	'''
	Allow other group user to access this folder
	'''
	os.system("chmod -R g=u "+self.working_dir)
    
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


