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
import core.file_handler as fh
from core.activity import *
from utilities import misc, parallel_run

from structures.structure import Structure
from structures.structure_collection import StructureCollection, get_collection
from copy import deepcopy

ui = user_input.get_config()

def main(input_structure):
	'''	
	It must take as arguments a 'Structure'.
	It must return a Structure object, and a status. 

	Status includes:
	"failed": structure failed current relaxation
	"rejected": structure rejected due to certain threshold
	"accepted": structure accepted and should be further processed
	'''
	sname = "FHI-aims"
	input_structure = deepcopy(input_structure)
	working_dir = os.path.join(fh.tmp_dir,ui.get_replica_name())
	control_list = ui.get_list(sname,"control_in_filelist")
	control_dir = ui.get(sname,"control_in_directory")
	abs_success = ui.get_boolean(sname,"absolute_success")
	if ui.get_boolean(sname,"monitor_execution"):
		monitor = True
		upt = ui.get_list(sname,"update_poll_times",eval=True)
		if len(upt)!=len(control_list):
			raise ValueError("Number of specified update poll interval must match that of control.in files")
		upi = ui.get_list(sname,"update_poll_interval",eval=True)
		if len(upi)!=len(control_list):
			raise ValueError("Number of specified update poll times must match that of control.in files")
	else:
		monitor = False
		upt = [None]*len(control_list)
		upi = [None]*len(control_list)

	if ui.has_option(sname,"absolute_energy_thresholds"):
		at = ui.get_list(sname,"absolute_energy_thresholds",eval=True)
		if len(at)!=len(control_list):
			raise ValueError("Number of absolute energy thresholds not matching number of control.in files specified")
	else:
		at = [None]*len(control_list)
	if ui.has_option(sname,"relative_energy_thresholds"):
		rt = ui.get_list(sname,"relative_energy_thresholds",eval=True)
		if len(rt)!=len(control_list):
			raise ValueError("Number of relative energy thresholds not matching number of control.in files specified")
	else:
		rt = [None]*len(control_list)

	if ui.has_option(sname,"reject_if_worst_energy"):
		worst_energy = ui.get_list_of_booleans(sname,"reject_if_worst_energy")
		if len(worst_energy)!=len(control_list):
			raise ValueError("Number of reject_if_worst_energy boolean flags not matching number of control.in files specified")
	else:
		worst_energy = [False]*len(control_list)

	if ui.has_option(sname,"store_energy_names"):
		sen = ui.get_list(sname,"store_energy_names")
		if len(sen)!=len(control_list):
			raise ValueError("Number of store energy names not matching number of control.in files specified")
	else:
		sen = ["energy"]*len(control_list)

	if ui.has_option(sname,"relative_energy_names"):
		ren = ui.get_list(sname,"relative_energy_names")
		if len(ren)!=len(control_list):
			raise ValueError("Number of relative energy names not matching number of control.in files specified")
	else:
		ren = ["energy"]*len(control_list)

	output.local_message("\n|-------------------------- FHI-aims evaluation --------------------------|")
	stoic = input_structure.get_stoic()

	#####################Beginning cascade#####################
	for i in range(len(control_list)):
		if rt[i]!=None or worst_energy[i]:
			struct_coll = get_collection(stoic,0)
			struct_coll.update_local()
			energies = []
			for key, struct in struct_coll:
				energies.append(struct.get_property(ren[i]))
		
		et = None #energy threshold
		if rt != None and rt[i]!=None: #Relative threshold
			if ren!=None:
				ekey = ren[i]
			else:
				ekey = "energy"
			et = min(energies)+rt[i]

		if at != None and at[i]!=None: #Absolte threshold
			if et == None:
				et = at[i]
			else:
				et = min(et,et[i])

		if worst_energy[i]: #Reject if worst energy:
			struct_coll = get_collection(stoic,0)
			struct_coll.update_local()
			if et == None:
				et = max(energies)
			else:
				et = min(et,max(energies))

		control_path = os.path.join(os.path.join(fh.cwd,control_dir),
					    control_list[i])
		control_string = fh.read_data(control_path)

		FHI = FHIAimsRelaxation(input_structure, 
					working_dir, 
					control_string,
					ui.get_replica_name())

		output.time_log("Executing FHI-aims using control file: "
				+ control_list[i])

		output.local_message("\nFHI-aims evaluation using control file: "
					+ control_list[i])
		output.local_message("-- Rejection energy threshold: "+str(et)+" eV")

		begin_time = time.time()
		FHI.execute(monitor,upi[i],upt[i])
		end_time = time.time()

		output.local_message("-- Job execution time: " + 
					str(end_time-begin_time) + " seconds")

		#Need to pass in update poll interval and update poll time
		if abs_success:
			execute_success = FHI.is_successful_absolute()
		else:
			execute_success = FHI.is_successful_relative()


		if not execute_success: 
			output.time_log("FHI-aims execution using control file %s failed to launch, hung, or failed" % (control_list[i]))
			output.local_message("-- Job failed to launch, hung, or failed")
			if ui.get_boolean(sname,"save_failed_calc"):
				path = os.path.abspath(os.path.join(fh.fail_dir,
				misc.get_random_index()))
				shutil.copytree(working_dir,path)
				output.time_log("Failed calc folder saved to "+path)
				output.local_message("-- Failed calc folder saved to "+path)
			return input_structure, "failed"


		output.time_log("Execution with control file " + control_list[i] + " completed.")
		extract_result = FHI.extract()
		if extract_result == "failed":
			output.time_log("FHI-aims execution using control file %s failed to yield energy" % (control_list[i]))
			output.local_message("-- Job failed to yield energy")
			if ui.get_boolean(sname,"save_failed_calc"):
				path = os.path.abspath(os.path.join(fh.fail_dir,
				misc.get_random_index()))
				shutil.copytree(working_dir,path)
				output.time_log("Calculation folder saved to "+path)
				output.local_message("Calculation folder saved to "+path)
			return input_structure, "failed"

			
		#Update energy
		energy = FHI.result_struct.properties["energy"]
		input_structure.set_property(sen[i],energy)	
		output.local_message("-- Energy extracted: "+str(energy)+ " eV")

		#Update geometry if necessary
		if extract_result == "relaxed": 
			output.local_message("-- Updated geometry retrieved")
			if ui.all_geo():
				output.local_message(FHI.result_struct.\
				get_geometry_atom_format())
			input_structure.geometry = FHI.result_struct.geometry
			properties = ["lattice_vector_a","lattice_vector_b",\
		"lattice_vector_c","cell_vol","a","b","c","alpha","beta","gamma"]
			for prop in properties:
				if prop in FHI.result_struct.properties:
					input_structure.properties[prop] = \
					FHI.result_struct.properties[prop]
		
		if ui.get_boolean(sname,"save_successful_calc"):
			path = os.path.abspath(os.path.join(fh.success_dir,
			misc.get_random_index()))
			shutil.copytree(working_dir,path)
			output.time_log("Successful calc folder saved to "+path)
			output.local_message("-- Successful calc folder saved to "+path)

		if et!=None and energy >= et: #Rejected
			output.local_message("-- Energy threshold not met; structure rejected")
			return input_structure,"rejected"
		else:
			output.local_message("-- Energy threshold met")
        output.local_message('\n|------------------------ End FHI-aims evaluation ------------------------|')			
	return input_structure,"accepted"

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
	fh.mkdir_p_clean(self.working_dir)
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

    def execute(self,enable_monitor=False,update_poll_interval=None,update_poll_times=None):
        '''
	Directly calls mpirun to run the aims executable        
	'''
	def end_of_execution_tasks():
		outfile.close()
		output.local_message("-- aims job exit status: "
					+ str(p.poll()),self.replica)
		output.time_log("aims job exited with status "
				+ str(p.poll()),self.replica)
		os.chdir(original_dir)

	if enable_monitor and (update_poll_interval==None or update_poll_times==None):
		raise ValueError("FHI-aims job monitoring enabled, but no poll interval or times specified")

	out_location = str(self.working_dir)
        ui=user_input.get_config()
        bin=ui.get('FHI-aims','path_to_aims_executable')
	self.bin=bin
	execute_command = ui.get('FHI-aims','execute_command')

#	output.local_message("Aims relaxation being called. out_location=%s" % (out_location),self.replica)
#	output.local_message("Binary location is"+bin,self.replica)
	sname = "parallel_settings"

	original_dir = os.getcwd()
	if execute_command == "mpirun":
		arglist = ["mpirun","-wdir",self.working_dir]
		if ui.has_option("parallel_settings","allocated_nodes"):
			arglist += ["-host",",".join(map(str,ui.get_eval("parallel_settings","allocated_nodes")))]
		if ui.has_option("parallel_settings","processes_per_replica"):
			arglist += ["-n",ui.get("parallel_settings","processes_per_replica")]	
		if ui.has_option("parallel_settings","additional_arguments"):
			arglist += ui.get_eval("parallel_settings","additional_arguments")
		arglist += [self.bin]


	elif execute_command == "srun":
		arglist = ["srun","-D",self.working_dir]
		if ui.has_option("parallel_settings","allocated_nodes"):
			arglist += ["-w",",".join(map(str,ui.get_eval("parallel_settings","allocated_nodes")))]
			arglist += ["-N",str(len(ui.get_eval("parallel_settings","allocated_nodes")))]
		if ui.has_option("parallel_settings","processes_per_replica"):
			np = ui.get("parallel_settings","processes_per_replica")
			arglist += ["-n",np]
			mem = int(np)*ui.get_eval(sname,"srun_memory_per_core")
			arglist += ["--mem",str(mem)]

		arglist += ["--gres",ui.get(sname,"srun_gres_name")+":1"]
			
		if ui.has_option(sname,"additional_arguments"):
			arglist += ui.get_eval(sname,"additional_arguments")
		arglist += [self.bin]


	elif execute_command == "runjob":
		block_size=ui.get_eval('parallel_settings','nodes_per_replica')
		modes = ui.get_eval("parallel_settings","runjob_processes_per_node")
		arglist = ["runjob","--np",str(modes*block_size),"-p",str(modes),"--cwd",self.working_dir,"--exe",self.bin]
		if ui.has_option("parallel_settings","runjob_block"):
			arglist += ["--block",ui.get("parallel_settings","runjob_block")]
		if ui.has_option("parallel_settings","runjob_corner"):
			arglist += ["--corner",ui.get("parallel_settings","runjob_corner")]
		if ui.has_option("parallel_settings","runjob_shape"):
			arglist += ["--shape",ui.get("parallel_settings","runjob_shape")]
		if ui.has_option("parallel_settings","additional_arguments"):
			arglist += ui.get_eval("parallel_settings","additional_arguments")

	elif execute_command == "shell":
		os.chdir(self.working_dir)
		arglist = [self.bin]

	else:
		raise ValueError("Unknown execute command: %s; supporting mpirun, srun, runjob, and shell" % execute_command)

	aimsout=os.path.join(self.working_dir,"aims.out")
	aimserr=os.path.join(self.working_dir,"aims.err")

	if execute_command == "srun":
		#Requires special implementation
		stat = parallel_run.srun_call(arglist,aimsout,aimserr,self.replica)
		output.local_message("-- aims job exit status: "
					+ str(stat),self.replica)
		output.time_log("aims job exited with status "
				+ str(stat),self.replica)
		return stat

	if not enable_monitor:
		outfile = open(aimsout,"w")
		errfile = open(aimserr,"w")
		get_execute_clearance(request_folder=self.working_dir)
		output.time_log("aims job execute clearance acquired",self.replica)
		output.time_log("Aims execution with arguments: "+" ".join(map(str,arglist)),self.replica)
		p=subprocess.Popen(arglist,stdout=outfile,stderr=errfile)
		p.wait()
		end_of_execution_tasks()
		return True

	for i in range (10): #Allow 10 times for the job to successfully launch
		outfile=open(aimsout,"w")
		errfile = open(aimserr,"w")
		get_execute_clearance(request_folder=self.working_dir)
		output.time_log("aims job execute clearance acquired",self.replica)
		output.time_log("Aims execution with arguments: "+" ".join(map(str,arglist)),self.replica)
		p=subprocess.Popen(arglist,stdout=outfile,stderr=errfile)
		time.sleep(1)
		try:
			status=p.poll()
		except: #OSError Errno 3 Process does not exist
			output.time_log("Nodes failure ; replica will pass out from now on")
			time.sleep(86400)
			os.chdir(original_dir)
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
			os.chdir(original_dir)
			return False
	counter=0; last=os.stat(aimsout).st_size
	while counter<update_poll_times and p.poll()==None: #The output file needs to update at least once in every 5 minutes
		write_active(self.working_dir)
		self.set_permission()
		time.sleep(update_poll_interval)
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
	
	end_of_execution_tasks()
	

    def output(self, message): output.local_message(message, self.replica)

    def is_successful_absolute(self):
        '''
        checks if relaxation/optimization was successful
        '''
	aims_path = os.path.join(self.working_dir, 'aims.out')
	try:
	        aims_out = open(aims_path,"r")
	except:
		return False #File missing
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


    def is_successful_relative(self):
        '''
        checks if relaxation/optimization was relatively successful
	'''
        aims_path = os.path.join(self.working_dir, 'aims.out')
	try:
	        aims_out = open(aims_path,"r")
	except:
		return False
        counter = 0
        while True:
                line = aims_out.readline()
                if "Leaving FHI-aims" in line:
			return True
                elif line == '':
                        counter += 1
		else:
			counter = 0
                if counter > 10: #Allowing 10 empty lines in a row before determining eof
                        break
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
#	self.output("extracted energy" +str(energy))
        if energy == False: return "failed"  # did not converge
        self.result_struct.set_property('energy', energy)
       
	update_geometry = self.extract_geometry()
        if update_geometry == False: return "SPE" #No geometry changes
	
	latA = np.asarray(self.result_struct.properties["lattice_vector_a"])
	latB = np.asarray(self.result_struct.properties["lattice_vector_b"])
	latC = np.asarray(self.result_struct.properties["lattice_vector_c"])

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
        self.result_struct.set_property('b',b)
        self.result_struct.set_property('c',c)
	
	# extracts time from aims output
	wall_time = self.extract_time()
	self.result_struct.set_property('relax_time', wall_time)
        # extract the spin moment
        if self.ui.get('FHI-aims','initial_moment') == 'custom': self.extract_spin_moment()
        
        return "relaxed"

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
	counter = 0 
        while counter<10: #Allow 10 empty lines before determining EOF
            line = aims_out.readline()
            if not line: return False  # not converged
            if 'Final atomic structure' in line: break
            elif line == '':
                counter += 1
            else:
                counter = 0
	if counter==10:
            return False

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
	counter = 0 
        while counter<10: #Allow 10 empty lines before determining EOF
            line = aims_out.readline()
            if not line: return False  # not converged
            if 'Final atomic structure' in line: break
            elif line == '':
                counter += 1
            else:
                counter = 0
	if counter==10:
            return False
   
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


