'''
@authors: patrick/farren
'''
import os
import subprocess
import time
import numpy as np
import shutil
import json
from copy import deepcopy
import core.file_handler as fh
from core import user_input, output
import core.file_handler as fh
from core.activity import *
from utilities import misc, parallel_run
from structures.structure import Structure
from structures.structure_collection import StructureCollection, get_collection


def main(input_struct):
    ''' 
    Inputs:
    A Structure object

    Returns
    A Structure object and a status

    Status includes:
    "failed": structure failed current relaxation
    "rejected": structure rejected due to certain threshold
    "accepted": structure accepted and should be further processed
    '''
    # Get user-defined settings and replica name
    ui = user_input.get_config()
    replica = ui.get_replica_name()

    # Make a copy of the input Structure
    input_struct = deepcopy(input_struct)

    # Set paths for DFT bin, control files, command
    control_dir = ui.get("FHI-aims","control_in_directory")
    control_list = ui.get_list("FHI-aims","control_in_filelist")
    DFT_bin = ui.get('FHI-aims','path_to_aims_executable')
    execute_command = ui.get('FHI-aims','execute_command')
    num_control_files = len(control_list)

    # Get names to store energy's for diff control files
    if ui.has_option("FHI-aims","store_energy_names"):
        energy_names = ui.get_list("FHI-aims","store_energy_names")
        if len(energy_names) != num_control_files:
            message = "Number of store energy names does not match "+\
                      "number of control.in files specified"
            raise ValueError(messsage)
    else:
        energy_names = ["energy"] * num_control_files

    # Set working directory path for aims calculation
    working_dir = os.path.join(fh.tmp_dir, ui.get_replica_name())

    header ="\n|--------------------------" +\
             " FHI-aims evaluation" + \
             "--------------------------|"
    output.local_message(header)

    # Setup and run aims for each user-defined control file
    for i in range(num_control_files):
        current_control_name = control_list[i]
        current_control_path = os.path.join(os.path.join(fh.cwd, control_dir),
                                    current_control_name)
        current_control_string = fh.read_data(current_control_path)
        FHI = FHIAimsEvaluation(i,
                                input_struct, 
                                working_dir, 
                                num_control_files,
                                current_control_name,
                                current_control_string,
                                DFT_bin,
                                execute_command,
                                replica)
        # Launch and run aims job
        FHI.execute()

        # Check if aims calculation was successful
        if not FHI.check_if_successful():
            FHI.save_failed_calculation()
            return input_struct, "failed"

        # Extract energy and check thresholds
        energy = FHI.extract_energy()
        if not energy:
            output.time_log("Total energy not found in aims.out")
            FHI.save_failed_calculation()
            return input_struct, "failed"

        # Add if energy is false, possibly b/c SCF cycle didnt converge
        if not FHI.is_energy_acceptable(energy):
            return input_struct, "rejected"

        # Extract energy and geometry from aims output
        result_struct = FHI.extract_and_build_structure()
        result_struct.set_property("energy", energy)
        result_struct.set_property(energy_names[i], energy)    
    
        # Print to output end of current aims evaluation
        output.local_message("\n"+
                             "|------------------------ " +
                             "End of FHI-aims evaluation" + 
                             "------------------------|")
        input_struct = deepcopy(result_struct) 

    return input_struct, "accepted" 




class FHIAimsEvaluation():
    '''
    Creates an optimizer object using FHI-Aims to execute energy
    evaluations and relaxations.
    
    This is an example of using a class to define a state, in this case the
    state of the relaxation process through setup, execution, and extraction
    '''
    def __init__(self, 
                 i,
                 input_struct, 
                 working_dir, 
                 num_control_files,
                 current_control_name,
                 current_control_string,
                 DFT_bin,
                 execute_command,
                 replica):
        self.iter = i
        self.working_dir = working_dir
        fh.mkdir_p_clean(self.working_dir)
        self.ui = user_input.get_config()   
        self.bin = DFT_bin
        self.execute_command = execute_command
        self.replica = replica
        self.input_struct = input_struct
        self.control_in_string = current_control_string
        self.control_len = num_control_files
        self.control_name = current_control_name
        self.set_monitor_execution()
        self.set_energy_thresholds()
        self.output_beginning_messages()
        self.create_geometry_in()
        self.create_geometry_orig()
        self.create_control_in()


    def output(self, message): output.local_message(message, self.replica)

    def output_time_log(self, message): output.time_log(message, self.replica)

    def output_beginning_messages(self):
        self.output_time_log("Executing FHI-aims using control file: %s"
                             % (self.control_name))
        self.output("\nFHI-aims evaluation using control file: %s"
                             % (self.control_name))
        self.output("-- Rejection energy threshold: %s eV"
                             % (self.energy_threshold))

    def create_geometry_in(self):
        '''
        Writes the geometry.in file called by FHI-aims. 
        Calls from FHI-aims control settings found in res directory.
        Returns: None
        '''
        geometry_file = open(os.path.join(self.working_dir, 'geometry.in'), 'w')
        geometry_file.write(self.input_struct.get_geometry_atom_format())
        geometry_file.close()
        json_file = open(os.path.join(self.working_dir, "struct.json"),"w")
        json_file.write(self.input_struct.dumps())
        json_file.close()

    def create_geometry_orig(self):
        '''
        Saves the original/unrelaxed aims geometry as geometry.orig
        '''
        orig_path = os.path.join(self.working_dir, 'geometry.orig')
        if not os.path.isfile(orig_path):
            geometry_file = open(orig_path, 'w')
            geometry_file.write(self.input_struct.get_geometry_atom_format())
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
        Run the aims executable with user-defined options
        '''
        begin_time = time.time()
        self.original_dir = os.getcwd()

        # Set argument list for system-specific commands to run executable
        arglist = self.return_execution_arglist()

        # Set aims.out and aims.err paths
        self.aimsout = os.path.join(self.working_dir,"aims.out")
        self.aimserr = os.path.join(self.working_dir,"aims.err")

        # Check proper monitoring keywords set in ui.conf
        if self.enable_monitor and (self.update_poll_interval is None or
                                    self.update_poll_time is None):
            raise ValueError("FHI-aims job monitoring enabled, " +
                             "but no poll interval or times specified")

        #Srun command requires special implementation
        if self.execute_command == "srun":
            #Requires special implementation
            stat = parallel_run.srun_call(arglist, 
                                     self.aimsout, 
                                     self.aimserr, 
                                     self.replica)
            output.local_message("-- aims job exit status: "
                        + str(stat),self.replica)
            output.time_log("aims job exited with status "
                    + str(stat),self.replica)
            return stat

        # Monitor aims job with advanced polling
        if self.enable_monitor:

            # Get clearance to run in directory
            self.outfile = open(self.aimsout,"w")
            self.errfile = open(self.aimserr,"w")
            get_execute_clearance(request_folder=self.working_dir)
            self.output_time_log("Aims job execute clearance acquired")
            self.output_time_log("Aims execution with arguments: "+
                                 " ".join(map(str,arglist)))

            # Use subprocessing to launch aims executable
            p = subprocess.Popen(arglist,
                               stdout=self.outfile,
                               stderr=self.errfile,
                               cwd=self.working_dir)
            self.monitor_launch(p)
            self.monitor_aims_output(p)
            self.end_of_execution_outputs(p.poll())

        # If job not monitored, don't need advanced job polling
        else:
            outfile = open(aimsout,"w")
            errfile = open(aimserr,"w")
            get_execute_clearance(request_folder=self.working_dir)
            output.time_log("Aims job execute clearance acquired", self.replica)
            output.time_log("Aims execution with arguments: "+
                            " ".join(map(str,arglist)),self.replica)
            p=subprocess.Popen(arglist, 
                               stdout=outfile, 
                               stderr=errfile, 
                               cwd=self.working_dir)
            p.wait()
            self.end_of_execution_outputs(p.poll())
            return True
        
        end_time = time.time()
        self.output("-- Job execution time: %s seconds"  
                            % str(end_time-begin_time))
        self.output_time_log("Job execution time: %s seconds"
                            % str(end_time-begin_time))
        return

    def check_if_successful(self):
        '''
        Looks for "Leaving FHI-aims" in aims.out
        Indicates calculation finished sucessfully
        '''
        aims_path = os.path.join(self.working_dir, 'aims.out')
        try:
            aims_out = open(aims_path,"r")
        except:
            return False #File missing
        counter = 0
        while True:
                line = aims_out.readline()
                if "Leaving FHI-aims" in line:
                    output.time_log("Execution with control file " + 
                                     self.control_name  + " completed.",
                                     self.replica)
                    return True
                elif line == '':
                    counter += 1
                else:
                    counter = 0
                #Allows 10 empty lines in a row before determining eof
                if counter > 10: 
                        break
        return False

    def extract_and_build_structure(self):
        '''
        Reads the output files from relaxation and 
        specifies properties of the structure (e.g. energy, geometry)
        Returns: status= "relaxed", "SPE" or False
        '''
        self.result_struct = Structure()
        # Update geometry if relaxed, if not return old structure
        update_geometry = self.extract_geometry()
        if update_geometry == False:
            return deepcopy(self.input_struct)
        
        # Else create new structure for relaxed geometry
        output.local_message("-- Updated geometry retrieved")

        # Lattice vectors and magnitudes
        latA = self.result_struct.get_property("lattice_vector_a")
        latB = self.result_struct.get_property("lattice_vector_b")
        latC = self.result_struct.get_property("lattice_vector_c")
        a = np.linalg.norm(np.asarray(latA))
        b = np.linalg.norm(np.asarray(latB))
        c = np.linalg.norm(np.asarray(latC))

        # Unit cell volume and angles
        temp_vol = np.dot(np.cross(latA, latB), latC)
        alpha = self.angle(latB, latC)
        beta = self.angle(latA, latC)
        gamma = self.angle(latA, latB)

        # Time for calculation
        wall_time = self.extract_time()

        # Set unit cell parameters and time as properties
        self.result_struct.set_property('cell_vol', temp_vol)
        self.result_struct.set_property('alpha',alpha)
        self.result_struct.set_property('beta', beta)
        self.result_struct.set_property('gamma', gamma)
        self.result_struct.set_property('a',a)
        self.result_struct.set_property('b',b)
        self.result_struct.set_property('c',c)
        self.result_struct.set_property('relax_time', wall_time)
        return self.result_struct


    def extract_energy(self):
        '''
        reads the resulting energy from the FHI-aims log
        Returns: float if successful, False if unsuccessful
        '''
        aims_out = open(os.path.join(self.working_dir, 'aims.out'))
        st = '| Total energy of the DFT / Hartree-Fock s.c.f. calculation      :'
        while True:
            line = aims_out.readline()
            if not line:
                return False  # energy not converged
            if st in line:
                tokens = line.split() 
                energy = float(tokens[11])  
                return energy

    def is_energy_acceptable(self, energy):
        ''' 
        Checks if energy of structure is within user-defined cutoffs 
        Returns: Boolean
        '''
        self.output_time_log("Checking if structure passes energy cutoff")
        self.output("-- Current energy:   %s  eV" % energy)
        self.output("-- Cutoff energy:    %s  eV" % self.energy_threshold) 
        if energy <= self.energy_threshold:
            self.output("-- Structure passes energy cutoff")
            self.output_time_log("Structure passes energy cutoff") 
            return True
        elif energy > self.energy_threshold:
            self.output("-- Structure does not pass energy cutoff")
            self.output_time_log("Structure does not pass energy cutoff")
            self.output_time_log("Performing new selection")
            return False

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

    def angle(self, v1, v2):
        numdot = np.dot(v1, v2)
        anglerad = np.arccos(numdot/(np.linalg.norm(v1)*np.linalg.norm(v2)))
        angledeg = anglerad*180/np.pi
        return angledeg

    def return_execution_arglist(self):
        '''
        Returns the system-specific commands for running the aims executable
        '''
        sname = "parallel_settings"
        ui = self.ui
        original_dir = os.getcwd()

        if self.execute_command == "mpirun":
            arglist = ["mpirun","-wdir",self.working_dir]
            if ui.has_option("parallel_settings","allocated_nodes"):
                arglist += ["-host",",".join(map(str,ui.get_eval("parallel_settings","allocated_nodes")))]
            if ui.has_option("parallel_settings","processes_per_replica"):
                arglist += ["-n",ui.get("parallel_settings","processes_per_replica")]
            if ui.has_option("parallel_settings","additional_arguments"):
                arglist += ui.get_eval("parallel_settings","additional_arguments")
            arglist += [self.bin]

        elif self.execute_command == "srun":
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

        elif self.execute_command == "aprun":
            arglist = ["aprun"]
            if ui.has_option("parallel_settings","processes_per_replica"):
                    arglist += ["-n",ui.get("parallel_settings","aims_processes_per_replica")]
            arglist+=["-e","OMP_NUM_THREADS=1"]
            if ui.has_option("parallel_settings","additional_arguments"):
                    arglist += ui.get_eval("parallel_settings","additional_arguments")
            arglist += [self.bin]

        elif self.execute_command == "runjob":
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

        elif self.execute_command == "shell":
            os.chdir(self.working_dir)
            arglist = [self.bin]

        else:
            raise ValueError("Unknown execute command: %s; supporting mpirun,"+
                                 " srun, runjob, and shell" % execute_command)
        return arglist

    def set_monitor_execution(self):
        sname = "FHI-aims"
        if self.ui.get_boolean(sname,"monitor_execution"):
            self.enable_monitor = True
            # Check monitor parameters match number of control files
            self.update_poll_times = self.ui.get_list(sname,"update_poll_times",eval=True)
            if len(self.update_poll_times)!=self.control_len:
                raise ValueError("Number of specified update poll times" + 
                                 "must match that of control.in files")
            self.update_poll_intervals = self.ui.get_list(sname,"update_poll_intervals",eval=True)
            if len(self.update_poll_intervals)!=self.control_len:
                raise ValueError("Number of specified update poll intervals" + 
                                 "must match that of control.in files")
            # Set monitor parameters for current control file
            self.update_poll_time = self.update_poll_times[self.iter]
            self.update_poll_interval = self.update_poll_intervals[self.iter]

        else:
            self.enable_monitor = False
            self.update_poll_time = None
            self.update_poll_interval = None

    def set_energy_thresholds(self):
        # Read Cutoff information from control file
        sname = "FHI-aims"
        if self.ui.has_option(sname,"store_energy_names"):
            self.sen = self.ui.get_list(sname,"store_energy_names")
            if len(self.sen)!=self.control_len:
                raise ValueError("Number of store energy names not matching number of control.in files specified")
        else:
            self.sen = ["energy"]*self.control_len
        if self.ui.has_option(sname,"absolute_energy_thresholds"):
            self.at=self.ui.get_list(sname,"absolute_energy_thresholds",eval=True)
            if len(self.at) != self.control_len:
                raise ValueError("Number of absolute energy thresholds not" + 
                              "matching number of control.in files specified")
        else:
            self.at = [None] * self.control_len
        if self.ui.has_option(sname,"relative_energy_thresholds"):
            self.rt=self.ui.get_list(sname,"relative_energy_thresholds",eval=True)
            if len(self.rt) != self.control_len:
                raise ValueError("Number of relative energy thresholds not" + 
                              "matching number of control.in files specified")
        else:
            self.rt = [None] * self.control_len
        if self.ui.has_option(sname,"reject_if_worst_energy"):
            self.worst_energy = self.ui.get_list_of_booleans(sname,"reject_if_worst_energy")
            if len(self.worst_energy) != self.control_len:
                raise ValueError("Number of reject_if_worst_energy boolean" + 
                    "flags not matching number of control.in files specified")
        else:
            self.worst_energy = [False]*self.control_len

        # Find minimum energy in collection according to keyword
        if self.rt[self.iter]!=None or self.worst_energy[self.iter]:
            struct_coll = get_collection(self.input_struct.get_stoic(), 0)
            struct_coll.update_local()
            energies = []
            energy_name = self.sen[self.iter]
            if energy_name is None:
                energy_name = "energy"
            for key, struct in struct_coll:
                energies.append(struct.get_property(energy_name))
            if energies[0] is None: 
                energies = []
                for key, struct in struct_coll:
                    energies.append(struct.get_property("energy"))

        # Relative Threshold
        self.energy_threshold = None
        if self.rt[self.iter] != None and self.rt[self.iter]!=None:
            if self.sen!=None:
                ekey = self.sen[self.iter]
            else:
                ekey = "energy"
            self.energy_threshold = min(energies) + self.rt[self.iter]

        # Absolute THreshold
        if self.at != None and self.at[self.iter]!=None: 
            if self.energy_threshold == None:
                self.energy_threshold = self.at[self.iter]
            else:
                self.energy_threshold = min(self.energy_threshold,
                                        self.energy_threshold[self.iter])
        # Worst energy threshold
        if self.worst_energy[self.iter]:
            struct_coll = get_collection(stoic,0)
            struct_coll.update_local()
            if self.energy_threshold == None:
                self.energy_threshold = max(energies)
            else:
                self.energy_threshold = min(self.energy_threshold, max(energies))

    def monitor_launch(self, p):
        '''
        Monitors the job launching process 
        Checks for aims output and/or failures
        '''
        # Allow 10 attempts for aims to start outputting
        for i in range (10):
            time.sleep(1)

            # Check subprocess exists
            try:
                status=p.poll()
            except: #OSError Errno 3 Process does not exist
                self.output_time_log("Nodes failure; replica will sleep from now on")
                time.sleep(86400)
                os.chdir(original_dir)
                return False
            
            #Allow 90 seconds for aims to start outputting
            time_limit = 90
            for j in range (time_limit): 
                if (p.poll()!=None) or (os.stat(self.aimsout).st_size>512):
                    break
                write_active(self.working_dir)
                time.sleep(1)
            if (os.stat(self.aimsout).st_size>512):
                self.output_time_log("aims.out begins output")
                break

            # If aims hasn't output anything, kill process and call failure
            self.output_time_log("aims job launch failure")
            try:
                p.send_signal(2)
            except:
                output.time_log("Unable to kill process; possible node failures")
                time.sleep(86400)
            active_sleep(60,self.working_dir)

            # If launch failure keeps happening, kill jobs
            if i==9:
                self.output_time_log("WARNING: Repeated launch failure; exiting")
                os.chdir(self.original_dir)
                return False

    def monitor_aims_output(self, p):
        ''' 
        Monitors the aims job as it runs
        Makes sure job doesn't get hung
        '''
        counter = 0
        last = os.stat(self.aimsout).st_size

        #The output file needs to update at least once in every update_poll_time
        while counter < self.update_poll_times and p.poll() == None: 
            write_active(self.working_dir)
            time.sleep(self.update_poll_interval)
            if os.stat(self.aimsout).st_size>last:
                last=os.stat(self.aimsout).st_size
                counter=0
            else:
                counter+=1
        #If hasn't updated in a while, consider job hung
        if counter == 60:
            output.time_log("aims job hung",self.replica)
            try:
                p.send_signal(2)
            except:
                output.time_log("Unable to kill process; possible node failures")
                time.sleep(86400)
            active_sleep(60, self.working_dir)


    def save_failed_calculation(self):
        '''
        Checks if aims calculation failed and optionally
        saves failed realxation folder
        '''
        self.output_time_log("FHI-aims execution using control file %s" + 
                            "failed to launch, hung, or failed" % (self.control_name))
        self.output("-- Job failed to launch, hung, or failed")
        if self.ui.get_boolean("FHI-aims","save_failed_calc"):
            path = os.path.abspath(os.path.join(fh.fail_dir,
                   misc.get_random_index()))
            shutil.copytree(self.working_dir,path)
            self.output_time_log("Failed calc folder saved to %s" % (path))
            self.output("-- Failed calc folder saved to %s" % (path))

    def end_of_execution_outputs(self, status):
        output.local_message("-- aims job exit status: "
                    + str(status),self.replica)
        output.time_log("aims job exited with status "
                + str(status),self.replica)
        os.chdir(self.original_dir)



