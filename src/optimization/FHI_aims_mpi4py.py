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
import time
import numpy as np
import shutil
import json
from copy import deepcopy
import core.file_handler as fh
from core import user_input, output
from core.activity import *
from utilities import misc
from structures.structure import Structure
from structures.structure_collection import StructureCollection, get_collection

from mpi4py import MPI
from external_libs import aims_w

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

def main(input_struct, replica, comm):
    ''' 
    Inputs:
    A Structure object, stoichiometry of crystal
    A Structure object and a status

    Status includes:
    "failed": structure failed current relaxation
    "rejected": structure rejected due to certain threshold
    "accepted": structure accepted and should be further processed
    '''
    # Get user-defined settings and replica name
    ui = user_input.get_config()

    # Set paths for DFT bin, control files, command
    control_dir = ui.get("FHI-aims","control_in_directory")
    control_list = ui.get_list("FHI-aims","control_in_filelist")
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
    working_dir = os.path.join(fh.tmp_dir, replica)

    header ="\n|--------------------------" +\
             " FHI-aims evaluation" + \
             "--------------------------|"
    if comm.Get_rank() == 0:
        output.local_message(header, replica)

    # Setup and run aims for each user-defined control file
    for i in range(num_control_files):
        # Get control file for current run
        current_control_name = control_list[i]
        current_control_path = os.path.join(os.path.join(fh.cwd, control_dir),
                                    current_control_name)
        current_control_string = fh.read_data(current_control_path)
 
        comm.Barrier()

        # Construct FHI instance
        FHI = FHIAimsEvaluation(i,
                                input_struct, 
                                working_dir, 
                                num_control_files,
                                current_control_name,
                                current_control_string,
                                replica,
                                comm)

        # Launch and run aims job
        FHI.execute()
        comm.Barrier()

        # Check if aims job was sucessful
        if not FHI.check_if_successful():
            FHI.save_failed_calculation()
            return input_struct, "failed"

        # Extract energy and check thresholds
        energy = FHI.extract_energy()
        if not energy:
            energy = FHI.extract_energy_no_scf()
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
        if comm.Get_rank() == 0:
            output.local_message("\n"+
                                 "|------------------------ " +
                                 "End of FHI-aims evaluation" + 
                                 "------------------------|", replica)
        input_struct = result_struct 
        comm.Barrier()
    
    comm.Barrier()
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
                 replica,
                 comm):
        self.iter = i
        self.comm = comm
        self.commf = comm.py2f()
        self.working_dir = working_dir
        self.ui = user_input.get_config()   
        self.replica = replica
        self.input_struct = input_struct
        self.control_in_string = current_control_string
        self.control_len = num_control_files
        self.control_name = current_control_name
        self.set_energy_thresholds()
        self.output_beginning_messages()
        self.create_geometry_in()
        self.create_geometry_orig()
        self.create_control_in()
        self.original_dir = fh.cwd


    def output(self, message): output.local_message(message, self.replica)

    def output_time_log(self, message): output.time_log(message, self.replica)

    def output_beginning_messages(self):
        if self.comm.Get_rank() == 0:
            self.output_time_log("Executing FHI-aims using control file: %s"
                                 % (self.control_name))
            self.output("\nFHI-aims evaluation using control file: %s"
                                 % (self.control_name))
            self.output("-- Rejection energy threshold: %s eV"
                                 % (self.energy_threshold))
        self.comm.Barrier() 

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
        Run the aims calculation using the aims wrapper
        aims_w
        '''
        from external_libs import aims_w 
        
        # Record Begging time
        begin_time = time.time()
        #self.original_dir = os.getcwd()

        # Set aims.out and aims.err 
        self.aimsout = os.path.join(self.working_dir,"aims.out")
        self.aimserr = os.path.join(self.working_dir,"aims.err")
        self.outfile = open(self.aimsout,"w")
        self.errfile = open(self.aimserr,"w")

        # Change into the directory where aims output is produced
        os.chdir(self.working_dir)

        # Write folder is active (used for restarts)
        write_active(self.working_dir)

        # Run aims job with fortran communicator commf
        self.comm.Barrier()
        aims_w.aims_w(self.commf)
        self.comm.Barrier()

        # End tasks
        end_time = time.time()
        if self.comm.Get_rank() == 0:
            self.output("-- Job execution time: %s seconds"  
                                % str(end_time-begin_time))
            self.output_time_log("Job execution time: %s seconds"
                                % str(end_time-begin_time))

    def check_if_successful(self):
        """
        Looks for "Leaving FHI-aims" in aims.out
        """
        aims_path = os.path.join(tmp_dir, self.replica, 'aims.out')
        try:
            aims_out = open(aims_path,"r")
        except:
            return False #File missing
        counter = 0
        while True:
            line = aims_out.readline()
            if "Leaving FHI-aims" in line:
                if self.comm.Get_rank() == 0:
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
        """
        Reads the output files from relaxation and 
        specifies properties of the structure (e.g. energy, geometry)
        Returns: status= "relaxed", "SPE" or False
        """
        self.result_struct = Structure()
        # Update geometry if relaxed, if not return old structure
        update_geometry = self.extract_geometry()
        if update_geometry == False:
            return deepcopy(self.input_struct)
        
        # Else create new structure for relaxed geometry
        output.local_message("-- Updated geometry retrieved")

        # Save any previous energies stored
        for name in self.sen:
            prev_energy = self.input_struct.get_property(name)
            if prev_energy is not None:
                self.result_struct.set_property(name, prev_energy)

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

    def extract_energy_no_scf(self):
        '''
        reads the resulting energy from the FHI-aims log
        looks for | Total energy corrected        : in cases
        where SCF cycle is not converged
        Returns: float if successful, False if unsuccessful
        '''
        st = '| Total energy corrected        :' 
        aims_out = open(os.path.join(self.working_dir, 'aims.out'))
        aims_data = reversed(aims_out.readlines())
        for line in aims_data:
            if st in line:
                tokens = line.split()
                energy = float(tokens[5])
                self.output(energy)
                return energy
        return False # energy not found

    def is_energy_acceptable(self, energy):
        ''' 
        Checks if energy of structure is within user-defined cutoffs 
        Returns: Boolean
        '''
        if self.comm.Get_rank() == 0:
            self.output_time_log("Checking if structure passes energy cutoff")
            self.output("-- Current energy:   %s  eV" % energy)
            self.output("-- Cutoff energy:    %s  eV" % self.energy_threshold) 
        if energy <= self.energy_threshold:
            if self.comm.Get_rank() == 0:
                self.output("-- Structure passes energy cutoff")
                self.output_time_log("Structure passes energy cutoff") 
            return True
        elif energy > self.energy_threshold:
            if self.comm.Get_rank() == 0:
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
        '''
        Computes unit cell angles from lattice vectors
        '''
        numdot = np.dot(v1, v2)
        anglerad = np.arccos(numdot/(np.linalg.norm(v1)*np.linalg.norm(v2)))
        angledeg = anglerad*180/np.pi
        return angledeg

    def set_energy_thresholds(self):
        '''
        Set energy cutoffs for different control files
        as specified by the user in ui.conf
        '''
        # Read info from ui.conf 
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
            #if energies[0] is None:
            if None in energies: 
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
            struct_coll = get_collection(self.input_struct.get_stoic(),0)
            struct_coll.update_local()
            if self.energy_threshold == None:
                self.energy_threshold = max(energies)
            else:
                self.energy_threshold = min(self.energy_threshold, max(energies))

    def save_failed_calculation(self):
        '''
        Checks if aims calculation failed and optionally
        saves failed realxation folder
        '''
        self.output_time_log("%s FHI-aims execution using control file" + 
                            "failed to launch, hung, or failed" %(self.comm.Get_rank()))
        self.output_time_log("-- Job failed to launch, hung, or failed")
        try:
            if self.ui.get_boolean("FHI-aims","save_failed_calc"):
                path = os.path.abspath(os.path.join(fh.fail_dir,
                       misc.get_random_index()))
        
                shutil.copytree(self.working_dir,path)
                self.output_time_log("Failed calc folder saved to %s" % (path))
                self.output("-- Failed calc folder saved to %s" % (path))
        except: pass

    def end_of_execution_outputs(self, status):                                
        output.local_message("-- aims job exit status: "                       
                    + str(status),self.replica)                                
        output.time_log("aims job exited with status "                         
                + str(status),self.replica)                                    
        os.chdir(self.original_dir) 

