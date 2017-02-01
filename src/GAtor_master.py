'''
New master script for running of GAtor

Created on June 26th, 2016
'''

import os, subprocess, shutil
import sys,socket
import core.file_handler as fh
from core import user_input, kill, output
from utilities import parallel_run,stoic_model, misc
try:
	import imp.reload as reload #Python 3.0-3.3
except:
	try:
		import importlib.reload as reload #Python 3.4 above
	except:
		pass #Python 2

def main():
#	main_process = GAtor(sys.argv[-1])
	src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
	sys.path.append(src_dir)
	main_process = GAtor()

class GAtor():
	'''
	This is the master class of GAtorr
	Takes the path to a configuration file as an necessary input
	'''
	def __init__(self):
		sname = "GAtor_master"
		self.ui = user_input.get_config()
		if self.ui.has_option(sname,"working_directory"):
			os.chdir(self.ui.get(sname,"working_directory"))
			reload_modules()
			
		if self.ui.get_boolean(sname,"testing_mode"):
			print("Testing mode starting")
			self.testing_mode()
			return #Under test mode, no further comamnd is read

		if self.ui.get_boolean(sname,"bundled_ga"):
			parallel_run.launch_bundled()
			return

		if self.ui.get_boolean(sname,"kill_ga") and\
					self.ui.is_master_process():
			self.kill_GA()

		if self.ui.get_boolean(sname,"clean_folder") and\
					self.ui.is_master_process():
			self.clean_folder()

		if self.ui.get_boolean(sname,"fill_initial_pool") and\
					self.ui.is_master_process():
			self.fill_initial_pool()

		if self.ui.get_boolean(sname,"run_ga"):
			self.run_ga()		

	def testing_mode(self):
		from utilities import test_and_debug
		if not self.ui.has_option("test_and_debug","testing_procedure"):
			raise ValueError("testing mode is enabled, but testing_procedure parameter is missing from test_and_debug section")

		test_procedure = self.ui.get("test_and_debug","testing_procedure")
		if self.ui.is_master_process():
			output.time_log("Testing and debugging mode enabled")
			output.time_log("Testing procedure used: "+test_procedure)
		
		getattr(test_and_debug,test_procedure)()
		return

	def kill_GA(self):
		kill.set_kill()
		return

	def fill_initial_pool(self):
		IP_module = fh.my_import(self.ui.get("modules","initial_pool_module"),package="initial_pool")
		fh.mkdir_p(fh.out_tmp_dir)
		IP_module.main()
		

	def clean_folder(self):
		sname = "clean_folder"
		output.time_log("Cleaning folder is activied; pending user confirmation",sname)
		user_answer = user_input.keyboard_input("Are you sure you want to clean the working directory (y/n):",allowed_answers=["y","n"],time_out=30)
		if user_answer=="y":
			output.time_log("User confirmation received. Begins cleaning.",sname)
			clean()
			#output.time_log("Cleaning completed.",sname)
		else:
			output.time_log("User denied cleaning or no user response is received within 30 s.",sname)
			output.time_log("Moving on",sname)

	def run_ga(self):
		from core import run_GA
		sname = "parallel_settings"
		if self.ui.get(sname,"parallelization_method") != "serial":
		#Launch parallelism
			parallel_run.launch_parallel()
			return

		if self.ui.get_replica_name() == "master":
			#Need to assign new name
			self.ui.set(sname,"replica_name",misc.get_random_index())
			conf_path = os.path.join(fh.conf_tmp_dir,
						 self.ui.get_replica_name()+".conf")

			f = open(conf_path,"w")
			self.ui.write(f)
			f.close()
			sys.argv.append(conf_path)
			reload_modules()

		message = "GAtor instance reporting from " + \
				socket.gethostname()
		if self.ui.has_option(sname,"allocated_nodes"):
			nodes = self.ui.get_eval(sname,"allocated_nodes")
			message += "; controlling node(s): " + \
					", ".join(map(str,nodes))

		if self.ui.has_option(sname,"processes_per_replica"):
			ppr = self.ui.get_eval(sname, "processes_per_replica")
			message += "; controlling %i process(es)" % ppr

		if self.ui.has_option(sname,"runjob_block"):
			message += "; block: " + self.ui.get(sname,"runjob_block")

		if self.ui.has_option(sname,"runjob_corner"):
			message += "; corner: " + self.ui.get(sname,"runjob_corner")

		output.time_log(message)
	
		stoic = stoic_model.determine_stoic()
		ga = run_GA.RunGA(self.ui.get_replica_name(),stoic)
		ga.start()
		output.move_to_shared_output(self.ui.get_replica_name())
	


def reload_modules():
	'''
	These modules need to be reloaded after the configuration file is changed
	'''
	reload(fh)
	reload(user_input)
	reload(output)
	reload(stoic_model)
	reload(misc)
	reload(kill)

def clean():
    print 'resetting environment'
    sname = "clean"
    directory_to_clean = [fh.tmp_dir, fh.structure_dir, fh.conf_tmp_dir,
    fh.out_tmp_dir]
    directory_to_remove = [fh.success_dir, fh.scavenge_dir, fh.fail_dir]
    try:
        files_to_remove = [fh.output_file, fh.restart_relaxation_file, fh.restart_replica_file]
    except: pass
    output.time_log("Cleaning all the .out and .err files in: " + fh.cwd, sname)
    p = subprocess.Popen(['rm *.out'], cwd=fh.cwd, shell=True)
    p.wait()
    p = subprocess.Popen(['rm *.err'], cwd=fh.cwd, shell=True)
    p.wait()
    # tmp index is to keep track of replica number
    for directory in directory_to_clean:
        output.time_log("Cleaning directory: "+directory, sname)
        fh.mkdir_p_clean(directory)
    for directory in directory_to_remove:
        if os.path.exists(directory):
            output.time_log("Removing directory: "+directory, sname)
            shutil.rmtree(directory,ignore_errors=True)
    for rmfile in files_to_remove:
        if os.path.exists(rmfile):
            output.time_log("Removing file: "+rmfile, sname)
            os.remove(rmfile)
    p = subprocess.Popen(['rm *.log'], cwd=fh.cwd, shell=True)
    p.wait()
    return			

if __name__ == "__main__":
	main()
