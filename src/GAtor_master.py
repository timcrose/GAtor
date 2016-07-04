'''
New master script for running of GAtor

Created on June 26th, 2016
'''

import os, subprocess, shutil
import sys,socket
import core.file_handler as fh
from core import user_input, kill, output
from utilities import parallel_run,stoic_model

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
		if self.ui.get_boolean(sname,"testing_mode"):
			print("Testing mode starting")
			self.testing_mode()
			return #Under test mode, no further comamnd is read

		if self.ui.get_boolean(sname,"kill_ga") and self.ui.is_master_process():
			self.kill_GA()

		if self.ui.get_boolean(sname,"clean_folder") and self.ui.is_master_process():
			self.clean_folder()

		if self.ui.get_boolean(sname,"fill_initial_pool") and self.ui.is_master_process():
			self.fill_initial_pool()

		if self.ui.get_boolean(sname,"run_ga"):
			from core import run_GA
			if self.ui.get("parallel_settings","parallelization_method")!="serial":
			#Launch parallelism
				parallel_run.launch_parallel()
			else:
				sname = "parallel_settings"
				message = "GAtor instance reporting from "+socket.gethostname()
				if self.ui.has_option(sname,"allocated_nodes"):
					message += "; controlling nodes: "+", ".join(map(str,self.ui.get_eval(sname,"allocated_nodes")))
				if self.ui.has_option(sname,"processes_per_replica"):
					message += "; controlling %i processes" % (self.ui.get_eval(sname,"processes_per_replica"))
				output.time_log(message)

				
				stoic = stoic_model.determine_stoic()
				ga = run_GA.RunGA(self.ui.get_replica_name(),stoic)
				ga.start()
				output.move_to_shared_output(self.ui.get_replica_name())
				

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
                fh.mkdir_p(fh.tmp_dir)
		fh.mkdir_p(fh.structure_dir)
		IP_module = fh.my_import(self.ui.get("modules","initial_pool_module"),package="initial_pool")
		#fh.mkdir_p(fh.tmp_dir)
		#fh.mkdir_p(fh.structure_dir)
		IP_module.main()
		

	def clean_folder(self):
		sname = "clean_folder"
		output.time_log("Cleaning folder is activied; pending user confirmation",sname)
		user_answer = user_input.keyboard_input("Are you sure you want to clean the working directory (y/n):",allowed_answers=["y","n"],time_out=30)
		if user_answer=="y":
			output.time_log("User confirmation received. Begins cleaning.",sname)
			clean()
			output.time_log("Cleaning compleated.",sname)
		else:
			output.time_log("User denied cleaning or no user response is received within 30 s.",sname)
			output.time_log("Moving on",sname)

			


def clean():
	print 'resetting environment'
	sname = "clean"
	directory_to_remove = [fh.tmp_dir, fh.structure_dir, fh.conf_tmp_dir]
	try:
		files_to_remove = [fh.output_file, fh.restart_relaxation_file, fh.restart_replica_file]
	except: pass
	output.time_log("Cleaning all the .out and .err files in: " + fh.cwd, sname)
	p = subprocess.Popen(['rm *.out'], cwd=fh.cwd, shell=True)
	p.wait()
	p = subprocess.Popen(['rm *.err'], cwd=fh.cwd, shell=True)
	p.wait()
	# tmp index is to keep track of replica number
	for directory in directory_to_remove:
		if os.path.exists(directory): 
			output.time_log("Removing directory: "+directory, sname)
			shutil.rmtree(directory)

	for rmfile in files_to_remove:
		if os.path.exists(rmfile): 
			output.time_log("Removing file: "+rmfile, sname)
			os.remove(rmfile)
	return			

if __name__ == "__main__":
	main()
