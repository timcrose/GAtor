#!/usr/bin/env python
'''
Created on Aug 5, 2013

@author: newhouse
'''

# adding src directory to path at runtime
import os
import shutil
import subprocess
import multiprocessing
import sys
import time
# add source directory to python path
src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.append(src_dir)

#from core import user_input
from core.file_handler import *
from core.kill import set_unkill,set_kill
from core import output

def main(ui_name,reset_e,kill_e,data_e,run_e,fip_e):
    # read command line arguments
    if reset_e: #User requires reseting environment
        clean()
    if fip_e:
	from core import IP_filling
	mkdir_p(tmp_dir) #tmp stores the added_user_structures.dat
	mkdir_p(structure_dir) #structure_dir stores the common storage
	IP_filling.main()
    if data_e: #Unknown feature
	from core import data_tools
	data_tools.main(sys.argv)
	return
    if kill_e: #User requires terminating GA
	set_kill()
	return
    if run_e:
	from core import user_input
	ui=user_input.get_config()
	number_of_multi=ui.get_eval('parallel_settings','number_of_multiprocesses')
	environment=ui.get('parallel_settings','system')
	print "Setting up multiprocessing for %i processes on the %s system" % (number_of_multi,environment)
	mkdir_p(tmp_dir)
	mkdir_p(structure_dir)
	set_unkill()
	if environment=='Cypress' or environment=='cypress' or environment=='1':
		#Replica name will be a random string
		pool=multiprocessing.Pool(processes=number_of_multi)
		replica_name_list=[get_random_index() for i in range (number_of_multi)]
		pool.map(run_GA_master,replica_name_list)
	elif environment=='Cypress_login' or environment=='cypress_login' or environment=='CYpress-login' or environment=='cypress-login':
		for i in range (number_of_multi):
			replica_name=get_random_index()
			exe_string='''#!/bin/bash
#SBATCH --qos=normal
#SBATCH --job-name=fhi_aims

#SBATCH --time=%s
#SBATCH -e %s.err
#SBATCH --nodes=%i
#SBATCH --ntasks-per-node=20
#SBATCH --workdir=%s

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE

python %s -f %s --rn %s
''' % (ui.get('parallel_settings','replica_walltime'),replica_name,ui.get_eval('parallel_settings','nodes_per_replica'),os.getcwd(),os.path.join(src_dir,'core','run_GA.py'),ui_name,replica_name)
			ss=open('submit.ss','w')
			ss.write(exe_string)
			ss.close()
			os.system("sbatch submit.ss")    
        
def run_GA_master(replica):
	from utilities.stoic_model import determine_stoic
	from core import run_GA	
	stoic=determine_stoic()
	ga=run_GA.RunGA(replica,stoic)
	ga.start()
	output.move_to_shared_output(replica)
#	os.system('python '+os.path.join(src_dir,'core','run_GA.py')+' '+replica)
#    minor_version = sys.version_info[1]
#    ui = user_input.get_config()
#        if (minor_version == 7) or (minor_version == 6): os.system('python ' + os.path.join(src_dir, 'core', 'run_GA.py ') + stoic_filename)  # laptop
#        if minor_version == 6: os.system('qsub jobs')  # thnec        
    
def clean():
    print 'resetting environment'
    directory_to_remove = [tmp_dir, structure_dir]
    files_to_remove = [output_file]
    p = subprocess.Popen(['rm *.out'], cwd=cwd, shell=True)
    p.wait()
    # tmp index is to keep track of replica number
    for directory in directory_to_remove:
        if os.path.exists(directory): rmtree(directory)
    for rmfile in files_to_remove:
        if os.path.exists(rmfile): os.remove(rmfile)
    return


if __name__ == '__main__':
	try:#PYthon 2.6 still uses optparse
		(options,args)=argument_opt()
		ui_name=options.user_input
		reset_e=options.reset
		kill_e=options.kill		
		data_e=options.data
		run_e=options.run_e
		fip_e=options.fip_e
	except:
		pass
	main(ui_name,reset_e,kill_e,data_e,run_e,fip_e)
