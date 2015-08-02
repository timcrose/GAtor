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
	mkdir_p(fail_dir)
	ui=user_input.get_config()
	number_of_multi=ui.get_eval('parallel_settings','number_of_multiprocesses')
	restart_true = ui.get_eval('parallel_settings','restart_replicas')
	environment=ui.get('parallel_settings','system')
	print "Setting up multiprocessing for %i processes on the %s system" % (number_of_multi,environment)

	# Currently Not Used #
	#############################################################################
	mole_list=ui.get_eval('unit_cell_settings','molecule_list')
	for (molename,napm,occurance) in mole_list:
		if os.path.isfile(os.path.join(molecules_dir,molename+"_com_adjusted")):
			pass
		else:
			if not os.path.isfile(os.path.join(molecules_dir,molename)):
				raise RuntimeError("molecule geometry not found!")
				return
			atom_list=get_molecule_geo(molename,False)
			cm=[0,0,0]; tm=0
			for i in range (len(atom_list)):
				mass=ui.get_eval('molar_mass',atom_list[i][3])
				tm+=mass
				for j in range (3):
					cm[j]+=atom_list[i][j]*mass
			cm=[cm[j]/tm for j in range (3)]
			for atom in atom_list:
				for j in range (3):
					atom[j]-=cm[j]
#			print "this is atom_list",atom_list
			f=open(os.path.join(molecules_dir,molename+'_com_adjusted'),"w")
			for i in range (len(atom_list)):
				f.write('%f %f %f %s\n' % (atom_list[i][0],atom_list[i][1],atom_list[i][2],atom_list[i][3])) 
			#fix
			f.close()
	#######################################################################################

        if restart_true == True:
        	restart_files = [d for d in os.listdir(tmp_dir) if os.path.isdir(os.path.join(tmp_dir, d))]
                number_of_restarts = len(restart_files)
        else:
                restart_data_file = open(restart_replica_file, 'a')
                number_of_restarts = 0
	print "Number of restarts: "+str(number_of_restarts)
	if environment=='Cypress' or environment=='cypress' or environment=='1':
		#Replica name will be a random string
		pool=multiprocessing.Pool(processes=number_of_multi)
		replica_name_list=[get_random_index() for i in range (number_of_multi)]
		pool.map(run_GA_master,replica_name_list)
	elif environment=='Cypress_login' or environment=='cypress_login' or environment=='CYpress-login' or environment=='cypress-login':

		for i in range(number_of_multi):
			replica_name=get_random_index()
			print "In master.py, this is replica_name", replica_name
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
			os.remove('submit.ss')    
        #####CYPRESS SBATCH SUBMISSION SCRIPT#############
        elif environment=='Edison_login':
                for i in range (number_of_multi):
                        replica_name=get_random_index()
                        nodes = ui.get_eval('parallel_settings','nodes_per_replica')
                        time = ui.get('parallel_settings','replica_walltime')
                        ppn_edison = 24
                        width = ppn_edison*nodes
                        cwd = os.getcwd()
                        ex = os.path.join(src_dir,'core','run_GA.py')
                        exe_string="#!/bin/bash" + \
                                   "\n#PBS -q regular" + \
                                   "\n#PBS -A m2016" + \
                                   "\n#PBS -j eo" + \
                                   "\n#PBS -o " +str(replica_name)+".log" + \
                                   "\n#PBS -e " +str(replica_name)+".err" + \
                                   "\n#PBS -l mppwidth="+str(width) + \
                                   "\n#PBS -l walltime="+str(time) + \
				   "\n\ncd " +str(cwd)+ \
                                   "\npython " + str(ex)+ " -f " + str(ui_name)+ " --rn " +str(replica_name)
                        ss=open('submit.ss','w')
                        ss.write(exe_string)
                        ss.close()
                        os.system("qsub submit.ss")
	elif environment=="Cetus" or environment=="cetus" or environment=="mira" or environment=="Mira":
		block_size=ui.get_eval('parallel_settings','nodes_per_replica')
		from external_libs import bgqtools
		if block_size>=128:
		#Will get individual blocks
			
			partsize, partition, job_id = bgqtools.get_cobalt_info()
			blocks = bgqtools.get_bootable_blocks(partition,partsize,block_size)
			bgqtools.boot_blocks(blocks)
			number_of_multi=int(partsize/block_size)
#			print "In master.py, this is number_of_multi", number_of_multi
			pool=multiprocessing.Pool(processes=number_of_multi)
			pool.map(run_GA_master,blocks)
#			run_GA_master(blocks[0])
			
		else:
		#Will need to get corners; will boot blocks of 128
			partsize, partition, job_id = bgqtools.get_cobalt_info()
			if environment=="Cetus" or environment=="cetus":
				blocks=bgqtools.get_bootable_blocks(partition,partsize,128)
			else:
				blocks=bgqtools.get_bootable_blocks(partition,partsize,512)
			bgqtools.boot_blocks(blocks)
			corners = bgqtools.block_corner_iter(blocks,block_size)
			replica_name=[]
			while True: #converts the iterable into a list of replica names
				try:
					corner=next(corners)
					replica_name.append(corner[0]+'%'+corner[1]+'%'+corner[2]+"%"+get_random_index())
				except:
					break

			# Rename folders in tmp file to be the names of new replica names
			number_of_multi=int(partsize/block_size)
			pool=multiprocessing.Pool(processes=number_of_multi)
			pool.map(run_GA_master,replica_name)

	#mkdir_p(tmp_dir)
        #mkdir_p(structure_dir)
        #mkdir_p(fail_dir)
        #set_unkill()
	else:
		raise RuntimeError("parallel_settings.system not recognized!")
        
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
    try:
	files_to_remove = [output_file, restart_relaxation_file,restart_replica_file]
    except: pass
    p = subprocess.Popen(['rm *.out'], cwd=cwd, shell=True)
    p = subprocess.Popen(['rm *.err'], cwd=cwd, shell=True)	
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
