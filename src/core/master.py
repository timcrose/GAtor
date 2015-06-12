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
from core.kill import set_kill


def argument_opt():
	from optparse import OptionParser
	parser=OptionParser()
	parser.add_option('-c','--clean','-r','--reset',action='store_true',dest='reset',default=False,help='Enables resetting environment before run')
	parser.add_option('-k',action='store_true',dest='kill',default=False,help='Terminates the GA')
	parser.add_option('-f','--file',action='store',type='string',dest='user_input',default='ui.conf',help='User input file name (default="ui.conf"')
	parser.add_option('-d','--data',action='store_true',dest='data',default=False,help='Enables datatools')
	parser.add_option('-n','--dontrun',action='store_false',dest='run_e',default=True,help='Disables the actual running of the GA')
	return parser.parse_args()

def main(reset_e,kill_e,data_e,run_e):
    # read command line arguments
    if reset_e:
        clean()
    if data_e:
	from core import data_tools
	data_tools.main(sys.argv)
	return
    if kill_e:
	set_kill()
	return
    if run_e:
	from core import user_input
	ui=user_input.get_config()
	number_of_multi=ui.get_eval('parallel_settings','number_of_multiprocesses')
	environment=ui.get_eval('parallel_settings','system')
	system_dict={1:"Cypress"}
	print "Setting up multiprocessing for %i processes on the %s system" % (number_of_multi,system_dict[environment])
	if environment==1: #Meaning cypress
		#Replica name will be a random string
		pool=multiprocessing.Pool(processes=number_of_multi)
		replica_name_list=[get_random_index() for i in range (number_of_multi)]
		pool.map(run_GA,replica_name_list)
    
        
def run_GA(replica):
	os.system('python '+os.path.join(src_dir,'core','run_GA.py')+' '+replica)
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
	try:
		(options,args)=argument_opt()
		user_input=options.user_input
		reset_e=options.reset
		kill_e=options.kill		
		data_e=options.data
		run_e=options.run_e
	except:
		pass
	main(reset_e,kill_e,data_e,run_e)
