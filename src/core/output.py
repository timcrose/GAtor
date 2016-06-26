'''
Created on Dec 19, 2013

@author: newhouse
'''
import os

from core.file_handler import cwd, output_file, log_file
import time
from external_libs.filelock import FileLock

def restart_message(message):
    out_file = os.path.join(cwd, 'restart_relaxations.dat')
    data_file = open(out_file, 'a')
    data_file.write(str(message) + '\n')
    data_file.close()
    os.system("chmod g=u "+out_file)


def local_message(message, replica):
    out_file = os.path.join(cwd, str(replica) + '.out')
    data_file = open(out_file, 'a')
    data_file.write(str(message) + '\n')
    data_file.close()

def time_log(message,replica,file=log_file):
	message=time.strftime("%Y-%m-%d %H:%M:%S")+' '+replica+" : "+message+"\n"
	with FileLock(file,cwd,3600):
		if not os.path.exists(os.path.join(cwd,file)):
			f=open(os.path.join(cwd,file),"w")
			f.close()
			os.system("chmod g=u "+os.path.join(cwd,file))
		f=open(os.path.join(cwd,file),"a")
		f.write(message)
		f.close()

def error(message, replica=None):
    if replica == None: r = ''
    else: r = str(replica) + ' '
    out_file = os.path.join(cwd, 'error.out')
    data_file = open(out_file, 'a')
    data_file.write(r + str(message) + '\n')
    data_file.close()

def reset_local(replica):
    out_file = os.path.join(cwd, str(replica) + '.out')
    data_file = open(out_file, 'w')
    data_file.write(str('') + '\n')
    data_file.close()
    
def move_to_shared_output(replica):
    local_out_file = os.path.join(cwd, str(replica) + '.out')
    if not os.path.exists(local_out_file): pass
    else: 
        d_file = open(local_out_file, 'r')
        contents_string = d_file.read()
        d_file.close()
        
        data_file = open(output_file, 'a')
        data_file.write('Replica: ' + str(replica) + str(contents_string) + '\n')
        data_file.close()
	os.system("chmod g=u "+output_file)
    reset_local(replica)
