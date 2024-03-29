"""                                                                            
If any part of this module is used for a publication please cite:              
                                                                               
F. Curtis, X. Li, T. Rose, A. Vazquez-Mayagoitia, S. Bhattacharya,             
L. M. Ghiringhelli, and N. Marom "GAtor: A First-Principles Genetic            
Algorithm for Molecular Crystal Structure Prediction",                         
J. Chem. Theory Comput., DOI: 10.1021/acs.jctc.7b01152;                        
arXiv 1802.08602 (2018)                                                        
"""   
import os
from shutil import rmtree, copytree
import time
from external_libs.filelock import FileLock
from core import user_input
from core.file_handler import *
from core import output

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

cwd = os.getcwd()

def get_execute_clearance(request_folder="",buffer=3,max_wait=1000):
        '''
        Reads the execute.info in the working directory and gets clearance 
        for executing commands such as runjob and mpirun
        '''
        for i in range (max_wait):
                with FileLock("execute.info",tmp_dir,3600):
                        if not os.path.exists(os.path.join(tmp_dir,"execute.info")):
                                data_file=open(os.path.join(tmp_dir,"execute.info"),"w")
                                data_file.write(str(time.time()))
                                data_file.close()
                                return True

                        data_file=open(os.path.join(tmp_dir,"execute.info"))
                        last_time=float(data_file.read())
                        data_file.close()
                        current_time=time.time()
                        if (current_time-last_time>buffer) or (i==max_wait):
                                data_file=open(os.path.join(tmp_dir,"execute.info"),"w")
                                data_file.write(str(time.time()))
                                data_file.close()
                                if i==max_wait and current_time-last_time<=buffer:
                                        return False
                                return True
                if request_folder!="":
                        write_active(request_folder)
                time.sleep(buffer)

def request_folder_to_check():
    '''
	Reads folder.info and returns the last one
	Updates folder.info if the time has exceeded user defined
    '''
    ui = user_input.get_config()
    time_wait = 60
    result = False
    with FileLock("folder.info",tmp_dir,3600):
        need_update = False
        if not os.path.exists(os.path.join(tmp_dir,"folder.info")):
            need_update = True
        else:
            f = open(os.path.join(tmp_dir,"folder.info"),"r")
            folder_list = f.read().split("\n")
            if folder_list[0] =="" or time.time()-float(folder_list[0])>time_wait:
                need_update=True
            f.close()
        if need_update:
            folder_list = [str(time.time())] #Time stamp comes first
            restart_wait = 30
            list_of_folders = ([name for name in os.listdir(tmp_dir) if 
                                os.path.isdir(os.path.join(tmp_dir,name)) and 
                                os.path.exists(os.path.join(tmp_dir,name,"active.info"))])
            for folder in list_of_folders:
                time_stamp = read_active(os.path.join(tmp_dir,folder))
                if time_stamp != False and time.time()-time_stamp>restart_wait:
                    folder_list.append(folder)
        f = open(os.path.join(tmp_dir,"folder.info"),"w")
        restart_wait = 30
        while len(folder_list)>1 and result==False: #Take the last on the list
            result=folder_list.pop()
            if (result=="") or not os.path.exists(os.path.join(tmp_dir,result,"active.info")):
                result = False
            else:
                time_stamp = read_active(os.path.join(tmp_dir,result))
                if time_stamp == False or time.time()-time_stamp < restart_wait: 
                    result = False
        for info in folder_list:
            f.write(info+"\n")
        f.close()
    return result

def bk_folder(fdir,folder,bk_path,naming_scheme="original",nname=get_random_index()):
        '''
        Backs up a folder in the bk_folder
        '''
        if naming_scheme=="original":
                nname=folder
        elif naming_scheme=="random":
                nname=get_random_index(seed=folder)
        elif naming_scheme=="const":
                nname=nname
        else:
            raise ValueError("Unknown naming_scheme in bk_folder")
        try:
            copytree(os.path.join(fdir,folder),os.path.join(bk_path,nname))
            output.time_log("bk_folder success: from %s to %s" % (os.path.join(fdir,folder),os.path.join(bk_path,nname)),"utility")
            return True
        except:
            output.time_log("bk_folder failure: from %s to %s" % (os.path.join(fdir,folder),os.path.join(bk_path,nname)),"utility!!!")
            return False

def write_active(fdir):
	'''
	Creates or overwrite a active.info file
	Recording current time stamp
	'''
	f=open(os.path.join(fdir,"active.info"),"w")
	f.write(str(time.time()))
	f.close()

def read_active(folder):
        '''
        Reads the active.info in the folder and returns the time
        '''
        try:
                f=open(os.path.join(folder,"active.info"))
                last_time=float(f.read())
                return last_time
        except:
                return False
def active_sleep(total_time,folder,write_active_interval=10):
        '''
        Sleep for the total time but continues to write_active in the folder
        '''
        times = int(total_time/write_active_interval)
        for i in range (times):
                write_active(folder)
                time.sleep(write_active_interval)
        write_active(folder)
        try:
                time.sleep(total_time-write_active_interval*times)
                write_active(folder)
        except:
                pass
					
if __name__ == '__main__':
    print(cwd)
