import errno
import errno
import os,sys
from shutil import rmtree
import time
from hashlib import sha1

from external_libs.filelock import FileLock

cwd = os.getcwd()

def argument_opt():
    from optparse import OptionParser
    parser=OptionParser()
    parser.add_option('-c','--clean',
                      '-r','--reset',
                      action='store_true',
                      dest='reset',
                      default=False,
                      help='Enables resetting environment before run')
    parser.add_option('-k',
                      action='store_true',
                      dest='kill',
                      default=False,
                      help='Terminates the GA')
    parser.add_option('-f','--file',
                      action='store',
                      type='string',
                      dest='user_input',
                      default='ui.conf',
                      help='User input file name (default="ui.conf"')
    parser.add_option('-d','--data',
                      action='store_true',
                      dest='data',
                      default=False,
                      help='Enables datatools')
    parser.add_option('-n','--dontrun',
                      action='store_false',
                      dest='run_e',
                      default=True,
                      help='Disables the actual running of the GA')
    parser.add_option('-i','--fill_ip',
                      action='store_true',
                      dest='fip_e',
                      default=False,
                      help='Enables reading initial pool from user defined directory')
    parser.add_option('-t','--test',
                      action="store_true",
                      dest="test_e",
                      default=False,
                      help="Enables testing and debugging mode and calls the testing\
                            procedure specified by test_and_debug.testing_procedure")
    parser.add_option('--rn','--replica_name',
                      action='store',
                      type='string',
                      dest='replica',
                      default=None, 
                      help='Replica name for the run_GA; not used for master.py')
    return parser.parse_args()


# source directories
src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
GAtor_master_path = os.path.join(src_dir,"GAtor_master.py")
GA_dir = os.path.abspath(os.path.join(src_dir, os.pardir))
res_dir = os.path.join(src_dir, 'res')

# working directories TODO: make this movable
tmp_dir = os.path.join(cwd, 'tmp/')
conf_tmp_dir = os.path.join(tmp_dir,'conf_tmp/')
out_tmp_dir = os.path.join(tmp_dir,'replica_out/')

#back-up directory
fail_dir = os.path.join(tmp_dir, 'save_calc_failed')  
success_dir = os.path.join(tmp_dir,"save_calc_success")
scavenge_dir = os.path.join(tmp_dir,"save_calc_scavenged")

# filesystem storage
structure_dir = os.path.join(cwd, 'structures')

# molecule directory
molecules_dir = os.path.join(cwd,"molecules")

# database storage
db_file = os.path.join(cwd, 'structures.sqlite') 

# useful files
progress_file = os.path.join(cwd, tmp_dir , 'progress.dat')
default_config = os.path.join(res_dir, 'default.conf')
(options,argv)=argument_opt() #Retrieve the user_input from command line
ui_conf = os.path.abspath(sys.argv[-1])
replica_file = os.path.join(tmp_dir, 'replica_index.dat')
output_file = os.path.join(cwd, 'GAtor.out')
error_file = os.path.join(cwd,'GAtor.err')
sys.stderr = open(error_file,"a")
log_file = os.path.join(cwd,"GAtor.log")
restart_relaxation_file = os.path.join(cwd, 'restart_relaxations.dat')
restart_replica_file = os.path.join(cwd, 'restart_replicas.dat')
# constants
INITIAL_POOL_REFID = -1

def mkdir_p(path):
    '''
    makes full directory path
    '''
    try:
	if not os.path.isdir(os.path.dirname(path)):
		mkdir_p(os.path.dirname(path))
        os.makedirs(path)
	try:
		os.system("chmod -R 771 "+path)
	except: pass
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
        
def mkdir_p_clean(path):
    '''
    removes directory and recreates it
    '''
    if os.path.exists(path): rmtree(path)
    mkdir_p(path)

def rmdir_silence(path):
	'''
	Removes directory silently
	'''
	try:
		rmtree(path)
	except:
		pass
    
def my_import(name, package=''):
    '''
    dynamically (at runtime) imports modules portentially specified in the UI
    taken from http://function.name/in/Python/__import__ 
    '''
    name = package + '.' + name
    mod = __import__(name)
    components = name.split('.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod 

def write_data(filepath, filename, contents):
    if not filename == None: filepath = os.path.join(filepath, filename)
    d_file = open(filepath, 'w')
    d_file.write(str(contents))
    d_file.close()
    os.system("chmod 771 "+filepath)

def read_data(filepath, filename=None):
    if filename is not None: full_filepath = os.path.join(filepath, filename)
    else: full_filepath = filepath
    d_file = open(full_filepath, 'r')
    contents_string = d_file.read()
    d_file.close()
    return contents_string
    
def get_progress():
    with FileLock(progress_file):
        p_file = open(progress_file, 'r')
        progress_string = p_file.readline()
        p_file.close()
    return progress_string

def set_progress(progress_string):
    with FileLock(progress_file):
        p_file = open(progress_file, 'w')
        p_file.write(progress_string)
        p_file.close()
    return True

def get_index(index_path):
    if not os.path.exists(index_path):
        data_file = open(index_path, 'w')
        data_file.write(str(0))
        data_file.close()
        
    with FileLock(index_path):
        index_file = open(index_path, 'r')
        index = int(index_file.readline())
        index_file.close()
    return index

def get_and_increment_index(index_path):
    '''
    retrieves the next valid index and increments the index in the shared file
    uses FileLock to avoid conflicts 
    '''
    # create if does not exist
    if not os.path.exists(index_path):
        mkdir_p(os.path.abspath(os.path.join(index_path, os.pardir)))
        data_file = open(index_path, 'w')
        data_file.write(str(0))
        data_file.close()
        
    with FileLock(index_path):
        index_file = open(index_path, 'r')
        index = int(index_file.readline())
        index_file.close()
    
        data_file = open(index_path, 'w')
        data_file.write(str(index + 1))
        data_file.close()
    return index

def get_random_index(seed=None):
    LENGTH_OF_INDEX = 10
    return sha1(repr(time.time())+str(seed)).hexdigest()[:LENGTH_OF_INDEX]

def get_molecule_geo(molename,adjusted=True):
	'''
	Reads in the molecule geometry from the molecules_dir
	if adjusted=True, then the file opened will be molename_com_adjusted
	'''
	if adjusted:
		f=open(os.path.join(molecules_dir,molename+"_com_adjusted"),"r")
	else:
		f=open(os.path.join(molecules_dir,molename),"r")
	st=f.readline()
	llist=[]
	while st!='':
		atom=st.split()
		if len(atom)<4:
			st=f.readline()
			continue
		for j in range (3):
			atom[j]=float(atom[j])
		llist.append(atom)
		st=f.readline()
	f.close()
	return llist

def print_to_file(message):
    with FileLock(output_file):
        data_file = open(output_file, 'a')
        data_file.write(str(message) + '\n')
        data_file.close()
    try:
	os.system("chmod g=u "+output_file)
    except:
	pass


if __name__ == '__main__':
    print cwd    
mkdir_p(tmp_dir)
mkdir_p(conf_tmp_dir)

