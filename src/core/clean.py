from file_handler import *
import subprocess

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
    clean()
