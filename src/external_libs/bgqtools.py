import itertools
import logging
import os
import sys
import socket
import subprocess

_fqdn = socket.getfqdn()
if 'alcf.anl.gov' in _fqdn:
   if 'cetus' in _fqdn:
      sys.path.append('/soft/cobalt/cetus/bgq_hardware_mapper/')
   else:
      sys.path.append('/soft/cobalt/bgq_hardware_mapper/')
   import coord_functions

class Unbuffered:
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)

def set_unbuffered_stdout():
    sys.stdout = Unbuffered(sys.stdout)

def get_corners(min_coords, max_coords, extents):
    corners = []
    for A in range(min_coords[0], max_coords[0]+1, extents[0]):
        for B in range(min_coords[1], max_coords[1]+1, extents[1]):
            for C in range(min_coords[2], max_coords[2]+1, extents[2]):
                for D in range(min_coords[3], max_coords[3]+1, extents[3]):
                    for E in range(min_coords[4], max_coords[4]+1, extents[4]):
                        corners.append([A,B,C,D,E])
    return corners

def get_block_data(name):
    name_tok = name.split('-')
    return ([int(m, 16) for m in name_tok[1]],
            [int(m, 16) for m in name_tok[2]],
            int(name_tok[len(name_tok)-1]))

def get_bootable_blocks(partition, nodes, block_size=None):
    if block_size is None:
        block_size = min(max(min_bootable_block_from_hostname(), nodes), 512) # booted blocks limited to 512 nodes size
    if not is_power2(block_size):
        block_size = closest_lower_power2(block_size)
    job = subprocess.Popen(['get-bootable-blocks', '--size=%d' % block_size, partition],
                           stdout=subprocess.PIPE)
    stdout, stderr = job.communicate()
    blocks = stdout.split()
    return blocks

def boot_blocks(blocks):
    logging.info('Booting blocks...')
    jobs = []
    for block in blocks:
        logging.info('  block %s' % block)
        jobs.append(subprocess.Popen(['boot-block', '--block', block]))
    for job in jobs:
        job.wait()
    logging.info('done booting blocks.')


def nodes_to_extents(nodes):
    """
    Return a 5-tuple of the optimal extent for a given number of nodes
    """
    _nodes_to_extents = {
       1:   (1, 1, 1, 1, 1),
       4:   (1, 1, 2, 2, 1),
       8:   (2, 1, 2, 2, 1),
       16:  (2, 2, 2, 2, 1),
       32:  (2, 2, 2, 2, 2),
       64:  (2, 2, 4, 2, 2),
       128: (2, 2, 4, 4, 2),
       256: (4, 2, 4, 4, 2),
       512: (4, 4, 4, 4, 2),
       1024: (4, 4, 4, 8, 2),
       2048: (4, 4, 4, 16, 2),
       4096: (4, 4, 8, 16, 2),
       8192: (4, 4, 16, 16, 2),
       12288: (8, 4, 12, 16, 2),
       16384: (4, 8, 16, 16, 2),
       24576: (4, 12, 16, 16, 2),
       32768: (8, 8, 16, 16, 2),
       49152: (8, 12, 16, 16, 2)
       }
   
    try:
       return _nodes_to_extents[nodes]
    except KeyError:
       raise ValueError('No extent layout for %d nodes' % nodes)

def block_corner_iter(blocks, nodes_per_subblock):
    """
    Iterate over the corners in `blocks` of sizes `nodes_per_subblock`.

    Yields tuples (block, corner, shape) suitable for passsing to
    runjob. nodes_per_sublock can be a scalar int or a list of ints.
    If it's a scalar, all sub-block partitions will have the same
    shape. If it's a list, each job will reserve a sub-block of the
    maximum requested size but smaller jobs will get a reduced 
    shape argument (i.e. some nodes will be unoccupied).
    """
    # check if nodes_per_subblock is scalar or list
    try:
       max_nodes = max(nodes_per_subblock)
       nodes_iter = iter(nodes_per_subblock)
    except TypeError:
       # fallback to constant nodes per job (int)
       max_nodes = nodes_per_subblock 
       nodes_iter = itertools.repeat(nodes_per_subblock)
    extents = nodes_to_extents(max_nodes)
    for block in blocks:
        min_coords, max_coords, size = get_block_data(block)
        if size > 512:
            # Size of sub-block must be 512 nodes or less, so just yield the whole block
            yield (block, None, None)
            continue

        for c in get_corners(min_coords, max_coords, extents):
            corner = coord_functions.get_compute_node_from_global_coords(c)
            job_nodes = next(nodes_iter)
            job_extents = nodes_to_extents(job_nodes)
            shape = 'x'.join(str(x) for x in job_extents)
            yield (block, corner, shape)

def get_hostname_ip():
    hostname = socket.gethostname()

    if 'alcf.anl.gov' in hostname:
       if 'vesta' in hostname:
           network='.vesta.itd.alcf.anl.gov'
       elif 'mira' in hostname:
           network='.mira.i2b.alcf.anl.gov'
       elif 'cetus' in hostname:
           network='.cetus.i2b.alcf.anl.gov'
       else:
          raise ValueError('Unknown login node hostname %s' % hostname)

       hostname_data = hostname + network
       ip = socket.gethostbyname(hostname_data)
    else:
       output = os.popen('ip addr show ib0').read()
       lines = output.split('\n')
       hostname = lines[2].split()[1].split('/')[0]
       ip = socket.gethostbyname(hostname)

    return (hostname, ip)

def get_cobalt_info():
    partsize = int(os.environ['COBALT_PARTSIZE'])
    partition = os.environ['COBALT_PARTNAME']
    job_id = int(os.environ['COBALT_JOBID'])
    
    return (partsize, partition, job_id)

def get_ll_info():
    partsize = int(os.environ['LOADL_BG_SIZE'])
    partition = os.environ['LOADL_BG_BLOCK']
    job_name = os.environ['LOADL_STEP_ID']

    return (partsize, partition, job_name)

def min_bootable_block_from_hostname():
    hostname = socket.gethostname()
    if 'vesta' in hostname:
       minblock = 64
    elif 'mira' in hostname:
       minblock = 512
    elif 'cetus' in hostname:
       minblock = 128
    else:
        raise ValueError('Function min_bootable_block_from_hostname() has unknown login node hostname %s' % hostname)
    return minblock

def is_power2(num):
    'States if a number is a power of two'
    return num != 0 and ((num & (num - 1)) == 0)

def closest_lower_power2(num):
    'Returns the closest power of two that is less than num. Limited to results between 32 (Vesta minimum) and 512 for our use.'
    return max([num & 2**n for n in range(5,10)])
