'''
Handles parallelization of the GA

Created by Patrick Kilecdi on June 26th, 2016
'''
import subprocess, os, math
import core.file_handler as fh
from core import user_input, output
from core.activity import *
from copy import deepcopy
from utilities import misc
from external_libs import bgqtools
from external_libs.filelock import FileLock
import time
ui = user_input.get_config()

def launch_parallel():
    # Get user-defined parallization method
    sname = "parallel_settings"
    spawn_method = ui.get(sname,"parallelization_method")
    output.time_log("Parallelization method: "+spawn_method)
    if spawn_method == "serial":
        if ui.is_master_process():
            output.time_log("No parallelization is called")
        return 
    fh.mkdir_p(fh.conf_tmp_dir)

    # Different subroutines for different parallelization methods	
    if spawn_method == "subprocess":
        launch_parallel_subprocess()
    elif spawn_method == "mira" or spawn_method == "cetus":
        launch_parallel_bgq()
    elif spawn_method == "mpirun":
		launch_parallel_mpirun()
    elif spawn_method == "srun":
		launch_parallel_mpirun(use_srun=True)
    elif spawn_method == "aprun":
        launch_parallel_aprun()
    elif spawn_method == "theta":
        launch_parallel_theta()
    else:
        raise ValueError("Unknown parallelization method: "+spawn_method)


def launch_bundled():
	sname = "bundled_run_settings"
	spawn_method = ui.get(sname,"parallelization_method")
	output.time_log("Parallelization method: "+spawn_method)

	if spawn_method == "mira" or spawn_method == "cetus":
		launch_bundled_bgq()
	else:
		raise ValueError("Unsupported bundling method: "+spawn_method)

def launch_parallel_subprocess():
	'''
	Launches new replicas using subprocessing
	'''
	sname = "parallel_settings"
	nor = ui.get_eval(sname,"number_of_replicas")
	python_command = ui.get(sname,"python_command")

	output.time_log("Beginning to launch new replicas using subprocessing")
	new_ui = deepcopy(ui)
	new_ui.set(sname,"parallelization_method","serial")
	new_ui.set(sname,"im_not_master_process","TRUE")
	new_ui.set(sname,"number_of_replicas","1")
	processes = []

	for i in range(nor):
		new_ui.set(sname,"replica_name",misc.get_random_index())		
		conf_path = os.path.join(fh.conf_tmp_dir,new_ui.get(sname,"replica_name")+".conf")
		f = open(conf_path,"w")
		new_ui.write(f)
		f.close()
		arglist = [python_command,fh.GAtor_master_path,conf_path]
		output.time_log("Running subprocess: "+" ".join(map(str,arglist)))
		p = subprocess.Popen(arglist)
		processes.append(p)
	
	for p in processes:
		p.wait()

def launch_parallel_bgq():
	'''
	Launches new replicas using subprocessing
	with each replica assigned an appropriate block and corner
	on ALCF's Mira and Cetus machine
	'''
	sname = "parallel_settings"
	spawn_method = ui.get(sname,"parallelization_method")
	python_command = ui.get(sname,"python_command")
	npr = ui.get_eval(sname,"nodes_per_replica")
	if ui.has_option(sname,"bgq_block_size"):
		block_size = ui.get_eval(sname,"bgq_block_size")
		if spawn_method == "mira" and block_size < 512:
			message = "bgq_block_size specified as %i\n" % (block_size)
			message+= "On mira, a block is at least 512 nodes!"
			raise ValueError(message)
	else:
		if spawn_method == "mira":
			block_size = 512
		else: #Cetus
			block_size = 128


	partsize, partition, job_id = bgqtools.get_cobalt_info()
	if block_size > partsize:
		message = "block size larger than the total number of nodes; "
		message += "modify bgq_block_size or submission script"
		raise ValueError(message)
	if npr > block_size:
		message = "Number of nodes per replica (%i) " % npr
		message += "larger than block size (%i); " % block_size
		message += "Modify nodes_per_replica or bgq_block_size"
		raise ValueError(message)

	otl = output.time_log
	otl("Launching parallel replicas on machines: " + spawn_method)
	otl("Total number of nodes: " + str(partsize))
	otl("Blocks of %i nodes are created" % block_size)
	blocks = bgqtools.get_bootable_blocks(partition, partsize, block_size)
	bgqtools.boot_blocks(blocks)
	otl("All blocks are booted")

	otl("Corners of %i nodes are created for each replica" % (npr))
	corners = bgqtools.block_corner_iter(blocks, npr)

	new_ui = deepcopy(ui)
	new_ui.set(sname,"parallelization_method","serial")
	new_ui.set(sname,"im_not_master_process","TRUE")

	processes = []
	for corner in corners:
		new_ui.set(sname,"replica_name",misc.get_random_index())
		new_ui.set(sname,"runjob_block", corner[0])
		new_ui.set(sname,"runjob_corner", corner[1])
		new_ui.set(sname,"runjob_shape",corner[2])

		conf_path = os.path.join(fh.conf_tmp_dir,
					 new_ui.get(sname,"replica_name")+".conf")
		f = open(conf_path,"w")
		new_ui.write(f)
		f.close()
	
		arglist = [python_command,fh.GAtor_master_path,conf_path]
		otl("Running subprocess: "+" ".join(map(str,arglist)))
		p = subprocess.Popen(arglist)
		processes.append(p)
	
	for p in processes:
		p.wait()


def launch_bundled_bgq():
	'''
	Launch multiple runs of GAtor bundled together
	'''
	sname = "bundled_run_settings"
	spawn_method = ui.get(sname,"parallelization_method")
	python_command = ui.get(sname,"python_command")
#	npr = ui.get_eval(sname,"nodes_per_replica")
	if ui.has_option(sname,"bgq_block_size"):
		block_size = ui.get_eval(sname,"bgq_block_size")
		if spawn_method == "mira" and block_size < 512:
			message = "bgq_block_size specified as %i\n" % (block_size)
			message+= "On mira, a block is at least 512 nodes!"
			raise ValueError(message)
	else:
		if spawn_method == "mira":
			block_size = 512
		else: #Cetus
			block_size = 128

	partsize, partition, job_id = bgqtools.get_cobalt_info()
	if block_size > partsize:
		message = "block size larger than the total number of nodes; "
		message += "modify bgq_block_size or submission script"
		raise ValueError(message)
	otl = output.time_log
	otl("Launching parallel replicas on machines: " + spawn_method)
	otl("Total number of nodes: " + str(partsize))
	otl("Blocks of %i nodes are created" % block_size)
	blocks = bgqtools.get_bootable_blocks(partition, partsize, block_size)
	bgqtools.boot_blocks(blocks)
	otl("All blocks are booted")

	run_names = ui.get_list(sname,"run_names")
	run_info = ui.get_bundled_run_info()
	sum_of_blocks = sum([int(run["number_of_blocks"]) for run in run_info])

	if sum_of_blocks > len(blocks):
		message = "Runs subscribe to %i blocks; " % sum_of_blocks
		message += "Only %i are available" % len(blocks)
		otl(message)
		raise ValueError(message)

	elif sum_of_blocks < len(blocks):
		message = "Under-subscription of blocks has occured!"
		otl(message)

	blocks_used = 0
	processes = []
	sname = "parallel_settings"

	for run in run_info:
		new_ui = user_input.get_config_with_path(run["config_file_path"])
		if new_ui.get_boolean("GAtor_master","fill_initial_pool"):
			#Need to fill the initial pool for each one
			otl("Setting up initial pool for run " + run["run_name"])
			new_ui.remove_section("GAtor_master")
			new_ui.add_section("GAtor_master")
			new_ui.set_working_dir(run["working_directory"])
			new_ui.set_true("GAtor_master","fill_initial_pool")
			new_ui.set(sname,"replica_name",misc.get_random_index())

			conf_path = os.path.join(fh.conf_tmp_dir,
					 new_ui.get(sname,"replica_name")+".conf")
			f = open(conf_path,"w")
			new_ui.write(f)
			f.close()
	
			arglist = [python_command,fh.GAtor_master_path,conf_path]
			otl("Running subprocess: "+" ".join(map(str,arglist)))
			p = subprocess.Popen(arglist)
			processes.append(p)

	for p in processes:
		p.wait()

	processes = []	
	for run in run_info:
		if run["nodes_per_replica"] > block_size:
			message = "Run %s nodes_per_replica exceeds block_size"\
					% run["run_name"]
			otl(message)
			raise ValueError(message)

		alloted = blocks[blocks_used : blocks_used+run["number_of_blocks"]]
		blocks_used += run["number_of_blocks"]


		corners = bgqtools.block_corner_iter(alloted,
						     run["nodes_per_replica"])
			
		new_ui = user_input.get_config_with_path(run["config_file_path"])
		new_ui.set(sname,"parallelization_method","serial")
		new_ui.set(sname,"im_not_master_process","TRUE")
		new_ui.set(sname,"nodes_per_replica",str(run["nodes_per_replica"]))
		new_ui.set_working_dir(run["working_directory"])	
		new_ui.set_false("GAtor_master","bundled_ga")

		otl("Run %s is allocated %i blocks" % (run["run_name"],
						       run["number_of_blocks"]))
		otl("Running with %i nodes per replica" % run["nodes_per_replica"])
		otl("Working directory: " + run["working_directory"])
		otl("Using configuration file: " + run["config_file_path"])

		for corner in corners:
			new_ui.set(sname,"replica_name",misc.get_random_index())
			new_ui.set(sname,"runjob_block", corner[0])
			new_ui.set(sname,"runjob_corner", corner[1])
			new_ui.set(sname,"runjob_shape",corner[2])
	
			conf_path = os.path.join(fh.conf_tmp_dir,
					 new_ui.get(sname,"replica_name")+".conf")
			f = open(conf_path,"w")
			new_ui.write(f)
			f.close()
	
			arglist = [python_command,fh.GAtor_master_path,conf_path]
			otl("Running subprocess: "+" ".join(map(str,arglist)))
			p = subprocess.Popen(arglist)
			processes.append(p)
		
	for p in processes:
		p.wait()

def launch_parallel_aprun():
    sname = "parallel_settings"
    output.time_log("aprun parallelization method is called")

    all_processes = get_all_processes("aprun")
    all_nodes = list(set(all_processes))

    output.time_log("All nodes: "+", ".join(map(str,all_nodes)))
    nop = len(all_processes) #Number of processes
    non = len(all_nodes) #Number of nodes
    ppn = int(nop/non) #Processes per node
    npr = [] #Nodes per replica
    ppr = [] #Processes per replica
    nor = 0 #Number of replicas

    if ui.has_option(sname,"processes_per_replica"):
        ppr = ui.get_eval(sname,"processes_per_replica")
        if ppr >= ppn:
            npr = int(math.ceil(ppr/(ppn+0.0)))
            if ui.has_option(sname,"number_of_replicas"):
                nor = ui.get_eval(sname,"number_of_replicas")
            else:
                nor = non/npr
            npr = [npr]*nor
            ppr = [ppr]*nor
        else:
            rpn = ppn/ppr
            if ui.has_option(sname,"number_of_replicas"):
                nor = ui.get_eval(sname,"number_of_replicas")
            else:
                nor = non*rpn
            npr = [int(0+(x%3)==0) for x in range(1,nor+1)]
            ppr = [ppr]*nor
    else:
        message = "aprun parallelization method requires the setting "
        message += "setting the following parameter within parallel_settings: "
        message += "nodes_per_replica, processes_per_replica"
        raise KeyError(message)


    # Setup individual ui.conf for each replica
    processes = []
    new_ui = deepcopy(ui)
    new_ui.set(sname,"parallelization_method","serial")
    if not ui.has_option(sname,"processes_per_node"):
        new_ui.set(sname,"processes_per_node",str(ppn))
    python_command = ui.get(sname,"python_command")
    new_ui.set(sname,"im_not_master_process","TRUE")

    # Main loop to launch parallel replicas
    node_count = 0

    for i in range(nor):
        if node_count == non:
            node_count = 0
        new_ui.set(sname,"replica_name",misc.get_random_index())
        new_ui.set(sname,"processes_per_replica",str(ppr[i]))
        new_ui.set(sname,"allocated_nodes","['"+str(all_nodes[node_count])+"']")
        conf_path = os.path.join(fh.conf_tmp_dir,new_ui.get(sname,"replica_name")+".conf")
        exe_path = os.path.join(fh.conf_tmp_dir,new_ui.get(sname,"replica_name")+".exe")
        f = open(conf_path,"w")
        new_ui.write(f)
        f.close()
        arglist = [python_command,fh.GAtor_master_path,conf_path]
        output.time_log("Running: "+" ".join(map(str,arglist)))
        p = subprocess.Popen(arglist)
        processes.append(p)
        node_count +=1
    for p in processes:
        p.wait()

def launch_parallel_theta():
    sname = "parallel_settings"
    output.time_log("aprun parallelization method is called")

    all_processes = get_all_processes("aprun")
    all_nodes = list(set(all_processes))

    output.time_log("All nodes: "+", ".join(map(str,all_nodes)))
    nop = len(all_processes) #Number of processes
    non = len(all_nodes) #Number of nodes
    ppn = int(nop/non) #Processes per node
    npr = [] #Nodes per replica
    ppr = [] #Processes per replica
    nor = 0 #Number of replicas

    if ui.has_option(sname,"processes_per_replica"):
        ppr = ui.get_eval(sname,"processes_per_replica")
        if ppr >= ppn:
            npr = int(math.ceil(ppr/(ppn+0.0)))
            if ui.has_option(sname,"number_of_replicas"):
                nor = ui.get_eval(sname,"number_of_replicas")
            else:
                nor = non/npr
            npr = [npr]*nor
            ppr = [ppr]*nor
        else:
            rpn = ppn/ppr
            if ui.has_option(sname,"number_of_replicas"):
                nor = ui.get_eval(sname,"number_of_replicas")
            else:
                nor = non*rpn
            npr = [int(0+(x%3)==0) for x in range(1,nor+1)]
            ppr = [ppr]*nor
    else:
        message = "aprun parallelization method requires the setting "
        message += "setting the following parameter within parallel_settings: "
        message += "nodes_per_replica, processes_per_replica"
        raise KeyError(message)


    # Setup individual ui.conf for each replica
    processes = []
    new_ui = deepcopy(ui)
    new_ui.set(sname,"parallelization_method","serial")
    if not ui.has_option(sname,"processes_per_node"):
        new_ui.set(sname,"processes_per_node",str(ppn))
    python_command = ui.get(sname,"python_command")
    new_ui.set(sname,"im_not_master_process","TRUE")

    # Main loop to launch parallel replicas
    node_count = 0

    for i in range(nor):
        if node_count == non:
            node_count = 0
        new_ui.set(sname,"replica_name",misc.get_random_index())
        new_ui.set(sname,"processes_per_replica",str(ppr[i]))
        new_ui.set(sname,"allocated_nodes","['"+str(all_nodes[node_count])+"']")
        conf_path = os.path.join(fh.conf_tmp_dir,new_ui.get(sname,"replica_name")+".conf")
        exe_path = os.path.join(fh.conf_tmp_dir,new_ui.get(sname,"replica_name")+".exe")
        f = open(conf_path,"w")
        new_ui.write(f)
        f.close()
        arglist = [python_command,fh.GAtor_master_path,conf_path]
        output.time_log("Running: "+" ".join(map(str,arglist)))
        p = subprocess.Popen(arglist)
        processes.append(p)
        node_count +=1
    for p in processes:
        p.wait()
def launch_parallel_mpirun(use_srun=False):
    sname = "parallel_settings"
    output.time_log("mpirun/srun parallelization method is called")
    if ui.has_option(sname,"allocated_nodes"):
    #Unsual case where a wrapper script pre-determines the nodes allocated to a master replica
        all_nodes = ui.get_eval(sname,"allocated_nodes")
        if not use_srun:
            all_processes = get_all_processes("mpirun",all_nodes)
        else:
            all_processes = get_all_processes("srun",all_nodes)
    else: #Obtains all nodes and processes available for this job
        if not use_srun:
            all_processes = get_all_processes("mpirun")
        else:
            all_processes = get_all_processes("srun")
        all_nodes = list(set(all_processes))

    output.time_log("All nodes: "+", ".join(map(str,all_nodes)))
    nop = len(all_processes) #Number of processes
    non = len(all_nodes) #Number of nodes
    ppn = int(nop/non) #Processes per node
    npr = [] #Nodes per replica
    ppr = [] #Processes per replica
    nor = 0 #Number of replicas

    if ui.has_option(sname,"nodes_per_replica"):
		npr = ui.get_eval(sname,"nodes_per_replica")
		if ui.has_option(sname,"processes_per_replica"):
			ppr = ui.get_eval(sname,"processes_per_replica")
		else:
			ppr = npr*ppn
		if ui.has_option(sname,"number_of_replicas"):
			nor = ui.get_eval(sname,"number_of_replicas")
		else:
			nor = non/npr
		npr = [npr]*nor
		ppr = [ppr]*nor
    elif ui.has_option(sname,"processes_per_replica"):
        ppr = ui.get_eval(sname,"processes_per_replica")
        if ppr >= ppn:
            npr = int(math.ceil(ppr/(ppn+0.0)))
            if ui.has_option(sname,"number_of_replicas"):
                nor = ui.get_eval(sname,"number_of_replicas")
            else:
                nor = non/npr
            npr = [npr]*nor
            ppr = [ppr]*nor
        else:
			rpn = ppn/ppr
			if ui.has_option(sname,"number_of_replicas"):
				nor = ui.get_eval(sname,"number_of_replicas")
			else:
				nor = non*rpn
			npr = [int(0+(x%3)==0) for x in range(1,nor+1)]
			ppr = [ppr]*nor
    elif ui.has_option(sname,"number_of_replicas"):
        nor = ui.get_eval(sname,"number_of_replicas")
        if nor > non: #Multiple replicas have to share same node
            rpn = nor / non
            add = nor % non #Certain nodes may have to handle more replicas
            npr = ([0]*rpn+[1])*add + ([0]*(rpn-1)+[1])*(non-add) 
            #Only the last replica on the node gets a 1
			#Later, when assigning nodes, this is when the next node is accessed
			#First assign the number of processes for each replica on the extra-loaded nodes
            ppr = ([ppn/(rpn+1)+1]*(ppn%(rpn+1))\
				+[ppn/(rpn+1)]*(rpn+1-(ppn%(rpn+1))))*add
            ppr+= ([ppn/rpn+1]*(ppn%rpn)+[ppn/rpn]*(rpn-(ppn%rpn)))*(non-add)
        else:
			npr = non / nor
			add = non % nor
			npr = [npr+1]*add + [npr]*(nor-add)
			ppr = [ppn*x for x in npr]
    else:
        message = "mpirun/srun parallelization method requires the setting " 
        message += "of at least one of the following parameter within parallel_settings: "
        message += "nodes_per_replica, processes_per_replica, and number_of_replicas; None is found"
        raise KeyError(message)    

    # Time Log outputs for Nodes
    message = " Number of parallel replicas: %s \n" % (nor)   
    message += "Nodes assigned to each replica: %s \n" % (npr)
    message += "(0 indicates that a replica is assigned a fraction of a node) \n"               
    message += "Processes assigned to each replica: %s \n" %(ppr)
    message += "Total available nodes: %i; assigned nodes: %s" % (non,sum(npr))
    output.time_log(message)

    # Check for Oversubscription or underutilization
    if sum(npr) > non:
        message = "Oversubscription of node has occurred;"
        message +="Check the compatibility of number_of_replicas "
        message +="and nodes_per_replica/processes_per_node; aborting..."
        output.time_log(message)
        raise ValueError(message)
    elif sum(npr) < non:
        message = "Not all nodes are utilized; try setting only number_of_replicas "
        message += "or setting a replica_per_node that divides the total number of nodes"
        output.time_log(message)

    # Time log outputs for processes
    message = "Total available processes: %i; assigned processes: %i" % (nop,sum(ppr))
    output.time_log(message)
    if sum(ppr)>nop:
        message = "Oversubscription of processes has occurred;"
        message += "This should be avoided in general,"
        message += "but GAtor will proceed for now."
        output.time_log(message)
    elif sum(ppr)<nop:
        message = "Undersubscription of processes has occurred;"
        message += "If wish to utilize all processes, try setting only number_of_replicas "
        message += "or setting a processes_per_replica that fits the processes_per_node"
        output.time_log(message)

    # Setup individual ui.conf for each replica
    processes = [] 
    new_ui = deepcopy(ui)
    new_ui.set(sname,"parallelization_method","serial")
    if not ui.has_option(sname,"processes_per_node"):
        new_ui.set(sname,"processes_per_node",str(ppn))
    python_command = ui.get(sname,"python_command")
    new_ui.set(sname,"im_not_master_process","TRUE")

    # Main loop to launch parallel replicas
    for i in range(nor):
        new_ui.set(sname,"replica_name",misc.get_random_index())
        new_ui.set(sname,"processes_per_replica",str(ppr[i]))
        new_ui.set(sname,"allocated_nodes",str(all_nodes[:max(1,npr[i])]))
        conf_path = os.path.join(fh.conf_tmp_dir,new_ui.get(sname,"replica_name")+".conf")
        exe_path = os.path.join(fh.conf_tmp_dir,new_ui.get(sname,"replica_name")+".exe")
        f = open(conf_path,"w")
        new_ui.write(f)
        f.close()
		#Launch 1 instance of GAtor on the first node allocated
        if not use_srun:
            p = subprocess.Popen(["mpirun","-n","1","-host",all_nodes[0],python_command,fh.GAtor_master_path,conf_path])
            processes.append(p)
        else:
            arguments = ["srun","-n","1",
				     "-w",all_nodes[0],
				     "-c",new_ui.get(sname,"processes_per_node"),
				     "--mem",ui.get(sname,"srun_gator_memory"),
				     "--gres",
				     ui.get(sname,"srun_gres_name")+":1",
				     python_command,
				     fh.GAtor_master_path,
				     conf_path]
            p = subprocess.Popen(arguments)
            processes.append(p)
            output.time_log("Running command: " + 
					" ".join(map(str,arguments)))
        all_nodes = all_nodes[npr[i]:]
    monitor_srun(processes)


def get_all_processes(command,hostlist=None):
    '''
    This function returns all the processors available for the current job
    Specify hostlist if only wishes to access certain nodes
    '''
    if command!="mpirun" and command!="srun" and command!="aprun":
        raise ValueError("Unsupported command for get_all_hosts; only supporting mpirun and srun")
    arglist = [command]
    if hostlist!=None:
        if command == "mpirun":
            arglist += ["--host",",".join(map(str,hostlist))]
        elif command == "srun":
            arglist += ["-w",",".join(map(str,hostlist))]
    elif command == "aprun":
        sname = "parallel_settings"
        npr = ui.get_eval(sname,"number_of_replicas")
        ppr = ui.get_eval(sname,"processes_per_replica")
        total_processes = npr*ppr
        arglist += ["-n",str(total_processes)]
        print "HERE"

    arglist += ["python",os.path.join(fh.src_dir,"utilities","print_host.py")]
    print arglist
    output.time_log("Acquiring all available processes using this command: "+" ".join(map(str,arglist)))
    p = subprocess.Popen(arglist,stdout=subprocess.PIPE)
    time.sleep(2)
    p.wait()
    out , err = p.communicate()
    try:
        out = str(out,"utf-8") #Python 3
    except:
        pass
    if command == "aprun":
        hosts = out.split("\n")
        hosts = hosts[:-2]
    else:
        hosts = out.split("\n")
        hosts.pop() #Last line empty
    output.time_log("hosts \n" + str(hosts))

    output.time_log("Number of processes acquired: "+str(len(hosts)))


    return hosts

def monitor_srun(processes):
	'''
	Takes all the currently running replicas
	Monitors the srun command file that submits the srun job
	Until all the processes exits
	'''
	sname = "parallel_settings"
	cfile = ui.get(sname, "srun_command_file")
	fname = os.path.basename(cfile)
	fdir = os.path.dirname(cfile)
	sfile = ui.get(sname,"srun_submitted_file")
	ssname = os.path.basename(sfile)
	sdir = os.path.dirname(sfile)
	rfile = ui.get(sname,"srun_completed_file")
	rname = os.path.basename(rfile)
	rdir = os.path.dirname(rfile)
	calls = []
	runtime = ui.get_eval(sname,"srun_max_runtime")
	otl = output.time_log
	received_id = []
	while True:
		if os.path.exists(cfile):
			with FileLock(fname,fdir):
				f = open(cfile,"r")
				commands = f.read()
				f.close()
				os.remove(cfile)
			commands = commands.split("\n")
			for command in commands:
				if len(command)==0:
					continue

				k = command.split("  ")
				arglist = eval(k[0])
				stdout = k[1]
				stderr = k[2]
				replica = k[3]
				job_id = k[4]
				if job_id in received_id:
				#Job already submitted
					continue

				otl("Execution with arguments: " + 
				    " ".join(map(str,arglist)),
				    replica)

				with FileLock(ssname,sdir,timeout=60):
				#Record that this job is received
					f = open(sfile,"a")
					f.write(job_id+"\n")
					f.close()

				outfile = open(stdout,"w")
				errfile = open(stderr,"w")
				p = subprocess.Popen(arglist,
						     stdout=outfile,
						     stderr=errfile)
				calls.append([p,time.time(),job_id])
				received_id.append(job_id)

		completed = []
		for call in calls:
			p = call[0]
			if p.poll()!=None: #Job completed
				with FileLock(rname,rdir):
					f = open(rfile,"a")
					f.write(call[2] + " " + str(p.poll())+"\n")
					f.close()
				completed.append(calls.index(call))
			
			if (time.time()-call[1]) > runtime:
				try:
					p.send_signal(2)
				except:
					pass

		calls = [call for call in calls if calls.index(call) not in completed]
		

		still_running = False
		for process in processes:
			if process.poll()==None:
				still_running = True
				break
		if not still_running: #All replicas exited
			break
				
		time.sleep(5)
				

def srun_call (arglist,stdout,stderr,replica):
	'''
	Submit a srun call to be picked up by monitor_srun in the master process
	'''
	sname = "parallel_settings"
	cfile = ui.get(sname, "srun_command_file")
	fname = os.path.basename(cfile)
	fdir = os.path.dirname(cfile)
	sfile = ui.get(sname,"srun_submitted_file")
	ssname = os.path.basename(sfile)
	sdir = os.path.dirname(sfile)
	rfile = ui.get(sname,"srun_completed_file")
	rname = os.path.basename(rfile)
	rdir = os.path.dirname(rfile)

	job_id = misc.get_random_index()
	output.local_message("-- srun execution internal job id: " + job_id)
	while True:
		with FileLock(fname,fdir):
			f = open(cfile,"a")
			f.write("  ".join(map(str,[arglist,stdout,stderr,replica,job_id]))+"\n")
			f.close()
		time.sleep(10)
		with FileLock(ssname,sdir):
			f = open(sfile,"r")
			received_id = f.read()
			f.close()
		if job_id in received_id.split("\n"):
		#Command picked up by master process
			break

	
	while True:
		if os.path.exists(rfile):
			with FileLock(rname,rdir,timeout=600):
				f = open(rfile,"r")
				completed = f.read()
				f.close()

			completed = completed.split("\n")
			for i in range(len(completed)):
				if job_id in completed[i]:
					stat=eval(completed[i].split(" ")[1])
					return stat
		time.sleep(5)

	return stat

def launch_and_monitor(working_dir,arglist,stdout,stderr,enable_monitor=False,update_poll_interval=None,update_poll_times=None,execute_srun=False,replica=ui.get_replica_name()):
	'''
	Conduct a binary call with stdout and stderr for a replica
	Originally adapted from FHI-aims binary call in the FHI-aims module
	If enable_monitor is set to TRUE, update_poll_interval and update_poll_times cannot be None
	update_poll_interval is the interval in seconds between polls of the process status
	update_poll_times is the number of polls without output before determining process hung
	execute_srun is set to True when this function is called by a master process

	NOTE: THIS FUNCTION'S ADAPTATION IS NOT COMPLETE YET 
	'''
	if arglist[0] == "srun" and not execute_srun:
		#Requires special implementation
		#Call is not directly made by the replica
		#But rather the master process that first spawned all of them
		return srun_call(arglist,
				 stdout,
				 stderr,
				 enable_monitor,
				 update_poll_interval,
				 update_poll_time,
				 replica)


	otl = output.time_log
	if not enable_monitor:
		outfile = open(stdout,"w")
		errfile = open(stderr,"w")

		get_execute_clearance(request_folder=working_dir)
		otl("Job execution clearance acquired",self.replica)
		otl("Execution without monitoring using arguments: " + 
		    " ".join(map(str,arglist)),
		    self.replica)

		p=subprocess.Popen(arglist,stdout=outfile,stderr=errfile)
		p.wait()
		return True

	for i in range (10): #Allow 10 times for the job to successfully launch
		outfile = open(stdout,"w")
		errfile = open(stderr,"w")
		get_execute_clearance(request_folder=working_dir)
		otl("Job execution clearance acquired",self.replica)
		otl("Job (with monitoring) launch attempt " + str(i) +
		    " with arguments: "+" ".join(map(str,arglist)),
		    self.replica)


		#Incomplete adaptation!!!!!!!
		p=subprocess.Popen(arglist,stdout=outfile,stderr=errfile)
		time.sleep(1)
		try:
			status=p.poll()
		except: #OSError Errno 3 Process does not exist
			otl("Nodes failure ; replica will pass out from now on")
			time.sleep(86400)
			os.chdir(original_dir)
			return False
			
		time_limit=60
		for j in range (time_limit): #Allow 60 seconds for aims to start outputting
			if (p.poll()!=None) or (os.stat(aimsout).st_size>512):
				break
			write_active(self.working_dir)
			self.set_permission()
			time.sleep(1)
		if (os.stat(aimsout).st_size>512):
			output.time_log("aims.out begins output", self.replica)
			break
		outfile.close()
		output.time_log("aims job launch failure",self.replica)
		try:
			p.send_signal(2)
		except:
			output.time_log("Unable to kill process ; possible node failures", self.replica)
			time.sleep(86400)
		active_sleep(60,self.working_dir)
		try:
			self.set_permission()
		except:
			pass


		if i==9:
			output.time_log("WARNING: Repeated launch failure ; exiting",self.replica)
			os.chdir(original_dir)
			return False
	counter=0; last=os.stat(aimsout).st_size
	while counter<update_poll_times and p.poll()==None: #The output file needs to update at least once in every 5 minutes
		write_active(self.working_dir)
		self.set_permission()
		time.sleep(update_poll_interval)
		if os.stat(aimsout).st_size>last:
			last=os.stat(aimsout).st_size
			counter=0
		else:
			counter+=1
	if counter==60:
		output.time_log("aims job hung",self.replica)
		try:
			p.send_signal(2)
		except:
			output.time_log("Unable to kill process ; possible node failures", self.replica)
			time.sleep(86400)
		active_sleep(60,self.working_dir)
		try:
			self.set_permission()
		except:
			pass
	 
						
