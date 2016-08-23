'''
Handles parallelization of the GA

Created by Patrick Kilecdi on June 26th, 2016
'''
import subprocess, os, math
import core.file_handler as fh
from core import user_input, output
from copy import deepcopy
from utilities import misc
from external_libs import bgqtools
import time
ui = user_input.get_config()

def launch_parallel():
	sname = "parallel_settings"
#	nor = ui.get_eval(sname,"number_of_replicas")
#	if nor == 1:
#		output.time_log("launch_replica is called by number_of_replicas is set to 1; exiting...",replica)
#		return None

	spawn_method = ui.get(sname,"parallelization_method")
	output.time_log("Parallelization method: "+spawn_method)
	if spawn_method == "serial":
		if ui.is_master_process():
			output.time_log("No parallelization is called")
			
		return 

	fh.mkdir_p(fh.conf_tmp_dir)
	
	if spawn_method == "subprocess":
		launch_parallel_subprocess()

	elif spawn_method == "mira" or spawn_method == "cetus":
		launch_parallel_bgq()
	
	#At this point, all is left is regular mpirun and srun spawn
	elif spawn_method == "mpirun":
		launch_parallel_mpirun()

	elif spawn_method == "srun":
		launch_parallel_mpirun(use_srun=True)
	
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
#	if npr > block_size:
#		message = "Number of nodes per replica (%i) " % npr
#		message += "larger than block size (%i); " % block_size
#		message += "Modify nodes_per_replica or bgq_block_size"
#		raise ValueError(message)

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


def launch_parallel_mpirun(use_srun=False):
	
	sname = "parallel_settings"
	output.time_log("mpirun/srun parallelization method is called")

#	nor = ui.get_eval(sname,"number_of_replicas")
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
		raise KeyError("mpirun/srun parallelization method requires the setting of at least one of the following parameter within parallel_settings: nodes_per_replica, processes_per_replica, and number_of_replicas; None is found")

	output.time_log("Number of parallel replicas: "+str(nor))
	output.time_log("Nodes assigned to each replica (0 indicates that this replica is assigned a fraction of a node and is not the last replica on the node: " + " ".join(map(str,npr)))
	output.time_log("Processes assigned to each replica: "+" ".join(map(str,ppr)))

	output.time_log("Total available nodes: %i; assigned nodes: %i" % (non,sum(npr)))

	if sum(npr) > non:
		output.time_log("Oversubscription of node has occured; Check the compatibility of number_of_replicas and nodes_per_replica/processes_per_node; aborting...")
		raise ValueError("Oversubscription of node has occured; Check the compatibility of number_of_replicas and nodes_per_replica/processes_per_node")
	elif sum(npr) < non:
		output.time_log("Not all nodes are utilized; try setting only number_of_replicas or setting a replica_per_node that divides the total number of nodes")

	output.time_log("Total available processes: %i; assigned processes: %i" % (nop,sum(ppr)))
	if sum(ppr)>nop:
		output.time_log("Oversubscription of processes has occured; This should be avoided in general, but GAtor will proceed for now.")
	elif sum(ppr)<nop:
		output.time_log("Undersubscription of processes has occured; If wish to utilize all processes, try setting only number_of_replicas or setting a processes_per_replica that fits the processes_per_node")

	processes = [] 
	new_ui = deepcopy(ui)
	new_ui.set(sname,"parallelization_method","serial")
	if not ui.has_option(sname,"processes_per_node"):
		new_ui.set(sname,"processes_per_node",str(ppn))
	
	python_command = ui.get(sname,"python_command")
	new_ui.set(sname,"im_not_master_process","TRUE")
	for i in range(nor):
		new_ui.set(sname,"replica_name",misc.get_random_index())
		new_ui.set(sname,"processes_per_replica",str(ppr[i]))
		new_ui.set(sname,"allocated_nodes",str(all_nodes[:max(1,npr[i])]))
		conf_path = os.path.join(fh.conf_tmp_dir,new_ui.get(sname,"replica_name")+".conf")
		f = open(conf_path,"w")
		new_ui.write(f)
		f.close()
		#Launch 1 instance of GAtor on the first node allocated
		if not use_srun:
			p = subprocess.Popen(["mpirun","-n","1","-host",all_nodes[0],python_command,fh.GAtor_master_path,conf_path]) 
		else:
			p = subprocess.Popen(["srun","-n","1","-w",all_nodes[0],python_command,fh.GAtor_master_path,conf_path])
		processes.append(p)
	
		all_nodes = all_nodes[npr[i]:]

	for p in processes:
		p.wait()


def get_all_processes(command,hostlist=None):
	'''
	This function returns all the processors available for the current job
	Specify hostlist if only wishes to access certain nodes
	'''
	if command!="mpirun" and command!="srun":
		raise ValueError("Unsupported command for get_all_hosts; only supporting mpirun and srun")
	arglist = [command]
	if hostlist!=None:
		if command == "mpirun":
			arglist += ["--host",",".join(map(str,hostlist))]
		elif command == "srun":
			arglist += ["-w",",".join(map(str,hostlist))]
	
	arglist += ["python",os.path.join(fh.src_dir,"utilities","print_host.py")]
	output.time_log("Acquiring all available processes using this command: "+" ".join(map(str,arglist)))

	p = subprocess.Popen(arglist,stdout=subprocess.PIPE)
	time.sleep(2)
	try:
		p.kill()
	except:
		pass
	out , err = p.communicate()
	try:
		out = str(out,"utf-8") #Python 3
	except:
		pass
	hosts = out.split("\n")
	hosts.pop() #Last line empty
	output.time_log("Number of processes acquired: "+str(len(hosts)))
	return hosts

def allocate_nodes(number_of_replicas,all_nodes,all_hosts,conf=ui):
	nor = number_of_replicas
	
