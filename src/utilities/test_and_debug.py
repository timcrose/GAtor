'''
Module containing various testing and debugging procedures

Created by Patrick Kilecdi on June 25th, 2016
'''

from core import user_input,output
from utilities import parallel_run
import socket

ui = user_input.get_config()
def test_parallelism():
	'''
	Tests the parallelization method on different platforms
	'''
	sname = "parallel_settings"
	if ui.get_boolean(sname,"im_not_master_process"):
		message = "Replica %s reporting from %s" % (ui.get(sname,"replica_name"),socket.gethostname())
		if ui.has_option(sname,"allocated_nodes"):
			message += "; controlling nodes: "+", ".join(map(str,ui.get_eval(sname,"allocated_nodes")))
		if ui.has_option(sname,"processes_per_replica"):
			message += "; controlling %i processes" % (ui.get_eval(sname,"processes_per_replica"))
		output.time_log(message)
		return
	message = "test_parallelism launched to test parallelization method: "+ui.get(sname,"parallelization_method")
	output.time_log(message)
	parallel_run.launch_parallel()

