'''
@author: newhouse
This module contains the default user input values.
Values are overridden by textual user input
'''

from core.file_handler import default_config, ui_conf
import sys
from select import select
try:
        from ConfigParser import SafeConfigParser
except ImportError:
        from configparser import SafeConfigParser
import ast, os

try:
        from StringIO import StringIO
except ImportError:
        from io import StringIO


DEFAULT_CONFIG_REPLICA = -1

class ListSafeConfigParser(SafeConfigParser):
	'''Inherits SafeConfigParser and provides list parsing with json'''
	
	# TODO: maybe i could use literaleval(super.get()) instead, so to always return lists and ints
	def get_list(self, section, option,eval=False):
		'''provides list parsing with json'''
		if not eval:
			return self.get(section, option).split()
		else:
			return [ast.literal_eval(x) for x in self.get(section,option).split()]
	
	def get_eval(self, section, option):
		return ast.literal_eval(self.get(section, option))

	def get_boolean(self,section,option):
		'''
		Check and see if the section and option is specified
		if specified, has to set to "TRUE", else an error will be raised
		'''
		if self.has_option(section,option):
			if self.get(section,option)!="TRUE":
				raise ValueError("Optional boolean flag has to be set to TRUE if present")
			return True
		return False

	def get_list_of_booleans(self,section,option):
		'''
		Allows only TRUE and FALSE in the list
		'''
		l = self.get_list(section,option)
		result = []
		for k in l:
			if k == "TRUE":
				l.append(True)
			elif k == "FALSE":
				l.append(False)
			else:
				raise ValueError("List of boolean parameter must contain only TRUE and FALSE, all caps")
		
		

	def get_section_as_dict(self,section,eval=False):
		'''
		Return all the option under a section as a dictionary
		'''
		dicti = {}
		for key in self.options(section):
			if eval:
				dicti[key] = self.get_eval(section,key)
			else:
				dicti[key] = self.get(section,key)
		return dicti

	def get_replica_name(self):
		return self.get("parallel_settings","replica_name")
	
	def get_property_to_optimize(self):
		return self.get("run_settings","property_to_optimize")
	
	def verbose(self):
		return self.get_boolean("run_settings","verbose")

	def all_geo(self):
		return self.get_boolean("run_settings","output_all_geometries")
	
	def ortho(self):
		return self.get_boolean("unit_cell_settings","orthogonalize_unit_cell")

	def is_master_process(self):
		return not self.get_boolean("parallel_settings","im_not_master_process")
			

	def __deepcopy__(self,memo):
		'''
		Due to the inability to deepcopy a configuration file
		Will generate a temporary config file and read it back in
		'''
		config_string = StringIO()
		self.write(config_string)
		config_string.seek(0)
		copied = ListSafeConfigParser()
		copied.readfp(config_string)
		return copied





	
def get_config():
	'''
	Reads in default and user defined UI from the filesystem
	'''
	config = ListSafeConfigParser()

	default_config_file = open(default_config, 'r')
	config.readfp(default_config_file)
	default_config_file.close()

	local_config_file = open(ui_conf, 'r')
	config.readfp(local_config_file)
	local_config_file.close()
	return config

def keyboard_input(prompt,allowed_answers=None,time_out=86400, attempts=10):
	'''
	Allows interactive user input
	'''
	user_answer = None
	while attempts>0:
		sys.stdout.write(prompt+" ")
		sys.stdout.flush()
		rlist, _, _ = select([sys.stdin], [], [], time_out)
		if rlist:
			user_answer = sys.stdin.readline()
		else:
			return None
		if allowed_answers==None or user_answer[0:-1] in list(allowed_answers):
			return user_answer[0:-1]
		attempts -= 1
		
