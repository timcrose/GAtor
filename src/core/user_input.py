'''
@author: newhouse
This module contains the default user input values.
Values are overridden by textual user input
'''
from ConfigParser import SafeConfigParser
import ast

from core.file_handler import default_config, ui_conf


DEFAULT_CONFIG_REPLICA = -1

class ListSafeConfigParser(SafeConfigParser):
	'''Inherits SafeConfigParser and provides list parsing with json'''
	
	# TODO: maybe i could use literaleval(super.get()) instead, so to always return lists and ints
	def get_list(self, section, option):
		'''provides list parsing with json'''
		return self.get(section, option).split()
	
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

	def get_section_as_dict(self,section,eval=False):
		'''
		Return all the option under a section as a dictionary
		'''
		dicti = {}
		for key in self.options:
			if eval:
				dicti[key] = inst.get_eval(section,key)
			else:
				dicti[key] = inst.get(section,key)
		return dicti
			

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
