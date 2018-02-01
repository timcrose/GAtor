import os
import shutil
import user_input

class ConfigurationChecker():
    def __init__(self):
        self.ui = user_input.get_config()
        self.control_dir = self.ui.get("FHI-aims","control_in_directory")
        self.control_files = self.ui.get_list("FHI-aims","control_in_filelist")
        self.aims_x = self.ui.get("FHI-aims","path_to_aims_executable")
        self.initial_pool_dir = self.ui.get("initial_pool","user_structures_dir")
        self.selection_module = self.ui.get("modules","selection_module")

    def run_checks(self):
        self.check_init_pool_paths()
        self.check_control_paths()
        self.check_aims_executable()
        self.check_selection_parameters()

    def check_init_pool_paths(self):
        if os.path.isdir(self.initial_pool_dir):
            if os.listdir(self.initial_pool_dir) == []:
                msg = ('Initial pool directory %s is empty.' 
                                  % (self.initial_pool_dir))
                raise Exception(msg)
        if not os.path.isdir(self.initial_pool_dir):
            msg = ('Initial pool directory %s doesnt exist.' 
                                  % (self.initial_pool_dir))
            raise IOError(msg)

    def check_control_paths(self):
        if not os.path.isdir(self.control_dir):
            msg = 'Control directory %s doesnt exist.' % (self.control_dir)
            raise IOError(msg)
        for f in self.control_files:
            path = os.path.join(self.control_dir, f)
            if not os.path.isfile(path):
                msg = 'Control file %s doesnt exist.' % (path)
                raise IOError(msg)

    def check_aims_executable(self):
        if not os.path.isfile(self.aims_x):
            msg = "Aims executable %s does not exist" % (self.aims_x)
            raise IOError(msg)

    def check_selection_parameters(self):
        if self.selection_module == "tournament_selection":
            if not self.ui.has_option("selection","tournament_size"):
                msg = "Tournament selection called but [selection]/"
                msg +="tournament_size not set in .conf file"
                raise ValueError(msg)
