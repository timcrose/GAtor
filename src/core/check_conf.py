import os
import shutil
import user_input

class ConfigurationChecker():
    def __init__(self):
        self.ui = user_input.get_config()
        self.control_dir = self.ui.get("FHI-aims","control_in_directory")
        self.control_files = self.ui.get_list('FHI-aims','control_in_filelist')

    def run_checks(self):
        self.check_init_pool_paths()
        self.check_control_paths()


    def check_init_pool_paths():
        pass

    def check_control_paths(self):
        if not os.path.isdir(self.control_dir):
            msg = 'Control directory %s doesnt exist' % (self.control_dir)
            raise IOError(msg)
        for f in self.control_files:
            path = os.path.join(self.control_dir, f)
            if not os.path.isfile(path):
                msg = 'Control file %s doesnt exist' % (path)
                raise IOError(msg)

