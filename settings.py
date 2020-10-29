from misc import getValues
from sympyplus import UserDefinedFunction
import yaml

# Change the filename to choose a different settings file

####################################
settings_filename = 'settings.yml'
####################################

# Open the file and load it using YAML

settings_file = open(f'settings/{settings_filename}')
settings_dict = yaml.load(settings_file, Loader=yaml.FullLoader)

# Settings

class const:
	L0, L1, L2, L3, A, B, C, ep, W0, W1, W2 = getValues(settings_dict['const'],'L0, L1, L2, L3, A, B, C, ep, W0, W1, W2')
class meshdata: # Perhaps I could make a class that constains both the mesh and the meshdata later on
	file_path, numnodes_init, numnodes_max = getValues(settings_dict['meshdata'],'file_path, numnodes_init, numnodes_max')
class options:
	omit_init_printoff, visualize, manufactured = getValues(settings_dict['options'],'omit_init_printoff, visualize, manufactured')
class paraview: # I could also make a class for paraview stuff
	file_path = settings_dict['paraview']['file_path']
class solverdata:
	ksp_type, pc_type = getValues(settings_dict['solverdata'],'ksp_type, pc_type')
class timedata:
	time_step, end_time = getValues(settings_dict['timedata'],'time_step, end_time')

# Add 'dt' entry to 'const' dictionary

const.dt = timedata.time_step

# END OF CODE