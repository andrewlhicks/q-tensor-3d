from misc import getValues, color
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
	from numpy import sqrt
	L1, L2, L3, q0, A, B, C, ep, W0, W1, W2 = getValues(settings_dict['const'],'L1, L2, L3, q0, A, B, C, ep, W0, W1, W2')
	L0 = 2*(A+B**2/C) if settings_dict['const']['L0'] == 'auto' else settings_dict['const']['L0']
	S0 = (B + sqrt(B**2 + 24.0*A*C))/(4.0*C)
class meshdata: # Perhaps I could make a class that constains both the mesh and the meshdata later on
	file_path, numnodes_init, numnodes_max = getValues(settings_dict['meshdata'],'file_path, numnodes_init, numnodes_max')
class options:
	omit_init_printoff, visualize, manufactured = getValues(settings_dict['options'],'omit_init_printoff, visualize, manufactured')
class visdata: # I could also make a class for paraview stuff
	file_path = settings_dict['visdata']['file_path']
class solverdata:
	ksp_type, pc_type = getValues(settings_dict['solverdata'],'ksp_type, pc_type')
class timedata:
	time_step, end_time = getValues(settings_dict['timedata'],'time_step, end_time')

# Add 'dt' entry to 'const' dictionary

const.dt = timedata.time_step

# Check if L1, L2, and L3 are set properly

if not isinstance(options.manufactured,bool):
    raise ValueError("Variable 'manufactured' must be a boolean.")

if not isinstance(options.omit_init_printoff,bool):
    raise ValueError("Variable 'omit_init_printoff' must be a boolean.")

if not isinstance(options.visualize,bool):
    raise ValueError("Variable 'visualize' must be a boolean.")

if 0>=const.L1 or -const.L1>=const.L3 or const.L3>=2*const.L1 or -3/5*const.L1-1/10*const.L3>=const.L2:
    print()
    print(f"{color.warning}WARNING: L1, L2, and L3 do not satisfy the proper inequalities{color.end}")

# END OF CODE