from misc import getValues
import yaml
import os

# Functions that deal with the file

def _load_file(file_path):
	# Open the file and load it using YAML

	with open(f'settings/{file_path}') as settings_file:
		settings_dict = yaml.load(settings_file, Loader=yaml.FullLoader)

	# Settings
	global mesh, options, visdata, saves, solverdata, timedata, constants

	class mesh: # Perhaps I could make a class that constains both the mesh and the meshdata later on
		name, refs = getValues(settings_dict['mesh'],'name, refs')
	class options:
		visualize, manufactured, weak_boundary, strong_boundary = getValues(settings_dict['options'],'visualize, manufactured, weak_boundary, strong_boundary')
	class saves:
		save, mode, name = getValues(settings_dict['saves'],'save, mode, name')
	class solverdata:
		ksp_type, pc_type = getValues(settings_dict['solverdata'],'ksp_type, pc_type')
	class timedata:
		time_step, end_time = getValues(settings_dict['timedata'],'time_step, end_time')
	class constants:
		file_path = settings_dict['constants']['file_path']

	if not isinstance(options.manufactured,bool):
	    raise ValueError("Variable 'manufactured' must be a boolean.")

	if not isinstance(options.visualize,bool):
	    raise ValueError("Variable 'visualize' must be a boolean.")

# END OF CODE
