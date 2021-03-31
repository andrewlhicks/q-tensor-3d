from misc import getValues
import yaml
import os

# Functions that deal with the file

def _load_file(file_path):
	# Open the file and load it using YAML

	with open(f'settings/{file_path}') as settings_file:
		settings_dict = yaml.load(settings_file, Loader=yaml.FullLoader)

	# Settings
	global mesh, options, visdata, saves, solver, time, constants, vis

	class mesh: # Perhaps I could make a class that constains both the mesh and the meshdata later on
		name, refs = getValues(settings_dict['mesh'],'name, refs')
	class options:
		visualize, manufactured, weak_boundary, strong_boundary = getValues(settings_dict['options'],'visualize, manufactured, weak_boundary, strong_boundary')
	class saves:
		save, mode, name = getValues(settings_dict['saves'],'save, mode, name')
	class solver:
		grad_desc, ksp_type, pc_type = getValues(settings_dict['solver'],'grad_desc, ksp_type, pc_type')
	class time:
		step, end = getValues(settings_dict['time'],'step, end')
	class constants:
		file_path = settings_dict['constants']['file_path']
	class vis:
		normal = settings_dict['vis']['normal']

	if not isinstance(options.manufactured,bool):
	    raise ValueError("Variable 'manufactured' must be a boolean.")

	if not isinstance(options.visualize,bool):
	    raise ValueError("Variable 'visualize' must be a boolean.")

# END OF CODE
