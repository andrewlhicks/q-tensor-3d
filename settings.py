from misc import getValues
import yaml
import os

# Functions that deal with the file

def _load_file(file_path):
	# Open the file and load it using YAML

	with open(file_path) as settings_file:
		settings_dict = yaml.load(settings_file, Loader=yaml.Loader)

	# Settings
	global mesh, options, visdata, saves, solver, time, constants, vis

	class mesh: # Perhaps I could make a class that constains both the mesh and the meshdata later on
		name, refs = getValues(settings_dict['mesh'],'name, refs')
	class options:
		manufactured, weak_boundary, strong_boundary = getValues(settings_dict['options'],'manufactured, weak_boundary, strong_boundary')
	class solver:
		grad_desc, ksp_type, ls_type, pc_type = getValues(settings_dict['solver'],'grad_desc, ksp_type, ls_type, pc_type')
	class time:
		save_every, step, num = getValues(settings_dict['time'],'save_every, step, num')
	class vis:
		normal = settings_dict['vis']['normal']

	if not isinstance(options.manufactured,bool):
	    raise TypeError("Variable 'manufactured' must be a boolean.")

# END OF CODE
