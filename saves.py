import os
import settings

def _choose_directory_name(directory_protoname):
	""" Chooses a directory name and returns it. More specifically, takes an
	arbitrary name, and then places a digit behind it to create a unique
	directory name which hasn't yet been used.  """

	ii = 0
	while True:
		name = f'{directory_protoname}{ii}'
		path = f'saves/{name}'
		if not os.path.exists(path):
			return name
		ii += 1

def _set_current_directory():
	""" Sets the global variable 'current_directory' to an appropriate path.
	More specifically, if the save mode is 'new', creates a new directory, using
	the save name as the protoname. If the save mode is 'overwrite' or 'resume',
	uses the save name as the explicit directory name. """

	global current_directory

	if settings.saves.mode == 'new':
		path = f'saves/{_choose_directory_name(settings.saves.name)}'
		if not os.path.exists(path):
			current_directory = path
			_create_directory(path)
		else:
			raise OSError(f'Directory {path} already exists.')
	elif settings.saves.mode == 'overwrite' or settings.saves.mode == 'resume':
		path = f'saves/{settings.saves.name}'
		if os.path.exists(path):
			current_directory = path
		else:
			raise OSError(f'Directory {path} does not exist.')
	else:
		raise ValueError('Save mode must be \'new\', \'overwrite\', or \'resume\'.')

def _create_directory(path):
	""" Creates a directory in the path 'path'. """

	os.makedirs(path)
	os.makedirs(path+'/chk') # stores checkpoints
	os.makedirs(path+'/energy') # contains the plot of the energy decrease
	os.makedirs(path+'/log') # contains the log of what is printed in the terminal
	os.makedirs(path+'/vis') # contains paraview stuff

###

def load_checkpoint(vector_space):
	from firedrake import Function, DumbCheckpoint, FILE_READ

	q_dump = Function(vector_space,name='dump')

	with DumbCheckpoint(f'{current_directory}/chk/dump',mode=FILE_READ) as chk:
		chk.load(q_dump)

	return q_dump

def save_checkpoint(q_dump):
	from firedrake import Function, DumbCheckpoint, FILE_CREATE

	if not isinstance(q_dump,Function):
		raise TypeError('Must be a Firedrake Function.')

	q_dump.rename('dump')

	with DumbCheckpoint(f'{current_directory}/chk/dump',mode=FILE_CREATE) as chk:
	    chk.store(q_dump)

def load_energies():
	import yaml

	with open(f'{current_directory}/energy/energies.yml') as energies_file:
		energies = yaml.load(energies_file, Loader=yaml.FullLoader)

	return energies

def save_energies(energies):
	import yaml

	with open(f'{current_directory}/energy/energies.yml','w') as energies_file:
		energies_file.write(yaml.dump(energies))

def load_times():
	import yaml

	with open(f'{current_directory}/energy/times.yml') as times_file:
		times = yaml.load(times_file, Loader=yaml.FullLoader)

	return times

def save_times(times):
	import yaml

	with open(f'{current_directory}/energy/times.yml','w') as times_file:
		times_file.write(yaml.dump(times))
