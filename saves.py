import os

def _choose_directory_name(directory_protoname):
	""" Chooses a directory name and returns it. More specifically, takes an
arbitrary name, and then places a digit behind it to create a unique directory
name which hasn't yet been used.  """
	ii = 0
	while True:
		path = f'saves/{directory_protoname}{ii}'
		if not os.path.exists(path):
			return path
		ii += 1

def _set_current_directory(directory_name):
	""" Sets the global variable 'current_directory' to the directory name you 
specify. """
	global current_directory
	current_directory = directory_name

def _create_directory(path):
	""" Creates a directory in the path 'path'. """
	os.makedirs(path)
	os.makedirs(path+'/log')
	os.makedirs(path+'/vis') # Directory that contains paraview stuff
	os.makedirs(path+'/energy')