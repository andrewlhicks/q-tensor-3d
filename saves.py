import os

def _create_directory(name):
	ii = 0
	while True:
		path = f'saves/{name}{ii}'
		if not os.path.exists(path):
			os.makedirs(path)
			os.makedirs(path+'/log')
			os.makedirs(path+'/vis') # Directory that contains paraview stuff
			os.makedirs(path+'/energy')
			global current_directory
			current_directory = path
			break
		ii += 1