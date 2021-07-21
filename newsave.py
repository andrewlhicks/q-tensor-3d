import getopt
import os
import sys
import yaml

def nondimensionalize(const,R):
	""" For a constants dictionary const and a length R (for a sphere-like mesh,
	the radius), returns the non-dimensionalized version of these constants """

	from math import sqrt

	nd_const = {}

	# Define auxiliary constants

	AM = max(const['A'],const['B'],const['C'])
	LM = max(const['L1'],const['L2'],const['L3'])

	# Define non-dimensionalized constants

	nd_const['L0'] = const['L0']
	nd_const['A'] = const['A']/AM
	nd_const['B'] = const['B']/AM
	nd_const['C'] = const['C']/AM
	nd_const['L1'] = const['L1']/LM
	nd_const['L2'] = const['L2']/LM
	nd_const['L3'] = const['L3']/LM
	nd_const['W0'] = const['W0']*R/(LM)
	nd_const['W1'] = const['W1']*R/(LM)
	nd_const['W2'] = const['W2']*R/(LM)
	nd_const['q0'] = const['q0']*R
	nd_const['ep'] = sqrt(LM/AM)/R
	try:
		nd_const['beta'] = const['beta']
	except:
		nd_const['beta'] = 1

	# Return the non-dimensionalized constants

	return nd_const

def main():
	if not os.path.exists('saves'):
		os.makedirs('saves')

	help_text = 'python newsave.py <name> <settings-file>.yml <constants-file>.yml -n <radius>'
	radius = None

	# Process the required arguments

	if len(sys.argv[1:]) < 3:
		print(help_text)
		sys.exit()
	if not sys.argv[2].endswith('.yml') or not sys.argv[3].endswith('.yml'):
		print(help_text)
		sys.exit()

	# Now process the optional arguments

	try:
		opts, args = getopt.getopt(sys.argv[4:], "n:", ['nondim='])
	except getopt.GetoptError as err:
		# print help information and exit:
		print(err)  # will print something like "option -a not recognized"
		print(help_text)
		sys.exit()
	output = None
	verbose = False
	for o, a in opts:
		if o in ("-n", "--nondim"):
			radius = float(a)
		else:
			assert False, "unhandled option"

	# First, load the settings and the constants files

	settings_name = sys.argv[2]
	settings_path = 'settings/' + settings_name

	with open(settings_path) as settings_file:
		settings_dict = yaml.load(settings_file, Loader=yaml.Loader)

	constants_name = sys.argv[3]
	constants_path = 'constants/' + constants_name

	with open(constants_path) as constants_file:
		constants_dict = yaml.load(constants_file, Loader=yaml.Loader)

	if radius is not None:
		print(f"Creating non-dimensionalized version of '{sys.argv[3]}' with radius {radius}")
		constants_dict = nondimensionalize(constants_dict,radius)

	# Next, create the new save and put the settings and constants file therein

	save_name = sys.argv[1]
	save_path = 'saves/' + save_name

	if os.path.exists(save_path):
		i = 1
		while True:
			if not os.path.exists(f'{save_path}{i}'):
				new_save_name = f'{save_path}{i}'
				new_save_path = f'{save_path}{i}'
				break
			i += 1
		while True:
			answer = input(f"Save '{save_name}' already exists. Create save at '{new_save_name}' instead? (y/n) ")
			if answer in ('y','Y'):
				save_name = new_save_name
				save_path = new_save_path
				break
			if answer in ('n','N'):
				print("No new save created.")
				sys.exit()

	os.makedirs(save_path)
	os.makedirs(save_path+'/chk') # stores checkpoints
	os.makedirs(save_path+'/energy') # contains the plot of the energy decrease
	os.makedirs(save_path+'/vis') # contains paraview stuff

	with open(save_path + '/settings.yml','w') as settings_file:
		settings_file.write(yaml.dump(settings_dict))

	with open(save_path + '/constants.yml','w') as constants_file:
		constants_file.write(yaml.dump(constants_dict))

	print(f"New save '{save_name}' successfully created.")

if __name__ == '__main__':
	main()
