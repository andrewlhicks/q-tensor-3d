import getopt
import os
import sys
import yaml

usage_str = """usage: python newsave.py -b [-n <radius>] <savename> <settingsfile>.yml <constantsfile>.yml
                         -c <old-savename> <new-savename>
"""

def usage():
	print(usage_str)

def getopts(arglist,shorthands,longhands):
	try:
		opts, args = getopt.getopt(arglist,shorthands,longhands)
	except getopt.GetoptError as err:
		# print help information and exit:
		print(err)  # will print something like "option -a not recognized"
		sys.exit()
	return opts

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

def build():
    argv = sys.argv[2:] # truncate

    radius = None

    opts = getopts(argv, 'n:', [])
    for o, a in opts:
        if o in ('-n'):
            radius = float(a)
            argv = argv[2:] # after processing, remove (o,a) pair from argv
        else:
        	assert False, "unhandled option"

    if len(argv) != 3:
        print("3 args required for -b.")
        sys.exit()
    if not argv[1].endswith('.yml') or not argv[2].endswith('.yml'):
        print("Arg 2 and 3 must end in '.yml'.")
        sys.exit()

    # Process arguments

    save_name = argv[0]
    settings_name = argv[1]
    constants_name = argv[2]

    # Get paths

    settings_path = 'settings/' + settings_name
    constants_path = 'constants/' + constants_name
    save_path = 'saves/' + save_name

    # Load settings and constants dicts

    with open(settings_path) as settings_file:
    	settings_dict = yaml.load(settings_file, Loader=yaml.Loader)
    with open(constants_path) as constants_file:
    	constants_dict = yaml.load(constants_file, Loader=yaml.Loader)

    if radius is not None:
    	print(f"Creating non-dimensionalized version of '{constants_name}' with radius {radius}")
    	constants_dict = nondimensionalize(constants_dict,radius)

    create_save(save_name,save_path,settings_dict,constants_dict)

def copy():
    argv = sys.argv[2:]

    if len(argv) != 2:
        print("2 args required for -c.")
        sys.exit()

    # Process arguments

    old_save_name = argv[0]
    save_name = argv[1]

    # Get paths

    settings_path = 'saves/' + old_save_name + '/settings.yml'
    constants_path = 'saves/' + old_save_name + '/constants.yml'
    save_path = 'saves/' + save_name

    # Load settings and constants dicts

    with open(settings_path) as settings_file:
    	settings_dict = yaml.load(settings_file, Loader=yaml.Loader)
    with open(constants_path) as constants_file:
    	constants_dict = yaml.load(constants_file, Loader=yaml.Loader)

    create_save(save_name,save_path,settings_dict,constants_dict)

def create_save(save_name,save_path,settings_dict,constants_dict):
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

def main():
    if not os.path.exists('saves'):
    	os.makedirs('saves')

    try:
    	opts, args = getopt.getopt([sys.argv[1]],'bc',['help'])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
        print("Must specify -b or -c. Use 'python newsave.py --help' for usage.")
        sys.exit()

    for o, a in opts:
        if o in ('-b'):
            build()
        elif o in ('-c'):
            copy()
        elif o in ('--help'):
            usage()
        else:
            sys.exit()

    # # Process the required arguments
    #
    # if len(sys.argv[1:]) < 3:
    # 	print(help_text)
    # 	sys.exit()
    # if not sys.argv[2].endswith('.yml') or not sys.argv[3].endswith('.yml'):
    # 	print(help_text)
    # 	sys.exit()
    #
    # # Now process the optional arguments
    #
    # opts = getopts(sys.argv[1], "bc", [])
    #
    # for o, a in opts:
    # 	if o in ('-b'):
    # 		if len(sys.argv[2:]) != 3:
    # 			print("Three options required for '-b'.")
    # 			sys.exit()
    # 		if not sys.argv[2].endswith('.yml') or not sys.argv[3].endswith('.yml'):
    # 		radius = float(a)
    # 	if o in ('-c'):
    # 	else:
    # 		assert False, "unhandled option"
    #
    # # First, load the settings and the constants files
    #
    # settings_name = sys.argv[2]
    # settings_path = 'settings/' + settings_name
    #
    # with open(settings_path) as settings_file:
    # 	settings_dict = yaml.load(settings_file, Loader=yaml.Loader)
    #
    # constants_name = sys.argv[3]
    # constants_path = 'constants/' + constants_name
    #
    # with open(constants_path) as constants_file:
    # 	constants_dict = yaml.load(constants_file, Loader=yaml.Loader)
    #
    # if radius is not None:
    # 	print(f"Creating non-dimensionalized version of '{sys.argv[3]}' with radius {radius}")
    # 	constants_dict = nondimensionalize(constants_dict,radius)
    #
    # # Next, create the new save and put the settings and constants file therein
    #
    # save_name = sys.argv[1]
    # save_path = 'saves/' + save_name
    #
    # if os.path.exists(save_path):
    # 	i = 1
    # 	while True:
    # 		if not os.path.exists(f'{save_path}{i}'):
    # 			new_save_name = f'{save_path}{i}'
    # 			new_save_path = f'{save_path}{i}'
    # 			break
    # 		i += 1
    # 	while True:
    # 		answer = input(f"Save '{save_name}' already exists. Create save at '{new_save_name}' instead? (y/n) ")
    # 		if answer in ('y','Y'):
    # 			save_name = new_save_name
    # 			save_path = new_save_path
    # 			break
    # 		if answer in ('n','N'):
    # 			print("No new save created.")
    # 			sys.exit()
    #
    # os.makedirs(save_path)
    # os.makedirs(save_path+'/chk') # stores checkpoints
    # os.makedirs(save_path+'/energy') # contains the plot of the energy decrease
    # os.makedirs(save_path+'/vis') # contains paraview stuff
    #
    # with open(save_path + '/settings.yml','w') as settings_file:
    # 	settings_file.write(yaml.dump(settings_dict))
    #
    # with open(save_path + '/constants.yml','w') as constants_file:
    # 	constants_file.write(yaml.dump(constants_dict))
    #
    # print(f"New save '{save_name}' successfully created.")

if __name__ == '__main__':
	main()
