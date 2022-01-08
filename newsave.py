import getopt
import os
import sys
import yaml
from loaddump import *
from userexpr import * # needed to process userexpr.yml

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
    """ Takes specified settings and constants files, along with potentially
    a nondimensionalization length, and creates a new save. Should be
    eventually deprecated in favor of a default build. """
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

    # load settings, constants, and userexpr dicts
    settings_dict = load_yml(settings_path)
    constants_dict = load_yml(constants_path)
    userexpr_txt = load_txt('defaults/userexpr.yml') #  <---- this is to be the behavior of settings/constants in the future

    if radius is not None:
        print(f"Creating non-dimensionalized version of '{constants_name}' with radius {radius}")
        constants_dict = nondimensionalize(constants_dict,radius)

    create_save(save_name,save_path,settings_dict,constants_dict,userexpr_txt)

def copy():
    argv = sys.argv[2:]

    if len(argv) != 2:
        print("2 args required for -c.")
        sys.exit()

    # process args
    old_save_name = argv[0]
    save_name = argv[1]
    save_path = 'saves/' + save_name

    # load settings, constants, and userexpr dicts
    settings_dict = load_yml(f'saves/{old_save_name}/settings.yml')
    constants_dict = load_yml(f'saves/{old_save_name}/constants.yml')
    userexpr_txt = load_txt(f'saves/{old_save_name}/userexpr.yml')

    create_save(save_name,save_path,settings_dict,constants_dict,userexpr_txt)

def create_save(save_name,save_path,settings_dict,constants_dict,userexpr_txt):
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

    # dump yaml files into save
    dump_yml(settings_dict,f'{save_path}/settings.yml')
    dump_yml(constants_dict,f'{save_path}/constants.yml')
    dump_txt(userexpr_txt,f'{save_path}/userexpr.yml')

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

if __name__ == '__main__':
	main()
