import getopt
import os
import sys
from loaddump import *
from userexpr import * # needed to process userexpr.yml

def usage():
    usage_str = """usage: python newsave.py (-b | -c) <save_name>
  b: builds new save <save_name> from yml's in ./defaults
  c: copies yml's from <save_name> to new save
"""
    print(usage_str)

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

def build(save_name:str):
    """ Builds new save <save_name> from yml's in ./defaults. """

    # load settings, constants, and userexpr text
    settings_txt = load_txt('defaults/settings.yml')
    constants_txt = load_txt('defaults/constants.yml')
    userexpr_txt = load_txt('defaults/userexpr.yml')

    # create save
    create_save(save_name,settings_txt,constants_txt,userexpr_txt)

def copy(save_name:str):
    """ Copies yml's from <save_name> to new save. """

    # load settings, constants, and userexpr text
    try:
        settings_txt = load_txt(f'saves/{save_name}/settings.yml')
    except FileNotFoundError:
        print('No settings file found. Reverting to default...')
        settings_txt = load_txt('defaults/settings.yml')
    try:
        constants_txt = load_txt(f'saves/{save_name}/constants.yml')
    except FileNotFoundError:
        print('No constants file found. Reverting to default...')
        constants_txt = load_txt('defaults/constants.yml')
    try:
        userexpr_txt = load_txt(f'saves/{save_name}/userexpr.yml')
    except FileNotFoundError:
        print('No userepxr file found. Reverting to default...')
        userexpr_txt = load_txt('defaults/userexpr.yml')

    create_save(save_name,settings_txt,constants_txt,userexpr_txt)

def create_save(save_name,settings_txt,constants_txt,userexpr_txt):
    """ Creates a new save <save_name> (unless there is a naming conflict) with
    specifed settings, constants, and userexpr text. """

    # get path
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

    # dump yaml files into save
    dump_txt(settings_txt,f'{save_path}/settings.yml')
    dump_txt(constants_txt,f'{save_path}/constants.yml')
    dump_txt(userexpr_txt,f'{save_path}/userexpr.yml')

    print(f"New save '{save_name}' successfully created.")

def main():
    if not os.path.exists('saves'):
        os.makedirs('saves')

    try:
        opts, args = getopt.getopt(sys.argv[1:],'b:c:',['help'])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)
        print('use --help for usage')
        sys.exit()

    for o, a in opts:
        if o in ('-b'):
            build(a)
        elif o in ('-c'):
            copy(a)
        elif o in ('--help'):
            usage()
        else:
            sys.exit()

if __name__ == '__main__':
	main()
