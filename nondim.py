""" Converts a YAML file with constants to a new file using their non-dimensionalized counterparts """

import yaml
from math import sqrt
from misc import getValues

settings_filename = input('Enter settings file name: ')

# Open the file and load it using YAML

with open(f'settings/{settings_filename}.yml') as settings_file:
	settings_dict = yaml.load(settings_file, Loader=yaml.FullLoader)

const = settings_dict['const']
nd_const = {}

# Define auxiliary constants

R = float(input('Enter length scale R: '))
AM = max(const['A'],const['B'],const['C'])
# AM = const['A']
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

# Replace old constants with new ones in the settings_dict

settings_dict['const'] = nd_const

with open(f'settings/{settings_filename}_nd.yml','w') as settings_file:
	settings_file.write(yaml.dump(settings_dict))