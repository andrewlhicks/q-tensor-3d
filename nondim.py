""" Converts a YAML file with constants to a new file using their non-dimensionalized counterparts """

import yaml
from math import sqrt
from misc import getValues

file_path = input('Enter constants file name: ')

# Open the file and load it using YAML

with open(f'constants/{file_path}.yml') as constants_file:
	constants_dict = yaml.load(constants_file, Loader=yaml.FullLoader)

const = constants_dict
nd_const = {}

# Define auxiliary constants

R = 3.0e-5 # Zoomed Slab for Lagerwall
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

# Replace old constants with new ones in the constants_dict

constants_dict = nd_const

with open(f'constants/{file_path}_nd.yml','w') as constants_file:
	constants_file.write(yaml.dump(constants_dict))
