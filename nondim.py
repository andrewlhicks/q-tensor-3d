""" This is now the legacy version of the non-dimensionalization. Will be phased
out in favor of newsave.py. Converts a YAML file with constants to a new file
using their non-dimensionalized counterparts. """

import os
import sys
import yaml
from newsave import nondimensionalize

file_name = input("Enter constants file name: ")

# Open the file and load it using YAML

with open(f'constants/{file_name}.yml') as constants_file:
	constants_dict = yaml.load(constants_file, Loader=yaml.FullLoader)

# R = 3.0e-5 # Zoomed Slab for Lagerwall
R = 0.42e-6 # Lavrentovich shell radius

constants_dict = nondimensionalize(constants_dict,R)

if os.path.exists(f'constants/{file_name}_nd.yml'):
	while True:
		answer = input(f"Constants file '{file_name}_nd.yml' already exists. Overwrite? (y/n)")
		if answer in ('y','Y'):
			break
		if answer in ('n','N'):
			print("No file overwritten.")
			sys.exit()

with open(f'constants/{file_name}_nd.yml','w') as constants_file:
	constants_file.write(yaml.dump(constants_dict))

print(f"Successfully created '{file_name}_nd.yml'")
