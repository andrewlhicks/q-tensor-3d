""" Converts a YAML file with constants to a new file using their
non-dimensionalized counterparts. """

import os
import sys
import getopt
import yaml
from newsave import nondimensionalize

def usage():
	print("python nondim.py <constants-file>.yml <radius> [-s]")

def main():
	save_mode = False

	try:
		opts, args = getopt.getopt(sys.argv[3:], 'hs', [''])
	except getopt.GetoptError as err:
		# print help information and exit:
		print(err)  # will print something like "option -a not recognized"
		usage()
		sys.exit()

	for o, a in opts:
		if o in ('-h'):
			usage()
			sys.exit()
		if o in ('-s'):
			save_mode = True
		else:
			assert False, "unhandled option"

	file = sys.argv[1]
	file_name = file.split('.yml')[0]
	R = float(sys.argv[2])

	if not file.endswith('.yml'):
		print("Constants file must end in '.yml'.")
		sys.exit()
	if not R > 0:
		print("Radius must be positive.")
		sys.exit()

	with open(f'constants/{file}') as constants_file:
		constants_dict = yaml.load(constants_file, Loader=yaml.Loader)
	constants_dict = nondimensionalize(constants_dict,R)

	print("")
	print(yaml.dump(constants_dict))

	if save_mode:
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
		print(f"File successfully saved as '{file_name}_nd.yml'.")

if __name__ == '__main__':
	main()
