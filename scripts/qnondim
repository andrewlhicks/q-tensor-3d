""" Converts a YAML file with constants to a new file using their
non-dimensionalized counterparts. """

import os
import sys
import getopt
import yaml

usage_str = """
usage: python nondim.py [-s] <constants-file>.yml <radius>
"""

def usage():
	print(usage_str)

def nondimensionalize(const,R):
	""" For a constants dictionary const and a length R (for a sphere-like mesh,
	the radius), returns the non-dimensionalized version of these constants """

	from math import sqrt

	nd_const = {}

	# Define auxiliary constants

	S0 = (const['B'] + sqrt(const['B']**2 + 24*const['A']*const['C']))/(4*const['C'])
	C0 = - const['A']/3*S0**2 - 2*const['B']/27*S0**3 + const['C']/9*S0**4
	# AM = max(const['A'],const['B'],const['C'])
	LM = max(const['L1'],const['L2'],const['L3'])

	# Define non-dimensionalized constants

	nd_const['L0'] = const['L0']
	nd_const['A'] = const['A']/C0
	nd_const['B'] = const['B']/C0
	nd_const['C'] = const['C']/C0
	nd_const['L1'] = const['L1']/LM
	nd_const['L2'] = const['L2']/LM
	nd_const['L3'] = const['L3']/LM
	nd_const['W0'] = const['W0']*R/(LM)
	nd_const['W1'] = const['W1']*R/(LM)
	nd_const['W2'] = const['W2']*R/(LM)
	nd_const['q0'] = const['q0']*R
	nd_const['ep'] = sqrt(LM/(C0*R**2))
	try:
		nd_const['beta'] = const['beta']
	except:
		nd_const['beta'] = 1

	# Return the non-dimensionalized constants

	return nd_const

def main():
	save_mode = False

	try:
		opts, args = getopt.getopt(sys.argv[1:], 's', ['help'])
	except getopt.GetoptError as err:
		# print help information and exit:
		print(err)  # will print something like "option -a not recognized"
		usage()
		sys.exit()

	for o, a in opts:
		if o in ('--help'):
			usage()
			sys.exit()
		if o in ('-s'):
			save_mode = True
		else:
			assert False, "unhandled option"

	try:
		float(sys.argv[-1])
	except:
		print("Radius must be float.")
		sys.exit()
	R = float(sys.argv[-1])
	if not R > 0:
		print("Radius must be positive.")
		sys.exit()

	file = sys.argv[-2]
	if not file.endswith('.yml'):
		print("Constants file must end in '.yml'.")
		sys.exit()
	file_name = file.split('.yml')[0]


	with open(f'constants/{file}') as constants_file:
		constants_dict = yaml.load(constants_file, Loader=yaml.Loader)
	constants_dict = nondimensionalize(constants_dict,R)

	print("")
	print(yaml.dump(constants_dict))

	if save_mode:
		if os.path.exists(f'constants/{file_name}_nd.yml'):
			while True:
				answer = input(f"Constants file '{file_name}_nd.yml' already exists. Overwrite? (y/n) ")
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
