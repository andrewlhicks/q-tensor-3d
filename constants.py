from misc import getValues
import yaml

def _load_file(file_path):
	# Open the file and load it using YAML

	with open(f'constants/{file_path}') as constants_file:
		constants_dict = yaml.load(constants_file, Loader=yaml.FullLoader)

	global const

	class const:
		from numpy import sqrt
		from math import ceil
		L1, L2, L3, q0, A, B, C, ep, W0, W1, W2 = getValues(constants_dict['const'],'L1, L2, L3, q0, A, B, C, ep, W0, W1, W2')
		L0 = ceil(2*(A+B**2/C)) if constants_dict['const']['L0'] == 'auto' else constants_dict['const']['L0']
		S0 = (B + sqrt(B**2 + 24.0*A*C))/(4.0*C)

	# Check if L1, L2, and L3 are set properly

	if 0>=const.L1 or -const.L1>=const.L3 or const.L3>=2*const.L1 or -3/5*const.L1-1/10*const.L3>=const.L2:
		print('Warning: L1, L2, and L3 do not satisfy the proper inequalities')