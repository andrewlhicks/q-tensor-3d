from misc import getValues
import yaml
from numpy import sqrt
from math import ceil
import settings

def _load_file(file_path):
	# Open the file and load it using YAML

	with open(f'constants/{file_path}') as constants_file:
		constants_dict = yaml.load(constants_file, Loader=yaml.FullLoader)

	global L1, L2, L3, q0, A, B, C, ep, W0, W1, W2, beta, L0, S0, dt

	L1, L2, L3, q0, A, B, C, ep, W0, W1, W2, beta = getValues(constants_dict,'L1, L2, L3, q0, A, B, C, ep, W0, W1, W2, beta')
	L0 = ceil(2*(A+B**2/C)) if constants_dict['L0'] == 'auto' else constants_dict['L0']
	S0 = (B + sqrt(B**2 + 24.0*A*C))/(4.0*C)
	dt = settings.time.step
