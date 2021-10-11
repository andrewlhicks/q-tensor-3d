""" Takes the setup data given and compiles the file eqndata.json, which
controls the PDE to be solved. """

import json
import yaml
from sympyplus import *

def usage():
    print('')

def load_json(path):
    with open(path) as file:
        dict = json.loads(file.read())
    return dict

def dump_json(dump,path):
    with open(path,'w') as file:
        file.write(json.dumps(dump))

def load_yml(path):
    with open(path) as file:
        dict = yaml.load(file, Loader=yaml.Loader)
    return dict

def dump_yml(dump,path):
    with open(path,'w') as file:
        file.write(yaml.dump(dump))

def process_constants(constants):
    from math import sqrt

    A = constants['A']
    B = constants['B']
    C = constants['C']

    if constants['L0'] == 'auto':
        constants['L0'] = ceil(2*(A+B**2/C))

    constants['S0'] = (B + sqrt(B**2 + 24.0*A*C))/(4.0*C)

pathl = 'saves/experiment/'

constants = load_yml(f'{pathl}setup/constants.yml')
process_constants(constants)

settings = load_yml(f'{pathl}setup/settings.yml')

eqndata = {'constants':constants,'settings':settings}

A = load_yml(f'{pathl}setup/initial_q.yml')
print(A)

# dump_json(eqndata,f'{pathl}eqndata.json')
