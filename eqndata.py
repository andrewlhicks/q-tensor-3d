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

def initialq_getufl():
    dict = load_yml(f'{pathl}setup/initial_q.yml')
    if dict is None:
        raise ValueError('Empty dict not accepted for "initial_q.yml".')
    if len(dict) > 1:
        raise KeyError('Only one dict key accepted for "initial_q.yml".')
    for key in dict.keys():
        if key not in ('ic_simple','ic_piecewise'):
            raise KeyError(f'Dict keys for "initial_q.yml" may only be "ic_simple" or "ic_piecewise", not "{key}".')

    if dict.keys() == ['ic_simple']:
        return dict['ic_simple'].uflfy()

    # If nothing is returned, assume ic_piecewise

    dict = dict['ic_piecewise']

    if set(dict.keys()) != {'if','then','else'}:
        raise KeyError('Keys for a piecewise ic must be "if", "then", and "else".')

    condition = dict['if']
    ufl_a = dict['then'].uflfy()
    ufl_b = dict['else'].uflfy()

    return f'conditional({condition},{ufl_a},{ufl_b})'


pathl = 'saves/experiment/'

constants = load_yml(f'{pathl}setup/constants.yml')
process_constants(constants)

settings = load_yml(f'{pathl}setup/settings.yml')

eqndata = {'constants':constants,'settings':settings}

print(initialq_getufl())

# dump_json(eqndata,f'{pathl}eqndata.json')
