""" Takes the setup data given and compiles the file eqndata.json, which
controls the PDE to be solved. """

from loaddump import *

def usage():
    print('')

def process_constants(constants):
    from math import sqrt

    A = constants['A']
    B = constants['B']
    C = constants['C']

    if constants['L0'] == 'auto':
        constants['L0'] = ceil(2*(A+B**2/C))

    constants['S0'] = (B + sqrt(B**2 + 24.0*A*C))/(4.0*C)

def initcond_getufl(path):
    dict = load_yml(path)
    if dict is None:
        raise ValueError('Empty dict not accepted for user expression.')
    if len(dict) > 1:
        raise KeyError('Only one dict key accepted for user expression.')
    for key in dict.keys():
        if key not in ('simple','piecewise'):
            raise KeyError(f'Dict keys for user expression may only be "simple" or "piecewise", not "{key}".')

    if 'simple' in dict.keys():
        return dict['simple'].uflfy()

    # If nothing is returned, assume ic_piecewise

    dict = dict['piecewise']

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
dump_json(eqndata,f'{pathl}eqndata.json')

from userexpression import *

print(initcond_getufl(f'{pathl}setup/initial_q.yml'))

# dump_json(eqndata,f'{pathl}eqndata.json')
