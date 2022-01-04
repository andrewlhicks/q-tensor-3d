""" Takes the setup data given and compiles the file uflcache.json, which
controls the PDE to be solved. """

from loaddump import *
from userexpression import *

def usage():
    print('')

def initcond_getufl(dict):
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

def main():
    import config
    config.initialize('saves/ver/settings.yml','saves/ver/constants.yml')
    userexpr = load_yml('saves/ver/userexpr.yml')
    initcond = userexpr['initcond']
    uflcache = {'initcond':initcond_getufl(initcond)}
    dump_json(uflcache,f'saves/ver/uflcache.json')

if __name__ == '__main__':
    main()
