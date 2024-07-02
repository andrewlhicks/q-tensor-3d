""" Takes the setup data given and compiles the file uflcache.json, which
controls the PDE to be solved. """

import sys
from q3d.loaddump import *
from q3d.uflcache.userexpr import *
from q3d.uflplus import *
from ufl.core.expr import Expr
import re

add_sympy_constructors()
add_ufl_constructors()

userexpr_types = ('initcond','w_bdy_nu','s_bdy','manu_q','forcing_f','forcing_g')

def usage():
    usage = """usage: python uflcache.py (-l | -r) <save-name>
  l: file in './saves'
  r: file in './saves-remote'"""
    print(usage)

def process_condition(condition):
    # if not dictionary, return simple condition
    if not isinstance(condition,dict):
        return simple_condition(condition)
    
    # if dictionary, then go through rigorous checking process
    if not condition:
        raise ValueError('Empty dict not accepted for user expression.')
    if len(condition) > 1:
        raise KeyError('Only one dict key accepted for user expression.')
    for key in condition.keys():
        if key not in ('simple','piecewise'):
            raise KeyError(f'Dict keys for user expression may only be "simple" or "piecewise", not "{key}".')

    # if 'simple' key is given, return simple condition
    if 'simple' in condition.keys():
        return simple_condition(condition['simple'])

    # If nothing is returned, return piecewise
    return piecewise_condition(condition['piecewise'])

def simple_condition(condition: Expr | FromSympy) -> str:
    if isinstance(condition, Expr):
        string = repr(condition)
        string = re.sub(r'<ufl.domain.AbstractDomain object at 0[xX][0-9a-fA-F]+?>', 'mesh', string)
        return string
    elif isinstance(condition, FromSympy):
        return condition.uflfy()
    else:
        raise TypeError('Conditions must be an instance of ufl.core.expr.Expr or FromSympy')

def piecewise_condition(condition):
    if set(condition.keys()) != {'if','then','else'}:
        raise KeyError('Keys for a piecewise condition must be "if", "then", and "else".')

    condition = condition['if']
    ufl_a = condition['then'].uflfy()
    ufl_b = condition['else'].uflfy()

    return f'conditional({condition},{ufl_a},{ufl_b})'

def load_userexpr_yml(path_head):
    try:
        userexpr_dict = load_yml(f'{path_head}/userexpr.yml')
    except FileNotFoundError:
        print(f'File "{path_head}/userexpr.yml" not found, please provide.')
        sys.exit()
    for key in userexpr_dict.keys():
        if key not in userexpr_types:
            types = '("' + '","'.join(userexpr_types) + '")'
            raise KeyError(f'Dict keys for user expression may only be {types}, not "{key}".')
    return userexpr_dict

def dump_uflcache_json(userexpr_dict,path_head):
    uflcache = {key:process_condition(val) for key, val in userexpr_dict.items()}
    dump_json(uflcache,f'{path_head}/uflcache.json')

def build_uflcache(path_head):
    userexpr_dict = load_userexpr_yml(path_head)
    dump_uflcache_json(userexpr_dict,path_head)

def main():
    import getopt
    import sys
    import q3d.config as config

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'l:r:', ['help'])
    except getopt.GetoptError as err:
        print(err)  # will print something like "option -a not recognized"
        sys.exit()

    for o, a in opts:
        if o in ('--help'):
            usage()
            sys.exit()
        elif o in ('-l'):
            remote = ''
            save_name = a
        elif o in ('-r'):
            remote = '-remote'
            save_name = a
        else:
            assert False, "unhandled option"

    try:
        save_name
    except NameError:
        print("Must choose local (-l) or remote (-r).")
        sys.exit()
    
    path_head = f'saves{remote}/{save_name}'

    config.initialize(f'{path_head}/settings.yml',f'{path_head}/constants.yml') # So userexpr.py has constants.S0

    build_uflcache(path_head)
    print(f"Done writing uflcache.json at '{path_head}'.")

if __name__ == '__main__':
    main()
