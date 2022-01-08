""" Takes the setup data given and compiles the file uflcache.json, which
controls the PDE to be solved. """

from loaddump import *
from userexpr import *

def usage():
    usage = """usage: python uflcache.py (-l | -r) <save-name>
  l: file in './saves'
  r: file in './saves-remote'"""
    print(usage)

def process_initcond(dict):
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

def load_userexpr(path_head):
    try:
        userexpr_dict = load_yml(f'{path_head}/userexpr.yml')
    except FileNotFoundError:
        print(f'File "{path_head}/userexpr.yml" not found, falling back to default.')
        userexpr_dict = load_yml('defaults/userexpr.yml')
    return userexpr_dict

def dump_uflcache(userexpr_dict,path_head):
    initcond = userexpr_dict['initcond']
    uflcache = {'initcond':process_initcond(initcond)}
    dump_json(uflcache,f'{path_head}/uflcache.json')

def build_uflcache(path_head):
    userexpr_dict = load_userexpr(path_head)
    dump_uflcache(userexpr_dict,path_head)

def main():
    import getopt
    import sys
    import config

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
