""" Here, most of the printing is controlled. Most of it is done through
'plogging', i.e. print-logging, where the information is printed to the
console and then put into a log file.

My intent is to deprecate this in favor of using the built-in logging
module, but this will take time to migrate. """

from multiprocessing.sharedctypes import Value
import saves
import functools

def main():
    if saves.SaveMode:
        from datetime import datetime
        from firedrake import COMM_WORLD
        if COMM_WORLD.rank == 0:
            mode = 'a' if saves.SaveMode == 'resume' else 'w'
            with open(f'{saves.current_directory}/log.txt',mode) as file:
                now = datetime.now()
                file.write(now.strftime('%c') + '\n')
                file.write('\n')
                Print(now.strftime('%c'))
                Print()

# Decorators

def Print(string='', color=None, *args, **kwargs):
    from firedrake.petsc import PETSc

    colors = {'header' : '\033[95m%s\033[0m',
        'blue' : '\033[94m%s\033[0m',
        'green' : '\033[92m%s\033[0m',
        'warning' : '\033[93m%s\033[0m',
        'fail' : '\033[91m%s\033[0m',
        'bold' : '\033[1m%s\033[0m',
        'uline' : '\033[4m%s\033[0m'}

    if color is not None:
        string = colors[color] % string

    PETSc.Sys.Print(string, *args, **kwargs)

def plogger(func):
    """ A decorator that defines a plog function every time it is called. The
    reason for this is to ensure that the 'with' statement is used correctly.
    """
    @functools.wraps(func)
    def wrapper_plogger(*args, **kwargs):
        from firedrake import COMM_WORLD
        global plog

        if COMM_WORLD.rank != 0:
            return

        if not saves.SaveMode:
            def plog(string='', color=None, *args, **kwargs):
                Print(string, color, *args, **kwargs)
            value = func(*args, **kwargs)
            return value

        with open(f'{saves.current_directory}/log.txt','a') as file:
            def plog(string='', color=None, *args, **kwargs):
                file.write(string+'\n')
                Print(string, color, *args, **kwargs)
            value = func(*args, **kwargs)
        return value

    return wrapper_plogger

# Functions that print lines

def print_lines(*args):
    """ Prints a lines, each line containing a tile and text. Args are
    dictionaries with 'title' and 'text' attributes. """

    def print_line(title,text):
        if not isinstance(title,str):
            raise TypeError('Title must be a string.')

        title_len = len(title)
        spaces_len = int(indent_len - 1 - title_len)
        if spaces_len < 0: spaces_len = 0

        spaces = ''.join([' ' for _ in range(spaces_len)])
        plog(f'{title}:{spaces}{text}')

    indent_len = int(0)

    for arg in args:
        if len(arg['title']) + 2 > indent_len:
            indent_len = len(arg['title']) + 2

    plog()
    for arg in args:
        print_line(arg['title'],arg['text'])
    plog()

@plogger
def iter_info(*strings: str, i: int, p: str='-', b: str='()', show_numbering: bool=True):
    for string in strings:
        if not isinstance(string,str):
            raise TypeError('Must be str')
    if len(b) != 2:
        raise ValueError('Must be str of len 2')
    
    plog(''.join([p for _ in range(15)]))

    s = ''.join([' ' for _ in range(len(str(i)))])
    if show_numbering:
        plog(f'{p} {b[0]}{i}{b[1]} ' + strings[0])
    else:
        plog(f'{p}  {s}  ' + strings[0])
    for string in strings[1:]:
        plog(f'{p}  {s}  ' + string)
    
    plog(''.join([p for _ in range(15)]))

# Plogger functions

@plogger
def constants_info():
    from config import constants as c
    plog()
    plog('CONSTANTS:',color='uline')
    plog()
    for key, val in c.as_dict().items():
        plog(f'{key} = {val}')
    plog()

@plogger
def settings_info():
    from config import settings

    plog('SETTINGS:',color='uline')
    plog()
    for key, val in settings.as_dict().items():
        plog(f'{key}:')
        for key, val in val.items():
            plog(f'  {key} : {val}')
    plog()

@plogger
def pde_solve_info(**kwargs):
    """ Prints desired information for the PDE solving after it is complete. """

    # Assemble dicts

    dicts = []

    for kw in kwargs:
        if kw == 'mesh_numnodes':
            dicts.append({'title':'Mesh node struc','text':f'{kwargs[kw]} x {kwargs[kw]} x {kwargs[kw]}'})
        elif kw == 'h1_error':
            dicts.append({'title':'H1 error','text':f'{kwargs[kw]:0.15f}'})
        elif kw == 'l2_error':
            dicts.append({'title':'L2 error','text':f'{kwargs[kw]:0.15f}'})
        elif kw == 'time_elapsed':
            dicts.append({'title':'Time elapsed','text':kwargs[kw]})
        elif kw == 'refinement_level':
            dicts.append({'title':'Refinement level','text':kwargs[kw]})
        elif kw == 'energy':
            dicts.append({'title':'Energy','text':kwargs[kw]})
        elif kw == 'custom':
            dicts.append(kwargs[kw])

    print_lines(*dicts)

@plogger
def info(text, spaced=False, color=None, *args, **kwargs):
    if not isinstance(text,str):
        raise TypeError('Text must be composed of a string.')
    plog(text, color=color, *args, **kwargs)
    if spaced: plog('')

def green(text, spaced=False, *args, **kwargs):
    info(text, spaced=spaced, color='green', *args, **kwargs)

def blue(text, spaced=False, *args, **kwargs):
    info(text, spaced=spaced, color='blue', *args, **kwargs)

def warning(text, spaced=False, *args, **kwargs):
    info(f'warning: {text}', spaced=spaced, color='warning', *args, **kwargs)

def fail(text, spaced=False, *args, **kwargs):
    info(f'fatal: {text}', spaced=spaced, color='fail', *args, **kwargs)

def sinfo(text, color=None, *args, **kwargs):
    info(text, spaced=True, color=color, *args, **kwargs)

def sgreen(text, *args, **kwargs):
    green(text, spaced=True, *args, **kwargs)

def sblue(text, *args, **kwargs):
    blue(text, spaced=True, *args, **kwargs)

def swarning(text, *args, **kwargs):
    warning(text, spaced=True, *args, **kwargs)

def sfail(text, *args, **kwargs):
    fail(text, spaced=True, *args, **kwargs)

# MAIN

main()

# END OF CODE
