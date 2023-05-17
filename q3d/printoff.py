""" This module prints and logs information. This is done through
'plogging', i.e. print-logging, where the information is printed to the
console and then put into a log file. """

import functools
from datetime import datetime

from firedrake import COMM_WORLD
from firedrake.petsc import PETSc

import q3d.saves as saves

# utility functions

def rstr(string,number):
    return ''.join([string for _ in range(number)])

def Print(string='', *, color=None, **kwargs):
    colors = {'header' : '\033[95m%s\033[0m',
        'blue' : '\033[94m%s\033[0m',
        'green' : '\033[92m%s\033[0m',
        'warning' : '\033[93m%s\033[0m',
        'fail' : '\033[91m%s\033[0m',
        'bold' : '\033[1m%s\033[0m',
        'uline' : '\033[4m%s\033[0m'}

    if color is not None:
        string = colors[color] % string

    PETSc.Sys.Print(string, **kwargs)

# decorators

def plogger(func):
    """ A decorator that defines a plog function every time it is called. The
    reason for this is to ensure that the 'with' statement is used correctly.
    """
    @functools.wraps(func)
    def wrapper_plogger(*args, **kwargs):
        # exit function if not on processor 0
        if COMM_WORLD.rank != 0:
            return

        # must make plog global to reach @plogger functions' scopes
        global plog

        # allow user to specify write mode, otherwise give 'a'
        mode = kwargs.pop('mode', 'a')

        # if SaveMode not specified, plog is the same as Print
        if not saves.SaveMode:
            def plog(string='', *, color=None, **plog_kwargs):
                Print(string, color=color, **plog_kwargs)
            value = func(*args, **kwargs)
            del plog # important to ensure plog is not accessed by non-plogger functions
            return value

        # if SaveMode is specified, plog will Print and write to log file
        with open(f'{saves.SavePath}/log.txt', mode) as file:
            def plog(string='', *, color=None, **plog_kwargs):
                file.write(string+'\n')
                Print(string, color=color, **plog_kwargs)
            value = func(*args, **kwargs)
            del plog # important to ensure plog is not accessed by non-plogger functions
        return value
    return wrapper_plogger

# functions that print lines

def print_lines(*args):
    """ Prints a lines, each line containing a tile and text. Args are
    dictionaries with 'title' and 'text' attributes. MUST be called
    within a @plogger function. """

    try:
        plog
    except NameError:
        fail('must be called within @plogger function')
        raise

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

@plogger
def iter_info_verbose(*strings: str, i: int, j: int=None, b: str='()', spaced=False):
    now = datetime.now().strftime('%c')
    for string in strings:
        if not isinstance(string,str):
            raise TypeError('Must be str')
    if len(b) != 2:
        raise ValueError('Must be str of len 2')
    if j is None:
        j = ''
    else:
        j = f'-{j}'
    
    now_str = f'[{now}]'
    it_str = f'{b[0]}{i}{j}{b[1]}'
    spacing = rstr(' ',len(now_str + it_str) + 2)
    plog(f'{now_str} {it_str} {strings[0]}')
    for string in strings[1:]:
        plog(spacing + string)
    
    if spaced: plog('')


# plogger functions

@plogger
def constants_info():
    from q3d.config import constants as c
    plog()
    plog('CONSTANTS:',color='uline')
    plog()
    for key, val in c.as_dict().items():
        plog(f'{key} = {val}')
    plog()

@plogger
def settings_info():
    from q3d.config import settings

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
def text(string, *, spaced=False, color=None, **kwargs):
    if not isinstance(string,str):
        raise TypeError('Text must be composed of a string.')
    plog(string, color=color, **kwargs)
    if spaced: plog('')

def stext(string, *, color=None, **kwargs):
    text(string, spaced=True, color=color, **kwargs)

@plogger
def info(string, *, spaced=False, color=None, **kwargs):
    now = datetime.now().strftime('%c')
    if not isinstance(string,str):
        raise TypeError('Text must be composed of a string.')
    plog(f'[{now}] {string}', color=color, **kwargs)
    if spaced: plog('')

def green(string, *, spaced=False, **kwargs):
    info(string, spaced=spaced, color='green', **kwargs)

def blue(string, *, spaced=False, **kwargs):
    info(string, spaced=spaced, color='blue', **kwargs)

def warning(string, *, spaced=False, **kwargs):
    info(f'warning: {string}', spaced=spaced, color='warning', **kwargs)

def fail(string, *, spaced=False, **kwargs):
    info(f'fatal: {string}', spaced=spaced, color='fail', **kwargs)

def sinfo(string, *, color=None, **kwargs):
    info(string, spaced=True, color=color, **kwargs)

def sgreen(string, **kwargs):
    green(string, spaced=True, **kwargs)

def sblue(string, **kwargs):
    blue(string, spaced=True, **kwargs)

def swarning(string, **kwargs):
    warning(string, spaced=True, **kwargs)

def sfail(string, **kwargs):
    fail(string, spaced=True, **kwargs)

# END OF CODE
