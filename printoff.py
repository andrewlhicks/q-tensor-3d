""" Here, most of the printing is controlled. Most of it is done through
'plogging', i.e. print-logging, where the information is printed to the
console and then put into a log file. """

from time import sleep
import saves
import settings
import const
import functools

if settings.saves.save:
    from datetime import datetime
    from firedrake import COMM_WORLD
    if COMM_WORLD.rank == 0:
        mode = 'a' if settings.saves.mode == 'resume' else 'w'
        with open(f'{saves.current_directory}/log/log.txt',mode) as file:
            now = datetime.now()
            file.write(now.strftime('%c') + '\n')

# Decorators

def color_print(string,color=None):
    from firedrake.petsc import PETSc

    colors = {'header' : '\033[95m',
        'blue' : '\033[94m',
        'green' : '\033[92m',
        'warning' : '\033[93m',
        'fail' : '\033[91m',
        'end' : '\033[0m',
        'bold' : '\033[1m',
        'uline' : '\033[4m'}

    if color is not None:
        string = colors[color] + string + colors['end']

    PETSc.Sys.Print(string)

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

        if not settings.saves.save:
            def plog(string='',color=None):
                color_print(string,color)
            value = func(*args, **kwargs)
            return value

        with open(f'{saves.current_directory}/log/log.txt','a') as file:
            def plog(string='',color=None):
                file.write(string+'\n')
                color_print(string,color)
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

# Plogger functions

@plogger
def constants_info():
    """ Plogs information for the constants """

    plog()
    plog('CONSTANTS:',color='uline')

    plog()
    plog(f'     L1 = {const.L1},')
    plog(f'     L2 = {const.L2},')
    plog(f'     L3 = {const.L3},')
    plog(f'     q0 = {const.q0}')
    plog(f'      A = {const.A},')
    plog(f'      B = {const.B},')
    plog(f'      C = {const.C},')
    plog(f'     S0 = {const.S0},')
    plog(f'     W0 = {const.W0},')
    plog(f'     W1 = {const.W1},')
    plog(f'     W2 = {const.W2},')
    plog(f'epsilon = {const.ep},')
    plog(f'     L0 = {const.L0}')
    plog()

@plogger
def mesh_info():
    """ Prints the mesh information, namely, the mesh name and any other desired properties. """

    plog('MESH:',color='uline')

    # Assemble dicts

    dicts = []

    dicts.append({'title':'Mesh','text':settings.mesh.name})
    dicts.append({'title':'Refinements','text':settings.mesh.refs})

    # Print lines

    print_lines(*dicts)

@plogger
def options_info():
    plog('OPTIONS:',color='uline')

    # Assemble dicts

    dicts = []

    dicts.append({'title':'Visualize','text':settings.options.visualize})
    dicts.append({'title':'Manufactured','text':settings.options.manufactured})
    dicts.append({'title':'Weak boundary','text':settings.options.weak_boundary})
    dicts.append({'title':'Strong boundary','text':settings.options.strong_boundary})

    # Print lines

    print_lines(*dicts)

@plogger
def saves_info():
    plog('SAVES:',color='uline')

    # Assemble dicts

    dicts = []

    dicts.append({'title':'Save','text':settings.saves.save})
    if settings.saves.save:
        dicts.append({'title':'Mode','text':settings.saves.mode})
        dicts.append({'title':'Save directory','text':saves.current_directory})

    # Print lines

    print_lines(*dicts)

@plogger
def solver_info():
    plog('SOLVER:',color='uline')

    # Assemble dicts

    dicts = []

    dicts.append({'title':'Gadient descent','text':settings.solver.grad_desc})
    dicts.append({'title':'KSP Type','text':settings.solver.ksp_type})
    dicts.append({'title':'PC Type','text':settings.solver.pc_type})

    # Print lines

    print_lines(*dicts)

@plogger
def time_info():
    plog('TIME:',color='uline')

    # Assemble dicts

    dicts = []

    dicts.append({'title':'Time step','text':settings.time.step})
    # dicts.append({'title':'End time','text':settings.time.end})
    dicts.append({'title':'No. time steps','text':f'{settings.time.num}'})

    # Print lines

    print_lines(*dicts)

@plogger
def vis_info():
    plog('VIS:',color='uline')

    # Assemble dicts

    dicts = []

    dicts.append({'title':'Normal vector','text':settings.vis.normal})
    dicts.append({'title':'Save every','text':settings.vis.save_every})

    # Print lines

    print_lines(*dicts)

@plogger
def prelimCompTitle():
    """ Prints title for the preliminary computations. """

    sleep(1)

    plog(f'PRELIMINARY COMPUTATIONS:',color='uline')
    plog()

@plogger
def prelimCompInfo(time_elapsed):
    """ Prints information after the preliminary computations are finished. """

    plog(f'Finished preliminary computations in {time_elapsed:0.2f} seconds.')
    plog()

@plogger
def pdeSolveTitle():
    """ Prints title for the PDE solving. """

    sleep(1)

    plog(f'PDE SOLVE:',color='uline')
    plog()

@plogger
def pdeSolveInfo(**kwargs):
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
            dicts.append({'title':'Time elapsed','text':f'{kwargs[kw]:0.2f} seconds'})
        elif kw == 'refinement_level':
            dicts.append({'title':'Refinement level','text':kwargs[kw]})
        elif kw == 'energy':
            dicts.append({'title':'Energy','text':kwargs[kw]})
        elif kw == 'custom':
            dicts.append(kwargs[kw])

    print_lines(*dicts)

@plogger
def warning(text,spaced=True):
    """ Plogs a warning. """

    if not isinstance(text,str):
        raise TypeError('Warnings must be composed of a string.')
    plog(f'Warning: {text}',color='warning')
    if spaced: plog('')

@plogger
def info(text,spaced=True):
    if not isinstance(text,str):
        raise TypeError('Info must be composed of a string.')
    plog(text,color='green')
    if spaced: plog('')

@plogger
def text(text,spaced=True):
    if not isinstance(text,str):
        raise TypeError('Text must be composed of a string.')
    plog(text)
    if spaced: plog('')

# END OF CODE
